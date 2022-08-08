/** **************************************************************************** 
*
* @name       POD
* @brief      Computation/manipulation of POD basis
*             independently of problem dimension
*
* @author     Swetlana Giere
* @date       08.03.2017 (start of implementation)
*
*******************************************************************************/

#include <algorithm>
#include <cmath>
#include <POD.h>
#include <SnapshotsCollector.h>
#include <Solver.h>

#ifdef _MPI
#include <mpi.h>
#endif

// Extern declaration of Lapack method.
extern "C" {
//  computes the eigenvalues and eigenvectors
void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w,
            double* work, int* lwork, int* info);
}

/** ***************************************************************************/
ParameterDatabase POD::default_pod_database()
{
  ParameterDatabase db("Default ParMooN parameter database for "
                       "POD-basis computation");

  db.add("pod_directory", ".",
         "This directory is where the POD basis and Co. are written. This "
         "directory will be created, if it does not exist already. Files in "
         "this directory will be overwritten without any warning.");

  db.add("pod_basename", "parmoon_pod",
         "Basename for pod basis and related files. When writing the POD "
         "basis, the basis elements will be written into pod_basename.pod, "
         "the average (if needed) into pod_basename.mean. When reading the "
         "basis, the program expect to find files ending with .pod and .mean");

  db.add("pod_rank", (size_t) 0,
         "This integer specifies the dimension of the POD space to be computed."
         " If pod_rank <= 0, then all possible POD modes will be computed.");

  db.add("pod_fluctuations_only", true,
         "This is the flag whether the POD basis should be computed only "
         "fluctuation part (without average) of the snapshots "
         "(also central trajectory method).",
         {true,false});

  db.add("pod_inner_product", "euclidean",
         "Specification of the inner product which is used to compute POD "
         "basis. Besides default value, 'L2' (and 'H1' for Navier Stokes) are "
         "possible at the moment.",
         {"euclidean", "L2", "H1"});

  // Merge with snapshots database needed to get the path of the snapshots
  db.merge(SnapshotsCollector::default_snapshots_database(), true);

  return db;
}

/** ***************************************************************************/
POD::POD(const ParameterDatabase& param_db, std::string filename_suffix) :
  db(POD::default_pod_database()),
  rank(0),
  length(0),
  length_snaps(0),
  number_snaps(0),
  eigen_threshold(1.0e-10),
  valid_eigs(0),
  eigs(nullptr),
  gramian_ptr(nullptr),
  projection_solver_db(Solver<>::default_solver_database()),
  snaps_mat(0, 0),
  snaps_mean(0),
  pod_basis(0, 0)
{
  this->db.merge(param_db, false);
  POD::check_parameters();
  rank = db["pod_rank"];

  set_basename(filename_suffix);
}

/** ***************************************************************************/
void POD::set_basename(std::string filename_suffix)
{
  this->filename_suffix = filename_suffix;
  this->basename =  this->db["pod_directory"].get<std::string>() + "/"
                  + this->db["pod_basename"].get<std::string>() + ".";

  if (!filename_suffix.empty())
  {
    this->basename += filename_suffix + ".";
  }

  this->basename += this->db["pod_inner_product"].get<std::string>() + ".";

  if (this->db["pod_fluctuations_only"])
  {
    basename += "fluc.";
  }
}

/** ***************************************************************************/
POD::~POD()
{
  if (eigs != nullptr)
  {
    delete[] eigs;
  }
}

/** ***************************************************************************/
void POD::check_parameters()
{
  if (!this->db["pod_inner_product"].is("L2")
    && !this->db["pod_inner_product"].is("H1")
    && !this->db["pod_inner_product"].is("euclidean"))
  {
    ErrThrow("Error: POD::check_parameters(), the parameter ",
             "'pod_inner_product' does not currently support the value ",
             this->db["pod_inner_product"].get<std::string>(), ".");
  }
}

/** ***************************************************************************/
void POD::set_gramian_ptr(std::shared_ptr<TMatrix> mat)
{
  this->gramian_ptr = mat;
}

void POD::set_projection(const BlockFEMatrix& mat, const ParameterDatabase& solver_db)
{
  projection_matrix = mat;
  projection_solver_db.merge(solver_db, false);
}

/** ***************************************************************************/
void POD::read_snapshots()
{
  if (length <= 0)
  {
    ErrThrow("POD read_snapshots(): length is not available!");
  }

  // read snapshots from file
  bool bin = false;

  std::string snap_filename = SnapshotsCollector::get_snapshot_filename(this->db, this->filename_suffix);

  if (db["write_solution_binary"])
  {
    bin = true;
  }

  Output::print<1>("Reading snapshots from file: ", snap_filename, "...");

  std::vector < std::vector<double> > tmp_snaps;
  SnapshotsCollector::read_data(snap_filename, bin, length, tmp_snaps);

  this->length_snaps = tmp_snaps[0].size(); //total dof
  this->number_snaps = tmp_snaps.size();

  // resize snaps_mat
  DenseMatrix snaps_mat(this->number_snaps, this->length_snaps);
  this->snaps_mat = snaps_mat;

  // store snapshots into member matrix snaps_mat 
  for (int i = 0; i < this->length_snaps; i++)
  {
    for (int j = 0; j < this->number_snaps; j++)
    {
      this->snaps_mat.setEntry(j, i, tmp_snaps[j][i]);
    }
  }

  Output::print<1>("... snapshots read.");
  Output::print<1>("  * Length of snapshots : ", this->length_snaps);
  Output::print<1>("  * Number of snapshots : ", this->number_snaps);
}

/** ***************************************************************************/
void POD::decompose_snaps()
{
  // transforms raw snaps into mean + fluctuations.
  // this is allowed to be MPI-agnostic because there's no interaction

  this->snaps_mean.resize(this->length_snaps);

  if (this->length_snaps == 0)
  {
    ErrThrow("Error: Snapshots are not available. Before computing POD basis "
             "read_snapshots() has to be called!");
  }

  bool with_time_derivatives = db["snaps_time_derivative"];

  if (with_time_derivatives)
  {
    Output::root_info("POD", "Snapshots include time derivatives - averaging every second snapshot only.");
  }

  for (int i = 0; i < length_snaps; i++)
  {
    // calculate sum of the ith column of snapshot matrix
    double sum_i = 0.0;

    if (with_time_derivatives)
    {
      // no time derivatives are written for the first entry

      sum_i = snaps_mat.getEntry(0, i);

      for (int j = 1; j < this->number_snaps; j += 2)
      {
        sum_i += snaps_mat.getEntry(j, i);
      }

      this->snaps_mean[i] = sum_i / (this->number_snaps / 2 + 1);
    }
    else
    {
      for (int j = 0; j < this->number_snaps; j++)
      {
        sum_i += snaps_mat.getEntry(j, i);
      }

      this->snaps_mean[i] = sum_i / this->number_snaps;
    }
  }

  if (with_time_derivatives)
  {
    snaps_mat.addToRow(&snaps_mean, 0, -1.0);

    for (int j = 1; j < this->number_snaps; j += 2)
    {
      snaps_mat.addToRow(&snaps_mean, j, -1.0);
    }
  }
  else
  {
    // subtract the mean from snapshots
    bool transpose = true;
    snaps_mat.add(&snaps_mean, -1.0, transpose);
  }
}

/** ***************************************************************************/
void POD::compute_basis(int l)
{
  bool i_am_root = true;

#ifdef _MPI
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  i_am_root = mpi_rank == 0;
#endif

  this->length = l;
  read_snapshots();

  std::string name = "";
  if (!filename_suffix.empty())
  {
    name += filename_suffix + "-";
  }

  if (db["pod_fluctuations_only"])
  {
    Output::print<1>(name, "POD will be computed from fluctuating part of "
                     "snapshots...");
    POD::decompose_snaps();
  }
  else
  {
    Output::print<1>(name, "POD will be computed from full/raw snapshots...");
  }

  if (db["snaps_time_derivative"])
  {
    double tau = db["snaps_time_derivative_tau"];

    Output::root_info("POD", "Scaling time derivatives by tau = ", tau);

    for (int j = 2; j < this->number_snaps; j += 2)
    {
      snaps_mat.scaleRow(j, tau);
    }

    if (db["snaps_time_derivative_projection"])
    {
      Output::root_info("POD", "Projecting time derivatives back to FE "
        "space...");

      auto solver = Solver<BlockFEMatrix, BlockVector>(projection_solver_db);
      solver.update_matrix(projection_matrix);

      BlockVector dt_u_star(projection_matrix, true);
      BlockVector dt_u(projection_matrix, false);

      for (int j = 2; j < this->number_snaps; j += 2)
      {
        Output::root_info("POD", "Snapshot ", j, "/", this->number_snaps, "...");

        snaps_mat.getRow(j, dt_u_star.get_entries());

        solver.solve(dt_u_star, dt_u);

    #ifdef _MPI
        // TODO: sync
    #endif

        snaps_mat.setRow(j, dt_u.get_entries());
      }
    }
  }

  DenseMatrix autocorr_mat(this->number_snaps, this->number_snaps);
  compute_autocorr_mat(autocorr_mat);

  //----------------------------------------------------------------------------
  // compute eigen-values and eigen-vectors with LAPACK
  //----------------------------------------------------------------------------
  /* parameters for LAPACK dsyev */
  double** autocorr_to_eig;

  if (this->eigs != nullptr)
  {
    delete[] this->eigs;
  }

  /* array for eigenvalues */
  eigs = new double[this->number_snaps];

  autocorr_to_eig = new double*[this->number_snaps];
  autocorr_to_eig[0] = new double[this->number_snaps * number_snaps];

  for (int i = 0; i < this->number_snaps; i++)
  {
    autocorr_to_eig[i] = autocorr_to_eig[0] + i * this->number_snaps;
  }

  for (int i = 0; i < this->number_snaps; i++)
  {
    for (int j = 0; j < this->number_snaps; j++)
    {
      autocorr_to_eig[j][i] = autocorr_mat.getEntry(j, i);
    }
  }

  if (i_am_root)
  {
    // eigenvalue decomposition happens on root node only

    int n = (int)this->number_snaps;
    double* work;
    double wkopt;
    char arg1 = 'V'; /* compute eigenvalues and eigenvectors */
    char arg2 = 'U'; /* for upper triangular matrix */
    int info = 0, lwork;

    /* workspace query: the routine only calculates the optimal size of work */
    lwork = -1;
    dsyev_( &arg1, &arg2, &n, *autocorr_to_eig, &n, eigs, &wkopt, &lwork, &info );

    lwork = (int)wkopt;
    work = (double*)malloc(lwork * sizeof(double));

    /* Solve eigenproblem, the eigenvalues are stored in ascending order */
    dsyev_( &arg1, &arg2, &n, *autocorr_to_eig, &n, eigs, work, &lwork, &info );

    /* Check for convergence */
    if (info > 0)
    {
      ErrThrow("The algorithm failed to compute eigenvalues.");
    }

    //----------------------------------------------------------------------------
    // check eigen-values and rank
    //----------------------------------------------------------------------------
    /* count valid eigenvalues */
    for (int i = this->number_snaps - 1; i >= 0; i--)
    {
      if (this->eigs[i] > this->eigen_threshold)
      {
        this->valid_eigs++;
      }
      else
      {
        break;
      }
    }

    delete[] work;
  }

#ifdef _MPI
  // distribute eigenvalues
  MPI_Bcast(&(this->eigs), this->number_snaps, MPI_DOUBLE,
    0, MPI_COMM_WORLD);

  // distribute eigenvalue count
  MPI_Bcast(&(this->valid_eigs), 1, MPI_INT,
    0, MPI_COMM_WORLD);
#endif

  /* check the rank of the pod basis */
  if (this->rank <= 0)
  {
    this->rank = valid_eigs; // take all
  }

  if (this->rank > this->valid_eigs)
  {
    Output::warn<1>("POD basis", "Rank of POD basis can't be greater than ",
                    this->valid_eigs, "! Parameter 'pod_rank' changed from ",
                    this->rank, " to ", this->valid_eigs, ".");
    this->rank = this->valid_eigs;
  }

  Output::print<1>("  Number of valid ", name, "POD eigenvalues: ",
                   this->valid_eigs);
  Output::print<1>("  Dimension of ", name, "POD basis         : ", this->rank);

  /* Compute the sum of the resting eigenvalues: SUM_i=rank^valid_eigs eigs[i]
   * eigen-values are in ascending order */
  if (rank < valid_eigs)
  {
    double sum_rest_eigs = 0.;
    int nb_rest = std::max(0, valid_eigs - rank);
    for (int i = 0; i < nb_rest; i++)
    {
      int idx = this->number_snaps - 1 - rank - i;
      sum_rest_eigs += eigs[idx];
    }

    Output::print<1>("  Sum of the resting ", name, "eigenvalues : ",
                     sum_rest_eigs);
  }

  //----------------------------------------------------------------------------
  // compute pod_basis
  //----------------------------------------------------------------------------

  /* resize pod_basis */
  DenseMatrix pod_basis(this->rank, this->length_snaps);
  this->pod_basis = pod_basis;

  std::vector<double> vec_eig(this->number_snaps);
  std::vector<double> vec_Seig(this->length_snaps);

  // this is allowed to be MPI-agnostic again, since we're
  // assembling linear combinations of snapshots

  for (int i = 0; i < this->rank; i++)
  {
    /* eigen-values and eigen-vectors are in ascending order */
    int idx = this->number_snaps - 1 - i;

    for (int j = 0; j < this->number_snaps; j++)
    {
      vec_eig[j] = 1. / std::sqrt(this->number_snaps * this->eigs[idx])
                   * autocorr_to_eig[idx][j];
    }

    bool transpose = true;
    vec_Seig = this->snaps_mat.multiply(&vec_eig, 1., transpose);

    for (int j = 0; j < this->length_snaps; j++)
    {
      this->pod_basis.setEntry(i, j, vec_Seig[j]);
    }
  }

  /* reverse order of eigs (descending order) */
  std::reverse(this->eigs, this->eigs + this->number_snaps);

  this->length = this->length_snaps;
  delete[] *autocorr_to_eig;
  delete[] autocorr_to_eig;

  write_pod();
}

/** ***************************************************************************/
void POD::write_pod()
{
  /* writing into basis into file */
  std::string pod_filename = basename + "pod";
  std::string name = "";

#ifdef _MPI
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  pod_filename += "_r" + std::to_string(mpi_rank);
#endif

  if (!filename_suffix.empty())
  {
    name += filename_suffix + "-";
  }

  std::ofstream ofile;
  ofile.open( pod_filename.c_str() , std::ios::out | std::ios::trunc );
  if (!ofile.is_open())
  {
    ErrThrow( "Error: File for POD basis ", pod_filename,
              " could not be created." );
  }
  ofile << setprecision( 12 );

  Output::print<1>("Writing ", name, "POD basis into file: ", pod_filename);

  if (db["pod_fluctuations_only"])
  {
    /* write averages of snapshots into file */
    write_averages();
  }

  for (int i = 0; i< rank; i++)
  {
    for (int j = 0; j < length_snaps; j++)
    {
      ofile << pod_basis.getEntry(i, j) << " ";
    }
    ofile << "\n";
  }

  ofile.close();

  /*write pod eigenvalues into file */
  write_eigenvalues();
}

/** ***************************************************************************/
void POD::read_basis(int rom_rank)
{
  this->rank = rom_rank;

  std::string filename = basename + "pod";

#ifdef _MPI
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  filename += "_r" + std::to_string(mpi_rank);
#endif

  std::vector<std::vector<double>> tmp_basis;
  std::ifstream podfile(filename.c_str());
  std::string line;
  if (!podfile)
  {
    Output::print<1>("Error: POD file ", filename ," could not be opened. "
                     "Check if the parameters in the input file are consistent "
                     "with the parameters used to compute the POD basis");
    exit(4711);
  }

  if (this->db["pod_fluctuations_only"])
  {
    Output::print<1>("POD basis is computed out of fluctuating part of "
                     "snapshots!");
    read_averages();
  }

  std::string name = "";
  if (!filename_suffix.empty())
  {
    name += filename_suffix + "-";
  }

  Output::print<1>("Reading ", name, "POD basis from file ", filename);
  while (getline(podfile, line))
  {
    std::vector<double> data;
    double value;
    std::istringstream iss(line);

    while (iss >> value)
    {
      data.push_back(value);
    }

    tmp_basis.push_back(data);
  }
  podfile.close();

  if (this->rank <= 0)
  {
    this->rank = (int)tmp_basis.size();
  }

  if (this->rank > (int)tmp_basis.size())
  {
    Output::warn<1>("POD basis", "Rank of POD basis can't be greater than ",
                    tmp_basis.size(), "! Parameter 'pod_rank' setting the "
                    "ROM dimension changed to ",
                    tmp_basis.size());

    this->rank = (int)tmp_basis.size();
  }

  this->length = (int)tmp_basis[0].size(); //total dof
  DenseMatrix pod_basis(this->rank, this->length);
  this->pod_basis = pod_basis;

  Output::print<1>("  Length of ", name, "POD-basis vector : ", this->length);
  Output::print<1>("  Dimension of ", name, "POD basis     : ", this->rank);

  //pod basis (for all components)
  for (int i = 0; i < this->length; i++)
  {
    for (int j = 0; j < this->rank; j++)
    {
      this->pod_basis.setEntry(j,i, tmp_basis[j][i]);
    }
  }

  if (rank < (int)tmp_basis.size())
  {
    read_eigenvalues();
  }
}

/** ***************************************************************************/
void POD::read_eigenvalues()
{
  std::string eig_filename = basename + "eigs";

  Output::print<1>("Reading eigenvalues from file ", eig_filename);

  std::vector<double> eigenvalues;
  read_data(eig_filename, eigenvalues, 2);
  this->valid_eigs = eigenvalues.size();

  double sum_rest_eigs = 0.;
  for (int i = rank; i < valid_eigs; i++)
  {
    sum_rest_eigs += eigenvalues[i];
  }
  Output::print<1>("  Sum of the resting eigenvalues : ", sum_rest_eigs);
}

/** ***************************************************************************/
void POD::write_averages()
{
  std::string avr_filename = basename + "mean";

#ifdef _MPI
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  avr_filename += "_r" + std::to_string(mpi_rank);
#endif

  std::ofstream ofile;
  ofile.open(avr_filename.c_str(), std::ios::out | std::ios::trunc);

  if (!ofile.is_open())
  {
    ErrThrow("Error: File for averages of snapshots ", avr_filename,
             " could not be created.");
  }

  ofile << setprecision(12);

  if ((int)this->snaps_mean.size() != this->length)
  {
    ErrThrow( "Error: Vector for averages of snapshots has the wrong length!" );
  }

  // write elements
  for (int i = 0; i < (int)this->snaps_mean.size(); ++i)
  {
      ofile << this->snaps_mean[i] << " ";
  }
  ofile.close();
}

/** ***************************************************************************/
void POD::write_eigenvalues()
{
  std::string eig_filename = basename + "eigs";

  int valid_eigen = 0;
  int no_eigen = this->number_snaps;
  double sum_eigen = 0.0, cumulative=0.0;

  if (this->eigs == nullptr)
  {
    ErrThrow("Error: Eigenvalues not found. POD not created?");
  }

  for (int i = 0; i < no_eigen; i++)
  {
    if (this->eigs[i] >= this->eigen_threshold)
    {
      valid_eigen++;
    }
  }

  std::ofstream ofile;
  ofile.open(eig_filename.c_str(), std::ios::out | std::ios::trunc);
  ofile << setprecision(12);

  for (int i = 0; i < no_eigen; i++)
  {
    sum_eigen += this->eigs[i];
  }

  for (int i = 0; i < std::min(no_eigen,valid_eigen); i++)
  {
    cumulative += this->eigs[i] / sum_eigen;
    ofile << i + 1 << " \t" << this->eigs[i] << " \t" << cumulative << endl;
  }

  ofile.close();
}

/** ***************************************************************************/
void POD::read_averages()
{
  std::string avg_filename = basename + "mean";

#ifdef _MPI
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  avg_filename += "_r" + std::to_string(mpi_rank);
#endif

  Output::print<1>("Reading averages of snapshots from file ", avg_filename);

  std::vector<double> data;

  read_data(avg_filename , data);
  this->snaps_mean.resize(data.size());

  for (int i = 0; i < (int)this->snaps_mean.size(); ++i)
  {
    this->snaps_mean[i] = data[i];
  }
}

/** ***************************************************************************/
void POD::read_data(std::string filename, std::vector<double>& data, int format)
{
  std::string line;
  std::ifstream file(filename.c_str());
  if (!file)
  {
    ErrThrow("Error: File ", filename, " could not be opened.");
  }

  while (getline(file, line))
  {
    double value;
    int column_i = 1;

    std::istringstream iss(line);
    while (iss >> value)
    {
      if (format == 0 || format == column_i)
      {
        data.push_back(value);
      }

      ++column_i;

      if (format != 0 && format < column_i)
      {
        break;
      }
    }
  }

  file.close();
}

/** ***************************************************************************/
void POD::compute_autocorr_mat(DenseMatrix& corr_mat)
{
  Output::print<1>("Computing autocorrelation matrix...");

  if (this->length_snaps == 0)
  {
    ErrThrow("Error: Snapshots are not available! Before computing POD "
             "basis read_snapshots() must be called!");
  }

  if (db["pod_inner_product"].is("euclidean"))
  {
  #ifdef _MPI
    Output::root_warn("POD", "POD with the Euclidean inner product will not "
      "give correct results in the MPI case since the POD class does not know "
      "what an FE space is. The inner product computed in this case is"
      "\n\n"
      "(phi_i, phi_j) = [number of ranks aware of both DOF #i and DOF #j]"
      "\n\n"
      "and so is dependent on the domain decomposition.");
  #endif

    bool transp_A = false;
    bool transp_B = true;

    corr_mat = *(snaps_mat.multiply( &this->snaps_mat, transp_A, transp_B ));
  }
  else
  {
    if (gramian_ptr == nullptr)
    {
      ErrThrow("POD with ", db["pod_inner_product"], " must be supplied with a "
        "Gramian matrix before calling compute_autocorr_mat.");
    }

    if ((int)this->gramian_ptr->get_n_rows() != this->length_snaps)
    {
      ErrThrow("ERROR: Dimension of inner product matrix and length of "
               "snapshots must coincide!\nCurrently: dim of matrix : ",
               this->gramian_ptr->get_n_rows(), " length of snapshots : ",
               this->length_snaps);
    }

    bool transpose = true;
    std::shared_ptr<DenseMatrix> tmp_mat = gramian_ptr->multiply(&snaps_mat,
                                                                 transpose);

    corr_mat = *(this->snaps_mat.multiply(tmp_mat.get()));
  } // end else

#ifdef _MPI
  MPI_Allreduce(MPI_IN_PLACE, corr_mat.get_entries(),
    this->number_snaps * this->number_snaps, MPI_DOUBLE,
    MPI_SUM, MPI_COMM_WORLD);
#endif

  // scale autocorrelation matrix with the number of snapshots according to
  //
  corr_mat.multiply(1. / this->number_snaps);
  Output::print<1>("... computation of autocorrelation matrix completed.");
}
