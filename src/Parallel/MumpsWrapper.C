/**
 * Implementation of class MumpsWrapper, declared in MumpsWrapper.h
 *
 * @author Clemens Bartsch
 * @date 2016/03/14
 */
#ifdef _MPI

#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <MumpsWrapper.h>
#include <ParFECommunicator3D.h>
#include "MooNMD_Io.h"

#include <mpi.h>
#include <memory>

#ifdef PARMOON_WITH_MUMPS

// two of the used mumps job codes
#define JOB_INIT -1
#define JOB_END -2
#define JOB_ANALYZE 1
#define JOB_FACTORIZE 2
#define JOB_SOLVE 3

// macros doing the fortran shift for mumps control and info params
#define ICNTL(I) icntl[(I)-1]
#define INFO(I) info[(I)-1]
#define INFOG(I) infog[(I)-1]

MumpsWrapper::MumpsWrapper(const BlockFEMatrix& bmatrix,
                           std::vector<double> pres0)
: analyzed_and_factorized(false), distributed_solution(false),
  distributed_solution_setup_done(false),
  distributed_rhs_setup_done(false),
  original_matrix(bmatrix)
{
  // check the input
  check_input_matrix(bmatrix);

  // copy the BlockFEMatrices communicators and store them
  comms_ = bmatrix.get_communicators();

  // initialize the mumps entity
  // set those "hard" parameters which must be set before initializing
  id_.par = 1; // root will take part in computation
  id_.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);
  id_.sym = 0; // non-symmetric matrix
  id_.job = JOB_INIT;

  dmumps_c(&id_);

  // the "softer" parameters must be set after initializing
  set_mumps_parameters();

  // transform and store the matrix distributed in coordinate format
  store_in_distributed_coordinate_form(bmatrix, pres0);
}

MumpsWrapper::MumpsWrapper( const BlockMatrix& bmatrix,
                            std::vector<const TParFECommunicator3D*> comms,
                            std::vector<double> pres0,
                            std::vector<std::vector<int>> loc_to_seq)
: analyzed_and_factorized(false), distributed_solution(false),
  distributed_solution_setup_done(false),
  distributed_rhs_setup_done(false),
  original_matrix(bmatrix)
{
  // check the input
  check_input_matrix(bmatrix);

  // copy the BlockFEMatrices communicators and store them
  comms_ = comms; 

  // initialize the mumps entity
  // set those "hard" parameters which must be set before initializing
  id_.par = 1; // root will take part in computation
  id_.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);
  id_.sym = 0; // non-symmetric matrix
  id_.job = JOB_INIT;

  dmumps_c(&id_);

  // the "softer" parameters must be set after initializing
  set_mumps_parameters();

  // transform and store the matrix distributed in coordinate format
  store_in_distributed_coordinate_form(bmatrix, pres0, loc_to_seq);
}

void MumpsWrapper::enable_distributed_solution(int max_iterations, double residual_tolerance)
{
  if (analyzed_and_factorized && !distributed_solution)
  {
    ErrThrow("MumpsWrapper cannot be set distributed after the first call to solve()!");
  }

  distributed_solution = true;
  distributed_max_iterations = max_iterations;
  distributed_residual_tolerance = residual_tolerance;

  set_mumps_parameters();
}

void MumpsWrapper::setup_distributed_solution()
{
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  size_t n_blocks = comms_.size();

  distr_send_count.clear();
  distr_send_count.resize(mpi_size, 0);
  distr_send_displ.resize(mpi_size + 1, 0);

  distr_recv_count.clear();
  distr_recv_count.resize(mpi_size, 0);
  distr_recv_displ.resize(mpi_size + 1, 0);

  // count master dofs

  int n_loc_masters = 0;

  for (size_t bl = 0; bl < n_blocks; ++bl)
  {
    int n_loc_masters_block = comms_.at(bl)->GetN_Master();
    n_loc_masters += n_loc_masters_block;
  }

  // communicate dof counts etc. rank to rank

  int max_send_total;

  MPI_Allreduce(&id_.lsol_loc, &max_send_total, 1,
    MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  std::vector<int> n_rank_masters(mpi_size, 0);

  MPI_Allgather(&n_loc_masters, 1, MPI_INT,
    n_rank_masters.data(), 1, MPI_INT,
    MPI_COMM_WORLD);

  // figure out index distribution

  std::vector<std::vector<int>> send_idx(mpi_size);
  std::vector<int> send_idx_flat(max_send_total * mpi_size);
  std::vector<int> recv_idx_flat(max_send_total * mpi_size);

  for (int i = 0; i < id_.lsol_loc; i++)
  {
    // we're getting our information from fortran land,
    // so the indices in isol_loc are 1-based
    int glob_master_index = id_.isol_loc[i] - 1;

    int dof_rank = -1;
    int loc_master_index = glob_master_index;

    for (int r = 0; r < mpi_size; r++)
    {
      if (loc_master_index >= 0 && loc_master_index < n_rank_masters[r])
      {
        dof_rank = r;
        break;
      }

      loc_master_index -= n_rank_masters[r];
    }

    if (dof_rank < 0)
    {
      Output::warn("MumpsWrapper", "Global DOF ", glob_master_index,
        " does not seem to belong to any rank!");
      continue;
    }

    int j = (int)send_idx[dof_rank].size();

    send_idx[dof_rank].push_back(i);
    send_idx_flat[dof_rank * max_send_total + j] = loc_master_index;
  }

  for (int r = 0; r < mpi_size; r++)
  {
    distr_send_count[r] = (int)send_idx[r].size();
  }

  // communicate indices sent to each rank

  MPI_Alltoall(distr_send_count.data(), 1, MPI_INT,
    distr_recv_count.data(), 1, MPI_INT,
    MPI_COMM_WORLD);

  MPI_Alltoall(send_idx_flat.data(), max_send_total, MPI_INT,
    recv_idx_flat.data(), max_send_total, MPI_INT,
    MPI_COMM_WORLD);

  // count block sizes and figure out master dofs

  std::vector<int> block_master_count(n_blocks);
  std::vector<std::vector<int>> block_masters(n_blocks);

  int block_offset = 0;

  for (size_t bl = 0; bl < n_blocks; bl++)
  {
    int master_count = comms_[bl]->GetN_Master();
    int dof_count = comms_[bl]->GetNDof();
    const int* masters = comms_[bl]->GetMaster();

    block_masters[bl].resize(master_count);
    block_master_count[bl] = master_count;

    int master_index = 0;

    for (int i = 0; i < dof_count; i++)
    {
      if (masters[i] == mpi_rank)
      {
        block_masters[bl][master_index++] = block_offset + i;
      }
    }

    block_offset += dof_count;
  }

  // figure out received indices

  std::vector<std::vector<int>> recv_idx(mpi_size);

  int recv_total = 0;

  for (int r = 0; r < mpi_size; r++)
  {
    int n = distr_recv_count[r];

    for (int i = 0; i < n; i++)
    {
      int master_index = recv_idx_flat[r * max_send_total + i];
      int block_master_index = master_index;
      int loc_dof_index = -1;

      for (size_t bl = 0; bl < n_blocks; bl++)
      {
        int n = block_master_count[bl];

        if (block_master_index >= 0 && block_master_index < n)
        {
          loc_dof_index = block_masters[bl][block_master_index];
          break;
        }

        block_master_index -= n;
      }

      if (loc_dof_index < 0)
      {
        Output::warn("MumpsWrapper", "DOF ", master_index, " on rank ", mpi_rank,
          " does not seem to belong to any block!");
        continue;
      }

      recv_idx[r].push_back(loc_dof_index);
    }

    recv_total += (int)recv_idx[r].size();
  }

  // allocate index arrays

  distr_send_index.resize(id_.lsol_loc);
  distr_recv_index.resize(recv_total);

  // allocate data buffers

  distr_send_buf.resize(id_.lsol_loc);
  distr_recv_buf.resize(recv_total);

  // build displ and index arrays

  int sdispl = 0;
  int rdispl = 0;

  for (int r = 0; r < mpi_size; r++)
  {
    distr_send_displ[r] = sdispl;
    distr_recv_displ[r] = rdispl;

    int n_send = distr_send_count[r];
    int n_recv = distr_recv_count[r];

    for (int i = 0; i < n_send; i++)
    {
      distr_send_index[sdispl + i] = send_idx[r][i];
    }

    for (int i = 0; i < n_recv; i++)
    {
      distr_recv_index[rdispl + i] = recv_idx[r][i];
    }

    sdispl += n_send;
    rdispl += n_recv;
  }

  distr_send_displ[mpi_size] = sdispl;
  distr_recv_displ[mpi_size] = rdispl;
}

void MumpsWrapper::setup_distributed_rhs()
{
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int n_loc_masters = 0;

  for (size_t bl = 0; bl < comms_.size(); bl++)
  {
    n_loc_masters += comms_[bl]->GetN_Master();
  }

  MPI_Allreduce(&n_loc_masters, &distr_lrhs_loc,
    1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  // each process contributes exactly its master dofs

  distr_rhs_loc.resize(distr_lrhs_loc);
  distr_irhs_loc.resize(n_loc_masters);
  distr_rhs_index.resize(n_loc_masters);

  // figure out master offset by summing up earlier ranks' master counts

  std::vector<int> masters_per_rank(mpi_size);

  MPI_Allgather(&n_loc_masters, 1, MPI_INT,
    masters_per_rank.data(), 1, MPI_INT,
    MPI_COMM_WORLD);

  // we're sending these indices to fortran land, so add 1
  int master_offset = 1;

  for (int r = 0; r < mpi_rank; r++)
  {
    master_offset += masters_per_rank[r];
  }

  // set up rhs indices

  int master_index = 0;

  size_t n_blocks = comms_.size();
  int dof_offset = 0;

  for (size_t bl = 0; bl < n_blocks; bl++)
  {
    int dof_count = comms_[bl]->GetNDof();
    const int* masters = comms_[bl]->GetMaster();

    for (int i = 0; i < dof_count; i++)
    {
      if (masters[i] == mpi_rank)
      {
        distr_rhs_index[master_index] = dof_offset + i;
        distr_irhs_loc[master_index] = master_offset + master_index;

        ++master_index;
      }
    }

    dof_offset += dof_count;
  }
}

void MumpsWrapper::prepare_rhs_distributed(const BlockVector& rhs, BlockVector& solution)
{
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int n_masters_local_comms;
  int n_dofs_global_comms;
  check_input_solve(rhs, solution, n_masters_local_comms, n_dofs_global_comms);

  size_t n = distr_rhs_index.size();
  for (size_t i = 0; i < n; i++)
  {
    distr_rhs_loc[i] = rhs[distr_rhs_index[i]];
  }
}

void MumpsWrapper::redistribute_solution_distributed(BlockVector& solution)
{
  int nsend = distr_send_index.size();
  for (int i = 0; i < nsend; i++)
  {
    distr_send_buf[i] = id_.sol_loc[distr_send_index[i]];
  }

  MPI_Alltoallv(
    distr_send_buf.data(), distr_send_count.data(), distr_send_displ.data(), MPI_DOUBLE,
    distr_recv_buf.data(), distr_recv_count.data(), distr_recv_displ.data(), MPI_DOUBLE,
    MPI_COMM_WORLD);

  int nrecv = distr_recv_index.size();
  for (int i = 0; i < nrecv; i++)
  {
    solution[distr_recv_index[i]] = distr_recv_buf[i];
  }

  for (size_t bl = 0; bl < comms_.size(); bl++)
  {
    comms_.at(bl)->queue_consistency_update(solution.block(bl), 3);
  }

  TParFECommunicator3D::flush_consistency_updates();
}

void MumpsWrapper::prepare_rhs_centralized(const BlockVector& rhs, BlockVector& solution)
{
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int root_rank = 0;
  bool i_am_root = (mpi_rank == root_rank);

  // 0) do input check and gather two important numbers
  int n_masters_local_comms;
  int n_dofs_global_comms;
  check_input_solve(rhs, solution, n_masters_local_comms, n_dofs_global_comms);

  // 1) set up the global right hand side in root
  master_values.clear(); // note that clear() does not deallocate, so no heap churn
  master_values.reserve(n_masters_local_comms); // space for loc rhs and solution

  if (i_am_root) // put up an array for the global right hand side
  {
    rhs_global.resize(n_dofs_global_comms, 0.0);
  }

  // three different shifts in this loop, all due to arrangement in blocks
  // - "loc_dof_shift", accounts for all local dofs (master, slave, halo)
  // - "loc_master_shift", which accounts for master dofs only
  // - "glob_dof_shift", which accounts for the global number of dofs per block
  //   (determined by adding up master dofs)
  size_t loc_master_shift = 0;
  size_t loc_dof_shift = 0;
  size_t glob_dof_shift = 0;

  for (size_t index = 0; index < comms_.size(); ++index) // loop over blocks
  {
    // fill the local right hand side with master dof rows only
    const int* masters = comms_.at(index)->GetMaster(); // TODO this should be a vector (in ParFECommunicator)!

    size_t n_loc_dofs_block = comms_.at(index)->GetNDof();

    for (size_t i = 0; i < n_loc_dofs_block; ++i)
    {
      // push rhs values for master dofs on this rank and block to master_values

      if (masters[i] == mpi_rank)
      {
        master_values.push_back(rhs.at(loc_dof_shift + i));
      }
    }

    // ...and gather these local right hand sides globally
    int n_loc_masters_block = comms_.at(index)->GetN_Master();
    double* global_rhs_dummy = nullptr;
    if (i_am_root)
    {
      global_rhs_dummy = &rhs_global.at(glob_dof_shift);
    }

    // BUGFIX: make sure that master_values.at(loc_master_shift) does not throw,
    // even if the number of masters on this rank in this block was 0
    if (n_loc_masters_block == 0)
    {
      master_values.push_back(0);
    }

    gather_vector(
        global_rhs_dummy,                                          // receive
        &master_values.at(loc_master_shift), n_loc_masters_block,  // send
        root_rank);                                                // control

    if (n_loc_masters_block == 0)
    {
      master_values.pop_back();  // revert action due to BUGFIX above
    }

    // one block treated - count up the shifts
    loc_master_shift += n_loc_masters_block;
    loc_dof_shift += n_loc_dofs_block;

    int current_n_dofs_global;
    MPI_Allreduce(&n_loc_masters_block, &current_n_dofs_global, 1,
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    glob_dof_shift += current_n_dofs_global;
  }
}

void MumpsWrapper::redistribute_solution_centralized(BlockVector& solution)
{
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  int root_rank = 0;
  bool i_am_root = (mpi_rank == root_rank);

  size_t loc_dof_shift = 0;
  size_t loc_master_shift = 0;
  size_t glob_dof_shift = 0;

  for (size_t index = 0; index < comms_.size(); ++index) // loop over blocks
  {
    // receive all master values from root and store in master_values
    int n_loc_masters_block = comms_.at(index)->GetN_Master();
    int n_glob_dofs_block;

    MPI_Allreduce(&n_loc_masters_block, &n_glob_dofs_block,
                  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    double* global_rhs_dummy = nullptr;
    if (i_am_root)
    {
      global_rhs_dummy = &rhs_global.at(glob_dof_shift);
    }

    // BUGFIX: make sure that master_values.at(loc_master_shift) does not throw,
    // even if the number of masters on this rank in this block was 0
    if (n_loc_masters_block == 0)
    {
      master_values.push_back(0);
    }

    scatter_vector(
        global_rhs_dummy,                                         // send
        &master_values.at(loc_master_shift), n_loc_masters_block, // receive
        root_rank);                                               // control

    // BUGFIX: make sure that master_values.at(loc_master_shift) does not throw,
    // even if the number of masters on this rank in this block was 0
    if (n_loc_masters_block == 0)
    {
      master_values.pop_back();
    }

    // write master values into local solution vector
    const int* masters = comms_.at(index)->GetMaster(); // TODO this should be a vector (in ParFECommunicator)!
    int n_local_dofs_block = comms_.at(index)->GetNDof();

    int i_master = 0;
    for (int i_dof = 0; i_dof < n_local_dofs_block; ++i_dof)
    {
      // write solution values for master dofs on this rank into solution vector

      if (masters[i_dof] == mpi_rank)
      {
        solution.at(loc_dof_shift + i_dof)
            = master_values.at(loc_master_shift + i_master);
        ++i_master;
      }
    }

    // Update all non-master values in solution vector. Big fire!
    comms_.at(index)->consistency_update(solution.block(index), 3);

    // count up the block shifts
    int n_loc_dofs_block = comms_.at(index)->GetNDof();
    loc_dof_shift += n_loc_dofs_block;
    loc_master_shift += n_loc_masters_block;
    glob_dof_shift += n_glob_dofs_block;
  }
}

void MumpsWrapper::solve(const BlockVector& rhs, BlockVector& solution, int iteration)
{
  double rhs_max = rhs.norm_infty(comms_);

  if (rhs_max == 0.0)
  {
    // RHS is zero - don't bother doing anything else, just return zero

    Output::root_info("MumpsWrapper", "RHS is zero. Bless.");

    solution.reset();

    return;
  }

  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  int root_rank = 0;
  bool i_am_root = (mpi_rank == root_rank);

  if (distributed_solution)
  {
    if (!distributed_rhs_setup_done)
    {
      setup_distributed_rhs();

      distributed_rhs_setup_done = true;
    }

    prepare_rhs_distributed(rhs, solution);
  }
  else
  {
    prepare_rhs_centralized(rhs, solution);
  }

  // 2) let mumps do its jobs
  if (!analyzed_and_factorized)
  {
    id_.nz_loc = matrix_.nz_loc;
    id_.irn_loc = matrix_.irn_loc.data();
    id_.jcn_loc = matrix_.jcn_loc.data();
    id_.a_loc = matrix_.a_loc.data();
  }

  id_.nrhs = 1;
  id_.n = matrix_.n;

  if (distributed_solution)
  {
    id_.nloc_rhs = (int)distr_irhs_loc.size();
    id_.lrhs_loc = distr_lrhs_loc;
    id_.irhs_loc = distr_irhs_loc.data();
    id_.rhs_loc = distr_rhs_loc.data();
  }
  else
  {
    if (i_am_root)
    {
      id_.rhs = rhs_global.data();
      id_.lrhs = matrix_.n;
    }
  }

  if (!analyzed_and_factorized)
  {
    kick_off_job("analyze");

    long real_estimate = id_.INFOG(3);
    long int_estimate = id_.INFOG(4);

    if (real_estimate < 0)
    {
      real_estimate *= -1000 * 1000;
    }

    if (int_estimate < 0)
    {
      int_estimate *= -1000 * 1000;
    }

    int mbytes_estimate_max = id_.INFOG(16);
    int mbytes_estimate_sum = id_.INFOG(17);

    Output::root_info("MumpsWrapper", "Analysis space estimates: "
      "\nReal workspace: ", real_estimate, " entries"
      "\nInteger workspace: ", int_estimate, " entries"
      "\nTotal memory: ", mbytes_estimate_sum, " MB total, ",
      mbytes_estimate_max, " MB max.");

    kick_off_job("factorize");

    long real_space = id_.INFOG(9);
    long int_space = id_.INFOG(10);

    if (real_space < 0)
    {
      real_space *= -1000 * 1000;
    }

    if (int_space < 0)
    {
      int_space *= -1000 * 1000;
    }

    int mbytes_space_max = id_.INFOG(18);
    int mbytes_space_sum = id_.INFOG(19);

    Output::root_info("MumpsWrapper", "Factorize space requirements: "
      "\nReal workspace: ", real_space, " entries"
      "\nInteger workspace: ", int_space, " entries"
      "\nTotal memory: ", mbytes_space_sum, " MB total, ",
      mbytes_space_max, " MB max.");

    if (distributed_solution)
    {
      // set up distributed solution arrays

      int lsol_loc = id_.INFO(23);

      distr_sol_loc.resize(lsol_loc);
      distr_isol_loc.resize(lsol_loc);

      id_.lsol_loc = lsol_loc;
      id_.sol_loc = distr_sol_loc.data();
      id_.isol_loc = distr_isol_loc.data();

      distributed_solution_setup_done = false;
    }

    analyzed_and_factorized = true;
  }

  kick_off_job("solve");

  // 3) distribute solution among processes

  if (distributed_solution)
  {
    // write the sol_loc over into solution

    if (!distributed_solution_setup_done)
    {
      // set up the buffers for sending/receiving distributed solution data

      setup_distributed_solution();

      distributed_solution_setup_done = true;
    }

    redistribute_solution_distributed(solution);

    if (iteration == 0)
    {
      BlockVector iteration_residual(rhs);

      // solution is at consistency level 3 after redistribution
      original_matrix.apply(solution, iteration_residual);
      iteration_residual -= rhs;

      double residual = iteration_residual.norm_infty(comms_) / rhs_max;

      Output::root_info<2>("MumpsWrapper", "Iteration ", iteration, " residual: ", residual,
        ", abs: ", residual * rhs_max);

      if (residual > distributed_residual_tolerance && distributed_max_iterations > 1)
      {
        BlockVector iteration_solution(rhs);

        for (int i = 1; i < distributed_max_iterations; i++)
        {
          // r = Ax - y
          // => A^-1 y = x - A^-1 r

          solve(iteration_residual, iteration_solution, i);
          solution -= iteration_solution;

          // iteration_solution is at consistency level 3 after solve(),
          // so solution does not lose consistency
          original_matrix.apply(solution, iteration_residual);
          iteration_residual -= rhs;

          residual = iteration_residual.norm_infty(comms_) / rhs_max;

          Output::root_info<2>("MumpsWrapper", "Iteration ", i, " residual: ", residual,
            ", abs: ", residual * rhs_max);

          if (residual <= distributed_residual_tolerance)
          {
            break;
          }
        }
      }
    }
  }
  else
  {
    redistribute_solution_centralized(solution);
  }
}

void MumpsWrapper::write_matrix_distributed(const std::string& filename) const
{
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  std::ofstream matrixfile;
  std::string file = filename + std::to_string(mpi_rank);
  matrixfile.open(file.c_str());

  // write the header line - coordinate format, real values, no symmetry used
  matrixfile << "%%MatrixMarket matrix coordinate real general \n";

  // write general matrix information
  matrixfile << matrix_.n << "\t" << matrix_.n << "\t" << matrix_.nz_loc << "\n";

  // loop and write info to file
  for (size_t index = 0; index < matrix_.nz_loc; ++index)
  {
    matrixfile << matrix_.irn_loc.at(index) << "\t" << matrix_.jcn_loc.at(index) << "\t"
        << matrix_.a_loc.at(index) << "\n";
  }
}

// Special member function.
MumpsWrapper::~MumpsWrapper()
{
  Output::root_info<2>("MumpsWrapper", "Cleaning up MUMPS instance.");

  id_.job = JOB_END;
  dmumps_c(&id_);
}

// Private functions.
void MumpsWrapper::check_input_solve(
    const BlockVector& rhs, const BlockVector& solution,
    int& n_masters_local_comms, int& n_dofs_global_comms)
{
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  // 0) input checks
  n_masters_local_comms = 0;
  n_dofs_global_comms = 0;
  for (size_t i = 0; i < comms_.size(); ++i)
  {
    if (comms_.at(i)->get_n_dim() != 1 )
    {
      ErrThrow("MumpsWrapper can only be used with "
        "ParFECommunicators of dimension 1");
    }

    size_t local_n_dofs = comms_.at(i)->GetNDof();
    size_t block_length_rhs = rhs.length(i);
    size_t block_length_sol = solution.length(i);

    if (local_n_dofs != block_length_rhs)
    {
      ErrThrow("Length of rhs block ", i, " does not"
        " fit n of local dofs in given communicator.");
    }

    if (local_n_dofs != block_length_sol)
    {
      ErrThrow("Length of sol block ", i, " does not fit n "
        "of local dofs in given communicator.");
    }

    // add up all masters
    n_masters_local_comms += comms_.at(i)->GetN_Master();
  }

  MPI_Allreduce(&n_masters_local_comms, &n_dofs_global_comms,
    1, MPI_INT,MPI_SUM, MPI_COMM_WORLD);

  if (n_dofs_global_comms != (int) matrix_.n)
  {
    ErrThrow("Total number of masters of given communicators "
      "not equal global stored matrix order.");
  }
}

void MumpsWrapper::kick_off_job(const std::string& job)
{
  int job_id = 0;

  // determine job
  if (job.compare("analyze") == 0)
  {
    job_id = JOB_ANALYZE;
  }
  else if (job.compare("factorize") == 0)
  {
    job_id = JOB_FACTORIZE;
  }
  else if (job.compare("solve") == 0)
  {
    job_id = JOB_SOLVE;
  }
  else
  {
    ErrThrow("The string '",  job ,"' does not describe a MUMPS job!");
  }

  // run the chosen job
  id_.job = job_id;
  dmumps_c(&id_);

  int err = id_.INFOG(1);

  // handle workspace size errors by allowing more space
  while ((err == -8
    || err == -9
    || err == -11
    || err == -12
    || err == -14
    || err == -15
    || err == -17
    || err == -20)
    && id_.ICNTL(14) < 2000)
  {
    id_.ICNTL(14) += 100;

    Output::root_warn("MumpsWrapper", "Insufficient workarray size. "
      "Increasing ICNTL(14) to ", id_.ICNTL(14), "%.");

    dmumps_c(&id_);
    err = id_.INFOG(1);
  }

  // inform about eventual errors
  if (err < 0)
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    bool i_am_root = (rank == 0);
    if (i_am_root)
    {
      long val = id_.INFOG(2);
      if (val < 0)
      {
        val *= -1000 * 1000;
      }

      Output::info("MumpsWrapper", "MUMPS job ", job ," failed INFOG(1) = ", err);
      Output::info("MumpsWrapper", "MUMPS job ", job ," failed INFOG(2) = ", val);

      switch (err)
      {
        // error codes as of MUMPS 5.4.1

        case -2:
        {
          Output::info("MumpsWrapper", "NNZ out of range:\nNNZ = ", val);
          break;
        }

        case -3:
        {
          Output::info("MumpsWrapper", "JOB out of range or unavailable:"
            "\nJOB = ", val);
          break;
        }

        case -4:
        {
          Output::info("MumpsWrapper", "Error in PERM_IN at position i: (1-based)"
            "\ni = ", val);
          break;
        }

        case -5:
        {
          Output::info("MumpsWrapper", "Error allocating real workspace for analysis:"
            "\ntried to allocate n entries, n = ", val);
          break;
        }

        case -6:
        {
          Output::info("MumpsWrapper", "Matrix is structurally singular:"
            "\nstructural rank = ", val);
          break;
        }

        case -7:
        {
          Output::info("MumpsWrapper", "Error allocating integer workspace for analysis:"
            "\ntried to allocate n entries, n = ", val);
          break;
        }

        case -8:
        {
          Output::info("MumpsWrapper", "Integer workarray too small.");
          break;
        }

        case -9:
        {
          Output::info("MumpsWrapper", "Real workarray too small:"
            "\nMissing n entries, n = ", val);
          break;
        }

        case -10:
        {
          Output::info("MumpsWrapper", "Matrix is numerically singular:"
            "\nEliminated n pivots, n = ", val);
          break;
        }

        case -11:
        {
          Output::info("MumpsWrapper", "Real workarray too small to hold solution:"
            "\nMissing n entries, n = ", val);
          break;
        }

        case -12:
        {
          Output::info("MumpsWrapper", "Real workarray too small for iterative refinement:"
            "\nMissing n entries, n = ", val);
          break;
        }

        case -13:
        {
          Output::info("MumpsWrapper", "Error allocating real workspace:"
            "\ntried to allocate n entries, n = ", val);
          break;
        }

        case -14:
        {
          Output::info("MumpsWrapper", "Integer workarray too small to hold solution:"
            "\nMissing n entries, n = ", val);
          break;
        }

        case -15:
        {
          Output::info("MumpsWrapper", "Integer workarray too small for iterative refinement:"
            "\nMissing n entries, n = ", val);
          break;
        }

        case -16:
        {
          Output::info("MumpsWrapper", "N out of range:\nN = ", val);
          break;
        }

        case -17:
        {
          Output::info("MumpsWrapper", "Internal send buffer too small.");
          break;
        }

        case -19:
        {
          Output::info("MumpsWrapper", "Maximum working memory size too small:"
            "\nMissing n entries, n = ", val);
          break;
        }

        case -20:
        {
          Output::info("MumpsWrapper", "Internal receive buffer too small:"
            "\nShould have n bytes, n = ", val);
          break;
        }

        default:
        {
          Output::info("MumpsWrapper", "Some other issue has occurred. Please check the "
            "MUMPS user guide for details.");
          break;
        }
      }
    }

    MPI_Finalize();
    exit(-1);
  }
}

void MumpsWrapper::set_mumps_parameters()
{
  // input format parameters
  id_.ICNTL(5)  = 0; // matrices in "assembled" format
  id_.ICNTL(18) = 3; // structure AND entries are distributed among processes

  if (distributed_solution)
  {
    id_.ICNTL(20) = 11; // distributed right hand side
    id_.ICNTL(21) = 1; // distributed solution

    // iterative refinement is only available when the solution is centralized
    id_.ICNTL(10) = 0;
  }
  else
  {
    id_.ICNTL(20) = 0; // dense right hand side (stored in nrhs and lrhs)
    id_.ICNTL(21) = 0; // centralized solution

    // This parameter forces MUMPS to repeat the solution 2 times
    // using a simple iterative scheme. This (usually) makes the solution
    // more exact by several orders of magnitude.
    id_.ICNTL(10) = -2;
  }

  // parameters collected from old MumpsSolver.C
  id_.ICNTL(1) = 6;    // standard outstream for errors
  id_.ICNTL(2) = 6;    // standard warnings output
  id_.ICNTL(3) = 0;    // global info output - surpressed ("max trans not allowed...")
  id_.ICNTL(4) = 1;    // verbosity level
  id_.ICNTL(14) = 100;  // estimated working space increase (%)

  // the following block is for choice of ordering tools in analysis phase
  // FIXME parellel ordering with parmetis is segfaulting!
  id_.ICNTL(28) = 1; // request seq(1)/par(2) ordering in analysis phase
  // request mumps to independently choose a sequential ordering tool if id_.ICNTL(28)=1
  // in the current setup METIS and PORD are available plus some built-in orderings
  id_.ICNTL(7) = 7;
  // if id_.ICNTL(28)=2 ordering is done in parallel, if a tool is available
  // FIXME in the current setup of the libraries, parmetis is segfaulting when
  // called - that's why id_.ICNTL(28) is set to 1 (sequential ordering) so far
  id_.ICNTL(29) = 2; // request parmetis to do the ordering if id.ICNTL(28)=2
}

void MumpsWrapper::store_in_distributed_coordinate_form(
    const BlockMatrix& bmatrix,
    std::vector<double> pres0,
    std::vector<std::vector<int>> loc_to_seq
)
{
  // check what kind of block matrix we deal with
  bool is_block_fe_matrix = true;

  try
  {
    dynamic_cast<const BlockFEMatrix&>(bmatrix);
  }
  catch (const std::bad_cast& e)
  {
    // cast did not work
    is_block_fe_matrix = false;
  }

  // no input checks - assume that this method is only called from
  // the constructor, which already performed checks
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  // 0) preparation
  // 0.0) determine global number of dofs by adding up the masters
  size_t n_dofs_global = 0; // sum over all local masters is number of global dofs
  {
    size_t n_masters_local = 0;
    for (auto comm: comms_)
    {
      n_masters_local += comm->GetN_Master();
    }

    MPI_Allreduce(&n_masters_local, &n_dofs_global,
                  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }

  // 0.1)
  // determine an array which holds index-shifts
  size_t nComms = comms_.size();
  std::vector<size_t> shifts(nComms, 0);

  // for each space, find out how many dofs there are in total
  // (over all processors) - adding up those will give the shift
  for (size_t index = 0; index < nComms - 1; ++index)
  {
    int n_masters_local = comms_.at(index)->GetN_Master();

    // sum up all nLocalMasters over all processes and store
    MPI_Allreduce(&n_masters_local, &shifts.at(index + 1),
                  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // add up the shift with the previous one
    shifts.at(index + 1) += shifts.at(index);
  }

  // 0.2)
  // peel the master arrays from the Communicators -make copies!
  std::vector<std::vector<int>> dof_masters(nComms);
  for (size_t index = 0; index < nComms; ++index)
  {
    const int* master_array = comms_.at(index)->GetMaster();

    size_t master_array_length = comms_.at(index)->GetNDof();

    dof_masters.at(index).resize(master_array_length);

    std::copy(master_array,
              master_array + master_array_length, // pointer arithmetic!
              dof_masters.at(index).begin());
  }

  // 0.3)
  // peel the local2global vectors from the communicators - make copies!
  std::vector<std::vector<int>> local2globals(nComms);
  for (size_t index = 0; index < nComms; ++index)
  {
    const int* l2g = comms_.at(index)->Get_Local2Global();

    size_t array_length = comms_.at(index)->GetNDof();

    local2globals.at(index).resize(array_length);

    std::copy(l2g,
              l2g + array_length, // pointer arithmetic!
              local2globals.at(index).begin());
  }

  // 1) setting up the matrix
  // get a rough overview and reserve some space
  matrix_.n = n_dofs_global;
  matrix_.nz_loc = 0;
  {
    size_t expected_entries = bmatrix.get_n_total_entries();
    matrix_.irn_loc.reserve(expected_entries);
    matrix_.jcn_loc.reserve(expected_entries);
    matrix_.a_loc.reserve(expected_entries);
  }

  // loop through the block matrix
  // will always hold the transposed state of the last treated block
  bool transp;

  // will always hold the last treated block
  std::shared_ptr<const TMatrix> block = bmatrix.get_block(0, 0, transp);

  for (size_t cell_row = 0; cell_row < nComms; ++cell_row) // loop over rows
  {
    const std::vector<int>& masters = dof_masters.at(cell_row); // fix masters
    size_t row_shift = shifts.at(cell_row);                     // fix shift
    const std::vector<int>& row_l2g = local2globals.at(cell_row); // fix local-to-global mapping

    for (size_t cell_col = 0; cell_col < nComms; ++cell_col) // loop over columns
    {
      size_t col_shift = shifts.at(cell_col); // fix shift
      const std::vector<int>& col_l2g = local2globals.at(cell_col); // fix local-to-global mapping

      // fetch all the info from the current block
      block = bmatrix.get_block(cell_row, cell_col, transp);
      bool is_bbt_type = (block->get_sparse_type() == SparsityType::B_TIMES_BT);
      const int* row_ptr = block->get_row_ptr();
      const int* k_col = block->get_vector_columns();
      const double* entries = block->GetEntries();

      if (!transp) // block is stored non-transposed
      {
        for (size_t row = 0; (int)row < block->get_n_rows(); ++row)
        {
          // active rows

          size_t row_start = row_ptr[row];
          size_t row_end = row_ptr[row + 1];

          for (size_t k = row_start; k < row_end; ++k)
          {
            size_t col = k_col[k];
            double entry = entries[k];

            // gather the conditions under which (AND) we will to put the entry
            // into the sparsity structure on this process
            bool is_nonzero = (entry != 0.0);

            bool is_in_master_row = masters.at(row) == mpi_rank;

            bool is_in_active_row = true; // true for standard block matrix

            if (is_block_fe_matrix)
            {
              is_in_active_row = (row < static_cast<const BlockFEMatrix&>(bmatrix).get_n_row_actives(cell_row));
            }

            // If this is BB^T case, the matrix is assumed to be
            // in additive storage and all entries are made count.
            if (is_bbt_type)
            {
              is_in_master_row = true;
            }

            if (is_nonzero && is_in_master_row && is_in_active_row)
            {
              // put entry into the new mumps matrix

              if (loc_to_seq.empty())
              {
                matrix_.irn_loc.push_back(row_shift + row_l2g.at(row) + 1);
                matrix_.jcn_loc.push_back(col_shift + col_l2g.at(col) + 1);
              }
              else
              {
                matrix_.irn_loc.push_back(row_shift + loc_to_seq.at(cell_row).at(row) + 1);
                matrix_.jcn_loc.push_back(col_shift + loc_to_seq.at(cell_col).at(col) + 1);
              }

              matrix_.a_loc.push_back(entry);
              ++matrix_.nz_loc;
            }
          }
        }
      }// end non-transposed case
      else if (transp) // block is stored transposed
      {
        if (is_block_fe_matrix)
        {
          // static cast is okay, because we are sure to deal with a BlockFEMatrix
          if (std::static_pointer_cast<const FEMatrix>(block)->get_n_active_rows() != block->get_n_rows()) // check for security
          {
            ErrThrow("This block (", cell_row, ",", cell_col, ") has test space "
                     "non-actives, it should never have been stored in transposed state!");
          }
        }

        for (int row = 0; row < block->get_n_rows(); ++row)
        {
          size_t row_start = row_ptr[row];
          size_t row_end = row_ptr[row+1];

          for (size_t k = row_start; k < row_end; ++k)
          {
            size_t col = k_col[k];
            double entry = entries[k];

            // gather the conditions under which (AND) we will to put the entry
            // into the sparsity structure on this process
            bool is_nonzero = (entry != 0.0);

            bool is_in_master_row = masters.at(col) == mpi_rank;

            bool is_in_active_row = true; // true for standard block matrix

            if(is_block_fe_matrix)
            {
              is_in_active_row = (col < static_cast<const BlockFEMatrix&>(bmatrix).get_n_row_actives(cell_row));
            }

            if (is_nonzero && is_in_master_row && is_in_active_row)
            {
              // fill entry in the sparsity structure

              if (loc_to_seq.empty())
              {
                matrix_.irn_loc.push_back(row_shift + row_l2g.at(col) + 1);
                matrix_.jcn_loc.push_back(col_shift + col_l2g.at(row) + 1);
              }
              else
              {
                matrix_.irn_loc.push_back(row_shift + loc_to_seq.at(cell_row).at(col) + 1);
                matrix_.jcn_loc.push_back(col_shift + loc_to_seq.at(cell_col).at(row) + 1);
              }

              matrix_.a_loc.push_back(entry);
              ++matrix_.nz_loc;
            }
          }
        }
      }

      if (is_block_fe_matrix)
      {
        // on diagonal blocks, set the ones on diagonals of non-active rows
        // (same code both cases, but loop does only make iterations for
        // non-transp case, because transposed storage of blocks with testspace
        // non-actives is not allowed in BlockFEMatrix)
        if (cell_row == cell_col)
        {
          for (int row = std::static_pointer_cast<const FEMatrix>(block)->get_n_active_rows();
            row < block->get_n_rows(); ++row)
          {
            if (masters.at(row) == mpi_rank)
            {
              // only add the entry if the current process is master of its row
              // put a 1 on the diagonal

              if (loc_to_seq.empty())
              {
                matrix_.irn_loc.push_back(row_shift + row_l2g.at(row) + 1);
                matrix_.jcn_loc.push_back(row_shift + col_l2g.at(row) + 1);
              }
              else
              {
                matrix_.irn_loc.push_back(row_shift + loc_to_seq.at(cell_row).at(row) + 1);
                matrix_.jcn_loc.push_back(row_shift + loc_to_seq.at(cell_col).at(row) + 1);
              }

              matrix_.a_loc.push_back(1);
              ++matrix_.nz_loc;
            }
          }
        }// end treating dirichlet rows
      }
    }// end loop over columns
  }// end loop over rows

  if (is_block_fe_matrix)
  {
    // if the matrix comes from an enclosed flow problem,
    // an internal pressure row correction is necessary
    if (static_cast<const BlockFEMatrix&>(bmatrix).pressure_projection_enabled())
    {
      Output::root_info<5>("Pressure Projection", "MumpsWrapper applying pressure row correction.");
      pressure_row_correction(pres0);
    }
  }
}

void MumpsWrapper::pressure_row_correction(
    std::vector<double> pres0)
{
  if (!pres0.empty())
  {
    const TFESpace3D* pres_space = this->comms_.at(3)->get_fe_space();
    double p0_x = pres0.at(0);
    double p0_y = pres0.at(1);
    double p0_z = pres0.at(2);

    int pressure_dof_to_correct = -1;
    for (int d = 0; d < pres_space->get_n_dof(); ++d)
    {
      double x,y,z;
      double tol = 1e-6;

      pres_space->GetDOFPosition(d, x, y, z);
      if(    std::abs(p0_x - x) < tol
          && std::abs(p0_y - y) < tol
          && std::abs(p0_z - z) < tol)
      {
        pressure_dof_to_correct = d;
        break;
      }
    }

    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int n_velo_dofs_local = comms_.at(0)->GetN_Master(); // master dofs
    int sendbuf[1] = {n_velo_dofs_local};
    int recvbuf[1] = {0};

    // gather total number of velo dofs in root 0
    MPI_Allreduce(sendbuf, recvbuf, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    bool i_am_master = false;
    if (pressure_dof_to_correct != -1)
    {
      i_am_master = (comms_.at(3)->GetMaster()[pressure_dof_to_correct] == rank);
    }

    if (i_am_master)
    {
      // correct that entire row on "my" process

      int n_velo_dofs_global = 3 * recvbuf[0]; // we're in 3D, thus multiply by 3

      // go through the matrix and remove all entries
      bool found_diagonal = false;
      for (size_t i = 0; i < matrix_.nz_loc; ++i)
      {
        if (matrix_.irn_loc.at(i) == n_velo_dofs_global + pressure_dof_to_correct + 1)
        {
          if (matrix_.jcn_loc.at(i) != n_velo_dofs_global + pressure_dof_to_correct + 1)
          {
            // off-diagonal entry - set to zero!
            matrix_.a_loc.at(i) = 0;
          }
          else
          {
            // diagonal entry - put to one!
            matrix_.a_loc.at(i) = 1;
            found_diagonal = true;
          }
        }
      }

      if (!found_diagonal)
      {
        matrix_.irn_loc.push_back(n_velo_dofs_global + pressure_dof_to_correct + 1);
        matrix_.jcn_loc.push_back(n_velo_dofs_global + pressure_dof_to_correct + 1);
        matrix_.a_loc.push_back(1);
        matrix_.nz_loc = matrix_.nz_loc + 1;
      }
    }
  }
  else
  {
    int size, my_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int n_velo_dofs_local = comms_.at(0)->GetN_Master(); // master dofs
    int sendbuf[1] = {n_velo_dofs_local};
    int recvbuf[1] = {0};

    // gather total number of velo dofs in root 0
    MPI_Reduce(sendbuf, recvbuf, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // due to the used global ordering it is process 0 that holds the
    // globally first pressure row. that row should to be set to unit vector
    if (my_rank == 0)
    {
      int n_velo_dofs_global = 3 * recvbuf[0]; // we're in 3D, thus multiply by 3
      int first_p_rw = n_velo_dofs_global;

      // go through the matrix and remove all entries
      for (size_t i = 0; i < matrix_.nz_loc; ++i)
      {
        if (matrix_.irn_loc.at(i) == first_p_rw)
        {
          if (matrix_.jcn_loc.at(i) != first_p_rw)
          {
            // off-diagonal entry - set to zero!
            matrix_.a_loc.at(i) = 0;
          }
          else
          {
            // diagonal entry - put to one!
            matrix_.a_loc.at(i) = 1;
          }
        }
      }
    }
  }

}

void MumpsWrapper::gather_vector(
    double* GlobalArray, double *LocalArray, int LocalSize, int root) const
{
  MPI_Comm Comm = MPI_COMM_WORLD;
  int size;
  MPI_Comm_size(Comm, &size);

  int *displ         = new int[size];
  int *N_ElementsAll = new int[size];

  // determine how many elements to receive per process - N_ElementsAll
  MPI_Allgather(&LocalSize, 1, MPI_INT, N_ElementsAll, 1, MPI_INT, Comm);

  displ[0] = 0;
  for (int i = 1; i < size; i++)
  {
    displ[i] = displ[i - 1] + N_ElementsAll[i - 1];
  }

  MPI_Gatherv(LocalArray, LocalSize, MPI_DOUBLE, // send
              GlobalArray, N_ElementsAll, displ, MPI_DOUBLE, // receive
              root, Comm); // control

  delete[] displ;
  delete[] N_ElementsAll;
}

void MumpsWrapper::scatter_vector(
    double *GlobalArray, double *LocalArray, int LocalSize, int root) const
{
  MPI_Comm Comm = MPI_COMM_WORLD;
  int size;
  MPI_Comm_size(Comm, &size);

  int *displ         = new int[size];
  int *N_ElementsAll = new int[size];
  MPI_Allgather(&LocalSize, 1, MPI_INT, N_ElementsAll, 1, MPI_INT, Comm);

  displ[0] = 0;
  for (int i = 1; i < size; i++)
  {
    displ[i] = displ[i - 1] + N_ElementsAll[i - 1];
  }

  MPI_Scatterv(GlobalArray, N_ElementsAll, displ, MPI_DOUBLE, // send
               LocalArray, LocalSize, MPI_DOUBLE,             // receive
               root, Comm);                                   // control

  delete[] displ;
  delete[] N_ElementsAll;
}

void MumpsWrapper::check_input_matrix(const BlockMatrix& bmatrix)
{
  // Check that, if there is a B*B^T-matrix contained, it
  // is the only matrix and it is non-transposed.
  for (size_t row = 0; row < bmatrix.get_n_cell_rows(); ++row)
  {
    for(size_t col = 0; col < bmatrix.get_n_cell_rows(); ++col)
    {
      bool transp;

      if (bmatrix.get_block(row,col,transp)->get_sparse_type() == SparsityType::B_TIMES_BT)
      {
        if (bmatrix.get_n_cell_rows() != 1 || bmatrix.get_n_cell_rows() != 1 || transp)
        {
          // currently I do not feel like thinking about this case
          ErrThrow("SparsityType::B_TIMES_BT matrix found in "
              "non-1x1 BlockMatrix (or in transposed state).");
        }
        else
        {
          Output::root_info("MumpsWrapper", "SparsityType::B_TIMES_BT matrix found.");
        }
      }
    }
  }
}

#else // PARMOON_WITH_MUMPS

MumpsWrapper::MumpsWrapper(const BlockFEMatrix&, std::vector<double>)
{
  ErrThrow("ParMooN has been compiled without mumps, therefore a MumpsWrapper "
           "can not be created.");
}

MumpsWrapper::MumpsWrapper(const BlockMatrix&,
                           std::vector<const TParFECommunicator3D*>,
                           std::vector<double>,
                           std::vector<std::vector<int>>)
{
  ErrThrow("ParMooN has been compiled without mumps, therefore a MumpsWrapper "
           "can not be created.");
}

MumpsWrapper::~MumpsWrapper()
{
}

void MumpsWrapper::enable_distributed_solution(int, double)
{
  ErrThrow("ParMooN has been compiled without mumps, therefore "
           "MumpsWrapper::enable_distributed_solution can not be called.");
}

void MumpsWrapper::write_matrix_distributed(const std::string&)
{
  ErrThrow("ParMooN has been compiled without mumps, therefore "
           "MumpsWrapper::write_matrix_distributed can not be called.");
}

void MumpsWrapper::solve(const BlockVector&, BlockVector&, int)
{
  ErrThrow("ParMooN has been compiled without mumps, therefore "
           "MumpsWrapper::solve can not be called.");
}

#endif // PARMOON_WITH_MUMPS

#endif // _MPI