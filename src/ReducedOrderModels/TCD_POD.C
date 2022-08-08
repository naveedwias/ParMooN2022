/** **************************************************************************** 
*
* @name       TCD_POD
* @brief      Computation of POD basis - specific routines for TCD problems
*
*******************************************************************************/

#include "TCD_POD.h"
#include "Database.h"
#include "MainUtilities.h"
#include <sys/stat.h>

#ifdef __2D__
#include "Assemble2D.h"
#include "SquareMatrix2D.h"
#else
#include "Assemble3D.h"
#include "SquareMatrix3D.h"
#endif

#ifdef _MPI
#include "ParFECommunicator3D.h"
#endif

/* ************************************************************************** */
template <int d>
ParameterDatabase TCD_POD<d>::set_pod_basis_database(
                                              const ParameterDatabase& param_db)
{
  ParameterDatabase db("Parameter database for POD basis computation");

  db.merge(LocalAssembling<d>::default_local_assembling_database(), true);

  db.merge(ParameterDatabase::default_output_database(), true);
  db.merge(POD::default_pod_database(), true);
  db.merge(param_db, false);

  // change name to avoid overwritting standard FEM output
  // TODO: unification between naming system (outfile, output_basename,
  //       snaps_basename, pod_basename)
  std::string output_pod_basename = db["output_basename"].get<std::string>()
                                    + "_pod";
  db["output_basename"].set(output_pod_basename, false);

  return db;
}

/* ************************************************************************** */
template <int d>
TCD_POD<d>::TCD_POD(TCollection& coll, const ParameterDatabase& param_db)
  : TCD_POD<d>(coll, param_db, Example_TimeCD(param_db))
{
}

/* ************************************************************************** */
template <int d>
TCD_POD<d>::TCD_POD(TCollection& coll,
                    const ParameterDatabase& param_db,
                    const Example_TimeCD&    ex)
  : fe_space(new FESpace(&coll, "space", ex.get_bc(0),
                         TDatabase::ParamDB->ANSATZ_ORDER)),
    pod_c(param_db),
    db(TCD_POD<d>::set_pod_basis_database(param_db)),
    example(ex),
    outputWriter(db)
{
  // create directory db["pod_directory"]
  std::string directory_name = this->db["pod_directory"].value_as_string();
  mkdir(directory_name.c_str(), 0777);

#ifdef __3D__
  this->gramian_matrix = BlockFEMatrix::CD3D(fe_space);
#else
  this->gramian_matrix = BlockFEMatrix::CD2D(fe_space);
#endif
  this->pod_mode = BlockVector(this->gramian_matrix, false);
  this->fe_function = FEFunction(this->fe_space, "c",
                                 this->pod_mode.get_entries());

  this->set_parameters();

  this->output_problem_size_info();

  this->outputWriter.add_fe_function(&this->fe_function);
}

/* ************************************************************************** */
template <int d>
void TCD_POD<d>::set_parameters()
{
  //set problem_type to Time_CD if not yet set
  if(!db["problem_type"].is(2))
  {
    if (db["problem_type"].is(0))
    {
      db["problem_type"] = 2;
    }
    else
    {
      Output::warn<2>("The parameter problem_type doesn't correspond to "
                      "Time_CD. It is now reset to the correct value for "
                      "Time_CD (=2).");
      db["problem_type"] = 2;
    }
  }

  // an error when using ansatz order 0
  if(TDatabase::ParamDB->ANSATZ_ORDER == 0)
  {
    throw std::runtime_error("Ansatz order 0 is no use in convection diffusion "
        "reaction problems! (Vanishing convection and diffusion term).");
  }

  // For assembling of the Gramian matrix for the computation of
  // the POD basis only the "Galerkin" part of the matrices is needed.
  // TDatabase::ParamDB->DISCTYPE = GALERKIN;
  if(!db["space_discretization_type"].is("galerkin"))
  {
    db["space_discretization_type"] = "galerkin";
  }

  if(db["pod_inner_product"].is("H1"))
  {
    ErrThrow("For Convection-Diffusion Problems, the parameter "
             "'pod_inner_product' does not support the value 'H1'.");
  }
}

/* ************************************************************************** */
template <int d>
void TCD_POD<d>::assemble_gramian()
{
  Output::print<1>("[pod_inner_product = L2]: "
                   "Assembling the gramian matrix for POD computation...");

  std::vector<const FEFunction*> feFunction = {&this->fe_function};
  LocalAssembling<d> la(this->db, LocalAssembling_type::TCDMassOnly,
                        feFunction, example.get_coeffs());

  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  int nFESpaces = 1;
  const FESpace *fe_space = this->fe_space.get();

  auto block = this->gramian_matrix.get_blocks_uniquely()[0].get();
  /**
     @attention
     we get block [1] because the TCD Mass Matrix functions
     assembles Mass as block[1] (block[0] is reserved for stiffness)
  **/
  const int nsqMat = 2;
  SquareMatrixD *sqMat[nsqMat];
  sqMat[0] = nullptr;
  sqMat[1] = reinterpret_cast<SquareMatrixD*>(block);
  sqMat[1]->reset();
  int nRhs = 0;

  auto * bound_cond = this->fe_space->get_boundary_condition();
  auto * bound_val = this->example.get_bd(0);
#ifdef __2D__
  Assemble2D(
#else
  Assemble3D(
#endif
    nFESpaces, &fe_space, nsqMat, sqMat, 0, nullptr, nRhs, nullptr,
    &fe_space, &bound_cond, &bound_val, la);

  Output::print<1>("... assembling done.");
}

/* ************************************************************************** */
template <int d>
void TCD_POD<d>::assemble_stiff(BlockFEMatrix& stiff_matrix)
{
  Output::print<1>("Assembling the stiffness matrix of the POD basis...");

  std::vector<const FEFunction*> feFunction = {&this->fe_function};
  LocalAssembling<d> la(this->db, LocalAssembling_type::TCDStiffOnly,
                        feFunction, example.get_coeffs());

  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  int nFESpaces = 1;
  const FESpace *fe_space = this->fe_space.get();

  auto block = stiff_matrix.get_blocks_uniquely()[0].get();
  /**
     @attention
     we get block [0] because the TCD Stiff Matrix functions
     assembles Stiff as block[0]
  **/
  const int nsqMat = 2;
  SquareMatrixD *sqMat[nsqMat];
  sqMat[0] = reinterpret_cast<SquareMatrixD*>(block);
  sqMat[0]->reset();
  sqMat[1] = nullptr;
  int nRhs = 0;

  auto * bound_cond = this->fe_space->get_boundary_condition();
  auto * bound_val = this->example.get_bd(0);
#ifdef __2D__
  Assemble2D(
#else
  Assemble3D(
#endif
    nFESpaces, &fe_space, nsqMat, sqMat, 0, nullptr, nRhs, nullptr,
    &fe_space, &bound_cond, &bound_val, la);

  Output::print<1>("... assembling done.");
}

/* ************************************************************************** */
template <int d> 
void TCD_POD<d>::write_S2()
{
  BlockFEMatrix stiff_matrix;
#ifdef __3D__
  stiff_matrix   = BlockFEMatrix::CD3D(fe_space);
#else
  stiff_matrix   = BlockFEMatrix::CD2D(fe_space);
#endif
  this->assemble_stiff(stiff_matrix);

  bool transpose = true;
  std::shared_ptr<DenseMatrix> tmp_mat = stiff_matrix.get_combined_matrix()->
                                         multiply(pod_c.get_basis(), transpose);
  DenseMatrix stiff_matrix_pod = *(pod_c.get_basis()->multiply(tmp_mat.get()));

  double s_norm2 = stiff_matrix_pod.norm2();
  Output::print<1>("Spectral norm of the stiffness matrix of the POD-basis: ",
                   s_norm2);

  std::string s2_filename = pod_c.get_basename() + "s2";
  std::ofstream ofile;
  ofile.open( s2_filename.c_str() , std::ios::out | std::ios::trunc );
  ofile << setprecision( 12 );
  ofile << s_norm2 << endl;
  ofile.close();
}

/* ************************************************************************** */
template <int d> 
void TCD_POD<d>::output_problem_size_info() const
{
  const TCollection *coll = this->fe_space->GetCollection();
  double hMin, hMax;
  coll->GetHminHmax(&hMin, &hMax);
  int n_cells = coll->GetN_Cells();
  int n_dof = this->fe_space->get_n_dof();
  Output::stat("TimeConvectionDiffusion",
                "Mesh data and problem size on finest grid");
  Output::dash("n cells     : ", setw(13), n_cells);
  Output::dash("h(min, max) : ", setw(13), hMin, " ", setw(13), hMax);
  Output::dash("n dofs      : ", setw(13), n_dof);
  Output::dash("n dof active: ", setw(13), this->fe_space->get_n_active());
}

/* ************************************************************************** */
template <int d>
void TCD_POD<d>::compute_pod_basis()
{
  double t_pod_start = GetTime();

  if(this->db["pod_inner_product"].is("L2"))
  {
    this->assemble_gramian();
    pod_c.set_gramian_ptr(this->gramian_matrix.get_combined_matrix());
  }

  pod_c.compute_basis( this->fe_space->get_n_dof());
  Output::print<1>("POD basis computed in ", GetTime()-t_pod_start, "s.");
}

/* ************************************************************************** */
template <int d>
void TCD_POD<d>::output()
{
  if(this->fe_space->get_n_dof() != pod_c.get_basis()->getNColumns())
  {
    ErrThrow("Current FE space dimension does not coincide with the number "
             "of dof of POD basis.\nDimension of FE space : ",
             this->fe_space->get_n_dof(),
             "\nDOF of POD basis      : ", pod_c.get_basis()->getNColumns());
  }

  // write spectral norm of the stiffness matrix of the POD basis
  this->write_S2();

  int shift = 0;
  /* write averages of snapshots into a vtk-file */
  if(this->db["pod_fluctuations_only"])
  {
    for(int j=0; j< this->fe_space->get_n_dof(); ++j)
    {
      this->pod_mode.get_entries()[j] = pod_c.get_snaps_avr()[j];
    }
    Output::print<1>("Writing vtk for snapshots average into ",
                      this->db["output_basename"], 0 ,".vtk");
    this->outputWriter.write(0);
    shift = 1;
  }

  for(int i=0; i < std::min(pod_c.get_rank(), 10); ++i)
  {
    for(int j=0; j<this->fe_space->get_n_dof(); ++j)
    {
      this->pod_mode.get_entries()[j] = pod_c.get_basis()->getEntry(i,j);
    }
    Output::print<1>("Writing vtk for POD basis into ",
                     this->db["output_basename"], i+shift ,".vtk");
    this->outputWriter.write(i+shift);
  }
}

#ifdef __3D__
template class TCD_POD<3>;
#else
template class TCD_POD<2>;
#endif
