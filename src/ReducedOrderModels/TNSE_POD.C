/** **************************************************************************** 
*
* @name       TNSE_POD
* @brief      Computation of POD basis - specific routines for TNSED problems
*
*******************************************************************************/

#include "CD_local_assembling_routines.h"
#include "Database.h"
#include "MainUtilities.h"
#include "NSE_local_assembling_routines.h"
#include "ROM_local_assembling_routines.h"
#include "Time_NSE_local_assembling_routines.h"
#include "TNSE_POD.h"
#include "Solver.h"

#ifdef __2D__
#include "Assemble2D.h"
#include "SquareMatrix2D.h"
#else
#include "Assemble3D.h"
#include "SquareMatrix3D.h"
#endif

#include <sys/stat.h>

/* ************************************************************************** */
template <int d>
ParameterDatabase TNSE_POD<d>::set_pod_basis_database(
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
TNSE_POD<d>::TNSE_POD(TCollection& coll, const ParameterDatabase& param_db)
  : TNSE_POD<d>(coll, param_db, Example_TimeNSE(param_db))
{
}

/* ************************************************************************** */
template <int d>
TNSE_POD<d>::TNSE_POD(TCollection&             coll,
                      const ParameterDatabase& param_db,
                      const Example_TimeNSE&   ex)
  : velocity_space(new FESpace(&coll, "velocity space", ex.get_bc(0),
                               TDatabase::ParamDB->VELOCITY_SPACE)),
    pressure_space(new FESpace(&coll, "pressure space", ex.get_bc(d),
                               TDatabase::ParamDB->PRESSURE_SPACE)),
    db(TNSE_POD<d>::set_pod_basis_database(param_db)),
    example(ex),
    outputWriter(db)
{
  // create directory db["pod_directory"]
  std::string directory_name = this->db["pod_directory"].value_as_string();
  mkdir(directory_name.c_str(), 0777);

  BlockFEMatrix diag_b_matrix;

#ifdef __3D__
  diag_b_matrix = NSE3D_Diag(velocity_space, pressure_space);
#else
  diag_b_matrix = NSE2D_Diag(velocity_space, pressure_space);
#endif

  this->pod_mode = BlockVector(diag_b_matrix, false);

  velocity_function = FEVectFunct(velocity_space, "u",
                                  pod_mode.block(0), d);
  pressure_function = FEFunction(pressure_space, "p",
                                 pod_mode.block(d));

  this->set_parameters();

  this->output_problem_size_info();

  this->outputWriter.add_fe_vector_function(&this->velocity_function);
  this->outputWriter.add_fe_function(&this->pressure_function);

}

/* ************************************************************************** */
template <int d>
void TNSE_POD<d>::set_parameters()
{
  // set problem_type to Time_NSE if not yet set
  if (!db["problem_type"].is(6))
  {
    if (db["problem_type"].is(0))
    {
      db["problem_type"] = 6;
    }
    else
    {
      Output::warn<2>("The parameter problem_type doesn't correspond to "
                      "TimeNavierStokes. It is now reset to the correct value "
                      "for TimeNavierStokes (=6).");
      db["problem_type"] = 6;
    }
  }

  // For assembling of the Gramian matrix for the computation of
  // the POD basis only the "Galerkin" part of the matrices is needed.
  if (!db["space_discretization_type"].is("galerkin"))
  {
    db["space_discretization_type"] = "galerkin";
  }
}

/* ************************************************************************** */
template <int d>
void TNSE_POD<d>::assemble_gramian(BlockFEMatrix* mat_ptr, bool is_velocity, std::string inner_product)
{
  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  using BdValFunct = typename Template_names<d>::BoundaryValuesFunction;
  using BdCondFunct = typename Template_names<d>::BoundaryConditionFunction;
  using namespace std::placeholders;

  LocalAssembling_type type;
  std::vector<const FESpace*> spaces_mat;
  std::vector<const FESpace*> spaces_rhs;
  std::vector<const FEFunction*> fefunctions;
  std::vector<SquareMatrixD*> sqMatrices;
  std::vector<double*> rhs_array;
  std::vector<BdCondFunct*> boundCondition;
  std::vector<BdValFunct*> boundValues;

  if (is_velocity)
  {
    type = LocalAssembling_type::TimeNavierStokesMass;
    fefunctions.resize(d);
    for (int i = 0; i < d; ++i)
    {
      fefunctions[i] = velocity_function.GetComponent(i).release();
    }
  }
  else
  {
    type = LocalAssembling_type::TCDMassOnly;
    fefunctions.resize(1);
    fefunctions[0] = &pressure_function;
  }

  LocalAssembling<d> la(this->db, type, fefunctions, example.get_coeffs());

  std::vector<AssembleFctParam> la_modif;
  auto block = mat_ptr->get_blocks_uniquely()[0].get();
  double coeff_one = 1.;

  if (is_velocity)
  {
    spaces_mat.resize(2);
    spaces_mat[0] = velocity_space.get();
    spaces_mat[1] = pressure_space.get();
    spaces_rhs.resize(d + 1);

    for (int i = 0; i < d; ++i)
    {
      spaces_rhs[i] = velocity_space.get();
    }
    spaces_rhs[d] = pressure_space.get();

    sqMatrices.resize(1);
    sqMatrices[0] = reinterpret_cast<SquareMatrixD*>(block);

    rhs_array.resize(d + 1, nullptr);
    boundCondition.resize(d);
    for (int i = 0; i < d; ++i)
    {
      boundCondition[i] = velocity_space->get_boundary_condition();
    }

    boundValues.resize(d);
    for (int i = 0; i < d; ++i)
    {
      boundValues[i] = this->example.get_bd(i);
    }

    if (inner_product == "L2")
    {
      la_modif.push_back(std::bind(NSMassMatrixSingle<d>,
                                   _1, _2, _3, _4, _5, _6, _7, _8 ) );
    }
    else if (inner_product == "H1")
    {
      la_modif.push_back(std::bind(NSLaplaceGradGradSingle<d>,
                                   _1, &coeff_one, _3, _4, _5, _6, _7, _8 ) );
    }
    else
    {
      ErrThrow("pod_inner_product: ", inner_product, " not supported");
    }
  }
  else
  {
    spaces_mat.resize(1);
    spaces_mat[0] = pressure_space.get();
    sqMatrices.resize(2);
    sqMatrices[0] = nullptr;
    sqMatrices[1] = reinterpret_cast<SquareMatrixD*>(block);
    boundCondition.resize(1);
    boundCondition[0] = pressure_space->get_boundary_condition();
    boundValues.resize(1);
    boundValues[0] = this->example.get_bd(d);

    if (inner_product == "L2")
    {
      la_modif.push_back(std::bind(TCDMass<d>,
                                   _1, _2, _3, _4, _5, _6, _7, _8 ) );
    }
    else if (inner_product == "H1")
    {
      la_modif.push_back(std::bind(TCDDiffOnlyScale<d>,
                                   _1, _2, _3, 1., _5, _6, _7, _8 ) );
    }
    else
    {
      ErrThrow("pod_inner_product: ", inner_product, " not supported");
    }
  }

  la.replace_local_assembling(la_modif);
  for (auto m : sqMatrices)
  {
    if (m != nullptr)
    {
      m->reset();
    }
  }

#ifdef __2D__
    Assemble2D(
#else
    Assemble3D(
#endif
               spaces_mat.size(), spaces_mat.data(), sqMatrices.size(),
               sqMatrices.data(), 0, nullptr,
               rhs_array.size(), rhs_array.data(), spaces_rhs.data(),
               boundCondition.data(), boundValues.data(), la);
}

/* ************************************************************************** */
template <int d> 
void TNSE_POD<d>::output_problem_size_info() const
{
  double hMin, hMax;
  auto coll   = velocity_space->GetCollection();
  int n_cells = coll->GetN_Cells();
  coll->GetHminHmax(&hMin, &hMax);

  size_t n_dof_velocity = d*velocity_space->get_n_dof();
  size_t n_dof_pressure = pressure_space->get_n_dof();
  size_t n_Total        = n_dof_velocity + n_dof_pressure;
  size_t n_Active       = d*velocity_space->get_n_active();

  Output::stat("TimeNavierStokes",
                "Mesh data and problem size on finest grid");
  Output::print("n cells     : ", setw(10), n_cells);
  Output::print("h(min,max)  : ", setw(12), hMin ," ", setw(12), hMax);
  Output::print("dof Velocity: ", setw(10), n_dof_velocity);
  Output::print("dof Pressure: ", setw(10), n_dof_pressure);
  Output::print("dof Total   : ", setw(10), n_Total);
  Output::print("active dof  : ", setw(10), n_Active);
}

/* ************************************************************************** */
template <int d>
void TNSE_POD<d>::compute_pod_basis(const ParameterDatabase& param_db)
{
  // container to save values for output
  std::vector<std::vector<double>> snaps_avr(2);
  std::vector<DenseMatrix> basis(2, DenseMatrix(0,0));

  auto vs_db = Solver<>::default_solver_database();
  auto ps_db = Solver<>::default_solver_database();

  BlockFEMatrix vs_mat;
  BlockFEMatrix ps_mat;

  if (param_db["snaps_time_derivative_projection"])
  {
    std::string v_db_name = "POD Database - velocity projection Solver";
    std::string p_db_name = "POD Database - pressure projection Solver";

    if (param_db.has_nested_database(v_db_name))
    {
      vs_db.merge(param_db.get_nested_database(v_db_name), false);
    }
    else
    {
      Output::warn("POD", "Missing velocity solver "
                   "database. You should provide a nested database named ",
                   v_db_name, ". Now a default database is used. Is this "
                   "intended?");
    }

    if (param_db.has_nested_database(p_db_name))
    {
      ps_db.merge(param_db.get_nested_database(p_db_name), false);
    }
    else
    {
      Output::warn("POD", "Missing pressure solver "
                   "database. You should provide a nested database named ",
                   p_db_name, ". Now a default database is used. Is this "
                   "intended?");
    }

#ifdef __3D__
    vs_mat = BlockFEMatrix::Mass_NSE3D(velocity_space);
    ps_mat = BlockFEMatrix::CD3D(pressure_space);
#else
    vs_mat = BlockFEMatrix::Mass_NSE2D(velocity_space);
    ps_mat = BlockFEMatrix::CD2D(pressure_space);
#endif

    assemble_gramian(&vs_mat, true, "L2");
    assemble_gramian(&ps_mat, false, "L2");

    vs_mat = vs_mat.get_sub_blockfematrix(0, d - 1);
  }

  for (auto& compo : data_suffix)
  {
    bool is_velocity = (compo == "u") ? true : false;
    std::string name = compo + "-";
    int length = is_velocity ? d * velocity_space->get_n_dof()
                             : pressure_space->get_n_dof();

    POD pod_tmp(param_db, compo);

    if (!this->db["pod_inner_product"].is("euclidean"))
    {
#ifdef __3D__
      BlockFEMatrix gramian_matrix = is_velocity
                                     ? BlockFEMatrix::Mass_NSE3D(velocity_space)
                                     : BlockFEMatrix::CD3D(pressure_space);
#else
      BlockFEMatrix gramian_matrix = is_velocity
                                     ? BlockFEMatrix::Mass_NSE2D(velocity_space)
                                     : BlockFEMatrix::CD2D(pressure_space);
#endif
      assemble_gramian(&gramian_matrix, is_velocity, db["pod_inner_product"]);
      pod_tmp.set_gramian_ptr(gramian_matrix.get_combined_matrix());
    }

    if (db["snaps_time_derivative_projection"])
    {
      if (is_velocity)
      {
        pod_tmp.set_projection(vs_mat, vs_db);
      }
      else
      {
        pod_tmp.set_projection(ps_mat, ps_db);
      }
    }

    double t_pod_start = GetTime();
    pod_tmp.compute_basis(length);
    Output::print<1>(name, "POD-basis computed in ",
                     GetTime()-t_pod_start, "s.");

    // save value for output
    if (length != pod_tmp.get_basis()->getNColumns())
    {
      ErrThrow("Current FE space dimension does not coincide with the number "
               "of dof of POD basis.\nDimension of FE space : ", length,
               "\nDOF of POD basis      : ",pod_tmp.get_basis()->getNColumns());
    }

    int i = is_velocity ? 0 : 1;
    snaps_avr[i] = pod_tmp.get_snaps_avr();
    basis[i]     = *pod_tmp.get_basis();
  }

  output(snaps_avr, basis);
}

/* ************************************************************************** */
template <int d>
void TNSE_POD<d>::output(const std::vector<std::vector<double>>& snaps_avr,
                         const std::vector<DenseMatrix>&         basis)
{
  int shift = 0;
  // write averages of snapshots into a vtk-file
  if (this->db["pod_fluctuations_only"])
  {
    int index = 0;
    for (int i = 0; i < (int)snaps_avr.size(); ++i)
    {
      for (int j = 0; j < (int)snaps_avr[i].size(); ++j)
      {
        this->pod_mode.get_entries()[index] = snaps_avr[i][j];
        index++;
      }
    }
    Output::print<1>("Writing vtk for snapshots average into ",
                      this->db["output_basename"], 0, ".vtk");
    this->outputWriter.write(0);
    shift = 1;
  }

  // calculate maximal rank
  int rank_max = 0;
  for (int k = 0; k < (int)basis.size(); ++k)
  {
    rank_max = std::max(rank_max, basis[k].getNRows());
  }

  // write the first min(rank_max, 10) POD-basis functions
  for (int i = 0; i < std::min(rank_max, 10); ++i)
  {
    int index = 0;

    for (int k = 0; k < (int)basis.size(); ++k)
    {
      for (int j = 0; j < basis[k].getNColumns(); ++j)
      {
        if (i < basis[k].getNRows())
        {
          this->pod_mode.get_entries()[index] = basis[k].getEntry(i,j);
        }
        else
        {
          this->pod_mode.get_entries()[index] = 0.;
        }

        index++;
      }
    }
    Output::print<1>("Writing vtk for POD basis into ",
                     this->db["output_basename"], i+shift ,".vtk");
    this->outputWriter.write(i+shift);
  }
}



#ifdef __2D__
/**
  * Named constructor for a matrix taking the block structure
  *
  * ( A  0  0 )
  * ( 0  A  0 )
  * ( 0  0  C )
  *
  * @param velocity The velocity finite element space.
  * @param pressure The pressure finite element space.
  * @return A newly constructed BlockFEMatrix for NSE3D POD-ROM problems
  */
BlockFEMatrix NSE2D_Diag(std::shared_ptr<const TFESpace2D> velocity,
                         std::shared_ptr<const TFESpace2D> pressure)
{
  BlockFEMatrix my_matrix({velocity, velocity, pressure});

  //create new blocks with correct structures filled with 0
  FEMatrix velo_velo(velocity, velocity);                // A block
  FEMatrix velo_velo_zero(velocity, velocity, true);     // velocity zero block
  FEMatrix pressure_velo_zero(pressure, velocity, true); // p-v zero block
  FEMatrix pressure_pressure(pressure, pressure);        // C block

  // fill in the velo-velo blocks
  my_matrix.replace_blocks(velo_velo, {{0,0}, {1,1}}, {false, false});
  my_matrix.replace_blocks(velo_velo_zero, {{1,0}, {0,1}}, {false, false});

  // fill in the pressure_velo_zero blocks
  my_matrix.replace_blocks(pressure_velo_zero, {{2,0}, {0,2}, {2,1}, {1,2}},
                           {false, true, false, true});

  // fill in the pressure_pressure block
  my_matrix.replace_blocks(pressure_pressure, {{2,2}}, {false});

  return my_matrix;
}
#else
/**
  * Named constructor for a matrix taking the block structure
  *
  * ( A  0  0  0 )
  * ( 0  A  0  0 )
  * ( 0  0  A  0 )
  * ( 0  0  0  C )
  *
  * @param velocity The velocity finite element space.
  * @param pressure The pressure finite element space.
  * @return A newly constructed BlockFEMatrix for NSE3D POD-ROM problems
  */
BlockFEMatrix NSE3D_Diag(std::shared_ptr<const TFESpace3D> velocity,
                             std::shared_ptr<const TFESpace3D> pressure)
{
  BlockFEMatrix my_matrix({velocity, velocity, velocity, pressure});

  //create new blocks with correct structures filled with 0
  FEMatrix velo_velo(velocity, velocity);                // A block
  FEMatrix velo_velo_zero(velocity, velocity, true);     // velocity zero block
  FEMatrix pressure_velo_zero(pressure, velocity, true); // p-v zero block
  FEMatrix pressure_pressure(pressure, pressure);        // C block

  // fill in the velo-velo blocks
  my_matrix.replace_blocks(velo_velo,
                           {{0,0}, {1,1}, {2,2}},
                           {false, false, false});
  my_matrix.replace_blocks(velo_velo_zero,
                           {{1,0}, {0,1}, {2,0}, {0,2}, {2,1}, {1,2}},
                           {false, false, false, false, false, false});

  // fill in the pressure_velo_zero blocks
  my_matrix.replace_blocks(pressure_velo_zero,
                           {{3,0}, {0,3}, {3,1}, {1,3}, {3,2}, {2,3}},
                           {false, true, false, true, false, true});

  // fill in the pressure_pressure block
  my_matrix.replace_blocks(pressure_pressure, {{3,3}}, {false});

  return my_matrix;
}
#endif



#ifdef __3D__
template class TNSE_POD<3>;
#else
template class TNSE_POD<2>;
#endif
