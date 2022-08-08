/** ****************************************************************************
*
* @name      TCD_ROM
* @brief     solve ROM for TCD (only time-independent BC)
*
*******************************************************************************/

#include "ConvDiff.h"
#include "Database.h"
#include "LinAlg.h"
#include "MainUtilities.h"
#include "TCD_ROM.h"

#ifdef __2D__
#include "Assemble2D.h"
#include "AuxParam2D.h"
#include "SquareMatrix2D.h"
#else
#include "Assemble3D.h"
#include "AuxParam3D.h"
#include "SquareMatrix3D.h"
#endif

#include <algorithm>

/** ***************************************************************************/
template <int d>
ParameterDatabase TCD_ROM<d>::set_databaseROM(const ParameterDatabase& param_db)
{
  ParameterDatabase db("Parameter database for TCD_ROM");

  db.add("rom_rank", (size_t) 0,
         "This integer specifies the dimension of POD-basis to be used for "
         "the ROM space. "
         "If = 0, then all possible POD modes will be used.");

  db.merge(LocalAssembling<d>::default_local_assembling_database(), true);

  db.merge(ParameterDatabase::default_output_database(), true);
  db.merge(TimeDiscretization::default_TimeDiscretization_database(), true);
  db.merge(ROM::default_rom_database(), true);
  db.merge(param_db, false);

  // change name to avoid overwritting standard FEM output or POD-basis
  // TODO: unification between naming system (outfile, output_basename,
  //       snaps_basename, pod_basename)
  std::string output_rom_basename = db["output_basename"].get<std::string>()
                                    + "_rom";
  db["output_basename"].set(output_rom_basename, false);

  return db;
}

/** ***************************************************************************/
template<int d>
TCD_ROM<d>::TCD_ROM(const TDomain&           domain,
                    const ParameterDatabase& param_db)
  : TCD_ROM<d>(domain, param_db, Example_TimeCD(param_db))
{
}

/** ***************************************************************************/
template<int d>
TCD_ROM<d>::TCD_ROM(const TDomain&           domain,
                    const ParameterDatabase& param_db,
                    const Example_TimeCD&    ex)
  : TimeConvectionDiffusion<d>(domain,
                               TCD_ROM<d>::set_databaseROM(param_db), // output
                               ex),
    ROM(param_db,
       TCD_ROM<d>::set_databaseROM(param_db)["rom_rank"].Parameter::get<size_t>()),
    db(TCD_ROM<d>::set_databaseROM(param_db)),
    mat_time_dependent(db["mat_time_dependent"]),
    sys_mat_r(POD::rank, POD::rank),
    sys_rhs_r(POD::rank),
    sol_r(POD::rank),
    sol_m1_r(POD::rank),
    mass_mat_r(POD::rank, POD::rank),
    cdr_mat_r(POD::rank, POD::rank)
{
  check_and_set_parameters();

  if(db["time_discretization"].is("bdf_two"))
  {
    sol_m2_r.resize(POD::rank);
    bdf2_coeffs = this->get_time_stepping_scheme().get_bdf_coefficients();
  }

  bool usingMultigrid = this->solver.is_using_multigrid();
  auto collections = domain.get_grid_collections();
  int ansatz_order = TDatabase::ParamDB->ANSATZ_ORDER;

  if(!usingMultigrid)
  {
    // the given collection for particular cell
    TCollection& cellCollection = *collections.front();
    this->systems.emplace_back(this->example, cellCollection, ansatz_order);

    // initial condition on the solution
    this->systems.front().fe_function.Interpolate(this->example.
                                                  get_initial_cond(0));
  }

  // print useful information
  this->output_problem_size_info();
  this->outputWriter.add_fe_function(&this->systems.front().fe_function);
//    // initialize L_inf error to some negative number
//   errors[TCD_ROM::n_errors-1] = -1.;
}

/** ***************************************************************************/
template<int d>
void TCD_ROM<d>::check_and_set_parameters()
{
  if(this->solver.is_using_multigrid())
  {
    ErrThrow("Multigrid for ROM is not yet enabled in TimeConvectionDiffusion");
  }

  if(db["space_discretization_type"].is("supg"))
  {
    if(mat_time_dependent != true)
      Output::print<1>("WARNING: SUPG-ROM computation with time-independent "
                       "matrices, be sure that parameters for convection, "
                       "diffusion and reaction are time-independent.");
  }

  if(db["time_discretization"].is("bdf_two"))
  {
    pre_stage_bdf = true;
  }

#ifdef _MPI
  ErrThrow("MPI is not yet enabled for TimeConvectionDiffusion with ROM");
#endif

}

/** ***************************************************************************/
template<int d>
void TCD_ROM<d>::compute_initial_condition()
{
  // assemble and set gramian matrix
  if(db["pod_inner_product"].is("L2"))
  {
    assemble_and_reduce(LocalAssembling_type::TCDMassOnly);
    POD::set_gramian_ptr(this->systems.front().
                         mass_matrix.get_combined_matrix());
  }

  if(db["rom_init_regularized"])
  {
    /* see S.Giere, PhD Thesis 2016 */
    Output::print<1>("Type of ROM initial condition: regularized");

    // reduce initial solution using gramian matrix
    ROM::reduce_init_solution(this->systems.front().solution.get_entries(),
                              sol_r);

    // assemble mass and stiffness matrices and Rhs of the second order
    // Helmholtz equation: -mu*mu*delta(u_) + u_ = u
    assemble_and_reduce(LocalAssembling_type::Custom);

    double mu = db["differential_filter_width"];

    // set system matrix
    sys_mat_r = mass_mat_r;
    sys_mat_r.add(cdr_mat_r, mu*mu);

    // set system rhs
    sys_rhs_r = mass_mat_r.multiply(&sol_r);

    if(POD::db["pod_fluctuations_only"])
    {
      Daxpy(POD::rank, -mu * mu, &cdr_mat_mean_r[0], &sys_rhs_r[0]);
    }

    // solve Helmholtz equation
    solve();
    // clear LU values corresponding to the second order Helmholtz equation
    sys_mat_r.clearLU();
  }
  else
  {
    // reduce initial solution using gramian matrix
    ROM::reduce_init_solution(this->systems.front().solution.get_entries(),
                              sol_r);
    Output::print<1>("Reducing initial condition done.");
  }

  // reset gramian matrix which is no longer needed
  if(!db["pod_inner_product"].is("euclidean"))
  {
    POD::set_gramian_ptr(nullptr);
  }

  sol_m1_r = sol_r;
}

/** ***************************************************************************/
template<int d>
void TCD_ROM<d>::assemble_initial_time()
{
  if(mat_time_dependent)
  {
    Output::print<1>("ROM computation will be proceeded with time dependent "
                     "stiffness and mass matrices.");
  }
  else
  {
    Output::print<1>("ROM computation will be proceeded with time independent "
                     "stiffness and mass matrices.");
  }

  compute_initial_condition();

  // assemble Mass and Stiffness matrices and Rhs
  assemble_and_reduce(LocalAssembling_type::TCDStiffMassRhs);

  // prepare system_matrix (only needed in case of mat_time_dependent = false)
  if(!mat_time_dependent)
  {
    prepare_system_matrix_rhs(LocalAssembling_type::Custom);
  }
}

/** ***************************************************************************/
template<int d>
void TCD_ROM<d>::assemble()
{
  LocalAssembling_type type;

  // in case of BDF2 time scheme, the matrices have to be reassemble at the
  // second time step to take into account the actual time discretisation (BDF2)
  if(mat_time_dependent || pre_stage_bdf)
  {
    if(db["space_discretization_type"].is("supg"))
    {
      // assemble the mass matrix
      type = LocalAssembling_type::TCDStiffMassRhs;
    }
    else
    {
    // assemble the stiffness matrix and right hand side
    type = LocalAssembling_type::TCDStiffRhs;
    }
  }
  else
  {
    // assemble the right hand side only
    type = LocalAssembling_type::TCDRhsOnly;
  }

  // in case of BDF2 time scheme, check if the first time step is already done
  // and update before assembling (but after deciding what has to be assembled)
  if(pre_stage_bdf == true)
  {
    pre_stage_bdf = this->get_time_stepping_scheme().pre_stage_bdf;
  }

  prepare_system_matrix_rhs(type);
}

/** ***************************************************************************/
template<int d>
void TCD_ROM<d>::prepare_system_matrix_rhs(LocalAssembling_type type)
{
  double tau = this->time_stepping_scheme.get_step_length();

  // ===========================================================================
  // Set system_rhs
  // ===========================================================================
  if(  type == LocalAssembling_type::TCDRhsOnly
    || type == LocalAssembling_type::TCDStiffRhs
    || type == LocalAssembling_type::TCDStiffMassRhs)
  {
    rhs_m1_r = std::move(rhs_r);

    // -------------------------------------------------------------------------
    // The next line updates the source terms (rhs and if needed stiffness and
    // mass and matrices, depending on 'type')
    // NOTE: Online (= within the time loop) assembling involving FE dimension
    // is not efficient in the ROM context. To avoid it, one could alternatively
    // store all FE coefficients of the source term and reduce them offline, and
    // use here the corresponding reduced-order vectors. For source terms that
    // are formulated in the separated time-space form one can pre-assemble and
    // reduce the space part before the time loop and multiply it in every time
    // step with the appropriate temporal factor.
    // -------------------------------------------------------------------------
    assemble_and_reduce(type);

    // -------------------------------------------------------------------------
    // add the old/new terms according to the time discretization scheme
    // -------------------------------------------------------------------------
    if(db["time_discretization"].is("backward_euler")
      || (pre_stage_bdf && db["bdf_two_first_step"].is("backward_euler")))
    {
      sys_rhs_r = mass_mat_r.multiply(&sol_m1_r);
      Daxpy(POD::rank, tau, &rhs_r[0], &sys_rhs_r[0]);
      if(db["pod_fluctuations_only"])
      {
        Daxpy(POD::rank, -tau, &cdr_mat_mean_r[0], &sys_rhs_r[0]);
      }
      if(db["time_discretization"].is("bdf_two") && pre_stage_bdf)
      {
        Output::print<3>("First step in BDF2 scheme is performed by the "
                         "backward-Euler scheme");
      }
    }
    else if(db["time_discretization"].is("crank_nicolson")
      || (pre_stage_bdf && db["bdf_two_first_step"].is("crank_nicolson")))
    {
      sys_rhs_r = mass_mat_r.multiply(&sol_m1_r);
      Daxpy(POD::rank, -tau * 0.5, &(cdr_mat_r.multiply(&sol_m1_r))[0],
            &sys_rhs_r[0]);
      Daxpy(POD::rank, tau * 0.5, &rhs_m1_r[0], &sys_rhs_r[0]);
      Daxpy(POD::rank, tau * 0.5, &rhs_r[0], &sys_rhs_r[0]);
      if(db["pod_fluctuations_only"])
      {
        Daxpy(POD::rank, -tau, &cdr_mat_mean_r[0], &sys_rhs_r[0]);
      }
      if(db["time_discretization"].is("bdf_two") && pre_stage_bdf)
      {
        Output::print<3>("First step in BDF2 scheme is performed by the "
                         "Crank-Nicolson scheme");
      }
    }
    else if(db["time_discretization"].is("bdf_two") && !pre_stage_bdf)
    {
      sys_rhs_r = mass_mat_r.multiply(&sol_m1_r, bdf2_coeffs[0]);
      Daxpy(POD::rank, bdf2_coeffs[1], &(mass_mat_r.multiply(&sol_m2_r))[0],
            &sys_rhs_r[0]);
      Daxpy(POD::rank, tau * bdf2_coeffs[2], &rhs_r[0], &sys_rhs_r[0]);
      if(db["pod_fluctuations_only"])
      {
        Daxpy(POD::rank, -tau * bdf2_coeffs[2], &cdr_mat_mean_r[0],
              &sys_rhs_r[0]);
      }
    }
    else
    {
      ErrThrow("time discretization: ", db["time_discretization"],
              " is not supported");
    }
  }

  // ===========================================================================
  // Set system_matrix
  // ===========================================================================
  if(  type == LocalAssembling_type::TCDStiffMassRhs
    || type == LocalAssembling_type::TCDStiffRhs
    || type == LocalAssembling_type::Custom)
  {
    sys_mat_r = mass_mat_r;

    if(db["time_discretization"].is("backward_euler")
      || (pre_stage_bdf && db["bdf_two_first_step"].is("backward_euler")))
    {
      sys_mat_r.add(cdr_mat_r, tau);
    }
    else if(db["time_discretization"].is("crank_nicolson")
      || (pre_stage_bdf && db["bdf_two_first_step"].is("crank_nicolson")))
    {
      sys_mat_r.add(cdr_mat_r, tau * 0.5);
    }
    else if(db["time_discretization"].is("bdf_two") && !pre_stage_bdf)
    {
      sys_mat_r.add(cdr_mat_r, tau * bdf2_coeffs[2]);
    }
    else
    {
      ErrThrow("time discretization: ", db["time_discretization"],
               " is not supported");
    }
  }
}

/** ***************************************************************************/
template<int d>
void TCD_ROM<d>::solve()
{
  if(sys_mat_r.getNRows() != (int)sys_rhs_r.size())
  {
    ErrThrow("Dimension mismatch during A*x = rhs system resolution: ",
             sys_mat_r.getNRows(), " rows in the system matrix and ",
             sys_rhs_r.size(), " elements in the rhs.");
  }

  if( !sys_mat_r.getEntriesLU() || !sys_mat_r.getPivotsLU()
    || mat_time_dependent )
  {
    sys_mat_r.decomposeLU();
  }

  sys_mat_r.solve(&sys_rhs_r[0]); // solution is stored in sys_rhs_r
  sol_r = std::move(sys_rhs_r);

  // update old solutions for the next time step
  if(db["time_discretization"].is("bdf_two"))
  {
    sol_m2_r = std::move(sol_m1_r);
  }
  sol_m1_r = sol_r;
}

/** ***************************************************************************/
template<int d>
void Helmholtz_coeffs(int n_points, const double*, const double*,
#ifdef __3D__
                                                      const double*,
#endif
                      const double*const*, double** coeffs)
{
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = 1.; // diffusion coefficient epsilon
    coeffs[i][1] = 0.; // convection coefficient b1
    coeffs[i][2] = 0.; // convection coefficient b2
    if(d==3)
    {
      coeffs[i][d] = 0.; // convection coefficient b3
    }
    coeffs[i][d+1] = 0.; // reaction coefficient c
    coeffs[i][d+2] = 0.;
  }
}

/** ***************************************************************************/
template<int d>
void TCD_ROM<d>::assemble_and_reduce(LocalAssembling_type type)
{

  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  std::vector<std::shared_ptr<FEMatrix>> block;

  int n_fespaces = 1;
  int nRhs = 0;
  const int nsqMatrices = 2; // maximum number of matrices

  auto& system = this->systems.front();
  const FESpace* fe_space = system.fe_space.get();
  double *rhsEntries = nullptr;
  SquareMatrixD *sqMatrices[nsqMatrices] = {nullptr, nullptr};

  // ===========================================================================
  // Assemble
  // ===========================================================================
  if(  type == LocalAssembling_type::TCDStiffMassRhs
    || type == LocalAssembling_type::Custom
    || type == LocalAssembling_type::TCDMassOnly)
  {
    sqMatrices[1] = reinterpret_cast<SquareMatrixD*>(system.mass_matrix.
                                                get_blocks_uniquely()[0].get());
    sqMatrices[1]->reset();
  }
  if(  type == LocalAssembling_type::TCDStiffRhs
    || type == LocalAssembling_type::TCDStiffMassRhs
    || type == LocalAssembling_type::Custom)
  {
    sqMatrices[0] = reinterpret_cast<SquareMatrixD*>(system.stiffness_matrix.
                                                get_blocks_uniquely()[0].get());
    sqMatrices[0]->reset();
  }
  if(  type == LocalAssembling_type::TCDRhsOnly
    || type == LocalAssembling_type::TCDStiffRhs
    || type == LocalAssembling_type::TCDStiffMassRhs
    || type == LocalAssembling_type::Custom)
  {
    nRhs = 1;
    rhsEntries = system.rhs.get_entries();
    system.rhs.reset();
  }
  if(  type != LocalAssembling_type::TCDRhsOnly
    && type != LocalAssembling_type::TCDStiffRhs
    && type != LocalAssembling_type::TCDStiffMassRhs
    && type != LocalAssembling_type::Custom
    && type != LocalAssembling_type::TCDMassOnly)
  {
    ErrThrow("LocalAssembling_type ", type, " not supported");
  }

  std::vector<const FEFunction*> pointer_to_function = {&system.fe_function};
  ParameterDatabase db_assembling(db);

  // flag for standard assembling (true) or special assembling needed for the
  // computation of the initial solution
  bool standard = (type != LocalAssembling_type::Custom) ? true : false;

  // in case of special assembling, insure that no SUPG terms will be added
  if(!standard)
  {
    type = LocalAssembling_type::TCDStiffMassRhs;
    db_assembling["space_discretization_type"].set("galerkin", false);
  }

  // local assembling
  LocalAssembling<d> la(db_assembling, type, pointer_to_function,
                        standard ? this->example.get_coeffs()
                                 : Helmholtz_coeffs<d>);

  // boundary conditions and boundary values
  auto* boundary_conditions = fe_space->get_boundary_condition();
  auto* boundary_value      = this->example.get_bd(0);

#ifdef __2D__
  Assemble2D(
#else
  Assemble3D(
#endif
             n_fespaces, &fe_space, nsqMatrices, sqMatrices,
             0, nullptr, nRhs, &rhsEntries,
             &fe_space, &boundary_conditions, &boundary_value, la);

  // ===========================================================================
  // Reduce
  // ===========================================================================
  if(type == LocalAssembling_type::TCDStiffMassRhs)
  {
    Output::print<5>("Reducing finite element mass matrix...");
    ROM::reduce(system.mass_matrix.get_combined_matrix(), mass_mat_r);
  }

  if(  type == LocalAssembling_type::TCDStiffRhs
    || type == LocalAssembling_type::TCDStiffMassRhs)
  {
    Output::print<5>("Reducing finite element stiffness matrix...");
    ROM::reduce(system.stiffness_matrix.get_combined_matrix(), cdr_mat_r);

    if(db["pod_fluctuations_only"])
    {
      ROM::reduce_mat_mean(system.stiffness_matrix.get_combined_matrix(),
                           cdr_mat_mean_r);
    }
  }
  if(  type == LocalAssembling_type::TCDRhsOnly
    || type == LocalAssembling_type::TCDStiffRhs
    || type == LocalAssembling_type::TCDStiffMassRhs)
  {
    Output::print<5>("Reducing finite element source term...");
    ROM::reduce(system.rhs.get_vector_entries(), rhs_r);
  }
}

/** ***************************************************************************/
template<int d>
void TCD_ROM<d>::output()
{
  ROM::get_full_solution(this->sol_r,
                         this->systems.front().solution.get_entries());

  TimeConvectionDiffusion<d>::output();
}

#ifdef __3D__
template class TCD_ROM<3>;
#else
template class TCD_ROM<2>;
#endif
