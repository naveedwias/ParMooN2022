/** ****************************************************************************
*
* @name      TNSE_ROM
* @brief     solve ROM for TNSE
*
*******************************************************************************/

#include "CD_local_assembling_routines.h" // for P-ROM gramian
#include "ConvDiff.h"
#include "Database.h"
#include "LinAlg.h"
#include "MainUtilities.h"
#include "ROM_local_assembling_routines.h"
#include "TNSE_ROM.h"

#ifdef __2D__
#include "Assemble2D.h"
#include "AuxParam2D.h"
#include "Matrix2D.h"
#include "SquareMatrix2D.h"
#else
#include "Assemble3D.h"
#include "AuxParam3D.h"
#include "Matrix3D.h"
#include "SquareMatrix3D.h"
#endif

#include <algorithm>

/** ***************************************************************************/
template <int d>
ParameterDatabase TNSE_ROM<d>::set_databaseROM(
                                              const ParameterDatabase& param_db,
                                              bool is_full)
{
  ParameterDatabase db("Parameter database for TNSE_ROM");
  db.add("vp_rom_rank", std::vector<size_t>(2, 0),
         "These integers specify the dimensions of POD-basis to be used for "
         "the u-ROM and p-ROM spaces, respectively. "
         "If = 0, then all possible POD modes will be used.");

//   db.add("print_fem_residuals", false,
//          "This is the flag indicating if the corresponding FEM residuals have "
//          " to be computed and printed out.",
//          {true,false});

  db.merge(LocalAssembling<d>::default_local_assembling_database(), true);

  db.merge(ParameterDatabase::default_output_database(), true);
  db.merge(TimeDiscretization::default_TimeDiscretization_database(), true);
  db.merge(ROM::default_rom_database(), true);

  if (is_full)
  {
    db.merge(param_db, true);
  }
  else
  {
    db.merge(param_db, false);
  }

  // change name to avoid overwritting standard FEM output or POD-basis
  // when passed to constructor TimeNavierStokes()
  // TODO: unification between naming system (outfile, output_basename,
  //       snaps_basename, pod_basename)
  std::string output_rom_basename = db["output_basename"].get<std::string>()
                                    + "_rom";
  db["output_basename"].set(output_rom_basename, false);

  return db;
}

/** ***************************************************************************/
template<int d>
TNSE_ROM<d>::TNSE_ROM(const TDomain&           domain,
                      const ParameterDatabase& param_db)
  : TNSE_ROM<d>(domain, param_db, Example_TimeNSE(param_db))
{
}

/** ***************************************************************************/
template<int d>
TNSE_ROM<d>::TNSE_ROM(const TDomain&           domain,
                      const ParameterDatabase& param_db,
                      const Example_TimeNSE&   ex)
  : TimeNavierStokes<d>(domain, TNSE_ROM<d>::set_databaseROM(param_db,true),ex),
    ROM_u(param_db,
          TNSE_ROM<d>::set_databaseROM(param_db)["vp_rom_rank"].
            Parameter::get<std::vector<size_t>>()[0], "u"),
    ROM_p(param_db,
          TNSE_ROM<d>::set_databaseROM(param_db)["vp_rom_rank"].
            Parameter::get<std::vector<size_t>>()[1], "p"),
    db(TNSE_ROM<d>::set_databaseROM(param_db)),
    rank_u(ROM_u.get_rank()),
    sys_mat_u_r(rank_u, rank_u),
    lhs_mat_u_r(rank_u, rank_u),
    sys_rhs_u_r(rank_u),
    rhs_mat_u_r(rank_u, rank_u),
    mass_mat_u_r(rank_u, rank_u),
    diff_mat_u_r(rank_u, rank_u),
    diff_mean_u_r(rank_u),
    conv_lin_mat_u_r(rank_u, rank_u),
    conv_lin_mean_u_r(rank_u),
    conv_non_lin_mat_u_r(rank_u, conv_lin_mat_u_r),
    conv_non_lin_mean_u_r(rank_u, sys_rhs_u_r),
    rhs_src_u_r(rank_u),
    sol_u_r(rank_u),
    sol_m1_u_r(rank_u),
    rank_p(ROM_p.get_rank()),
    sys_mat_p_r(rank_p, rank_p),
    sys_rhs_p_r(rank_p),
    diff_mat_p_r(rank_p, rank_p),
    div_mat_p_r(rank_p, rank_u),
    div_mean_p_r(rank_p),
    conv_lin_mat_p_r(rank_p, rank_u),
    conv_lin_mean_p_r(rank_p),
    conv_non_lin_mat_p_r(rank_u, conv_lin_mat_p_r),
    conv_non_lin_mean_p_r(rank_u, sys_rhs_p_r),
    rhs_src_p_r(rank_p),
    sol_p_r(rank_p)
{
  check_and_set_parameters();

  if (db["time_discretization"].is("crank_nicolson"))
  {
    rhs_src_m1_u_r.resize(rank_u);
  }

  if (db["time_discretization"].is("bdf_two"))
  {
    sol_m2_u_r.resize(rank_u);
    bdf2_coeffs = this->get_time_stepping_scheme().get_bdf_coefficients();
  }
}

/** ***************************************************************************/
template<int d>
void TNSE_ROM<d>::check_and_set_parameters()
{
  if (db["pod_fluctuations_only"].is(false))
  {
    ErrThrow("The parameter 'pod_fluctuations_only' has to be activated "
             "otherwise only problems with homogeneous boundary conditions "
             "can be solved.");
  }

  if (this->solver.is_using_multigrid())
  {
    ErrThrow("Multigrid for ROM is not yet enabled in TimeNavierStokes");
  }

  if (!db["space_discretization_type"].is("galerkin"))
  {
    ErrThrow("The space discretization ", db["space_discretization_type"],
             " is not yet enabled in TimeNavierStokes for ROM computation.");
  }

  if (this->get_db()["time_discretization_nonlinear_term"].is("fully_explicit"))
  {
    ErrThrow("The fully explicit time discretization for the nonlinear term, "
             "is not yet enabled in TimeNavierStokes for ROM computation.");
  }

  // if (db["nse_nonlinear_form"].is("divergence"))
  if (!db["nse_nonlinear_form"].is("convective"))
  {
    ErrThrow("The nonlinear form ", db["nse_nonlinear_form"],
             " is not yet enabled in TimeNavierStokes for ROM computation.");
  }

  if (!db["graddiv_stab"].is(0.))
  {
    Output::warn("Stabilization parameters are not yet supported for ROM-TNSE. "
                 "The parameter 'graddiv_stab' is reset to 0 ");
    db["graddiv_stab"] = 0.;
  }

  if (db["mat_time_dependent"].is(true))
  {
    Output::warn("The parameter 'mat_time_dependent' is not yet implemented "
                 "for ROM-TNSE");
  }

  if (db["rom_init_regularized"].is(true))
  {
    Output::warn("The parameter 'rom_init_regularized' is not yet implemented "
                 "for ROM-TNSE");
  }

#ifdef _MPI
  ErrThrow("MPI is not yet enabled for TimeNavierStokes with ROM");
#endif

//   if (db["print_fem_residuals"])
//   {
//     fem_residuals = true;
//   }

  if (db["time_discretization"].is("bdf_two"))
  {
    pre_stage_bdf = true;
  }

}

/** ***************************************************************************/
template<int d>
void TNSE_ROM<d>::compute_initial_condition()
{
  // assemble and set gramian matrices
  if (!db["pod_inner_product"].is("euclidean"))
  {
    ROM_u.set_gramian_ptr(assemble_gramian());

    auto pressure_space = this->get_pressure_space();
#ifdef __3D__
    BlockFEMatrix pressure_matrix = BlockFEMatrix::CD3D(pressure_space);
#else
    BlockFEMatrix pressure_matrix = BlockFEMatrix::CD2D(pressure_space);
#endif
    ROM_p.set_gramian_ptr(assemble_gramian(&pressure_matrix));
  }

  auto& s = this->systems.front();
  ROM_u.reduce_init_solution(s.solution.block(0), sol_u_r);
  ROM_p.reduce_init_solution(s.solution.block(d), sol_p_r);
  Output::print<1>("Reducing initial condition done.");

  // reset gramian matrix which is no longer needed
  if (!db["pod_inner_product"].is("euclidean"))
  {
    ROM_u.set_gramian_ptr(nullptr);
    ROM_p.set_gramian_ptr(nullptr);
  }

  sol_m1_u_r = sol_u_r;
}

/** ***************************************************************************/
template<int d>
void TNSE_ROM<d>::assemble_initial_time()
{
  compute_initial_condition();

  // assemble Mass, Diff, Conv matrices and Rhs source
  bool is_velocity = true;
  assemble_system(LocalAssembling_type::TimeNavierStokesAll, is_velocity);
  assemble_system(LocalAssembling_type::TimeNavierStokesAll, !is_velocity);

// for FEM residual computation
//   if (fem_residuals)
//     TimeNavierStokes<d>::assemble_initial_time();
}

/** ***************************************************************************/
template<int d>
void TNSE_ROM<d>::assemble_matrices_rhs(unsigned int it_counter)
{
  bool first_nonlinear_step = (it_counter == 0);

  // save old values and assemble Rhs only at each first nonlinear step
  if (first_nonlinear_step)
  {
    // update old rhs and solutions for the current time step
    if (db["time_discretization"].is("crank_nicolson"))
    {
      rhs_src_m1_u_r = rhs_src_u_r;
    }
    if (db["time_discretization"].is("bdf_two"))
    {
      sol_m2_u_r = std::move(sol_m1_u_r);
    }
    sol_m1_u_r = sol_u_r;

    bool is_velocity = true;
    assemble_system(LocalAssembling_type::TimeNavierStokesRhs, is_velocity);
    assemble_system(LocalAssembling_type::TimeNavierStokesRhs, !is_velocity);
  }

  // in case of BDF2 time scheme, check if the first time step is already done
  if (pre_stage_bdf == true)
  {
    pre_stage_bdf = this->get_time_stepping_scheme().pre_stage_bdf;
  }

  prepare_system_velocity(first_nonlinear_step);

// for FEM residual computation
//   if (fem_residuals)
//     TimeNavierStokes<d>::assemble_matrices_rhs(it_counter);
}

/** ***************************************************************************/
template<int d>
void TNSE_ROM<d>::prepare_system_velocity(bool first_nonlinear_step)
{
  double tau = this->time_stepping_scheme.get_step_length();

  // compute linear parts (to be updated at each first nonlinear step
  // if tau is not constant)
  if (first_nonlinear_step && (tau_old != tau))
  {
    prepare_offline_velocity();
    tau_old = tau;
  }

  // compute nonlinear parts of the convective term
  DenseMatrix aux_conv_mat(rank_u, rank_u);
  std::vector<double>  aux_conv_mean_vect(rank_u);
  for (unsigned int r = 0; r < (unsigned int)rank_u; r++)
  {
    aux_conv_mat.add(conv_non_lin_mat_u_r[r], sol_u_r[r]);
    if (db["pod_fluctuations_only"])
    {
      Daxpy(rank_u, sol_u_r[r], &(conv_non_lin_mean_u_r[r])[0],
            &aux_conv_mean_vect[0]);
    }
  }

  // compute system matrix and system rhs
  sys_mat_u_r = lhs_mat_u_r;
  if (db["time_discretization"].is("backward_euler")
    || (pre_stage_bdf && db["bdf_two_first_step"].is("backward_euler")))
  {
    sys_mat_u_r.add(aux_conv_mat, tau);
    sys_rhs_u_r = rhs_mat_u_r.multiply(&sol_m1_u_r);
    Daxpy(rank_u, tau, &rhs_src_u_r[0], &sys_rhs_u_r[0]); ///add source

    if (db["pod_fluctuations_only"])
    {
      /// T*A*u_mean
      Daxpy(rank_u, -tau, &diff_mean_u_r[0], &sys_rhs_u_r[0]);
      /// Phi^T*N(u_mean)*u_mean
      Daxpy(rank_u, -tau, &conv_lin_mean_u_r[0], &sys_rhs_u_r[0]); 
      Daxpy(rank_u, -tau, &aux_conv_mean_vect[0], &sys_rhs_u_r[0]);
    }
  }
  else if (db["time_discretization"].is("crank_nicolson")
    || (pre_stage_bdf && db["bdf_two_first_step"].is("crank_nicolson")))
  {
    sys_mat_u_r.add(aux_conv_mat, tau * 0.5);
    sys_rhs_u_r = rhs_mat_u_r.multiply(&sol_m1_u_r);
    Daxpy(rank_u, tau * 0.5, &rhs_src_m1_u_r[0], &sys_rhs_u_r[0]); ///add source
    Daxpy(rank_u, tau * 0.5, &rhs_src_u_r[0], &sys_rhs_u_r[0]); ///add source
    if (db["pod_fluctuations_only"])
    {
      Daxpy(rank_u, -tau, &diff_mean_u_r[0], &sys_rhs_u_r[0]);
      Daxpy(rank_u, -tau, &conv_lin_mean_u_r[0], &sys_rhs_u_r[0]);
      Daxpy(rank_u, -tau, &aux_conv_mean_vect[0], &sys_rhs_u_r[0]);
    }
  }
  else if (db["time_discretization"].is("bdf_two") && !pre_stage_bdf)
  {
    sys_mat_u_r.add(aux_conv_mat, tau * bdf2_coeffs[2]);
    sys_rhs_u_r = rhs_mat_u_r.multiply(&sol_m1_u_r, bdf2_coeffs[0]);
    Daxpy(rank_u, bdf2_coeffs[1],
          &(rhs_mat_u_r.multiply(&sol_m2_u_r))[0], &sys_rhs_u_r[0]);
    if (db["pod_fluctuations_only"])
    {
      Daxpy(rank_u, -tau * bdf2_coeffs[2], &diff_mean_u_r[0], &sys_rhs_u_r[0]);
      Daxpy(rank_u, -tau * bdf2_coeffs[2], &conv_lin_mean_u_r[0],
                                           &sys_rhs_u_r[0]);
      Daxpy(rank_u, -tau * bdf2_coeffs[2], &aux_conv_mean_vect[0],
                                           &sys_rhs_u_r[0]);
    }
  }
}

/** ***************************************************************************/
template<int d>
void TNSE_ROM<d>::prepare_offline_velocity()
{
  double tau = this->time_stepping_scheme.get_step_length();

  lhs_mat_u_r = mass_mat_u_r;
  rhs_mat_u_r = mass_mat_u_r;

  if (db["time_discretization"].is("backward_euler")
    || (pre_stage_bdf && db["bdf_two_first_step"].is("backward_euler")))
  {
    lhs_mat_u_r.add(diff_mat_u_r, tau);
    lhs_mat_u_r.add(conv_lin_mat_u_r, tau);
  }
  else if (db["time_discretization"].is("crank_nicolson")
    || (pre_stage_bdf && db["bdf_two_first_step"].is("crank_nicolson")))
  {
    lhs_mat_u_r.add(diff_mat_u_r, tau * 0.5);
    lhs_mat_u_r.add(conv_lin_mat_u_r, tau * 0.5);
    rhs_mat_u_r.add(diff_mat_u_r, -tau * 0.5);
    rhs_mat_u_r.add(conv_lin_mat_u_r, -tau * 0.5);
  }
  else if (db["time_discretization"].is("bdf_two") && !pre_stage_bdf)
  {
    lhs_mat_u_r.add(diff_mat_u_r, tau * bdf2_coeffs[2]);
    lhs_mat_u_r.add(conv_lin_mat_u_r, tau * bdf2_coeffs[2]);
  }
  else
  {
    ErrThrow("time discretization: ", db["time_discretization"],
              " is not supported");
  }
}

/** ***************************************************************************/
template<int d>
void TNSE_ROM<d>::prepare_system_pressure(void)
{
  // only assemble once (sys_mat_p_r is time independent) for this purpose
  // check if LU decomposition was not performed yet:
  if (!sys_mat_p_r.getEntriesLU())
  {
    sys_mat_p_r = diff_mat_p_r;
  }

  double tau = this->time_stepping_scheme.get_step_length();
  std::vector<double> u_dt(sol_u_r.size());
  // u_dt = sol_u_r - sol_m1_u_r:
  std::transform(sol_u_r.begin(), sol_u_r.end(), sol_m1_u_r.begin(),
                 u_dt.begin(), std::minus<double>());

  // compute nonlinear parts of the convective term
  DenseMatrix aux_conv_mat(rank_p, rank_u);
  std::vector<double>  aux_conv_mean_vect(rank_p);
  for (unsigned int r = 0; r < (unsigned int)rank_u; r++)
  {
    aux_conv_mat.add(conv_non_lin_mat_p_r[r], sol_u_r[r]);
    if (db["pod_fluctuations_only"])
    {
      Daxpy(rank_p, sol_u_r[r], &(conv_non_lin_mean_p_r[r])[0],
            &aux_conv_mean_vect[0]);
    }
  }

  sys_rhs_p_r = div_mat_p_r.multiply(&u_dt, 1./tau);
  Daxpy(rank_p, 1., &(conv_lin_mat_p_r.multiply(&sol_u_r))[0], &sys_rhs_p_r[0]);
  Daxpy(rank_p, 1., &(aux_conv_mat.multiply(&sol_u_r))[0], &sys_rhs_p_r[0]);
  Daxpy(rank_p, 1., &rhs_src_p_r[0], &sys_rhs_p_r[0]);
  if (db["pod_fluctuations_only"])
  {
    Daxpy(rank_p, 1., &conv_lin_mean_p_r[0], &sys_rhs_p_r[0]);
    Daxpy(rank_p, 1., &aux_conv_mean_vect[0], &sys_rhs_p_r[0]);
  }
}

/** ***************************************************************************/
template<int d>
std::shared_ptr<TMatrix> TNSE_ROM<d>::assemble_gramian(BlockFEMatrix* mat_ptr)
{
  bool is_velocity = (mat_ptr==nullptr) ? true : false;
  bool is_reduced = false;
  LocalAssembling_type type;

  bool is_L2 = false;
  if (db["pod_inner_product"].is("L2"))
  {
    is_L2 = true;
    type = is_velocity ? LocalAssembling_type::TimeNavierStokesMass
                       : LocalAssembling_type::TCDMassOnly;
  }
  else if (db["pod_inner_product"].is("H1"))
  {
    type = is_velocity ? LocalAssembling_type::NavierStokesLinear
                       : LocalAssembling_type::TCDDiffOnly;
  }
  else
  {
    ErrThrow("pod_inner_product: ",db["pod_inner_product"]," not supported");
  }

  assemble_and_reduce(is_reduced, type, is_velocity, mat_ptr);

  if (is_velocity)
  {
    auto& s = this->systems.front();
    size_t n_row = s.mass_matrix.get_n_cell_rows();
    size_t n_col = s.mass_matrix.get_n_cell_columns();

    return (is_L2) ? s.mass_matrix.get_combined_submatrix({0,0},
                                                          {n_row-2,n_col-2})
                   : s.matrix.get_combined_submatrix({0,0},
                                                     {n_row-2,n_col-2});
  }
  else
  {
    return mat_ptr->get_combined_matrix();
  }
}

/** ***************************************************************************/
template<int d>
void TNSE_ROM<d>::assemble_system(LocalAssembling_type type, bool is_velocity)
{
  bool is_reduced = true;

  if (type == LocalAssembling_type::TimeNavierStokesAll)
  {
    if (is_velocity)
    {
      assemble_and_reduce(is_reduced,LocalAssembling_type::TimeNavierStokesMass,
                          is_velocity);
      assemble_and_reduce(is_reduced, LocalAssembling_type::NavierStokesLinear,
                          is_velocity);
      assemble_and_reduce(is_reduced, LocalAssembling_type::TimeNavierStokesNL,
                          is_velocity);
      assemble_and_reduce(is_reduced, LocalAssembling_type::TimeNavierStokesRhs,
                          is_velocity);
    }
    else
    {
      auto pressure_space = this->get_pressure_space();
#ifdef __3D__
      BlockFEMatrix pressure_matrix = BlockFEMatrix::CD3D(pressure_space);
#else
      BlockFEMatrix pressure_matrix = BlockFEMatrix::CD2D(pressure_space);
#endif
      assemble_and_reduce(is_reduced, LocalAssembling_type::Custom,
                          is_velocity, &pressure_matrix);
      assemble_and_reduce(is_reduced, LocalAssembling_type::NavierStokesLinear,
                          is_velocity);
      assemble_and_reduce(is_reduced, LocalAssembling_type::TimeNavierStokesNL,
                          is_velocity);
      assemble_and_reduce(is_reduced, LocalAssembling_type::TimeNavierStokesRhs,
                          is_velocity);
    }
  }
  else if (type == LocalAssembling_type::TimeNavierStokesRhs)
  {
    assemble_and_reduce(is_reduced, type, is_velocity);
  }
}

/** ***************************************************************************/
template<int d>
void TNSE_ROM<d>::assemble_and_reduce(bool is_reduced,
                                      LocalAssembling_type type,
                                      bool is_velocity,
                                      BlockFEMatrix* pressure_matrix)
{
  std::vector<const FESpace*> spaces_mat;
  std::vector<const FESpace*> spaces_rhs;
  std::vector<const FEFunction*> fefunctions;

  std::vector<SquareMatrixD*> sqMatrices;
  std::vector<MatrixD*> rectMatrices;
  std::vector<double*> rhs_array;
  std::vector<BdCondFunct*> boundCondition;
  std::vector<BdValFunct*> boundValues;

  // is_new_term and new_term_idx are usefull for TimeNavierStokesNL:
  // only TimeNavierStokesNL has to be assembled more than one time, it is
  // splitted up into 1 assemble related to the mean velocity function
  // and rank_u assembles related to the velocity POD basis functions
  bool is_new_term = true;
  int new_term_idx = -1;
  while(is_new_term)
  {
    set_arrays(spaces_mat, spaces_rhs, fefunctions, boundCondition, boundValues,
               type, new_term_idx);
    set_matrices_rhs(sqMatrices, rectMatrices, rhs_array,
                     type, is_velocity, pressure_matrix);

    LocalAssembling<d> la(this->db, adapt_type(is_velocity, type), fefunctions,
                        this->example.get_coeffs(),
                        this->space_disc_global);

    adapt_local_assembling(is_velocity, type, la);

#ifdef __3D__
    Assemble3D(
#else
    Assemble2D(
#endif
               spaces_mat.size(), spaces_mat.data(), sqMatrices.size(),
               sqMatrices.data(), rectMatrices.size(), rectMatrices.data(),
               rhs_array.size(), rhs_array.data(), spaces_rhs.data(),
               boundCondition.data(), boundValues.data(), la);

    if (!is_reduced)
    {
      return;
    }

    is_new_term = reduce(type, is_velocity, new_term_idx, pressure_matrix);
  } // while(is_new_term)
}

/** ***************************************************************************/
template<int d>
void TNSE_ROM<d>::set_arrays(std::vector<const FESpace*>&    spaces,
                             std::vector<const FESpace*>&    spaces_rhs,
                             std::vector<const FEFunction*>& functions,
                             std::vector<BdCondFunct*>&      bdCond,
                             std::vector<BdValFunct*>&       bdVal,
                             LocalAssembling_type            type,
                             int                             new_term_idx)
{
  auto& s = this->systems.front();
  auto& velocity_space = s.velocity_space;
  auto& pressure_space = s.pressure_space;

  // for Pressure assembles with TCD structure
  if (type == LocalAssembling_type::TCDMassOnly
    || type == LocalAssembling_type::TCDDiffOnly
    || type == LocalAssembling_type::Custom)
  {
    spaces.resize(1);
    spaces[0] = pressure_space.get();
    spaces_rhs.resize(1);
    spaces_rhs[0] = pressure_space.get();
    functions.resize(1);
    functions[0] = &s.p;
    bdCond.resize(1);
    bdCond[0] = pressure_space->get_boundary_condition();
    bdVal.resize(1);
    bdVal[0] = this->example.get_bd(d);
  }
  // for every assembles (Velocity and Pressure) with TNSE structure:
  // LocalAssembling_type::TimeNavierStokesMass,
  // LocalAssembling_type::NavierStokesLinear,
  // LocalAssembling_type::TimeNavierStokesNL,
  // LocalAssembling_type::TimeNavierStokesRhs
  else
  {
    spaces.resize(2);
    spaces[0] = velocity_space.get();
    spaces[1] = pressure_space.get();
    spaces_rhs.resize(d + 1);
    for (int i = 0; i < d; ++i)
    {
      spaces_rhs[i] = velocity_space.get();
    }
    spaces_rhs[d] = pressure_space.get();
    functions.resize(d + 1);
    FEVectFunct fe_vec_fct;
    if (type == LocalAssembling_type::TimeNavierStokesNL && new_term_idx < 0)
    {
      fe_vec_fct = FEVectFunct(velocity_space,
                               "u_mean",
                               ROM_u.get_snaps_avr_ptr(),
                               d);
    }
    else if (type == LocalAssembling_type::TimeNavierStokesNL && new_term_idx >= 0)
    {
      fe_vec_fct = FEVectFunct(velocity_space,
                               "pod_basis",
                               ROM_u.get_basis_vect(new_term_idx),
                               d);
    }
    else
    {
      fe_vec_fct = this->get_velocity();
    }
    for (int i = 0; i < d; ++i)
    {
      functions[i] = fe_vec_fct.GetComponent(i).release();
    }
    functions[d] = &s.p;
    bdCond.resize(d + 1);
    for (int i = 0; i < d; ++i)
    {
      bdCond[i] = velocity_space->get_boundary_condition();
    }
    bdCond[d] = pressure_space->get_boundary_condition();
    bdVal.resize(d + 1);
    for (int i = 0; i < d + 1; ++i)
    {
      bdVal[i] = this->example.get_bd(i);
    }
  }
}


/** ***************************************************************************/
template<int d>
void TNSE_ROM<d>::set_matrices_rhs(std::vector<SquareMatrixD*>& sqMat,
                                   std::vector<MatrixD*>&       reMat,
                                   std::vector<double*>&        rhs_array,
                                   LocalAssembling_type         type,
                                   bool                         is_velocity,
                                   BlockFEMatrix*               pressure_matrix)
{
  if (is_velocity)
  {
    auto& s = this->systems.front();
    TimeNavierStokes<d>::set_matrices_rhs(s, type, sqMat,
                                          reMat, rhs_array);
  }
  // for Pressure assembles with TCD structure
  else if (type == LocalAssembling_type::TCDMassOnly
    || type == LocalAssembling_type::TCDDiffOnly
    || type == LocalAssembling_type::Custom )
  {
    sqMat.resize(2, nullptr);
    reMat.resize(0, nullptr);
    rhs_array.resize(0, nullptr);
    auto pressureBlock = pressure_matrix->get_blocks_uniquely()[0].get();
    sqMat[1] = reinterpret_cast<SquareMatrixD*>(pressureBlock);
    sqMat[1]->reset();
  }
  // for Pressure assembles with TNSE structure:
  // LocalAssembling_type::NavierStokesLinear,
  // LocalAssembling_type::TimeNavierStokesNL,
  // LocalAssembling_type::TimeNavierStokesRhs
  else
  {
    sqMat.resize(d*d+1, nullptr);
    reMat.resize(2*d, nullptr);
    rhs_array.resize(d+1, nullptr);
    auto& s = this->systems.front();

    if (type == LocalAssembling_type::TimeNavierStokesRhs)
    {
      rhs_array[d] = s.rhs.block(d);
      s.rhs.reset();
    }
    else
    {
      auto blocks = s.matrix.get_blocks_uniquely();
      int nstype = TDatabase::ParamDB->NSTYPE;
      switch(nstype)
      {
        case 1:
          for (int i = 0; i < d; ++i)
          {
            reMat[i] = reinterpret_cast<MatrixD*>(blocks[1 + i].get());
          }
          break;
        case 2:
          for (int i = 0; i < d; ++i)
          {
            reMat[i] = reinterpret_cast<MatrixD*>(blocks[d + 1 + i].get());
          }
          break;
        case 3:
          for (int i = 0; i < d; ++i)
          {
            reMat[i] = reinterpret_cast<MatrixD*>(blocks[d + (d + 1) * i].get());
          }
          break;
        case 4:
          for (int i = 0; i < d; ++i)
          {
            reMat[i] = reinterpret_cast<MatrixD*>(blocks[d * (d + 1) + i].get());
          }
          break;
        case 14:
          for (int i = 0; i < d; ++i)
          {
            reMat[i] = reinterpret_cast<MatrixD*>(blocks[d * (d + 1) + i].get());
          }
          break;
      }

      for (auto rm: reMat)
      {
        if (rm != nullptr)
        {
          rm->reset();
        }
      }
    }
  }
}

/** ***************************************************************************/
template<int d>
LocalAssembling_type TNSE_ROM<d>::adapt_type(bool is_velocity,
                                             LocalAssembling_type type)
{
  // the LocalAssembling_type is changed only for Pressure assembles
  if (type == LocalAssembling_type::TCDDiffOnly
    || type == LocalAssembling_type::Custom )
  {
    return LocalAssembling_type::TCDMassOnly;
  }
  else if (!is_velocity && (type == LocalAssembling_type::NavierStokesLinear
    || type == LocalAssembling_type::TimeNavierStokesNL))
  {
    return LocalAssembling_type::TimeNavierStokesMass;
  }
  else
  {
    return type;
  }
}

/** ***************************************************************************/
template<int d>
void TNSE_ROM<d>::adapt_local_assembling(bool                 is_velocity,
                                         LocalAssembling_type type,
                                         LocalAssembling<d>&  la)
{
  using namespace std::placeholders;
  std::vector<AssembleFctParam> la_modif;

  if (is_velocity)
  {
    // in LocalAssembling::set_parameters_for_tnse() with TimeNavierStokesNL:
    // assembling routine for diffusion is added at the beginning and must be
    // removed here
    if (type == LocalAssembling_type::TimeNavierStokesNL)
    {
      la_modif = la.getAssemblingRoutines();
      la_modif.erase(la_modif.begin());
      la.replace_local_assembling(la_modif);
      return;
    }
    else
    {
      return;
    }
  }

  // modifications for Pressure assembles
  if (type == LocalAssembling_type::TCDMassOnly)
  {
    return;
  }
  else if (type == LocalAssembling_type::TCDDiffOnly)
  {
    la_modif.push_back(std::bind(TCDDiffOnlyScale<d>,
                                 _1, _2, _3, 1., _5, _6, _7, _8) );
  }
  else if (type == LocalAssembling_type::Custom)
  {
    la_modif.push_back(std::bind(TCDDiffOnlyScale<d>,
                                 _1, _2, _3, _4, _5, _6, _7, _8) );
  }
  else
  {
    // adapt for the pressure derivatives
    using MultiIndex_vector = typename Template_names<d>::MultiIndex_vector;
#ifdef __3D__
    MultiIndex_vector p_deriv {MultiIndex3D::D100, MultiIndex3D::D010,
                               MultiIndex3D::D001};
#else
    MultiIndex_vector p_deriv {MultiIndex2D::D10, MultiIndex2D::D01};
#endif
    std::vector<int> FESpaceNumber(d, 1);
    la.add_derivatives(p_deriv, FESpaceNumber);

    if (type == LocalAssembling_type::NavierStokesLinear)
    {
      la_modif.push_back(std::bind(NSPROM_Divergence<d>,
                                   _1, _2, _3, _4, _5, _6, _7, _8) );
    }
    if (type == LocalAssembling_type::TimeNavierStokesNL)
    {
      la_modif.push_back(std::bind(NSPROM_NonlinearTerm<d>,
                                  _1, _2, _3, _4, _5, _6, _7, _8) );
    }
    if (type == LocalAssembling_type::TimeNavierStokesRhs)
    {
      la_modif.push_back(std::bind(NSPROM_RightHandSide<d>,
                                  _1, _2, _3, _4, _5, _6, _7, _8) );
    }
  }
  la.replace_local_assembling(la_modif);
}

/** ***************************************************************************/
template<int d>
bool TNSE_ROM<d>::reduce(LocalAssembling_type type,
                         bool                 is_velocity,
                         int&                 new_term_idx,
                         BlockFEMatrix*       pressure_matrix)
{
  bool is_new_term = false;
  auto& s = this->systems.front();
  std::pair<size_t,size_t> upper_left;
  std::pair<size_t,size_t> lower_right;
  size_t n_row = s.mass_matrix.get_n_cell_rows();
  size_t n_col = s.mass_matrix.get_n_cell_columns();

  if (is_velocity)
  {
    upper_left = std::make_pair(0, 0);
    lower_right = std::make_pair(n_row-2, n_col-2);
  }
  else
  {
    upper_left = std::make_pair(n_row-1, 0);
    lower_right = std::make_pair(n_row-1, n_col-2);
  }

  if (type == LocalAssembling_type::TimeNavierStokesRhs)
  {
    if (is_velocity)
    {
      ROM_u.reduce(s.rhs.get_vector_entries()->begin(),
                   s.rhs.get_vector_entries()->end() - s.rhs.length(d),
                   rhs_src_u_r);
    }
    else
    {
      ROM_p.reduce(s.rhs.get_vector_entries()->end() - s.rhs.length(d),
                   s.rhs.get_vector_entries()->end(),
                   rhs_src_p_r);
    }
  }
  else if (type == LocalAssembling_type::TimeNavierStokesMass)
  {
    ROM_u.reduce(s.mass_matrix.get_combined_submatrix(upper_left, lower_right),
                 mass_mat_u_r);
  }
  else if (type == LocalAssembling_type::Custom)
  {
    ROM_p.reduce(pressure_matrix->get_combined_matrix(), diff_mat_p_r);
  }
  else
  {
    DenseMatrix* matrix_ptr = nullptr;
    std::vector<double>* vector_ptr = nullptr;
    if (type == LocalAssembling_type::NavierStokesLinear)
    {
      matrix_ptr = is_velocity ? &diff_mat_u_r : &div_mat_p_r;
      vector_ptr = is_velocity ? &diff_mean_u_r : &div_mean_p_r;
    }
    else if (type == LocalAssembling_type::TimeNavierStokesNL)
    {
      // related to the mean velocity function
      if (new_term_idx < 0)
      {
        matrix_ptr = is_velocity ? &conv_lin_mat_u_r : &conv_lin_mat_p_r;
        vector_ptr = is_velocity ? &conv_lin_mean_u_r : &conv_lin_mean_p_r;
      }
      // related to the rank_u velocity POD basis functions
      else
      {
        matrix_ptr = is_velocity ? &conv_non_lin_mat_u_r[new_term_idx]
                                : &conv_non_lin_mat_p_r[new_term_idx];
        vector_ptr = is_velocity ? &conv_non_lin_mean_u_r[new_term_idx]
                                : &conv_non_lin_mean_p_r[new_term_idx];
      }
      new_term_idx++;
      // check if the first term related to the mean velocity function
      // and the rank_u terms related to the velocity POD basis functions
      // have been computed or not yet
      is_new_term = (new_term_idx==rank_u) ? false : true;
    }

    if (is_velocity)
    {
      ROM_u.reduce(s.matrix.get_combined_submatrix(upper_left, lower_right),
                   *matrix_ptr);
    }
    else
    {
      ROM_p.reduce(s.matrix.get_combined_submatrix(upper_left, lower_right),
                   *matrix_ptr, ROM_u.get_basis());
    }

    if (db["pod_fluctuations_only"])
    {
      if (is_velocity)
      {
        ROM_u.reduce_mat_mean(s.matrix.get_combined_submatrix(upper_left,
                                                              lower_right),
                              *vector_ptr);
      }
      else
      {
        ROM_p.reduce_mat_mean(s.matrix.get_combined_submatrix(upper_left,
                                                              lower_right),
                              *vector_ptr, &ROM_u.get_snaps_avr()[0]);
      }
    }
  }
  return is_new_term;
}


/** ***************************************************************************/
template<int d>
void TNSE_ROM<d>::solve(bool is_velocity)
{
  if (is_velocity)
  {
    if (sys_mat_u_r.getNRows() != (int)sys_rhs_u_r.size())
    {
      ErrThrow("Dimension mismatch during A*x = rhs for velocity resolution: ",
              sys_mat_u_r.getNRows(), " rows in the system matrix and ",
              sys_rhs_u_r.size(), " elements in the rhs.");
    }

    sys_mat_u_r.decomposeLU();

    // copy sys_rhs_r and keep it for compute_residuals()
    std::vector<double> sys_rhs_u_r_tmp(sys_rhs_u_r);
    sys_mat_u_r.solve(&sys_rhs_u_r_tmp[0]); //solution stored in sys_rhs_u_r_tmp
    sol_u_r = std::move(sys_rhs_u_r_tmp);
  }
  else
  {
    prepare_system_pressure();

    if (sys_mat_p_r.getNRows() != (int)sys_rhs_p_r.size())
    {
      ErrThrow("Dimension mismatch during A*x = rhs for pressure resolution: ",
              sys_mat_p_r.getNRows(), " rows in the system matrix and ",
              sys_rhs_p_r.size(), " elements in the rhs.");
    }

    if (!sys_mat_p_r.getEntriesLU() || !sys_mat_p_r.getPivotsLU())
    {
      sys_mat_p_r.decomposeLU();
    }

    sys_mat_p_r.solve(&sys_rhs_p_r[0]); //solution stored in sys_rhs_p_r_tmp
    sol_p_r = std::move(sys_rhs_p_r);

    ROM_p.get_full_solution(this->sol_p_r,
                            this->get_solution().get_entries()
                            + d * this->get_velocity_space()->get_n_dof());
  }
}

/* ************************************************************************** */
template <int d>
bool TNSE_ROM<d>::stop_it(unsigned int it_counter)
{
  bool is_velocity = true;
  compute_residuals();

  // check if minimum number of iterations was performed already
  size_t min_it = TimeNavierStokes<d>::db["nonlinloop_minit"];
  if (it_counter < min_it)
  {
    return false;
  }

  size_t max_it = TimeNavierStokes<d>::db["nonlinloop_maxit"];
  double limit = TimeNavierStokes<d>::db["nonlinloop_epsilon"];
  bool linearized_scheme = this->linearized();
  if (TimeNavierStokes<d>::db["nonlinloop_scale_epsilon_with_size"])
  {
    limit *= std::sqrt(sol_u_r.size());
    Output::print("stopping tolerance for nonlinear iteration ", limit);
  }

  if (residual_u_r <= limit || it_counter == max_it
    || (linearized_scheme && it_counter > 0))
  {
    ROM_u.get_full_solution(this->sol_u_r, this->get_solution().get_entries());

    solve(!is_velocity);

// FEM residual computation
//     if (fem_residuals)
//       TimeNavierStokes<d>::compute_residuals();

    auto& s = this->systems.front();
    s.solution_m2 = s.solution_m1;
    s.solution_m1 = s.solution;

    this->adjust_pressure();
    return true;
  }
  else
  {
    return false;
  }
}

/* ************************************************************************** */
template <int d>
void TNSE_ROM<d>::compute_residuals(void)
{
  // compute reduced velocity defect = sys_rhs_u_r - sys_mat_u_r * sol_u_r
  std::vector<double> defect(sys_rhs_u_r.size());
  std::transform(sys_rhs_u_r.begin(), sys_rhs_u_r.end(),
                  sys_mat_u_r.multiply(&sol_u_r).begin(),
                  defect.begin(), std::minus<double>());
  residual_u_r = std::sqrt(std::inner_product(defect.begin(),
                                              defect.end(),
                                              defect.begin(), 0.));
}


/** ***************************************************************************/
template<int d>
void TNSE_ROM<d>::print_fem_residuals(int loop_index)
{
  if (!fem_residuals)
  {
    return;
  }

  using namespace std;
  std::stringstream s;
  s << "time step " << loop_index << " corresponding fem ";
  s <<left << setprecision(10) << setw(15);
  s << this->get_residuals();
  Output::print<1>(s.str());
}


#ifdef __3D__
template class TNSE_ROM<3>;
#else
template class TNSE_ROM<2>;
#endif
