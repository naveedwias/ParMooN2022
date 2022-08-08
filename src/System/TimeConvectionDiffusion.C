#include "TimeConvectionDiffusion.h"
#include "Database.h"
#include "ConvDiff.h"
#include "Domain.h"
#include "Multigrid.h"
#ifdef __2D__
 #include "Assemble2D.h"
 #include "SquareMatrix2D.h"
 #include "AuxParam2D.h"
#else
 #include "Assemble3D.h"
 #include "SquareMatrix3D.h"
 #include "AuxParam3D.h"
#endif
#include "LocalProjection.h"
#ifdef _MPI
 #include "ParFECommunicator3D.h"
 #include "MumpsWrapper.h"
#endif
#include "MainUtilities.h"
#include <cmath>
#include "CD_local_assembling_routines.h"


/* ************************************************************************* */
template <int d>
ParameterDatabase TimeConvectionDiffusion<d>::default_tcd_database(
  bool complete)
{
  Output::print<5>("creating a default TCD parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default TCD2D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("TimeConvectionDiffusion parameter database");

  if(complete)
  {
    db.merge(TDomain::default_domain_parameters());
    if(d == 3)
      db.add_nested_database(TDomain::default_sandwich_grid_parameters());
    db.merge(Example_TimeCD::default_example_database());
  }
  db.merge(ParameterDatabase::default_output_database());
  db.merge(default_afc_tcd_database());
  db.merge(LocalAssembling<d>::default_local_assembling_database());
  db.merge(TimeDiscretization::default_TimeDiscretization_database());
  if(complete)
  {
    db.merge(CheckpointIO::default_checkpoint_io_database());
    db.merge(Solver<>::default_solver_database());
  }
  db["problem_type"] = 2; // means time-dependent convection-diffusion
  return db;
}

/* ************************************************************************* */
template <int d>
TimeConvectionDiffusion<d>::SystemPerGrid::SystemPerGrid(
  const Example_TimeCD& example, const TCollection& coll, int ansatz_order)
: fe_space(new FESpace(&coll, "space", example.get_bc(0), ansatz_order))
{
#ifdef _MPI
  fe_space->get_communicator().print_info();
#endif
#ifdef __3D__
  stiffness_matrix = BlockFEMatrix::CD3D(fe_space);
  mass_matrix = BlockFEMatrix::CD3D(fe_space);
#else
  stiffness_matrix = BlockFEMatrix::CD2D(fe_space);
  mass_matrix = BlockFEMatrix::CD2D(fe_space);
#endif

  rhs = BlockVector(stiffness_matrix, true);
  solution = BlockVector(stiffness_matrix, false);
  fe_function = FEFunction(fe_space, "c", solution.get_entries());
  solution_m1 = BlockVector(stiffness_matrix, false);
  u_m1 = FEFunction(fe_space,"um1", solution_m1.get_entries());
  solution_m2 = BlockVector(stiffness_matrix, false);
  u_m2 = FEFunction(fe_space,"um2", solution_m2.get_entries());
}

/* ************************************************************************* */
template <int d>
TimeConvectionDiffusion<d>::TimeConvectionDiffusion(
  const TDomain& domain, const ParameterDatabase& param_db)
 : TimeConvectionDiffusion<d>(domain, param_db, Example_TimeCD(param_db))
{
}

/* ************************************************************************* */
template <int d>
TimeConvectionDiffusion<d>::TimeConvectionDiffusion(
  const TDomain& domain, const ParameterDatabase &param_db,
  const Example_TimeCD& ex)
 : db(default_tcd_database()), solver(param_db), systems(), example(ex),
   old_rhs(), errors({}), outputWriter(param_db),
   time_stepping_scheme(param_db), rhs_from_time_disc(),
   checkpoint_io(param_db)
{
  this->db.merge(param_db, false); // update this database with given values
  this->check_and_set_parameters();
  
  bool usingMultigrid = this->solver.is_using_multigrid();
  auto collections = domain.get_grid_collections();
  int ansatz_order = TDatabase::ParamDB->ANSATZ_ORDER;
  if(!usingMultigrid)
  {
    // the given collection for particular cell
    TCollection& cellCollection = *collections.front();
    systems.emplace_back(example, cellCollection, ansatz_order);
    // initial condition on the solution
    this->systems.front().fe_function.Interpolate(example.get_initial_cond(0));
  }
  else
  {
    auto multigrid = this->solver.get_multigrid();
    size_t nMgLevels = multigrid->get_n_geometric_levels();
    size_t n_grids = collections.size();
    if(n_grids < nMgLevels)
    {
      ErrThrow("Multigrid: expected ", nMgLevels, " collections, ", n_grids,
               " provided.");
    }
    // remove not needed coarser grid from list of collections
    for(size_t i = nMgLevels; i < n_grids; ++i)
    {
      collections.pop_back();
    }
    std::list<BlockFEMatrix*> matrices;
    // construct all SystemPerGrid and store them
    for(auto it : collections)
    {
      systems.emplace_back(example, *it, ansatz_order);
      systems.front().fe_function.Interpolate(example.get_initial_cond(0));
      matrices.push_front(&systems.back().stiffness_matrix);
    }
    multigrid->initialize(matrices);
  }// multigrid case
  
  // print useful information
  this->output_problem_size_info();
  outputWriter.add_fe_function(&this->systems.front().fe_function);
  this->rhs_from_time_disc.copy_structure(this->systems.front().rhs);
  errors.fill(0.); // initialize the array
   // initialize L_inf errors to some negative number
  errors[4] = -1.;
  errors[5] = -1.;

  checkpoint_io.read(systems.front().solution);
}

/* ************************************************************************* */
template <int d>
void TimeConvectionDiffusion<d>::output_problem_size_info() const
{
  // print some useful information
  const FESpace& space = *this->systems.front().fe_space;
  auto coll = space.GetCollection();
#ifndef _MPI
  double hMin, hMax;
  coll->GetHminHmax(&hMin, &hMax);
  int n_cells = coll->GetN_Cells();
  int n_dof = space.get_n_dof();
#else // _MPI
  int root = 0; // root process number
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  
  int n_local_master_cells = coll->GetN_OwnCells();
  int n_cells;
  MPI_Reduce(&n_local_master_cells, &n_cells, 1, MPI_DOUBLE, MPI_SUM, root,
             MPI_COMM_WORLD);
  
  double local_hmin, local_hmax;
  coll->GetHminHmax(&local_hmin, &local_hmax);
  double hMin, hMax;
  MPI_Reduce(&local_hmin, &hMin, 1, MPI_DOUBLE, MPI_MIN, root, MPI_COMM_WORLD);
  MPI_Reduce(&local_hmax, &hMax, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
  
  auto par_comm = systems.front().fe_space->get_communicator();
  int n_dof  = par_comm.get_n_global_dof();
  if(my_rank == root)
#endif
  {
    Output::stat("TimeConvectionDiffusion",
                 "Mesh data and problem size on finest grid");
    Output::dash("n cells     : ", setw(13), n_cells);
    Output::dash("h(min, max) : ", setw(13), hMin, " ", setw(13), hMax);
    Output::dash("n dofs      : ", setw(13), n_dof);
#ifndef _MPI
    Output::dash("n dof active: ", setw(13), space.get_n_active());
#endif
  }
}

/* ************************************************************************* */
template <int d>
void TimeConvectionDiffusion<d>::check_and_set_parameters()
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
      Output::warn<2>("The parameter problem_type doesn't correspond to Time_CD."
          "It is now reset to the correct value for Time_CD (=2).");
      db["problem_type"] = 2;
    }
  }

  // an error when using ansatz order 0
  if(TDatabase::ParamDB->ANSATZ_ORDER == 0)
  {
    throw std::runtime_error("Ansatz order 0 is no use in convection diffusion "
        "reaction problems! (Vanishing convection and diffusion term).");
  }

  //////////////// Algebraic flux correction ////////////
  if(!db["algebraic_flux_correction"].is("none"))
  {//some kind of afc enabled
    if(!db["algebraic_flux_correction"].is("fem-fct-cn"))
    {
      db["algebraic_flux_correction"].set("fem-fct-cn");
      Output::warn("TimeConvectionDiffusion::check_and_set_parameters()",
                   "Parameter 'algebraic_flux_correction' changed to "
                   "'fem-fct-cn', which is the only implemented kind of "
                   "algebraic flux correction for TCD problems.");
    }

    if(TDatabase::TimeDB->TIME_DISC != 2) /// @todo are we still using this?
    {
      ErrThrow("Algebraic flux correction with FEM-FCT is implemented for "
               "Crank-Nicolson time stepping scheme only. Set "
               "TDatabase::TimeDB->TIME_DISC to 2.")
    }

    if(TDatabase::ParamDB->ANSATZ_ORDER != 1)
    {
      ErrThrow("Algebraic flux correction with FEM-FCT does only work for"
          "linear elements. Change ANSATZ_ORDER to 1!");
    }

    //make sure that galerkin discretization is used
    if (!db["space_discretization_type"].is("galerkin"))
    {//some other disctype than galerkin
      db["space_discretization_type"] = "galerkin";
      Output::warn("TimeConvectionDiffusion::check_and_set_parameters()",
                   "Parameter 'space_discretization_type' changed to 'galerkin'"
                   " because Algebraic Flux Correction is enabled.");
    }
  }
}

/* ************************************************************************* */
template <int d>
void TimeConvectionDiffusion<d>::assemble_initial_time(FEFunction *velo1, FEFunction *velo2, FEFunction *velo3)
{
  LocalAssembling_type allMatrices = LocalAssembling_type::TCDStiffMassRhs;
  for(auto &s : this->systems)
  {
    std::vector<const FEFunction*> feFunction = {&s.fe_function};
    LocalAssembling<d> la(this->db, allMatrices, feFunction,
                          example.get_coeffs());
    
    // Assemble stiffness, mass matrices and the rhs
    if (velo1 && velo2 && velo3)
      modify_and_call_assembling_routine(s, la, true, velo1, velo2, velo3);
    else
    call_assembling_routine(s, la, true);
        // apply local projection stabilization method on stiffness matrix only!
    if(db["space_discretization_type"].is("local_projection")
        && TDatabase::ParamDB->LP_FULL_GRADIENT>0)
    {
      if(TDatabase::ParamDB->LP_FULL_GRADIENT==1)
      {
        if(d == 3)
          ErrThrow("local_projection for TimeConvectionDiffusion<3> is not "
                   "implemented"); /// @todo
#ifdef __2D__
        //fetch stiffness matrix as block
        std::vector<std::shared_ptr<FEMatrix>> stiff_blocks =
            s.stiffness_matrix.get_blocks_uniquely();
        auto stiff_block = stiff_blocks.at(0).get();
        //call ultra local projection
        UltraLocalProjection((void *)stiff_block, false);
#endif // 2D
      }
      else
      {
        ErrThrow("LP_FULL_GRADIENT needs to be 1 to use LOCAL_PROJECTION");
      }

    }
  }
  
  if (!db["algebraic_flux_correction"].is("none") )
  {
    AFC_TCD_params<d> afcTcdObjectParameters;
    setAfcTcdObjectParameters(afcTcdObjectParameters);
    create_afcTcdObject(afcTcdObjectParameters);
  }

  SystemPerGrid& s = this->systems.front();
  old_rhs = s.rhs;
  s.solution_m1 = s.solution;
  s.solution_m2 = s.solution_m1;  
}

/* ************************************************************************* */
template <int d>
void TimeConvectionDiffusion<d>::assemble(FEFunction *velo1, FEFunction *velo2, FEFunction *velo3, FEFunction *sources_and_sinks)
{
  // In the case of SUPG: local assemble function itself take care of the 
  // number of matrices. One have to assemble also the weighted mass matrix 
  // which comes from the time discretization of the SUPG method. We have 
  // to assemble the Mass matrix also for each time step due to the convection 
  // field which might also depend on time as well.
  
  if(afcTcdPtr)
  {
    if(afcTcdPtr->getComplexity() == EXPLICIT)
    {
      AFC_TCD_explicit<d>* afcTcdExplicitObjectPtr;
      afcTcdExplicitObjectPtr = static_cast<AFC_TCD_explicit<d>*>(afcTcdPtr);
      
       if(afcTcdExplicitObjectPtr->doAssembling == false)
         return;
    }
    else
    {
      SystemPerGrid& s = this->systems.front();
      
      old_rhs = s.rhs;
      s.solution_m1 = s.solution;
    }
  }
  
  LocalAssembling_type la_type = LocalAssembling_type::TCDStiffRhs;
  bool assemble_both = false;
  for(auto &s : this->systems)
  {
    // call assembling routine 
    if(db["space_discretization_type"].is("supg"))
    {
      la_type = LocalAssembling_type::TCDStiffMassRhs;
      assemble_both = true;
    }
    
    std::vector<const FEFunction*> feFunction = {&s.fe_function};
    LocalAssembling<d> la(this->db, la_type, feFunction, example.get_coeffs());
    
    if (velo1 && velo2 && velo3 && sources_and_sinks)
      modify_and_call_assembling_routine(s, la, assemble_both, velo1, velo2, velo3, sources_and_sinks);
    else if (velo1 && velo2 && velo3)
      modify_and_call_assembling_routine(s, la, assemble_both, velo1, velo2, velo3);
    else
        call_assembling_routine(s, la, assemble_both);
  }  
  
  if(afcTcdPtr)
  {
    if( afcTcdPtr->getComplexity() == LINEAR )
    {
      AFC_TCD_implicit<d>* afcTcdImplicitPtr;
      afcTcdImplicitPtr = static_cast<AFC_TCD_implicit<d>*>(afcTcdPtr);
      afcTcdImplicitPtr->calculateMatAndRhs();
    }
    
    return;
  }
  // preparing the right hand side discretized by the used time
  // stepping scheme
  SystemPerGrid& s = this->systems.front();
  rhs_from_time_disc.reset();
  rhs_from_time_disc = s.rhs;
  // all matrices are available
  unsigned int n_sols = time_stepping_scheme.n_old_solutions();
  std::vector<BlockVector> old_sols(n_sols);
  old_sols[0] = s.solution_m1;
  if(old_sols.size() == 2)
    old_sols[1] = s.solution_m2;
  std::vector<BlockVector> rhs(2);
  rhs[0] = rhs_from_time_disc;
  rhs[1] = old_rhs;
  // prepare the right hand side from the previous time step
  time_stepping_scheme.prepare_rhs_from_time_disc(s.stiffness_matrix, 
                                                  s.mass_matrix, rhs, old_sols);
  rhs_from_time_disc = rhs[0];
  old_rhs = s.rhs;
  rhs_from_time_disc.copy_nonactive(s.rhs);
  
  for(auto &s : this->systems)
    time_stepping_scheme.prepare_system_matrix(s.stiffness_matrix,
                                               s.mass_matrix);
  s.solution.copy_nonactive(s.rhs);
}

/* ************************************************************************* */
template <int d>
void TimeConvectionDiffusion<d>::solve()
{
  if(afcTcdPtr)
  {
    if( afcTcdPtr->getComplexity() == NONLINEAR )
    {
      fem_fct_nonlinear_loop();
      return;
    }
    else if( afcTcdPtr->getComplexity() == EXPLICIT )
    {
      calculateNewSolutionExplicitly();
      return;
    }
  }
  
  SystemPerGrid& s = systems.front();
#ifndef _MPI
  solver.solve(s.stiffness_matrix, rhs_from_time_disc, s.solution);
#else
  if(solver.get_db()["solver_type"].is("direct"))
  {
    MumpsWrapper mumps_wrapper(s.stiffness_matrix);
    mumps_wrapper.solve(rhs_from_time_disc, s.solution);
  }
  else
    // same as sequential
    solver.solve(s.stiffness_matrix, rhs_from_time_disc, s.solution);
#endif

  compute_residuals();

  if(afcTcdPtr)
    if(afcTcdPtr->projectionOntoAdmissibleValues)
      afcTcdPtr->projectOntoAdmissibleValues(s.solution);
  
  // restore stiffness matrix
  if(db["algebraic_flux_correction"].is("none"))
  {
    for(auto &s : this->systems)
      time_stepping_scheme.reset_linear_matrices(s.stiffness_matrix,
                                                 s.mass_matrix);
  }

  s.solution_m2 = s.solution_m1;
  s.solution_m1 = s.solution;
}

/* ************************************************************************** */
template <int d>
void TimeConvectionDiffusion<d>::compute_residuals()
{
  if(Output::getVerbosity()<2)
    return;

  SystemPerGrid& s = this->systems.front();
  bool i_am_root = true;
#ifdef _MPI
  int my_rank;
  int root = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  i_am_root = (my_rank == root);
  const TParFECommunicator3D& comm = s.fe_space->get_communicator();
  // the residuals before and after solving will be computed and
  // require consistency level of 2
  comm.consistency_update(s.solution_m1.get_entries(), 2);
  comm.consistency_update(s.solution.get_entries(), 2);
#endif

  // compute residuals before and after solving
  BlockVector defect_init = rhs_from_time_disc;
  BlockVector defect_fin  = rhs_from_time_disc;
  s.stiffness_matrix.apply_scaled_add(s.solution_m1, defect_init, -1.);
  s.stiffness_matrix.apply_scaled_add(s.solution, defect_fin, -1.);

  double residual_init;
  double residual_fin;
#ifdef _MPI
  residual_init = defect_init.norm({&comm});
  residual_fin  = defect_fin.norm({&comm});
#else
  residual_init = defect_init.norm();
  residual_fin  = defect_fin.norm();
#endif

  if(i_am_root)
  {
    Output::print("initial residual: ", std::setprecision(14), residual_init);
    Output::print("final residual  : ", std::setprecision(14), residual_fin);
  }
}

/* ************************************************************************** */
template <int d>
void TimeConvectionDiffusion<d>::fem_fct_nonlinear_loop()
{
  SystemPerGrid& sys = systems.front();
  
#ifdef _MPI
  const TParFECommunicator3D& comm = sys.fe_space->get_communicator();
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  bool i_am_root = (my_rank == 0);
  int numberOfDofs = comm.get_n_global_dof();
#else
  bool i_am_root = true;
  int numberOfDofs = sys.solution.length();
#endif
  
  double toleranceForResidual = db["afc_nonlinloop_epsilon"];
  toleranceForResidual *= std::sqrt(numberOfDofs);
  if(i_am_root)
    Output::print<1>("Stopping tolerance for residual in nonlinear iteration is ",
                  toleranceForResidual);
  
  unsigned int maxNumberOfInterations = get_db()["afc_nonlinloop_maxit"];
  if(i_am_root)
    Output::print<1>("Maximum number of nonlinear iterations is ",
                  maxNumberOfInterations);
  
  BlockVector residuals;
  double residual;
  
  TimeDiscretization& tss = get_time_stepping_scheme();
  
  
  for(tss.nonlin_iteration = 0; ; tss.nonlin_iteration++)
  {
    AFC_TCD_implicit<d>* afcTcdImplicitPtr;
    afcTcdImplicitPtr = static_cast<AFC_TCD_implicit<d>*>(afcTcdPtr);
    
    afcTcdImplicitPtr->readSolPrevNonlin(sys.solution);
    afcTcdImplicitPtr->calculateMatAndRhs();
    
    residuals = rhs_from_time_disc;
    sys.stiffness_matrix.apply_scaled_add(sys.solution , residuals ,-1.0);
    
    #ifdef _MPI
    residual = residuals.norm({&comm});
    #else
    residual = residuals.norm();
    #endif
    
    if(i_am_root)
      Output::print<1>("Residual in step ", tss.nonlin_iteration, " is ", residual);
    
    if(residual < toleranceForResidual && tss.nonlin_iteration > 0 )
      break;
    else 
      if(tss.nonlin_iteration >= maxNumberOfInterations)
      {
        Output::print<1>("Maximum number of iterations reached");
        break;
      }
      
    #ifndef _MPI
    solver.solve(sys.stiffness_matrix, rhs_from_time_disc, sys.solution);
    #else
    if(solver.get_db()["solver_type"].is("direct"))
    {
      MumpsWrapper mumps_wrapper(sys.stiffness_matrix);
      mumps_wrapper.solve(rhs_from_time_disc, sys.solution);
    }
    else
      // same as sequential
      solver.solve(sys.stiffness_matrix, rhs_from_time_disc, sys.solution);
    #endif
    
    if(afcTcdPtr->projectionOntoAdmissibleValues)
      afcTcdPtr->projectOntoAdmissibleValues(sys.solution);
  }
  
}


template <int d>
void TimeConvectionDiffusion<d>::calculateNewSolutionExplicitly()
{
  SystemPerGrid& sys = systems.front();
  TimeDiscretization& tss = time_stepping_scheme;
  
  #ifdef _MPI
  const TParFECommunicator3D& comm = sys.fe_space->get_communicator();
  #endif
  
  AFC_TCD_explicit<d>* afcTcdExplicitPtr;
  afcTcdExplicitPtr = static_cast<AFC_TCD_explicit<d>*>(afcTcdPtr);
  
  // same for both coefficients
  if(afcTcdPtr->getLimiter() == ZALESAK)
    afcTcdExplicitPtr->setZalesakSolutionTime(tss.current_time_);
  
  // 1. Runge-Kutta coefficient
  double oldTime = tss.current_time_ - tss.get_step_length();
  TDatabase::TimeDB->CURRENTTIME = oldTime;
  tss.current_time_ = oldTime;
  BlockVector solToPass = sys.solution;
  
  afcTcdExplicitPtr->doAssembling = true;
  afcTcdExplicitPtr->readSolForRungeKutta(solToPass);
  afcTcdExplicitPtr->readTimeStep(tss.get_step_length());
  
  assemble();
  BlockVector rk1Coeff = afcTcdExplicitPtr->getRungeKuttaCoeff();
  
  
  // 2. Runge-Kutta coefficient
  double newTime = oldTime + tss.get_step_length();
  TDatabase::TimeDB->CURRENTTIME = newTime;
  tss.current_time_ = newTime;
  
  solToPass.operator+=(rk1Coeff);
  
  // prepare boundary values for 2. coefficient and for new solution
  BlockVector boundaryValues;
  boundaryValues.copy_structure(rk1Coeff);
  
  std::vector<double> boundaryValuesFEFunctionEntries(rk1Coeff.length());
  FEFunction boundaryValuesFEFunction (sys.fe_space,
                                       "boundary_values",
                                       boundaryValuesFEFunctionEntries.data());
  boundaryValuesFEFunction.SetDirichletBC(example.get_bc(0), example.get_bd(0) );
  memcpy( boundaryValues.get_entries(),
          boundaryValuesFEFunction.GetValues(),
          sizeof(double)*boundaryValuesFEFunction.GetLength() );
  
  
  solToPass.copy_nonactive(boundaryValues);
  
  if(afcTcdPtr->projectionOntoAdmissibleValues)
    afcTcdPtr->projectOntoAdmissibleValues(solToPass);
  
  #ifdef _MPI
  comm.consistency_update(solToPass.get_entries(), 2);
  #endif
  
  assemble();
  afcTcdExplicitPtr->readSolForRungeKutta(solToPass);
  BlockVector rk2Coeff = afcTcdExplicitPtr->getRungeKuttaCoeff();
  
  
  // new solution
  for(unsigned int i = 0; i < rk1Coeff.length(); i++)
  {
    sys.solution[i] = sys.solution[i] + 0.5*( rk1Coeff[i] + rk2Coeff[i] );
  }
  sys.solution.copy_nonactive(boundaryValues);
  
  if(afcTcdPtr->projectionOntoAdmissibleValues)
    afcTcdPtr->projectOntoAdmissibleValues(sys.solution);
  
  #ifdef _MPI
  comm.consistency_update(sys.solution.get_entries(), 2);
  #endif
  
  
  // preparation for new interation
  afcTcdExplicitPtr->doAssembling = false;
}


template <int d>
void TimeConvectionDiffusion<d>::setAfcTcdObjectParameters
(AFC_TCD_params<d>& params)
{
  params.checkParamsConsistency(db);
  
  SystemPerGrid& s = systems.front();
      
  params.stiffMatPtr = s.stiffness_matrix.get_blocks_uniquely().at(0).get();
  params.massMatPtr = s.mass_matrix.get_blocks().at(0).get();
  params.rhsPtr = &s.rhs;
  params.rhsOldPtr = &old_rhs;
  params.rhsFinalPtr = &rhs_from_time_disc;
  params.solOldPtr = &s.solution_m1;
  params.interationNumberPtr = &time_stepping_scheme.nonlin_iteration;
  params.solFESpacePtr = systems.front().fe_space;
  
  params.oldTime = time_stepping_scheme.current_time_;
  params.timeStep = time_stepping_scheme.get_step_length();
  
  params.limiter = db["afc_limiter"].is("zalesak") ? ZALESAK : MONOLITHIC;
  
  params.solDotUsedInZalesak = db["solDotUsedInZalesak"];
  
  params.boundaryConditionFunction = example.get_bc(0);
  params.boundaryValuesFunction = example.get_bd(0);
  
  const int positionOfBoundaryTimeDerivatives = 1;
  
  if(example.get_n_bd_fcts() > positionOfBoundaryTimeDerivatives)
    params.boundaryTimeDerivativeFct = example.get_bd(positionOfBoundaryTimeDerivatives);
  else
  {
    if(params.limiter == MONOLITHIC ||
      (params.limiter == ZALESAK && params.solDotUsedInZalesak))
    Output::print("WARNING: No boudary values of time derivative provided."
                  " If Dirichlet boundary conditions are prescribed,"
                  " these derivatives will be calculated"
                  " via finite difference.");
    params.boundaryTimeDerivativeFct = nullptr;
  }
  
  if(db["afc_fct_scheme"].is("linear"))
    params.complexity = LINEAR;
  else if(db["afc_fct_scheme"].is("non-linear"))
    params.complexity = NONLINEAR;
  else if(db["afc_fct_scheme"].is("explicit"))
    params.complexity = EXPLICIT;
  else
    ErrThrow("Unknow type of afc_fct_scheme.",
             "Possible types: linear, non-linear, explicit.");
  
  if(db["time_discretization"].is("crank_nicolson"))
  {
    params.timeSteppingScheme = CRANK_NICOLSON;
    params.theta = 0.5;
  }
  else if(db["time_discretization"].is("backward_euler"))
  {
    params.timeSteppingScheme = BACKWARD_EULER;
    params.theta = 1.0;
  }
  else if(db["time_discretization"].is("Runge_Kutta_Heun"))
    params.timeSteppingScheme = RUNGE_KUTTA_HEUN;
  
  switch((int) db["afc_prelimiter"])
  {
    case 1:
      params.prelimiter = AFC_TCD_Prelimiter::MIN_MOD;
      break;
    case 2:
      params.prelimiter = AFC_TCD_Prelimiter::GRAD_DIRECTION;
      break;
    case 3:
      params.prelimiter = AFC_TCD_Prelimiter::BOTH;
      break;
    default:
      params.prelimiter = AFC_TCD_Prelimiter::NONE;
  }
  
  params.toleranceForZero = db["grad_direction_tolerance"];
  
  if(db["projection_onto_admissible_interval"].is("yes"))
  {
    params.projectionOntoAdmissibleValues = true;
    params.admissibleMinimum = db["admissible_minimum"];
    params.admissibleMaximum = db["admissible_maximum"];
  }
  else
    params.projectionOntoAdmissibleValues = false;
}


template <int d>
void TimeConvectionDiffusion<d>::create_afcTcdObject
(AFC_TCD_params<d>& passedParameters)
{ 
  rowsColsRecordPtr = new RowsColsRecord<d>(passedParameters);
  bcSetterPtr = new TCD_BC_setter<d>(passedParameters);
  
  passedParameters.RowsColsRecordPtr = rowsColsRecordPtr;
  passedParameters.tcdBcSetterPtr = bcSetterPtr;
  
  switch(passedParameters.complexity)
  {
    case EXPLICIT:
        afcTcdPtr = new AFC_TCD_explicit<d>(passedParameters);
      break;
    
    default:
      if(passedParameters.complexity == LINEAR && passedParameters.limiter == MONOLITHIC)
        ErrThrow("Monolithic limiter with complexity 'linear' is not implemented!");
      
        afcTcdPtr = new AFC_TCD_implicit<d>(passedParameters);
  }
}

template <int d>
void TimeConvectionDiffusion<d>::destroy_afcTcdObject()
{
  delete rowsColsRecordPtr;
  delete bcSetterPtr;
  
  switch(afcTcdPtr->getComplexity())
  {
    case EXPLICIT:
      AFC_TCD_explicit<d>* afcTcd2ExplicitPtr;
      afcTcd2ExplicitPtr = static_cast<AFC_TCD_explicit<d>*>(afcTcdPtr);
      delete afcTcd2ExplicitPtr;
      break;
    
    default:
      AFC_TCD_implicit<d>* afcTcd2ImplicitPtr;
      afcTcd2ImplicitPtr = static_cast<AFC_TCD_implicit<d>*>(afcTcdPtr);
      delete afcTcd2ImplicitPtr;
  }
}


template<int d>
void TimeConvectionDiffusion<d>::printValuesAtPoint
(double x, double y, double z, const std::string& fileName)
{
  #ifdef _MPI
  int my_rank;
  int root = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  #endif // _MPI
  
  FEFunction& solutionFct = get_function();
  
  #ifdef __2D__
  double values[3];
  solutionFct.FindGradient(x, y, values);
  #else
  std::vector<double> values(4);  
  solutionFct.FindGradient(x, y, z, values);
  #endif
  
  double value = values[0];
  
  #ifdef _MPI    
  double finalValue;
  
  MPI_Reduce(&value, &finalValue, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
  
  value = finalValue;
  #endif
  
  #ifdef _MPI
  if(my_rank == root)
  #endif
  {
    if(!fileName.empty())
    {
      std::ofstream file;
      file.open(fileName, std::ifstream::app);
        
      file << get_time_stepping_scheme().current_time_ << " " << value << std::endl;
      
      file.close();
    }
    else
      Output::print<1>("Solution at point (", x, ", ", y, ", ", z,"): ", value);
  }
}


/* ************************************************************************* */
template <int d>
void TimeConvectionDiffusion<d>::output()
{
  SystemPerGrid &s = this->systems.front();
  if(db["output_compute_minmax"])
  {
    s.fe_function.PrintMinMax();
  }
  bool i_am_root = true;
#ifdef _MPI
  int my_rank;
  int root = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  i_am_root = (my_rank == root);
  // computing errors as well as writing vtk files requires a minimum 
  // consistency level of 1
  s.fe_space->get_communicator().consistency_update(s.solution.get_entries(), 1);
#endif // _MPI

  //write solution for visualization 
  outputWriter.write(time_stepping_scheme.current_time_);

  // compute errors
  if(db["output_compute_errors"])
  {
    const int n_errors = 5;
    std::array<double, n_errors> locError;
#ifdef __3D__
    TAuxParam3D aux;
    MultiIndex3D all_derivatives[4] = { MultiIndex3D::D000, MultiIndex3D::D100,
                                        MultiIndex3D::D010, MultiIndex3D::D001 };
#else
    TAuxParam2D aux;
    MultiIndex2D all_derivatives[3] = { MultiIndex2D::D00, MultiIndex2D::D10,
                                        MultiIndex2D::D01 };
#endif
    const FESpace* space = s.fe_space.get();

#ifdef __3D__
    s.fe_function.GetErrors(example.get_exact(0), d+1, all_derivatives,
                            n_errors, conv_diff_l2_h1_linf_error<d>,
                            example.get_coeffs(), &aux, 1, &space,
                            locError.data());
#else
    s.fe_function.GetErrors(example.get_exact(0), d+1, all_derivatives,
                            n_errors, conv_diff_l2_h1_linf_error<d>,
                            example.get_coeffs(), &aux, 1, &space,
                            locError.data(), false,[](const TBaseCell*, int)
                            {return false;}, db);
#endif

#ifdef _MPI
    /// @todo the GetErrors method in TFEFunction3D should already do the 
    /// communication, it's surprising that in the mpi case the errors are 
    /// squared, while the square root has been taken already in the sequential 
    /// case.
    // global (across all processes) errors:
    // L2, H1, SD-error (for SUPG), DG-error
    std::vector<double> errorsReduced(n_errors);
    MPI_Reduce(locError.data(), errorsReduced.data(), n_errors-1, MPI_DOUBLE,
               MPI_SUM, root, MPI_COMM_WORLD);
    // Linf
    MPI_Reduce(&locError[n_errors-1], &errorsReduced[n_errors-1], 1,
               MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
    // correct values only on the root process!
    for(int i = 0; i < n_errors-1; ++i)
      locError[i] = std::sqrt(errorsReduced[i]);
    locError[n_errors-1] = errorsReduced[n_errors-1];
#endif // _MPI
    
// Computation of FEM-FCT norm deleted.
// If you miss this functionality, tell Ondrej Partl.
// He will add it ASAP.

    if(i_am_root)
    {
      Output::print<1>(std::setprecision(10), time_stepping_scheme.current_time_,"  L2   : ", std::setprecision(14), locError[0]);
      Output::print<1>(std::setprecision(10), time_stepping_scheme.current_time_,"  H1   : ", std::setprecision(14), locError[1]);
      Output::print<1>(std::setprecision(10), time_stepping_scheme.current_time_,"  L_inf: ", std::setprecision(14), locError[4]);
      if(db["space_discretization_type"].is("supg"))
      {
        Output::print<1>(std::setprecision(10), time_stepping_scheme.current_time_,"  SUPG : ", std::setprecision(14), locError[2]);
      }
    }

    double tau = time_stepping_scheme.get_step_length();

    // errors[1] is the previous L2 error squared
    errors[0] += (locError[0] * locError[0] + errors[1] * errors[1]) * tau*0.5;
    errors[1] = locError[0];
    // errors[3] is the previous H1 error squared
    errors[2] += (locError[1] * locError[1] + errors[3] * errors[3]) * tau*0.5;
    errors[3] = locError[1];
    if(db["space_discretization_type"].is("supg"))
    {
      // errors[i] is the previous SUPG error squared
      errors[6] += (locError[2]*locError[2] + errors[7]*errors[7]) * tau*0.5;
      errors[7] = locError[2];
    }
    // L_inf(0,T;L_inf)
    if(errors[4] < locError[4])
    {
      errors[4] = locError[4];
      t_Linf = time_stepping_scheme.current_time_;
    }
    // L_inf(0,T;L2)
    if(errors[5] < locError[0])
    {
      errors[5] = locError[0];
      t_L2 = time_stepping_scheme.current_time_;
    }

    if(i_am_root)
    {
      Output::print<1>(std::setprecision(10), time_stepping_scheme.current_time_,std::setprecision(14), "  L2(0,T;L2)      : ", std::sqrt(errors[0]));
      Output::print<1>(std::setprecision(10), time_stepping_scheme.current_time_,std::setprecision(14), "  L2(0,T;H1)      : ", std::sqrt(errors[2]));
      Output::print<1>(std::setprecision(10), time_stepping_scheme.current_time_,std::setprecision(14), "  L_inf(0,T;L_inf): ", errors[4], " at time ", t_Linf);
      Output::print<1>(std::setprecision(10), time_stepping_scheme.current_time_,std::setprecision(14), "  L_inf(0,T;L2)   : ", errors[5], " at time ", t_L2);
      if(db["space_discretization_type"].is("supg"))
      {
        Output::print<1>(std::setprecision(10), time_stepping_scheme.current_time_,"  L2(0,T;SUPG)    : ", std::setprecision(14),
                         std::sqrt(errors[6]));
      }
    }
  }
  if(this->db["compute_cut_lines"].is("yes"))
    printValuesAtPoint (1.0, 7.0/16.0, 9.0/16.0, "");
  
  checkpoint_io.write(s.solution, time_stepping_scheme);
}

/* ************************************************************************* */
template <int d>
void TimeConvectionDiffusion<d>::call_assembling_routine(
  TimeConvectionDiffusion<d>::SystemPerGrid& system, LocalAssembling<d>& la,
  bool assemble_both)
{
  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  
  int nFESpaces = 1;
  const FESpace * fe_space = system.fe_space.get();
  
  auto stiffBlock = system.stiffness_matrix.get_blocks_uniquely()[0].get();
  auto massBlock = system.mass_matrix.get_blocks_uniquely()[0].get();
  
  int nsqMatrices = assemble_both ? 2 : 1;
  std::vector<SquareMatrixD*> sqMatrices(nsqMatrices, nullptr);
  
  sqMatrices[0] = reinterpret_cast<SquareMatrixD*>(stiffBlock);
  if(assemble_both)
  {
    sqMatrices[1] = reinterpret_cast<SquareMatrixD*>(massBlock);
    sqMatrices[1]->reset();
  }
  
  if(afcTcdPtr)
  {
    if(afcTcdPtr->getLimiter() == MONOLITHIC)
    {
      FEMatrix* diffPartMatPtr = &afcTcdPtr->getDiffPartMat();
      
      nsqMatrices++;
      sqMatrices.push_back(reinterpret_cast<SquareMatrixD*>(diffPartMatPtr));
      sqMatrices.back()->reset();
      
      std::vector<AssembleFctParam> laf = la.getAssemblingRoutines();
      
      using namespace std::placeholders;
      
      laf[0] = TCDStiffWithoutDiffAndLinearSource<d>;
      laf.push_back( std::bind( TCDDiffAndLinearSourceOnly<d>,
                                _1, _2, _3, _4, _5, _6, _7, _8,
                                nsqMatrices-1 ) );
      
      la.replace_local_assembling(laf);
    }
  }
  sqMatrices[0]->reset();
  system.rhs.reset();
  
  int nRectMatrices = 0;
  int nRhs = 1;
  double *rhsEntries = system.rhs.get_entries();
  
  auto * boundary_conditions = fe_space->get_boundary_condition();
  auto * boundary_value = example.get_bd(0);
  
  // In case of FEM-FCT, no resetting of boundary conditions is needed
  // because Assemble2D/Assemble3D does not modify matrix
  // because of Dirichlet boundary conditions
#ifdef __2D__
    Assemble2D(
#else      
    Assemble3D(
#endif
               nFESpaces, &fe_space, nsqMatrices, sqMatrices.data(),  nRectMatrices,
               nullptr, nRhs, &rhsEntries, &fe_space, &boundary_conditions,
               &boundary_value, la, true);
}

// Used as a "ParamFunction" (a MooNMD specific oddity
// of the assembling process) in the following.
// ...shears away the first three parameters (x,y,z) and passes on the three next parameters
void ThreeFEParametersFunction3D(const double *in, double *out)
{
  out[0] = in[3];
  out[1] = in[4];
  out[2] = in[5];
}
// ...shears away the first three parameters (x,y,z) and passes on the four next parameters
void FourFEParametersFunction3D(const double *in, double *out)
{
  out[0] = in[3];
  out[1] = in[4];
  out[2] = in[5];
  out[3] = in[6];
}

template <int d>
void TimeConvectionDiffusion<d>::modify_and_call_assembling_routine(
    SystemPerGrid& s,
    LocalAssembling<d>& la,
    bool assemble_both,
    FEFunction* velo1,
    FEFunction* velo2,
    FEFunction* velo3,
    FEFunction* sources_and_sinks)
{
  // NOTE: velo1 to velo3 and sources_and_sinks must be defined on the right grid
  // - but currently there is no way to check that...

  // set up the input...
  std::vector<int> beginParameter = {0};

  std::vector<const FEFunction*> fe_funct(4); //fill up the new fe function array (4th entry is optional, see below)
  fe_funct[0] = &s.fe_function;
  fe_funct[1] = velo1;
  fe_funct[2] = velo2;
  fe_funct[3] = velo3;

  std::vector<int> feValueFctIndex = {1,2,3}; // to produce first fe value use fe function 1,
                                            // for second fe value use function 2 and for third value function 3




   //for all three fe value use 0th derivative,
  int N_parameters = 3; // three parameters (AFTER application of parameterFct...)
  int N_feValues = 3;   // ...all of them stem from the evaluation of fe fcts
  int N_paramFct = 1;   // dealing with them is performed by 1 ParamFct

  // chose the parameter function ("in-out function") which shears away
  // the first three "in" values (x,y,z) and passes u_x, u_y and u_z
  std::vector<ParamFct*> parameterFct = {ThreeFEParametersFunction3D};

  if(sources_and_sinks) // rhs source and sink terms are given
  {
    // NOTE: sources and sinks must be defined on the right grid
    // - but currently there is no way to check that...

    fe_funct.push_back(sources_and_sinks);
    feValueFctIndex = {1,2,3,4}; // to produce first fe value use fe function 1,
                               // for second fe value use function 2,
                               // for third fe value use function 3 and for fourth value use function 4

#ifdef __2D__
    MultiIndex_vector feValueMultiIndex(d+1, MultiIndex2D::D00);
#else
    MultiIndex_vector feValueMultiIndex(d+1, MultiIndex3D::D000);
#endif
    N_parameters = 4; // four parameters (AFTER application of parameterFct...)
    N_feValues = 4;   // ...all four of them stem from the evaluation of fe fcts
    N_paramFct = 1;   // dealing with them is performed by 1 ParamFct
    // chose the parameter function ("in-out function") which shears away
    // the first three "in" values (x,y,z) and passes u_x, u_y, u_z and f
    parameterFct = {FourFEParametersFunction3D};


    //THIS IS DOUBLE CODE, but this stuff must be performed before
    // interpolated_sources_and_sinks and entries_source_and_sinks
    // go out of scope...
    // ...and call the corresponding setters
    la.setBeginParameter(beginParameter);
    la.setFeFunctions(fe_funct); //reset - now velo comp included
    la.setFeValueFctIndex(feValueFctIndex);
    la.setFeValueMultiIndex(feValueMultiIndex);
    la.SetN_Parameters(N_parameters);
    la.setN_FeValues(N_feValues);
    la.setN_ParamFct(N_paramFct);
    la.setParameterFct(parameterFct);
    //...I expect that to do the trick.

    // step 4 - the assembling must be done before the velo functions
    // run out of scope
    call_assembling_routine(s, la, assemble_both);
  }
  else
  {

  // ...and call the corresponding setters
#ifdef __2D__
    MultiIndex_vector feValueMultiIndex(d, MultiIndex2D::D00);
#else
    MultiIndex_vector feValueMultiIndex(d, MultiIndex3D::D000);
#endif
    la.setBeginParameter(beginParameter);
    la.setFeFunctions(fe_funct); //reset - now velo comp included
    la.setFeValueFctIndex(feValueFctIndex);
    la.setFeValueMultiIndex(feValueMultiIndex);
    la.SetN_Parameters(N_parameters);
    la.setN_FeValues(N_feValues);
    la.setN_ParamFct(N_paramFct);
    la.setParameterFct(parameterFct);
  //...I expect that to do the trick.

  // step 4 - the assembling must be done before the velo functions
  // run out of scope
  call_assembling_routine(s, la, assemble_both);
  }
}



/* ************************************************************************* */
template <int d>
std::array<double, 3> TimeConvectionDiffusion<d>::get_errors() const
{
  std::array<double, 3> error_at_time_points;
  error_at_time_points.at(0) = std::sqrt(errors.at(0));
  error_at_time_points.at(1) = std::sqrt(errors.at(2));
  error_at_time_points.at(2) = errors.at(4);
//   error_at_time_points.at(3) = errors.at(5); // SUPG
  return error_at_time_points;
}

#ifdef __3D__
template class TimeConvectionDiffusion<3>;
#else
template class TimeConvectionDiffusion<2>;
#endif
