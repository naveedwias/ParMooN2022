#include <Database.h>
#include <NavierStokes_DG.h>
#include <LocalAssembling.h>
#ifdef __2D__
#include <FEFunction2D.h>
#include <Assemble2D.h>
#include <AuxParam2D.h>
#else
#include <FEFunction3D.h>
#include <Assemble3D.h>
#include "AuxParam3D.h"
#endif
#include <MainUtilities.h>
#include "Multigrid.h"
#include "Saddle_point_preconditioner.h"
#include <Assemble_DG.h>

//added 21.02.2022
#include "GridTransfer.h"
#include "Upwind.h"

//added 22.02.2002
#include "Matrix2D.h"
#include "SquareMatrix2D.h"

template <int d>
ParameterDatabase NavierStokes_DG<d>::default_NavierStokes_DG_database(bool complete)
{
  Output::print<5>("creating a default NavierStokes_DG parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default NavierStokes_DG database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("NavierStokes_DG parameter database");

  if(complete)
  {
    db.merge(TDomain::default_domain_parameters(), true);
    if(d == 3)
      db.add_nested_database(TDomain::default_sandwich_grid_parameters());
    db.merge(Example_NSE::default_example_database());
  }

  db.merge(ParameterDatabase::default_output_database(), true);
  db.merge(LocalAssembling<d>::default_local_assembling_database(), true);
  db.merge(Solver<>::default_solver_database(), true);
  if(complete)
  {
    db.merge(CheckpointIO::default_checkpoint_io_database());
    db.add_nested_database(ParameterDatabase(
      Saddle_point_preconditioner::database_name_velocity_solver));
    db.add_nested_database(ParameterDatabase(
      Saddle_point_preconditioner::database_name_pressure_solver));
  }

  db.merge(default_dg_database(), true);

  return db;
}


/** ************************************************************************ */
template <int d>
NavierStokes_DG<d>::System_per_grid::System_per_grid(
  const Example_NSE& example, const TCollection& coll,
  std::pair<int,int> velocity_pressure_orders)
 : velocity_space(new FESpace(&coll, "u", example.get_bc(0),
                              velocity_pressure_orders.first)),
   pressure_space(new FESpace(&coll, "p", example.get_bc(1),
                              velocity_pressure_orders.second))
{
  matrix = BlockFEMatrix::NSE_DG(this->velocity_space, this->pressure_space);

  rhs = BlockVector(this->matrix, true);
  solution = BlockVector(this->matrix, false);

  u = FEFunction(this->velocity_space, "u", this->solution.block(0));
  p = FEFunction(this->pressure_space, "p", this->solution.block(1));
}

/** ************************************************************************ */
template <int d>
NavierStokes_DG<d>::NavierStokes_DG(const TDomain& domain, const ParameterDatabase& param_db)
 : NavierStokes_DG<d>(domain, param_db, Example_NSE(param_db))
{
}

/** ************************************************************************ */
template <int d>
NavierStokes_DG<d>::NavierStokes_DG(const TDomain& domain, const ParameterDatabase& param_db,
                 const Example_NSE ex)
 : systems(), example(ex), db(default_NavierStokes_DG_database(false)),
   outputWriter(param_db), checkpoint_io(param_db), solver(param_db), errors()
{
  // get the parameters to control the behavior of this class
  this->db.merge(param_db, false);
  // make sure all parameters in the database are set consistently
  this->set_parameters();
  std::pair <int,int>
      velocity_pressure_orders(TDatabase::ParamDB->VELOCITY_SPACE,
                               TDatabase::ParamDB->PRESSURE_SPACE);
  // The construction of the members differs, depending on whether a multigrid
  // solver will be used or not.
  bool usingMultigrid = solver.is_using_multigrid();

  // a collection is basically only an array of cells, which is needed to create
  // a finite element space
  auto collections = domain.get_grid_collections();
  if(!usingMultigrid)
  {
    // Get the collection on the finest grid
    auto coll = collections.front();
    // create finite element spaces and functions, a matrix, rhs, and solution
    this->systems.emplace_back(example, *coll, velocity_pressure_orders);
  }
  else
  {
    ErrThrow("The multigrid implementation for NavierStokes_DG type problems is not yet "
             "finished. You have to pick a different solver/preconditioner.");
    // Construct multigrid object
    auto mg = solver.get_multigrid();
    if(mg->is_using_mdml())
    {
      ErrThrow("mdml not supported for NavierStokes_DG type problems");
    }

    //Check whether number of given grids is alright
    size_t n_geo_multigrid_levels = mg->get_n_geometric_levels();
    size_t n_grids = collections.size();
    if(n_geo_multigrid_levels > n_grids )
      ErrThrow("Wrong number of grids for multigrid! I was expecting ",
               n_geo_multigrid_levels, " geometric grids but only got ", n_grids,".");
    // remove not needed coarser grid from list of collections
    for(size_t i = n_geo_multigrid_levels; i < n_grids; ++i)
    {
      collections.pop_back();
    }


    // Construct systems per grid and store them, finest level first
    std::list<BlockFEMatrix*> matrices;
    for(auto coll : collections)
    {
      systems.emplace_back(example, *coll, velocity_pressure_orders);
      //prepare input argument for multigrid object
      matrices.push_front(&systems.back().matrix);
    }
    mg->initialize(matrices);
  }

  outputWriter.add_fe_function(&this->get_velocity());
  outputWriter.add_fe_function(&this->get_pressure());
  checkpoint_io.read(this->get_solution());

  // print out some information on the finite element space
  auto& v_space = *this->systems.front().velocity_space;
  auto& p_space = *this->systems.front().pressure_space;
  int n_u = v_space.get_n_dof();
  int n_u_active = v_space.get_n_active();
  int n_p = p_space.get_n_dof();
  int n_dof = n_u + n_p;
  int n_cells = v_space.GetCollection()->GetN_Cells();

  Output::stat<1>("NavierStokes_DG", "Mesh data and problem size in ", d, "D");
  Output::dash<1>("cells                        : ", setw(5), n_cells);
  Output::dash<1>("dof velocity (vector-valued) : ", setw(5), n_u);
  Output::dash<1>("active dof velocity          : ", setw(5), n_u_active);
  Output::dash<1>("dof pressure                 : ", setw(5), n_p);
  Output::dash<1>("dof all                      : ", setw(5), n_dof);
}

/** ************************************************************************ */
template <int d>
void NavierStokes_DG<d>::set_parameters()
{
  // check if given velocity space is supported
  switch(TDatabase::ParamDB->VELOCITY_SPACE)
  {
    case 1000: case 1001: case 1002: case 1003: case 1011: case 1012: case 1013:
      break;
    default:
      ErrThrow("unknown velocity space for NavierStokes_DG");
      break;
  }
  if(TDatabase::ParamDB->PRESSURE_SPACE == -4711)
  {
    switch(TDatabase::ParamDB->VELOCITY_SPACE)
    {
      case 1000: // Raviart-Thomas, order 0
        TDatabase::ParamDB->PRESSURE_SPACE = 0;
        break;
      case 1001: // Raviart-Thomas, order 1
        TDatabase::ParamDB->PRESSURE_SPACE = -11;
        break;
      case 1002: // Raviart-Thomas, order 2
        TDatabase::ParamDB->PRESSURE_SPACE = -12;
        break;
      case 1003: // Raviart-Thomas, order 3
        TDatabase::ParamDB->PRESSURE_SPACE = -13;
        break;
      case 1011: // Brezzi-Douglas-Marini, order 1
        TDatabase::ParamDB->PRESSURE_SPACE = 0;
        break;
      case 1012: // Brezzi-Douglas-Marini, order 2
        TDatabase::ParamDB->PRESSURE_SPACE = -110;
        break;
      case 1013: // Brezzi-Douglas-Marini, order 3
        TDatabase::ParamDB->PRESSURE_SPACE = -120;
        break;
      default:
        ErrThrow("unknown velocity space for NavierStokes_DG");
        break;
    }
  }
}

/** ************************************************************************ */
template <int d>
void NavierStokes_DG<d>::assemble()
{
  using MatrixD = typename Template_names<d>::MatrixD;
  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  using BoundaryConditionFunction
    = typename Template_names<d>::BoundaryConditionFunction;
  using BoundaryValuesFunction
    = typename Template_names<d>::BoundaryValuesFunction;

  for(System_per_grid& s : this->systems)
  {
    // the class LocalAssembling3D which we will need next, requires an array of
    // pointers to finite element functions, i.e. TFEFunction3D **.
    std::vector<const FEFunction*> fe_functions{ &s.u, &s.p };
    // create a local assembling object which is needed to assemble the matrices
    LocalAssembling<d> la(this->db, LocalAssembling_type::NavierStokesLinear, fe_functions,
                          this->example.get_coeffs());

    // everything which follows within this for loop only has the goal to call
    // Assemble3D_VectFE at the end.
    const FESpace * v_space = s.velocity_space.get();
    const FESpace * p_space = s.pressure_space.get();

    const size_t n_fe_spaces = 2;

    const FESpace *fespmat[2] = {v_space, p_space};
    const size_t n_sq_mat = 2;
    SquareMatrixD *sq_matrices[n_sq_mat]{};

    const size_t n_rect_mat = 2;
    MatrixD *rect_matrices[n_rect_mat]{};

    const size_t n_rhs = 2;
    double *RHSs[n_rhs] = {s.rhs.block(0), s.rhs.block(1)};

    BoundaryConditionFunction * boundary_conditions[n_fe_spaces] = {
      v_space->get_boundary_condition(), p_space->get_boundary_condition() };
    std::array<BoundaryValuesFunction*, n_fe_spaces+1> non_const_bound_values;


    non_const_bound_values[0] = example.get_bd()[3];  // u * n
    non_const_bound_values[1] = example.get_bd()[2];  // p
    non_const_bound_values[2] = example.get_bd()[4];  // u * t

    std::vector<std::shared_ptr<FEMatrix>> blocks
      = s.matrix.get_blocks_uniquely();
    sq_matrices[0] = reinterpret_cast<SquareMatrixD*>(blocks.at(0).get());
    sq_matrices[1] = reinterpret_cast<SquareMatrixD*>(blocks.at(3).get());
    
    rect_matrices[0] = reinterpret_cast<MatrixD*>(blocks.at(1).get());
    rect_matrices[1] = reinterpret_cast<MatrixD*>(blocks.at(2).get());
#ifdef __3D__
    Assemble3D_mixed(
#else
    Assemble2D_VectFE(
#endif
                      n_fe_spaces, fespmat, n_sq_mat, sq_matrices, n_rect_mat,
                      rect_matrices, n_rhs, RHSs, fespmat, la,
                      boundary_conditions, non_const_bound_values.data());

    
    Assemble_DG<d>(la.GetCoeffFct(), 1, &v_space, 1, sq_matrices,
                0, nullptr, 1, RHSs, boundary_conditions,
                &non_const_bound_values[2], db, "NSE");
//     s.rhs.print("rhs");
//     blocks.at(1)->PrintFull("BT");
//     blocks.at(2)->PrintFull("B");
  }
  // copy Dirichlet values from rhs to solution vector (this is not really
  // necessary in case of a direct solver)
  this->systems.front().solution.copy_nonactive(this->systems.front().rhs);
} // void NavierStokes_DG3D::Assemble

//added 21.02.22
/* ************************************************************************* */
template <int d>
void NavierStokes_DG<d>::assemble_nonlinear_term()
{
  using MatrixD = typename Template_names<d>::MatrixD;
  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  
  using BoundaryValuesFunction
    = typename Template_names<d>::BoundaryValuesFunction;
  using BoundaryConditionFunction 
    = typename Template_names<d>::BoundaryConditionFunction;
  
  const size_t n_fe_spaces = 2; // space needed for assembling matrices
 
  const size_t n_sq_mat = 2; // no of square matrices
  //std::vector<SquareMatrixD*> sq_matrices(n_square_matrices, nullptr);
 
  const size_t n_rect_mat = 2; // no of rectangular matrices
  //std::vector<MatrixD*> re_matrices(n_rectangular_matrices, nullptr);
 
  const size_t n_rhs = 2; // number of right hand sides
  // std::vector<double*> rhs_array(n_rhs, nullptr); // right hand side 
  
 
  // boundary conditions and boundary values
 // std::array<BoundaryConditionFunction*, n_fe_spaces+1> boundCondition;
  //std::array<BoundaryValuesFunction*, n_fe_spaces+1> non_const_bound_values;
 // for(int i = 0; i < d+1; ++i)
 //   boundValues[i] = example.get_bd()[i];
  
  //std::array<const FESpace*, d+1> rhs_spaces;
  
  //Nonlinear assembling requires an approximate velocity solution on every grid!
  if(systems.size() > 1)
  {
    for( int block = 0; block < d ;++block)
    {
      std::vector<const FESpace*> spaces;
      std::vector<double*> u_entries;
      std::vector<size_t> u_ns_dofs;
      for(auto &s : systems )
      {
        spaces.push_back(s.velocity_space.get());
        u_entries.push_back(s.solution.block(block));
        u_ns_dofs.push_back(s.solution.length(block));
      }
      GridTransfer::RestrictFunctionRepeatedly(spaces, u_entries, u_ns_dofs);
    }
  }
  
  bool mdml =  this->solver.is_using_multigrid() 
            && this->solver.get_multigrid()->is_using_mdml();
            
  bool is_stokes = this->db["problem_type"].is(3); // otherwise Navier-Stokes
  if ((mdml && !is_stokes)|| db["space_discretization_type"].is("upwind"))
  {
    // in case of upwinding we only assemble the linear terms. The nonlinear
    // term is not assembled but replaced by a call to the upwind method.
    // Note that we assemble the same terms over and over again here. Not 
    // nice, but otherwise we would have to store the linear parts in a 
    // separate BlockFEMatrix.
    this->assemble();
  }

  for(auto &s : this->systems)
  {
    // the class LocalAssembling3D which we will need next, requires an array of
    // pointers to finite element functions, i.e. TFEFunction3D **.
    std::vector<const FEFunction*> fe_functions{&s.u, &s.p};
    
    //hold the velocity space, we'll need it...
    const FESpace * v_space = s.velocity_space.get();
    const FESpace * p_space = s.pressure_space.get();
    

    // spaces for matrices
    const FESpace *fespmat[2] = {v_space, p_space};
    
    SquareMatrixD *sq_matrices[n_sq_mat]{};
    MatrixD *rect_matrices[n_rect_mat]{};

    double *RHSs[n_rhs] = {s.rhs.block(0), s.rhs.block(1)};
    
    
   /* for(int i = 0; i < d; ++i)
      rhs_spaces[i] = v_space;
    rhs_spaces[d] = p_space;*/
    
    BoundaryConditionFunction * boundary_conditions[n_fe_spaces] = {
      v_space->get_boundary_condition(), p_space->get_boundary_condition() };
    std::array<BoundaryValuesFunction*, n_fe_spaces+1> non_const_bound_values;


    non_const_bound_values[0] = example.get_bd()[3];  // u * n
    non_const_bound_values[1] = example.get_bd()[2];  // p
    non_const_bound_values[2] = example.get_bd()[4];  // u * t

    std::vector<std::shared_ptr<FEMatrix>> blocks
      = s.matrix.get_blocks_uniquely();
    sq_matrices[0] = reinterpret_cast<SquareMatrixD*>(blocks.at(0).get());
    sq_matrices[1] = reinterpret_cast<SquareMatrixD*>(blocks.at(3).get());
    
    rect_matrices[0] = reinterpret_cast<MatrixD*>(blocks.at(1).get());
    rect_matrices[1] = reinterpret_cast<MatrixD*>(blocks.at(2).get());
 
    
    // What is this?
    /*std::vector<std::vector<size_t>> cells;
    for(size_t i = 0; i < d; ++i)
    {
      for(size_t j = 0; j < d; ++j)
        cells.push_back({{i, j}});
    }
    auto blocks = s.matrix.get_blocks_uniquely(cells);
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
      case 2:
        sq_matrices[0] = reinterpret_cast<SquareMatrixD*>(blocks[0].get());
        break;
      case 3: 
      case 4:
      case 14:
        for(int i = 0; i < d*d; ++i)
          sq_matrices[i] = reinterpret_cast<SquareMatrixD*>(blocks[i].get());
        if(db["space_discretization_type"].is("supg"))
        {
          auto blocks = s.matrix.get_blocks_uniquely();
          for(int i = 0; i < d; ++i)
            re_matrices[d+i] = reinterpret_cast<MatrixD*>(blocks[d+(d+1)*i].get());
          for(auto rm : re_matrices)
          {
            if(rm != nullptr)
              rm->reset();
          }
          for(int i = 0; i < d+1; ++i)
            rhs_array[i] = s.rhs.block(i);
          s.rhs.reset();
        }
        break;
    }// endswitch nstype

    for(int i = 0; i < d; ++i)
      boundCondition[i] = spaces[0]->get_boundary_condition();
    boundCondition[d] = spaces[1]->get_boundary_condition();

    // finite element functions
    std::array<std::unique_ptr<FEFunction>, d+1> velocity_components;
    for(int i = 0; i < d; ++i)
    {
      velocity_components[i] = s.u.GetComponent(i);
      feFunction[i] = velocity_components[i].get();
      s.velocity_space.get()
    }
    feFunction[d] = &s.p;
    */

    //decide wether to assemble by upwinding or not
    bool finest_grid = (&s == &systems.at(0));
    bool do_upwinding = (db["space_discretization_type"].is("upwind")
                        || (mdml && !finest_grid))
                        && !is_stokes;
    // local assembling object    
    LocalAssembling<d> la(this->db, LocalAssembling_type::NavierStokesNL, 
                          fe_functions, this->example.get_coeffs());
    if(!do_upwinding)
    {
      for(auto* mat : sq_matrices)
      {
        if(mat != nullptr)
          mat->reset();
      }

      // assemble now the matrices and right hand side 
#ifdef __3D__
      Assemble3D(
#else
      Assemble2D_VectFE(
#endif
                 n_fe_spaces, fespmat, n_sq_mat, sq_matrices,
                 n_rect_mat, rect_matrices, n_rhs,
                 RHSs, fespmat, la, boundary_conditions,
                 non_const_bound_values.data());
     
     Assemble_DG<d>(la.GetCoeffFct(), 1, &v_space, 1, sq_matrices,
                0, nullptr, 1, RHSs, boundary_conditions,
                &non_const_bound_values[2], db, "NSE");
    }
    else
    {
#ifdef __3D__
      // the inverse of the example's diffusion coefficient
      double one_over_nu = 1./example.get_nu();
      for(auto mat : sq_matrices)
      {
        UpwindForNavierStokes3D(mat, feFunction[0], feFunction[1],
                                feFunction[2], one_over_nu);
      }
#else // 3D -> 2D
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
        case 2:
          // do upwinding with one matrix
          UpwindForNavierStokes(la.GetCoeffFct(), sq_matrices[0],
                                la.get_fe_function(0), la.get_fe_function(1));
          break;

        case 3:
        case 4:
        case 14:
          // do upwinding with two matrices
          UpwindForNavierStokes(la.GetCoeffFct(), sq_matrices[0],
                                la.get_fe_function(0), la.get_fe_function(1));
          UpwindForNavierStokes(la.GetCoeffFct(), sq_matrices[1],
                                la.get_fe_function(0), la.get_fe_function(1));
          break;
      } // endswitch
#endif
    }
    
    }
  }// endfor auto grid*/



/* ************************************************************************* */
/*
// added 22.02.22
template <int d>
bool NavierStokes_DG<d>::stop_it(unsigned int iteration_counter)
{
  int my_rank = 0;
#ifdef _MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
  //compute and update defect and residuals
  compute_residuals();
  
  // the current norm of the residual
  const double normOfResidual = this->get_full_residual();
  // store initial residual, so later we can print the overall reduction
  if(iteration_counter == 0)
  {
    initial_residual = normOfResidual;
#ifdef _MPI
    auto comms = this->get_matrix().get_communicators();
    initial_rhs_norm = this->systems.front().rhs.norm(comms);
#else
    initial_rhs_norm = this->systems.front().rhs.norm();
#endif
  }
  
  // check if minimum number of iterations was performed already
  size_t min_it = db["nonlinloop_minit"];
  if(iteration_counter < min_it)
	  return false;

  // hold the residual from 10 iterations ago
  const double oldNormOfResidual = this->old_residuals.front().fullResidual;

  size_t max_it = db["nonlinloop_maxit"];
  double convergence_speed = db["nonlinloop_slowfactor"];
  bool slow_conv = false;

  if(normOfResidual >= convergence_speed*oldNormOfResidual)
    slow_conv = true;

  double limit = db["nonlinloop_epsilon"];
  if(db["nonlinloop_scale_epsilon_with_size"])
  {
    limit *= std::sqrt(this->get_size());
    if(my_rank==0)
     Output::print<1>("stopping tolerance for nonlinear iteration ", limit);
  }
  
  //check residual relative to initial right hand side
  if(db["nonlinloop_residual_relative_to_rhs"])
    limit *= initial_rhs_norm;

  // check if the iteration has converged, or reached the maximum number of
  // iterations or if convergence is too slow. Then return true otherwise false
  if((normOfResidual <= limit) || (iteration_counter == max_it) || (slow_conv))
  {
    if(slow_conv && my_rank==0)
      Output::print<1>(" SLOW !!! ", normOfResidual/oldNormOfResidual);

    // stop iteration
    adjust_pressure();
    return true;
  }
  else
    return false;
}

*/


/** ************************************************************************ */
template <int d>
void NavierStokes_DG<d>::solve()
{
  double t = GetTime();
  System_per_grid& s = this->systems.front();

  this->solver.solve(s.matrix, s.rhs, s.solution);

  if(s.matrix.pressure_projection_enabled())
    s.p.project_into_L20();
  
  t = GetTime() - t;
  Output::print<2>(" solving a NavierStokes_DG", d, "D problem done in ", t, " seconds");
}

/** ************************************************************************ */
template <int d>
void NavierStokes_DG<d>::output(int i)
{
  using ErrorMethod = typename FEFunction::ErrorMethod;

  System_per_grid & s = this->systems.front();
  if((size_t)db["verbosity"]> 1)
  {
    //s.u.PrintMinMax();
    //s.p.PrintMinMax();
  }

  if(i < 0)
    outputWriter.write();
  else
    outputWriter.write(i);
  checkpoint_io.write(this->get_solution());
  if(db["output_compute_errors"])
  {
    auto *const *Exact = &(example.get_exact())[0];
    ErrorMethod *L2DivH1 = L2DivH1Errors;
    // unfortunatly this needs to be longer than one would expect. This is
    // because TFEFunction3D::GetErrors needs this array to be of size
    // N_Errors+1 (here 2+1).
    std::array<double, 6> errors;
    s.u.GetErrorsForVectorValuedFunction(Exact, L2DivH1, errors.data());

#ifdef __3D__
    TAuxParam3D aux;
    MultiIndex3D AllDerivatives[4] = { MultiIndex3D::D000, MultiIndex3D::D100,
                                       MultiIndex3D::D010, MultiIndex3D::D001 };
#else
    TAuxParam2D aux;
    MultiIndex2D AllDerivatives[3] = { MultiIndex2D::D00, MultiIndex2D::D10,
                                       MultiIndex2D::D01 };
#endif
    const FESpace * pointer_to_p_space = s.pressure_space.get();
    // the 'd' for the number of errors is a dirty hack, we need to refactor
    // the GetErrors methods.
    s.p.GetErrors(example.get_exact(d), d+1, AllDerivatives, d, L2H1Errors,
                  nullptr, &aux, 1, &pointer_to_p_space, errors.data() + 3);

    Output::print<1>(" L2(u):      ", setprecision(14), errors[0]);
    Output::print<1>(" L2(div(u)): ", setprecision(14), errors[1]);
    Output::print<1>(" H1-semi(u): ", setprecision(14), errors[2]);
    Output::print<1>(" L2(p):      ", setprecision(14), errors[3]);
    Output::print<1>(" H1-semi(p): ", setprecision(14), errors[4]);

    // copy
    std::copy(errors.begin(), errors.end()-1, this->errors.begin());
  } // if(TDatabase::ParamDB->MEASURE_ERRORS)
  
//    s.solution.print("sol_computed");
  BlockVector r(s.rhs);
//   s.matrix.apply_scaled_add(s.solution, r, -1.0);
//   r.print("res");
  
  std::vector<double> val_u(s.u.GetLength());
  auto values = s.u.GetValues();
  for (size_t i = 0; i < val_u.size(); ++i)
  {
    val_u[i] = values[i];
  }
  std::vector<double> val_p(s.p.GetLength());
  values = s.p.GetValues();
  for (size_t i = 0; i < val_p.size(); ++i)
  {
    val_p[i] = values[i];
  }

  s.u.Interpolate(example.get_exact()[0], example.get_exact()[1]);
  s.p.Interpolate(example.get_exact()[2]);
//   s.solution.print("sol_exact");
  r = s.rhs;
  s.matrix.apply_scaled_add(s.solution, r, -1.0);
  //r.print("res");

  for (size_t i = 0; i < val_p.size(); ++i)
  {
    values[i] = val_p[i];
  }
  values = s.u.GetValues();
  for (size_t i = 0; i < val_u.size(); ++i)
  {
    values[i] = val_u[i];
  }

  // Output::print(">> norm of residual: ", r.norm());
}

/** ************************************************************************ */
template <int d>
const BlockVector & NavierStokes_DG<d>::get_solution() const
{
  return this->systems.front().solution;
}

/** ************************************************************************ */
template <int d>
BlockVector & NavierStokes_DG<d>::get_solution()
{
  return this->systems.front().solution;
}

/** ************************************************************************ */
template <int d>
double NavierStokes_DG<d>::getL2VelocityError() const
{
  return this->errors[0];
}

/** ************************************************************************ */
template <int d>
double NavierStokes_DG<d>::getL2DivergenceError() const
{
  return this->errors[1];
}

/** ************************************************************************ */
template <int d>
double NavierStokes_DG<d>::getH1SemiVelocityError() const
{
  return this->errors[2];
}

/** ************************************************************************ */
template <int d>
double NavierStokes_DG<d>::getL2PressureError() const
{
  return this->errors[3];
}

/** ************************************************************************ */
template <int d>
double NavierStokes_DG<d>::getH1SemiPressureError() const
{
  return this->errors[4];
}

#ifdef __3D__
template class NavierStokes_DG<3>;
#else
template class NavierStokes_DG<2>;
#endif



