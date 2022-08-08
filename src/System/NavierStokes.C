#include "NavierStokes.h"
#include "LocalAssembling.h"
#include "Database.h"
#include "Multigrid.h"
#include "GridTransfer.h"
#include "NSE_local_assembling_routines.h"
#include "LPS_scott_zhang.h"
#include "BaseCell.h"
#ifdef __3D__
#include "Assemble3D.h"
#include "SquareMatrix3D.h"
#include "Upwind3D.h"
#include "AuxParam3D.h"
#include "BoundaryAssembling3D.h"
#include "Matrix3D.h"
#else
#include "Assemble2D.h"
#include "SquareMatrix2D.h"
#include "Upwind.h"
#include "AuxParam2D.h"
#include "BoundaryAssembling2D.h"
#include "Matrix2D.h"
#endif
#ifdef _MPI
#include "ParFECommunicator3D.h"
#include "MumpsWrapper.h"
#endif
#include "Saddle_point_preconditioner.h"
#include <cmath>

#include <AddGadientJump.hpp>
#include "Assemble_DG.h"

template <int d>
ParameterDatabase NavierStokes<d>::default_nse_database(bool complete)
{
  Output::print<5>("creating a default NSE parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default NSE3D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("NSE parameter database");
  
  if(complete)
  {
    db.merge(TDomain::default_domain_parameters(), true);
    if(d == 3)
      db.add_nested_database(TDomain::default_sandwich_grid_parameters());
    db.merge(Example_NSE::default_example_database());
  }
  //NSE3D requires a nonlinear iteration, set up a nonlinit_database and merge
  db.merge(ParameterDatabase::default_nonlinit_database(), true);
  // a default output database - needed here as long as there's no class handling the output
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
  db["problem_type"] = 5; // means stationary Navier-Stokes
  return db;
}

/* ************************************************************************* */
template <int d>
NavierStokes<d>::System_per_grid::System_per_grid (
  const Example_NSE& example, const TCollection& coll,
  std::pair<int,int> velocity_pressure_orders)
 : velocity_space(new FESpace(&coll, "u", example.get_bc(0), 
                              velocity_pressure_orders.first)),
   pressure_space(new FESpace(&coll, "p", example.get_bc(d),
                              velocity_pressure_orders.second))
{
  // build the matrix due to NSE type
  switch (TDatabase::ParamDB->NSTYPE)
    {
#ifdef __2D__
    case 1:
      matrix = BlockFEMatrix::NSE2D_Type1(velocity_space, pressure_space);
      break;
    case 2:
      matrix = BlockFEMatrix::NSE2D_Type2(velocity_space, pressure_space);
      break;
    case 3:
      matrix = BlockFEMatrix::NSE2D_Type3(velocity_space, pressure_space);
      break;
    case 4:
      matrix = BlockFEMatrix::NSE2D_Type4(velocity_space, pressure_space);
      break;
    case 14:
      matrix = BlockFEMatrix::NSE2D_Type14(velocity_space, pressure_space);
      break;
#else // 3D
    case 1:
      matrix = BlockFEMatrix::NSE3D_Type1(velocity_space, pressure_space);
      break;                                               
    case 2:                                                
      matrix = BlockFEMatrix::NSE3D_Type2(velocity_space, pressure_space);
      break;                                               
    case 3:                                                
      matrix = BlockFEMatrix::NSE3D_Type3(velocity_space, pressure_space);
      break;                                               
    case 4:                                                
      matrix = BlockFEMatrix::NSE3D_Type4(velocity_space, pressure_space);
      break;                                               
    case 14:                                               
      matrix = BlockFEMatrix::NSE3D_Type14(velocity_space, pressure_space);
      break;
#endif
    default:
      ErrThrow("NSTYPE: ", TDatabase::ParamDB->NSTYPE, " is not known");
  }
  rhs = BlockVector(matrix, true);
  solution = BlockVector(matrix, false);

  u = FEVectFunct(velocity_space, "u", solution.block(0), d);
  p = FEFunction(pressure_space, "p", solution.block(d));
  
#ifdef _MPI
  //print some information
  velocity_space->get_communicator().print_info();
  pressure_space->get_communicator().print_info();
#endif
}

/* ************************************************************************* */
template <int d>
NavierStokes<d>::System_per_grid::System_per_grid(const System_per_grid& other)
 : velocity_space(other.velocity_space), pressure_space(other.pressure_space),
   matrix(other.matrix), rhs(other.rhs), solution(other.solution)
{
  // the fe functions must be newly created, because copying would mean 
  // referencing the BlockVectors in 'other'.
  u = FEVectFunct(velocity_space, "u", solution.block(0), d);
  p = FEFunction(pressure_space, "p", solution.block(d));
}

/* ************************************************************************* */
template <int d>
NavierStokes<d>::NavierStokes(const TDomain& domain,
                              const ParameterDatabase& param_db)
 : NavierStokes<d>(domain, param_db, Example_NSE(param_db))
{
}

/* ************************************************************************* */
template <int d>
NavierStokes<d>::NavierStokes(const TDomain& domain,
                              const ParameterDatabase& param_db,
                              Example_NSE example_in)
  : systems(), example(example_in), db(default_nse_database()),
    outputWriter(param_db), checkpoint_io(param_db), solver(param_db), defect(),
    old_residuals(), initial_residual(1e10), errors()
{
  this->db.merge(param_db, false);
  auto collections = domain.get_grid_collections();
  const TCollection *coll = collections.front(); //the finest grid collection
  this->check_parameters(coll);

  std::pair <int,int> 
      velocity_pressure_orders(TDatabase::ParamDB->VELOCITY_SPACE, 
                               TDatabase::ParamDB->PRESSURE_SPACE);
  // set the velocity and pressure spaces
  // this function returns a pair which consists of 
  // velocity and pressure order
  this->get_velocity_pressure_orders(velocity_pressure_orders);

  bool usingMultigrid = solver.is_using_multigrid();
  // create finite element space and function, a matrix, rhs, and solution
  systems.emplace_back(example, *coll, velocity_pressure_orders);
  
  if(usingMultigrid)
  {
    // Construct multigrid object
    auto mg = solver.get_multigrid();
    bool mdml = mg->is_using_mdml();

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

    if(mdml)
    {
      // change the discretization on the coarse grids to lowest order 
      // non-conforming(-1). The pressure space is chosen automatically(-4711).
      velocity_pressure_orders = {-1, -4711};
      this->get_velocity_pressure_orders(velocity_pressure_orders);
    }
    else
    {
      // for standard multigrid, pop the finest collection - it was already
      // used to construct a space before the "if(usingMultigrid)" clause
      // and will not (as in mdml) be used a second time with a different discretization
      collections.pop_front();
    }
    
     // Construct systems per grid and store them, finest level first
    std::list<BlockFEMatrix*> matrices;
    // matrix on finest grid is already constructed
    matrices.push_back(&systems.back().matrix);

    for(auto coll : collections) // initialize the coarse grid space hierarchy
    {
      systems.emplace_back(example, *coll, velocity_pressure_orders);
      // prepare input argument for multigrid object
      matrices.push_front(&systems.back().matrix);
    }
    // initialize the multigrid object with all the matrices on all levels
    mg->initialize(matrices);
  }
  
  outputWriter.add_fe_vector_function(&this->get_velocity());
  outputWriter.add_fe_function(&this->get_pressure());
  checkpoint_io.read(this->get_solution());

  if(db["output_write_exact_solution"])
  {
    // initialize variables for storing exact solution (FE functions)
    solution_exact = BlockVector(this->get_solution());
    u_exact = FEVectFunct(this->systems.front().velocity_space, "u_exact",
                          solution_exact.block(0), d);
    p_exact = FEFunction(this->systems.front().pressure_space, "p_exact",
                         solution_exact.block(d));
    // Interpolating the exact solution
    for(int i = 0; i < d; ++i)
    {
      auto ui = u_exact.GetComponent(i);
      ui->Interpolate(example.get_exact(i));
    }
    p_exact.Interpolate(example.get_exact(d));
    p_exact.PrintMinMax(std::string("p_exact"));
    outputWriter.add_fe_vector_function(&u_exact);
    outputWriter.add_fe_function(&p_exact);
  }
  
  output_problem_size_info();
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::check_parameters(const TCollection* collection)
{
  if(!db["problem_type"].is(3) && !db["problem_type"].is(5) &&
     !db["problem_type"].is(7))
  {
    Output::warn<2>("The parameter problem_type doesn't correspond neither to "
        "NSE nor to Stokes or Brinkman. It is now reset to the default value "
        "for NSE (=5).");
    db["problem_type"] = 5;
  }
  bool mat_c_enabled = (TDatabase::ParamDB->NSTYPE == 14);
  bool has_hanging_nodes = collection->includes_hanging_vertices();
  if(has_hanging_nodes && !mat_c_enabled)
  {
    ErrThrow("In order for hanging nodes to work, the NSTYPE hast to be 14. "
             "Currently it is set to ", TDatabase::ParamDB->NSTYPE, ".");
  }
  // more check needed
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::get_velocity_pressure_orders(
  std::pair<int, int>& velocity_pressure_orders)
{
  int velocity_order = velocity_pressure_orders.first;
  int pressure_order = velocity_pressure_orders.second;
  int order = 0;
  switch(velocity_order)
  {
    case 1: case 2: case 3: case 4: case 5: // P_k/Q_k
      order = velocity_order;
      break;
    case 12: case 13: case 14: case 15:
      if(d == 3 && velocity_order == 12)
      {
        ErrThrow("Scott-Vogelius P2/P1disc is not stable in 3D even if the final "
                 "refinement step is barycentric.");
      }
      // P2/P1disc and P3/P2disc elements are in general not inf-sup stable on 
      // triangles. If the last refinement step was barycentric refinement,
      // then these elements are inf-sup stable.  
      Output::warn("You chose to use Scott-Vogelius finite elements, make sure "
                   "that the final refinement step is barycentric, see the "
                   "parameter 'refinement_final_step_barycentric'.");
      order = velocity_order-10;
      break;
    case -1: case -2: case -3: case -4: case -5:
    case -101:
      order = velocity_order;
      break;
    case 101:
      order = velocity_order;
    // conforming fe spaces with bubbles on triangles
    case 22: case 23: case 24: case 25:
      order = velocity_order;
      if(d == 3 && velocity_order == 25)
      {
        ErrThrow("Velocity space 25 not supported in 3D");
      }
      break;
      // discontinuous spaces 
    case -11: case -12: case -13:
      order = velocity_order*10;
      break;
    default:
      ErrThrow("velocity space order ", velocity_order, " is not supported.");
  }
  TDatabase::ParamDB->VELOCITY_SPACE = order;
  velocity_pressure_orders.first = order;
  switch(pressure_order)
  {
    case -4711:
      switch(velocity_order)
      {
        case -1:
        case -2:
        case -3:
        case -4:
          // nonconforming pw (bi)linear velo/ pw constant pressure
          // conforming pw (bi)linear velo/ pw constant pressure (not stable !!!)
          pressure_order = -velocity_order-1;
          break; 
        case 1: // discontinuous space 
          pressure_order = 0;
          Output::warn("NSE3D", "The P1/P0 element pair (Q1/Q0 on hexa) is "
                       " not stable. Make sure to use stabilization!");
          break;
        case 2: case 3: case 4: case 5:
        // standard conforming velo and continuous pressure
          pressure_order = velocity_order-1;
          break;
          // Scott-Vogelius: discontinuous pressure spaces with standard 
          // conforming velocity space. This is not stable on general triangles,
          // be sure to use barycentric refinement.
        case 12: case 13: case 14: case 15:
          pressure_order = -velocity_order+1;
          break;
        case 22: case 23: case 24: case 25:
          pressure_order = -(velocity_order-11)*10;
          break;
        case 101:
          pressure_order = 1;
      }
      break;
    // continuous pressure spaces
    case 1: case 2: case 3: case 4: case 5:
      // nothing to do
      break;
    // discontinuous spaces
    case -11: case -12: case -13: case -14:
      pressure_order = pressure_order*10;
      break;
  }
  TDatabase::ParamDB->PRESSURE_SPACE  = pressure_order;
  velocity_pressure_orders.second = pressure_order;
  
  Output::print("velocity space", setw(10), TDatabase::ParamDB->VELOCITY_SPACE);
  Output::print("pressure space", setw(10), TDatabase::ParamDB->PRESSURE_SPACE);
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::output_problem_size_info() const
{
  int my_rank = 0;
#ifndef _MPI
  auto & velocity_space = *this->systems.front().velocity_space;
  auto & pressure_space = *this->systems.front().pressure_space;

  size_t nDofu  = velocity_space.get_n_dof();
  size_t nDofp  = pressure_space.get_n_dof();
  size_t nTotal = d*nDofu + nDofp;
  size_t nActive= d*velocity_space.get_n_active();

  auto coll = velocity_space.GetCollection();

  double hmin, hmax;
  coll->GetHminHmax(&hmin, &hmax);
#else
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  auto velocity_comm = systems.front().velocity_space->get_communicator();
  auto pressure_comm = systems.front().pressure_space->get_communicator();
  int nDofu  = velocity_comm.get_n_global_dof();
  int nDofp  = pressure_comm.get_n_global_dof();
  int nTotal = d*nDofu + nDofp;
#endif
  if(my_rank ==0)
  {
    Output::stat("NavierStokes", "Mesh data and problem size on finest grid");
#ifndef _MPI
    Output::dash("N_Cells            :  ", setw(10), coll->GetN_Cells());
    Output::dash("h(min, max)        :  ", setw(10), hmin, setw(10), " ", hmax);
#endif
    Output::dash("dof velocity       :  ", setw(10), d*nDofu);
#ifndef _MPI
    Output::dash("dof velocity active:  ", setw(10), nActive);
#endif
    Output::dash("dof pressure       :  ", setw(10), nDofp);
    Output::dash("dof total          :  ", setw(10), nTotal);
  }
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::assemble_linear_terms(
  std::function<void(LocalAssembling<d>& la)> modify_la)
{
  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  using MatrixD = typename Template_names<d>::MatrixD;
  using BoundaryValuesFunction
    = typename Template_names<d>::BoundaryValuesFunction;
  using BoundaryConditionFunction 
    = typename Template_names<d>::BoundaryConditionFunction;
  
  int n_fe_spaces = 2; // spaces used for assembling matrices
  int n_square_matrices=d*d + 1; // maximum no of square matrices (type 14)
  std::vector<SquareMatrixD*> sq_matrices(n_square_matrices, nullptr);
  int n_rectangular_matrices = 2*d; // maximum no of rectangular matrices
  std::vector<MatrixD*> re_matrices(n_rectangular_matrices, nullptr);
  constexpr int n_rhs = d+1; // maximum number of right hand sides
  std::array<double*, n_rhs> rhs_array; // right hand side 
  // finite element function used for nonlinear term
  std::vector<const FEFunction*> feFunction(d+1, nullptr);
  // boundary conditions and boundary values
  std::array<BoundaryConditionFunction*, d+1> boundCondition;
  std::array<BoundaryValuesFunction*, d+1> boundValues;
  for(int i = 0; i < d+1; ++i)
    boundValues[i] = example.get_bd()[i];
  
  std::array<const FESpace*, d+1> rhs_spaces;
  
  for(auto &s : this->systems)
  {
    s.rhs.reset(); // right hand side reset (is that necessary?)
    s.matrix.reset(); // reset matrix (needed for mdml where this is called)
    
    const FESpace *v_space = s.velocity_space.get();
    const FESpace *p_space = s.pressure_space.get();

    // spaces for matrices
    const FESpace *spaces[2] = {v_space, p_space};
    for(int i = 0; i < d; ++i)
      rhs_spaces[i] = v_space;
    rhs_spaces[d] = p_space;

    // spaces for right hand side
    for(int i = 0; i <= d; ++i)
      rhs_array[i] = s.rhs.block(i);
    
    std::vector<std::shared_ptr<FEMatrix>> blocks 
      = s.matrix.get_blocks_uniquely();
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        sq_matrices[0] = reinterpret_cast<SquareMatrixD*>(blocks[0].get());
        for(int i = 0; i < d; ++i)
          re_matrices[i] = reinterpret_cast<MatrixD*>(blocks[1+i].get());
        break;
      case 2:
        sq_matrices[0] = reinterpret_cast<SquareMatrixD*>(blocks[0].get());
        for(int i = 0; i < d; ++i)
          re_matrices[i] = reinterpret_cast<MatrixD*>(blocks[d+1+i].get());
        for(int i = 0; i < d; ++i)
          re_matrices[i+d] = reinterpret_cast<MatrixD*>(blocks[1+i].get());
        break;
      case 3:
        for(int i = 0, j = 0; i < d*d; ++i, ++j)
        {
          if(i%d == 0 && i > 0)
            j++;
          sq_matrices[i] = reinterpret_cast<SquareMatrixD*>(blocks[j].get());
        }
        for(int i = 0; i < d; ++i)
          re_matrices[i] = reinterpret_cast<MatrixD*>(blocks[d+(d+1)*i].get());
        break;
      case 4:
        for(int i = 0, j = 0; i < d*d; ++i, ++j)
        {
          if(i%d == 0 && i > 0)
            j++;
          sq_matrices[i] = reinterpret_cast<SquareMatrixD*>(blocks[j].get());
        }
        for(int i = 0; i < d; ++i)
          re_matrices[i] = reinterpret_cast<MatrixD*>(blocks[d*(d+1)+i].get());
        for(int i = 0; i < d; ++i)
          re_matrices[d+i] = reinterpret_cast<MatrixD*>(blocks[d+(d+1)*i].get());
        break;
      case 14:        
        for(int i = 0, j = 0; i < d*d; ++i, ++j)
        {
          if(i%d == 0 && i > 0)
            j++;
          sq_matrices[i] = reinterpret_cast<SquareMatrixD*>(blocks[j].get());
        }
        sq_matrices[d*d] = reinterpret_cast<SquareMatrixD*>(blocks[(d+1)*(d+1)-1].get());
        for(int i = 0; i < d; ++i)
          re_matrices[i] = reinterpret_cast<MatrixD*>(blocks[d*(d+1)+i].get());
        for(int i = 0; i < d; ++i)
          re_matrices[d+i] = reinterpret_cast<MatrixD*>(blocks[d+(d+1)*i].get());
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
    }
    feFunction[d] = &s.p;

    // local assembling object    
    LocalAssembling<d> la(this->db, LocalAssembling_type::NavierStokesAll, 
                          feFunction, example.get_coeffs());
    
    if(modify_la)
      modify_la(la);
    
    // assemble now the matrices and right hand side 
#ifdef __3D__
    Assemble3D(
#else
    Assemble2D(
#endif
               n_fe_spaces, spaces, 
               n_square_matrices, sq_matrices.data(),
               n_rectangular_matrices, re_matrices.data(), 
               n_rhs, rhs_array.data(), rhs_spaces.data(),
               boundCondition.data(), boundValues.data(), la);    
#ifdef __2D__
//     AddGadientJump(n_fe_spaces, spaces, 
//                     n_square_matrices, sq_matrices.data(), 
//                     0, nullptr, n_rhs,
//                     nullptr, nullptr, 
//                     boundCondition.data(), 
//                     boundValues.data(), la);
//     Assemble_DG<2>(la.GetCoeffFct(), 1, spaces, 1, &sq_matrices[0], 0, nullptr, 0, nullptr, boundCondition.data(), boundValues.data(), db, "Oseen", 1, &sq_matrices[3]);
    
    AddJumpStabilizationCIP(1, spaces, 1, sq_matrices.data(), 
            boundCondition.data(), boundValues.data(), db, la);
    // sq_matrices[2]->PrintFull("A12");exit(0);
//     Assemble_DG<2>(la.GetCoeffFct(), 1, spaces, 1, &sq_matrices[3], 0, nullptr, 0, nullptr, boundCondition.data(), boundValues.data(), db, "Oseen", 3);
//     TAuxParam2D aux;
//     Assemble2D_CIP(la.GetCoeffFct(), 1, spaces, 1, &sq_matrices[0], 0, nullptr, 0, nullptr, 
//                    nullptr,boundCondition.data(), boundValues.data(), &aux);
//     ComputeAddStab(spaces, la, &sq_matrices[0], &sq_matrices[3]);
//     sq_matrices[0]->Print();
// exit(0);
#endif
    // assemble on the boundary if needed
    assemble_boundary_terms();
    
    // copy Dirichlet values from right hand side into solution
    s.solution.copy_nonactive(s.rhs);
    
    if(db["space_discretization_type"].is("local_projection"))
    {
      if(d == 3)
        ErrThrow("local_projection stabilization is not implemented in 3D");
      LPS_parameter_set lp{db["lps_coeff_type"], db["lps_delta0"],
                           db["lps_delta1"]};
      auto C = LPS_for_pressure_Scott_Zhang(blocks.at(8), false,
                                            this->example.get_nu(), lp);
      s.matrix.replace_blocks(*C, {{d,d}}, {false}); // creates a copy
    }
    
    if(v_space->get_n_hanging() != 0 || p_space->get_n_hanging() != 0)
    {
      if(TDatabase::ParamDB->NSTYPE)
      {
        // remove entries in hanging and Dirichlet rows in non-diagonal blocks
        for(int i = 0; i < d+1; ++i)
        {
          for(int j = 0; j < d+1; ++j)
          {
            if(i == j)
              continue; // skip diagonal blocks
            blocks[(d+1)*i+j]->resetNonActive();
          }
        }
      }
      else
      {
        ErrThrow("hanging nodes and NSTYPE ", TDatabase::ParamDB->NSTYPE,
                 " is not supported.");
      }
    }
  }// endfor auto grid

/** When we call copy_nonactive in MPI-case, we have to remember the following:
   * it can happen that some slave ACTTIVE DoFs are placed in the block of
   * NON-ACTIVE DoFs (because they are at the interface between processors).
   * Doing copy_nonactive changes then the value of these DOFs,although they are
   * actually active.
   * That's why we have to update the values so that the vector becomes consistent again.
   * This is done here.
   */
#ifdef _MPI
  auto& s = this->systems.front();
  for(int i = 0; i < d; ++i)
  {
    auto ui = s.solution.block(i);
    s.velocity_space->get_communicator().consistency_update(ui, 3);
  }
  auto p = s.solution.block(3);
  s.pressure_space->get_communicator().consistency_update(p, 3);
#endif
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::assemble_nonlinear_term()
{
  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  using MatrixD = typename Template_names<d>::MatrixD;
  using BoundaryValuesFunction
    = typename Template_names<d>::BoundaryValuesFunction;
  using BoundaryConditionFunction 
    = typename Template_names<d>::BoundaryConditionFunction;
  
  size_t n_fe_spaces = 2; // space needed for assembling matrices
  size_t n_square_matrices = d*d + 1; // no of square matrices
  std::vector<SquareMatrixD*> sq_matrices(n_square_matrices, nullptr);
  int n_rectangular_matrices = 2*d; // maximum no of rectangular matrices
  std::vector<MatrixD*> re_matrices(n_rectangular_matrices, nullptr);
  int n_rhs = d+1; // maximum number of right hand sides
  std::vector<double*> rhs_array(n_rhs, nullptr); // right hand side 
  
  std::vector<const FEFunction*> feFunction(d+1, nullptr);
  // boundary conditions and boundary values
  std::array<BoundaryConditionFunction*, d+1> boundCondition;
  std::array<BoundaryValuesFunction*, d+1> boundValues;
  for(int i = 0; i < d+1; ++i)
    boundValues[i] = example.get_bd()[i];
  
  std::array<const FESpace*, d+1> rhs_spaces;
  
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
    this->assemble_linear_terms();
  }

  for(auto &s : this->systems)
  {
    //hold the velocity space, we'll need it...
    const FESpace * v_space = s.velocity_space.get();
    const FESpace * p_space = s.pressure_space.get();
    
#ifdef _MPI
    //MPI: solution in consistency level 3 (TODO: this might be superfluous here)
    for (size_t bl = 0; bl < s.solution.n_blocks(); ++bl)
    {
      s.matrix.get_communicators()[bl]->consistency_update(s.solution.block(bl),
                                                           3);
    }
#endif

    // spaces for matrices
    const FESpace *spaces[2] = {v_space, p_space};
    for(int i = 0; i < d; ++i)
      rhs_spaces[i] = v_space;
    rhs_spaces[d] = p_space;
    std::vector<std::vector<size_t>> cells;
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
    }
    feFunction[d] = &s.p;

    //decide wether to assemble by upwinding or not
    bool finest_grid = (&s == &systems.at(0));
    bool do_upwinding = (db["space_discretization_type"].is("upwind")
                        || (mdml && !finest_grid))
                        && !is_stokes;
    // local assembling object    
    LocalAssembling<d> la(this->db, LocalAssembling_type::NavierStokesNL, 
                          feFunction, example.get_coeffs());
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
      Assemble2D(
#endif
                 n_fe_spaces, spaces, n_square_matrices, sq_matrices.data(),
                 n_rectangular_matrices, re_matrices.data(), n_rhs,
                 rhs_array.data(), rhs_spaces.data(), boundCondition.data(),
                 boundValues.data(), la);
#ifdef __2D__
//     AddGadientJump(n_fe_spaces, spaces, 
//                     n_square_matrices, sq_matrices.data(), 
//                     0, nullptr, n_rhs,
//                     nullptr, nullptr, 
//                     boundCondition.data(), 
//                     boundValues.data(), la);
    
    Assemble_DG<2>(la.GetCoeffFct(), 1, spaces, 1, &sq_matrices[0], 0, nullptr, 0, nullptr, boundCondition.data(), boundValues.data(), db, "Oseen", 1);
#endif      
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
    
    // assemble on the boundary if needed (only the one needed for non-linear iteration: A-block)
    assemble_nonlinear_boundary_terms();

    //TODO: Copying non-actives??
    if(db["space_discretization_type"].is("supg"))
      s.solution.copy_nonactive(s.rhs);
    
    if(v_space->get_n_hanging() != 0 || p_space->get_n_hanging() != 0)
    {
      if(TDatabase::ParamDB->NSTYPE)
      {
        // remove entries in hanging and Dirichlet rows in non-diagonal blocks
        for(int i = 0; i < d+1; ++i)
        {
          for(int j = 0; j < d+1; ++j)
          {
            if(i == j)
              continue; // skip diagonal blocks
            blocks[(d+1)*i+j]->resetNonActive();
          }
        }
      }
      else
      {
        ErrThrow("hanging nodes and NSTYPE ", TDatabase::ParamDB->NSTYPE,
                 " is not supported.");
      }
    }
  }// endfor auto grid
}


/* ************************************************************************* */
template <int d>
bool NavierStokes<d>::stop_it(unsigned int iteration_counter)
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

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::compute_residuals()
{
  System_per_grid& s = this->systems.front();
  unsigned int n_u_dof = s.solution.length(0);
  unsigned int n_p_dof = s.solution.length(d);

  // copy rhs to defect and compute defect
#ifdef _MPI
    //MPI: solution in consistency level 3 (TODO: maybe this is superfluous here
    // (because solution might be in level 3 consistency already)!)
    auto comms = s.matrix.get_communicators();
    for (size_t bl = 0; bl < comms.size() ;++bl)
    {
      comms[bl]->consistency_update(s.solution.block(bl), 3);
    }
#endif

  defect = s.rhs;
  s.matrix.apply_scaled_add(s.solution, defect, -1.);

  if(s.matrix.pressure_projection_enabled())
  {
    FEFunction defect_fctn(s.pressure_space, "p_defect", &defect[d*n_u_dof]);
    defect_fctn.project_into_L20();
  }
  
  std::vector<unsigned int> velocity_blocks(d, 0);
  std::iota(std::begin(velocity_blocks), std::end(velocity_blocks), 0);
#ifdef _MPI
  std::vector<const TParFECommunicator3D*> velocity_comms(d, nullptr);
  for (int i = 0; i < d; ++i)
  {
    velocity_comms[i] = comms[i];
  }

  double momentum_residual_square = defect.norm(velocity_blocks, velocity_comms);
  double mass_residual_square = defect.norm({d}, {comms[d]});

  double momentum_residual_max = defect.norm_infty(velocity_blocks, velocity_comms);
  double mass_residual_max = defect.norm_infty({d}, {comms[d]});

#else
  double momentum_residual_square = defect.norm(velocity_blocks);
  double mass_residual_square = defect.norm({d});

  double momentum_residual_max = defect.norm_infty(velocity_blocks);
  double mass_residual_max = defect.norm_infty({d});
#endif

  // the struct 'Residuals' takes squares:
  momentum_residual_square *= momentum_residual_square;
  mass_residual_square *= mass_residual_square;

  Residuals current_residuals(momentum_residual_square, mass_residual_square,
    momentum_residual_max, mass_residual_max);

  old_residuals.add(current_residuals);
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::solve()
{
  System_per_grid& s = this->systems.front();
  double damping = this->db["nonlinloop_damping_factor"];
  // store previous solution for damping, it is a pointer so that we can avoid
  // the copy in case of no damping
  std::shared_ptr<BlockVector> old_solution(nullptr);
  if(damping != 1.0)
    old_solution = std::make_shared<BlockVector>(s.solution);  
  // solving:
#ifdef _MPI
  if(this->solver.get_db()["solver_type"].is("direct"))
  {
    if(damping != 1.0)
      Output::warn("NSE3D::solve", "damping in an MPI context is not tested");

    //set up a MUMPS wrapper
    MumpsWrapper mumps_wrapper(s.matrix);

    //kick off the solving process
    mumps_wrapper.solve(s.rhs, s.solution);
  }
  else
#endif
  this->solver.solve(s.matrix, s.rhs, s.solution);
  
  if(damping != 1.0)
  {
    s.solution.scale(damping);
    s.solution.add_scaled(*old_solution, 1-damping);
  }

  // project pressure if necessary
  if(s.matrix.pressure_projection_enabled())
    s.p.project_into_L20();

  

}

/* ************************************************************************* */
// this is a bad implementation, we need to find ways to pass parameters to the
// TFEFunction2D::GetErrors method
double delta0 = 0.;
template <int d>
void natural_error_norm_infsup_stabilizations(
  int N_Points, std::array<const double*, d> xyz, const double *AbsDetjk,
  const double *Weights, double hK, const double *const* Der,
  const double *const* Exact, const double *const* coeffs, double *LocError);
template <int d>
void parameter_function_for_errors(const double *in, double *out)
{
  // d=2: u1, u2, u1x u1y, u2x, u2y
  // d=3: u1, u2, u3, u1x u1y, u1z, u2x, u2y, u2z, u3x, u3y, u3z
  for(int i = 0; i < d*(d+1); ++i)
    out[i] = in[2+i];
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::output(int i)
{
  int my_rank =0;
#ifdef _MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

  System_per_grid& s = this->systems.front();
  std::array<std::unique_ptr<FEFunction>, d> u_components;
  std::array<FEFunction*, d> velocity_components;
  for(int i = 0; i < d; ++i)
  {
    u_components[i] = s.u.GetComponent(i);
    velocity_components[i] = u_components[i].get();
  }

  if(db["output_compute_minmax"])
  {
    for(int i = 0; i < d; ++i)
      velocity_components[i]->PrintMinMax("u" + std::to_string(i));
    s.p.PrintMinMax(std::string("p"));
  }
  
  if(i < 0)
    outputWriter.write();
  else
    outputWriter.write(i);
  checkpoint_io.write(this->get_solution());
  
#ifdef __2D__
 // velocity_over_line({0.5, 0.}, {0.5, 1.}, 81, velocity_components);

  // geothermal tests flow rate
 /*velocity_over_line({0.15, 0.4}, {0.15, 0.6}, 30, velocity_components, "velocity_channel_sigmae5_mueffe-1_x0.15.txt");
 velocity_over_line({0.3, 0.4}, {0.3, 0.6}, 30, velocity_components, "velocity_channel_sigmae5_mueffe-1_x0.3.txt");
 velocity_over_line({0.45, 0.4}, {0.45, 0.6}, 3, velocity_components, "velocity_channel_sigmae5_mueffe-1_x0.45.txt");
 velocity_over_line({0.6, 0.4}, {0.6, 0.6}, 30, velocity_components, "velocity_channel_sigmae5_mueffe-1_x0.6.txt");
 velocity_over_line({0.75, 0.4}, {0.75, 0.4}, 30, velocity_components, "velocity_channel_sigmae5_mueffe-1_x0.75.txt");
*/
#endif

  // measure errors to known solution
  // If an exact solution is not known, it is usually set to be zero, so that
  // in such a case here only integrals of the solution are computed.
  if(db["output_compute_errors"])
  {
    std::vector<std::array<double, 4>> computed_errors(d+1);
#ifdef __3D__
    TAuxParam3D aux;
    MultiIndex3D nsAllDerivs[d+1] = {MultiIndex3D::D000, MultiIndex3D::D100,
                                     MultiIndex3D::D010, MultiIndex3D::D001};
#else
    TAuxParam2D aux, aux2;
    MultiIndex2D nsAllDerivs[d+1] = {MultiIndex2D::D00, MultiIndex2D::D10,
                                     MultiIndex2D::D01};
#endif
    const FESpace *velocity_space = this->systems.front().velocity_space.get();
    const FESpace *pressure_space = this->systems.front().pressure_space.get();
    
    // errors in the velocity components
    for(int i = 0; i < d; ++i)
    {
      auto ui = velocity_components[i];
      ui->GetErrors(example.get_exact(i), d+1, nsAllDerivs, d, 
                    L2H1Errors,
                    example.get_coeffs(), &aux, 1, &velocity_space,
                    computed_errors[i].data());
    }
    // error in divergence
    double div_error = s.u.GetL2NormDivergenceError(example.get_exact(0),
                                                    example.get_exact(1)
#ifdef __3D__
                                                    , example.get_exact(2)
#endif
                                                   );
    // errors in pressure
    s.p.GetErrors(example.get_exact(d), d+1, nsAllDerivs, d, L2H1Errors,
                  nullptr, &aux, 1, &pressure_space, computed_errors[d].data());


#ifdef __2D__
    int boundary_component_id;
    double un_boundary_error = 0.;
    double boundary_error_l2_squared = 0.;
    double boundary_error_l2 = 0;
    const ParameterDatabase e_db = example.get_database();
    int n_nitsche_bd = e_db["n_nitsche_bd"];
  if (n_nitsche_bd)
  {
    std::vector<size_t> nitsche_id = e_db["nitsche_id"];
      std::vector<double> nitsche_penalty = e_db["nitsche_penalty"];
      for(unsigned int k = 0; k < nitsche_id.size(); k++)
      {
        boundary_component_id = nitsche_id[k];

        double boundary_error_l2_u0[1], boundary_error_l2_u1[1];

        velocity_components[0]->GetL2BoundaryError(example.get_bd(0),
                &aux2, 1, &velocity_space,
                boundary_error_l2_u0, boundary_component_id);
        velocity_components[1]->GetL2BoundaryError(example.get_bd(1),
                &aux2, 1, &velocity_space,
                boundary_error_l2_u1, boundary_component_id);

        boundary_error_l2_squared += boundary_error_l2_u0[0] + boundary_error_l2_u1[0];

        // compute the L2-norm of the normal velocity error at the Nitsche boundaries
        un_boundary_error += s.u.GetL2NormNormalComponentError(example.get_bd(0),
                example.get_bd(1), boundary_component_id);
      }
      un_boundary_error = std::sqrt(un_boundary_error);
      boundary_error_l2 = std::sqrt(boundary_error_l2_squared);
  }
#endif


#ifdef _MPI
    constexpr int n_send = 2*(d+1)+1;
    double err_red[n_send]; //memory for global (across all processes) error
    double err_send[n_send]; //fill send buffer
    err_send[0] = computed_errors[0][0];
    err_send[1] = computed_errors[0][1];
    err_send[2] = computed_errors[1][0];
    err_send[3] = computed_errors[1][1];
    err_send[4] = computed_errors[2][0];
    err_send[5] = computed_errors[2][1];
    err_send[6] = div_error*div_error;
    err_send[7] = computed_errors[d][0];
    err_send[8] = computed_errors[d][1];

    MPI_Allreduce(err_send, err_red, n_send, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for(i = 0; i < n_send; i++)
    { //MPI: std::sqrt was skipped in GetErrors function - do it here globally!
      err_red[i] = std::sqrt(err_red[i]);
    }
    //fill the reduced errors back where they belong
    computed_errors[0][0] = err_red[0];
    computed_errors[0][1] = err_red[1];
    computed_errors[1][0] = err_red[2];
    computed_errors[1][1] = err_red[3];
    computed_errors[2][0] = err_red[4];
    computed_errors[2][1] = err_red[5];
    div_error = err_red[6];
    computed_errors[d][0] = err_red[7];
    computed_errors[d][1] = err_red[8];
#endif

    double l2_error = 0;
    double h1_error = 0;
    for(int i = 0; i < d; ++i)
    {
      l2_error += computed_errors[i][0] * computed_errors[i][0];
      h1_error += computed_errors[i][1] * computed_errors[i][1];
    }

    errors.at(0) = std::sqrt(l2_error);
    errors.at(1) = div_error;
    errors.at(2) = std::sqrt(h1_error);
    errors.at(3) = computed_errors[d][0];
    errors.at(4) = computed_errors[d][1];
#ifdef __2D__
    errors.at(6) = boundary_error_l2;
    errors.at(7) = un_boundary_error;
#endif

    //print errors
    if(my_rank == 0)
    {
      Output::print("--------------------------------------------");
      Output::stat("NavierStokes", "Measured errors");
      Output::dash("L2(u)     : ", setprecision(14), errors.at(0));
      Output::dash("L2(div(u)): ", setprecision(14), errors.at(1));
      Output::dash("H1-semi(u): ", setprecision(14), errors.at(2));
      Output::dash("L2(p)     : ", setprecision(14), errors.at(3));
      Output::dash("H1-semi(p): ", setprecision(14), errors.at(4));
#ifdef __2D__
      Output::dash("L2(u)_boundary: ", setprecision(14), errors.at(6));
      Output::dash("L2(u.n)_boundary: ", setprecision(14), errors.at(7));
#endif
      Output::print("--------------------------------------------");
    }

    auto sdt(db["space_discretization_type"]);
    bool pspg = sdt.is("pspg");
    bool symm_gls = sdt.is("symm_gls");
    bool nonsymm_gls = sdt.is("nonsymm_gls");
    if(pspg || symm_gls || nonsymm_gls)
    {
      double nu = example.get_nu();
      double error_in_natural_norm = nu * errors[2]*errors[2];
      if(symm_gls)
        error_in_natural_norm += (1./nu) * errors[3]*errors[3];
      delta0 = db["pspg_delta0"];
      int beginparameter = 0;
      ParamFct * parameter_function = &parameter_function_for_errors<d>;
#ifdef __3D__
      int fevalue_fctindex[12] = {0, 1, 2, 0, 0, 0, 1, 1, 1, 2, 2, 2};
      MultiIndex3D fevalue_multiindex[] = {
            MultiIndex3D::D000, MultiIndex3D::D000, MultiIndex3D::D000,
            MultiIndex3D::D200, MultiIndex3D::D020, MultiIndex3D::D002,
            MultiIndex3D::D200, MultiIndex3D::D020, MultiIndex3D::D002,
            MultiIndex3D::D200, MultiIndex3D::D020, MultiIndex3D::D002};
      TAuxParam3D NSE_aux(1, 3, 1, d*(d+1), nullptr, velocity_components.data(),
                          &parameter_function, fevalue_fctindex,
                          fevalue_multiindex, d*(d+1), &beginparameter);
#else
      int fevalue_fctindex[6] = {0, 1, 0, 0, 1, 1};
      MultiIndex2D fevalue_multiindex[6] = 
         { MultiIndex2D::D00, MultiIndex2D::D00, MultiIndex2D::D20,
           MultiIndex2D::D02, MultiIndex2D::D20, MultiIndex2D::D02 };
      TAuxParam2D NSE_aux(1, d*(d+1), velocity_components.data(),
                          &parameter_function, fevalue_fctindex,
                          fevalue_multiindex, d*(d+1), &beginparameter);
#endif
      double err[4];
      s.p.GetErrors(example.get_exact(2), d+1, nsAllDerivs, 1,
                    natural_error_norm_infsup_stabilizations<d>,
                    get_example().get_coeffs(), &NSE_aux, 1,
                    &pressure_space, err);
      error_in_natural_norm += err[0]*err[0];
      error_in_natural_norm = std::sqrt(error_in_natural_norm);
      errors[5] = error_in_natural_norm;
      Output::print("Error in natural(", sdt, ") norm: ", std::setprecision(14),
                    errors[5]);
    }
  } // if(this->db["compute_errors"])

  bool write_matrix = false;
  if (write_matrix)
  {
    // create an output file containing the whole FE matrix. This can be read into Matlab using the Matlab function mmread.m
    std::stringstream matrix_name;
    matrix_name << "Coeff_Matrix_matrixmarket";
    s.matrix.get_combined_matrix()->write(matrix_name.str());
  }

  //do postprocessing step depending on what the example implements
  example.do_post_processing(*this);
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::adjust_pressure()
{
  System_per_grid& s = this->systems.front();
  if(db["problem_type"].is(5)) // Navier--Stokes
  {
    int sign = 0;
    if(db["nse_nonlinear_form"].is("rotational"))
      sign = -1;
    if(db["nse_nonlinear_form"].is("emac"))
      sign = 1;
    if(sign)
    {
      std::array<std::unique_ptr<FEFunction>, d> u_components;
      std::array<FEFunction*, d> velocity_components;
      for(int i = 0; i < d; ++i)
      {
        u_components[i] = s.u.GetComponent(i);
        velocity_components[i] = u_components[i].get();
      }
      typename FEFunction::AnalyticFunction 
      f = [&velocity_components, sign](const TBaseCell* cell, int i, 
                                       std::array<double, d> xyz)
          {
            double val = 0;
            for(int c = 0; c < d; ++c)
            {
              double val_ui;
#ifdef __3D__
              velocity_components[c]->FindValueLocal(cell, i, xyz[0], xyz[1],
                                                     xyz[2], &val_ui);
#else
              velocity_components[c]->FindValueLocal(cell, i, xyz[0], xyz[1],
                                                     &val_ui);
#endif
              val += val_ui*val_ui;
            }
            return sign * 0.5 * val;
          };
      s.p.add(f);
      // project pressure if necessary
      if(s.matrix.pressure_projection_enabled())
        s.p.project_into_L20();
    }
  }
}

/* ************************************************************************* */
template <int d>
std::unique_ptr<typename NavierStokes<d>::FEFunction>
  NavierStokes<d>::get_velocity_component(int i)
{
  if(i >= 0 && i < d)
    return this->systems.front().u.GetComponent(i);
  else
    ErrThrow("There are only ", d, " velocity components!");
}

/* ************************************************************************* */
template <int d>
const Residuals& NavierStokes<d>::get_residuals() const
{
  return old_residuals.back();
}

/* ************************************************************************* */
template <int d>
double NavierStokes<d>::get_impuls_residual() const
{
  return old_residuals.back().momentumResidual;
}

/* ************************************************************************* */
template <int d>
double NavierStokes<d>::get_mass_residual() const
{
  return old_residuals.back().massResidual;
}

/* ************************************************************************* */
template <int d>
double NavierStokes<d>::get_full_residual() const
{
  return old_residuals.back().fullResidual;
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::reset_residuals()
{
  this->old_residuals = FixedSizeQueue<10, Residuals>();
}

/* ************************************************************************* */
template<int d> std::array<double, 8> NavierStokes<d>::get_errors() const
{
  return errors;
}

/* ************************************************************************* */
template <int d>
void natural_error_norm_infsup_stabilizations(int N_Points,
                                              std::array<const double*, d>,
                                              const double *AbsDetjk,
                                              const double *Weights, double hK,
                                              const double *const* Der,
                                              const double *const* Exact,
                                              const double *const* coeffs,
                                              double *LocError)
{
  LocError[0] = 0.0;
  for(int i=0;i<N_Points;i++)
  {
    double nu = coeffs[i][0];
    double delta = compute_PSPG_delta(delta0, hK, nu);
    const double *deriv = Der[i];
    const double *exactval = Exact[i];
    double w = delta*Weights[i]*AbsDetjk[i];

    double t = deriv[1]-exactval[1];
    LocError[0] += w*t*t;
      
    t = deriv[2]-exactval[2];
    LocError[0] += w*t*t;
    if(d == 3)
    {
      t = deriv[3]-exactval[3];
      LocError[0] += w*t*t;
    }
  }
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::assemble_nonlinear_boundary_terms()
{
  const ParameterDatabase e_db = example.get_database();
  int n_nitsche_bd = e_db["n_nitsche_bd"];
  int n_windkessel_bd = e_db.try_get_value("n_windkessel_bd", 0);

  int my_rank = 0;
#ifdef _MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

  for(System_per_grid& s : this->systems)
  {
    const FESpace* v_space = s.velocity_space.get();

    // Windkessel BC
    if(n_windkessel_bd)
    {
      std::vector<size_t> windkessel_id = e_db["windkessel_id"];
      std::vector<double> windkessel_Rp = e_db["windkessel_Rp"];
      double windkessel_p0              = e_db["windkessel_distal_pressure"];

      for(size_t k = 0; k < windkessel_id.size(); k++)
      {
        double q_out = s.u.compute_flux(windkessel_id[k]);
        double q_exp = flux_face(k);
        double q_rel = (q_out-q_exp)/q_exp;
        double tolerance = 1e-3;
        // stabilize
        if(std::abs(q_rel) < tolerance)
        {
          if(my_rank == 0)
          {
            Output::print(" Windkessel BC on boundary: ", windkessel_id[k],
                          ", resistance = ",windkessel_Rp[k],
                          ", flux = ", q_out, " (Q_rel < ", tolerance, ")",
                          ", set p_new = p_old = ", p_out_old[2*k+1]);
          }
          p_out_old[2*k] = p_out_old[2*k+1];
          continue;
        }
        // solve windkessel model
        double p_out = windkessel_Rp[k]*q_out + windkessel_p0;
        double p_update = p_out - p_out_old[2*k+1];
        double p_update_m1 = p_out_old[2*k+1] - p_out_old[2*k];
        double sgn = (p_update == 0.) ? 0. : ((p_update > 0.) ? 1. : -1.);
        double sgn_m1 = (p_update_m1 == 0.) ? 0.
                                            : ((p_update_m1 > 0.) ? 1. : -1.);
        if(sgn * sgn_m1 < 0.)
        {
          windkessel_damping_limiter[k] = 0.5 * windkessel_damping_limiter[k]
                                         / (1. + windkessel_damping_limiter[k]);
        }
        // limit p_out
        if(windkessel_damping_limiter[k] != 0.)
        {
          p_update = sgn * std::min(std::abs(p_update),
                      std::abs(windkessel_damping_limiter[k]*p_out_old[2*k+1]));
          p_out = p_out_old[2*k+1] + p_update;
        }

        if(my_rank == 0)
        {
          Output::print(" Windkessel BC on boundary: ", windkessel_id[k],
                        ", resistance = ",windkessel_Rp[k],
                        ", flux = ", q_out, " (Q_rel: ", 100*q_rel, "%)",
                        ", p_old = ", p_out_old[2*k+1], " p_new = ", p_out,
                        " (P_step: ", p_update, ")");
        }
#ifdef __2D__
      BoundaryAssembling2D::rhs_g_v_n(s.rhs,
                                      v_space,
                                      nullptr,
                                      windkessel_id[k],
                                      -1.*p_update);
#else
      std::vector<TBoundFace*> boundaryFaceList;
      boundaryFaceList.clear();
      v_space->GetCollection()->get_face_list_on_component(windkessel_id[k],
                                                           boundaryFaceList);
      BoundaryAssembling3D ba;
      ba.rhs_g_v_n(s.rhs,
                   v_space,
                   nullptr,
                   boundaryFaceList,
                   windkessel_id[k],
                   -1.*p_update);
#endif
      p_out_old[2*k] = p_out_old[2*k+1];
      p_out_old[2*k+1] = p_out;
      }
    } /* if(n_windkessel_bd) */

    if (n_nitsche_bd)
    {
      std::vector<size_t> nitsche_id = e_db["nitsche_id"];
      std::vector<double> nitsche_penalty = e_db["nitsche_penalty"];
      double effective_viscosity = this->example.get_nu();

#ifdef __3D__
      auto coll = s.velocity_space.get()->GetCollection();
      std::vector<TBoundFace*> boundaryFaceList;
      boundaryFaceList.clear();
#endif
      
      for (unsigned int k = 0; k < nitsche_id.size(); k++)
      {
	int sym_u = e_db["symmetric_nitsche_u"];
	Output::print<2>(" (nonlinear) Nitsche BC on boundary: ", nitsche_id[k], ", nitsche penalty: ", nitsche_penalty[k]);


#ifdef __2D__
  double sigma = this->example.get_inverse_permeability();
  double L_0 = db["L_0"];
  
	if (e_db["example"].is(15))
	  Output::print(" ** ERROR: Example 15 should onlybe used in a Stokes problem ");
        else
        {
          

	  BoundaryAssembling2D::nitsche_bc_nonlinear_iteration(s.matrix, 
							       v_space, 
							       nitsche_id[k], nitsche_penalty[k],
							       effective_viscosity, sigma, L_0,
							       sym_u);
        }


	///@todo: corner stabilization for Brinkman + penalty free
#else
	
  int sym_p = e_db["symmetric_nitsche_p"];

	// Nitsche penalty for weak essential BC
	coll->get_face_list_on_component(nitsche_id[k], boundaryFaceList);
        const TFESpace3D * p_space = s.pressure_space.get();
	BoundaryAssembling3D ba;
	ba.nitsche_bc_nonlinear_iteration(s.matrix,
					  v_space, p_space,
					  boundaryFaceList,
					  nitsche_id[k], nitsche_penalty[k],
					  effective_viscosity,sym_u, sym_p);

#endif
    
      }
    }
  }
}

/* ************************************************************************* */
template <int d>
double NavierStokes<d>::flux_face(int face_i)
{
  const ParameterDatabase e_db = example.get_database();
  std::vector<size_t> windkessel_id = e_db["windkessel_id"];
  return systems.front().u.compute_flux(windkessel_id[face_i]);
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::assemble_boundary_terms()
{
  const ParameterDatabase e_db = example.get_database();
  int n_windkessel_bd = e_db.try_get_value("n_windkessel_bd", 0);
  int n_neumann_bd = e_db["n_neumann_bd"];
  int n_nitsche_bd = e_db["n_nitsche_bd"];

  int my_rank = 0;
#ifdef _MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

  ///@todo this part of the code needs still to be implemented dimension-independent
  for(System_per_grid& s : this->systems)
  {
    const FESpace* v_space = s.velocity_space.get();

    if(n_windkessel_bd)
    {
      // Windkessel BC
      std::vector<size_t> windkessel_id = e_db["windkessel_id"];
      std::vector<double> windkessel_Rp = e_db["windkessel_Rp"];
      double windkessel_p0              = e_db["windkessel_distal_pressure"];
      double w_d_l = e_db["windkessel_damping_limiter"];

      for(unsigned int k = 0; k <  windkessel_id.size(); k++)
      {
//         double q_out = s.u.compute_flux(windkessel_id[k]);
        double q_out = flux_face(k);
        // solve windkessel model
        double p_out = windkessel_Rp[k]*q_out + windkessel_p0;
        p_out_old.push_back(p_out);
        p_out_old.push_back(p_out);
        windkessel_damping_limiter.push_back(w_d_l);
        if(my_rank == 0)
        {
          Output::print<2>(" Windkessel BC on boundary: ", windkessel_id[k],
                           ", resistance = ", windkessel_Rp[k], "\n",
                           " == Q( ",windkessel_id[k],") = ",q_out,
                           " == P = ", p_out
                          );
        }

#ifdef __2D__
        BoundaryAssembling2D::rhs_g_v_n(s.rhs,
                                        v_space,
                                        nullptr,
                                        windkessel_id[k],
                                        -1.*p_out);
#else
        std::vector<TBoundFace*> boundaryFaceList;
        boundaryFaceList.clear();
        v_space->GetCollection()->get_face_list_on_component(windkessel_id[k],
                                                             boundaryFaceList);
//         BoundaryAssembling3D::rhs_g_v_n(s.rhs,
        BoundaryAssembling3D ba;
        ba.rhs_g_v_n(s.rhs,
                                        v_space,
                                        nullptr,
                                        boundaryFaceList,
                                        windkessel_id[k],
                                        -1.*p_out);
#endif
      }
    } /* if(n_windkessel_bd) */

    if(n_neumann_bd)
    {
      // Neumann BC
      std::vector<size_t> neumann_id = e_db["neumann_id"];
      std::vector<double> neumann_value = e_db["neumann_value"];

      for(unsigned int k = 0; k < neumann_id.size(); k++)
      {
        if(my_rank == 0)
        {
          Output::print<2>(" Neumann BC on boundary: ", neumann_id[k],
                           ", value = ",neumann_value[k] );
        }

#ifdef __2D__
        BoundaryAssembling2D::rhs_g_v_n(s.rhs,
                                        v_space,
                                        nullptr,
                                        neumann_id[k],
                                        -1.*neumann_value[k]);
#else
        std::vector<TBoundFace*> boundaryFaceList;
        boundaryFaceList.clear();
        v_space->GetCollection()->get_face_list_on_component(neumann_id[k],
                                                             boundaryFaceList);
//         BoundaryAssembling3D::rhs_g_v_n(s.rhs,
        BoundaryAssembling3D ba;
        ba.rhs_g_v_n(s.rhs,
                                        v_space,
                                        nullptr,
                                        boundaryFaceList,
                                        neumann_id[k],
                                        -1.*neumann_value[k]);
#endif
      }
    } /* if(n_neumann_bd) */

    if(n_nitsche_bd)
    {
      // Nitsche penalty for weak essential BC
      std::vector<size_t> nitsche_id = e_db["nitsche_id"];
      std::vector<double> nitsche_penalty = e_db["nitsche_penalty"];
      double effective_viscosity = this->example.get_nu();

#ifdef __2D__
      for(unsigned int k = 0; k < nitsche_id.size(); k++)
      {
        const FESpace * p_space = s.pressure_space.get();

        int sym_u = e_db["symmetric_nitsche_u"];
        int sym_p = e_db["symmetric_nitsche_p"];
        double sigma = this->example.get_inverse_permeability();
        double L_0 = db["L_0"];

        if(e_db["example"].is(15))
        {
          if(nitsche_id.size() == 2 || nitsche_id[0] == 1 || nitsche_id[1] == 4)
          {
            double effective_viscosity_porous = 0;
            double effective_viscosity_stokes = 0.001;
            double sigma_porous =  1.e7;
            double sigma_stokes =  0.;

            if(my_rank == 0)
            {
              Output::print<2>(" Nitsche BC on boundary: ", nitsche_id[0],
                               ", nitsche penalty: ", nitsche_penalty[0]);
            }
            BoundaryAssembling2D::nitsche_bc(s.matrix, s.rhs,
                                             v_space, p_space,
                                             this->example.get_bd(0),
                                             this->example.get_bd(1),
                                             nitsche_id[0], nitsche_penalty[0],
                                             effective_viscosity_porous,
                                             sigma_porous, L_0,
                                             sym_u, sym_p);

            if(my_rank == 0)
            {
              Output::print<2>(" Nitsche BC on boundary: ", nitsche_id[1],
                               ", nitsche penalty: ", nitsche_penalty[1]);
            }
            BoundaryAssembling2D::nitsche_bc(s.matrix, s.rhs,
                                             v_space, p_space,
                                             this->example.get_bd(0),
                                             this->example.get_bd(1),
                                             nitsche_id[1], nitsche_penalty[1],
                                             effective_viscosity_stokes,
                                             sigma_stokes, L_0,
                                             sym_u, sym_p);
            break;
          }
          else
          {
            ErrThrow("The riverbed example (15) is only hard coded for "
                     "the Nitsche method on top and bottom boundaries."
                     "The reason is that variations in coefficients defined "
                     "in the example file are not used in the Nitsche "
                     "functions.!!!");
            break;
          }
        } /* if(e_db["example"].is(15)) */
        else
        {
          if(my_rank == 0)
          {
            Output::print<2>(" Nitsche BC on boundary: ", nitsche_id[k], 
                             ", nitsche penalty: ", nitsche_penalty[k]);
          }
          BoundaryAssembling2D::nitsche_bc(s.matrix, s.rhs,
                                           v_space, p_space,
                                           this->example.get_bd(0),
                                           this->example.get_bd(1),
                                           nitsche_id[k], nitsche_penalty[k],
                                           effective_viscosity,
                                           sigma, L_0,
                                           sym_u, sym_p);
        }
      } /* for(unsigned int k = 0; k < nitsche_id.size(); k++) */

      double corner_stab = e_db["corner_stab"];
      if(corner_stab)
      {
        if(my_rank == 0)
        {
          Output::print<2>(" Corner stabilization is applied, corner_stab = ",
                           corner_stab);
        }
        double sigma = this->example.get_inverse_permeability();
        double L_0 = db["L_0"];
        corner_stab = corner_stab * (effective_viscosity + sigma * L_0* L_0);

        if(e_db["example"].is(15))
        {
          ErrThrow("The riverbed example (15) does not support the corner "
                   "stabilization currently.");
          break;
        }
        else
        {
          BoundaryAssembling2D::matrix_and_rhs_corner_stabilization(
                                                        s.matrix, s.rhs,
                                                        v_space,
                                                        this->example.get_bd(0),
                                                        this->example.get_bd(1),
                                                        nitsche_id,
                                                        corner_stab);
        }
      } /* if(corner_stab) */
#else
      const FESpace * p_space = s.pressure_space.get();
      std::vector<TBoundFace*> boundaryFaceList;
      boundaryFaceList.clear();
      int sym_u = e_db["symmetric_nitsche_u"];
      int sym_p = e_db["symmetric_nitsche_p"];
      //double sigma = this->example.get_inverse_permeability();
      //double L_0 = db["L_0"];

      for (size_t k = 0; k < nitsche_id.size(); k++)
      {
        v_space->GetCollection()->get_face_list_on_component(nitsche_id[k],
                                                             boundaryFaceList);
        if(my_rank == 0)
        {
          Output::print<2>(" Nitsche BC on boundary: ", nitsche_id[k],
                           ", nitsche penalty: ", nitsche_penalty[k]);
          Output::print<5>("boundaryFaceList.size(): ",
                            boundaryFaceList.size() );
        }
//         BoundaryAssembling3D::nitsche_bc(s.matrix, s.rhs,
        BoundaryAssembling3D ba;
        ba.nitsche_bc(s.matrix, s.rhs,
                                         v_space, p_space,
                                         this->example.get_bd(0),
                                         this->example.get_bd(1),
                                         this->example.get_bd(2),
                                         nullptr,
                                         boundaryFaceList,
                                         nitsche_id[k], nitsche_penalty[k],
                                         effective_viscosity,
                                         // sigma, L_0,
                                         sym_u, sym_p);
      }
#endif
    } /* if(n_nitsche_bd) */
  } /* for(System_per_grid& s : this->systems) */


//     auto coll = s.velocity_space.get()->GetCollection();
//     BoundaryAssembling3D ba;
// 
// 
//     if (n_windkessel_bd)
//     {
//       // Windkessel BC
//       std::vector<TBoundFace*> boundaryFaceList;
// 
//       std::vector<size_t> windkessel_id = e_db["windkessel_id"];
//       std::vector<double> windkessel_Rp = e_db["windkessel_Rp"];
// 
//       for (size_t k = 0; k < windkessel_id.size(); k++)
//       {
//         double q_out = s.u.compute_flux(windkessel_id[k]);
//         //double q_out = 0.6;
// 
//         if(my_rank == 0)
//         {
//           Output::print(" == Q( ",windkessel_id[k],") = ",q_out);
//         }
//         // solve windkessel model
//         double p_out = windkessel_Rp[k]*q_out;
// 
//         if(my_rank == 0)
//         {
//           Output::print<2>(" Windkessel BC on boundary: ", windkessel_id[k]);
//         }
//         coll->get_face_list_on_component(windkessel_id[k], boundaryFaceList);
//         const TFESpace3D * v_space = s.velocity_space.get();
//         ba.rhs_g_v_n(s.rhs,
//                      v_space,
//                      nullptr,
//                      boundaryFaceList,
//                      (int) windkessel_id[k],
//                      -1.*p_out);
//       }
//     }
// 
//     if (n_neumann_bd)
//     {
//       // Neumann BC
//       std::vector<TBoundFace*> boundaryFaceList;
//       boundaryFaceList.clear();
//       std::vector<size_t> neumann_id = e_db["neumann_id"];
//       std::vector<double> neumann_value = e_db["neumann_value"];
// 
//       std::vector<TBaseCell*> dummy;
//       for (size_t k = 0; k < neumann_id.size(); k++)
//       {
//         if(my_rank == 0)
//         {
//           Output::print<2>(" Neumann BC on boundary: ", neumann_id[k],
//                            " value = ", neumann_value[k]);
//         }
//         coll->get_face_list_on_component(neumann_id[k], boundaryFaceList);
//         const TFESpace3D * v_space = s.velocity_space.get();
//         ba.rhs_g_v_n(s.rhs,
//                      v_space,
//                      nullptr,
//                      boundaryFaceList,
//                      (int) neumann_id[k],
//                      -1.*neumann_value[k]);
//       }
//     }
// 
//     if (n_nitsche_bd)
//     {
//       // Nitsche penalty for weak essential BC
//       std::vector<TBoundFace*> boundaryFaceList;
//       boundaryFaceList.clear();
//       std::vector<size_t> nitsche_id = e_db["nitsche_id"];
//       std::vector<double> nitsche_penalty = e_db["nitsche_penalty"];
// 
//       for (size_t k = 0; k < nitsche_id.size(); k++)
//       {
//         if(my_rank == 0)
//         {
//           Output::print<2>(" Nitsche BC on boundary: ", nitsche_id[k],
//                            ", nitsche penalty: ", nitsche_penalty[k]);
//         }
//         coll->get_face_list_on_component(nitsche_id[k], boundaryFaceList);
// 
//         if(my_rank == 0)
//         {
//           Output::print<5>("boundaryFaceList.size(): ",
//                            boundaryFaceList.size() );
//         }
//         const TFESpace3D * v_space = s.velocity_space.get();
//         const TFESpace3D * p_space = s.pressure_space.get();
// 
//         double effective_viscosity = this->example.get_nu();
//         int sym_u = e_db["symmetric_nitsche_u"];
//         int sym_p = e_db["symmetric_nitsche_p"];
//         //double sigma = this->example.get_inverse_permeability();
//         //double L_0 = db["L_0"];
// 
//         ba.nitsche_bc(s.matrix, s.rhs, v_space, p_space,
//                       this->example.get_bd(0),
//                       this->example.get_bd(1),
//                       this->example.get_bd(2),
//                       nullptr,
//                       boundaryFaceList,
//                       nitsche_id[k], nitsche_penalty[k],
//                       effective_viscosity,// sigma, L_0,
//                       sym_u, sym_p);
//       }
//     }
//   }
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::velocity_over_line(
  const std::vector<double>& start_point, const std::vector<double>& end_point,
  size_t number_of_points, std::array<FEFunction*, d> velocity_components, std::string name_of_file)
{
  //----------------------------------------------------------------------------------
  // The following output is made for the geometry channel.mesh or channel_simple.mesh
  Output::print("The values the solution u1, u2 (, u3) takes (and the velocity magnitude) at the line between [x1, y1] and [x2, y2] are saved in", name_of_file);
  std::ostringstream oss;
  oss << name_of_file; //"u_values_over_line.txt";

  std::string var = oss.str();
  std::ofstream velfile(var);

  double values_u1[3] = {};
  double values_u2[3] = {};

  for (unsigned int k = 0; k < number_of_points; k++)
  {
    double Y = start_point[1] + k * (end_point[1]-start_point[1])/(number_of_points-1);
    double X = start_point[0] + k * (end_point[0]-start_point[0])/(number_of_points-1); //x[0];

    // for (int i = 0; i < d; i++)
    //{
#ifdef __2D__
    velocity_components[0]->FindGradient(X, Y, values_u1);
    velocity_components[1]->FindGradient(X, Y, values_u2);
#else
      (void) velocity_components;
#endif
    //}
    velfile << "(X, Y) = (" << X << ", " << Y << "); " << " " << " " << " " <<  " (u0, u1, |u|) = ( " << values_u1[0] << ", " << values_u2[0] << ", " << std::sqrt( values_u1[0] * values_u1[0] + values_u2[0] * values_u2[0]) << ")" << endl;
  }
  velfile.close();
}



/////////////// Routines for periodic boundary conditions /////////////////
/* ************************************************************************* */
template <int d>
void NavierStokes<d>::findPeriodicDOFs(double x_left, double x_right)
{
  if (d == 3)
    ErrThrow("findPeriodicDOFs() currently only works in 2D.");

#ifdef __2D__
  // threshold for two doubles to be equal
  const double eps = 1e-8;
  auto * fespace = this->systems.front().velocity_space.get();
  auto coll = this->systems.front().velocity_space.get()->GetCollection();
  int n_u_active = this->systems.front().velocity_space.get()->get_n_active();

  int N_Cells = coll->GetN_Cells();
  // first loop over cells
  for (int cell = 0; cell < N_Cells; cell++)
  {
    TBaseCell *cell1 = coll->GetCell(cell);
    // check if face on boundary
    for (int j1 = 0; j1 < cell1->GetN_Edges(); j1++)
    {
      auto *joint1 = cell1->GetJoint(j1);
      //not on boundary
      if (joint1->GetType() != BoundaryEdge)
        continue;

      const int n_Vert = cell1->GetN_Vertices();
      double x11, x12, y11, y12;
      // compute coordinates of vertices
      cell1->GetVertex(j1)->GetCoords(x11, y11);
      cell1->GetVertex((j1 + 1) % n_Vert)->GetCoords(x12, y12);

      // Todo: change implementation such that periodic dofs in one
      // and the same cell are kept as well

      // check if vertex x = x_left
      bool left = false;
      if (std::abs(x11 - x_left) > eps || std::abs(x12 - x_left) > eps)
      { // one vertex does not lie on the left boundary
        if (std::abs(x11 - x_right) > eps || std::abs(x12 - x_right) > eps)
          continue; // one vertex does not lie on the right boundary
      }
      else
      {
        left = true;
      }

      const FE_type FEid1 = fespace->get_fe_type(cell);
      const FiniteElement & FE1 = fespace->get_fe(cell);

      // global indices of all degrees of freedom in this cell
      const int *globalDOF1 = fespace->GetGlobalDOF(cell);
      // local degrees of freedom which correspond to this edge
      const int* localDOF1 = FE1.GetFEDesc()->GetJointDOF(j1);
      const int N_localDOF1 = FE1.GetFEDesc()->GetN_JointDOF();

      // midpoint of edge (specific for periodicity in x direction)
      const double y1 = (y11 + y12) / 2;

      // find cell which should be coupled to this cell
      // inner loop over the cells
      for (int i2_cell = cell + 1; i2_cell < N_Cells; i2_cell++)
      {
        auto *cell2 = coll->GetCell(i2_cell);
        // check if face on boundary
        for (int j2 = 0; j2 < cell2->GetN_Edges(); j2++)
        {
          const TJoint *joint2 = cell2->GetJoint(j2);
          //not on boundary
          if (joint2->GetType() != BoundaryEdge)
            continue;

          double x21, x22, y21, y22;
          cell2->GetVertex(j2)->GetCoords(x21, y21);
          cell2->GetVertex((j2 + 1) % n_Vert)->GetCoords(x22, y22);

          if (std::abs(x21 - (left ? x_right : x_left)) > eps || std::abs(x22 - (left ? x_right : x_left)) > eps)
            continue; // one vertex does not lie on the correct boundary
          // the two edges are on the correct boundary parts
          // check if their midpoints have the same y-coordinate
          const double y2 = (y21 + y22) / 2;
          if (std::abs(y1 - y2) > eps)
            continue;

          // found two edges will should be identified

          const FE_type FEid2 = fespace->get_fe_type(i2_cell);
          const FiniteElement & FE2 = fespace->get_fe(i2_cell);

          // global indices of all degrees of freedom in this cell
          const int *globalDOF2 = fespace->GetGlobalDOF(i2_cell);
          // local degrees of freedom which correspond to this edge
          const int* localDOF2 = FE2.GetFEDesc()->GetJointDOF(j2);
          const int N_localDOF2 = FE2.GetFEDesc()->GetN_JointDOF();

          if(FEid1 != FEid2)
          {
            ErrThrow("Error in making periodic boundary. ",
                    "Two different finite elements");
          }
          if(N_localDOF1 != N_localDOF2)
          {
            ErrThrow("Error in making periodic boundary. ",
                    "Different numbers of dofs on the periodic boundary");
          }

          /*if(TDatabase::ParamDB->SC_VERBOSE > 2)
            Output::print(" creating a vertical periodic boundary at y=(", y21, ",", y22, ")");
           */

          for (int edge_dof = 0; edge_dof < N_localDOF1; edge_dof++)
          {
            // due to counterclockwise numbering in each cell we have to go
            // through one edge the opposite way:
            const int dof1 = globalDOF1[localDOF1[N_localDOF1 - 1 - edge_dof]];
            const int dof2 = globalDOF2[localDOF2[edge_dof]];

            if (dof1 >= n_u_active && dof2 >= n_u_active)
              continue;
            else if ((dof1 >= n_u_active && dof2 < n_u_active) || (dof1 < n_u_active && dof2 >= n_u_active))
              ErrThrow("Error in findPeriodicDOFs() - active dofs.");

            if (left)
              periodic_dofs[dof1] = dof2;
            else
              periodic_dofs[dof2] = dof1;
          }
        }
      }
    }
  }
  Output::print(" There are ", periodic_dofs.size(),
                " periodic degrees of freedom.");
#else // 3D
  // avoid compiler warnings:
  (void) x_left;
  (void) x_right;
#endif
}


/* ************************************************************************* */
template <int d>
void NavierStokes<d>::makePeriodicBoundary(std::shared_ptr<TMatrix> mat,
    bool stokesMat, bool p )
{
  if (periodic_dofs.empty())
    ErrThrow("called NavierStokes::makePeriodicBoundary with map ",
        "'periodic_dofs' not yet set");

  int do_once = 0;

  if (!mat)
  {
    auto blocks = this->get_matrix().get_blocks_uniquely();
    makePeriodicBoundary(blocks[0], true, true);  // A11
    do_once +=1 ;
    makePeriodicBoundary(blocks[1], false, true); // A12
    makePeriodicBoundary(blocks[2], false, true); // B1T
    makePeriodicBoundary(blocks[3], false, true); // A21
    makePeriodicBoundary(blocks[4], true, true);  // A22
    makePeriodicBoundary(blocks[5], false, true); // B2T
    //makePeriodicBoundary(blocks[6], true, false); // B1
    //makePeriodicBoundary(blocks[7], true, false); // B2
    //makePeriodicBoundary(blocks[8], false, true); // C
    return;
  }

  int const * const rowPtr = mat->get_row_ptr();
  int const * const colPtr = mat->get_vector_columns();
  double const * const entries = mat->GetEntries();
  const int n_rows = mat->get_n_rows();
  // this map will be passed to "mat->changeRows(new_rows)"·
  std::map<int, std::map<int, double> > new_rows;

  for (std::map<int, int>::iterator ii = periodic_dofs.begin();
      ii != periodic_dofs.end() && ii->first < n_rows; ++ii)
  {
    // the row with number "ii->second" is added to the row number "ii->first".
    // Then row number "ii->second" is replaced by a row with two entries, 1·
    // on the diagonal and -1 on the "ii->first"-th entry. Also the right hand·
    // side is changed to be zero

    // loop over all entries in row "ii->first" of A-matrices
    for (int i = rowPtr[ii->first]; i < rowPtr[1 + ii->first]; i++)
    {
      if (entries[i] != 0.0 || p)
        (new_rows[ii->first])[colPtr[i]] += entries[i];
    } // up to here this row would simply be copied.
    // loop over all entries in row "ii->second" of A-matrices
    for (int i = rowPtr[ii->second]; i < rowPtr[1 + ii->second]; i++)
    {
      if (entries[i] != 0.0 || p)
         (new_rows[ii->first])[colPtr[i]] += entries[i];
    }
    if (stokesMat && p)
    {
      if (do_once == 0)
      {
      this->get_rhs()[ii->first] += this->get_rhs()[ii->second];
      this->get_rhs()[ii->first + this->get_rhs().length(0)] += this->get_rhs()[ii->second + this->get_rhs().length(0)];
      }
    }

    if(stokesMat)
    {
      (new_rows[ii->second])[ii->second] = 1.0; // diagonal
      (new_rows[ii->second])[ii->first] = -1.0; // coupling
      // here two entries for the right hand side are handled.
      // One would be enough, but which one depends on which matrix this is
      // (A11 or A22).

      this->get_rhs()[ii->second] = 0;
      this->get_rhs()[ii->second + this->get_rhs().length(0)] = 0;
    }
   else
      (new_rows[ii->second]); // sets the entire row to zero
  }
  mat->changeRows(new_rows);
}

/* ************************************************************************* */
template <int d>
void NavierStokes<d>::checkPeriodicDOFs()
{
  std::map<int, int>::iterator it;
  double const * const vals = this->get_solution().get_entries();
  for(it = periodic_dofs.begin(); it != periodic_dofs.end(); ++it)
  {
    double diff = vals[it->first] - vals[it->second];
    if(std::abs(diff) > 1e-6)
      Output::warn("NavierStokes<d>::checkPeriodicDOFs", "ERROR, dofs provided "
                   "by NavierStokes<d>::periodic_dofs ar not periodic. The "
                   "first problematic dof pair and the difference of function "
                   "values are ", it->first, "\t", it->second, "\t",
                   setprecision(10), diff);
  }
}

/* ************************************************************************* */
#ifdef __3D__
template class NavierStokes<3>;
#else
template class NavierStokes<2>;
#endif
