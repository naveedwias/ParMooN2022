#include <Database.h>
#include <Darcy.h>
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

template <int d>
ParameterDatabase Darcy<d>::default_darcy_database(bool complete)
{
  Output::print<5>("creating a default Darcy parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default Darcy database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("Darcy parameter database");
  
  if(complete)
  {
    db.merge(TDomain::default_domain_parameters(), true);
    if(d == 3)
      db.add_nested_database(TDomain::default_sandwich_grid_parameters());
    db.merge(Example_Darcy::default_example_database());
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

  return db;
}


/** ************************************************************************ */
template <int d>
Darcy<d>::System_per_grid::System_per_grid(
  const Example_Darcy& example, const TCollection& coll,
  std::pair<int,int> velocity_pressure_orders)
 : velocity_space(new FESpace(&coll, "u", example.get_bc(0),
                              velocity_pressure_orders.first)),
   pressure_space(new FESpace(&coll, "p", example.get_bc(1),
                              velocity_pressure_orders.second))
{
#ifdef __3D__
  matrix = BlockFEMatrix::Darcy3D(this->velocity_space, this->pressure_space);
#else
  matrix = BlockFEMatrix::Darcy2D(this->velocity_space, this->pressure_space);
#endif

  rhs = BlockVector(this->matrix, true);
  solution = BlockVector(this->matrix, false);

  u = FEFunction(this->velocity_space, "u", this->solution.block(0));
  p = FEFunction(this->pressure_space, "p", this->solution.block(1));
}

/** ************************************************************************ */
template <int d>
Darcy<d>::Darcy(const TDomain& domain, const ParameterDatabase& param_db)
 : Darcy<d>(domain, param_db, Example_Darcy(param_db))
{
}

/** ************************************************************************ */
template <int d>
Darcy<d>::Darcy(const TDomain& domain, const ParameterDatabase& param_db,
                 const Example_Darcy ex)
 : systems(), example(ex), db(default_darcy_database(false)),
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
    ErrThrow("The multigrid implementation for Darcy type problems is not yet "
             "finished. You have to pick a different solver/preconditioner.");
    // Construct multigrid object
    auto mg = solver.get_multigrid();
    if(mg->is_using_mdml())
    {
      ErrThrow("mdml not supported for Darcy type problems");
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
  
  Output::stat<1>("Darcy", "Mesh data and problem size in ", d, "D");
  Output::dash<1>("cells                        : ", setw(5), n_cells);
  Output::dash<1>("dof velocity (vector-valued) : ", setw(5), n_u);
  Output::dash<1>("active dof velocity          : ", setw(5), n_u_active);
  Output::dash<1>("dof pressure                 : ", setw(5), n_p);
  Output::dash<1>("dof all                      : ", setw(5), n_dof);
}

/** ************************************************************************ */
template <int d>
void Darcy<d>::set_parameters()
{
  // check if given velocity space is supported
  switch(TDatabase::ParamDB->VELOCITY_SPACE)
  {
    case 1000: case 1001: case 1002: case 1003: case 1011: case 1012: case 1013:
      break;
    default:
      ErrThrow("unknown velocity space for Darcy");
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
        ErrThrow("unknown velocity space for Darcy");
        break;
    }
  }
}

/** ************************************************************************ */
template <int d>
void Darcy<d>::assemble()
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
    LocalAssembling<d> la(this->db, LocalAssembling_type::Darcy, fe_functions,
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
    std::array<BoundaryValuesFunction*, n_fe_spaces> non_const_bound_values;
    non_const_bound_values[0] = example.get_bd()[0];
    non_const_bound_values[1] = example.get_bd()[1];
    
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
  }
  // copy Dirichlet values from rhs to solution vector (this is not really 
  // necessary in case of a direct solver)
  this->systems.front().solution.copy_nonactive(this->systems.front().rhs);
} // void Darcy3D::Assemble

/** ************************************************************************ */
template <int d>
void Darcy<d>::solve()
{
  double t = GetTime();
  System_per_grid& s = this->systems.front();
  
  this->solver.solve(s.matrix, s.rhs, s.solution);
  
  if(s.matrix.pressure_projection_enabled())
    s.p.project_into_L20();
  
  t = GetTime() - t;
  Output::print<2>(" solving a Darcy", d, "D problem done in ", t, " seconds");
}

/** ************************************************************************ */
template <int d>
void Darcy<d>::output(int i)
{
  using ErrorMethod = typename FEFunction::ErrorMethod;
  
  System_per_grid & s = this->systems.front();
  
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
}

/** ************************************************************************ */
template <int d>
const BlockVector & Darcy<d>::get_solution() const
{
  return this->systems.front().solution;
}
    
/** ************************************************************************ */
template <int d>
BlockVector & Darcy<d>::get_solution()
{
  return this->systems.front().solution;
}

/** ************************************************************************ */
template <int d>
double Darcy<d>::getL2VelocityError() const
{
  return this->errors[0];
}

/** ************************************************************************ */
template <int d>
double Darcy<d>::getL2DivergenceError() const
{
  return this->errors[1];
}

/** ************************************************************************ */
template <int d>
double Darcy<d>::getH1SemiVelocityError() const
{
  return this->errors[2];
}

/** ************************************************************************ */
template <int d>
double Darcy<d>::getL2PressureError() const
{
  return this->errors[3];
}

/** ************************************************************************ */
template <int d>
double Darcy<d>::getH1SemiPressureError() const
{
  return this->errors[4];
}

#ifdef __3D__
template class Darcy<3>;
#else
template class Darcy<2>;
#endif


