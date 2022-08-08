#include "ConvectionDiffusion.h"
#include "ConvDiff.h"
#include "Database.h"
#include "Domain.h"
#include "LocalAssembling.h"
#include "Multigrid.h"
#include "MainUtilities.h"
#include "BaseCell.h"
#include "Assemble_DG.h"
#include "Point.h"
#ifdef __2D__
 #include "Assemble2D.h"
 #include "SquareMatrix2D.h"
 #include "AuxParam2D.h"
 #include "Example_CD2D.h"
#else
 #include "Assemble3D.h"
 #include "SquareMatrix3D.h"
 #include "AuxParam3D.h"
 #include "Example_CD3D.h"
#endif
#ifdef _MPI
 #include "ParFECommunicator3D.h"
 #include "MumpsWrapper.h"
 #include <cmath>
#endif
#include <algorithm>
#include <Utilities.h>

/* ************************************************************************* */
template<int d>
ParameterDatabase ConvectionDiffusion<d>::default_cd_database(bool complete)
{
  Output::print<5>("creating a default ConvectionDiffusion parameter database");
  // we use a ParMooN default database because this way these parameters are
  // available in the default CD2D database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("ConvectionDiffusion parameter database");

  if(complete)
  {
    db.merge(TDomain::default_domain_parameters());
    if(d == 3)
      db.add_nested_database(TDomain::default_sandwich_grid_parameters());
    db.merge(Example_CD::default_example_database());
  }
  // a default output database - needed here as long as there's no class handling the output
  db.merge(ParameterDatabase::default_output_database());
  // default local assembling database
  db.merge(LocalAssembling<d>::default_local_assembling_database());
  if(complete)
  {
    db.merge(CheckpointIO::default_checkpoint_io_database());
    db.merge(Solver<>::default_solver_database());
#ifdef __2D__
    db.merge(TSlopeLimiter::default_slope_limiter_database());
#endif
  }
  // Note that the DG database cannot be in the complete case, since the DG
  // method is not a class yet. The reason for this behaviour is that the CD
  // class only saves the parameters needed by itself and forwards the other
  // parameters to other classes. The rest is deleted if not saved explicitly.
  // Hence, it cannot be in the complete case yet.
  db.merge(default_dg_database());

  db["problem_type"] = 1; // means stationary convection-diffusion
  return db;
}

/* ************************************************************************* */
template<int d>
ConvectionDiffusion<d>::SystemPerGrid::SystemPerGrid(const Example_CD& example,
    const TCollection& coll, int ansatz_order)
: fe_space(new FESpace(&coll, "space", example.get_bc(0), ansatz_order))
{
#ifdef _MPI
  fe_space->get_communicator().print_info();
#endif // _MPI
#ifdef __3D__
  matrix = BlockFEMatrix::CD3D(fe_space);
#else
  matrix = BlockFEMatrix::CD2D(fe_space);
#endif

  rhs = BlockVector(this->matrix, true);
  solution = BlockVector(this->matrix, false);
  fe_function = FEFunction(fe_space, "c", solution.get_entries());
}

/* ************************************************************************* */
template<int d>
ConvectionDiffusion<d>::SystemPerGrid::SystemPerGrid(const SystemPerGrid& other)
 : fe_space(other.fe_space), matrix(other.matrix), rhs(other.rhs),
   solution(other.solution)
{
  // the fe functions must be newly created, because copying would mean 
  // referencing the BlockVector in 'other'.
  fe_function = FEFunction(fe_space, "c", solution.get_entries());
}


/* ************************************************************************* */
template<int d>
ConvectionDiffusion<d>::ConvectionDiffusion(const TDomain& domain,
                                            const ParameterDatabase& param_db)
 : ConvectionDiffusion<d>(domain, param_db, Example_CD(param_db))
{
}

/* ************************************************************************* */
template<int d>
ConvectionDiffusion<d>::ConvectionDiffusion(const TDomain& domain,
                                            const ParameterDatabase& param_db,
                                            const Example_CD& example_cd)
 : systems(), example(example_cd), db(default_cd_database()),
   outputWriter(param_db), checkpoint_io(param_db), solver(param_db),
   limiter(param_db), errors()
{
  db.merge(param_db, false); // update this database with given values
  set_parameters();
  // The construction of the members differ, depending on whether
  // a multigrid solver will be used or not.
  bool usingMultigrid = solver.is_using_multigrid();
  int ansatz_order = TDatabase::ParamDB->ANSATZ_ORDER;
  
  auto collections = domain.get_grid_collections();
  if(!usingMultigrid)
  {
    // Get the collection on the finest grid
    TCollection& cellCollection = *collections.front();
    // create finite element space and function, a matrix, rhs, and solution
    systems.emplace_back(example, cellCollection, ansatz_order);
  }
  else
  {
    // we are using multigrid
    size_t n_levels = collections.size();
    size_t desirec_multigrid_levels = param_db["multigrid_n_levels"];
    if(desirec_multigrid_levels > n_levels)
      ErrThrow("Not enough collections (", n_levels, ") to use ",
               desirec_multigrid_levels, " multigrid levels!");
    // remove not needed coarser grid from list of collections
    for(size_t i = desirec_multigrid_levels; i < n_levels; ++i)
    {
      collections.pop_back();
    }

    auto mg = this->solver.get_multigrid();
    // Construct systems per grid and store them, finest level first
    std::list<BlockFEMatrix*> matrices;
    for(auto coll : collections)
    {
      Output::print("creating SystemPerGrid object, n_cells = ", coll->GetN_Cells());
      systems.emplace_back(example, *coll, ansatz_order);
      //prepare input argument for multigrid object
      matrices.push_front(&systems.back().matrix);
    }
    mg->initialize(matrices);
  }
  
  outputWriter.add_fe_function(&this->get_function());
  checkpoint_io.read(this->get_solution());
  output_problem_size_info();
  if (db["space_discretization_type"].is("dg") &&
      !this->get_space()->is_discontinuous())
  {
    // To me, it is not clear if this even makes sense at all. Since it is not
    // tested yet, we just exclude it in general even though it might work. But
    // you can try it out if you want to.
    ErrThrow("A DG discretization with continuous elements is not allowed.");
  }
  if (db["space_discretization_type"].is("galerkin") &&
      this->get_space()->is_discontinuous())
  {
    Output::warn("ConvectionDiffusion.C", "Your space is not continuous but ",
        "you choose the standard Galerkin discretization. If you want to use a",
        " DG discretization choose in your database:",
        "\tspace_discretization_type: dg");
  }
}

/* ************************************************************************* */
template<int d>
void ConvectionDiffusion<d>::set_parameters()
{
  //set problem_type to CD if not yet set
  if(!db["problem_type"].is(1))
  {
    if (db["problem_type"].is(0))
    {
      db["problem_type"] = 1;
    }
    else
    {
      Output::warn<2>("The parameter problem_type doesn't correspond to CD."
          "It is now reset to the correct value for CD (=1).");
      db["problem_type"] = 1;
    }
  }
  
  if(TDatabase::ParamDB->ANSATZ_ORDER == 0)
  {
    ErrThrow("Ansatz order 0 is not used in convection diffusion "
             "reaction problems! (Vanishing convection and diffusion term).");
  }
}

/* ************************************************************************* */
template<int d>
void ConvectionDiffusion<d>::output_problem_size_info() const
{
  // print some useful information
  auto& space = *this->systems.front().fe_space;
  double hMin, hMax;
  auto coll = space.GetCollection();
  coll->GetHminHmax(&hMin, &hMax);
  bool delaunay;
  int hanging_nodes = space.get_n_hanging();
  
  delaunay = coll -> IsDelaunay();
  Output::print<1>("N_Cells    : ", setw(13), coll->GetN_Cells());
  Output::print<1>("h(min, max): ", setw(13), hMin, " ", setw(13), hMax);
  Output::print<1>("dofs all   : ", setw(13), space.get_n_dof());
  Output::print<1>("dof active : ", setw(13), space.get_n_active());
  Output::print<1>("hanging nodes: ", setw(13), hanging_nodes);

  if(hanging_nodes != 0)
      Output::print<1>("NOTE: If the initial mesh is Delaunay then so would",
                        " be the refinements");
  if(delaunay)
    Output::print<1>("The mesh is Delaunay");
  else
    Output::print<1>("The mesh is not Delaunay");
}

/* ************************************************************************* */
template<int d>
void ConvectionDiffusion<d>::assemble(
  std::function<void(LocalAssembling<d>& la)> modify_la)
{  
  LocalAssembling_type laType = LocalAssembling_type::ConvDiff;
  // this loop has more than one iteration only in case of multigrid
  for(auto & s : systems)
  {
    std::vector<const FEFunction*> feFunctionPtr = {&s.fe_function};

    // create a local assembling object which is needed to assemble the matrix
    LocalAssembling<d> laObject(this->db, laType, feFunctionPtr,
                                example.get_coeffs());
    if(modify_la)
      modify_la(laObject);
    call_assembling_routine(s, laObject, false);
  }
}

/* ************************************************************************* */
template<int d>
void ConvectionDiffusion<d>::solve()
{
  double t = GetTime();
  SystemPerGrid& s = this->systems.front();
#ifndef _MPI
  this->solver.solve(s.matrix, s.rhs, s.solution);
#else
  if(this->solver.get_db()["solver_type"].is("direct"))
  {
    MumpsWrapper mumps_wrapper(s.matrix);
    mumps_wrapper.solve(s.rhs, s.solution);
  }
  else
    this->solver.solve(s.matrix, s.rhs, s.solution);
#endif
  t = GetTime() - t;
  Output::print<3>("solving time of a ConvectionDiffusion<", d, "> problem: ",
                   t, " seconds");

  // Apply limiter
#ifdef __2D__
  limiter.limit_function(&this->get_function());
#endif
}

/* ************************************************************************* */
template<int d>
void ConvectionDiffusion<d>::output(int i)
{
  int my_rank = 0;
#ifdef _MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  // computing errors as well as writing vtk files requires a minimum 
  // consistency level of 1
  this->systems.front().fe_space->get_communicator().consistency_update(
    this->systems.front().solution.get_entries(), 1);
#endif

  const FEFunction & fe_function = this->systems.front().fe_function;
  if(db["output_compute_minmax"])
  {
    fe_function.PrintMinMax();
  }

  // write solution
  if(i < 0)
    outputWriter.write();
  else
    outputWriter.write(i);
  checkpoint_io.write(this->get_solution());

  // measure errors to known solution
  // If an exact solution is not known, it is usually set to be zero, so that
  // in such a case here only integrals of the solution are computed.
  if(db["output_compute_errors"])
  {
    std::array<double, ConvectionDiffusion<d>::n_errors> errors;
#ifdef __3D__
    TAuxParam3D aux;
    MultiIndex3D AllDerivatives[4] = { MultiIndex3D::D000, MultiIndex3D::D100,
                                       MultiIndex3D::D010, MultiIndex3D::D001 };
#else
    TAuxParam2D aux;
    MultiIndex2D AllDerivatives[3] = { MultiIndex2D::D00, MultiIndex2D::D10,
                                       MultiIndex2D::D01 };
#endif
    const FESpace* space = this->systems.front().fe_space.get();

#ifdef __3D__
    fe_function.GetErrors(example.get_exact(0), d+1, AllDerivatives,
                          ConvectionDiffusion<d>::n_errors,
                          conv_diff_l2_h1_linf_error<d>, example.get_coeffs(),
                          &aux, 1, &space, errors.data());
#else
    fe_function.GetErrors(example.get_exact(0), d+1, AllDerivatives,
                          ConvectionDiffusion<d>::n_errors,
                          conv_diff_l2_h1_linf_error<d>, example.get_coeffs(),
                          &aux, 1, &space, errors.data(), false,[](const TBaseCell*, int){return false;}, db);
#endif

#ifdef _MPI
    // memory for global (across all processes) error
    std::array<double, ConvectionDiffusion<d>::n_errors> errorsReduced;
    MPI_Reduce(errors.data(), errorsReduced.data(),
               ConvectionDiffusion<d>::n_errors-1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&errors[ConvectionDiffusion<d>::n_errors-1],
               &errorsReduced[ConvectionDiffusion<d>::n_errors-1], 1,
               MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    for(int i = 0; i < ConvectionDiffusion<d>::n_errors-1; i++) // exclude l_inf
      errors[i] = std::sqrt(errorsReduced[i]);
    errors[ConvectionDiffusion<d>::n_errors-1] 
      = errorsReduced[ConvectionDiffusion<d>::n_errors-1];
#endif
    // copy local variable to member variable
    std::copy(errors.begin(), errors.end(), this->errors.begin());


    //print errors
    if(my_rank == 0)
    {
      Output::print<1>("L2     : ", setprecision(14), errors[0]);
      Output::print<1>("H1-semi: ", setprecision(14), errors[1]);
      Output::print<1>("SD     : ", setprecision(14), errors[2]);
      if(db["space_discretization_type"].is("dg"))
        Output::print<1>("DG     : ", setprecision(14), errors[3]);
      Output::print<1>("L_inf  : ", setprecision(14), errors[4]);
    }
  } // if(this->db["compute_errors"])

  // do postprocessing
  example.do_post_processing(*this);
}

/** ************************************************************************ */
template<int d>
double ConvectionDiffusion<d>::get_L2_error() const
{
  return this->errors[0];
}

/** ************************************************************************ */
template<int d>
double ConvectionDiffusion<d>::get_H1_semi_error() const
{
  return this->errors[1];
}

/** ************************************************************************ */
template<int d>
double ConvectionDiffusion<d>::get_SD_error() const
{
  return this->errors[2];
}

/** ************************************************************************ */
template<int d>
double ConvectionDiffusion<d>::get_DG_error() const
{
  return this->errors[3];
}

/** ************************************************************************ */
template<int d>
double ConvectionDiffusion<d>::get_L_inf_error() const
{
  return this->errors[4];
}

/** ************************************************************************ */
template<int d>
void ConvectionDiffusion<d>::call_assembling_routine(SystemPerGrid& s,
						     LocalAssembling<d>& local_assem, bool assemble_dirichlet_rows)
{
  // assemble the system matrix with given local assembling, solution and rhs
  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  int n_fe_spaces = 1;
  const FESpace * fe_space = s.fe_space.get();
  auto * boundary_conditions = fe_space->get_boundary_condition();
  int n_square_matrices = 1;
  int n_rect_matrices = 0;
  double * rhs_entries = s.rhs.get_entries();
  int n_rhs = 1;
  auto * non_const_bound_value = example.get_bd(0);

  //fetch stiffness matrix as block
  auto blocks = s.matrix.get_blocks_uniquely();
  SquareMatrixD* block = reinterpret_cast<SquareMatrixD*>(blocks.at(0).get());

  // reset right hand side and matrix to zero
  s.rhs.reset();
  block->reset();
  // and call the assembling method
#ifdef __3D__
  Assemble3D(
#else
  Assemble2D(
#endif
    n_fe_spaces, &fe_space, n_square_matrices, &block,
    n_rect_matrices, nullptr, n_rhs, &rhs_entries, &fe_space,
    &boundary_conditions, &non_const_bound_value, local_assem,
    assemble_dirichlet_rows);
  // copy Dirichlet values from rhs to solution vector (this is not really
  // necessary in case of a direct solver)
  s.solution.copy_nonactive(s.rhs);

  // Assemble DG boundary terms
  if ( db["space_discretization_type"].is("dg") )
  {
    Assemble_DG<d>( local_assem.GetCoeffFct(), n_fe_spaces, &fe_space,
      n_square_matrices, &block, n_rect_matrices, nullptr, n_rhs,
      &rhs_entries, &boundary_conditions, &non_const_bound_value, db, "CD" );
  }
}


#ifdef __3D__
template class ConvectionDiffusion<3>;
#else
template class ConvectionDiffusion<2>;
#endif
