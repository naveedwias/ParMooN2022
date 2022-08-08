/**
 * @file New implementation of a multigrid object, which holds the necessary
 * grid information for executing a multigrid iteration.
 *
 * @date 2016/05/10
 * @author Clemens Bartsch
 */

#include <BlockVector.h>
#include <GridTransfer.h>
#include <Multigrid.h>
#include <MooNMD_Io.h>
#include <ParameterDatabase.h>

#ifdef _MPI
#include <ParFECommunicator3D.h>
#endif

SmootherCode string_to_smoother_code(const std::string& code)
{
  if(code == std::string("direct_solve"))
    return SmootherCode::DIRECT_SOLVE;
  else if(code == std::string("jacobi"))
    return SmootherCode::JACOBI;
  else if(code == std::string("sor"))
    return SmootherCode::SOR;
  else if(code == std::string("ssor"))
    return SmootherCode::SSOR;
  else if(code == std::string("nodal_vanka"))
    return SmootherCode::NODAL_VANKA;
  else if(code == std::string("cell_vanka"))
    return SmootherCode::CELL_VANKA;
  else if(code == std::string("patch_vanka"))
    return SmootherCode::PATCH_VANKA;
  else if(code == std::string("cell_vanka_jacobi"))
    return SmootherCode::CELL_VANKA_JACOBI;
  else if(code == std::string("nodal_vanka_store"))
    return SmootherCode::NODAL_VANKA_STORE;
  else if(code == std::string("cell_vanka_store"))
    return SmootherCode::CELL_VANKA_STORE;
  else if(code == std::string("patch_vanka_store"))
    return SmootherCode::PATCH_VANKA_STORE;
  else
  {
    Output::warn("SmootherCode", "The string ", code,
                 " does not equal a smoother code. "
                 "Defaulting to DIRECT_SOLVE");
    return SmootherCode::DIRECT_SOLVE;
  }

}

MultigridType string_to_multigrid_type(const std::string& code)
{
  if(code == std::string("standard"))
      return MultigridType::STANDARD;
  else if(code == std::string("mdml"))
      return MultigridType::MDML;
  else
  {
    Output::warn("MultigridType", "The string ", code,
                   " does not equal a multigrid type. "
                   "Defaulting to STANDARD");
    return MultigridType::STANDARD;
  }
}

Multigrid::Multigrid(const ParameterDatabase& param_db)
: db(default_multigrid_database()), damp_factors_(), 
  type_(MultigridType::STANDARD)
{
  this->db.merge(param_db, false);
  this->type_ = string_to_multigrid_type(db["multigrid_type"]);
  
  size_t n_geo_levels = this->get_n_geometric_levels();

  //Determine the number of algebraic levels (depth of space hierarchy)
  if(type_ == MultigridType::MDML)
    n_algebraic_levels_ = n_geo_levels + 1;
  else
    n_algebraic_levels_ = n_geo_levels;

  if(n_algebraic_levels_ <= 1)
    ErrThrow("for multigrid you need at least two multigrid levels, ",
             n_algebraic_levels_, " is plain not enough.");

  damp_factors_.resize(n_algebraic_levels_, db["multigrid_correction_damp_factor"]);
  
  // Set up the cycle control.
  std::string cycle_str = db["multigrid_cycle_type"];
  control_ = CycleControl(cycle_str, n_algebraic_levels_);

  n_pre_smooths_ = db["multigrid_n_pre_smooth"];

  coarse_n_maxit = db["multigrid_coarse_max_n_iterations"];

  coarse_epsilon = db["multigrid_coarse_residual"];

  n_post_smooths_ = db["multigrid_n_post_smooth"];

  //Send the coarse grid stopwatch to sleep.
  coarse_grid_timer.reset();

}

void Multigrid::initialize(const std::list<BlockFEMatrix*>& matrices)
{
  // Create the levels and collect them in a list
  if(matrices.size() != n_algebraic_levels_)
  {
    ErrThrow("Multigrid::initialize", " The number of algebraic multigrid levels "
                    "was set to ", n_algebraic_levels_, ", but there are ",
                    matrices.size(), " matrices given. Please check your program!");
  }
  auto coarsest = matrices.front();
  for(auto mat : matrices)
  {
    SmootherCode sm;
    if(mat == coarsest)
      sm = string_to_smoother_code(db["multigrid_smoother_coarse"]);
    else
      sm = string_to_smoother_code(db["multigrid_smoother"]);

    // A database, can be filled with parameters the level might need.
    ParameterDatabase level_db(std::string("multigrid level database"));
    //TODO only the parameters required by sm should go in that level database!
    level_db.add(Parameter(db["multigrid_vanka_damp_factor"]));
    level_db.add(Parameter(db["jacobi_damp"]));
    level_db.add(Parameter(db["sor_omega"]));
#ifdef _MPI
    level_db.add(Parameter(db["sor_parallel_type"]));
#endif

    levels_.push_back(MultigridLevel(mat, sm, level_db));
  }
}


void Multigrid::set_finest_sol(const BlockVector& bv)
{
  //TODO Eventually here a call to "GridTransfer::RestrictFunctionRepeatedly"
  // is necessary - so that for step length control or Braess-Sarazin an initial
  // residual can be calculated on each level
  levels_.back().solution_ = bv;
}

void Multigrid::set_finest_rhs(const BlockVector& bv)
{
  levels_.back().rhs_ = bv;
}

const BlockVector& Multigrid::get_finest_sol()
{
  return levels_.back().solution_;
}

void Multigrid::cycle()
{
  Output::info<4>("Multigrid", "Starting multigrid cycle of type ", (int) control_.get_type());
  size_t n_steps = control_.get_n_steps();
  size_t finest = levels_.size() - 1;

  //Start with 0 solution TODO enable other start solution!
  levels_.at(finest).solution_ = 0.0;
  //...but copy non-actives
  levels_.at(finest).solution_.copy_nonactive(levels_.at(finest).rhs_);

  int next_level = finest;
  for(size_t i = 1; i <= n_steps; ++i)
  {
    next_level = cycle_step(i, next_level);
  }

}

void Multigrid::update()
{
  const auto& coarsest = levels_.front();
  for( auto& lvl : levels_)
  {
    if(&lvl == &coarsest) //measure time spent on coarsest grid
    {
      coarse_grid_timer.start();
    }

    lvl.update_smoother();

    if(&lvl == &coarsest) //stop measure time spent on coarsest grid
    {
      coarse_grid_timer.stop();
    }
  }
}

void Multigrid::print_coarse_grid_time_total() const
{
  coarse_grid_timer.print_total_time("coarse grid (total)");
}

int Multigrid::cycle_step(size_t step, size_t level)
{
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#else
  int my_rank = 0;
#endif
  if (my_rank == 0)
    Output::info<4>("Multigrid", "Doing cycle step no ", step, " on level ", level);

  //catch the coarsest level case
  size_t coarsest = 0;
  if(level == coarsest)
  {//we're on coarsest level
    double res = 1e10;
    if ((my_rank == 0) && (Output::getVerbosity()>=4))
    {
      Output::info<4>("COARSE SOLVE", "Level ", level);
      levels_.at(level).calculate_defect();
      res = levels_.at(level).residual_;
      Output::info<4>("COARSE GRID bef.", "Level ", level, " res: ", res);
    }
    else
      levels_.at(level).calculate_defect();
    
    //start coarse grid time measurement
    coarse_grid_timer.start();

    for(size_t i = 0; i < this->coarse_n_maxit && res > coarse_epsilon ; ++i)
    {
      levels_.at(coarsest).apply_smoother();
      levels_.at(coarsest).calculate_defect();
      res = levels_.at(coarsest).residual_;
      if (my_rank == 0)
        Output::dash<4>("Coarse Grid Iteration ", i+1, " res: ", res);
    }

    //stop coarse grid time measurement
    coarse_grid_timer.stop();

    update_solution_in_finer_grid(level);
    return level + 1;
  }
  else
  {
    bool coming_from_below = (control_.get_direction(step-1) == MGDirection::Up);
    bool going_up = (control_.get_direction(step) == MGDirection::Up);
    bool going_down = (control_.get_direction(step) == MGDirection::Down);

    if(coming_from_below)
    {
      // compute and print more information, for debugging
      if ((my_rank == 0) && (Output::getVerbosity()>=4))
      {   
	double res;
	levels_.at(level).calculate_defect();
	res = levels_.at(level).residual_;
        Output::info<4>("POST SMOOTH bef.", "Level ", level, " res: ", res);
      }
      for(size_t i = 0; i < n_post_smooths_; ++i)
        levels_.at(level).apply_smoother(); //post smoothing
      // compute and print more information, for debugging
      if ((my_rank == 0) && (Output::getVerbosity()>=4))
      {   
	double res;	
	levels_.at(level).calculate_defect();
	res = levels_.at(level).residual_;
        Output::info<4>("POST SMOOTH done", "Level ", level, " res: ", res);
      }
    }
    if(going_down)
    {
      // compute and print more information, for debugging
      if ((my_rank == 0) && (Output::getVerbosity()>=4))
      {
        double res;
	levels_.at(level).calculate_defect();
	res = levels_.at(level).residual_;
        Output::info<4>("PRE SMOOTH bef.", "Level ", level, " res: ", res);
      }
      for(size_t i = 0; i < n_pre_smooths_; ++i)
        levels_.at(level).apply_smoother(); //pre smoothing
      levels_.at(level).calculate_defect(); //defect calculation
      set_solution_in_coarser_grid_to_zero(level); //set start iterate on coarser level to 0
      update_rhs_in_coarser_grid(level);
      // compute and print more information, for debugging
      if ((my_rank == 0) && (Output::getVerbosity()>=4))
      {
        double res;
	res = levels_.at(level).residual_;
        Output::info<4>("PRE SMOOTH done", "Level ", level, " res: ", res);
      }
      return level - 1;
    }
    else if(going_up)
    {
      update_solution_in_finer_grid(level);
      return level + 1;
    }
    else
    {//this is the last step - just return!
      if (my_rank == 0)
        Output::info<4>("Multigrid", "Cycle finished");
      return -1;
    }

  }
}

void Multigrid::update_rhs_in_coarser_grid(size_t lvl)
{
  MultigridLevel& level_current = levels_.at(lvl);
  MultigridLevel& level_coarse = levels_.at(lvl - 1);

  const BlockFEMatrix& matrix_current = *level_current.matrix_;
  BlockVector& defect_current = level_current.defect_;

  const BlockFEMatrix& matrix_coarse = *level_coarse.matrix_;
  BlockVector& rhs_coarse = level_coarse.rhs_;

  for(size_t i = 0 ; i < matrix_current.get_n_cell_rows() ; ++i)
  {
#ifdef __2D__
    const TFESpace2D& space_current = *matrix_current.get_row_space(i);
    const TFESpace2D& space_coarse = *matrix_coarse.get_row_space(i);

#endif
#ifdef __3D__
    const TFESpace3D& space_current = *matrix_current.get_row_space(i);
    const TFESpace3D& space_coarse = *matrix_coarse.get_row_space(i);
#endif
    double* defect_current_entries = defect_current.block(i);
    int size_current_defect = defect_current.length(i);

    double* rhs_coarse_entries = rhs_coarse.block(i);
    int size_coarse_rhs = rhs_coarse.length(i);

    // Restrict current level's defect to coarser grid and write
    // that restriction into coarser level's right hand side.
#ifdef _MPI
  const TParFECommunicator3D& comm_fine = space_current.get_communicator();
  comm_fine.consistency_update(defect_current_entries, 1); //restore level 1 consistency
#endif
    GridTransfer::DefectRestriction(space_coarse, space_current,
                                    rhs_coarse_entries, size_coarse_rhs,
                                    defect_current_entries, size_current_defect);
#ifdef _MPI
    const TParFECommunicator3D& comm_coarse = space_coarse.get_communicator();
    comm_coarse.CommUpdateReduce(rhs_coarse_entries); //call CommUpdateReduce TODO Why exactly???
#endif
  }

  // Nuke all non-actives.
  level_coarse.rhs_.ResetNonActive();
}

void Multigrid::update_solution_in_finer_grid(size_t lvl)
{
  MultigridLevel& level_fine = levels_.at(lvl + 1);
  MultigridLevel& level_current = levels_.at(lvl);

  const BlockFEMatrix& matrix_fine = *level_fine.matrix_;
  BlockVector& solution_fine = level_fine.solution_;

  const BlockFEMatrix& matrix_current = *level_current.matrix_;
  BlockVector& solution_current = level_current.solution_;

  BlockVector copy_sol_fine = solution_fine; //working copy!

  for(size_t i = 0 ; i < matrix_current.get_n_cell_rows() ; ++i)
  {
#ifdef __2D__
    const TFESpace2D& space_fine = *matrix_fine.get_row_space(i);
    const TFESpace2D& space_current = *matrix_current.get_row_space(i);
#endif
#ifdef __3D__
    const TFESpace3D& space_fine = *matrix_fine.get_row_space(i);
    const TFESpace3D& space_current = *matrix_current.get_row_space(i);
#endif

    double* solution_prolongation = copy_sol_fine.block(i);
    int size_solution_prolongation = copy_sol_fine.length(i);

    double* sol_cur_entries = solution_current.block(i);
    int size_sol_cur = solution_current.length(i);

#ifdef _MPI
    const TParFECommunicator3D& comm_current = space_current.get_communicator();
    comm_current.consistency_update(sol_cur_entries, 1); //restore level 1 consistency
#endif

    //Do the prolongation and write its result into solution_prolongation
    GridTransfer::Prolongate(
        space_current, space_fine,
        sol_cur_entries, size_sol_cur,
        solution_prolongation, size_solution_prolongation);

    //Update the actual solution_fine_entries by damped addition
    double* sol_fine_entries = solution_fine.block(i);
    double damp = this->damp_factors_[lvl];

    int n_non_actives = solution_fine.n_non_actives(i);
    for( int j=0 ; j < size_solution_prolongation - n_non_actives; ++j )
    {
      sol_fine_entries[j] += damp * solution_prolongation[j];
    }

  }

}

void Multigrid::set_solution_in_coarser_grid_to_zero(size_t lvl)
{
  levels_.at(lvl-1).solution_ = 0.0;
}

ParameterDatabase Multigrid::default_multigrid_database()
{
  Output::print<5>("creating a default multigrid parameter database");
  ParameterDatabase db("default multigrid database");

  db.add<size_t>("multigrid_n_levels", 2,
         "Determine how many levels to use in the grid hierarchy. Note that "
         "this determines the number of geometrical grids to use, so if you use mdml "
         "the multigrid cycle will be multigrid_n_levels + 1 levels deep, since the "
         "finest geometrical grid is used twice.", 0, 5);

  db.add("multigrid_type", "standard",
         "The type of multigrid algorithm to apply. Besides "
         "the standard approach a 'multiple discretization multilevel'"
         " approach is available (for those systems where it is "
         "already implemented). See e.g. John et al. 2002. ",
         {"standard", "mdml"});

  db.add("multigrid_cycle_type", "V",
         "The recursion type how to traverse the multigrid levels. "
         "So far the three standard cycle V, W and F are implemented."
         , {"V", "W", "F"});

  db.add("multigrid_smoother", "nodal_vanka",
         "The smoother to use on all but the coarsest level. You should take "
         "care, that the smoother you chose fits your problem type, e.g. Vanka "
         "smoothers are best fitted for saddle point problems.",
         {"jacobi", "sor",
          "nodal_vanka", "cell_vanka", "patch_vanka",
          "cell_vanka_jacobi",
          "nodal_vanka_store", "cell_vanka_store", "patch_vanka_store", "no_smoother"});

  db.add("multigrid_smoother_coarse", "direct_solve",
         "The smoother to use on the coarsest level. You should take care, "
         "that the smoother you chose fits your problem type, e.g. Vanka "
         "smoothers are best fitted for saddle point problems.",
         {"direct_solve", "jacobi", "sor",
          "nodal_vanka", "cell_vanka", "patch_vanka",
          "cell_vanka_jacobi",
          "nodal_vanka_store", "cell_vanka_store", "patch_vanka_store", "no_smoother"});

  db.add("multigrid_correction_damp_factor", 1.0,
         "The damping factor which is used when applying the coarse grid "
         "correction to a finer grid. A factor of 1.0 means no damping, a "
         "factor of 0.0 means nothing changes.",
         0.0,1.0);

  db.add<size_t>("multigrid_n_pre_smooth", 1,
                 "The number of smoothing steps to apply per level before "
                 "going down to the next coarsest level.",
                 0, 10);

  db.add<size_t>("multigrid_n_post_smooth", 1,
                 "The number of smoothing steps to apply per level after "
                 "coming up from the next coarsest level.",
                 0, 10);

  db.add("multigrid_coarse_residual", 1.0e-1,
         "The target residual on the coarsest grid. When this residual is "
         "reached on the coarsest grid by solving or smoothing, the coarse "
         "level will return and the process continue on the next finest level.",
         0.0, 1.0);

  db.add<size_t>("multigrid_coarse_max_n_iterations", 10,
         "The maximal number of solver/smoother iterations to be performed "
         "whenever working on the coarsest level.",
         1, 100);

  db.add("multigrid_vanka_damp_factor", 1.0,
         "A damping factor relevant for Vanka type smoothers only. It is "
         "responsible for a damping when adding the solution of the local "
         "defect equation onto the global solution. Vanka smoothers tend to be "
         "quite responsive to this value. Although it defaults to 1.0 (no "
         "damping), a value of 0.8 is often a good start.", 0.0, 1.0);

  db.add("jacobi_damp", 1.0,
         "Damping factor for the Jacobi smoother. If chosen below 1, turns "
         "the Jacobi smoother into a damped Jacobi smoother, which is expected "
         "to perform better.", 0.0, 1.0);


  db.add("sor_omega", 1.5, "The overrelaxation parameter (typically called "
         "omega). This is only used for the (symmetric) successive "
         "overrelaxation method, adn only if that one is used "
         "as a smoother here.", 0., 2.);

#ifdef _MPI
  db.add("sor_parallel_type", "all_cells", "This is in an experimental state. It is "
         "about the parallelization strategy of the SOR (when used as a smoother).",
         {"all_cells", "own_cells"});
#endif
  return db;
}
