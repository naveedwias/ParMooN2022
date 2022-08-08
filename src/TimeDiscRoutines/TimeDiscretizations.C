
#include <cmath>
#include <TimeDiscretizations.h>
#include <Database.h>
#include "MooNMD_Io.h"

#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */

/* ****************************************************************************/
ParameterDatabase TimeDiscretization::default_TimeDiscretization_database()
{
  ParameterDatabase db("default ParMooN time parameters database");

  db.add("time_start", 0.,
         "This is the start time of the simulation. The total simulated time "
         "interval is [start_time, end_time].", -1000.,1000.);

  db.add("time_end", 1.,
         "This is the end time of the simulation. The total simulated time "
         "interval is [start_time, end_time].", -1000.,1000.);

  db.add("time_step_length", 0.05,
         "This is the time step length. Without time adaptivity, it "
         "is the same constant value for all time steps, whereas for "
         "time adaptivity, it only corresponds to the initial value.",
         0., 0.5);

  db.add("time_discretization", "backward_euler",
         "Determine the time stepping scheme: currently supported are",
         {"backward_euler", "crank_nicolson", "bdf_two", "fractional_step", "Runge_Kutta_Heun"});

  db.add("bdf_two_first_step", "backward_euler",
         "Determine the time stepping scheme to be used for the first step of "
         "bdf_two: currently supported are",
         {"backward_euler", "crank_nicolson"});

  return db;
}

/* ****************************************************************************/
TimeDiscretization::TimeDiscretization(const ParameterDatabase & param_db)
: db(default_TimeDiscretization_database())
{
  db.merge(param_db);

  if(db["time_discretization"].is("bdf_two") )
  {
    bdf_coefficients = { 4./3., -1./3., 2./3. };
    // pre-stage computation is done either by the
    // backwar-Euler or Crank-Nicolson scheme
    // This parameter is set to False in the main code
    // for the second time step
    pre_stage_bdf = true;
  }
  else
  {
    // rk = std::make_shared<RungeKuttaTable>(db["time_discretization"]);
    pre_stage_bdf = false;
  }

  start_time = db["time_start"];
  current_time_step_length = db["time_step_length"];
  current_time_ = db["time_start"];
  end_time = db["time_end"];

  TDatabase::TimeDB->base_time_step_length = current_time_step_length;
  TDatabase::TimeDB->TIMESTEPLENGTH = current_time_step_length;
  TDatabase::TimeDB->CURRENTTIMESTEPLENGTH = current_time_step_length;

  n_scale_block = 4;

  b_bt_linear_nl = "linear";

  if (db["time_discretization"].is("backward_euler")
    || db["time_discretization"].is("crank_nicolson")
    || db["time_discretization"].is("bdf_two"))
  {
    last_scale_factor = get_scale_factor(current_time_step_length);
  }
  else
  {
    last_scale_factor = 1.0;
  }

  last_b_factor = 1.0;

  Output::print<5>("Time discretization is initialized with ",
                db["time_discretization"]);
}

/* ****************************************************************************/
void TimeDiscretization::prepare_timestep(BlockFEMatrix&,
                                          const BlockFEMatrix&,
                                          std::vector<BlockVector>&,
                                          const std::vector<BlockVector>&)
{
  //TODO: here the idea can be implemented of multisteps
}

/* ****************************************************************************/
void TimeDiscretization::prepare_rhs_from_time_disc(
     BlockFEMatrix& system_matrix, const BlockFEMatrix& mass_matrix,
     std::vector<BlockVector> & rhs,
     const std::vector<BlockVector>& old_solutions)
{
  int n_blocks;
  int n_blocks_mass;

  if (system_matrix.get_n_cell_rows() == 1)
  {
    n_blocks = 1;
    n_blocks_mass = 1;
  }

  if(system_matrix.get_n_cell_rows() == 3)
  {
    n_blocks = 2;
    n_blocks_mass = 3;
  }

  if (system_matrix.get_n_cell_rows() == 4)
  {
    n_blocks = 3;
    n_blocks_mass = 4;
  }

  // at the moment it's implemented in the way it was done in MooNMD
  // but have to be adopted to be fit according to the Butscher table
  if ((db["time_discretization"].is("bdf_two")) && !pre_stage_bdf)
  {
    // right hand side with tau*..
    rhs[0].scaleActive(bdf_coefficients[2] * current_time_step_length);

    mass_matrix.apply_scaled_submatrix(old_solutions[0], rhs[0], n_blocks_mass, n_blocks_mass,
                                           bdf_coefficients[0]);

    mass_matrix.apply_scaled_submatrix(old_solutions[1], rhs[0], n_blocks_mass, n_blocks_mass,
                                           bdf_coefficients[1]);

    Output::print<5>(RED, "bdf2_stage ", pre_stage_bdf, BLACK);
  }
  else if (db["time_discretization"].is("backward_euler")
    || (pre_stage_bdf && db["bdf_two_first_step"].is("backward_euler")))
  {
    rhs[0].scaleActive(current_time_step_length);

    // mass matrix times old solution goes to right hand side
    mass_matrix.apply_scaled_submatrix(old_solutions[0], rhs[0], n_blocks_mass, n_blocks_mass, 1.);

    if ((db["time_discretization"].is("bdf_two") ) && pre_stage_bdf)
    {
      Output::print<5>("First step in BDF2 scheme is performed by the BDF1 scheme");
    }
  }
  else if (db["time_discretization"].is("crank_nicolson")
    || (pre_stage_bdf && db["bdf_two_first_step"].is("crank_nicolson")))
  {
    rhs[0].scaleActive(0.5 * current_time_step_length);
    rhs[0].addScaledActive(rhs[1], 0.5 * current_time_step_length);

    // mass matrix times old solution goes to right hand side
    mass_matrix.apply_scaled_submatrix(old_solutions[0], rhs[0], n_blocks_mass, n_blocks_mass, 1.0);

    // A * u_old
    system_matrix.apply_scaled_submatrix(old_solutions[0], rhs[0], n_blocks, n_blocks,
                                         -0.5 * current_time_step_length);

    if (db["time_discretization"].is("bdf_two") && pre_stage_bdf)
    {
      Output::print<5>("First step in BDF2 scheme is performed by the "
                       "Crank-Nicolson scheme");
    }
  }
  else
  {
    ErrThrow("time discretization scheme ", db["time_discretization"], " is not supported yet");
  }
}

double TimeDiscretization::get_scale_factor(double tau)
{
  if (db["time_discretization"].is("backward_euler")
    || (pre_stage_bdf && db["bdf_two_first_step"].is("backward_euler")))
  {
    return tau;
  }
  else if( (db["time_discretization"].is("bdf_two") ) && !pre_stage_bdf)
  {
    return tau * bdf_coefficients[2];
  }
  else if(db["time_discretization"].is("crank_nicolson")
    || (pre_stage_bdf && db["bdf_two_first_step"].is("crank_nicolson")))
  {
    return tau * 0.5;
  }
  else
  {
    ErrThrow("Time stepping scheme ", db["time_discretization"], " is not supported");
  }
}

/* ****************************************************************************/
void TimeDiscretization::prepare_system_matrix(
  BlockFEMatrix& system_matrix, const BlockFEMatrix& mass_matrix)
{
  std::vector<std::vector<size_t>> cells;

  if (system_matrix.get_n_cell_rows() == 3)
  {
    cells = {{0,0}, {0,1},
             {1,0}, {1,1} };
  }
  else if (system_matrix.get_n_cell_rows() == 4)
  {
    cells = {{0,0}, {0,1}, {0,2},
             {1,0}, {1,1}, {1,2},
             {2,0}, {2,1}, {2,2} };
  }
  else if (system_matrix.get_n_cell_rows() == 1)
  {
    cells = {{0,0}};
  }

  double factor = get_scale_factor(current_time_step_length);

  // checl this also for comparison add_matrix_actives
  system_matrix.scale_blocks_actives(factor, cells);

  // add the scaled block matrices
  system_matrix.add_blockfe_matrix(mass_matrix);

  Output::print<5>("addition of the mass and stiffness matrix done ");

  last_scale_factor = factor;
}

/* ****************************************************************************/
void TimeDiscretization::reset_linear_matrices(BlockFEMatrix& matrix,
                                               const BlockFEMatrix& mass)
{
  // first subtract the mass matrix from the matrix
  matrix.add_blockfe_matrix(mass, -1.);

  // reset the matrices
  std::vector<std::vector<size_t>> cells;
  if (matrix.get_n_cell_rows() == 3)
  {
    cells = {{0,0}, {0,1},
             {1,0}, {1,1} };
  }
  else if (matrix.get_n_cell_rows() == 4)
  {
    cells = {{0,0}, {0,1}, {0,2},
             {1,0}, {1,1}, {1,2},
             {2,0}, {2,1}, {2,2} };
  }
  else if (matrix.get_n_cell_rows() == 1)
  {
    cells = {{0,0}};
  }

  matrix.scale_blocks_actives(1.0 / last_scale_factor, cells);
  // No need to reset the B, BT and C blocks because they are only scaled once during
  // the time step or re-assembled and scaled in the nonlinear case.

  last_scale_factor = 1.0;
}

/* ****************************************************************************/
void TimeDiscretization::scale_descale_all_b_blocks(BlockFEMatrix& matrix)
{
  /// @DETAILS: scaling of the B and BT's blocks (also "C")
  double factor;

  if(db["time_discretization"].is("backward_euler")
    || (pre_stage_bdf && db["bdf_two_first_step"].is("backward_euler")))
  {
    factor = current_time_step_length;
  }
  else if((db["time_discretization"].is("bdf_two")) && !pre_stage_bdf)
  {
    // if time step size changes, then multiply with 
    // bdf_coefficients[2]*current_time_step_length/old_time_step_length
    // here we assume that the time step size is constant. The blocks have been
    // multiplied by the current_time_step_length in the pre_stage_bdf already
    factor = bdf_coefficients[2];
  }
  else if(db["time_discretization"].is("crank_nicolson")
    || (pre_stage_bdf && db["bdf_two_first_step"].is("crank_nicolson")))
  {
    factor = current_time_step_length;
  }
  else
  {
    ErrThrow("Time stepping scheme ", db["time_discretization"], " is not supported");
  }

  std::vector<std::vector<size_t>> cells;
  if (n_scale_block == 5 || n_scale_block == 7)
  {
#ifdef __2D__
    cells = {{0,2}, {1,2},
             {2,0}, {2,1}, {2,2}};
#endif
#ifdef __3D__
    cells = {{0,3}, {1,3}, {2,3},
             {3,0}, {3,1}, {3,2}, {3,3}};
#endif
  }
  else
  {
#ifdef __2D__
    cells = {{0,2}, {1,2}};
    if (TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT > 0)
    {
      cells = {{0,2}, {1,2},
               {2,0}, {2,1}};
      factor *= TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT;
    }
#endif
#ifdef __3D__
    cells = {{0,3}, {1,3}, {2,3}};
    if (TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT > 0)
    {
      cells = {{0,3}, {1,3}, {2,3},
               {3,0}, {3,1}, {3,2}};
      factor *= TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT;
    }
#endif
  }

  // in the case of bdf2:
  // the first step is done by backward Euler, 
  // then at the second time step rescale with the BDF2 coefficients  
  if (current_step_ == 1 ||
    (db["time_discretization"].is("bdf_two") && current_step_ == 2) )
  {
    matrix.scale_blocks(factor, cells);
  }
  else if (factor != last_b_factor)
  {
    matrix.scale_blocks(factor / last_b_factor, cells);
  }

  last_b_factor = factor;
}

/* ****************************************************************************/
void TimeDiscretization::scale_nl_b_blocks(BlockFEMatrix& matrix)
{
  /// @DETAILS: scaling of the B and BT's blocks (also "C")
  double factor;

  if (db["time_discretization"].is("backward_euler")
    || (pre_stage_bdf && db["bdf_two_first_step"].is("backward_euler")))
  {
    factor = current_time_step_length;
  }
  else if ((db["time_discretization"].is("bdf_two") )&& !pre_stage_bdf)
  {
    factor = current_time_step_length * bdf_coefficients[2];
  }
  else if (db["time_discretization"].is("crank_nicolson")
    || (pre_stage_bdf && db["bdf_two_first_step"].is("crank_nicolson")))
  {
    factor = current_time_step_length;//this have to be adopted according to the stages
  }
  else
  {
    ErrThrow("Time stepping scheme ", db["time_discretization"], " is not supported");
  }

  std::vector<std::vector<size_t>> cells;// = {{0,2},{1,2}};
  if (b_bt_linear_nl.compare("nonlinear") == 0)
  {
#ifdef __2D__
    if (n_scale_block == 2)
    {
      cells = {{0,2}, {1,2}};
    }
    else if (n_scale_block == 5)
    {
      cells = {{0,2}, {1,2}, {2,0}, {2,1}, {2,2}};
    }
    else
    {
      ErrThrow("No of blocks are out of range");
    }
#endif
#ifdef __3D__
    if (n_scale_block == 3)
    {
      cells = {{0,3}, {1,3}, {2,3}};
    }
    else if (n_scale_block == 7)
    {
      cells = {{0,3}, {1,3}, {2,3}, {3,0}, {3,1}, {3,2}, {3,3}};
    }
    else
    {
      ErrThrow("No of blocks are out of range");
    }
#endif
  }

  matrix.scale_blocks(factor, cells);
}

/* ****************************************************************************/
void TimeDiscretization::set_time_disc_parameters()
{
  // the first step needs to be done by the backward Euler
  // or Crank-Nicolson scheme
  if(db["time_discretization"].is("backward_euler"))
  {
    this->pre_stage_bdf = false;
    TDatabase::TimeDB->THETA1 = 1.0;
    TDatabase::TimeDB->THETA2 = 0.0;
    TDatabase::TimeDB->THETA3 = 0.0;
    TDatabase::TimeDB->THETA4 = 1.0;
  }
  if((db["time_discretization"].is("bdf_two"))&& (current_step_<=1) )
  {
    if(db["bdf_two_first_step"].is("backward_euler"))
    {
      TDatabase::TimeDB->THETA1 = 1.0;
      TDatabase::TimeDB->THETA2 = 0.0;
      TDatabase::TimeDB->THETA3 = 0.0;
      TDatabase::TimeDB->THETA4 = 1.0;
    }
    else if(db["bdf_two_first_step"].is("crank_nicolson"))
    {
      TDatabase::TimeDB->THETA1 = 0.5;
      TDatabase::TimeDB->THETA2 = 0.5;
      TDatabase::TimeDB->THETA3 = 0.5;
      TDatabase::TimeDB->THETA4 = 0.5;
    }
  }
  if((db["time_discretization"].is("bdf_two") )&&
          (this->current_step_ >= 2))
  {
    this->pre_stage_bdf = false;
    TDatabase::TimeDB->THETA1 = 2./3;
    TDatabase::TimeDB->THETA2 = 0.0;
    TDatabase::TimeDB->THETA3 = 0.0;
    TDatabase::TimeDB->THETA4 = 2./3.;
    if(current_step_==2)
    {
      Output::print("BDF2 with new parameters : ", bdf_coefficients[0], "  ",
                    bdf_coefficients[1],"  ", bdf_coefficients[2]);
    }
  }
  if(db["time_discretization"].is("crank_nicolson"))
  {
    this->pre_stage_bdf = false;
    TDatabase::TimeDB->THETA1 = 0.5;
    TDatabase::TimeDB->THETA2 = 0.5;
    TDatabase::TimeDB->THETA3 = 0.5;
    TDatabase::TimeDB->THETA4 = 0.5;
  }
  // set the global parameters used in the local assembling routines 
  // and at some other places
  TDatabase::TimeDB->TIMESTEPLENGTH = current_time_step_length;
  TDatabase::TimeDB->CURRENTTIMESTEPLENGTH = current_time_step_length;
  Output::print<5>("BDF2 scheme with the steps perform ", pre_stage_bdf);
}

void TimeDiscretization::set_step_length(double tau)
{
  if (tau != current_time_step_length)
  {
    Output::root_warn("Time stepping scheme", "Changed time step! Make sure you're using a robust scheme.");
  }

  current_time_step_length = tau;
  TDatabase::TimeDB->TIMESTEPLENGTH = current_time_step_length;
  TDatabase::TimeDB->CURRENTTIMESTEPLENGTH = current_time_step_length;
}