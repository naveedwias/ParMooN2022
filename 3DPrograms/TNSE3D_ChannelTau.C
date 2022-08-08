/**
 * @brief Main program for solving a 3D time-dependent Navier Stokes equation
 *        in a periodic channel using ParMooN.
 *
 */
#include <Domain.h>
#include <Database.h>
#include <Chrono.h>
#include <LoopInfo.h>
#include <TimeDiscretizations.h>
#include <ChannelTauRoutines.h>
#include <TimeDiscRout.h>
#include "ParMooN.h"

// main program
// =============================================================================
int main(int argc, char* argv[])
{
  auto parmoon_db = parmoon::parmoon_initialize(argc, argv);
  {
#ifdef _MPI
  TDatabase::ParamDB->Comm = MPI_COMM_WORLD;
  // Hold mpi rank and size ready, check whether the current processor
  // is responsible for output (usually root, 0).
  int my_rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if(my_rank==0)
  {
    Output::print("<<<<< Running ParMooN: NSE3D Main Program >>>>>");
    Output::info("Time_NSE3D", "MPI, using ", size, " processes");
  }
#else
  int my_rank = 0;
  Output::print("<<<<< Running ParMooN: NSE3D Main Program >>>>>");
  Output::info("Time_NSE3D", "SEQUENTIAL (or OMP...)");
#endif
  // start a stopwatch which measures time spent in program parts
  Chrono timer;

  // Choose and construct example.
  Example_TimeNSE3D example(parmoon_db);

  // Do the parameter check of the Database.
  check_parameters_consistency_NSE(parmoon_db);
  // ===========================================================================
  // set the database values and generate mesh
  // ===========================================================================
  TDomain domain(parmoon_db);

  // open OUTFILE,
  // this is where all output is written to (additionally to console)
  if(parmoon_db["problem_type"].is(0))
  {
    parmoon_db["problem_type"] = 6;
  }

  // Initial refinement
  std::list<TCollection*> gridCollections =
    domain.refine_and_get_hierarchy_of_collections(parmoon_db,
                                                   example.get_bc(0));

  //print information on the mesh partition on the finest grid
  domain.print_info("TNSE3D domain");
  // set some parameters for time stepping
  SetTimeDiscParameters(0);
  // Construct an object of the ChannelFlowRoutines type.
  TNS_ChannelTau tnse3d(domain, parmoon_db);

  TimeDiscretization& tss = tnse3d.get_time_stepping_scheme();
  tss.current_step_ = 0;
  tss.set_time_disc_parameters();
  TDatabase::TimeDB->CURRENTTIME = tss.get_start_time();

  tnse3d.assemble_initial_time();
  tnse3d.output();

  int n_substeps = GetN_SubSteps();

  LoopInfo<double> loop_info_time("time loop");
  loop_info_time.print_time_every_step = true;
  loop_info_time.verbosity_threshold = 1; // full verbosity
  int linear_iterations = 0; 

  timer.restart_and_print("setting up spaces, matrices and initial assembling");
  //============================================================================
  // time iteration
  //============================================================================
  while(!tss.reached_final_time_step())
  {
    // time measuring during every time iteration
    Chrono timer_timeit;

    tss.current_step_++;
    // set the time parameters
    tss.set_time_disc_parameters();
    // loop over substeps in one time iteration
    for(int j = 0 ; j < n_substeps ; ++j)
    {
      // setting the time discretization parameters
      SetTimeDiscParameters(1);
      if(tss.current_step_ == 1 && my_rank==0)
      {
        Output::print<1>("Theta1: ", TDatabase::TimeDB->THETA1);
        Output::print<1>("Theta2: ", TDatabase::TimeDB->THETA2);
        Output::print<1>("Theta3: ", TDatabase::TimeDB->THETA3);
        Output::print<1>("Theta4: ", TDatabase::TimeDB->THETA4);
      }
      // tau may change depending on the time discretization (adaptive time)
      double tau = tss.get_step_length();
      tss.current_time_ += tss.get_step_length();
      // this is used at several places, e.g., in the example file etc.
      TDatabase::TimeDB->CURRENTTIME += tau;

      if(my_rank==0)
      {
        Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);
      }

      tnse3d.assemble_matrices_rhs(0);
      timer_timeit.restart_and_print("preparation of nonlinear iteration");

      // nonlinear iteration
      LoopInfo<Residuals> loop_info("nonlinear");
      loop_info.print_time_every_step = true;
      loop_info.verbosity_threshold = 1; // full verbosity
      for(unsigned int k=0 ; ; k++)
      {
        if(my_rank==0) // some output
        {
          Output::print<1>("\nNONLINEAR ITERATION :", setw(3), k);
        }
        // checking residuals and stop conditions
        if(tnse3d.stop_it(k))
        {
          loop_info.finish(k, tnse3d.get_residuals());
          linear_iterations+=k;
          break;
        }
        else
        {
          loop_info.print(k, tnse3d.get_residuals());
        }
 
        tnse3d.solve();

        tnse3d.assemble_matrices_rhs(k+1);

        timer_timeit.restart_and_print("solving and reassembling in the "
                                       "nonlinear iteration "
                                       + std::to_string(k));
      } // end of nonlinear loop

      loop_info_time.print(linear_iterations, tnse3d.get_full_residual());

      timer_timeit.restart_and_print("solving the time iteration "
                              + std::to_string(TDatabase::TimeDB->CURRENTTIME));

      tnse3d.output();

      timer_timeit.print_total_time("time step "
                              + std::to_string(TDatabase::TimeDB->CURRENTTIME));
    } // end of subtime loop
  } // end of time loop

  loop_info_time.finish(linear_iterations, tnse3d.get_full_residual());

  timer.print_total_time("whole solving procedure ");

  // ===========================================================================
  Output::print("MEMORY: ", setw(11), GetMemory()/(1048576.0), " MB");
  // ===========================================================================
  }
  parmoon::parmoon_finalize();
}
