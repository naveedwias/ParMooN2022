#include <Chrono.h>
#include <Domain.h>
#include <Database.h>
#include <SnapshotsCollector.h>
#include "TimeConvectionDiffusion.h"
#include <TCD_POD.h>
#include <TCD_ROM.h>
#include <TimeDiscRout.h>
#include <TimeDiscretizations.h>
#include "ParMooN.h"

int main(int argc, char *argv[])
{
  parmoon::parmoon_initialize(argc, argv);
  Chrono timer;

  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.merge(TimeDiscretization::default_TimeDiscretization_database());
  parmoon_db.read(argv[1]);

  bool i_am_root = true;
#ifdef _MPI
  int my_rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  i_am_root = (my_rank == 0);
#endif

  if(i_am_root)
  {
    Output::print("<<<<< Running ParMooN: Time_CD3D Main Program >>>>>");
  }
#ifdef _MPI
  Output::info("Time_CD3D", "MPI, using ", size, " processes");
#endif

  // ===========================================================================
  // set the database values and generate mesh
  // ===========================================================================
  TDomain domain(parmoon_db);

  // Intial refinement
  domain.refine_and_get_hierarchy_of_collections(parmoon_db);

  //print information on the mesh partition on the finest grid
  domain.print_info("TCD3D domain");

  // set some parameters for time stepping
  SetTimeDiscParameters(0);

  timer.restart_and_print("setup(domain, example, database)");

  // ===========================================================================
  // Solve the PDE using FEM and acquire snapshots
  // ===========================================================================
  bool standard_FEM =  !parmoon_db["write_snaps"]
                    && !parmoon_db["compute_POD_basis"]
                    && !parmoon_db["ROM_method"];

  if( standard_FEM || parmoon_db["write_snaps"] )
  {
    // create an object of the class Time_CD3D
    TimeConvectionDiffusion<3> tcd3d(domain, parmoon_db);
    timer.restart_and_print("constructing Time_CD3D object");

    // initialize snapshot writer (only needed if parmoon_db["write_snaps"])
    // the solutions are stored in a special collector which can be later read
    // by the computation of POD the basis
    SnapshotsCollector snaps(parmoon_db);

    TimeDiscretization& tss = tcd3d.get_time_stepping_scheme();
    tss.current_step_ = 0;
    tss.current_time_ = parmoon_db["time_start"];

    // assemble the matrices and right hand side at the start time
    tcd3d.assemble_initial_time();
    tcd3d.output();
    timer.restart_and_print("initial assembling");
    Chrono timer_solve;
    timer_solve.stop();

    double start_time = parmoon_db["time_start"];
    TDatabase::TimeDB->CURRENTTIME = start_time;


    // store initial condition as snapshot
    if(parmoon_db["write_snaps"])
    {
      snaps.write_data(tcd3d.get_solution().block(0),
                       tcd3d.get_solution().length(0));
    }
    
    if( tcd3d.get_db()["algebraic_flux_correction"].is("fem-fct-cn")
        && tcd3d.get_db()["afc_fct_scheme"].is("non-linear") )
    {
      tss.nonlin_iteration = 0;
    }
    
    while(!tss.reached_final_time_step())
    {
      tss.current_step_++;
      tss.set_time_disc_parameters();
      SetTimeDiscParameters(1);
      tss.current_time_ += tss.get_step_length();
      // this is used at several places, e.g., in the example file etc.
      TDatabase::TimeDB->CURRENTTIME += tss.get_step_length();
      if(i_am_root)
      {
        Output::print<1>("\nCURRENT TIME: ", tss.current_time_);
      }
      
      tcd3d.assemble();
      timer_solve.start();
      tcd3d.solve();
      timer_solve.restart_and_print("solving in step t=" 
                                    + std::to_string(tss.current_time_));
      timer_solve.stop();
      
      tcd3d.output();
      timer.restart_and_print("time step (t="
                              + std::to_string(tss.current_time_)+ ")");

      // write the snapshots
      if(parmoon_db["write_snaps"])
      {
        // solution_m1 and solution_m2 were already updated
        snaps.write_data(tcd3d.get_solution().block(0),
                         tcd3d.get_solution().length(0),
                         tss.current_step_,
                         tcd3d.get_solution_m2().block(0),
                         tss.get_step_length());
      }
    }
    timer_solve.print_total_time("accumulated solver time");
  } // if FEM + snapshots

  // ===========================================================================
  // Compute the POD basis from snapshots
  // ===========================================================================
  if( parmoon_db["compute_POD_basis"] )
  {
    auto collections            = domain.get_grid_collections();
    TCollection& cellCollection = *collections.front();

    TCD_POD<3> tcd_pod(cellCollection, parmoon_db); // TODO add TYPE (i.e. TCD3D, TNSE3D,...)

    tcd_pod.compute_pod_basis();

    tcd_pod.output();
  } // if POD

  // ===========================================================================
  // Solve the PDE using ROM instead of FEM
  // ===========================================================================
  if( parmoon_db["ROM_method"] )
  {
    double start_time = parmoon_db["time_start"];
    TDatabase::TimeDB->CURRENTTIME = start_time;

    TCD_ROM<3> tcd_rom(domain, parmoon_db);

    TimeDiscretization& tss = tcd_rom.get_time_stepping_scheme();
    tss.current_step_ = 0;
    tss.set_time_disc_parameters();

    // assemble matrices and right hand side at start time
    tcd_rom.assemble_initial_time();
    tcd_rom.output();

    // time iteration
    while(!tss.reached_final_time_step())
    {
      tss.current_step_++;
      tss.set_time_disc_parameters();
      SetTimeDiscParameters(1);
      tss.current_time_ += tss.get_step_length();
      TDatabase::TimeDB->CURRENTTIME += tss.get_step_length();
      if(i_am_root)
      {
        Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);
      }

      tcd_rom.assemble();
      tcd_rom.solve();

      tcd_rom.output();
    }
  } // if ROM

  // ===========================================================================
  if(i_am_root)
  {
    Output::print("<<<<< ParMooN Finished: TCD3D Main Program >>>>>");
  }
  timer.print_total_time("TCD3D_ParMooN program");
  // ===========================================================================
  parmoon::parmoon_finalize();
}
