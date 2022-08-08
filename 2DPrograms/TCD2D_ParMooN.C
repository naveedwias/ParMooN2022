#include <Database.h>
#include <Domain.h>
#include "MainUtilities.h"
#include <SnapshotsCollector.h>
#include <TimeConvectionDiffusion.h>
#include <TCD_POD.h>
#include <TCD_ROM.h>
#include <TimeDiscRout.h>
#include "ParMooN.h"

int main(int argc, char* argv[])
{
  auto parmoon_db = parmoon::parmoon_initialize(argc, argv);
  // ===========================================================================
  // set the database values and generate mesh
  // ===========================================================================
  TDomain Domain(parmoon_db);

  // refine grid
  Domain.refine_and_get_hierarchy_of_collections(parmoon_db);

  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
  {
    Domain.PS("Domain.ps", It_Finest, 0);
  }

  // ===========================================================================
  // Solve the PDE using FEM and acquire snapshots
  // ===========================================================================
  bool standard_FEM =  !parmoon_db["write_snaps"]
                    && !parmoon_db["compute_POD_basis"]
                    && !parmoon_db["ROM_method"];

  if( standard_FEM || parmoon_db["write_snaps"] )
  {
    TimeConvectionDiffusion<2> tcd(Domain, parmoon_db);

    // initialize snapshot writer (only needed if parmoon_db["write_snaps"])
    // the solutions are stored in a special collector which can be later read
    // by the computation of POD the basis
    SnapshotsCollector snaps(parmoon_db);

    TimeDiscretization& tss = tcd.get_time_stepping_scheme();
    tss.current_step_ = 0;
    tss.set_time_disc_parameters();

    // assemble matrices and right hand side at start time
    tcd.assemble_initial_time();

    double start_time = parmoon_db["time_start"];
    TDatabase::TimeDB->CURRENTTIME = start_time;
    tcd.output();

    // store initial condition as snapshot
    if(parmoon_db["write_snaps"])
    {
      snaps.write_data(tcd.get_solution().block(0),
                       tcd.get_solution().length(0));
    }
    
    if( tcd.get_db()["algebraic_flux_correction"].is("fem-fct-cn")
        && tcd.get_db()["afc_fct_scheme"].is("non-linear") )
    {
      tss.nonlin_iteration = 0;
    }
    
    // time iteration
    while(!tss.reached_final_time_step())
    {
      tss.current_step_++;
      // Output::print("mem before: ", GetMemory());
      tss.set_time_disc_parameters();
      tss.current_time_ += tss.get_step_length();
      double tau = parmoon_db["time_step_length"];

      TDatabase::TimeDB->CURRENTTIME += tau;
      Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);
      SetTimeDiscParameters(1);

      tcd.assemble();
      tcd.solve();
      
      tcd.output();

      // write the snapshots
      if(parmoon_db["write_snaps"])
      {
        // solution_m1 and solution_m2 were already updated
        snaps.write_data(tcd.get_solution().block(0),
                         tcd.get_solution().length(0), tss.current_step_,
                         tcd.get_solution_m2().block(0), tau);
      }
    }
  } // if FEM + snapshots

  // ===========================================================================
  // Compute the POD basis from snapshots
  // ===========================================================================
  if( parmoon_db["compute_POD_basis"] )
  {
    auto collections            = Domain.get_grid_collections();
    TCollection& cellCollection = *collections.front();

    TCD_POD<2> tcd_pod(cellCollection, parmoon_db); // TODO add TYPE (i.e. TCD2D, TNSE2D,...)

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

    TCD_ROM<2> tcd_rom(Domain, parmoon_db);

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
      // Output::print("mem before: ", GetMemory());
      tss.set_time_disc_parameters();
      tss.current_time_ += tss.get_step_length();
      double tau = parmoon_db["time_step_length"];

      TDatabase::TimeDB->CURRENTTIME += tau;
      Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);
      SetTimeDiscParameters(1);

      tcd_rom.assemble();
      tcd_rom.solve();

      tcd_rom.output();
    }
  } // if ROM

  // ===========================================================================
  Output::print("MEMORY: ", setw(10), GetMemory()/(1048576.0), " MB");
  // ===========================================================================
  parmoon::parmoon_finalize();
} // end main
