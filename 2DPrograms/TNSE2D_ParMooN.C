#include <Domain.h>
#include <Database.h>
#include <Example_TimeNSE2D.h>
#include <SnapshotsCollector.h>
#include "TimeNavierStokes.h"
#include <TimeDiscretizations.h>
#include <TimeDiscRout.h>
#include <TNSE_POD.h>
#include <TNSE_ROM.h>
#include <LoopInfo.h>
#include "ParMooN.h"

using namespace std;

int main(int argc, char* argv[])
{
  auto parmoon_db = parmoon::parmoon_initialize(argc, argv);
  unsigned int dim = 2;

  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================
  TDomain Domain(parmoon_db);

  parmoon_db.write(Output::get_outfile());
  TDatabase::WriteParamDB(argv[0]);
  TDatabase::WriteTimeDB();

  // refine grid
  Domain.refine_and_get_hierarchy_of_collections(parmoon_db);

  // write grid into an Postscript file
  if (parmoon_db["output_write_ps"])
  {
    Domain.PS("Domain.ps", It_Finest, 0);
  }

  // set some parameters for time stepping
  SetTimeDiscParameters(0);
  // ===========================================================================
  // Solve the PDE using FEM and acquire snapshots
  // ===========================================================================
  bool standard_FEM =  !parmoon_db["write_snaps"]
                    && !parmoon_db["compute_POD_basis"]
                    && !parmoon_db["ROM_method"];

  if (standard_FEM || parmoon_db["write_snaps"])
  {
    // create an object of TimeNavierStokes<2> class
    TimeNavierStokes<2> tnse2d(Domain, parmoon_db);

    // initialize snapshot writer (only needed if parmoon_db["write_snaps"])
    // the solutions are stored in a special collector which can be later read
    // by the computation of POD the basis
    SnapshotsCollector snaps_u(parmoon_db, "u");
    SnapshotsCollector snaps_p(parmoon_db, "p");

    TimeDiscretization& tss = tnse2d.get_time_stepping_scheme();
    tss.current_step_ = 0;
    tss.set_time_disc_parameters();

    // assemble everything at the start time
    // this includes assembling of all A's, B's
    // and M's blocks that are necessary
    tnse2d.assemble_initial_time();

    // store initial condition
    tnse2d.output();
    if (parmoon_db["write_snaps"])
    {
      snaps_u.write_data(tnse2d.get_solution().block(0),
                         tnse2d.get_solution().length(0) * dim);
      snaps_p.write_data(tnse2d.get_solution().block(dim),
                         tnse2d.get_solution().length(dim));
    }

    LoopInfo<double> loop_info_time("time loop");
    loop_info_time.print_time_every_step = true;
    loop_info_time.verbosity_threshold = 1;
    int linear_iteration = 0;

    TDatabase::TimeDB->CURRENTTIME = tss.get_start_time();
    while (!tss.reached_final_time_step())
    {
      tss.current_step_++;

      // set the time parameters
      tss.set_time_disc_parameters();
      // tau may change depending on the time discretization (adaptive time)
      double tau = tss.get_step_length();
      tss.current_time_ += tss.get_step_length();
      // this is used at several places, e.g., in the example file etc.
      TDatabase::TimeDB->CURRENTTIME += tau;
      Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);

      tnse2d.assemble_matrices_rhs(0);

      LoopInfo<Residuals> loop_info("nonlinear");
      loop_info.print_time_every_step = true;
      loop_info.verbosity_threshold = 1;
      for (unsigned int i = 0; ; i++)
      {
        if (tnse2d.stop_it(i))
        {
          loop_info.finish(i,tnse2d.get_residuals());
          linear_iteration += i;
          loop_info_time.print(linear_iteration, tnse2d.get_full_residual());
          break;
        }
        else
        {
          loop_info.print(i, tnse2d.get_residuals());
        }
        tnse2d.solve();
        tnse2d.assemble_matrices_rhs(i + 1);
      }
      tnse2d.output();

      if (parmoon_db["write_snaps"])
      {
        if (parmoon_db["snaps_time_derivative_projection"])
        {
          BlockVector dt_form(tnse2d.get_solution());

          tnse2d.get_time_derivative(dt_form);

          snaps_u.write_data(tnse2d.get_solution().block(0),
                             tnse2d.get_solution().length(0) * dim,
                             tss.current_step_,
                             dt_form.block(0));

          snaps_p.write_data(tnse2d.get_solution().block(dim),
                             tnse2d.get_solution().length(dim),
                             tss.current_step_,
                             dt_form.block(dim));
        }
        else
        {
          // solution_m1 and solution_m2 were already updated
          snaps_u.write_data(tnse2d.get_solution().block(0),
                             tnse2d.get_solution().length(0) * dim,
                             tss.current_step_,
                             tnse2d.get_solution_m2().block(0),
                             tau);

          snaps_p.write_data(tnse2d.get_solution().block(dim),
                             tnse2d.get_solution().length(dim),
                             tss.current_step_,
                             tnse2d.get_solution_m2().block(dim),
                             tau);
        }
      }
    }
    loop_info_time.finish(linear_iteration, tnse2d.get_full_residual());
  } // if FEM + snapshots

  // ===========================================================================
  // Compute the POD basis from snapshots
  // ===========================================================================
  if (parmoon_db["compute_POD_basis"])
  {
    auto collections            = Domain.get_grid_collections();
    TCollection& cellCollection = *collections.front();

    if (!parmoon_db["write_snaps"])
    {
      // create an object of TimeNavierStokes<2> class to modify
      // TDatabase::ParamDB->VELOCITY_SPACE and
      // TDatabase::ParamDB->PRESSURE_SPACE
      // according to get_velocity_pressure_orders()
      TimeNavierStokes<2> tnse2d(Domain, parmoon_db);
    }
    TNSE_POD<2> tnse_pod(cellCollection, parmoon_db);

    tnse_pod.compute_pod_basis(parmoon_db);
  } // if POD

  // ===========================================================================
  // Solve the PDE using ROM instead of FEM
  // ===========================================================================
  if (parmoon_db["ROM_method"])
  {
    TNSE_ROM<2> tnse_rom(Domain, parmoon_db);

    TimeDiscretization& tss = tnse_rom.get_time_stepping_scheme();
    tss.current_step_ = 0;
    tss.set_time_disc_parameters();

    // assemble everything at the start time
    // this includes all the offline computation
    tnse_rom.assemble_initial_time();

    // store initial condition
    tnse_rom.output();

    int linear_iteration=0;
    TDatabase::TimeDB->CURRENTTIME = tss.get_start_time();
    while(!tss.reached_final_time_step())
    {
      tss.current_step_++;

      // set the time parameters
      tss.set_time_disc_parameters();
      // tau may change depending on the time discretization (adaptive time)
      double tau = tss.get_step_length();
      tss.current_time_ += tss.get_step_length();
      // this is used at several places, e.g., in the example file etc.
      TDatabase::TimeDB->CURRENTTIME += tau;
      Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);

      bool is_velocity = true;
      LoopInfo<double> loop_info("velocity (nl)");
      loop_info.print_time_every_step = true;
      loop_info.verbosity_threshold = 1;

      for (unsigned int i = 0;; i++)
      {
        tnse_rom.assemble_matrices_rhs(i);

        if (tnse_rom.stop_it(i))
        {
          loop_info.finish(i, tnse_rom.get_velocity_residual()); //ROM velocity
          linear_iteration++;
          tnse_rom.print_fem_residuals(linear_iteration);
          break;
        }
        else
        {
          loop_info.print(i, tnse_rom.get_velocity_residual()); //ROM velocity
        }

        tnse_rom.solve(is_velocity);
      }

      tnse_rom.output();
    }
  } // if ROM

  parmoon::parmoon_finalize();
}
