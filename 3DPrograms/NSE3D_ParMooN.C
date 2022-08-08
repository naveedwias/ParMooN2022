#include <Domain.h>
#include <Database.h>
#include "NavierStokes.h"
#include <Chrono.h>
#include <LoopInfo.h>
#include "ParMooN.h"

int main(int argc, char* argv[])
{
  auto parmoon_db = parmoon::parmoon_initialize(argc, argv);
  {
#ifdef _MPI
    TDatabase::ParamDB->Comm = MPI_COMM_WORLD;
    int my_rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(my_rank==0)
    {
      Output::print("<<<<< Running ParMooN: NSE3D Main Program >>>>>");
      Output::info("NSE3D", "MPI, using ", size, " processes");
    }
#else
    int my_rank = 0;
    Output::print("<<<<< Running ParMooN: NSE3D Main Program >>>>>");
    Output::info("NSE3D", "SEQUENTIAL (or OMP...)");
#endif

    Chrono timer; //start a stopwatch which measures time spent in program parts

    bool linear_problem = (parmoon_db["problem_type"].is(3)
                         || parmoon_db["problem_type"].is(7));
    
    TDomain domain(parmoon_db);

    if(my_rank==0) //Only one process should do that.
    {
      parmoon_db.write(Output::get_outfile());
      TDatabase::WriteParamDB(argv[0]);
    }

    // Do the parameter check of the Database.
    check_parameters_consistency_NSE(parmoon_db);

    // Intial refinement
    domain.refine_and_get_hierarchy_of_collections(parmoon_db);

    //print information on the mesh partition on the finest grid
    domain.print_info("NSE3D domain");

    timer.restart_and_print("setup(domain, example, database)");
    // Construct an object of the NSE3D-problem type.
    NavierStokes<3> nse3d(domain, parmoon_db);
    timer.restart_and_print("constructing NSE3D object");
    
    // assemble all matrices and right hand side
    nse3d.assemble_linear_terms();
    nse3d.stop_it(0);  // check initial residuals

    LoopInfo<Residuals> loop_info("nonlinear");
    loop_info.print_time_every_step = true;
    loop_info.verbosity_threshold = 1; // full verbosity
    if(my_rank==0)
      loop_info.print(0, nse3d.get_residuals());

    timer.restart_and_print("assembling linear terms");

    //Timer which measures time spent in solve() method solely.
    Chrono timer_sol;
    timer_sol.stop();

    //======================================================================
    for(unsigned int k=1;; k++)
    {
      if(my_rank == 0)
        Output::print<3>(); // new line for a new nonlinear iteration
      // solve the system
      timer_sol.start();
      nse3d.solve();
      timer_sol.stop();

      //no nonlinear iteration for Stokes problem
      if(linear_problem)
        break;

      nse3d.assemble_nonlinear_term();

      // checking residuals
      if(nse3d.stop_it(k))
      {
        loop_info.finish(k, nse3d.get_residuals());
        break;
      }
      else
        loop_info.print(k, nse3d.get_residuals());

    } // end for k
    timer.restart_and_print("nonlinear loop");
    
    timer_sol.print_total_time("solver only");

    nse3d.output();
    timer.restart_and_print("output");
    timer.print_total_time("NSE3D_ParMooN program");

    if(my_rank==0)
      Output::print("<<<<< ParMooN Finished: NSE3D Main Program >>>>>");

    if(my_rank == 0)
      Output::close_file();
  }

  parmoon::parmoon_finalize();
}

