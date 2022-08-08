#include <Domain.h>
#include <Database.h>
#include "NavierStokes.h"
#include <Example_NSE2D.h>
#include <Chrono.h>
#include <LoopInfo.h>
#include <ParameterDatabase.h>
#include "ParMooN.h"

int main(int argc, char* argv[])
{
  Chrono timer; // start a stopwatch which measures time spent in program parts
  
  parmoon::parmoon_initialize(argc, argv);
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.merge(NavierStokes<2>::default_nse_database());
  parmoon_db.read(argv[1]);
  
  bool linear_problem = (parmoon_db["problem_type"].is(3)
                         || parmoon_db["problem_type"].is(7));
  TDatabase::ParamDB->INTERNAL_PROBLEM_LINEAR = linear_problem;
  if(linear_problem)
    TDatabase::ParamDB->FLOW_PROBLEM_TYPE = OSEEN;
  
  TDomain domain(parmoon_db);
  
  // possibly change parameters in the database, if they are not meaningful now
  check_parameters_consistency_NSE(parmoon_db);
  
  // Choose and construct example.
  Example_NSE2D example(parmoon_db);

  // refine grid
  domain.refine_and_get_hierarchy_of_collections(parmoon_db, example.get_bc(0));
  
  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    domain.PS("Domain.ps", It_Finest, 0);
  
  // create an object of the Navier-Stokes class
  NavierStokes<2> ns(domain, parmoon_db);
  ns.assemble_linear_terms();
  // if solution was not zero up to here, you should call 
  //ns.assemble_nonlinear_term();
  
  // Routines for periodic boundary conditions (for example Brinkman_Riverbed.h)
  if(parmoon_db["example"].is(15))
  {
    ns.findPeriodicDOFs(0.0, 2.0);
    ns.checkPeriodicDOFs();
    ns.makePeriodicBoundary();
  }

  ns.stop_it(0);
  
  LoopInfo<Residuals> loop_info("nonlinear");
  loop_info.print_time_every_step = true;
  loop_info.verbosity_threshold = 1; // full verbosity
  loop_info.print(0, ns.get_residuals());
  
  timer.restart_and_print("setting up spaces, matrices, linear assemble");
  
  //======================================================================
  // nonlinear loop
  // in function 'stopIt' termination condition is checked
  for(unsigned int k = 1;; k++)
  {
    Output::print<3>(); // new line for a new nonlinear iteration
    ns.solve();
    
    //no nonlinear iteration for Stokes or Brinkman problems
    if(linear_problem)
      break;
    
    ns.assemble_nonlinear_term();
    
    if(ns.stop_it(k))
    {
      loop_info.finish(k, ns.get_residuals());
      break;
    }
    else
      loop_info.print(k, ns.get_residuals());
  } // end for k
  
  timer.restart_and_print("solving procedure");
  
  ns.output();
  
  timer.print_total_time("complete NSE2D_ParMooN program");
  
  parmoon::parmoon_finalize();
}
