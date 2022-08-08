#include <Domain.h>
#include <Database.h>
#include "NavierStokes.h"
#include <Example_NSE2D.h>
#include <Chrono.h>
#include <LoopInfo.h>
#include <ParameterDatabase.h>
#include "AuxParam2D.h"
#include "ParMooN.h"

#include <RefinementStrategy.h>
#include "NSEErrorEstimator.h"

// =======================================================================
// main program
// =======================================================================
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

  /** set variables' value in TDatabase using argv[1] (*.dat file) */
  TDomain domain(parmoon_db);

  // possibly change parameters in the database, if they are not meaningful now
  check_parameters_consistency_NSE(parmoon_db);
  
  // refine grid
  domain.refine_and_get_hierarchy_of_collections(parmoon_db);

  std::string base_name = parmoon_db["output_basename"];

  //=========================================================================
  NSEErrorEstimator<2> estimator {parmoon_db};
  RefinementStrategy<2> refinementStrategy(parmoon_db);
  LoopInfo<double> adaptive_loop_info("adaptive", true, true, 1);
  BlockVector values;
  size_t n_adaptive_steps = domain.get_database()["refinement_max_n_adaptive_steps"];
  TAuxParam2D aux;

  for(size_t curr_level = 0; ; ++curr_level)
  {
    Output::print("\nadaptive loop ", curr_level);
    std::ostringstream ostr;
    ostr << base_name << "_" << std::setw(2) << std::setfill('0') << curr_level;
    parmoon_db["output_basename"].set(ostr.str(), false);

    //=========================================================================
    // solving (Navier-)Stokes on this level
    // create an object of the Navier-Stokes class
    NavierStokes<2> ns(domain, parmoon_db);
    ns.assemble_linear_terms();
    // if solution was not zero up to here, you should call
    //ns.assemble_nonlinear_term();

    ns.stop_it(0);
    LoopInfo<Residuals> nonlinear_loop_info("nonlinear");
    nonlinear_loop_info.print_time_every_step = true;
    nonlinear_loop_info.verbosity_threshold = 1; // full verbosity
    nonlinear_loop_info.print(0, ns.get_residuals());

    timer.restart_and_print("setting up spaces, matrices, linear assemble");

    //======================================================================
    // nonlinear loop
    // in function 'stopIt' termination condition is checked
    for(unsigned int k = 1;; k++)
    {
      ns.solve();

      //no nonlinear iteration for Stokes or Brinkman problems
      if(linear_problem)
        break;

      ns.assemble_nonlinear_term();

      if(ns.stop_it(k))
      {
        nonlinear_loop_info.finish(k, ns.get_residuals());
        break;
      }
      else
        nonlinear_loop_info.print(k, ns.get_residuals());
    } // end for k

    timer.restart_and_print("solving procedure");

    //=========================================================================
    // done solving nonlinear problem, next: estimate errors, refine grid

    estimator.estimate(ns.get_example(), ns.get_velocity(), ns.get_pressure(),
                       aux);
    estimator.info();

    NavierStokes<2>::FEFunction estimated_error(estimator, values);
    ns.add_to_output(&estimated_error);
    ns.output();

    if(estimator.get_estimated_global_error()<1.0e-10 || curr_level == n_adaptive_steps)
    {
      adaptive_loop_info.finish(curr_level,
                                estimator.get_estimated_global_error());
      break;
    }
    else
    {
      adaptive_loop_info.print(curr_level,
                               estimator.get_estimated_global_error());
      refinementStrategy.apply_estimator(estimator);
      domain.RefineByRefinementStrategy(refinementStrategy);
      if(domain.get_database()["write_adaptive_mesh"])
      {
        std::string name_after = std::string("domain_after")
                               + std::to_string(curr_level) + ".ps";
        domain.PS(name_after.c_str(), It_Finest, 0);
      }
    }
  }
  //=========================================================================
  parmoon::parmoon_finalize();
} // end main
