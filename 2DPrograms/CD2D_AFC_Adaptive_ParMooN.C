// =======================================================================
//
// Purpose:  main program for solving a stationary scalar equation using ParMooN
//
// Author:   Sashikumaar Ganesan
//
// History:  Implementation started on 22.08.2014
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include "ConvectionDiffusion_AFC.h"
#include "CDErrorEstimator_AFC.h"
#include "CDErrorEstimator.h"
#include "RefinementStrategy.h"
#include "LoopInfo.h"
#include "Chrono.h"
#include "LocalAssembling.h"
#include "Multigrid.h"
#include "ParMooN.h"

#include <sys/stat.h>
#include <sys/types.h>


// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  parmoon::parmoon_initialize(argc, argv);
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.merge(ErrorEstimator<2>::default_error_estimator_database());
  parmoon_db.read(argv[1]);
  
  TDomain domain(parmoon_db);
  // refine grid
  domain.refine_and_get_hierarchy_of_collections(parmoon_db);
   
  std::string base_name = parmoon_db["output_basename"];
  
  size_t afc_estimator_type = 1;
  if(parmoon_db.contains("afc_estimator_type"))
    afc_estimator_type = parmoon_db["afc_estimator_type"];
  
  //=========================================================================
  CDErrorEstimator<2> estimator_initial(parmoon_db);
  
  parmoon_db["estimator_type"] = afc_estimator_type;
  
  CDErrorEstimator_AFC<2> estimator(parmoon_db);
  RefinementStrategy<2> refinementStrategy(parmoon_db);
  LoopInfo<double> loop_info("adaptive", true, true, 1);
  BlockVector values;
  int n_adaptive_steps = (int)domain.get_database()["refinement_max_n_adaptive_steps"]
                           + (int)domain.get_database()["refinement_max_n_uniform_steps"];
  int n_uniform_steps = 0;
  bool adaptive_converged = false;
 
  /* The idea is to use the last solution of the nth adaptive grid as the initial
   * solution of the (n+1)th grid. This works for grids with hanging nodes but
   * not for conforming grids. Hence, this code is highlighted */
  //int smear_done = 0;
  //TFEFunction2D function_coarse;
  //std::vector<double> interpolated_solution;
  for(int curr_level = 0; ; ++curr_level)
  {
    Output::print("\nadaptive loop ", curr_level);
    std::ostringstream ostr;
    ostr << base_name << "_" << std::setw(2) << std::setfill('0') << curr_level;
    parmoon_db["output_basename"].set(ostr.str(), false);
    
    ConvectionDiffusion_AFC<2> cd2d(domain, parmoon_db);
    //Only works for grids with hanging nodes
    /*if(curr_level > 0 && cd2d.get_space()->get_n_hanging() != 0)
    {
      interpolated_solution.resize(cd2d.get_space()->get_n_dof(), 0.0);
      TFEFunction2D function_interpolate(cd2d.get_space(), "interpolant", interpolated_solution.data(),
                                         cd2d.get_space()->get_n_dof());
      function_interpolate.Interpolate(&function_coarse);
    }*/
     // assemble and solve poisson equation with right-hand side
    cd2d.assemble_poisson();
    cd2d.solve_poisson();
    cd2d.assemble(0);
    cd2d.solve(0);
    if( cd2d.get_db()["algebraic_flux_correction"].is("afc") )
    {//nonlinear loop necessary
      size_t Max_It = cd2d.get_db()["afc_nonlinloop_maxit"];
      for(unsigned int k = 1;; k++)
      {
        bool converged;
        converged = cd2d.solve(k);
        
        if ((converged)||(k>= Max_It))
          break;
      }
    }

    // restrict solution to coarser grids
    // TO BE DONE
    //function_coarse = cd2d.get_function();
    
    estimator.estimate(cd2d.get_example(), cd2d.get_function(), 
                       cd2d.get_afc_D_entries(), cd2d.get_afc_alphas(),
                       cd2d.get_initial_fe_function());
    if(afc_estimator_type)
    {
      estimator_initial.estimate(cd2d.get_example(), 
                               cd2d.get_initial_fe_function());
      estimator.add_eta_K(estimator_initial);    
    }
    estimator.info();
    
    
    ConvectionDiffusion_AFC<2>::FEFunction estimated_error(estimator, values);
    cd2d.add_to_output(&estimated_error);
    cd2d.output();
    
    auto fespace = cd2d.get_space();

    //Computation of cutline only when DOFs > 40000. This is specific for Hemker problem
    // NOTE: Will be removed later
    /*if(cd2d.get_space()->get_n_dof()>40000 && smear_done == 0)
    {
      cd2d.ComputeCutLineY();
      //cd2d.ComputeCutLineX();
      smear_done = 1;
    }*/
    
    //Putting three conditions for checking the end of adaptive algorithm
    // TODO: Put the DOF stopping criteria in the database
    if(estimator.get_estimated_global_error()<1.0e-3 || cd2d.get_space()->get_n_dof()>100000 ||
       curr_level == n_adaptive_steps)
    {
      adaptive_converged = true;
    }
    if(adaptive_converged)
    {
      loop_info.finish(curr_level, estimator.get_estimated_global_error());
      break;
    }
    else
    {
      loop_info.print(curr_level,
                      estimator.get_estimated_global_error());
      if(n_uniform_steps != (int)domain.get_database()["refinement_max_n_uniform_steps"])
      {
        domain.RegRefineAll();
        n_uniform_steps++;
      }
      else  
      {
        refinementStrategy.apply_estimator(estimator);
        domain.RefineByRefinementStrategy(refinementStrategy);
      }
      if(domain.get_database()["write_adaptive_mesh"])
      {
       std::string name_after = std::string("domain_after")
                                + std::to_string(curr_level) + ".ps";
       domain.PS(name_after.c_str(), It_Finest, 0);
      }
      // interpolate solution from next coarser grid
      // TO BE DONE
    }
  }
  //=========================================================================
  parmoon::parmoon_finalize();
} // end main
