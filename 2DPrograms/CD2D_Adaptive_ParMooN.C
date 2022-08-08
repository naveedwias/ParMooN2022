#include <Domain.h>
#include <Database.h>
#include "ConvectionDiffusion.h"
#include "CDErrorEstimator.h"
#include "RefinementStrategy.h"
#include "LoopInfo.h"
#include "Chrono.h"
#include "LocalAssembling.h"
#include "Multigrid.h"

#include <sys/stat.h>
#include <sys/types.h>


// =======================================================================
// main program
// =======================================================================
int main(int, char* argv[])
{
  //  declaration of database, you need this in every program
  TDatabase::create(argv[1]);
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.merge(ErrorEstimator<2>::default_error_estimator_database());
  parmoon_db.read(argv[1]);
  
  //open OUTFILE, this is where all output is written to (additionally to console)
  Output::set_outfile(parmoon_db["outfile"], parmoon_db["script_mode"]);
  Output::setVerbosity(parmoon_db["verbosity"]);
  
  TDomain domain(parmoon_db);
  
  // write all Parameters to the OUTFILE (not to console) for later reference
  parmoon_db.write(Output::get_outfile());
  TDatabase::WriteParamDB(argv[0]);
  
  // refine grid
  domain.refine_and_get_hierarchy_of_collections(parmoon_db);
   
  std::string base_name = parmoon_db["output_basename"];
  
  //=========================================================================
  CDErrorEstimator<2> estimator(parmoon_db);
  
  RefinementStrategy<2> refinementStrategy(parmoon_db);
  LoopInfo<double> loop_info("adaptive", true, true, 1);
  BlockVector values;
   size_t n_adaptive_steps = domain.get_database()["refinement_max_n_adaptive_steps"];
  bool adaptive_converged = false;
  for(size_t curr_level = 0; ; ++curr_level)
  {
    Output::print("\nadaptive loop ", curr_level);
    std::ostringstream ostr;
    ostr << base_name << "_" << std::setw(2) << std::setfill('0') << curr_level;
    parmoon_db["output_basename"].set(ostr.str(), false);
    
    ConvectionDiffusion<2> cd2d(domain, parmoon_db);
    cd2d.assemble();
    cd2d.solve();
    
    estimator.estimate(cd2d.get_example(), cd2d.get_function());
    
    estimator.info();
    
    
    ConvectionDiffusion<2>::FEFunction estimated_error(estimator, values);
    cd2d.add_to_output(&estimated_error);
    cd2d.output();
    
    auto fespace = cd2d.get_space();
    
    if(estimator.get_estimated_global_error()<1.0e-10 || curr_level == n_adaptive_steps)
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
  TDatabase::destroy();
  Output::close_file();
  return 0;
} // end main
