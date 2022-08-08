/**
 * @brief A test program for the solving of Algebraic Flux Correction (AFC) in CD2D problems
 * using a Posteriori error estimation. 
 * Reference : Jha2020, arXiv, https://arxiv.org/abs/2005.02938
 *
 * This serves as a test for the solving of a posteriroi error porblem with AFC schemes applied to 
 * CD2D problems. It is intended to perform CD2D calculations with different examples in different 
 * setups to test a wide variety of ParMooN core functionality.
 *
 * The norms of the solution are compared with reference norms.
 * If those are not approximated well enough (or something in the process goes wrong)
 * the test fails.
 *
 * Should this test fail, there are two possibilities: either you made a mistake
 * which broke the programs functionality. Then you must find the mistake.
 * Or you changed some program setup (e.g. changed the default solver). Then this tests
 * shows you how many other program parts are affected by your changes.
 * If you are not perfectly sure how to repair this, it is a good idea
 * to describe your changes in the forum and request support.
 * 
 * NOTE: Tests can be extended to Q_1 elements as well.
 * 
 * If any problems are faced contact the author.
 *
 * @date 2020/07/10
 * @author Abhinav Jha
 * @email jha.abhinav0207@gmail.com
 *
 */
#include <cmath>
#include "AlgebraicFluxCorrection.h"
#include "Domain.h"
#include "Database.h"
#include "LoopInfo.h"
#include "ConvectionDiffusion_AFC.h"
#include "CDErrorEstimator_AFC.h"
#include "CDErrorEstimator.h"
#include "RefinementStrategy.h"
#include "ParMooN.h"

#include "LocalAssembling.h"
#include <MainUtilities.h> //for error measuring

void compareErrors(const ConvectionDiffusion<2>& cd2d,
                   std::array<double, 4> errors)
{
  const double eps = 1e-12;

  // check the errors
  if( std::abs(cd2d.get_L2_error() - errors[0]) > eps )
  {
    ErrThrow("Program 1: L2 error not correct. ",
             cd2d.get_L2_error() - errors[0]);
  }
  if( std::abs(cd2d.get_H1_semi_error() - errors[1]) > eps )
  {
    ErrThrow("Program 1: H1-semi error not correct. ",
             cd2d.get_H1_semi_error() - errors[1]);
  }
  if( std::abs(cd2d.get_SD_error() - errors[2]) > eps )
  {
    ErrThrow("Program 1: sd error not correct.",
             cd2d.get_SD_error() - errors[2]);
  }
  if( std::abs(cd2d.get_L_inf_error() - errors[3]) > eps )
  {
    ErrThrow("Program 1: L_inf error not correct.",
             cd2d.get_L_inf_error() - errors[3]);
  }
}

// Here the actual computations take place
void check_cd2d(TDomain & domain, ParameterDatabase& parmoon_db, std::array<double,4> errors)
{        
  CDErrorEstimator<2> estimator_initial(parmoon_db);  
  
  if(parmoon_db["afc_estimator_type"])
      parmoon_db["estimator_type"] = 1;
  CDErrorEstimator_AFC<2> estimator(parmoon_db);
  RefinementStrategy<2> refinementStrategy(parmoon_db);
  LoopInfo<double> loop_info("adaptive", true, true, 1);
  BlockVector values;
  int n_adaptive_steps = (int)domain.get_database()["refinement_max_n_adaptive_steps"]
                           + (int)domain.get_database()["refinement_max_n_uniform_steps"];
  int n_uniform_steps = 0;
  bool adaptive_converged = false;
  for(int curr_level = 0; ; ++curr_level)
  {
    Output::print("\nadaptive loop ", curr_level);    
    ConvectionDiffusion_AFC<2> cd2d(domain, parmoon_db);
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
    estimator.estimate(cd2d.get_example(), cd2d.get_function(), 
                       cd2d.get_afc_D_entries(), cd2d.get_afc_alphas(),
                       cd2d.get_initial_fe_function());
    if(parmoon_db["afc_estimator_type"])
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
    
    if(estimator.get_estimated_global_error()<1.0e-3 || curr_level == n_adaptive_steps)
    {
      adaptive_converged = true;
    }
    if(adaptive_converged)
    {
      loop_info.finish(curr_level, estimator.get_estimated_global_error());
      compareErrors(cd2d, errors); // throws upon a difference
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
    }
  }
}

// =======================================================================
// main program
// =======================================================================
int main(int, char* argv[])
{
  //testall argument
  bool testall = false;
  if(argv[1])
    testall = (std::string(argv[1]).compare("testall") == 0);
  ParameterDatabase db = parmoon::parmoon_initialize();
  Output::setVerbosity(2);
  
  Output::print("\ntesting with algebraic flux correction");
  
  size_t nRefinements = 2; // initial refinements prior to adaptive refinements
  
  db.merge(Example2D::default_example_database());
  db.merge(Solver<>::default_solver_database());
  db.merge(Example2D::default_example_database());
  db.merge(ParameterDatabase::default_output_database());
  db.merge(RefinementStrategy<2>::default_refinement_strategy_database());
  db.merge(TDomain::default_domain_parameters());
  db.merge(AlgebraicFluxCorrection::default_afc_database(), true);
  db["example"] = 8; //Example 1 ABR17
  
  //only direct solver for CD2D class
  db["solver_type"].set<>("direct");
  db["refinement_n_initial_steps"] = nRefinements;
  
  //addition for AFC database
  db["algebraic_flux_correction"].set("afc");
  db["afc_iteration_scheme"].set<>("fixed_point_rhs");
  db["afc_initial_iterate"].set<>("supg");
  
  //Diffusion Coeffecient
  TDatabase::ParamDB->RE_NR = 1000;
  db["diffusion_coefficient"] = 1.0e-3;
  
  //Uniform and Adaptive steps 
  db["refinement_max_n_adaptive_steps"] = 3;
  db["refinement_max_n_uniform_steps"] = 1;
  
  //Estimator Type
  db.add("estimator_type", (size_t) 0, "", (size_t) 0,(size_t) 6);
  
  // default construct a domain object
  db["boundary_file"].set<>("Default_UnitSquare");
  db["geo_file"].set<>("TwoTriangles");
  TDomain domain(db);
  
  //AFC only applicable to P1 elements
  TDatabase::ParamDB->ANSATZ_ORDER = 1; //P1 elements
  
  //maximum number of iterations for non linear loop
  db["afc_nonlinloop_maxit"]=1000;
 
  // refine grid up to the coarsest level
  domain.refine_and_get_hierarchy_of_collections(db);
  
  //Computation for P1 elements
  //================================================================//
  //P1 elements, conforming closure, with AFC-Energy technique
  Output::print("\n\n --------- P1 + kuzmin + AFC Energy + Conforming Closure ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["conforming_closure"] = true;
  db["estimator_type"] = (size_t) 0;
  //AFC Estimator Type chooses between AFC Energy and AFC SUPG Energy
  db.add("afc_estimator_type", false, "");
  std::array<double,4> errors = {{0.012317705871771, 3.8201374122139,
                                  0.57090416534603, 0.19029305056368}};
  check_cd2d(domain, db, errors);
  
  //P1 elements, conforming closure, with AFC-SUPG-Energy technique
  Output::print("\n\n --------- P1 + kuzmin + AFC SUPG Energy + Conforming Closure ---------\n");
  // Reset Domain
  db["boundary_file"].set<>("Default_UnitSquare");
  db["geo_file"].set<>("TwoTriangles");
  TDomain domain_afcse_conform_tria_kuzmin(db);
   // refine grid up to the coarsest level
  domain_afcse_conform_tria_kuzmin.refine_and_get_hierarchy_of_collections(db);
  db["afc_limiter"].set<>("kuzmin");
  db["conforming_closure"] = true;
  db["estimator_type"] = (size_t) 5;
  db["afc_estimator_type"] = true;
  errors = {{0.012115462769907, 3.8204282913751, 0.55153400604398, 0.1902404414289}};
  check_cd2d(domain_afcse_conform_tria_kuzmin, db, errors);
  
  if(testall)
  {
    //P1 elements, conforming closure, with AFC-SUPG-Energy technique
    Output::print("\n\n --------- P1 + BJK17 + AFC Energy + Conforming Closure ---------\n");
    // Reset Domain
    db["boundary_file"].set<>("Default_UnitSquare");
    db["geo_file"].set<>("TwoTriangles");
    TDomain domain_afce_conform_tria_bjk(db);
    // refine grid up to the coarsest level
    domain_afce_conform_tria_bjk.refine_and_get_hierarchy_of_collections(db);
    db["afc_limiter"].set<>("BJK17");
    db["conforming_closure"] = true;
    db["estimator_type"] = (size_t) 0;
    db["afc_estimator_type"] = false;
    errors = {{0.012653172069456, 3.8328630940764, 0.61814150791135, 0.1904324387886}};
    check_cd2d(domain_afce_conform_tria_bjk, db, errors);
    
    //P1 elements, conforming closure, with AFC-SUPG-Energy technique
    Output::print("\n\n --------- P1 + BJK17 + AFC SUPG Energy + Conforming Closure ---------\n");
    // Reset Domain
    db["boundary_file"].set<>("Default_UnitSquare");
    db["geo_file"].set<>("TwoTriangles");
    TDomain domain_afcse_conform_tria_bjk(db);
    // refine grid up to the coarsest level
    domain_afcse_conform_tria_bjk.refine_and_get_hierarchy_of_collections(db);
    db["afc_limiter"].set<>("BJK17");
    db["conforming_closure"] = true;
    db["estimator_type"] = (size_t) 5;
    db["afc_estimator_type"] = true;
    errors = {{0.011902356837195, 3.824773331893, 0.55261097047288, 0.19034433873272}};
    check_cd2d(domain_afcse_conform_tria_bjk, db, errors);
  }
  
  //Computation for P1 elements
  //================================================================//
  //P1 elements, Hanging nodes, with AFC-SUPG-Energy technique
  Output::print("\n\n --------- P1 + kuzmin + AFC Energy + Hanging Nodes ---------\n");
  // Reset Domain
  db["boundary_file"].set<>("Default_UnitSquare");
  db["geo_file"].set<>("TwoTriangles");
  db["conforming_closure"] = false;
  TDomain domain_afce_hang_tria_kuzmin(db);
   // refine grid up to the coarsest level
  domain_afce_hang_tria_kuzmin.refine_and_get_hierarchy_of_collections(db);
  db["afc_limiter"].set<>("kuzmin");
  db["estimator_type"] = (size_t) 0;
  db["afc_estimator_type"] = false;
  errors = {{0.013997998707125, 3.82377524303029, 0.60144111077041, 0.1905898680719}};
  check_cd2d(domain_afce_hang_tria_kuzmin, db, errors);
  
  if(testall)
  {
   //P1 elements, Hanging nodes, with AFC-SUPG-Energy technique
   Output::print("\n\n --------- P1 + kuzmin + AFC SUPG Energy + Hanging Nodes ---------\n");
   // Reset Domain
   db["boundary_file"].set<>("Default_UnitSquare");
   db["geo_file"].set<>("TwoTriangles");
   db["conforming_closure"] = false;
   TDomain domain_afcse_hang_tria_kuzmin(db);
    // refine grid up to the coarsest level
   domain_afcse_hang_tria_kuzmin.refine_and_get_hierarchy_of_collections(db);
   db["afc_limiter"].set<>("kuzmin");
   db["estimator_type"] = (size_t) 5;
   db["afc_estimator_type"] = true;
   errors = {{0.0142978522048055,3.86797860606598,0.562221910102313,0.190556616233352}};
   check_cd2d(domain_afcse_hang_tria_kuzmin, db, errors);
   
   //P1 elements, Hanging Nodes, with AFC-SUPG-Energy technique
   Output::print("\n\n --------- P1 + BJK17 + AFC Energy + Hanging Nodes ---------\n");
   // Reset Domain
   db["boundary_file"].set<>("Default_UnitSquare");
   db["geo_file"].set<>("TwoTriangles");
   db["conforming_closure"] = false;
   TDomain domain_afce_hang_tria_bjk(db);
    // refine grid up to the coarsest level
   domain_afce_hang_tria_bjk.refine_and_get_hierarchy_of_collections(db);
   db["afc_limiter"].set<>("BJK17");
   db["estimator_type"] = (size_t) 0;
   db["afc_estimator_type"] = false;
   errors = {{0.0153759321928701,3.84720493700289,0.658724745253243,0.190534432621869}};
   check_cd2d(domain_afce_hang_tria_bjk, db, errors);
   
   //P1 elements, Hanging nodes, with AFC-SUPG-Energy technique
   Output::print("\n\n --------- P1 + BJK17 + AFC SUPG Energy + Hanging Nodes ---------\n");
   // Reset Domain
   db["boundary_file"].set<>("Default_UnitSquare");
   db["geo_file"].set<>("TwoTriangles");
   db["conforming_closure"] = false;
   TDomain domain_afcse_hang_tria_bjk(db);
    // refine grid up to the coarsest level
   domain_afcse_hang_tria_bjk.refine_and_get_hierarchy_of_collections(db);
   db["afc_limiter"].set<>("BJK17");
   db["estimator_type"] = (size_t) 5;
   db["afc_estimator_type"] = true;
   errors = {{0.0141131632314105,3.87516806724602,0.564300824376585,0.190455374377342}};
   check_cd2d(domain_afcse_hang_tria_bjk, db, errors);
  }
  parmoon::parmoon_finalize();
}
