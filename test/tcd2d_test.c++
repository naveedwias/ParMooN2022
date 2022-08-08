/**
 * @brief A test program to test Time_CD2D program
 */
#include <cmath>

#include "all_defines_external_libraries.h"
#include <Domain.h>
#include <Database.h>
#include "TimeConvectionDiffusion.h"
#include <TimeDiscRout.h>
#include <TimeDiscretizations.h>
#include <MainUtilities.h>
#include <ConvDiff.h>
#include <AuxParam2D.h>
#include "ParMooN.h"
#include <list>


void testBE()
{
  
}

void testCN(TimeConvectionDiffusion<2> &tcd, int m)
{
  ParameterDatabase db = tcd.get_db();
  double errors[5];
  errors[0]=errors[1]=errors[2]=errors[3]=errors[4]=0.;
  TAuxParam2D aux;
  MultiIndex2D AllDerivatives[3] = {MultiIndex2D::D00, MultiIndex2D::D10,
                                    MultiIndex2D::D01};
  const TFEFunction2D& function = tcd.get_function();
  
  tcd.output();
  const TFESpace2D* space = function.GetFESpace2D().get();

  function.GetErrors(tcd.get_example().get_exact(0), 3, AllDerivatives, 5,
                     conv_diff_l2_h1_linf_error<2>,
                     tcd.get_example().get_coeffs(), &aux, 1, &space, errors);
  Output::print("output");
  if(m==1){cout << errors[0] << "  "<<
  errors[1]<< "  " <<endl;}
  double eps = 1e-10;

  const int nb_test = 7;
  int iter_m[nb_test] = {0, 1, 2, 3, 18, 19, 20};

  double L2_expect[nb_test] = {0.};
  double H1_expect[nb_test] = {0.};

  double L2_galerkin[nb_test] = {0.0033601882323796, 0.0021831948624723,
                                 0.0041373036975832, 0.0051536220684535,
                                 0.00084798363062268, 0.0011485020220284,
                                 0.0017171969767027};
  double H1_galerkin[nb_test] = {0.25214923985069, 0.34294475156715,
                                 0.45320721802076, 0.56543341117647,
                                 0.13988043445394, 0.18484695459687,
                                 0.25175783286665};

  double L2_supg[nb_test] = {0.0033601882323796, 0.0022877555421121,
                             0.0042674492028728, 0.0052925302485841,
                             0.00090066884841505, 0.0012154988855681,
                             0.0018006391409988};
  double H1_supg[nb_test] = {0.25214923985069, 0.34297576782859,
                             0.45325543222406, 0.56548475252577,
                             0.13989349915631, 0.18486434460595,
                             0.25178153241465};

  if(db["space_discretization_type"].is("galerkin"))
  {
    std::copy(std::begin(L2_galerkin), std::end(L2_galerkin),
              std::begin(L2_expect));
    std::copy(std::begin(H1_galerkin), std::end(H1_galerkin),
              std::begin(H1_expect));
  }
  else if(db["space_discretization_type"].is("supg"))
  {
    std::copy(std::begin(L2_supg), std::end(L2_supg),
              std::begin(L2_expect));
    std::copy(std::begin(H1_supg), std::end(H1_supg),
              std::begin(H1_expect));
  }

  for(int i=0 ; i<nb_test ; i++)
  {
    if(m==iter_m[i])
    {
      cout << "M: " << m << endl;
      if( std::abs(errors[0] - L2_expect[i]) > eps )
      {
        ErrThrow("test tcd2d for ",
                 tcd.get_solver().get_db()["solver_type"].get<std::string>(),
                 " solver with ", db["time_discretization"].get<std::string>(),
                 " and ", db["space_discretization_type"].get<std::string>(),
                 ": L2 norm not correct, ", std::setprecision(14), errors[0],
                 " != ", std::setprecision(14), L2_expect[i]);
      }
      if( std::abs(errors[1] - H1_expect[i]) > eps )
      {
        ErrThrow("test tcd2d for ",
                 tcd.get_solver().get_db()["solver_type"].get<std::string>(),
                 " solver with ", db["time_discretization"].get<std::string>(),
                 " and ", db["space_discretization_type"].get<std::string>(),
                 ": H1 norm not correct, ", std::setprecision(14), errors[1],
                 " != ", std::setprecision(14), H1_expect[i]);
      }
    }
  }
}

void time_integration(TimeConvectionDiffusion<2>& tcd)
{
  TimeDiscretization& tss = tcd.get_time_stepping_scheme();
  tss.current_step_ = 0;
  tss.set_time_disc_parameters();
  
  TDatabase::TimeDB->CURRENTTIME = tss.get_start_time();
  
  tcd.assemble_initial_time();
  
  testCN(tcd, tss.current_step_);

  while(!tss.reached_final_time_step())
  {
    tss.current_step_++;
    SetTimeDiscParameters(1);

    tss.current_time_ += tss.get_step_length();;
    TDatabase::TimeDB->CURRENTTIME += tss.get_step_length();

    Output::print<1>("\nCURRENT TIME: ", tss.current_time_);
    tcd.assemble();
    tcd.solve();
    testCN(tcd, tss.current_step_);
  }
  
}

int main(int, char**)
{
  parmoon::parmoon_initialize();
  // test with Crank Nicolson
  {
    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db.merge(Example2D::default_example_database());
    db.merge(LocalAssembling2D::default_local_assembling_database());
    db.merge(TimeDiscretization::default_TimeDiscretization_database());
    db["example"] = 0;
    db["reynolds_number"] = 1;

    db["time_discretization"] = "crank_nicolson";
    db["time_step_length"] = 0.05;
    db["time_end"]=1.;
    TDatabase::ParamDB->ANSATZ_ORDER=1;

    // declaration of databases
    db.add("boundary_file", "Default_UnitSquare", "");
    db.add("geo_file", "UnitSquare", "", {"UnitSquare", "TwoTriangles"});
    db.add("refinement_n_initial_steps", (size_t) 5,"");
    TDomain domain(db);
    // refine grid
    domain.refine_and_get_hierarchy_of_collections(db);

    // test direct solver with Galerkin
    Output::print("\n\nTesting Galerkin\n");
    db.add("solver_type", "direct", "", {"direct", "petsc"});
    db["space_discretization_type"] = "galerkin";
    TimeConvectionDiffusion<2> tcd(domain, db);
    Output::print("A");
    time_integration(tcd);
    Output::print("B");
    
#ifdef PARMOON_WITH_PETSC
    // test PETSc solver
    Output::print("\n\nTesting PETSc\n");
    db["solver_type"] = "petsc";
    TDatabase::TimeDB->CURRENTTIME = db["time_start"];
    TimeConvectionDiffusion<2> tcd_petsc(domain, db);
    time_integration(tcd_petsc);
#endif // PARMOON_WITH_PETSC
    
    // test direct solver with SUPG
    Output::print("\n\nTesting SUPG\n");
    db["solver_type"] = "direct";
    db["space_discretization_type"] = "supg";
    TDatabase::TimeDB->CURRENTTIME = db["time_start"];
    TDatabase::ParamDB->DELTA0 = 1.;
    TimeConvectionDiffusion<2> tcd_supg(domain, db);
    time_integration(tcd_supg);
  }
  parmoon::parmoon_finalize();
}
 
