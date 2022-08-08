#include "BaseCell.h"
#include "BlockVector.h"
#include "Domain.h"
#include "Database.h"
#include <ParMooN_repository_info.h>
#include "MooNMD_Io.h"
#include "MainUtilities.h"
#include "Utilities.h"
#ifdef __3D__
#include "AuxParam3D.h"
#include "FEFunction3D.h"
using FESpace = TFESpace3D;
using FEFunction = TFEFunction3D;
auto zero_solution = unknown_solution_3d;
constexpr int d = 3;
#else
#include "AuxParam2D.h"
#include "FEFunction2D.h"
#include "FEVectFunct2D.h"
using FESpace = TFESpace2D;
using FEFunction = TFEFunction2D;
auto zero_solution = unknown_solution_2d;
constexpr int d = 2;
#endif
#include <random>
#include "DataWriter.h"
#include "ParMooN.h"

#include "GridTransfer.h"

#ifdef _MPI
#include <ParFECommunicator3D.h>
#endif

const std::string path_to_repo = parmoon::source_directory;
const std::string path_to_meshes = path_to_repo + "/data/mesh/";
constexpr double tolerance = 1.e-10;

// the TestObject and the methods prepare_test* are mostly a copy taken from
// collection_test.c++
struct TestObject
{
    std::string Geo_file;
    std::string Boundary_file;
    double l2norm;
    double h1norm;
    double min;
    double max;
    double osc_mean;
    double min_nodal_fctn;
    double max_nodal_fctn;
    double osc_mean_nodal_fctn;
    std::pair<double, double> minmax;
    std::pair<double, double> minmax_pk;
    std::pair<parmoon::Point, double> value; // function value at a given point
};

TestObject prepare_test0()
{
#ifdef __2D__
  std::string geo_file = path_to_meshes + "unit_square/unit_square_quad6.mesh";
  std::string bd_file = path_to_meshes + "unit_square/unit_square.PRM";
  std::pair<double, double> minmax = std::make_pair(0.45851245101479, 0.60861542083476);
  std::pair<double, double> minmax_pk = std::make_pair(0.45851245101479, 
0.60861542083476);
  auto point_eval = std::make_pair<parmoon::Point, double>(
    parmoon::Point(0.4, 0.5), 0.9);
#else
  std::string geo_file = "Default_UnitCube_Tetra";
  std::string bd_file = "Default_UnitCube";
  std::pair<double, double> minmax = std::make_pair(0, 1);
  std::pair<double, double> minmax_pk = std::make_pair(0, 1);
  auto point_eval = std::make_pair<parmoon::Point, double>(
    parmoon::Point(0.4, 0.5, 0.6), 0.9);
#endif
  double l2norm = 1.080123449734; // = std::sqrt(7./6.);
  double h1norm = 1.4142135623731; // = std::sqrt(2);
  double min = 0.;
  double max = 2.;
  double osc_mean = d == 2 ? 0.20732610381507 : 0.54166666666667;
  double min_nodal_fctn = 0.;
  double max_nodal_fctn = 2.;
  double osc_mean_nodal_fctn = d == 2 ? 0.20732610381507 : 0.54166666666667;
  return TestObject{geo_file, bd_file, l2norm, h1norm, min, max, osc_mean,
    min_nodal_fctn, max_nodal_fctn, osc_mean_nodal_fctn, minmax, minmax_pk,
    point_eval};
}
TestObject prepare_test1()
{
#ifdef __2D__
  std::string geo_file = path_to_meshes + "unit_square/unit_square_tria6.mesh";
  std::string bd_file = path_to_meshes + "unit_square/unit_square.PRM";
  std::pair<double, double> minmax = std::make_pair(1.3538728796251, 1.5000000000021);
  std::pair<double, double> minmax_pk = std::make_pair(1.3538728796251, 
1.5000000000021);
  auto point_eval = std::make_pair<parmoon::Point, double>(
    parmoon::Point(0.4, 0.5), 0.9);
#else
  std::string geo_file = "Default_UnitCube_Hexa";
  std::string bd_file = "Default_UnitCube";
  std::pair<double, double> minmax = std::make_pair(0, 1);
  std::pair<double, double> minmax_pk = std::make_pair(0, 1);
  auto point_eval = std::make_pair<parmoon::Point, double>(
    parmoon::Point(0.4, 0.5, 0.6), 0.9);
#endif
  double l2norm = 1.080123449734; // = std::sqrt(7./6.);
  double h1norm = 1.4142135623731; // = std::sqrt(2);
  double min = 0.;
  double max = 2.;
  double osc_mean = d == 2 ? 0.19820922566297 : 0.5;
  double min_nodal_fctn = 0.;
  double max_nodal_fctn = 2.;
  double osc_mean_nodal_fctn = d == 2 ? 0.19820922566297 : 0.5;
  return TestObject{geo_file, bd_file, l2norm, h1norm, min, max, osc_mean,
    min_nodal_fctn, max_nodal_fctn, osc_mean_nodal_fctn, minmax, minmax_pk,
    point_eval};
}
TestObject prepare_test2()
{
#ifdef __2D__
  std::string geo_file = path_to_meshes + "Hemker_tria.mesh";
  std::string bd_file = path_to_meshes + "Hemker.PRM";
  double l2norm = 41.550323206975;
  double h1norm = 11.735293970182;
  double min = -6.;
  double max = 12.;
  double osc_mean = 3.0270082488921;
  double min_nodal_fctn = -6;
  double max_nodal_fctn = 12;
  double osc_mean_nodal_fctn = 3.0270082488921;
  std::pair<double, double> minmax = std::make_pair(-6, -5.05);
  std::pair<double, double> minmax_pk = std::make_pair(-6, -5.05);
  auto point_eval = std::make_pair<parmoon::Point, double>(
    parmoon::Point(3.0, 1.5), 4.5);
#else
  std::string geo_file = path_to_meshes + "cylinder.3d.3K.mesh";
  std::string bd_file = "";
  double l2norm = 15.751706293257;
  double h1norm = 15.802393311929;
  double min = -2.828428;
  double max = 2.828428;
  double osc_mean = 0.92345979102845;
  double min_nodal_fctn = -2.828428;
  double max_nodal_fctn = 2.828428;
  double osc_mean_nodal_fctn = 0.92345979102845;
  std::pair<double, double> minmax = std::make_pair(-0.0804295, 0.168075);
  std::pair<double, double> minmax_pk = std::make_pair(-0.0804295, 0.168075);
  auto point_eval = std::make_pair<parmoon::Point, double>(
    parmoon::Point(0.2, 1.5, 7.0), 1.7);
#endif
  return TestObject{geo_file, bd_file, l2norm, h1norm, min, max, osc_mean,
    min_nodal_fctn, max_nodal_fctn, osc_mean_nodal_fctn, minmax, minmax_pk,
    point_eval};
}
TestObject prepare_test3()
{
#ifdef __2D__
  std::string geo_file = path_to_meshes + "Hemker_quad.mesh";
  std::string bd_file = path_to_meshes + "Hemker.PRM";
  double l2norm = 41.550323206975;
  double h1norm = 11.735293970182;
  double min = -6.;
  double max = 12.;
  double osc_mean = 3.1240751947804;
  double min_nodal_fctn = -6;
  double max_nodal_fctn = 12.;
  double osc_mean_nodal_fctn = 3.1240751947804;
  std::pair<double, double> minmax = std::make_pair(-6, -5.045);
  std::pair<double, double> minmax_pk = std::make_pair(-6, -5.045);
  auto point_eval = std::make_pair<parmoon::Point, double>(
    parmoon::Point(3.0, 1.5), 4.5);
#else
  std::string geo_file = path_to_meshes + "channel.3d.mesh";
  std::string bd_file = "";
  double l2norm = 5.7321647429222;
  double h1norm = 8.2736671434135;
  double min = -2.4;
  double max = 2.4;
  double osc_mean = 0.87176246410779;
  double min_nodal_fctn = -2.4;
  double max_nodal_fctn = 2.4;
  double osc_mean_nodal_fctn = 0.87176246410779;
  std::pair<double, double> minmax = std::make_pair(-2.4, -1.79999984375);
  std::pair<double, double> minmax_pk = std::make_pair(-2.4, -1.79999984375);
  auto point_eval = std::make_pair<parmoon::Point, double>(
    parmoon::Point(0.2, 0.5, 3.0), 0.7);
#endif
  return TestObject{geo_file, bd_file, l2norm, h1norm, min, max, osc_mean,
    min_nodal_fctn, max_nodal_fctn, osc_mean_nodal_fctn, minmax, minmax_pk,
    point_eval};
}

void fill_with_arbitrary_data(BlockVector& x)
{
  constexpr double min = -5.;
  constexpr double max = 5.;
  int length = x.length();
  // initialize seed, we want reproducibility rather than randomness.
  std::mt19937 e2(length);
  std::uniform_real_distribution<> dist(min, max);
  for(int i = 0; i < length; ++i)
  {
    x[i] = dist(e2);
  }
}

void compare_with_get_errors(const FEFunction& function)
{
  auto norms = function.get_L2_and_H1_norm();

#ifdef __3D__
  int n_derivatives = 4;
  TAuxParam3D aux;
  MultiIndex3D all_derivs[] = {MultiIndex3D::D000, MultiIndex3D::D100,
                               MultiIndex3D::D010, MultiIndex3D::D001};
  const FESpace* space_ptr = function.GetFESpace3D().get();
#else
  int n_derivatives = 3;
  TAuxParam2D aux;
  MultiIndex2D all_derivs[] = {MultiIndex2D::D00, MultiIndex2D::D10,
                               MultiIndex2D::D01};
  const FESpace* space_ptr = function.GetFESpace2D().get();
#endif
  double errors[4] = {};
  function.GetErrors(zero_solution, n_derivatives, all_derivs, 3, L2H1Errors,
                      nullptr, &aux, 1, &space_ptr, errors);
#ifdef _MPI
  // for some reason, the GetErrors method does not communicate
  double sendbuf[2] = {errors[0], errors[1]};
  double recvbuf[2] = {0.0, 0.0};
  MPI_Allreduce(sendbuf, recvbuf, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  errors[0] = std::sqrt(recvbuf[0]);
  errors[1] = std::sqrt(recvbuf[1]);
#endif
  if(!utilities::are_equal(norms.first, errors[0]))
  {
    ErrThrow("difference in L2 error: ", setprecision(20), norms.first, " ", errors[0], ", ", std::abs(norms.first - errors[0]));
  }
  if(!utilities::are_equal(norms.second, errors[1]))
  {
    ErrThrow("difference in H1 error: ", setprecision(20), norms.second, " ", errors[1], ", ", std::abs(norms.second - errors[1]));
    ErrThrow("difference in H1 error: ", norms.second - errors[1]);
  }
}

void analytic_function(double x, double y,
#ifdef __3D__
                       double /*z*/,
#endif
                       double *vals)
{
  vals[0] = x+y;
}

void analyticFunction2(double x, double y,
#ifdef __3D__
                       double z,
#endif                    
                       double* vals)
{
//   vals[0] = x + y;
  vals[0] = x*x + 2.0*y*y;
  
  #ifdef __3D__
  vals[0] += 3.0*z*z;
  #endif
}


void errorIfToDifferent(const BlockVector& computed, const std::vector<double>& expected)
{
  const double tolerance = 1.0e-7;
  
  if(computed.length() != expected.size())
    ErrThrow("BlockVectors computed and expected have different lengths!");
  
  for(unsigned int i = 0; i < expected.size(); i++)
  {
    double diff = std::abs(computed[i] - expected.at(i));
    if( diff > tolerance )
      ErrThrow("Component i = ", i, " is ", computed[i], ", but it should be ", expected.at(i));
  }
}

void errorIfToDifferent(const double computed, const double expected)
{
  const double tolerance = 1.0e-7;
  
    double diff = std::abs(computed - expected);
    if( diff > tolerance )
      ErrThrow("Computed number is ", computed, ", but it should be ", expected);
}

#ifdef __2D__
void scalarFunctionTest2D(std::shared_ptr<const TFESpace2D>& spaceP2,
                          std::shared_ptr<const TFESpace2D>& spaceQ0,
                          std::shared_ptr<const TFESpace2D>& spaceP1)
{
  Output::print("\nRunning Interpolation test 2D, scalar.\n");
  
  // preparation of P2 function
  
  BlockVector function2Values(spaceP2->get_n_dof());
  TFEFunction2D function2(spaceP2, "f2", function2Values.get_entries());
  function2.Interpolate(analyticFunction2);
  
  
  // test of projection from P2 onto Q0
  
  BlockVector function0Values(spaceQ0->get_n_dof());
  TFEFunction2D function0(spaceQ0, "f0", function0Values.get_entries());
  
  function0.Interpolate(&function2);
  
  std::vector<double> function0Expected = {0.1875, 0.6875, 1.6875, 1.1875};
  errorIfToDifferent(function0Values, function0Expected);
  
  
  // test of projection from Q0 onto P1
  
  BlockVector function1Values(spaceP1->get_n_dof());
  TFEFunction2D function1(spaceP1, "f1", function1Values.get_entries());
  
  function1.Interpolate(&function0);
  
  std::vector<double> function1Expected = {0.9375, 0.1875, 0.4375, 0.6875, 1.1875,
                                           1.6875, 0.6875, 1.4375, 1.1875};
  errorIfToDifferent(function1Values, function1Expected);
}

void vectorFunctionTest2D(std::shared_ptr<const TFESpace2D>& spaceP2,
                          std::shared_ptr<const TFESpace2D>& spaceQ0,
                          std::shared_ptr<const TFESpace2D>& spaceP1)
{
  Output::print("\nRunning Interpolation test 2D, vector.\n");
  
  // preparation of P2 function
  
  unsigned int dofsP2 = spaceP2->get_n_dof();
  BlockVector function2VecValues({dofsP2, dofsP2});
  TFEFunction2D function2Aux(spaceP2, "f2Aux",
                             function2VecValues.get_entries());
  function2Aux.Interpolate(analyticFunction2);
  
  TFEVectFunct2D function2Vec(spaceP2, "f2Vec", function2VecValues.get_entries(),
                              2);
  
  for(unsigned int i = 0; i < dofsP2; i++)
    function2VecValues[dofsP2 + i] = - function2VecValues[i];
  
  
  // test of projection from P2 onto Q0
  
  unsigned int dofsQ0 = spaceQ0->get_n_dof();
  BlockVector function0VecValues({dofsQ0, dofsQ0});
  TFEVectFunct2D function0Vec(spaceQ0, "f0Vec", 
                              function0VecValues.get_entries(), 2);
  
  function0Vec.Interpolate(&function2Vec);
  
  std::vector<double> function0VecExpected = {0.1875, 0.6875, 1.6875, 1.1875,
                                              - 0.1875, - 0.6875, - 1.6875, - 1.1875};
  errorIfToDifferent(function0VecValues, function0VecExpected);
  
  
  // test of projection from Q0 onto P1
  
  unsigned int dofsP1 = spaceP1->get_n_dof();
  BlockVector function1VecValues({dofsP1, dofsP1});
  TFEVectFunct2D function1Vec(spaceP1, "f1Vec", 
                              function1VecValues.get_entries(), 2);
  
  function1Vec.Interpolate(&function0Vec);
  
  std::vector<double> function1VecExpected = {0.9375, 0.1875, 0.4375, 0.6875, 1.1875,
                                              1.6875, 0.6875, 1.4375, 1.1875,
                                              - 0.9375, - 0.1875, - 0.4375, - 0.6875, - 1.1875,
                                              - 1.6875, - 0.6875, - 1.4375, - 1.1875};
  errorIfToDifferent(function1VecValues, function1VecExpected);
}



void interpolationTests2D(ParameterDatabase& parmoon_db)
{
  Output::print("\nRunning Interpolation test 2D.\n");
  
  // preparation of grid collection
  
  parmoon_db["boundary_file"].set("Default_UnitSquare", false);
  parmoon_db["geo_file"].set("UnitSquare", false);
  parmoon_db["refinement_n_initial_steps"] = 2;
  
  TDomain domainQ(parmoon_db);
  auto collectionsQ = domainQ.refine_and_get_hierarchy_of_collections(parmoon_db);
  std::list<TCollection*>::iterator itQ = collectionsQ.begin();
  
  std::shared_ptr<const TFESpace2D> spaceQ0Fine{ new TFESpace2D(*itQ, std::string("spaceQ0Fine"),
                                                 BoundConditionNoBoundCondition, 0) };
                                                 
  std::shared_ptr<const TFESpace2D> spaceQ1Fine{ new TFESpace2D(*itQ, std::string("spaceQ1Fine"),
                                                 BoundConditionNoBoundCondition, 1) };
  
  std::advance(itQ, 1);
  
  parmoon_db["boundary_file"].set("Default_UnitSquare", false);
  parmoon_db["geo_file"].set("TwoTriangles", false);
  parmoon_db["refinement_n_initial_steps"] = 2;
  
  TDomain domainP(parmoon_db);
  auto collectionsP = domainP.refine_and_get_hierarchy_of_collections(parmoon_db);
  std::list<TCollection*>::iterator itP = collectionsP.begin();
  std::advance(itP, 1);
  
  
  // preparation of spaces
  
  std::shared_ptr<const TFESpace2D> spaceP2{ new TFESpace2D(*itP, std::string("spaceP2"),
                                            BoundConditionNoBoundCondition, 2) };
  
  std::shared_ptr<const TFESpace2D> spaceQ0{ new TFESpace2D(*itQ, std::string("spaceQ0"),
                                         BoundConditionNoBoundCondition, 0) };
  
  std::shared_ptr<const TFESpace2D> spaceP1{ new TFESpace2D(*itP, std::string("spaceP1"),
                                         BoundConditionNoBoundCondition, 1) };
  
  std::shared_ptr<const TFESpace2D> spaceQ1{ new TFESpace2D(*itQ, std::string("spaceQ1"),
                                         BoundConditionNoBoundCondition, 1) };
                                         
  // list of tests
                                         
  scalarFunctionTest2D(spaceP2, spaceQ0, spaceP1);
  vectorFunctionTest2D(spaceP2, spaceQ0, spaceP1);
}
#endif


#ifdef __3D__
void scalarFunctionTest3D(std::shared_ptr<const TFESpace3D>& spaceP2,
                          std::shared_ptr<const TFESpace3D>& spaceQ0,
                          std::shared_ptr<const TFESpace3D>& spaceP1)
{
  Output::print("\nRunning Interpolation test 3D, scalar.\n");
  
  // preparation of P2 function
  
  BlockVector function2Values(spaceP2->get_n_dof());
  TFEFunction3D function2(spaceP2, "f2", function2Values.get_entries());
  function2.Interpolate(analyticFunction2);
  
  
  // test of projection from P2 onto Q0
  
  BlockVector function0Values(spaceQ0->get_n_dof());
  TFEFunction3D function0(spaceQ0, "f0", function0Values.get_entries());
  
  function0.Interpolate(&function2);
  
  
  double l2Norm = function0.get_L2_norm();
  double integral, measure;
  function0.compute_integral_and_measure(integral, measure);
  
  errorIfToDifferent(l2Norm, 2.09538182678);
  errorIfToDifferent(integral, 1.875);

#ifndef  _MPI  
  std::vector<double> function0Expected = {0.375, 0.875, 1.875, 1.375,
                                           1.875, 2.375, 3.375, 2.875};
  errorIfToDifferent(function0Values, function0Expected);
#endif  
  
  
  // test of projection from Q0 onto P1
  
  BlockVector function1Values(spaceP1->get_n_dof());
  TFEFunction3D function1(spaceP1, "f1", function1Values.get_entries());
  
  function1.Interpolate(&function0);
  
  
  l2Norm = function1.get_L2_norm();
  function1.compute_integral_and_measure(integral, measure);
  
  errorIfToDifferent(l2Norm, 1.9512282456613);
  errorIfToDifferent(integral, 1.875);

#ifndef  _MPI  
  std::vector<double> function1Expected = {1.875, 0.375, 0.625, 1.125, 0.875, 1.375, 2.125,
                                           1.875, 2.625, 3.375, 1.375, 2.375, 1.625, 2.875,
                                           0.875, 1.625, 1.375, 2.375, 1.625, 2.125, 2.875,
                                           3.125, 1.125, 1.875, 2.125, 2.625, 2.375};
  errorIfToDifferent(function1Values, function1Expected);
#endif  
}


void vectorFunctionTest3D(std::shared_ptr<const TFESpace3D>& spaceP2,
                          std::shared_ptr<const TFESpace3D>& spaceQ0,
                          std::shared_ptr<const TFESpace3D>& spaceP1)
{
  Output::print("\nRunning Interpolation test 3D, vector.\n");
  
  // preparation of P2 function
  
  unsigned int dofsP2 = spaceP2->get_n_dof();
  BlockVector function2VecValues({dofsP2, dofsP2, dofsP2});
  TFEFunction3D function2Aux(spaceP2, "f2Aux", 
                             function2VecValues.get_entries());
  function2Aux.Interpolate(analyticFunction2);
  
  TFEVectFunct3D function2Vec(spaceP2, "f2Vec", 
                              function2VecValues.get_entries(), 3);
  
  for(unsigned int i = 0; i < dofsP2; i++)
  {
    function2VecValues[dofsP2 + i] = - function2VecValues[i];
    function2VecValues[dofsP2*2 + i] = function2VecValues[i];
  }
  
  
  // test of projection from P2 onto Q0
  
  unsigned int dofsQ0 = spaceQ0->get_n_dof();
  BlockVector function0VecValues({dofsQ0, dofsQ0, dofsQ0});
  TFEVectFunct3D function0Vec(spaceQ0, "f0Vec", 
                              function0VecValues.get_entries(), 3);
  
  function0Vec.Interpolate(&function2Vec);
  
  std::vector<double> l2NormExpected = {2.09538182678, 2.09538182678,
                                        2.09538182678};
                              
  std::vector<double> integralExpected = {1.875, - 1.875, 1.875};
  
  for(int i = 0; i < 3; i++)
  {
    auto function0Vec_i = function0Vec.GetComponent(i);
    double l2Norm = function0Vec_i->get_L2_norm();
    
    double integral, measure;
    function0Vec_i->compute_integral_and_measure(integral, measure);
    
    errorIfToDifferent(l2Norm, l2NormExpected.at(i));
    errorIfToDifferent(integral, integralExpected.at(i));
  }
  
  
#ifndef  _MPI  
  std::vector<double> function0VecExpected = {0.375, 0.875, 1.875, 1.375,
                                              1.875, 2.375, 3.375, 2.875,
                                              - 0.375, - 0.875, - 1.875, - 1.375,
                                              - 1.875, - 2.375, - 3.375, - 2.875,
                                              0.375, 0.875, 1.875, 1.375,
                                              1.875, 2.375, 3.375, 2.875};
  errorIfToDifferent(function0VecValues, function0VecExpected);
#endif 
  
  
  // test of projection from Q0 onto P1
  
  unsigned int dofsP1 = spaceP1->get_n_dof();
  BlockVector function1VecValues({dofsP1, dofsP1, dofsP1});
  TFEVectFunct3D function1Vec(spaceP1, "f1Vec", 
                              function1VecValues.get_entries(), 3);
  
  function1Vec.Interpolate(&function0Vec);
  
  
  l2NormExpected = {1.9512282456613, 1.9512282456613, 1.9512282456613};
  integralExpected = {1.875, - 1.875, 1.875};
  
  for(int i = 0; i < 3; i++)
  {
    auto function1Vec_i = function1Vec.GetComponent(i);
    
    double l2Norm = function1Vec_i->get_L2_norm();
    double integral, measure;
    function1Vec_i->compute_integral_and_measure(integral, measure);
    
    errorIfToDifferent(l2Norm, l2NormExpected.at(i));
    errorIfToDifferent(integral, integralExpected.at(i));
  }
  
#ifndef  _MPI
  std::vector<double> oneBlockValues = {1.875, 0.375, 0.625, 1.125, 0.875, 1.375, 2.125,
                                        1.875, 2.625, 3.375, 1.375, 2.375, 1.625, 2.875,
                                        0.875, 1.625, 1.375, 2.375, 1.625, 2.125, 2.875,
                                        3.125, 1.125, 1.875, 2.125, 2.625, 2.375};
  
  std::vector<double> function1VecExpected(3*dofsP1);
  
  for(unsigned int i = 0; i < oneBlockValues.size(); i++)
  {
    function1VecExpected.at(i) = oneBlockValues.at(i);
    function1VecExpected.at(i + dofsP1) = - oneBlockValues.at(i);
    function1VecExpected.at(i + 2*dofsP1) = oneBlockValues.at(i);
  }
  
  errorIfToDifferent(function1VecValues, function1VecExpected);
#endif 
}


void interpolationTests3D(ParameterDatabase& parmoon_db)
{
  Output::print("\nRunning Interpolation test 3D.\n");
  
  // preparation of grid collection
  
  parmoon_db["boundary_file"].set("Default_UnitCube", false);
  parmoon_db["geo_file"].set("Default_UnitCube_Hexa", false);
  parmoon_db["refinement_n_initial_steps"] = 1;
  
  TDomain domainQ(parmoon_db);
  auto collectionsQ = domainQ.refine_and_get_hierarchy_of_collections(parmoon_db);
  std::list<TCollection*>::iterator itQ = collectionsQ.begin();
  
  parmoon_db["boundary_file"].set("Default_UnitCube", false);
  parmoon_db["geo_file"].set("Default_UnitCube_Tetra", false);
  parmoon_db["refinement_n_initial_steps"] = 1;
  
  TDomain domainP(parmoon_db);
  auto collectionsP = domainP.refine_and_get_hierarchy_of_collections(parmoon_db);
  std::list<TCollection*>::iterator itP = collectionsP.begin();
  
  
  // preparation of spaces
  
  std::shared_ptr<const TFESpace3D> spaceP2{ new TFESpace3D(*itP, std::string("spaceP2"),
                                            BoundConditionNoBoundCondition, 2) };
  
  std::shared_ptr<const TFESpace3D> spaceQ0{ new TFESpace3D(*itQ, std::string("spaceQ0"),
                                         BoundConditionNoBoundCondition, 0) };
  
  std::shared_ptr<const TFESpace3D> spaceP1{ new TFESpace3D(*itP, std::string("spaceP1"),
                                         BoundConditionNoBoundCondition, 1) };
  
  // list of tests
                                         
  scalarFunctionTest3D(spaceP2, spaceQ0, spaceP1);
  vectorFunctionTest3D(spaceP2, spaceQ0, spaceP1);
}
#endif


int main(int argc, char* argv[])
{
  ParameterDatabase parmoon_db = parmoon::parmoon_initialize();
  int my_rank = 0;
#ifdef _MPI
  TDatabase::ParamDB->Comm = MPI_COMM_WORLD;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

  bool testall = false;
  if (argv[1])
  {
    testall = (std::string(argv[1]).compare("testall") == 0);
  }
  parmoon_db.merge(TDomain::default_domain_parameters());
  
  std::vector<TestObject> TestObjects = {prepare_test0(), prepare_test1(),
#ifndef _MPI // for some reason this fails in MPI mode, I have no idea why
                                         prepare_test2(),
#endif
                                         prepare_test3()};
  for(auto &Test : TestObjects)
  {
    if(my_rank == 0)
      Output::print("\n\ntesting geo_file, ", Test.Geo_file);
    parmoon_db["boundary_file"].set(Test.Boundary_file, false);
    parmoon_db["geo_file"].set(Test.Geo_file, false);
    parmoon_db["refinement_n_initial_steps"] = 1;
    
    // Construct domain, thereby read in controls from the input file.
    TDomain domain(parmoon_db);
    domain.refine_and_get_hierarchy_of_collections(parmoon_db);
    
    auto coll = domain.GetCollection(It_Finest, 0);
    // Check P2 and discontinuous P2 elements
    auto max_cases = 1;
    max_cases += (testall) ? 1 : 0;
    for (auto i = 0; i < max_cases; ++i)
    {
      std::shared_ptr<const FESpace> space;
      if (i == 0)
      {
        space = std::shared_ptr<const FESpace>{
          new FESpace(coll, std::string("fe space"),
              BoundConditionNoBoundCondition, 2)};
      }
      else
      {
        space = std::shared_ptr<const FESpace>{
          new FESpace(coll, std::string("fe space"),
              BoundConditionNoBoundCondition, -12)};
      }


      BlockVector data(space->get_n_dof());
      FEFunction function(space, "f", data.get_entries());


      if(my_rank == 0)
        Output::print(" testing with fe space of size ", space->get_n_dof());
      fill_with_arbitrary_data(data);
      compare_with_get_errors(function);

      function.Interpolate(analytic_function);
      auto norms = function.get_L2_and_H1_norm();
      double l2_norm = function.get_L2_norm();
      if(!utilities::are_equal(norms.first, l2_norm))
      {
        ErrThrow("difference in L2 norm ", norms.first - l2_norm);
      }
      if(!utilities::are_equal(norms.first, Test.l2norm))
      {
        ErrThrow("difference in L2 norm ", norms.first - Test.l2norm);
      }
      if(!utilities::are_equal(norms.second, Test.h1norm, tolerance))
      {
        ErrThrow("difference in H1 norm ", norms.second - Test.h1norm);
      }
      double min, max;
      function.MinMax(min, max, false);
      double osc_mean = function.compute_mean_oscillation(0, 1, false);
      double min_nodal_fctn, max_nodal_fctn;
      function.MinMax(min_nodal_fctn, max_nodal_fctn, true);
      double osc_mean_nodal_fctn = function.compute_mean_oscillation(0, 1,
          true);
      if(!utilities::are_equal(min, Test.min) ||
          !utilities::are_equal(max, Test.max))
      {
        ErrThrow("difference in minimum/maximum ", min-Test.min, " ",
            max-Test.max);
      }
      if(!utilities::are_equal(osc_mean, Test.osc_mean))
      {
        ErrThrow("difference in mean oscillation ", osc_mean, Test.osc_mean);
      }
      if(!utilities::are_equal(min_nodal_fctn, Test.min_nodal_fctn) ||
          !utilities::are_equal(max_nodal_fctn, Test.max_nodal_fctn))
      {
        ErrThrow("difference in minimum/maximum ",
            min_nodal_fctn-Test.min_nodal_fctn, " ",
            max_nodal_fctn-Test.max_nodal_fctn);
      }
      if(!utilities::are_equal(osc_mean_nodal_fctn, Test.osc_mean_nodal_fctn))
      {
        ErrThrow("difference in mean oscillation ", osc_mean_nodal_fctn
            - Test.osc_mean_nodal_fctn);
      }
      // Test compute_cell_min_max method by comparing on a single cell
      // The test object stores the minimal and maximal value of the FEFunction
      // on the globally first cell. Determining globally the first cell in MPI
      // case is quite annoying and needs a loop over all cells on each process.
      auto cell_to_check = 0;
      for (auto cell_i = 0; cell_i < coll->GetN_Cells(); ++cell_i)
      {
        auto cell = coll->GetCell(cell_i);
        if (cell->GetGlobalCellNo() == cell_to_check)
        {
#ifdef _MPI
          if (cell->IsHaloCell())
          {
            break;
          }
#endif
          std::vector<double> bf_values;
          BaseFunction_type current_type = BF_C_L_P0_1D; // dummy type
          auto minmax = function.compute_cell_min_max(0, bf_values, current_type, false);
          auto minmax_pk = function.compute_cell_min_max(0, bf_values, current_type, true);
          if(!utilities::are_equal(minmax.first, Test.minmax.first))
          {
            ErrThrow("difference in minimum on first cell ", minmax.first
                - Test.minmax.first);
          }
          if(!utilities::are_equal(minmax.second, Test.minmax.second))
          {
            ErrThrow("difference in maximum on first cell ", minmax.second
                - Test.minmax.second);
          }
          if(!utilities::are_equal(minmax_pk.first, Test.minmax_pk.first))
          {
            ErrThrow("difference in pk minimum on first cell ", minmax_pk.first
                - Test.minmax_pk.first);
          }
          if(!utilities::are_equal(minmax_pk.second, Test.minmax_pk.second))
          {
            ErrThrow("difference in pk maximum on first cell ", minmax_pk.second
                - Test.minmax_pk.second);
          }
          break;
        }
      }
      
      parmoon::Point p = Test.value.first;
      std::vector<double> vals(4, 0.);
      bool found = true; // may become false in mpi on some process
#ifdef __2D__
      function.FindGradient(p.x, p.y, &vals[0]);
#else
      found = function.FindGradient(p.x, p.y, p.z, vals);
#endif
      if(found && !utilities::are_equal(vals[0], Test.value.second))
      {
        ErrThrow("function evaluation at point (", p.x, ",", p.y, ",", p.z,
                 ") is wrong: ", vals[0], " expected: ", Test.value.second);
      }
    }

    delete coll;
    if(my_rank == 0)
      Output::print("finished testing geo_file ", Test.Geo_file);
  }
  
  
#ifdef __2D__
  interpolationTests2D(parmoon_db);
#else
  interpolationTests3D(parmoon_db);
#endif
  
  parmoon::parmoon_finalize();
}
