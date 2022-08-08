#include <Example_NSE3D.h>
#include "NavierStokes.h"

#include "FEDatabase.h"
#include <Database.h>
#include <MainUtilities.h>
#include "BaseCell.h"
#include "BoundFace.h"

#ifdef _MPI
#include <mpi.h>
#endif

#include <array>
#include <cmath>

/* examples */


namespace ansatz_lin_const //0
{
#include "NSE_3D/AnsatzLinConst.h"
}
namespace ansatz_quad_lin //1
{
#include "NSE_3D/AnsatzQuadLin.h"
}
namespace cos_sin_simple //2
{
#include "NSE_3D/CosSin_simple.h"
}
namespace driven_cavity3d //3
{
#include "NSE_3D/DrivenCavity3D.h"
}
namespace flow_around_cylinder_stat
{
#include "NSE_3D/FlowAroundCylinder_stat.h"
}
namespace poiseuille// 5
{
#include "NSE_3D/Poiseuille.h"
}
namespace simple_coriolis // 6
{
#include "NSE_3D/Coriolis_simple.h"
}
namespace brinkman3d_poiseuille // 7
{
  #include "NSE_3D/Brinkman3D_Poiseuille.h"
}
namespace aorta // 8
{
  #include "NSE_3D/Aorta.h"
}
namespace tube // 9
{
  #include "NSE_3D/Tube.h"
}
namespace aortacharite // 10
{
  #include "NSE_3D/AortaCharite.h"
}


//test examples
namespace test_u_0_p_0 //-1
{
#include "NSE_3D/test_u_0_p_0.h"
}
namespace test_u_1_p_0 //-2
{
#include "NSE_3D/test_u_1_p_0.h"
}
namespace test_u_2_p_1 //-3
{
#include "NSE_3D/test_u_2_p_1.h"
}
namespace test_u_3_p_2 //-4
{
#include "NSE_3D/test_u_3_p_2.h"
}


//========================================

Example_NSE3D::Example_NSE3D(const ParameterDatabase& user_input_parameter_db)
 : Example3D(user_input_parameter_db)
{
  int example_code = this->example_database["example"];
  switch( example_code )
  {
    case 0:
    {
      /** exact_solution */
      exact_solution.push_back( ansatz_lin_const::ExactU1 );
      exact_solution.push_back( ansatz_lin_const::ExactU2 );
      exact_solution.push_back( ansatz_lin_const::ExactU3 );
      exact_solution.push_back( ansatz_lin_const::ExactP );

      /* boundary condition */
      boundary_conditions.push_back( ansatz_lin_const::BoundCondition );
      boundary_conditions.push_back( ansatz_lin_const::BoundCondition );
      boundary_conditions.push_back( ansatz_lin_const::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /* boundary values */
      boundary_data.push_back( ansatz_lin_const::U1BoundValue );
      boundary_data.push_back( ansatz_lin_const::U2BoundValue );
      boundary_data.push_back( ansatz_lin_const::U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /* coefficients */
      problem_coefficients = ansatz_lin_const::LinCoeffs;

      /** some variables to change values in the example */
      ansatz_lin_const::DIMENSIONLESS_VISCOSITY = this->get_nu();

      ansatz_lin_const::ExampleFile();
      break;
    }
    case 1:
    {
      /** exact_solution */
      exact_solution.push_back( ansatz_quad_lin::ExactU1 );
      exact_solution.push_back( ansatz_quad_lin::ExactU2 );
      exact_solution.push_back( ansatz_quad_lin::ExactU3 );
      exact_solution.push_back( ansatz_quad_lin::ExactP );

      /* boundary condition */
      boundary_conditions.push_back( ansatz_quad_lin::BoundCondition );
      boundary_conditions.push_back( ansatz_quad_lin::BoundCondition );
      boundary_conditions.push_back( ansatz_quad_lin::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /* boundary values */
      boundary_data.push_back( ansatz_quad_lin::U1BoundValue );
      boundary_data.push_back( ansatz_quad_lin::U2BoundValue );
      boundary_data.push_back( ansatz_quad_lin::U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /* coefficients */
      problem_coefficients = ansatz_quad_lin::LinCoeffs;

      /** some variables to change values in the example */
      ansatz_quad_lin::DIMENSIONLESS_VISCOSITY = this->get_nu();

      ansatz_quad_lin::ExampleFile();
      break;
    }
    case 2:
    {
      /** exact_solution */
      exact_solution.push_back( cos_sin_simple::ExactU1 );
      exact_solution.push_back( cos_sin_simple::ExactU2 );
      exact_solution.push_back( cos_sin_simple::ExactU3 );
      exact_solution.push_back( cos_sin_simple::ExactP );

      /* boundary condition */
      boundary_conditions.push_back( cos_sin_simple::BoundCondition );
      boundary_conditions.push_back( cos_sin_simple::BoundCondition );
      boundary_conditions.push_back( cos_sin_simple::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /* boundary values */
      boundary_data.push_back( cos_sin_simple::U1BoundValue );
      boundary_data.push_back( cos_sin_simple::U2BoundValue );
      boundary_data.push_back( cos_sin_simple::U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /* coefficients */
      problem_coefficients = cos_sin_simple::LinCoeffs;
      
      /** some variables to change values in the example */
      cos_sin_simple::DIMENSIONLESS_VISCOSITY = this->get_nu();

      cos_sin_simple::ExampleFile();
      break;
    }
    case 3:
    {
      /** exact_solution */
      exact_solution.push_back( driven_cavity3d::ExactU1 );
      exact_solution.push_back( driven_cavity3d::ExactU2 );
      exact_solution.push_back( driven_cavity3d::ExactU3 );
      exact_solution.push_back( driven_cavity3d::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( driven_cavity3d::BoundCondition );
      boundary_conditions.push_back( driven_cavity3d::BoundCondition );
      boundary_conditions.push_back( driven_cavity3d::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( driven_cavity3d::U1BoundValue );
      boundary_data.push_back( driven_cavity3d::U2BoundValue );
      boundary_data.push_back( driven_cavity3d::U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = driven_cavity3d::LinCoeffs;
      
      /** some variables to change values in the example */
      driven_cavity3d::DIMENSIONLESS_VISCOSITY = this->get_nu();

      driven_cavity3d::ExampleFile();
      break;
    }
    case 4:
    {
      using namespace flow_around_cylinder_stat;
      /** exact_solution */
      exact_solution.push_back( ExactU1 );
      exact_solution.push_back( ExactU2 );
      exact_solution.push_back( ExactU3 );
      exact_solution.push_back( ExactP );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( U1BoundValue );
      boundary_data.push_back( U2BoundValue );
      boundary_data.push_back( U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = LinCoeffs;

      /**post processing - drag and lift calculation and output */
      post_processing_stat = compute_drag_lift_pdiff;

      /** some variables to change values in the example */
      flow_around_cylinder_stat::DIMENSIONLESS_VISCOSITY = this->get_nu();

      ExampleFile();
      break;
    }
    case 5:
    {
      /** exact_solution */
      exact_solution.push_back( poiseuille::ExactU1 );
      exact_solution.push_back( poiseuille::ExactU2 );
      exact_solution.push_back( poiseuille::ExactU3 );
      exact_solution.push_back( poiseuille::ExactP );
      
      /* boundary condition */
      boundary_conditions.push_back( poiseuille::BoundCondition );
      boundary_conditions.push_back( poiseuille::BoundCondition );
      boundary_conditions.push_back( poiseuille::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /* boundary values */
      boundary_data.push_back( poiseuille::U1BoundValue );
      boundary_data.push_back( poiseuille::U2BoundValue );
      boundary_data.push_back( poiseuille::U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /* coefficients */
      problem_coefficients = poiseuille::LinCoeffs;
      
      /** some variables to change values in the example */
      poiseuille::DIMENSIONLESS_VISCOSITY = this->get_nu();
      
      poiseuille::ExampleFile();
      break;
    }
    case 6:
    {
      /** exact_solution */
      exact_solution.push_back(simple_coriolis::ExactU1 );
      exact_solution.push_back(simple_coriolis::ExactU2 );
      exact_solution.push_back(simple_coriolis::ExactU3 );
      exact_solution.push_back(simple_coriolis::ExactP );
      
      /* boundary condition */
      boundary_conditions.push_back(simple_coriolis::BoundCondition );
      boundary_conditions.push_back(simple_coriolis::BoundCondition );
      boundary_conditions.push_back(simple_coriolis::BoundCondition );
      boundary_conditions.push_back(BoundConditionNoBoundCondition );
      
      /* boundary values */
      boundary_data.push_back(simple_coriolis::U1BoundValue );
      boundary_data.push_back(simple_coriolis::U2BoundValue );
      boundary_data.push_back(simple_coriolis::U3BoundValue );
      boundary_data.push_back(BoundaryValueHomogenous );
      
      /* coefficients */
      problem_coefficients = simple_coriolis::LinCoeffs;
      
      /** some variables to change values in the example */
      simple_coriolis::DIMENSIONLESS_VISCOSITY = this->get_nu();
      simple_coriolis::alpha = 0.05;
      simple_coriolis::v_0 = 1.;
      simple_coriolis::omega = 1.;
      simple_coriolis::include_nonlinear_term = false;
      simple_coriolis::include_coriolis_term = true;
      simple_coriolis::ExampleFile();
      break;
    }
    case 7:
    {
      /** exact_solution */
      exact_solution.push_back( brinkman3d_poiseuille::ExactU1 );
      exact_solution.push_back( brinkman3d_poiseuille::ExactU2 );
      exact_solution.push_back( brinkman3d_poiseuille::ExactU3 );
      exact_solution.push_back( brinkman3d_poiseuille::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( brinkman3d_poiseuille::BoundCondition );
      boundary_conditions.push_back( brinkman3d_poiseuille::BoundCondition );
      boundary_conditions.push_back( brinkman3d_poiseuille::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( brinkman3d_poiseuille::U1BoundValue );
      boundary_data.push_back( brinkman3d_poiseuille::U2BoundValue );
      boundary_data.push_back( brinkman3d_poiseuille::U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = brinkman3d_poiseuille::LinCoeffs;

      // Set dimensionless viscosity
      brinkman3d_poiseuille::effective_viscosity = get_effective_viscosity();
      brinkman3d_poiseuille::sigma = get_inverse_permeability();
      brinkman3d_poiseuille::neumann_id = get_neumann_id();
      brinkman3d_poiseuille::nitsche_id = get_nitsche_id();

      brinkman3d_poiseuille::ExampleFile();
      break;
    }
    case 8:
    {
      /** exact_solution */
      exact_solution.push_back( aorta::ExactU1 );
      exact_solution.push_back( aorta::ExactU2 );
      exact_solution.push_back( aorta::ExactU3 );
      exact_solution.push_back( aorta::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( aorta::BoundCondition );
      boundary_conditions.push_back( aorta::BoundCondition );
      boundary_conditions.push_back( aorta::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( aorta::U1BoundValue );
      boundary_data.push_back( aorta::U2BoundValue );
      boundary_data.push_back( aorta::U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = aorta::LinCoeffs;

      // Set dimensionless viscosity
      aorta::effective_viscosity = get_effective_viscosity();
      aorta::neumann_id = get_neumann_id();
      aorta::nitsche_id = get_nitsche_id();

      aorta::ExampleFile();
      break;
    }
    case 9:
    {
      /** exact_solution */
      exact_solution.push_back( tube::ExactU1 );
      exact_solution.push_back( tube::ExactU2 );
      exact_solution.push_back( tube::ExactU3 );
      exact_solution.push_back( tube::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( tube::BoundCondition );
      boundary_conditions.push_back( tube::BoundCondition );
      boundary_conditions.push_back( tube::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( tube::U1BoundValue );
      boundary_data.push_back( tube::U2BoundValue );
      boundary_data.push_back( tube::U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = tube::LinCoeffs;

      // Set dimensionless viscosity
      tube::effective_viscosity = get_effective_viscosity();
      tube::neumann_id = get_neumann_id();
      tube::nitsche_id = get_nitsche_id();
      tube::windkessel_id = get_windkessel_id();

      tube::ExampleFile();
      break;
    }
    case 10:
    {
      /** exact_solution */
      exact_solution.push_back( aortacharite::ExactU1 );
      exact_solution.push_back( aortacharite::ExactU2 );
      exact_solution.push_back( aortacharite::ExactU3 );
      exact_solution.push_back( aortacharite::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( aortacharite::BoundCondition );
      boundary_conditions.push_back( aortacharite::BoundCondition );
      boundary_conditions.push_back( aortacharite::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( aortacharite::U1BoundValue );
      boundary_data.push_back( aortacharite::U2BoundValue );
      boundary_data.push_back( aortacharite::U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = aortacharite::LinCoeffs;

      // Set dimensionless viscosity
      aortacharite::effective_viscosity = get_effective_viscosity();
      aortacharite::neumann_id = get_neumann_id();
      aortacharite::nitsche_id = get_nitsche_id();
      aortacharite::windkessel_id = get_windkessel_id();

      aortacharite::ExampleFile();
      break;
    }
    case -1:
    {
      using namespace test_u_0_p_0;
      /** exact_solution */
      exact_solution.push_back( ExactU1 );
      exact_solution.push_back( ExactU2 );
      exact_solution.push_back( ExactU3 );
      exact_solution.push_back( ExactP );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( U1BoundValue );
      boundary_data.push_back( U2BoundValue );
      boundary_data.push_back( U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = LinCoeffs;

      /** some variables to change values in the example */
      test_u_0_p_0::DIMENSIONLESS_VISCOSITY = this->get_nu();

      ExampleFile();
      break;
    }
    case -2:
    {
      using namespace test_u_1_p_0;
      /** exact_solution */
      exact_solution.push_back( ExactU1 );
      exact_solution.push_back( ExactU2 );
      exact_solution.push_back( ExactU3 );
      exact_solution.push_back( ExactP );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( U1BoundValue );
      boundary_data.push_back( U2BoundValue );
      boundary_data.push_back( U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = LinCoeffs;

      /** some variables to change values in the example */
      test_u_1_p_0::DIMENSIONLESS_VISCOSITY = this->get_nu();

      ExampleFile();
      break;
    }
    case -3:
    {
      using namespace test_u_2_p_1;
      /** exact_solution */
      exact_solution.push_back( ExactU1 );
      exact_solution.push_back( ExactU2 );
      exact_solution.push_back( ExactU3 );
      exact_solution.push_back( ExactP );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( U1BoundValue );
      boundary_data.push_back( U2BoundValue );
      boundary_data.push_back( U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = LinCoeffs;

      /** some variables to change values in the example */
      test_u_2_p_1::DIMENSIONLESS_VISCOSITY = this->get_nu();

      ExampleFile();
      break;
    }
    case -4:
    {
      using namespace test_u_3_p_2;
      /** exact_solution */
      exact_solution.push_back( ExactU1 );
      exact_solution.push_back( ExactU2 );
      exact_solution.push_back( ExactU3 );
      exact_solution.push_back( ExactP );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( U1BoundValue );
      boundary_data.push_back( U2BoundValue );
      boundary_data.push_back( U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = LinCoeffs;

      /** some variables to change values in the example */
      test_u_3_p_2::DIMENSIONLESS_VISCOSITY = this->get_nu();

      ExampleFile();
      break;
    }
  }
}

void Example_NSE3D::do_post_processing(NavierStokes<3>& nse3d) const
{
  if(post_processing_stat)
  {
    post_processing_stat(nse3d);
  }
  else
  {
#ifdef _MPI
	  int my_rank;
	  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	  if (my_rank == 0)
#endif
	    Output::info<2>("Example_NSE3D","No post processing done for the current example.");
  }
}

double Example_NSE3D::get_nu() const
{
  double inverse_reynolds = this->example_database["reynolds_number"];
  inverse_reynolds = 1/inverse_reynolds;
  return inverse_reynolds;
}

double Example_NSE3D::get_effective_viscosity() const
{
double effective_viscosity = this->example_database["effective_viscosity"];
return effective_viscosity;
}

double Example_NSE3D::get_inverse_permeability() const
{
  return this->example_database["inverse_permeability"];
}

std::vector<size_t> Example_NSE3D::get_neumann_id() const
{
  std::vector<size_t> neumann_id;
  int n_neumann_bd = this->example_database["n_neumann_bd"];
  if (n_neumann_bd)
    return this->example_database["neumann_id"];
  else
    return neumann_id;
}


std::vector<size_t> Example_NSE3D::get_nitsche_id() const
{
  std::vector<size_t> nitsche_id;
  int n_nitsche_bd = this->example_database["n_nitsche_bd"];
  if (n_nitsche_bd) 
    return this->example_database["nitsche_id"];
  else
    return nitsche_id;
}

std::vector<size_t> Example_NSE3D::get_windkessel_id() const
{
  std::vector<size_t> windkessel_id;
  int n_windkessel_bd = this->example_database["n_windkessel_bd"];
  if (n_windkessel_bd) 
    return this->example_database["windkessel_id"];
  else
    return windkessel_id;
}
