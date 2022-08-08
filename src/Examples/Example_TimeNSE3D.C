#include <Example_TimeNSE3D.h>
#include "TimeNavierStokes.h"
#include "FEDatabase.h"
#include <Database.h>
#include <MainUtilities.h>
#include "BaseCell.h"
#include "BoundFace.h"
#include <array>
#include <cmath>

namespace lin_space_time
{
  #include "TNSE_3D/linear_space_time.h"  // 0
}
namespace AnsatzLinConst
{
  #include "TNSE_3D/AnsatzLinConst.h"     // 1
}
namespace Bsp0
{
  #include "TNSE_3D/Bsp0.h"   // 2
}
namespace Bsp1
{
  #include "TNSE_3D/Bsp1.h"   // 3
}
namespace Bsp2
{
 #include "TNSE_3D/Bsp2.h"    // 4
}
namespace Bsp3
{
  #include "TNSE_3D/Bsp3.h"   // 5
}

namespace flow_around_cylinder_instationary
{
#include "TNSE_3D/FlowAroundCylinder_instat.h"   // 6
}

namespace cylinder
{
#include "TNSE_3D/Cylinder.h"   // 7
}

namespace dairybuilding
{
#include "TNSE_3D/Dairybuilding.h"   // 8
}

namespace tube_instationary
{
#include "TNSE_3D/Tube.h"   // 9
}

namespace aortacharite_time
{
  #include "TNSE_3D/AortaChariteTime.h" // 10
}

#include "TNSE_3D/ChannelTau.h" // 11

namespace driven_cavity
{
#include "TNSE_3D/DrivenCavity.h"  
}


//=========================================================
Example_TimeNSE3D::Example_TimeNSE3D(
  const ParameterDatabase& user_input_parameter_db)
 : Example_NonStationary3D(user_input_parameter_db)
{
  int example_code = this->example_database["example"];
  switch(example_code)
  {
    case 0:
    {
      using namespace lin_space_time;
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

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );
      initialCondtion.push_back( InitialP );

      /** some variables to change values in the example */
      lin_space_time::DIMENSIONLESS_VISCOSITY = this->get_nu();

      ExampleFile();
      break;
    }
    case 1:
    {
      using namespace AnsatzLinConst;
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

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );
      initialCondtion.push_back( InitialP );

      /** some variables to change values in the example */
      AnsatzLinConst::DIMENSIONLESS_VISCOSITY = this->get_nu();

      ExampleFile();
      break;
    }
    case 2:
    {
      using namespace Bsp0;
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

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );
      initialCondtion.push_back( InitialP );

      /** some variables to change values in the example */
      Bsp0::DIMENSIONLESS_VISCOSITY = this->get_nu();

      ExampleFile();
      break;
    }
    case 3:
    {
      using namespace Bsp1;
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

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );
      initialCondtion.push_back( InitialP );

      /** some variables to change values in the example */
      Bsp1::DIMENSIONLESS_VISCOSITY = this->get_nu();

      ExampleFile();
      break;
    }    
    case 4:
    {
      using namespace Bsp2;
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

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );
      initialCondtion.push_back( InitialP );

      /** some variables to change values in the example */
      Bsp2::DIMENSIONLESS_VISCOSITY = this->get_nu();

      ExampleFile();
      break;
    }
    case 5:
    {
      using namespace Bsp3;
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

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );
      initialCondtion.push_back( InitialP );

      /** some variables to change values in the example */
      Bsp3::DIMENSIONLESS_VISCOSITY = this->get_nu();

      ExampleFile();
      break;
    }
    case 6:
    {
      using namespace flow_around_cylinder_instationary;
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

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );
      initialCondtion.push_back( InitialP );

      /**post processing - drag and lift calculation and output */
      post_processing_stat = flow_around_cylinder_instationary::compute_drag_lift_pdiff;

      /** some variables to change values in the example */
      flow_around_cylinder_instationary::DIMENSIONLESS_VISCOSITY = this->get_nu();

      ExampleFile();
      break;
    }
    case 7:
    {
      using namespace cylinder;
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

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );
      initialCondtion.push_back( InitialP );

 
      /** some variables to change values in the example */
      cylinder::DIMENSIONLESS_VISCOSITY = this->get_nu();

      ExampleFile();
      break;
    }
    case 8:
    {
      using namespace dairybuilding;
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

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );
      initialCondtion.push_back( InitialP );

 
      /** some variables to change values in the example */
      dairybuilding::DIMENSIONLESS_VISCOSITY = this->get_nu();

      ExampleFile();
      break;
    }
    case 9:
    {
      using namespace tube_instationary;
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

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );
      initialCondtion.push_back( InitialP );

 
      /** some variables to change values in the example */
      tube_instationary::effective_viscosity = this->get_nu();
      tube_instationary::neumann_id = this->get_neumann_id();
      tube_instationary::windkessel_id = this->get_windkessel_id();

      ExampleFile();
      break;
    }
    case 10:
    {
      using namespace aortacharite_time;
      /** exact_solution */
      exact_solution.push_back( aortacharite_time::ExactU1 );
      exact_solution.push_back( aortacharite_time::ExactU2 );
      exact_solution.push_back( aortacharite_time::ExactU3 );
      exact_solution.push_back( aortacharite_time::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( aortacharite_time::BoundCondition );
      boundary_conditions.push_back( aortacharite_time::BoundCondition );
      boundary_conditions.push_back( aortacharite_time::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( aortacharite_time::U1BoundValue );
      boundary_data.push_back( aortacharite_time::U2BoundValue );
      boundary_data.push_back( aortacharite_time::U3BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = aortacharite_time::LinCoeffs;

      /** initial conditions */
      initialCondtion.push_back( aortacharite_time::InitialU1 );
      initialCondtion.push_back( aortacharite_time::InitialU2 );
      initialCondtion.push_back( aortacharite_time::InitialU3 );
      initialCondtion.push_back( aortacharite_time::InitialP );

      /** some variables to change values in the example */
      aortacharite_time::effective_viscosity = get_effective_viscosity();//this->get_nu();
      aortacharite_time::neumann_id = this->get_neumann_id();
      aortacharite_time::windkessel_id = this->get_windkessel_id();

      ExampleFile();
      break;
    }
    case 11:
    {
      using namespace channel_tau;
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

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );
      initialCondtion.push_back( InitialP );
 
      /** some variables to change values in the example */
      channel_tau::RE = this->get_reynolds();
      ExampleFile(user_input_parameter_db);
      break;
    }
    case 12:
    {
      using namespace driven_cavity;
      exact_solution.push_back(ExactU1);
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

      /** initial conditions */
      initialCondtion.push_back( InitialU1 );
      initialCondtion.push_back( InitialU2 );
      initialCondtion.push_back( InitialU3 );
      initialCondtion.push_back( InitialP );
 
      /** some variables to change values in the example */
      driven_cavity::DIMENSIONLESS_VISCOSITY = this->get_nu();
      ExampleFile();
      break;
    }
  default:
      ErrThrow("Unknown Example_TimeNSE3D example!");
  }
}

void Example_TimeNSE3D::do_post_processing(TimeNavierStokes<3>& tnse3d, double&) const
{
  if(post_processing_stat)
  {
    post_processing_stat(tnse3d);
  }
}

double Example_TimeNSE3D::get_reynolds() const
{
  double reynolds = this->example_database["reynolds_number"];
  return reynolds;
}

double Example_TimeNSE3D::get_nu() const
{
  double inverse_reynolds = this->example_database["reynolds_number"];
  inverse_reynolds = 1/inverse_reynolds;
  return inverse_reynolds;
}

double Example_TimeNSE3D::get_effective_viscosity() const
{
  double effective_viscosity = this->example_database["effective_viscosity"];
  return effective_viscosity;
}


std::vector<size_t> Example_TimeNSE3D::get_neumann_id() const
{
  std::vector<size_t> neumann_id;
  int n_neumann_bd = this->example_database["n_neumann_bd"];
  if (n_neumann_bd)
    return this->example_database["neumann_id"];
  else
    return neumann_id;
}


std::vector<size_t> Example_TimeNSE3D::get_nitsche_id() const
{
  std::vector<size_t> nitsche_id;
  int n_nitsche_bd = this->example_database["n_nitsche_bd"];
  if (n_nitsche_bd) 
    return this->example_database["nitsche_id"];
  else
    return nitsche_id;
}

std::vector<size_t> Example_TimeNSE3D::get_windkessel_id() const
{
  std::vector<size_t> ids;
  int n = this->example_database.try_get_value("n_windkessel_bd", 0);
  if (n) 
    return this->example_database["windkessel_id"];
  else
    return ids;
}

