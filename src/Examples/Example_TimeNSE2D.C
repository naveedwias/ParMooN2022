#include <Example_TimeNSE2D.h>
#include <Database.h>
#include <MainUtilities.h>
#include <BoundEdge.h>
#include "FEDatabase.h"
#include "TimeNavierStokes.h"
#include "AuxParam2D.h" // used in MixingLayerSlipSmallSquares.h
#include "BaseCell.h"
#include <string>

namespace bsp1              // case 0
{
 #include "TNSE_2D/Bsp1.h"
}

namespace lin_space_time   // case 1
{
#include "TNSE_2D/linear_space_time.h"
}

namespace sincosexp        // case 2
{
#include "TNSE_2D/SinCosExp.h"
}

namespace flow_around_cylinder_steady_inflow     // case 3
{
#include "flow_around_cylinder_steady_inflow.h"
}

namespace backward_facing_step_time  // case 4
{
#include "TNSE_2D/backward_facing_step.h"
}

namespace driven_cavity_time         // case 5
{
#include "TNSE_2D/DrivenCavity.h"
}

namespace mixing_layer_us       // case 6
{
#include "TNSE_2D/MixingLayerSlipSmallSquares.h"
}

namespace flow_around_cylinder_transient_inflow
{
#include "TNSE_2D/flow_around_cylinder_transient_inflow.h"
}

namespace two_outlets
{
#include "TNSE_2D/Two_Outlets.h"
}

namespace poly_sin
{
#include "TNSE_2D/PolySin.h"
}

namespace sinsincospi
{
#include "TNSE_2D/SinSinCosPi.h"
}

Example_TimeNSE2D::Example_TimeNSE2D(
  const ParameterDatabase& user_input_parameter_db)
 : Example_NonStationary2D(user_input_parameter_db)
{
  int example_code = this->example_database["example"];
  switch(example_code)
  {
    case 0:
      /** exact_solution */
      exact_solution.push_back( bsp1::ExactU1 );
      exact_solution.push_back( bsp1::ExactU2 );
      exact_solution.push_back( bsp1::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( bsp1::BoundCondition );
      boundary_conditions.push_back( bsp1::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( bsp1::U1BoundValue );
      boundary_data.push_back( bsp1::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = bsp1::LinCoeffs;
      
      /** initial condition */
      initialCondition.push_back(bsp1::InitialU1);
      initialCondition.push_back(bsp1::InitialU2);
      initialCondition.push_back(bsp1::InitialP);
      bsp1::ExampleFile();
      break;
    case 1:
      /** exact_solution */
      exact_solution.push_back( lin_space_time::ExactU1 );
      exact_solution.push_back( lin_space_time::ExactU2 );
      exact_solution.push_back( lin_space_time::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( lin_space_time::BoundCondition );
      boundary_conditions.push_back( lin_space_time::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( lin_space_time::U1BoundValue );
      boundary_data.push_back( lin_space_time::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = lin_space_time::LinCoeffs;
      
      initialCondition.push_back(lin_space_time::InitialU1);
      initialCondition.push_back(lin_space_time::InitialU2);
      initialCondition.push_back(lin_space_time::InitialP);
      
      // Set dimensionless viscosity
      lin_space_time::DIMENSIONLESS_VISCOSITY = get_nu();
      
      lin_space_time::ExampleFile();
      break;
      
    case 2: // SinCosExp
      /** exact_solution */
      exact_solution.push_back(sincosexp::ExactU1 );
      exact_solution.push_back(sincosexp::ExactU2 );
      exact_solution.push_back(sincosexp::ExactP );

      /** boundary condition */
      boundary_conditions.push_back(sincosexp::BoundCondition );
      boundary_conditions.push_back(sincosexp::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back(sincosexp::U1BoundValue );
      boundary_data.push_back(sincosexp::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      // Set dimensionless viscosity
      sincosexp::DIMENSIONLESS_VISCOSITY = get_nu();

      /** coefficients */
      problem_coefficients =sincosexp::LinCoeffs;

      initialCondition.push_back(sincosexp::InitialU1);
      initialCondition.push_back(sincosexp::InitialU2);
      initialCondition.push_back(sincosexp::InitialP);

      sincosexp::ExampleFile();
      break;

    case 3:
      /** exact_solution */
      exact_solution.push_back( flow_around_cylinder_steady_inflow::ExactU1 );
      exact_solution.push_back( flow_around_cylinder_steady_inflow::ExactU2 );
      exact_solution.push_back( flow_around_cylinder_steady_inflow::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( flow_around_cylinder_steady_inflow::BoundCondition );
      boundary_conditions.push_back( flow_around_cylinder_steady_inflow::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( flow_around_cylinder_steady_inflow::U1BoundValue );
      boundary_data.push_back( flow_around_cylinder_steady_inflow::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = flow_around_cylinder_steady_inflow::LinCoeffs;
      
      initialCondition.push_back(flow_around_cylinder_steady_inflow::InitialU1);
      initialCondition.push_back(flow_around_cylinder_steady_inflow::InitialU2);
      initialCondition.push_back(flow_around_cylinder_steady_inflow::InitialP);
      
      // Set dimensionless viscosity
      flow_around_cylinder_steady_inflow::DIMENSIONLESS_VISCOSITY = get_nu();

      /**post processing - drag and lift calculation and output */
      post_processing_stat = flow_around_cylinder_steady_inflow::compute_drag_lift_pdiff;

      flow_around_cylinder_steady_inflow::ExampleFile();
      break;
    case 4:
      exact_solution.push_back( backward_facing_step_time::ExactU1 );
      exact_solution.push_back( backward_facing_step_time::ExactU2 );
      exact_solution.push_back( backward_facing_step_time::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( backward_facing_step_time::BoundCondition );
      boundary_conditions.push_back( backward_facing_step_time::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( backward_facing_step_time::U1BoundValue );
      boundary_data.push_back( backward_facing_step_time::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = backward_facing_step_time::LinCoeffs;
      
      initialCondition.push_back(backward_facing_step_time::InitialU1);
      initialCondition.push_back(backward_facing_step_time::InitialU2);
      initialCondition.push_back(backward_facing_step_time::InitialP);
      
      backward_facing_step_time::ExampleFile();
      backward_facing_step_time::DIMENSIONLESS_VISCOSITY = this->get_nu();
      break;
    case 5:
      exact_solution.push_back( driven_cavity_time::ExactU1 );
      exact_solution.push_back( driven_cavity_time::ExactU2 );
      exact_solution.push_back( driven_cavity_time::ExactP );

      /** boundary condition */
      boundary_conditions.push_back( driven_cavity_time::BoundCondition );
      boundary_conditions.push_back( driven_cavity_time::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );

      /** boundary values */
      boundary_data.push_back( driven_cavity_time::U1BoundValue );
      boundary_data.push_back( driven_cavity_time::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );

      /** coefficients */
      problem_coefficients = driven_cavity_time::LinCoeffs;

      initialCondition.push_back(driven_cavity_time::InitialU1);
      initialCondition.push_back(driven_cavity_time::InitialU2);
      initialCondition.push_back(driven_cavity_time::InitialP);

      driven_cavity_time::DIMENSIONLESS_VISCOSITY = this->get_nu();
      driven_cavity_time::ExampleFile();
      break;
    case 6:
      exact_solution.push_back( mixing_layer_us::ExactU1 );
      exact_solution.push_back( mixing_layer_us::ExactU2 );
      exact_solution.push_back( mixing_layer_us::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( mixing_layer_us::BoundCondition );
      boundary_conditions.push_back( mixing_layer_us::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( mixing_layer_us::U1BoundValue );
      boundary_data.push_back( mixing_layer_us::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = mixing_layer_us::LinCoeffs;
      
      initialCondition.push_back(mixing_layer_us::InitialU1);
      initialCondition.push_back(mixing_layer_us::InitialU2);
      initialCondition.push_back(mixing_layer_us::InitialP);
      
      mixing_layer_us::ExampleFile();
      
      mixing_layer_us::DIMENSIONLESS_VISCOSITY = this->get_nu();
      
      /**post processing - drag and lift calculation and output */
      post_processing_stat = mixing_layer_us::EvaluateSolution;
      break;
    case 7:
      /** exact_solution */
      exact_solution.push_back( flow_around_cylinder_transient_inflow::ExactU1 );
      exact_solution.push_back( flow_around_cylinder_transient_inflow::ExactU2 );
      exact_solution.push_back( flow_around_cylinder_transient_inflow::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( flow_around_cylinder_transient_inflow::BoundCondition );
      boundary_conditions.push_back( flow_around_cylinder_transient_inflow::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( flow_around_cylinder_transient_inflow::U1BoundValue );
      boundary_data.push_back( flow_around_cylinder_transient_inflow::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = flow_around_cylinder_transient_inflow::LinCoeffs;
      
      initialCondition.push_back(flow_around_cylinder_transient_inflow::InitialU1);
      initialCondition.push_back(flow_around_cylinder_transient_inflow::InitialU2);
      initialCondition.push_back(flow_around_cylinder_transient_inflow::InitialP);
      
      // Set dimensionless viscosity
      flow_around_cylinder_transient_inflow::DIMENSIONLESS_VISCOSITY = get_nu();

      /**post processing - drag and lift calculation and output */
      post_processing_stat = flow_around_cylinder_transient_inflow::compute_drag_lift_pdiff;

      flow_around_cylinder_transient_inflow::ExampleFile();
      break;
    case 8:
      /** exact_solution */
      exact_solution.push_back( poly_sin::ExactU1 );
      exact_solution.push_back( poly_sin::ExactU2 );
      exact_solution.push_back( poly_sin::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( poly_sin::BoundCondition );
      boundary_conditions.push_back( poly_sin::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( poly_sin::U1BoundValue );
      boundary_data.push_back( poly_sin::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = poly_sin::LinCoeffs;
      
      initialCondition.push_back(poly_sin::InitialU1);
      initialCondition.push_back(poly_sin::InitialU2);
      initialCondition.push_back(poly_sin::InitialP);
      
      // Set dimensionless viscosity
      poly_sin::DIMENSIONLESS_VISCOSITY = get_nu();


      poly_sin::ExampleFile();
      break;
    case 9:
      /** exact_solution */
      exact_solution.push_back( sinsincospi::ExactU1 );
      exact_solution.push_back( sinsincospi::ExactU2 );
      exact_solution.push_back( sinsincospi::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( sinsincospi::BoundCondition );
      boundary_conditions.push_back( sinsincospi::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( sinsincospi::U1BoundValue );
      boundary_data.push_back( sinsincospi::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = sinsincospi::LinCoeffs;
      
      initialCondition.push_back(sinsincospi::InitialU1);
      initialCondition.push_back(sinsincospi::InitialU2);
      initialCondition.push_back(sinsincospi::InitialP);
      
      // Set dimensionless viscosity
      sinsincospi::DIMENSIONLESS_VISCOSITY = get_nu();

      sinsincospi::ExampleFile();
      break;
      /* we use example 14 to be consistent with the steady NSE case */
  case 14:
      /** exact_solution */
      exact_solution.push_back( two_outlets::ExactU1 );
      exact_solution.push_back( two_outlets::ExactU2 );
      exact_solution.push_back( two_outlets::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( two_outlets::BoundCondition );
      boundary_conditions.push_back( two_outlets::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( two_outlets::U1BoundValue );
      boundary_data.push_back( two_outlets::U2BoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = two_outlets::LinCoeffs;
      
      initialCondition.push_back(two_outlets::InitialU1);
      initialCondition.push_back(two_outlets::InitialU2);
      initialCondition.push_back(two_outlets::InitialP);
      
      // Set dimensionless viscosity
      two_outlets::effective_viscosity = get_nu();

      // boundary conditions
      two_outlets::neumann_id = get_neumann_id();
      two_outlets::nitsche_id = get_nitsche_id();
      two_outlets::windkessel_id = get_windkessel_id();
     

      two_outlets::ExampleFile();
      break;
  default:
      ErrThrow("Unknown time-dependent Example_TimeNSE2D example!");
  }
}

Example_TimeNSE2D::Example_TimeNSE2D(
  const std::vector<DoubleFunct2D*>& exact,
  const std::vector<BoundCondFunct2D*>& bc,
  const std::vector<BoundValueFunct2D*>& bd, const CoeffFct2D& coeffs,
  bool timedependentrhs, bool timedependentcoeffs,
  const std::vector<DoubleFunct2D*>& init_cond)
  : Example_NonStationary2D(exact, bc, bd, coeffs, timedependentrhs,
                            timedependentcoeffs, init_cond)
  {

  }

void Example_TimeNSE2D::do_post_processing(TimeNavierStokes<2>& tnse2d,
                                           double& val) const
{
  if(post_processing_stat)
  {
    post_processing_stat(tnse2d, val);
  }
  else
  {
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0)
#endif
      Output::info<2>("Example_TimeNSE2D","No post processing done for the current example.");
  }
}
void Example_TimeNSE2D::do_post_processing(TimeNavierStokes<2>& tnse2d) const
{
  if(post_processing_stat_old)
  {
    post_processing_stat_old(tnse2d);
  }
  else
  {
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0)
#endif
      Output::info<2>("Example_TimeNSE2D","No post processing done for the current example.");
  }
}

double Example_TimeNSE2D::get_nu() const
{
  int example_code = this->example_database["example"];
  double inverse_reynolds;
  double nu=-1.;

  if (example_code == 8) {
    nu = this->example_database["effective_viscosity"];
  } else {
    inverse_reynolds = this->example_database["reynolds_number"];
    nu = 1/inverse_reynolds;
  }
  return nu;
}

double Example_TimeNSE2D::get_inverse_permeability() const
{
  return this->example_database["inverse_permeability"];
}


std::vector<size_t> Example_TimeNSE2D::get_neumann_id() const
{
  std::vector<size_t> empty_vector;
  int n_neumann_bd = this->example_database["n_neumann_bd"];
  if (n_neumann_bd){
    return this->example_database["neumann_id"];
  } else {
    return empty_vector;
  }
}


std::vector<size_t> Example_TimeNSE2D::get_nitsche_id() const
{
  std::vector<size_t> empty_vector;
  int n_nitsche_bd = this->example_database["n_nitsche_bd"];
  if (n_nitsche_bd){
    return this->example_database["nitsche_id"];
  } else {
    return empty_vector;
  }
}

std::vector<size_t> Example_TimeNSE2D::get_windkessel_id() const
{
  std::vector<size_t> empty_vector;
  int n_windkessel_bd = this->example_database["n_windkessel_bd"];
  if (n_windkessel_bd){
    return this->example_database["windkessel_id"];
  } else {
    return empty_vector;
  }
}
