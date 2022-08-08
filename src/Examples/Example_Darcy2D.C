#include <Example_Darcy2D.h>

#include <Database.h>
#include <MainUtilities.h>
#include <cmath>

/* examples */
namespace sine_simple
{
  #include "Darcy_2D/SinSinSolution.h"  
}
namespace simple_benchmark
{
  #include "Darcy_2D/Benchmark.h"
}
namespace cubic_pressure
{
  #include "Darcy_2D/CubicPressure.h"
}
namespace obstacle
{
  #include "Darcy_2D/Obstacle.h"
}
namespace five_spot
{
  #include "Darcy_2D/5SpotProblem.h"
}
namespace polynomial_darcy
{
  #include "Darcy_2D/Polynomial.h"
}


Example_Darcy2D::Example_Darcy2D(
  const ParameterDatabase& user_input_parameter_db)
 : Example2D(user_input_parameter_db)
{
  int example_code = this->example_database["example"];
  switch( example_code )
  {
    case 0:
      /** exact_solution */
      exact_solution.push_back( sine_simple::ExactU1 );
      exact_solution.push_back( sine_simple::ExactU2 );
      exact_solution.push_back( sine_simple::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( sine_simple::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( sine_simple::FluxBoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = sine_simple::LinCoeffs;
      
      sine_simple::ExampleFile();
      break;
    case 1:
      /** exact_solution */
      exact_solution.push_back( simple_benchmark::ExactU1 );
      exact_solution.push_back( simple_benchmark::ExactU2 );
      exact_solution.push_back( simple_benchmark::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( simple_benchmark::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( simple_benchmark::UNBoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = simple_benchmark::LinCoeffs;
      
      simple_benchmark::ExampleFile();
      break;
    case 2:
      /** exact_solution */
      exact_solution.push_back( cubic_pressure::ExactU1 );
      exact_solution.push_back( cubic_pressure::ExactU2 );
      exact_solution.push_back( cubic_pressure::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( cubic_pressure::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( cubic_pressure::FluxBoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = cubic_pressure::LinCoeffs;
      
      cubic_pressure::ExampleFile();
      break;
    case 3:
      /** exact_solution */
      exact_solution.push_back( obstacle::ExactU1 );
      exact_solution.push_back( obstacle::ExactU2 );
      exact_solution.push_back( obstacle::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( obstacle::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( obstacle::FluxBoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = obstacle::LinCoeffs;
      
      obstacle::factor = this->example_database["darcy_permeability_jump"];
      
      obstacle::ExampleFile();
      break;
    case 4:
      /** exact_solution */
      exact_solution.push_back( five_spot::ExactU1 );
      exact_solution.push_back( five_spot::ExactU2 );
      exact_solution.push_back( five_spot::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( five_spot::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( five_spot::FluxBoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = five_spot::LinCoeffs;
      
      five_spot::factor = this->example_database["darcy_permeability_jump"];
      // valid values for VELOCITY_SPACE are: 1000,1001,1002,1003,1011,1012,1013
      five_spot::polynomial_order = TDatabase::ParamDB->VELOCITY_SPACE % 10; 
      if(user_input_parameter_db.contains("refinement_n_initial_steps"))
      {
        // assume we use a uniformly refined mesh on the unit square:
        size_t nr = user_input_parameter_db["refinement_n_initial_steps"];
        five_spot::mesh_width_at_corners = 1./std::pow(2, nr);
      }
      
      five_spot::ExampleFile();
      break;
    case 5:
      /** exact_solution */
      exact_solution.push_back( polynomial_darcy::ExactU1 );
      exact_solution.push_back( polynomial_darcy::ExactU2 );
      exact_solution.push_back( polynomial_darcy::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( polynomial_darcy::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( polynomial_darcy::UNormalValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = polynomial_darcy::LinCoeffsDarcy;
      
      polynomial_darcy::deg = this->example_database["degree_polynomial"];
      polynomial_darcy::is_RT = this->example_database["is_RT_polynomial"];
      
      polynomial_darcy::ExampleFileDarcy();
      break;
    default:
      ErrThrow("Unknown name of the mixed Darcy example!");
      break;
  }
}

void Example_Darcy2D::do_post_processing(Darcy<2>& darcy2d) const
{
  if(post_processing_stat)
  {
    post_processing_stat(darcy2d);
  }
  else
  {
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0)
#endif
      Output::info<2>("Example_Darcy2D","No post processing done for the current example.");
  }
}

double Example_Darcy2D::get_nu() const
{
  double inverse_reynolds = this->example_database["reynolds_number"];
  inverse_reynolds = 1/inverse_reynolds;
  return inverse_reynolds;
}

