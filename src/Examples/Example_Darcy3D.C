#include <Example_Darcy3D.h>

#include <Database.h>
#include <MainUtilities.h>
#include <algorithm>
#include <cmath>


/* examples */

namespace sine_simple
{
  #include "Darcy_3D/SinCos.h"  
}
namespace simple_benchmark
{
  #include "Darcy_3D/Benchmark.h"
}
namespace nine_spot
{
  #include "Darcy_3D/9SpotProblem.h"
}


Example_Darcy3D::Example_Darcy3D(
  const ParameterDatabase& user_input_parameter_db)
 : Example3D(user_input_parameter_db)
{
  int example_code = this->example_database["example"];
  switch( example_code )
  {
    case 0:
      /** exact_solution */
      exact_solution.push_back( sine_simple::ExactU1 );
      exact_solution.push_back( sine_simple::ExactU2 );
      exact_solution.push_back( sine_simple::ExactU3 );
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
      exact_solution.push_back( simple_benchmark::ExactU3 );
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
      exact_solution.push_back( nine_spot::ExactU1 );
      exact_solution.push_back( nine_spot::ExactU2 );
      exact_solution.push_back( nine_spot::ExactU3 );
      exact_solution.push_back( nine_spot::ExactP );
      
      /** boundary condition */
      boundary_conditions.push_back( nine_spot::BoundCondition );
      boundary_conditions.push_back( BoundConditionNoBoundCondition );
      
      /** boundary values */
      boundary_data.push_back( nine_spot::FluxBoundValue );
      boundary_data.push_back( BoundaryValueHomogenous );
      
      /** coefficients */
      problem_coefficients = nine_spot::LinCoeffs;
      nine_spot::factor = this->example_database["darcy_permeability_jump"];
      if(user_input_parameter_db.contains("refinement_n_initial_steps"))
      {
        // assume we use a uniformly refined mesh on the unit cube:
        size_t nr = user_input_parameter_db["refinement_n_initial_steps"];
        nine_spot::mesh_width_at_corners = 1./std::pow(2, nr);
      }
      
      nine_spot::ExampleFile();
      break;
    default:
      ErrThrow("Unknown name of the mixed Darcy example!");
  }
}

void Example_Darcy3D::do_post_processing(Darcy<3>& darcy2d) const
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
      Output::info<2>("Example_Darcy3D","No post processing done for the current example.");
  }
}

double Example_Darcy3D::get_nu() const
{
  double inverse_reynolds = this->example_database["reynolds_number"];
  inverse_reynolds = 1/inverse_reynolds;
  return inverse_reynolds;
}

