#include <Example_TimeCD3D.h>
#include <Database.h>
#include "MainUtilities.h"
#include <cmath>
#include <string.h>

namespace linear_space_time
{
  #include "TCD_3D/linear_space_time.h"
}

namespace quad_space_time
{
#include "TCD_3D/quadratic_space_time.h"
}

namespace concentration
{
#include "TCD_3D/concentrationOfSpecies_3d.h"
}

namespace HeatChannel
{
#include "TCD_3D/HeatChanel.h"
}

namespace block_xyzt_dirichlet
{
#include "TCD_3D/block_xyzt_dirichlet.h"
}

namespace block_xyzt_neumann
{
#include "TCD_3D/block_xyzt_neumann.h"
}

namespace rotating_shapes
{
#include "TCD_3D/rotating_shapes.h"
}



Example_TimeCD3D::Example_TimeCD3D(
  const ParameterDatabase& user_input_parameter_db)
 : Example_NonStationary3D(user_input_parameter_db)
{
  int example_code = this->example_database["example"];
  switch(example_code)
  {
    case -2:
    {
      // linear space and time solution example 
      using namespace linear_space_time;
      exact_solution.push_back(Exact);
      boundary_conditions.push_back(BoundCondition);
      boundary_data.push_back(BoundValue);
      problem_coefficients = BilinearCoeffs;
      initialCondtion.push_back(InitialCondition);
      ExampleFile();
    }
      break;
    case -1:
    {
      using namespace quad_space_time;
      exact_solution.push_back(Exact);
      boundary_conditions.push_back(BoundCondition);
      boundary_data.push_back(BoundValue);
      problem_coefficients = BilinearCoeffs;
      initialCondtion.push_back(InitialCondition);
      ExampleFile();
    }
      break;
    case 0:
      using namespace concentration;
      exact_solution.push_back(Exact);
      boundary_conditions.push_back(BoundCondition);
      boundary_data.push_back(BoundValue);
      problem_coefficients = BilinearCoeffs;
      initialCondtion.push_back(InitialCondition);
      
      concentration::diffusionCoefficient = this->getDiffCoeff();
      ExampleFile();
      break;
    case 1:
      using namespace HeatChannel;
      exact_solution.push_back(HeatChannel::Exact);
      boundary_conditions.push_back(HeatChannel::BoundCondition);
      boundary_data.push_back(HeatChannel::BoundValue);
      problem_coefficients = HeatChannel::BilinearCoeffs;
      initialCondtion.push_back(HeatChannel::InitialCondition);
      HeatChannel::ExampleFile();
      break;
    case 2:
      using namespace block_xyzt_dirichlet;
      exact_solution.push_back(block_xyzt_dirichlet::get_c_exact_solution_values_no_t);
      boundary_conditions.push_back(block_xyzt_dirichlet::get_boundary_condition);
      boundary_data.push_back(block_xyzt_dirichlet::get_c_boundary_value_no_t);
      boundary_data.push_back(block_xyzt_dirichlet::get_c_boundary_time_derivative_no_t);
      problem_coefficients = block_xyzt_dirichlet::bilinear_coeffs;
      initialCondtion.push_back(block_xyzt_dirichlet::get_c_initial_values);
      block_xyzt_dirichlet::print_information_about_example();
      block_xyzt_dirichlet::diffCoeff = this->getDiffCoeff();
      break;
    case 3:
      using namespace block_xyzt_neumann;
      exact_solution.push_back(block_xyzt_neumann::get_c_exact_solution_values_no_t);
      boundary_conditions.push_back(block_xyzt_neumann::get_boundary_condition);
      boundary_data.push_back(block_xyzt_neumann::get_c_boundary_value_no_t);
      boundary_data.push_back(block_xyzt_neumann::get_c_boundary_time_derivative_no_t);
      problem_coefficients = block_xyzt_neumann::bilinear_coeffs;
      initialCondtion.push_back(block_xyzt_neumann::get_c_initial_values);
      block_xyzt_neumann::print_information_about_example();
      block_xyzt_neumann::diffCoeff = this->getDiffCoeff();
      break;
    case 4:
      using namespace rotating_shapes; 
      exact_solution.push_back(rotating_shapes::Exact);
      boundary_conditions.push_back(rotating_shapes::BoundCondition);
      boundary_data.push_back(rotating_shapes::BoundValue);
      problem_coefficients = rotating_shapes::BilinearCoeffs;
      initialCondtion.push_back(rotating_shapes::InitialCondition);
      rotating_shapes::print_information_about_example();
      rotating_shapes::diffCoeff = this->getDiffCoeff();
      break;
    default:
      ErrThrow("Unknown Example_TimeCD3D example!");
  }
}

void Example_TimeCD3D::do_post_processing(Time_CD3D& tcd3d) const
{
  if(post_processing_stat)
  {
    post_processing_stat(tcd3d);
  }
  else
  {
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0)
#endif
      Output::info<2>("Example_TimeCD3D","No post processing done for the current example.");
  }
}

double Example_TimeCD3D::getDiffCoeff() const
{
  return this->example_database["diffusion_coefficient"];
}
