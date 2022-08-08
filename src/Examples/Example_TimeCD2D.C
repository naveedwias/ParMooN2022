#include <Example_TimeCD2D.h>
#include <Database.h>
#include <FEFunction2D.h>
#include "MainUtilities.h"
#include <cmath>
#include<string>

namespace linear_space_time
{
#include "TCD_2D/linear_space_time.h"
}
namespace exp_sin_cos
{
#include "TCD_2D/exp.h"
}
namespace sin_sin_sin
{
#include "TCD_2D/Sin3.h"
}

namespace sin_cos
{
#include "TCD_2D/SinCos1.h"
}

namespace rotating_bodies_1
{
#include "TCD_2D/Rotating_Bodies.h"
}

namespace Geothermal_Energy_TCD2D
{
  #include "TCD_2D/Geothermal_Energy_TCD2D.h"
}

namespace Tube2D
{
  #include "TCD_2D/Tube2D.h"
}

namespace smooth_solution_time
{
#include "TCD_2D/smooth_solution_time.h"
}

namespace time_dominated
{
#include "TCD_2D/time_dominated.h"
}

namespace rectangle_xyt_dirichlet
{
#include "TCD_2D/rectangle_xyt_dirichlet.h"
}

namespace rectangle_xyt_neumann
{
#include "TCD_2D/rectangle_xyt_neumann.h"
}

#include "TCD_2D/traveling_wave_sqrt.h"  // 10


#include "TCD_2D/JohnMaubachTobiska1997inst.h"  // 11

#include "TCD_2D/Sin3_0_pi.h"  // 12

#include "TCD_2D/exp_sin3_0_pi.h"  // 13

//=======================================================

Example_TimeCD2D::Example_TimeCD2D(
  const ParameterDatabase& user_input_parameter_db)
 : Example_NonStationary2D(user_input_parameter_db)
{
  int example_code = this->example_database["example"];
  switch(example_code)
  {
    case -1:
    {
      /**Exact solution"**/
      exact_solution.push_back(linear_space_time::Exact);
      /** boundary condition */
      boundary_conditions.push_back( linear_space_time::BoundCondition );
      
      /** boundary values */
      boundary_data.push_back( linear_space_time::BoundValue );
      
      /** coefficients */
      problem_coefficients = linear_space_time::BilinearCoeffs;
      
      /** Initial condition*/
      initialCondition.push_back(linear_space_time::InitialCondition);
     linear_space_time::ExampleFile();

     this->timeDependentRhs = linear_space_time::rhs_depends_on_time;
     this->timeDependentCoeffs=linear_space_time::coefficients_depend_on_time;
    }
    break;
    case 0:
      /**Exact solution"**/
      exact_solution.push_back(exp_sin_cos::Exact);
      /** boundary condition */
      boundary_conditions.push_back( exp_sin_cos::BoundCondition );
      
      /** boundary values */
      boundary_data.push_back( exp_sin_cos::BoundValue );
      
      /** coefficients */
      problem_coefficients = exp_sin_cos::BilinearCoeffs;
      
      /** Initial condition*/
      initialCondition.push_back(exp_sin_cos::InitialCondition);
     exp_sin_cos::ExampleFile();
     
     this->timeDependentRhs = exp_sin_cos::rhs_depends_on_time;
     this->timeDependentCoeffs=exp_sin_cos::coefficients_depend_on_time;
     break;
    case 1:
      /**Exact solution"**/
       exact_solution.push_back(sin_sin_sin::Exact);
      /** boundary condition */
      boundary_conditions.push_back( sin_sin_sin::BoundCondition );
      
      /** boundary values */
      boundary_data.push_back( sin_sin_sin::BoundValue );
      
      /** coefficients */
      problem_coefficients = sin_sin_sin::BilinearCoeffs;
      
      /** Initial condition*/
      initialCondition.push_back(sin_sin_sin::InitialCondition);
      sin_sin_sin::ExampleFile();
      
      this->timeDependentRhs = sin_sin_sin::rhs_depends_on_time;
     this->timeDependentCoeffs=sin_sin_sin::coefficients_depend_on_time;
      break;
    case 2:
      /**Exact solution"**/
      exact_solution.push_back(sin_cos::Exact);

      /** boundary condition */
      boundary_conditions.push_back( sin_cos::BoundCondition );
      
      /** boundary values */
      boundary_data.push_back( sin_cos::BoundValue );
      
      /** coefficients */
      problem_coefficients = sin_cos::BilinearCoeffs;
      
      /** Initial condition*/
      initialCondition.push_back(sin_cos::InitialCondition);
      sin_cos::ExampleFile();      
      break;
    case 3:
      /**Exact solution"**/
      exact_solution.push_back(rotating_bodies_1::Exact);

      /** boundary condition */
      boundary_conditions.push_back( rotating_bodies_1::BoundCondition );

      /** boundary values */
      boundary_data.push_back( rotating_bodies_1::BoundValue );

      /** coefficients */
      problem_coefficients = rotating_bodies_1::BilinearCoeffs;

      /** Initial condition*/
      initialCondition.push_back(rotating_bodies_1::InitialCondition);

      // Print some example specific information.
      rotating_bodies_1::ExampleFile();
      break;
 case 4:
      /**Exact solution"**/
      exact_solution.push_back(Geothermal_Energy_TCD2D::Exact);

      /** boundary condition */
      boundary_conditions.push_back( Geothermal_Energy_TCD2D::BoundCondition );

      /** boundary values */
      boundary_data.push_back( Geothermal_Energy_TCD2D::BoundValue );

      /** coefficients */
      problem_coefficients = Geothermal_Energy_TCD2D::BilinearCoeffs;

      /** Initial condition*/
      initialCondition.push_back(Geothermal_Energy_TCD2D::InitialCondition);
 
      // Print some example specific information.
      Geothermal_Energy_TCD2D::ExampleFile();
      break;
      
  case 5:
      /**Exact solution"**/
      exact_solution.push_back(Tube2D::Exact);

      /** boundary condition */
      boundary_conditions.push_back( Tube2D::BoundCondition );

      /** boundary values */
      boundary_data.push_back( Tube2D::BoundValue );

      /** coefficients */
      problem_coefficients = Tube2D::BilinearCoeffs;

      /** Initial condition*/
      initialCondition.push_back(Tube2D::InitialCondition);

      Tube2D::diffusion_coefficient = getDiffCoeff();
      
      // Print some example specific information.
      Tube2D::ExampleFile();
      break;
      
    case 6:
      /**Exact solution"**/
      exact_solution.push_back(smooth_solution_time::Exact);
      /** boundary condition */
      boundary_conditions.push_back( smooth_solution_time::BoundCondition );
      
      /** boundary values */
      boundary_data.push_back( smooth_solution_time::BoundValue );
      
      /** coefficients */
      problem_coefficients = smooth_solution_time::BilinearCoeffs;
      
      /** Initial condition*/
      initialCondition.push_back(smooth_solution_time::InitialCondition);
       smooth_solution_time::ExampleFile();
     
     this->timeDependentRhs = smooth_solution_time::rhs_depends_on_time;
     this->timeDependentCoeffs = smooth_solution_time::coefficients_depend_on_time;
     break;
     
    case 7:
      /**Exact solution"**/
      exact_solution.push_back(time_dominated::Exact);
      /** boundary condition */
      boundary_conditions.push_back( time_dominated::BoundCondition );
      
      /** boundary values */
      boundary_data.push_back( time_dominated::BoundValue );
      
      /** coefficients */
      problem_coefficients = time_dominated::BilinearCoeffs;
      
      /** Initial condition*/
      initialCondition.push_back(time_dominated::InitialCondition);
       time_dominated::ExampleFile();
     
     this->timeDependentRhs = time_dominated::rhs_depends_on_time;
     this->timeDependentCoeffs = time_dominated::coefficients_depend_on_time;
     break;
     
    case 8:
      /**Exact solution"**/
      exact_solution.push_back(rectangle_xyt_dirichlet::get_c_exact_solution_values_no_t);
      /** boundary condition */
      boundary_conditions.push_back( rectangle_xyt_dirichlet::get_boundary_condition );
      
      /** boundary values */
      boundary_data.push_back( rectangle_xyt_dirichlet::get_c_boundary_value_no_t );
      boundary_data.push_back( rectangle_xyt_dirichlet::get_c_boundary_time_derivative_no_t );
      
      /** coefficients */
      problem_coefficients = rectangle_xyt_dirichlet::bilinear_coeffs;
      
      /** Initial condition*/
      initialCondition.push_back(rectangle_xyt_dirichlet::get_c_initial_values);
      
      this->timeDependentRhs = rectangle_xyt_dirichlet::rhs_depends_on_time;
      this->timeDependentCoeffs = rectangle_xyt_dirichlet::coefficients_depend_on_time;
      rectangle_xyt_dirichlet::diffCoeff = getDiffCoeff();
      
      rectangle_xyt_dirichlet::print_information_about_example();
     break;
     
    case 9:
      /**Exact solution"**/
      exact_solution.push_back(rectangle_xyt_neumann::get_c_exact_solution_values_no_t);
      /** boundary condition */
      boundary_conditions.push_back( rectangle_xyt_neumann::get_boundary_condition );
      
      /** boundary values */
      boundary_data.push_back( rectangle_xyt_neumann::get_c_boundary_value_no_t );
      boundary_data.push_back( rectangle_xyt_neumann::get_c_boundary_time_derivative_no_t );
      
      /** coefficients */
      problem_coefficients = rectangle_xyt_neumann::bilinear_coeffs;
      
      /** Initial condition*/
      initialCondition.push_back(rectangle_xyt_neumann::get_c_initial_values);
      
      this->timeDependentRhs = rectangle_xyt_neumann::rhs_depends_on_time;
      this->timeDependentCoeffs = rectangle_xyt_neumann::coefficients_depend_on_time;
      rectangle_xyt_neumann::diffCoeff = getDiffCoeff();
      
      rectangle_xyt_neumann::print_information_about_example();
      
     break;

    case 10:
    {
      using namespace traveling_wave;
      /**Exact solution"**/
      exact_solution.push_back( Exact );

      /** Initial condition*/
      initialCondition.push_back( InitialCondition );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );

      /** boundary values */
      boundary_data.push_back( BoundValue );

      /** coefficients */
      problem_coefficients = BilinearCoeffs;

      DiffCoeff = getDiffCoeff();
      ExampleFile();
     break;
    }

    case 11:
    {
      using namespace john_maubach_tobiska_inst;
      /**Exact solution"**/
      exact_solution.push_back( Exact );

      /** Initial condition*/
      initialCondition.push_back( InitialCondition );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );

      /** boundary values */
      boundary_data.push_back( BoundValue );

      /** coefficients */
      problem_coefficients = BilinearCoeffs;

      DiffCoeff = getDiffCoeff();
      ExampleFile();
     break;
    }

    case 12:
    {
      using namespace sin3_0_pi;
      /**Exact solution"**/
      exact_solution.push_back( Exact );

      /** Initial condition*/
      initialCondition.push_back( InitialCondition );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );

      /** boundary values */
      boundary_data.push_back( BoundValue );

      /** coefficients */
      problem_coefficients = BilinearCoeffs;

      DiffCoeff = getDiffCoeff();
      ExampleFile();
     break;
    }

    case 13:
    {
      using namespace exp_sin3_0_pi;
      /**Exact solution"**/
      exact_solution.push_back( Exact );

      /** Initial condition*/
      initialCondition.push_back( InitialCondition );

      /** boundary condition */
      boundary_conditions.push_back( BoundCondition );

      /** boundary values */
      boundary_data.push_back( BoundValue );

      /** coefficients */
      problem_coefficients = BilinearCoeffs;

      DiffCoeff = getDiffCoeff();
      ExampleFile();
     break;
    }

    default:
      ErrThrow("Unknown name of the transient-convection-diffusion (Time_CD2D) "
               "example!", example_code);
  }
}

void Example_TimeCD2D::do_post_processing(Time_CD2D& tcd2d) const
{
  if(post_processing_stat)
  {
    post_processing_stat(tcd2d);
  }
  else
  {
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0)
#endif
      Output::info<2>("Example_TimeCD2D","No post processing done for the current example.");
  }
}

double Example_TimeCD2D::getDiffCoeff() const
{
  return this->example_database["diffusion_coefficient"];
}

