#include <Example_CD2D.h>
#include <ConvectionDiffusion.h>

#include <Database.h>
#include <FEFunction2D.h>
#include "FEDatabase.h"
#include <SquareMatrix2D.h>
#include "MainUtilities.h"
#include "BaseCell.h"

#include <BoundEdge.h>
#include <IsoBoundEdge.h>
#include <BoundComp.h>


#include <algorithm>
#include <cmath>
#include <string.h>

//===========================================================
// examples for stationary convection-diffusion-reaction
// problems
//===========================================================

namespace sine_laplace
{
  #include "CD_2D/SineLaplace.h"
}

namespace two_interior_layers
{
  #include "CD_2D/TwoInteriorLayers.h"
}

namespace hemker_1996
{
  #include "CD_2D/Hemker1996.h"
}

namespace sharp_boundary_layer
{
  #include "CD_2D/SharpBoundaryLayer.h"
}

namespace smooth_solution
{
  #include "CD_2D/Smooth.h"
}

namespace HMM1986
{
  #include "CD_2D/HMM1986.h"
}

namespace hump
{
  #include "CD_2D/hump.h"
}

namespace boundary_layer_known
{
  #include "CD_2D/boundary_layer_known.h"
}

namespace example_1_ABR17
{
  #include "CD_2D/example_1_ABR17.h"
}

namespace rotating_convection_field
{
  #include "CD_2D/rotating_convection_field.h"
}

namespace polynomial
{
  #include "CD_2D/Polynomial.h"
}

namespace convection_reaction
{
  #include "CD_2D/convection_reaction.h"
}

namespace jump_polynomial
{
  #include "CD_2D/jump_polynomial.h"
}

namespace three_interior_layers_inlet_M
{
  #include "CD_2D/ThreeInteriorLayers_inlet_M.h"
}

namespace three_interior_layers_inlet_M_bump
{
  #include "CD_2D/ThreeInteriorLayers_inlet_M_bump.h"
}

namespace dominant_reaction
{
  #include "CD_2D/dominant_reaction.h"
}

namespace HMM1986_modified_bc
{
  #include "CD_2D/HMM1986_modified_bc.h"
}

namespace ParabolicLayers
{
  #include "CD_2D/ParabolicLayers.h"
}

namespace steady_circular_advection
{
  #include "CD_2D/steady_circular_advection.h"
}

Example_CD2D::Example_CD2D(const ParameterDatabase& user_input_parameter_db) 
 : Example2D(user_input_parameter_db)
{
  int example_code = this->example_database["example"];
  int do_poisson = TDatabase::ParamDB->INTERNAL_COERCIVITY;
  switch( example_code)
  {
    case 0:
      /** exact_solution */
      exact_solution.push_back( sine_laplace::Exact );

      /** boundary condition */
      boundary_conditions.push_back( sine_laplace::BoundCondition );

      /** boundary values */
      boundary_data.push_back( sine_laplace::BoundValue );

      /** coefficients */
      problem_coefficients = sine_laplace::BilinearCoeffs;

      sine_laplace::ExampleFile();

      sine_laplace::DIFFUSION = this->get_nu();

      break;

    case 1:
      /** exact_solution */
      exact_solution.push_back( two_interior_layers::Exact );

      /** boundary condition */
      boundary_conditions.push_back( two_interior_layers::BoundCondition );

      /** boundary values */
      boundary_data.push_back( two_interior_layers::BoundValue );

      /** coefficients */
      problem_coefficients = two_interior_layers::BilinearCoeffs;

      two_interior_layers::ExampleFile();

      two_interior_layers::DIFFUSION = this->get_nu();

      break;

    case 2:
    case 20:
      /** exact_solution */
      exact_solution.push_back( hemker_1996::Exact );

      /** boundary condition */
      boundary_conditions.push_back( hemker_1996::BoundCondition );

      /** boundary values */
      boundary_data.push_back( hemker_1996::BoundValue );

      /** coefficients */
      problem_coefficients = hemker_1996::BilinearCoeffs;

      post_processing_stat = hemker_1996::hemker_postprocessing;
      hemker_1996::DIFFUSION = this->get_nu();
      hemker_1996::modified_Hemker = (example_code == 20);
      hemker_1996::ExampleFile();
      break;

    case 3:
      /** exact_solution */
      exact_solution.push_back( sharp_boundary_layer::Exact );

      /** boundary condition */
      boundary_conditions.push_back( sharp_boundary_layer::BoundCondition );

      /** boundary values */
      boundary_data.push_back( sharp_boundary_layer::BoundValue );

      /** coefficients */
      problem_coefficients = sharp_boundary_layer::BilinearCoeffs;

      sharp_boundary_layer::ExampleFile();
      break;

    case 4:
      /** exact_solution */
      exact_solution.push_back( smooth_solution::Exact );

      if(do_poisson == 27)
      {
        /** boundary condition */
        boundary_conditions.push_back( smooth_solution::BoundCondition_Poisson );
        /** boundary values */
        boundary_data.push_back( smooth_solution::BoundValue_Poisson );
      }
      else
      {
        /** boundary condition */
        boundary_conditions.push_back( smooth_solution::BoundCondition );
        /** boundary values */
        boundary_data.push_back( smooth_solution::BoundValue );
      }

      /** coefficients */
      problem_coefficients = smooth_solution::BilinearCoeffs;

      smooth_solution::ExampleFile();

      smooth_solution::DIFFUSION = this->get_nu();

      break;

    case 5:
    case 13:
      /** exact_solution */
      exact_solution.push_back( HMM1986::Exact );

      /** boundary condition */
      boundary_conditions.push_back( HMM1986::BoundCondition );

      /** boundary values */
      boundary_data.push_back( HMM1986::BoundValue );

      /** coefficients */
      problem_coefficients = HMM1986::BilinearCoeffs;

      HMM1986::modified_HMM = (example_code == 5) ? false : true;

      post_processing_stat = HMM1986::hmm1986_postprocessing;
      HMM1986::DIFFUSION = this->get_nu();
      
      HMM1986::ExampleFile();
      break;

    case 6:
      /** exact_solution */
      exact_solution.push_back( hump::Exact );

      /** boundary condition */
      boundary_conditions.push_back( hump::BoundCondition );

      /** boundary values */
      boundary_data.push_back( hump::BoundValue );

      /** coefficients */
      problem_coefficients = hump::BilinearCoeffs;

      hump::ExampleFile();

      hump::DIFFUSION = this->get_nu();

      break;

    case 7:
      /** exact_solution */
      exact_solution.push_back( boundary_layer_known::Exact );

      /** boundary condition */
      boundary_conditions.push_back( boundary_layer_known::BoundCondition );

      /** boundary values */
      boundary_data.push_back( boundary_layer_known::BoundValue );

      /** coefficients */
      problem_coefficients = boundary_layer_known::BilinearCoeffs;

      boundary_layer_known::ExampleFile();

      boundary_layer_known::DIFFUSION = this->get_nu();

      break;

    case 8:
      /** exact_solution */
      exact_solution.push_back( example_1_ABR17::Exact );

      /** boundary condition */
      boundary_conditions.push_back( example_1_ABR17::BoundCondition );

      /** boundary values */
      boundary_data.push_back( example_1_ABR17::BoundValue );

      /** coefficients */
      problem_coefficients = example_1_ABR17::BilinearCoeffs;

      example_1_ABR17::ExampleFile();

      example_1_ABR17::DIFFUSION = this->get_nu();
      break;

    case 9:
      /** exact_solution */
      exact_solution.push_back( rotating_convection_field::Exact );

      /** boundary condition */
      boundary_conditions.push_back( rotating_convection_field::BoundCondition );

      /** boundary values */
      boundary_data.push_back( rotating_convection_field::BoundValue );

      /** coefficients */
      problem_coefficients = rotating_convection_field::BilinearCoeffs;

      rotating_convection_field::ExampleFile();

      rotating_convection_field::DIFFUSION = this->get_nu();

      break;

    case 10:  // polynomial solution
      /** exact_solution */
      exact_solution.push_back( polynomial::Exact );

      /** boundary condition */
      boundary_conditions.push_back( polynomial::BoundCondition );

      /** boundary values */
      boundary_data.push_back( polynomial::BoundValue );

      /** coefficients */
      problem_coefficients = polynomial::BilinearCoeffs;
      polynomial::deg = this->example_database["degree_polynomial"];
      polynomial::ExampleFile();

      polynomial::DIFFUSION = this->get_nu();

      break;

    case 11:  // pure convection-reaction
      /** exact_solution */
      exact_solution.push_back( convection_reaction::Exact );

      /** boundary condition */
      boundary_conditions.push_back( convection_reaction::BoundCondition );

      /** boundary values */
      boundary_data.push_back( convection_reaction::BoundValue );

      /** coefficients */
      problem_coefficients = convection_reaction::BilinearCoeffs;
      convection_reaction::ExampleFile();
      break;

    case 12:  // jump_polynomial
      /** exact_solution */
      exact_solution.push_back( jump_polynomial::Exact );

      /** boundary condition */
      boundary_conditions.push_back( jump_polynomial::BoundCondition );

      /** boundary values */
      boundary_data.push_back( jump_polynomial::BoundValue );

      /** coefficients */
      problem_coefficients = jump_polynomial::BilinearCoeffs;
      jump_polynomial::deg = this->example_database["degree_polynomial"];
      jump_polynomial::ExampleFile();
      break;

    case 14: // three_interior_layers_inlet_M
      /** exact_solution */
      exact_solution.push_back( three_interior_layers_inlet_M::Exact );

      /** boundary condition */
      boundary_conditions.push_back( three_interior_layers_inlet_M::BoundCondition );

      /** boundary values */
      boundary_data.push_back( three_interior_layers_inlet_M::BoundValue );

      /** coefficients */
      problem_coefficients = three_interior_layers_inlet_M::BilinearCoeffs;
      three_interior_layers_inlet_M::ExampleFile();
      three_interior_layers_inlet_M::DIFFUSION = this->get_nu();
      post_processing_stat = three_interior_layers_inlet_M::three_interior_layers_inlet_M_postprocessing; 
      break;

    case 15: // three_interior_layers_inlet_M_bump
      /** exact_solution */
      exact_solution.push_back( three_interior_layers_inlet_M_bump::Exact );

      /** boundary condition */
      boundary_conditions.push_back( three_interior_layers_inlet_M_bump::BoundCondition );

      /** boundary values */
      boundary_data.push_back( three_interior_layers_inlet_M_bump::BoundValue );

      /** coefficients */
      problem_coefficients = three_interior_layers_inlet_M_bump::BilinearCoeffs;
      three_interior_layers_inlet_M_bump::ExampleFile();
      three_interior_layers_inlet_M_bump::DIFFUSION = this->get_nu();
      post_processing_stat = three_interior_layers_inlet_M_bump::three_interior_layers_inlet_M_bump_postprocessing;

      break;
    case 16: // dominant_reaction
      /** exact_solution */
      exact_solution.push_back(dominant_reaction::Exact );

      /** boundary condition */
      boundary_conditions.push_back(dominant_reaction::BoundCondition );

      /** boundary values */
      boundary_data.push_back(dominant_reaction::BoundValue );

      /** coefficients */
      problem_coefficients = dominant_reaction::BilinearCoeffs;
      dominant_reaction::ExampleFile();
      dominant_reaction::DIFFUSION = this->get_nu();
      //post_processing_stat = dominant_reaction::dominant_reaction_postprocessing;

      break;

    case 17: // HMM1986_modified_bc
      /** exact_solution */
      exact_solution.push_back(HMM1986_modified_bc::Exact );

      /** boundary condition */
      boundary_conditions.push_back(HMM1986_modified_bc::BoundCondition );

      /** boundary values */
      boundary_data.push_back(HMM1986_modified_bc::BoundValue );

      /** coefficients */
      problem_coefficients = HMM1986_modified_bc::BilinearCoeffs;
      HMM1986_modified_bc::ExampleFile();
      HMM1986_modified_bc::DIFFUSION = this->get_nu();
      //post_processing_stat = HMM1986_modified_bc::HMM1986_modified_bc_postprocessing;

      break;

    case 18:
      /** exact_solution */
      exact_solution.push_back( ParabolicLayers::Exact );

      /** boundary condition */
      boundary_conditions.push_back( ParabolicLayers::BoundCondition );

      /** boundary values */
      boundary_data.push_back( ParabolicLayers::BoundValue );

      /** coefficients */
      problem_coefficients = ParabolicLayers::BilinearCoeffs;

      ParabolicLayers::ExampleFile();

      ParabolicLayers::DIFFUSION = this->get_nu();

      break;    
    case 19:
      /** exact_solution */
      exact_solution.push_back( steady_circular_advection::Exact );

     if(do_poisson == 27)
      {
        /** boundary condition */
        boundary_conditions.push_back( steady_circular_advection::BoundCondition_Poisson );
        /** boundary values */
        boundary_data.push_back( steady_circular_advection::BoundValue_Poisson );
      }
      else
      {
        /** boundary condition */
        boundary_conditions.push_back( steady_circular_advection::BoundCondition );
        /** boundary values */
        boundary_data.push_back( steady_circular_advection::BoundValue );
      }

      /** coefficients */
      problem_coefficients = steady_circular_advection::BilinearCoeffs;

      steady_circular_advection::ExampleFile();

      steady_circular_advection::DIFFUSION = this->get_nu();
      break;
    default:
      ErrThrow("Unknown name of the convection-diffusion (CD2D) example!", 
               example_code);
  }
}

Example_CD2D::Example_CD2D(const std::vector<DoubleFunct2D*>& exact,
                           const std::vector<BoundCondFunct2D*>& bc,
                           const std::vector<BoundValueFunct2D*>& bd,
                           const CoeffFct2D& coeffs, double nu)
: Example2D(exact, bc, bd, coeffs)
{
  this->example_database["diffusion_coefficient"] = nu;
}


void Example_CD2D::do_post_processing(ConvectionDiffusion<2>& cd2d) const
{
  if(post_processing_stat)
  {
    post_processing_stat(cd2d);
  }
  else
  {
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0)
#endif
      Output::info<2>("Example_CD2D","No post processing done for the current example.");
  }
}

double Example_CD2D::get_nu() const
{
  double diffusion_coefficient = this->example_database["diffusion_coefficient"];
  return diffusion_coefficient;
}

