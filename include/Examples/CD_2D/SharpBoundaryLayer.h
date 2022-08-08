/*!
 * Example file containing assembling information and boundary conditions for a 2D convection-
 * diffusion example taken from
 *
 * Kuzmin & Moeller (2005): Algebraic Flux Correction I. Scalar Conservation Laws. (Chapter 7.4).
 *
 * The example is used to test the implementation of the stationary FEM-TVD algebraic flux correction
 * implementation scheme as it exhibits a sharp layer near the x=1 boundary.
 *
 * Use UnitSquare.GEO and UnitSquare.PRM for the geometry.
 *
 * @author Clemens Bartsch
 * @date 2015/07/22
 */


void ExampleFile()
{
  Output::print<1>("Example: SharpBoundaryLayer.h. No exact solution known, ",
                   "put to 0. Interpret error values as norm of solution.");
}

// Unknown exact solution - put to zero.
auto& Exact = unknown_solution_2d;

// kind of boundary condition
void BoundCondition(int BdComp, double, BoundCond &cond)
{
	// numbering of boundary components starts with 0 at the y==0 bdry
	// and proceeds counterclockwise for UnitSquare.PRM
  if(BdComp==2)
    cond = NEUMANN;
  else
    cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
	//on one Dirichlet boundary: jump function
	if(BdComp == 3){
		if(Param < 0.5)
			value = 1;
		else
			value = 0;
	}
	//on the other boundaries zero
	else
		value = 0;
}

void BilinearCoeffs(int n_points, const double *, const double *,
                    const double *const*parameters, double **coeffs)
{
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = 0.001; //diffusion coefficient

    bool use_spatially_varying_convection = false;
    if (use_spatially_varying_convection)
    {
    coeffs[i][1] = parameters[i][0]; //std::cos(M_PI/18);//convection in x direction: std::cos(10 deg)
    coeffs[i][2] = parameters[i][1]; //std::sin(M_PI/18);//convection in y direction: std::sin(10 deg)
    }
    else
    {
    coeffs[i][1] = std::cos(M_PI/18);//convection in x direction: std::cos(10 deg)
    coeffs[i][2] = std::sin(M_PI/18);//convection in y direction: std::sin(10 deg)
    }

    coeffs[i][3] = 0; //reaction coefficient

    coeffs[i][4] = 0; //rhs, outer body force coefficients
  }
}
