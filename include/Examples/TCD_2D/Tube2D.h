/*!
 * A two-dimensional tube with given C = 1 on the left boundary, and neumann boundary condition elsewhere
 */

#ifndef __TUBE2D__
#define __TUBE2D__

double diffusion_coefficient = -1.;

/** Print out some information on the example file. */
void ExampleFile()
{
  Output::print("Example: Tube2D.h");

  // TODO This global variable is a total mess - rework!
  TDatabase::ParamDB->RE_NR = 1./diffusion_coefficient;
  Output::print("TDatabase::ParamDB->RE_NR was set to",
		TDatabase::ParamDB->RE_NR);
}

constexpr bool rhs_depends_on_time = false;
constexpr bool coefficients_depend_on_time = false;
/** The exact solution */
auto& Exact = unknown_solution_2d;


/** The type of the boundary condition - Dirichlet on all components. */
void BoundCondition(int component, double, BoundCond &cond)
{
  if (component==0) {
    cond = DIRICHLET;
  } else {
    cond = NEUMANN;
  }
}


/** The value of boundary condition - 0 everywhere*/
void BoundValue(int component, double, double &value)
{
  if (component==0) {
    value = 1.;
  } else {
    value = 0.;
  }
}


/** The initial condition - describes the initial shape of the
 *  bodies.*/
void InitialCondition(double, double, double *values)
{
  
  values[0] = 0;

}

/**
 * Coefficient function, used in the assembling process of the CDR problem.
 */
void BilinearCoeffs(int n_points, const double*, const double*,
                    const double*const*, double **coeffs)
{

  double eps = diffusion_coefficient;

  int i;
  double *coeff;                                  // *param;
  // double t = TDatabase::TimeDB->CURRENTTIME;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    // param = parameters[i];

    //x = X[i];
    //y = Y[i];

    // diffusion
    coeff[0] = eps;
    // convection in x direction
    coeff[1] = 0.;
    // convection in y direction
    coeff[2] = 1.;
    // reaction
    coeff[3] = 0;
    // rhs
    coeff[4] = 0;
    // rhs from previous time step
    coeff[5] = 0;
  }
}


#endif
