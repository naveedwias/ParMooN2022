// ==========================================================================
// instationary problem
// ==========================================================================

//===========================================================================
// example file
// =========================================================================
#include <MooNMD_Io.h>

double diffusionCoefficient;

void ExampleFile()
{
  Output::print("concentrationOfSpecies_3d.h");
}

// exact solution
auto& Exact = unknown_solution_3d;

// initial conditon
void InitialCondition(double, double, double, double *values)
{
    values[0] = 0;
}

// kind of boundary condition
void BoundCondition(int, double x, double y, double z, BoundCond &cond)
{
  double eps = 1e-6;
  
  if ((std::abs(x-1)<eps)&& (y>=3.0/8.0-eps) && (y<=0.5+eps)
    && (z>=0.5-eps) && (z<=5.0/8.0+eps))
    cond = NEUMANN;
  else
    cond = DIRICHLET;
}


// value of boundary condition
void BoundValue(int, double x, double y, double z, double &value)
{

  double t = TDatabase::TimeDB->CURRENTTIME;
  double eps = 1e-6;

  if ((x==0) && (y>=5.0/8.0-eps) && (y<=0.75+eps) &&
    (z>= 5.0/8.0-eps) && (z<= 0.75 +eps))
  {
    if ((t>=0) && t<=1)
      value =  std::sin(M_PI*t/2.0)*1.0;
    else
    {
      if ((t>1)&&t<=2)
	value = 1.0;
      else
      {
	if ((t>2)&&(t<=3))
	  value = std::sin(M_PI*(t-1)/2.0)*1.0;
	else
	{
	  value = 0;
	}
      }
    }
  }
  else
    value = 0.0;
}

void BilinearCoeffs(int n_points, const double *X, const double *Y,
                    const double *Z, const double *const*, double **coeffs)
{
  double eps = diffusionCoefficient;
  int i;
  double *coeff;                                  // *param;
  double x, y, z, c, a[3], b[3], s[3], h;
  // double t = TDatabase::TimeDB->CURRENTTIME;
  
  // norm of convection vector
  c = 1.0/16 + 1.0/64 + 1;
  c = std::sqrt(c);

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    // param = parameters[i];

    x = X[i];
    y = Y[i];
    z = Z[i];

    // diffusion
    coeff[0] = eps;
    // convection in x direction
    coeff[1] = 1;
    // convection in y direction
    coeff[2] = -0.25;
    // convection in z direction
    coeff[3] = -0.125;
    // reaction
    // vector from center of outflow (1,7/16,9/16) to (x,y,z)
    a[0] = x-1;
    a[1] = y - 7.0/16;
    a[2] = z - 9.0/ 16;
    // vector from center of inflow (0,11/16,11/16) to (x,y,z)
    b[0] = x;
    b[1] = y - 11.0/16;
    b[2] = z - 11.0/16;
    // cross product
    s[0] = a[1] * b[2] - a[2] * b[1];
    s[1] = a[2] * b[0] - a[0] * b[2];
    s[2] = a[0] * b[1] - a[1] * b[0];
    // norm of cross product = 2 * length of convection * height 
    // area of parallelogram = 2 * area of triangle = 2 * c * h /2 = c * h
    h = std::sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
    // 2 * height
    h /= c;

    if (h<=0.1)
    {
      coeff[4] = 1;
    }
    else
	coeff[4] = 0;
    // rhs
    coeff[5] = 0;
    // rhs from previous time step
    coeff[6] = 0;
  }
}
