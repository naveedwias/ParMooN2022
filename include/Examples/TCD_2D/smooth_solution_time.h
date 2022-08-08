// ======================================================================
// instationary problem
// ======================================================================

/// ========================================================================
// example file
// ========================================================================

#define __Smooth_Sol__

void ExampleFile()
{
  Output::print<1>("Example: smooth_solution_time.h");
}

constexpr bool rhs_depends_on_time = true;
constexpr bool coefficients_depend_on_time = false;

// exact solution
void Exact(double x, double y, double *values)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*100.0*x*x*(1.0-x)*(1.0-x)*y*(1.0-2.0*y)*(1.0-y); //original solution
  values[1] = -t*200.0*(1.0-x)*x*x*y*(1.0-2.0*y)*(1.0-y) + t*200.0*(1-x)*(1-x)*x*y*(1.0-2.0*y)*(1.0-y); //x-derivative
  values[2] = t*100.0*x*x*(1.0-x)*(1.0-x)*(1.0-2.0*y)*(1.0-y) - t*200.0*y*(1.0-y)*x*x*(1.0-x)*(1.0-x)-100.0*x*x*(1.0-x)*(1.0-x)*y*(1.0-2.0*y); //y-derivative
  values[3] = t*(200.0*(1.0-x)*(1.0-x)*y*(1.0-y)*(1.0-2*y)-800.0*x*(1.0-x)*y*(1.0-y)*(1.0-2.0*y)+ 200.0*x*x*y*(1.0-y)*(1.0-2.0*y) - 200.0*x*x*(1.0-x)*(1.0-x)*(1.0-2.0*y) - 400*x*x*(1.0-x)*(1.0-x)*(1.0-y)+400.0*x*x*(1.0-x)*(1.0-x)*y); //laplacian;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int, double, BoundCond &cond)
{
 cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int, double, double &value)
{
  value = 0;
}

// initial conditon
void InitialCondition(double x, double y, double *values)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = t*100.0*x*x*(1.0-x)*(1.0-x)*y*(1.0-2.0*y)*(1.0-y);
}


void BilinearCoeffs(int n_points, const double *X, const double *Y,
                    const double *const*, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->RE_NR;
  double b1=3, b2=2, c=1;
  int i;
  double *coeff;
  double x, y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    double exact[4];
    
    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = b1;
    coeff[2] = b2;
    coeff[3] = c;
    
    Exact(X[i], Y[i], exact);
    coeff[4] = -coeff[0]*exact[3]; // diffusion
    coeff[4] += coeff[1]*exact[1] + coeff[2]*exact[2]; // convection
    coeff[4] += coeff[3]*exact[0]; // reaction
    coeff[4] += 100.0*x*x*(1.0-x)*(1.0-x)*y*(1.0-2.0*y)*(1.0-y); //time_derivative
  }
}
