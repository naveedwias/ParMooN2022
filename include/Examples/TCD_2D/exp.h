// ======================================================================
// instationary problem
// ======================================================================

/// ========================================================================
// example file
// ========================================================================

#define __SIN3__

void ExampleFile()
{
  Output::print<1>("Example: exp.h");
}

constexpr bool rhs_depends_on_time = true;
constexpr bool coefficients_depend_on_time = false;

// exact solution
void Exact(double x, double y, double *values)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = std::exp(std::sin(2*M_PI*t))*(std::sin(2*M_PI*x)*std::sin(2*M_PI*y));
  values[1] = std::exp(std::sin(2*M_PI*t))*2*M_PI*std::cos(2*M_PI*x)*std::sin(2*M_PI*y);
  values[2] = std::exp(std::sin(2*M_PI*t))*2*M_PI*std::sin(2*M_PI*x)*std::cos(2*M_PI*y);
  values[3] = 0;
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
  values[0] = std::exp(std::sin(2*M_PI*t))*(std::sin(2*M_PI*x)*std::sin(2*M_PI*y));
}


void BilinearCoeffs(int n_points, const double *X, const double *Y,
                    const double *const*, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->RE_NR;
  double b1=1, b2=-1, c=1;
  int i;
  double *coeff;
  double x, y;
  double t = TDatabase::TimeDB->CURRENTTIME;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = b1;
    coeff[2] = b2;
    coeff[3] = c;

    coeff[4] = std::exp(std::sin(2*M_PI*t))*2*M_PI*std::cos(2*M_PI*t)*(std::sin(2*M_PI*x)*std::sin(2*M_PI*y))
	- eps * std::exp(std::sin(2*M_PI*t))*4*M_PI*M_PI*(-std::sin(2*M_PI*x)*std::sin(2*M_PI*y)-std::sin(2*M_PI*x)*std::sin(2*M_PI*y))
       + b1 * std::exp(std::sin(2*M_PI*t))*2*M_PI*std::cos(2*M_PI*x)*std::sin(2*M_PI*y)
       + b2 * std::exp(std::sin(2*M_PI*t))*2*M_PI*std::sin(2*M_PI*x)*std::cos(2*M_PI*y)
	+ c *  std::exp(std::sin(2*M_PI*t))*(std::sin(2*M_PI*x)*std::sin(2*M_PI*y));
  }
}
