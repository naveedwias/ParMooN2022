// ======================================================================
// instationary problem
// ======================================================================

/// ========================================================================
// example file
// ========================================================================

#define __SIN3__

void ExampleFile()
{
  Output::print<1>("Example: Sin3.h");
}

constexpr bool rhs_depends_on_time = true;
constexpr bool coefficients_depend_on_time = false;

// exact solution
void Exact(double x, double y, double *values)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = std::sin(t)*(std::sin(2*M_PI*x)*std::sin(2*M_PI*y)+1);
  values[1] = std::sin(t)*2*M_PI*std::cos(2*M_PI*x)*std::sin(2*M_PI*y);
  values[2] = std::sin(t)*2*M_PI*std::sin(2*M_PI*x)*std::cos(2*M_PI*y);
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
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;
    value = std::sin(t);
}

// initial conditon
void InitialCondition(double x, double y, double *values)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = std::sin(t)*(std::sin(2*M_PI*x)*std::sin(2*M_PI*y)+1);
}


void BilinearCoeffs(int n_points, const double *X, const double *Y,
                    const double *const*, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->PE_NR;
  double b1=1, b2=2, c=1;
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

    coeff[4] = std::cos(t)*(std::sin(2*M_PI*x)*std::sin(2*M_PI*y)+1)
	- eps * std::sin(t)*4*M_PI*M_PI*(-std::sin(2*M_PI*x)*std::sin(2*M_PI*y)-std::sin(2*M_PI*x)*std::sin(2*M_PI*y))
       + b1 * std::sin(t)*2*M_PI*std::cos(2*M_PI*x)*std::sin(2*M_PI*y)
       + b2 * std::sin(t)*2*M_PI*std::sin(2*M_PI*x)*std::cos(2*M_PI*y)
	+ c *  std::sin(t)*(std::sin(2*M_PI*x)*std::sin(2*M_PI*y)+1);
  }
}
