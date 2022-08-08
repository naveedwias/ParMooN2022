// ======================================================================
// instationary problem
// ======================================================================

#define __SINCOS1__

/// ========================================================================
// example file
// ========================================================================

void ExampleFile()
{
  Output::print<1>("Example: SinCos1.h");
}

constexpr bool rhs_depends_on_time = true;
constexpr bool coefficients_depend_on_time = false;

// exact solution
void Exact(double x, double y, double *values)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = std::sin(t*t)*std::cos(x*y*y);
  values[1] = -std::sin(t*t)*std::sin(x*y*y)*y*y;
  values[2] = -std::sin(t*t)*std::sin(x*y*y)*2*x*y;
  values[3] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;
  switch(BdComp)
  {
    case 0: value = std::sin(t*t);
            break;
    case 1: value = std::sin(t*t)*std::cos(Param*Param);
            break;
    case 2: value = std::sin(t*t)*std::cos(1-Param);
            break;
    case 3: value = std::sin(t*t);
            break;
  }
}

// initial conditon
void InitialCondition(double x, double y, double *values)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = std::sin(t*t)*std::cos(x*y*y);
}

void BilinearCoeffs(int n_points, const double *X, const double *Y,
                    const double *const*, double **coeffs)
{
  double eps=1./TDatabase::ParamDB->PE_NR;
  double b1=2., b2=-1., c=1.;
  int i;
  double *coeff;
  double x, y;
  double t = TDatabase::TimeDB->CURRENTTIME;
  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  // previous discrete time
  tau = t-tau;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];   

    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = b1;
    coeff[2] = b2;
    coeff[3] = c;

    coeff[4] = std::cos(t*t)*2*t*std::cos(x*y*y)
       - eps* std::sin(t*t)*(-std::cos(x*y*y)*(y*y*y*y+4*x*x*y*y) - std::sin(x*y*y)*2*x)  
       - b1*std::sin(t*t)*std::sin(x*y*y)*y*y - b2*std::sin(t*t)*std::sin(x*y*y)*2*x*y 
       + c*std::sin(t*t)*std::cos(x*y*y);
    // rhs from previous time step
    coeff[5] = std::cos(tau*tau)*2*tau*std::cos(x*y*y)
       - eps* std::sin(tau*tau)*(-std::cos(x*y*y)*(y*y*y*y+4*x*x*y*y) - std::sin(x*y*y)*2*x)  
       - b1*std::sin(tau*tau)*std::sin(x*y*y)*y*y - b2*std::sin(tau*tau)*std::sin(x*y*y)*2*x*y 
       + c*std::sin(tau*tau)*std::cos(x*y*y);
  }
}
