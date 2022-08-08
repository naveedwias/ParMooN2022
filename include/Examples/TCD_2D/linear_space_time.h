
void ExampleFile()
{
  Output::print("Example: linear_space_time.h");
}

double get_nu()
{
  return 1;
}

constexpr bool rhs_depends_on_time = true;
constexpr bool coefficients_depend_on_time = false;

// exact solution
void Exact(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 1+2*x+3*t*y;
  values[1] = 2;
  values[2] = 3*t;
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
  double t=TDatabase::TimeDB->CURRENTTIME;

  switch(BdComp)
  {
    case 0:
      value = 1+2*Param;
    break;
    case 1:
      value = 3+3*Param*t;
    break;
    case 2:
      value = 1+3*t+2*(1-Param);
    break;
    case 3:
      value = 1+3*t*(1-Param);
    break;
  } // endswitch
}

// initial conditon
void InitialCondition(double x,  double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 1+2*x+3*t*y;
}

void BilinearCoeffs(int n_points, const double *X, const double *Y,
                    const double *const*, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->RE_NR;
  double a=1, b=2, c=1;
  int i;
  double *coeff, x, y;
  double t=TDatabase::TimeDB->CURRENTTIME;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = a;
    coeff[2] = b;
    coeff[3] = c;

    coeff[4] = c*(1+2*x+3*t*y)
              +2*a
              +3*t*b
              +3*y;
  }
}

// exact solution
void Initial(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 1+2*x+3*t*y;
  values[1] = 2;
  values[2] = 3*t;
  values[3] = 0;
}

