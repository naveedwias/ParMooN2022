#define __time_domianted__

void ExampleFile()
{
  Output::print<1>("Example: time_dominated.h\n") ;
}

constexpr bool rhs_depends_on_time = true;
constexpr bool coefficients_depend_on_time = false;

void Exact(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  
  values[0] = std::sin(50*t)*x*(1-x)*y*(1-y);  //original solution
  values[1] = y*(1-y)*std::sin(50*t)*(1-2*x);  //x-derivative
  values[2] = x*(1-x)*std::sin(50*t)*(1-2*y);  //y-derivative
  values[3] = -2*y*(1-y)*std::sin(50*t)-2*x*(1-x)*std::sin(50*t); //Laplacian
}

void BoundCondition(int, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

// initial conditon
void InitialCondition(double , double , double *values)
{
  values[0] = 0;
}

// value of boundary condition
void BoundValue(int, double, double &value)
{
  value = 0;
}

void BilinearCoeffs(int n_points, const double *X, const double *Y,
                    const double *const*, double **coeffs)
{
  double eps = 1/TDatabase::ParamDB->RE_NR;
  double b1=3, b2=2, c=1;
  double t=TDatabase::TimeDB->CURRENTTIME;
  double *coeff;
  double x,y;

  for(int i=0;i<n_points;i++)
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
    coeff[4] += 50*std::cos(50*t)*x*(1-x)*y*(1-y); //time_derivative
  }
}
