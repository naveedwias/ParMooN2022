// ======================================================================
// Known solution with a boundary layer.
// Example 1 from ABR17.SISC (Allendes, Barrenchea, and Rankin)
// ======================================================================
#define __example_1_ABR17__

double DIFFUSION;

void ExampleFile()
{
  Output::print<1>("Example: example_1_ABR17.h");
  TDatabase::ParamDB->INTERNAL_QUAD_RULE = 99;
}

// exact solution
void Exact(double x, double y, double *values)
{
  double exp_x, exp_eps, u_0_x;
  double eps = DIFFUSION;
  exp_x = std::exp((x-1)/eps);
  exp_eps = 1./(1-std::exp(-1./eps));
  u_0_x=x-exp_eps*(exp_x-std::exp(-1./eps));
  
  
  values[0] = y*(1-y)*u_0_x; //original solution
  values[1] = y*(1-y)*(1-(exp_x*exp_eps)/eps); //x-derivative
  values[2] = u_0_x*(1-y)-u_0_x*y; //y-derivative
  values[3] = -y*(1-y)*(exp_x*exp_eps)/(eps*eps)-2*u_0_x; //laplacian;
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
   Exact(x,y,values);
}


void BilinearCoeffs(int n_points, const double *X, const double *Y,
                    const double *const*, double **coeffs)
{
  double eps = DIFFUSION;
  double *coeff;
  double exact[4];

  for(int i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    //param = parameters[i];

    coeff[0] = eps;
    coeff[1] = 2; //x-coordiante of convection b
    coeff[2] = 1; //y-coordiante of convection b
    coeff[3] = 1; //reaction coeffecient c
    
    Exact(X[i], Y[i], exact);
    coeff[4] = -coeff[0]*exact[3]; // diffusion
    coeff[4] += coeff[1]*exact[1] + coeff[2]*exact[2]; // convection
    coeff[4] += coeff[3]*exact[0]; // reaction
  }
}
