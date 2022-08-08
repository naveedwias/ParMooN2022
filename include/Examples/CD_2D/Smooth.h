// ======================================================================
// smooth solution problem
// ======================================================================
#define __smooth_solution__

double DIFFUSION;

void ExampleFile()
{
  Output::print<1>("Example: Smooth.h");
  TDatabase::ParamDB->INTERNAL_QUAD_RULE = 99;
}

// exact solution
void Exact(double x, double y, double *values)
{
  values[0]=100.0*x*x*(1.0-x)*(1.0-x)*y*(1.0-2.0*y)*(1.0-y); //original solution
  values[1]= -200.0*(1.0-x)*x*x*y*(1.0-2.0*y)*(1.0-y)+200.0*(1-x)*(1-x)*x*y*(1.0-2.0*y)*(1.0-y); //x-derivative
  values[2]= 100.0*x*x*(1.0-x)*(1.0-x)*(1.0-2.0*y)*(1.0-y)-200.0*y*(1.0-y)*x*x*(1.0-x)*(1.0-x)-100.0*x*x*(1.0-x)*(1.0-x)*y*(1.0-2.0*y); //y-derivative
  values[3] = 200.0*(1.0-x)*(1.0-x)*y*(1.0-y)*(1.0-2*y)-800.0*x*(1.0-x)*y*(1.0-y)*(1.0-2.0*y)+200.0*x*x*y*(1.0-y)*(1.0-2.0*y)
             -200.0*x*x*(1.0-x)*(1.0-x)*(1.0-2.0*y)-400*x*x*(1.0-x)*(1.0-x)*(1.0-y)+400.0*x*x*(1.0-x)*(1.0-x)*y; //laplacian;
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

// kind of boundary condition (for FE space needed)
void BoundCondition_Poisson(int BdComp, double Param, BoundCond &cond)
{
  cond = NEUMANN;
}

// value of boundary condition
void BoundValue_Poisson(int, double, double &value)
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
  double eps=DIFFUSION;
  double *coeff;
  double exact[4];

  for(int i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    //param = parameters[i];

    coeff[0] = eps;
    coeff[1] = 3; //x-coordiante of convection b
    coeff[2] = 2; //y-coordiante of convection b
    coeff[3] = 1; //reaction coeffecient c
    
    Exact(X[i], Y[i], exact);
    coeff[4] = -coeff[0]*exact[3]; // diffusion
    coeff[4] += coeff[1]*exact[1] + coeff[2]*exact[2]; // convection
    coeff[4] += coeff[3]*exact[0]; // reaction
  }
}
