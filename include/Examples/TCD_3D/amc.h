// #define __UREA__
// #define __SIMPATURS__
//   
// #include <Urea_3d4d.h>
// #include <MacroCell.h>

void ExampleFile()
{
  
  Output::print<1>("Example: amc.h");
}

// ========================================================================
// definitions for the temperature
// ========================================================================

void Exact(double x, double y, double z, double *values)
{
 double t = TDatabase::TimeDB->CURRENTTIME;
 double k = 0.1;
  
  values[0] = (std::exp(-k*t))*std::sin(Pi*x)*std::cos(Pi*y)*std::cos(Pi*z);
  values[1] = Pi*(std::exp(-k*t))*std::cos(Pi*x)*std::cos(Pi*y)*std::cos(Pi*z);
  values[2] = -Pi*(std::exp(-k*t))*std::sin(Pi*x)*std::sin(Pi*y)*std::cos(Pi*z);
  values[3] = -Pi*(std::exp(-k*t))*std::sin(Pi*x)*std::cos(Pi*y)*std::sin(Pi*z);
  values[4] = -3.*Pi*Pi*(std::exp(-k*t))*std::sin(Pi*x)*std::cos(Pi*y)*std::cos(Pi*z);
}

// initial conditon
void InitialCondition(double x, double y, double z, double *values)
{ 
 double t = 0;
 double k = 0.1;
  
  values[0] = (std::exp(-k*t))*std::sin(Pi*x)*std::cos(Pi*y)*std::cos(Pi*z);
 /* values[1] = Pi*(std::exp(-k*t))*std::cos(Pi*x)*std::cos(Pi*y)*std::cos(Pi*z);
  values[2] = -Pi*(std::exp(-k*t))*std::cos(Pi*x)*std::sin(Pi*y)*std::cos(Pi*z);
  values[3] = -Pi*(std::exp(-k*t))*std::cos(Pi*x)*std::cos(Pi*y)*std::sin(Pi*z);
  values[4] = -3.*Pi*Pi*(std::exp(-k*t))*std::sin(Pi*x)*std::cos(Pi*y)*std::cos(Pi*z);*/  
}

void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{
 double t = TDatabase::TimeDB->CURRENTTIME;
 double k = 0.1;
  
  value = (std::exp(-k*t))*std::sin(Pi*x)*std::cos(Pi*y)*std::cos(Pi*z);  
}

// ========================================================================
// BilinearCoeffs for Heat 
// ========================================================================
void BilinearCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
//   double eps= 1.0/TDatabase::ParamDB->RE_NR;
  double eps= 1.0;
  int i;
  double *coeff;                                  // *param;
  double x, y, z, c, a[3], b[3], s[3], h;
  double t = TDatabase::TimeDB->CURRENTTIME;
  
  b[0] = 0;
  b[1] = 0;
  b[2] = 0;  
  c = 0;

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
    coeff[1] = b[0];
    // convection in y direction
    coeff[2] = b[1];
    // convection in z direction
    coeff[3] = b[2];
    // reaction
    coeff[4] = c;
     // rhs
    coeff[5] = (3.*eps*Pi*Pi - 0.1)*(std::exp(-0.1*t))*std::sin(Pi*x)*std::cos(Pi*y)*std::cos(Pi*z);
    coeff[6] = 0;

  }
}
