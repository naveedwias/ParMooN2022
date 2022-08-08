#include "MooNMD_Io.h"

void ExampleFile()
{
  Output::print("Example: PolySin.h");
}
// =====================================================================
// exact solution
// =====================================================================
double pi = 3.14159265358979;
// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;

void ExactU1(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  
  double u1 = x*x*(1.-x)*(1.-x)*(2.*y*(1.-y)*(1.-y)
  -2.*y*y*(1.-y))*std::sin(10.*pi*t);
  
  double u1_x = 2.*x*(1.-x)*(1.-x)*(2.*y*(1.-y)*(1.-y)
    -2.*y*y*(1.-y))*std::sin(10.*pi*t)
    -2.*x*x*(1.-x)*(2.*y*(1.-y)*(1.-y)
    -2.*y*y*(1.-y))*std::sin(10.*pi*t);
  
  double u1_y = x*x*(1.-x)*(1.-x)*(2.*(1.-y)*(1.-y)
    -8*y*(1.-y)+2.*y*y)*std::sin(10.*pi*t);
    
  double u1_xx = 2.*(1.-x)*(1.-x)*(2.*y*(1.-y)*(1.-y) -2.*y*y*(1.-y))*std::sin(10.*pi*t)
            - 8*x*(1.-x)*(2.*y*(1.-y)*(1.-y)-2.*y*y*(1.-y))*std::sin(10.*pi*t)
            + 2.*x*x*(2.*y*(1.-y)*(1.-y)-2.*y*y*(1.-y))*std::sin(10.*pi*t);

  double u1_yy = x*x*(1.-x)*(1.-x)*(-12.+24.*y)*std::sin(10.*pi*t);
  
  values[0] = u1;
  values[1] = u1_x;
  values[2] = u1_y;
  values[3] = u1_xx + u1_yy; 
}

void ExactU2(double x, double y, double *values)
{
  
  double t =TDatabase::TimeDB->CURRENTTIME;
  
  double u2   = -y*y*(1.-y)*(1.-y)*(2.*x*(1.-x)*(1.-x) 
                -2.*x*x*(1.-x))*std::sin(10.*pi*t);
  double u2_x = -y*y*(1.-y)*(1.-y)*(2.*(1.-x)*(1.-x)-8*x*(1.-x)
                + 2.*x*x)*std::sin(10.*pi*t);
  double u2_y = -2.*y*(1.-y)*(1.-y)*(2.*x*(1.-x)*(1.-x)
                -2.*x*x*(1.-x))*std::sin(10.*pi*t) 
                +2.*y*y*(1.-y)*(2.*x*(1.-x)*(1.-x)
                -2.*x*x*(1.-x))*std::sin(10.*pi*t);

  double u2_xx = -y*y*(1.-y)*(1.-y)*(-12.+24.*x)*std::sin(10.*pi*t);
  
  double u2_yy = -2.*(1.-y)*(1.-y)*(2.*x*(1.-x)*(1.-x)-2.*x*x*(1.-x))*std::sin(10.*pi*t)
                 +8*y*(1.-y)*(2.*x*(1.-x)*(1.-x)-2.*x*x*(1.-x))*std::sin(10.*pi*t)
                 - 2.*y*y*(2.*x*(1.-x)*(1.-x)-2.*x*x*(1.-x))*std::sin(10.*pi*t);
  values[0] = u2;
  values[1] = u2_x;
  values[2] = u2_y;
  values[3] = u2_xx + u2_yy; 
}
void ExactP(double x, double y, double *values)
{
  double t, p, p_x, p_y;
  
  t = TDatabase::TimeDB->CURRENTTIME;
  
  p = -(x*x*x + y*y*y -0.5)*(1.5 + 0.5*std::sin(10.*pi*t));
  p_x = -3.*x*x*(1.5+.5*std::sin(10.*pi*t));
  p_y = -3.*y*y*(1.5+.5*std::sin(10.*pi*t));
  
  values[0] = p;
  values[1] = p_x;
  values[2] = p_y;
  values[3] = 0;
  
}
// =====================================================================
// initial solution
// =====================================================================

void InitialU1(double, double, double *values)
{
  values[0] = 0;
}

void InitialU2(double, double, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y, double *values)
{
   double t = TDatabase::TimeDB->CURRENTTIME;
   values[0] = -(x*x*x + y*y*y -0.5)*(1.5+0.5*std::sin(10.*pi*t));
}

// =====================================================================
// boundary conditions
// =====================================================================
void BoundCondition(int, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

void U1BoundValue(int, double, double &value)
{
  value = 0;
}

void U2BoundValue(int, double, double &value)
{
  value = 0;
}
// =====================================================================
// coefficients 
// =====================================================================
void LinCoeffs(int n_points, const double *X, const double *Y, 
               const double *const*, double **coeffs)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  double nu = DIMENSIONLESS_VISCOSITY;
  
  double u1[4], u2[4], p[4];
  
  for(int i=0;i<n_points;i++)
  {
    double x = X[i];
    double y = Y[i];

    double u1_t = 10.*x*x*(1.-x)*(1.-x)*(2.*y*(1.-y)*(1.-y)
             -2.*y*y*(1.-y))*std::cos(10.*pi*t)*pi;
    double u2_t = -10.*y*y*(1.-y)*(1.-y)*(2.*x*(1.-x)*(1.-x)
            -2.*x*x*(1.-x))*std::cos(10.*pi*t)*pi;
    ExactU1(x,y,u1);
    ExactU2(x,y,u2);
    ExactP(x,y,p);
    
    coeffs[i][0] = nu;
    coeffs[i][1] = u1_t -nu*u1[3] + u1[0]*u1[1] + u2[0]*u1[2] + p[1];
    coeffs[i][2] = u2_t -nu*u2[3] + u1[0]*u2[1] + u2[0]*u2[2] + p[2]; 
  }  
}
