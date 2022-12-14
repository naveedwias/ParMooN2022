// time dependent Navier-Stokes problem 3D, ansatz
//
//u1(t,x,y,z) = std::exp(-t)*(y^2+z)
//u2(t,x,y,z) = t^4*(x-2*z^2)
//u3(t,x,y,z) = std::exp(-50*t)*(2*x*y+y)
//p(t,x,y,z) = t^3*(x+2*y+3*z)-3*t^3

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::root_info<1>("Example","Bsp3.h");
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = std::exp(-t)*(y*y+z);
values[1] = 0;
values[2] = 2*std::exp(-t)*y;
values[3] = std::exp(-t);
values[4] = 2*std::exp(-t);
}

void InitialU2(double x, double, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = t*t*t*t*(x-2*z*z);
values[1] = t*t*t*t;
values[2] = 0;
values[3] = -4*t*t*t*t*z;
values[4] = -4*t*t*t*t;
}

void InitialU3(double x, double y, double, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = std::exp(-50*t)*(2*x*y+y);
values[1] = 2*std::exp(-50*t)*y;
values[2] = std::exp(-50*t)*(2*x+1);
values[3] = 0;
values[4] = 0;
}

void InitialP(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = t*t*t*(x+2*y+3*z)-3*t*t*t;
values[1] = t*t*t;
values[2] = 2*t*t*t;
values[3] = 3*t*t*t;
values[4] = 0;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = std::exp(-t)*(y*y+z);
values[1] = 0;
values[2] = 2*std::exp(-t)*y;
values[3] = std::exp(-t);
values[4] = 2*std::exp(-t);
}

void ExactU2(double x, double, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = t*t*t*t*(x-2*z*z);
values[1] = t*t*t*t;
values[2] = 0;
values[3] = -4*t*t*t*t*z;
values[4] = -4*t*t*t*t;
}

void ExactU3(double x, double y, double, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = std::exp(-50*t)*(2*x*y+y);
values[1] = 2*std::exp(-50*t)*y;
values[2] = std::exp(-50*t)*(2*x+1);
values[3] = 0;
values[4] = 0;
}

void ExactP(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = t*t*t*(x+2*y+3*z)-3*t*t*t;
values[1] = t*t*t;
values[2] = 2*t*t*t;
values[3] = 3*t*t*t;
values[4] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
// kind of boundary condition (for FE space needed)
void BoundCondition(int, double, double, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(int, double, double y, double z, double &value)
{
double t=TDatabase::TimeDB->CURRENTTIME;
value = std::exp(-t)*(y*y+z);
}

void U2BoundValue(int, double x, double, double z, double &value)
{
double t=TDatabase::TimeDB->CURRENTTIME;
value = t*t*t*t*(x-2*z*z);
}

void U3BoundValue(int, double x, double y, double, double &value)
{
double t=TDatabase::TimeDB->CURRENTTIME;
value = std::exp(-50*t)*(2*x*y+y);
}

void U1BoundValue_diff(int, double, double y, double z, double &value)
{
double t=TDatabase::TimeDB->CURRENTTIME;
value = -std::exp(-t)*(y*y+z);
}

void U2BoundValue_diff(int, double x, double, double z, double &value)
{
double t=TDatabase::TimeDB->CURRENTTIME;
value = 4*t*t*t*(x-2*z*z);
}

void U3BoundValue_diff(int, double x, double y, double, double &value)
{
double t=TDatabase::TimeDB->CURRENTTIME;
value = -50*std::exp(-50*t)*(2*x*y+y);
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *X, const double *Y, const double *Z,
               const double *const*, double **coeffs)
{
  static double eps = DIMENSIONLESS_VISCOSITY;
  int i;
  double *coeff, x, y, z;
  double t=TDatabase::TimeDB->CURRENTTIME;
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    x = X[i];
    y = Y[i];
    z = Z[i];
    coeff[0] = eps;
    coeff[1] = -std::exp(-t)*(y*y+z)-2*eps*std::exp(-t)+2*t*t*t*t*(x-2*z*z)*std::exp(-t)*y+std::exp(-50*t)*(2*x*y+y)*std::exp(-t)+t*t*t;
    coeff[2] = 4*t*t*t*(x-2*z*z)+4*eps*t*t*t*t+std::exp(-t)*(y*y+z)*t*t*t*t-4*std::exp(-50*t)*(2*x*y+y)*t*t*t*t*z+2*t*t*t;
    coeff[3] = -50*std::exp(-50*t)*(2*x*y+y)+2*std::exp(-t)*(y*y+z)*std::exp(-50*t)*y+t*t*t*t*(x-2*z*z)*std::exp(-50*t)*(2*x+1)+3*t*t*t;
    coeff[4] = std::exp(-t)*(y*y+z)+2*eps*std::exp(-t)+8*t*t*t*(x-2*z*z)*std::exp(-t)*y-2*t*t*t*t*(x-2*z*z)*std::exp(-t)*y-51*std::exp(-50*t)*(2*x*y+y)*std::exp(-t)+3*t*t;
    coeff[5] = 12*t*t*(x-2*z*z)+16*eps*t*t*t-std::exp(-t)*(y*y+z)*t*t*t*t+4*std::exp(-t)*(y*y+z)*t*t*t+200*std::exp(-50*t)*(2*x*y+y)*t*t*t*t*z-16*std::exp(-50*t)*(2*x*y+y)*t*t*t*z+6*t*t;
    coeff[6] = 2500*std::exp(-50*t)*(2*x*y+y)-102*std::exp(-t)*(y*y+z)*std::exp(-50*t)*y+4*t*t*t*(x-2*z*z)*std::exp(-50*t)*(2*x+1)-50*t*t*t*t*(x-2*z*z)*std::exp(-50*t)*(2*x+1)+9*t*t;
  }
}
