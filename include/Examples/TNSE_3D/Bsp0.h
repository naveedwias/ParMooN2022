// time dependent Navier-Stokes problem 3D, ansatz
//
//u1(t,x,y,z) = 2*t
//u2(t,x,y,z) = 3
//u3(t,x,y,z) = 4
//p(t,x,y,z) = 0

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::root_info<1>("Example", "Bsp0.h");
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double, double, double, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = 2*t;
values[1] = 0;
values[2] = 0;
values[3] = 0;
values[4] = 0;
}

void InitialU2(double, double, double, double *values)
{
values[0] = 3;
values[1] = 0;
values[2] = 0;
values[3] = 0;
values[4] = 0;
}

void InitialU3(double, double, double, double *values)
{
values[0] = 4;
values[1] = 0;
values[2] = 0;
values[3] = 0;
values[4] = 0;
}

void InitialP(double, double, double, double *values)
{
values[0] = 0;
values[1] = 0;
values[2] = 0;
values[3] = 0;
values[4] = 0;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double, double, double, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = 2*t;
values[1] = 0;
values[2] = 0;
values[3] = 0;
values[4] = 0;
}

void ExactU2(double, double, double, double *values)
{
values[0] = 3;
values[1] = 0;
values[2] = 0;
values[3] = 0;
values[4] = 0;
}

void ExactU3(double, double, double, double *values)
{
values[0] = 4;
values[1] = 0;
values[2] = 0;
values[3] = 0;
values[4] = 0;
}

void ExactP(double, double, double, double *values)
{
values[0] = 0;
values[1] = 0;
values[2] = 0;
values[3] = 0;
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
void U1BoundValue(int, double, double, double, double &value)
{
double t=TDatabase::TimeDB->CURRENTTIME;
value = 2*t;
}

void U2BoundValue(int, double, double, double, double &value)
{
value = 3;
}

void U3BoundValue(int, double, double, double, double &value)
{
value = 4;
}

void U1BoundValue_diff(int, double, double, double, double &value)
{
value = 2;
}

void U2BoundValue_diff(int, double, double, double, double &value)
{
value = 0;
}

void U3BoundValue_diff(int, double, double, double, double &value)
{
value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *, const double *, const double *,
               const double *const*, double **coeffs)
{
  const double eps = DIMENSIONLESS_VISCOSITY;
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = eps;
    coeffs[i][1] = 2;
    coeffs[i][2] = 0;
    coeffs[i][3] = 0;
    coeffs[i][4] = 0;
    coeffs[i][5] = 0;
    coeffs[i][6] = 0;
  }
}
