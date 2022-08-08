// Time-dependent 3D Navier-Stokes problem
//
// The velocity part of the solution is chosen to be such that 
//   1. the nonlinear term (in convective form) vanishes
//   2. the solution is linear in space and time.
// Therefore using a fully explicit temporal
// discretization for this term is sufficient, i.e., the computed solution 
// should be found exactly. Of course also imex and full implicit 
// discretizations work.
// If you set the following boolean parameter to true, a fully explicit temporal
// discretization should no longer work, but the imex scheme should.


// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::root_info<1>("Example","linear_space_time.h");
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double z, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = 6*t*(x + y + z);
  values[1] = 6*t;
  values[2] = 6*t;
  values[3] = 6*t;
  values[4] = 0;
}

void ExactU2(double x, double y, double z, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = -3*t*(x + y + z);
  values[1] = -3*t;
  values[2] = -3*t;
  values[3] = -3*t;
  values[4] = 0;
}

void ExactU3(double x, double y, double z, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = -3*t*(x + y + z);
  values[1] = -3*t;
  values[2] = -3*t;
  values[3] = -3*t;
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
// initial conditions
// ========================================================================
void InitialU1(double x, double y, double z, double *values)
{
  double val[5];
  ExactU1(x, y, z, val);
  values[0] = val[0];
}

void InitialU2(double x, double y, double z, double *values)
{
  double val[5];
  ExactU2(x, y, z, val);
  values[0] = val[0];
}

void InitialU3(double x, double y, double z, double *values)
{
  double val[5];
  ExactU3(x, y, z, val);
  values[0] = val[0];
}

void InitialP(double, double, double, double *values)
{
  values[0] = 0;
}


// ========================================================================
// kind of boundary condition (for FE space needed) and values
// ========================================================================
void BoundCondition(int, double, double, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(int, double x, double y, double z, double &value)
{
  double val[5];
  ExactU1(x, y, z, val);
  value = val[0];
}

// value of boundary condition
void U2BoundValue(int, double x, double y, double z, double &value)
{
  double val[5];
  ExactU2(x, y, z, val);
  value = val[0];
}

// value of boundary condition
void U3BoundValue(int, double x, double y, double z, double &value)
{
  double val[5];
  ExactU3(x, y, z, val);
  value = val[0];
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *X, const double *Y, const double *Z,
               const double *const*, double **coeffs)
{
  double u1[5];
  double u2[5];
  double u3[5];
  double p[5];
  for(int i = 0; i < n_points; i++)
  {
    double x = X[i];
    double y = Y[i];
    double z = Z[i];
    ExactU1(x, y, z, u1);
    ExactU2(x, y, z, u2);
    ExactU3(x, y, z, u3);
    ExactP(x, y, z, p);
    coeffs[i][0] = DIMENSIONLESS_VISCOSITY;
    coeffs[i][1] = -DIMENSIONLESS_VISCOSITY*u1[4] + p[1];
    coeffs[i][2] = -DIMENSIONLESS_VISCOSITY*u2[4] + p[2];
    coeffs[i][3] = -DIMENSIONLESS_VISCOSITY*u3[4] + p[3];
    if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE != STOKES)
    {
      coeffs[i][1] += u1[0]*u1[1] + u2[0]*u1[2] + u3[0]*u1[3];
      coeffs[i][2] += u1[0]*u2[1] + u2[0]*u2[2] + u3[0]*u2[3];
      coeffs[i][3] += u1[0]*u3[1] + u2[0]*u3[2] + u3[0]*u3[3];
    }
    // unfortunately the time derivative has to be done by hand:
    coeffs[i][1] += 6*(x + y + z);
    coeffs[i][2] += -3*(x + y + z);
    coeffs[i][3] += -3*(x + y + z);
    // divergence
    coeffs[i][4] = u1[1] + u2[2] + u3[3];
  }
}
