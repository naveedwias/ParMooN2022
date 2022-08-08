// Navier-Stokes problem.
// 
// The velocity part of the solution is chosen to be such that 
//   1. the nonlinear term (in convective form) vanishes
//   2. the solution is linear in space and time.
// Therefore using a fully explicit temporal
// discretization for this term is sufficient, i.e., the computed solution 
// should be found exactly. Of course also imex and fully implicit 
// discretizations work.
// If you set the following boolean parameter to true, a fully explicit temporal
// discretization should no longer work, but the imex scheme should.

bool test_no_fully_explicit = false;
double DIMENSIONLESS_VISCOSITY = 1.;

void ExampleFile()
{
  Output::print("Example: linear_space_time.h. The velocity solution is a "
                "first order polynomial in space and time. The pressure is "
                "zero. The nonlinear term (convective form) vanishes for this "
                "solution.");
}

void ExactU1(double x, double y, double *values);
void ExactU2(double x, double y, double *values);
void ExactP(double x, double y, double *values);

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  double val[4];
  ExactU1(x, y, val);
  values[0] = val[0];
}

void InitialU2(double x, double y, double *values)
{
  double val[4];
  ExactU2(x, y, val);
  values[0] = val[0];
}
void InitialP(double x, double y, double *values)
{
  double val[4];
  ExactP(x, y, val);
  values[0] = val[0];
}


// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 3*t*(x + y) + (test_no_fully_explicit ? 1. : 0.);
  values[1] = 3*t;
  values[2] = 3*t;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  // U2 = -U1:
  ExactU1(x, y, values);
  values[0] = -1. * values[0] + (test_no_fully_explicit ? 1. : 0.);
  for(int i = 1; i < 4; ++i)
    values[i] = -1. * values[i];
}

void ExactP(double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

void get_xy_of_boundary_parametrization(int boundary_component, double Param,
                                        double&x, double& y)
{
  switch(boundary_component)
  {
    case 0:
      x = Param;
      y = 0;
    break;
    case 1:
      x = 1;
      y = Param;
    break;
    case 2:
      x = 1-Param;
      y = 1;
    break;
    case 3:
      x = 1-Param;
      y = 0;
    break;
    default:
      ErrThrow("wrong boundary component ", boundary_component);
  }
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  double x, y;
  get_xy_of_boundary_parametrization(BdComp, Param, x, y);
  double val[4];
  ExactU1(x, y, val);
  value = val[0];
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  double x, y;
  get_xy_of_boundary_parametrization(BdComp, Param, x, y);
  double val[4];
  ExactU2(x, y, val);
  value = val[0];
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *X, const double *Y,
               const double *const*, double **coeffs)
{
  double u1[4];
  double u2[4];
  double p[4];
  for(int i = 0; i < n_points; ++i)
  {
    double x = X[i];
    double y = Y[i];
    ExactU1(x, y, u1);
    ExactU2(x, y, u2);
    ExactP(x, y, p);
    coeffs[i][0] = DIMENSIONLESS_VISCOSITY;
    coeffs[i][1] = -DIMENSIONLESS_VISCOSITY*u1[3] + p[1];
    coeffs[i][2] = -DIMENSIONLESS_VISCOSITY*u2[3] + p[2];
    if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE != STOKES)
    {
      coeffs[i][1] += u1[0]*u1[1] + u2[0]*u1[2];
      coeffs[i][2] += u1[0]*u2[1] + u2[0]*u2[2];
    }
    // unfortunately the time derivative has to be done by hand:
    coeffs[i][1] += 3*(x+y);
    coeffs[i][2] += -3*(x+y);
    // divergence
    coeffs[i][3] = u1[1] + u2[2];
  }
}
