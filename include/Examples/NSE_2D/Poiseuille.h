// Navier-Stokes problem, Poiseuille-Problem
// 
// u(x,y) = (4*y*(1-y), 0)
// p(x,y) = x-1/2

void ExampleFile()
{
  Output::print<1>("Example: Poiseuille.h");
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double, double y, double *values)
{
  values[0] = 4*y*(1-y);
  values[1] = 0;
  values[2] = 4-8*y;
  values[3] = -8;
}

void ExactU2(double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactP(double x, double, double *values)
{
  values[0] = x-0.5;
  values[1] = 1;
  values[2] = 0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=4*Param*(1-Param);
            break;
    case 2: value=0;
            break;
    case 3: value=4*Param*(1-Param);
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void U2BoundValue(int BdComp, double, double &value)
{
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=0;
            break;
    case 2: value=0;
            break;
    case 3: value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *, const double *,
               const double *const*, double **coeffs)
{
  static double eps = 1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 1+8*eps; // f1
    coeff[2] = 0; // f2
    coeff[3] = 0.;//g
    // additional coefficient (used only in the Brinkman problem)
    coeff[4] = 0.;
  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void NonLinCoeffs(int n_points, const double *, const double *,
                  const double *const*parameters, double **coeffs)
{
  static double eps = 1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;
  const double *param;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = eps;
    coeff[1] = param[0];
    coeff[2] = param[1];

  }
}

