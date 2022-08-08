// Navier-Stokes problem, Driven cavity
// 
// u(x,y) = unknown
// p(x,y) = unknown

//WIll be manipulated by the example class.
double DIMENSIONLESS_VISCOSITY = 10;

void ExampleFile()
{
  Output::print<1>("Example: DrivenCavity.h");
}

// ========================================================================
// exact solution
// ========================================================================
auto& ExactU1 = unknown_solution_2d;
auto& ExactU2 = unknown_solution_2d;
auto& ExactP = unknown_solution_2d;

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
    case 0: value = 0;
            break;
    case 1: value = 0;
            break;
    case 2: if(Param<0.00001 || Param>0.99999) 
              value = 0;
            else
              value = 1;
            break;
    case 3: value = 0;
            break;
    default: cout << "wrong boundary part number" << endl;
  }
}

void U2BoundValue(int BdComp, double, double &value)
{
  value = 0;
  if(BdComp>3) cout << "wrong boundary part number" << endl;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *, const double *,
               const double *const*, double **coeffs)
{
  double eps = DIMENSIONLESS_VISCOSITY;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
    coeff[3] = 0; // divergence
    // additional coefficient (used only in the Brinkman problem)
    coeff[4] = 0.;
  }
}

