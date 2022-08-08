// time dependent Navier-Stokes problem 
//
// Tafti, Computers & Fluids (25), p. 647 - 665, 1996
// 
// u_1(x,y) = - std::cos(n pi x) std::sin(n pi y) std::exp(-2 (n pi)^2 nu t)
// u_2(x,y) = std::sin(n pi x) std::cos(n pi y) std::exp(-2 (n pi)^2 nu t)
// p(x,y) = -1/4 (std::cos(2n Pi x) + std::cos(2n Pi y))  std::exp(-4 (n pi)^2 nu t)
//
// n - number of oscillations

#define N_OSCILLATIONS    4

void ExampleFile()
{
  Output::print("Example: SinCosExp.h Number of oscillations: ",
                N_OSCILLATIONS);
}

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  double pi = 3.14159265358979;
  double nu=1/TDatabase::ParamDB->RE_NR;
 
  fac = std::exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
  values[0] = -std::cos(N_OSCILLATIONS*pi*x)*std::sin(N_OSCILLATIONS*pi*y)*fac;
  //Output::print(values[0]);
  //values[0] = -std::cos(N_OSCILLATIONS*pi*x)*std::sin(N_OSCILLATIONS*pi*y);
}

void InitialU2(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  double pi = 3.14159265358979;
  double nu=1/TDatabase::ParamDB->RE_NR;
 
  fac = std::exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
  values[0] = std::sin(N_OSCILLATIONS*pi*x)*std::cos(N_OSCILLATIONS*pi*y)*fac;
}

void InitialP(double x, double y, double *values)
{
  // values[0] = -(std::cos(2*N_OSCILLATIONS*pi*x)+std::cos(2*N_OSCILLATIONS*pi*y))/4.0;
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  double pi = 3.14159265358979;
  double nu=1/TDatabase::ParamDB->RE_NR;
 
  fac = std::exp(-4*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
  values[0] = -(std::cos(2*N_OSCILLATIONS*pi*x)+std::cos(2*N_OSCILLATIONS*pi*y))*fac/4.0;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  double pi = 3.14159265358979;
  double nu=1/TDatabase::ParamDB->RE_NR;
 
  fac = std::exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
  values[0] = -std::cos(N_OSCILLATIONS*pi*x)*std::sin(N_OSCILLATIONS*pi*y)*fac;
  values[1] = N_OSCILLATIONS*pi*std::sin(N_OSCILLATIONS*pi*x)*std::sin(N_OSCILLATIONS*pi*y)*fac;
  values[2] = -N_OSCILLATIONS*pi*std::cos(N_OSCILLATIONS*pi*x)*std::cos(N_OSCILLATIONS*pi*y)*fac;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  double pi = 3.14159265358979;
  double nu=1/TDatabase::ParamDB->RE_NR;
 
  fac = std::exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
  values[0] = std::sin(N_OSCILLATIONS*pi*x)*std::cos(N_OSCILLATIONS*pi*y)*fac;
  values[1] = N_OSCILLATIONS*pi*std::cos(N_OSCILLATIONS*pi*x)*std::cos(N_OSCILLATIONS*pi*y)*fac;
  values[2] = -N_OSCILLATIONS*pi* std::sin(N_OSCILLATIONS*pi*x)*std::sin(N_OSCILLATIONS*pi*y)*fac;;
  values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  double pi = 3.14159265358979;
  double nu=1/TDatabase::ParamDB->RE_NR;
 
  fac = std::exp(-4*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
  values[0] = -(std::cos(2*N_OSCILLATIONS*pi*x)+std::cos(2*N_OSCILLATIONS*pi*y))*fac/4.0;
  values[1] = 2*N_OSCILLATIONS*pi*std::sin(2*N_OSCILLATIONS*pi*x)*fac/4.0;
  values[2] = 2*N_OSCILLATIONS*pi*std::sin(2*N_OSCILLATIONS*pi*y)*fac/4.0;
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
  double t = TDatabase::TimeDB->CURRENTTIME;
  double nu=1/TDatabase::ParamDB->RE_NR;
  double fac;
  double pi = 3.14159265358979;

  switch(BdComp)
  {
  case 0: 
  case 2:
    value=0;
    break;
    case 1: 
      fac = std::exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
      value = -std::cos(N_OSCILLATIONS*pi)*std::sin(N_OSCILLATIONS*pi*Param)*fac;
      break;
    case 3: 
      fac = std::exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
      value = -std::sin(N_OSCILLATIONS*pi*(1-Param))*fac;
      break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void U1BoundValue_diff(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  double nu=1/TDatabase::ParamDB->RE_NR;
  double fac;
  double pi = 3.14159265358979;

  switch(BdComp)
  {
  case 0: 
  case 2:
    value=0;
    break;
    case 1: 
      fac = -2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*std::exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
      value = -std::cos(N_OSCILLATIONS*pi)*std::sin(N_OSCILLATIONS*pi*Param)*fac;
      break;
    case 3: 
      fac = -2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*std::exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
      value = -std::sin(N_OSCILLATIONS*pi*(1-Param))*fac;
      break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  double nu=1/TDatabase::ParamDB->RE_NR;
  double fac;
  double pi = 3.14159265358979;
  value = 0;
  switch(BdComp)
  {
  case 1: 
  case 3:
    value=0;
    break;
  case 0: 
    fac = std::exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
    value = std::sin(N_OSCILLATIONS*pi*Param)*fac;
    break;
  case 2: 
    fac = std::exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
    value = std::sin(N_OSCILLATIONS*pi*(1-Param))*std::cos(N_OSCILLATIONS*pi)*fac;
    break;
  default: cout << "wrong boundary part number" << endl;
    break;
  }
}
void U2BoundValue_diff(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  double nu=1/TDatabase::ParamDB->RE_NR;
  double fac;
  double pi = 3.14159265358979;
  value = 0;
  switch(BdComp)
  {
  case 1: 
  case 3:
    value=0;
    break;
  case 0: 
    fac = -2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*std::exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
    value = std::sin(N_OSCILLATIONS*pi*Param)*fac;
    break;
  case 2: 
    fac = -2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*std::exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
    value = std::sin(N_OSCILLATIONS*pi*(1-Param))*std::cos(N_OSCILLATIONS*pi)*fac;
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
  int i;
  double *coeff;//, x, y;
  double eps = DIMENSIONLESS_VISCOSITY;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
//    x = X[i];
//    y = Y[i];

    coeff[0] = eps;

    coeff[1] = 0;
    coeff[2] = 0;
  }
}

void BoundConditionPressure(int, double, BoundCond &cond)
{
     cond = NEUMANN;
}

void BoundConditionPressureLaplace(int bdcomp, double, BoundCond &cond)
{
  switch(bdcomp)
  {
    case 1:
      cond = DIRICHLET;
      break;
    default:
      cond = NEUMANN;
      break;
  }      
}

void PressureBoundValue(int BdComp, double, double &value)
{
  switch(BdComp)
  {
    case 0: value =0;
            break;
    case 1: 
            value = 0;					
            break;
    case 2: value = 0;
            break;
    case 3: 
	    value = 0;
            break;
    default: cout << "wrong boundary part number: " << BdComp << endl;
  }
}

void PressureBoundValueLaplace(int BdComp, double, double &value)
{
  switch(BdComp)
  {
    case 0: value =0;
            break;
    case 1: 
            value = 0;					
            break;
    case 2: value = 0;
            break;
    case 3: 
	    value = 0;
            break;
    default: cout << "wrong boundary part number: " << BdComp << endl;
  }
}

void EvaluateSolution(TFEFunction2D **, TFEVectFunct2D  **, double *, int *)
{
  // so special evaluations in this example
  return;
}
