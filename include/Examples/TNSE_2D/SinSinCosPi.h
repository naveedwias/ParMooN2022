// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::print("Example: SinCosExp.h   viscosity: ", DIMENSIONLESS_VISCOSITY);
}

// ========================================================================
// initial solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  double Pi = 3.14159265358979;
  double scale = TDatabase::ParamDB->P5;
  t *= Pi;
 
  values[0] = 2*Pi*sin(t)*sin(Pi*x)*sin(Pi*x)*sin(Pi*y)*cos(Pi*y);
  values[1] = 4*Pi*Pi*sin(t)*sin(Pi*x)*cos(Pi*x)*sin(Pi*y)*cos(Pi*y);
  values[2] = 2*Pi*Pi*sin(t)*sin(Pi*x)*sin(Pi*x)*(cos(Pi*y)*cos(Pi*y)-sin(Pi*y)*sin(Pi*y));
  values[3] = 4*Pi*Pi*Pi*sin(t)*(cos(Pi*x)*cos(Pi*x)
              -sin(Pi*x)*sin(Pi*x))*sin(Pi*y)*cos(Pi*y)
              - 8*Pi*Pi*Pi*sin(t)*sin(Pi*x)*sin(Pi*x)*cos(Pi*y)*sin(Pi*y);
  
  values[0] *= scale;
  values[1] *= scale;
  values[2] *= scale;
  values[3] *= scale;
}

void ExactU2(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  double scale = TDatabase::ParamDB->P5;
  double Pi = 3.14159265358979;
  t*=Pi;

  values[0] = -2*Pi*sin(t)*sin(Pi*y)*sin(Pi*y)*sin(Pi*x)*cos(Pi*x);
  values[1] = -2*Pi*Pi*sin(t)*sin(Pi*y)*sin(Pi*y)*(cos(Pi*x)*cos(Pi*x)-sin(Pi*x)*sin(Pi*x));
  values[2] = -4*Pi*Pi*sin(t)*sin(Pi*y)*cos(Pi*y)*sin(Pi*x)*cos(Pi*x);
  values[3] = 8*Pi*Pi*Pi*sin(t)*sin(Pi*y)*sin(Pi*y)*cos(Pi*x)*sin(Pi*x)
  -4*Pi*Pi*Pi*sin(t)*(cos(Pi*y)*cos(Pi*y)-sin(Pi*y)*sin(Pi*y))*sin(Pi*x)*cos(Pi*x);
  
  values[0] *= scale;
  values[1] *= scale;
  values[2] *= scale;
  values[3] *= scale;
}

void ExactP(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  double Pi = 3.14159265358979;
  t *=Pi;
 
  values[0] = 20*sin(t)*(x*x*y-1.0/6.0);
  values[1] = 40*sin(t)*x*y;
  values[2] = 20*sin(t)*x*x;
  values[3] = 0;
}

// ========================================================================
// exact solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  double val[4];
  
  ExactU1(x,y,val);
  values[0] = val[0];
}

void InitialU2(double x, double y, double *values)
{
  double val[4];
  
  ExactU2(x,y,val);
  values[0] = val[0];
}

void InitialP(double x, double y, double *values)
{
  double val[4];
  
  ExactP(x,y,val);
  values[0] = val[0];
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
  double val[4], x, y;
  
  switch(BdComp)
  {
    case 0:
      x = Param; y = 0;
      break;
    case 1:
      x = 1; y = Param;
      break;
    case 2:
      x = 1-Param; y = 1;
      break;
    case 3:
      x = 0; y = 1-Param;
      break;
  }
  
  ExactU1(x,y,val);
  value = val[0]; 
}


void U1BoundValue_diff(int, double, double &value)
{
   value = 0.;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  double val[4], x, y;
  
  switch(BdComp)
  {
    case 0:
      x = Param; y = 0;
      break;
    case 1:
      x = 1; y = Param;
      break;
    case 2:
      x = 1-Param; y = 1;
      break;
    case 3:
      x = 0; y = 1-Param;
      break;
  }
  
  ExactU2(x,y,val);
  value = val[0]; 
}
void U2BoundValue_diff(int, double, double &value)
{
  value = 0.;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *X, const double *Y,
               const double *const*, double **coeffs)
{
  int i;
  double *coeff, x, y, values[15]; 
  double nu=DIMENSIONLESS_VISCOSITY;
  double alpha = TDatabase::ParamDB->OSEEN_ZERO_ORDER_COEFF;
  double t=TDatabase::TimeDB->CURRENTTIME;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    coeff[0] = nu;
    x = X[i];
    y = Y[i];
    ExactU1(x, y, values);
    ExactU2(x, y, values+5);
    ExactP (x, y, values+10);
    

    constexpr double Pi = 3.14159265358979;
    double dt_u1 = 2*Pi*Pi*cos(t)*sin(Pi*x)*sin(Pi*x)*sin(Pi*y)*cos(Pi*y);
    double dt_u2 = -2*Pi*Pi*cos(t)*sin(Pi*y)*sin(Pi*y)*sin(Pi*x)*cos(Pi*x);
    coeff[1] = dt_u1 -coeff[0] * values[3]  +  values[0] * values[1] + values[5] * values[2]  + alpha * values[0] +  values[11];
    coeff[2] = dt_u2 -coeff[0] * values[8] +  values[0] * values[6] + values[5] * values[7] + alpha * values[5] +  values[12];
    
    coeff[3] = 0; // divergence
    coeff[4] = alpha; // Brinkman term
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
