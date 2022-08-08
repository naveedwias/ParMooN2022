// Navier-Stokes problem, Poiseuille-Problem
// 
// u(x,y) = (std::sin(t)*std::sin(Pi*x)*std::sin(Pi*y), std::sin(t)*std::cos(Pi*x)*std::cos(Pi*y))
// p(x,y) = std::sin(t)*(std::sin(Pi*x)+std::cos(Pi*y)-2/Pi)

void ExampleFile()
{
  Output::print("Example: Bsp1.h");
  TDatabase::ParamDB->INPUT_QUAD_RULE = 99;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = std::sin(t)*std::sin(M_PI*x)*std::sin(M_PI*y);
}

void InitialU2(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = std::sin(t)*std::cos(M_PI*x)*std::cos(M_PI*y);
}

void InitialP(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = std::sin(t)*(std::sin(M_PI*x)+std::cos(M_PI*y)-2/M_PI);
}


// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double t1, t2, t4, t5, t6, t15;
  double t=TDatabase::TimeDB->CURRENTTIME;

  t1 = std::sin(t);
  t2 = 0.3141592653589793E1*x;
  t4 = t1*std::sin(t2);
  t5 = 0.3141592653589793E1*y;
  t6 = std::sin(t5);
  t15 = 0.3141592653589793E1*0.3141592653589793E1;

  values[0] = t4*t6;
  values[1] = t1*std::cos(t2)*0.3141592653589793E1*t6;
  values[2] = t4*std::cos(t5)*0.3141592653589793E1;
  values[3] = -2.0*t4*t15*t6;
}

void ExactU2(double x, double y, double *values)
{
  double t1, t2, t4, t5, t6, t15;
  double t=TDatabase::TimeDB->CURRENTTIME;

  t1 = std::sin(t);
  t2 = 0.3141592653589793E1*x;
  t4 = t1*std::cos(t2);
  t5 = 0.3141592653589793E1*y;
  t6 = std::cos(t5);
  t15 = 0.3141592653589793E1*0.3141592653589793E1;

  values[0] = t4*t6;
  values[1] = -t1*std::sin(t2)*t6*0.3141592653589793E1;
  values[2] = -t4*0.3141592653589793E1*std::sin(t5);
  values[3] = -2.0*t4*t15*t6;
}

void ExactP(double x, double y, double *values)
{
  double t1, t2, t4;
  double t = TDatabase::TimeDB->CURRENTTIME;

  t1 = std::sin(t);
  t2 = 0.3141592653589793E1*x;
  t4 = 0.3141592653589793E1*y;

  values[0] = t1*(std::sin(t2)+std::cos(t4)-2.0/0.3141592653589793E1);
  values[1] = t1*std::cos(t2)*0.3141592653589793E1;
  values[2] = -t1*std::sin(t4)*0.3141592653589793E1;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

void U1BoundValue(int, double, double &value)
{
  value = 0;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  switch(BdComp)
  {
    case 0: value=std::sin(t)*std::cos(M_PI*Param);
            break;
    case 1: value=-std::sin(t)*std::cos(M_PI*Param);
            break;
    case 2: value=-std::sin(t)*std::cos(M_PI*(1-Param));
            break;
    case 3: value=std::sin(t)*std::cos(M_PI*(1-Param));
            break;
    default: Output::print("wrong boundary part number");
            break;
  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *X, const double *Y,
               const double *const*, double **coeffs)
{
  static double nu = 1/TDatabase::ParamDB->RE_NR;
  double t = TDatabase::TimeDB->CURRENTTIME;
  double t1, t2, t3, t5, t6, t8, t9, t10, t14, t16;
  double t17, t22, t23, t27, t36, t39, t41;
  int i;
  double *coeff, x, y; 
  static double a=1;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];

    coeff[0] = nu;

    t1 = std::cos(t);
    t2 = 0.3141592653589793E1*x;
    t3 = std::sin(t2);
    t5 = 0.3141592653589793E1*y;
    t6 = std::sin(t5);
    t8 = std::sin(t);
    t9 = nu*t8;
    t10 = 0.3141592653589793E1*0.3141592653589793E1;
    t14 = t8*t8;
    t16 = t6*t6;
    t17 = std::cos(t2);
    t22 = std::cos(t5);
    t23 = t22*t22;
    t27 = a*t8;
    t36 = t3*t3;
    t39 = t6*t22*0.3141592653589793E1;
    t41 = t17*t17;

/*
    // Stokes
    coeff[1] = t1*t3*t6+2.0*t9*t3*t10*t6
              +t27*t17*0.3141592653589793E1;
    coeff[2] = t1*t17*t22+2.0*t9*t17*t10*t22
              -t27*t6*0.3141592653589793E1;
*/

    // Navier-Stokes
    coeff[1] = t1*t3*t6+2.0*t9*t3*t10*t6
              +t14*t3*t16*t17*0.3141592653589793E1
              +t14*t17*t23*t3*0.3141592653589793E1
              +t27*t17*0.3141592653589793E1;
    coeff[2] = t1*t17*t22+2.0*t9*t17*t10*t22
              -t14*t36*t39-t14*t41*t39
              -t27*t6*0.3141592653589793E1;
    coeff[3] = 0;
  }
}
