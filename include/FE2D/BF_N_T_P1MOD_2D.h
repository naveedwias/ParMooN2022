// ***********************************************************************
// P1MOD element, nonconforming, 2D
// ***********************************************************************

// base function values
static void N_T_P1MOD_2D_Funct(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
  double t15, t16, t20, t21, t23;

  t1 = 12.0*eta;
  t2 = xi*xi;
  t3 = t2*eta;
  t4 = 20.0*t3;
  t5 = eta*eta;
  t6 = xi*t5;
  t7 = 20.0*t6;
  t8 = 30.0*t5;
  t9 = t5*eta;
  t10 = 20.0*t9;
  t11 = xi*eta;
  t12 = 20.0*t11;
  t15 = 2.0*t11;
  t16 = t2*xi;
  t20 = 12.0*xi;
  t21 = 30.0*t2;
  t23 = 20.0*t16;

  values[0] = 1.0-t1-t4-t7+t8-t10+t12;
  values[1] = -xi+3.0*t2+t15-2.0*t16-3.0*t3-t6;
  values[2] = -1.0+t20+t1-t21-40.0*t11+t23+40.0*t3+40.0*t6-t8+t10;
  values[3] = -t3+t6;
  values[4] = 1.0-t20+t21+t12-t23-t4-t7;
  values[5] = -3.0*t5+3.0*t6+2.0*t9+eta-t15+t3;
}

// values of the derivatives in xi direction
static void N_T_P1MOD_2D_DeriveXi(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t8, t9, t13, t15, t19;

  t1 = xi*eta;
  t2 = 40.0*t1;
  t3 = eta*eta;
  t4 = 20.0*t3;
  t5 = 20.0*eta;
  t8 = 2.0*eta;
  t9 = xi*xi;
  t13 = 60.0*xi;
  t15 = 60.0*t9;
  t19 = 2.0*t1;

  values[0] = -t2-t4+t5;
  values[1] = -1.0+6.0*xi+t8-6.0*t9-6.0*t1-t3;
  values[2] = 12.0-t13-40.0*eta+t15+80.0*t1+40.0*t3;
  values[3] = -t19+t3;
  values[4] = -12.0+t13+t5-t15-t2-t4;
  values[5] = 3.0*t3-t8+t19;
}

// values of the derivatives in eta direction
static void N_T_P1MOD_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t8, t10, t12;

  t1 = xi*xi;
  t2 = 20.0*t1;
  t3 = xi*eta;
  t4 = 40.0*t3;
  t5 = 60.0*eta;
  t6 = eta*eta;
  t7 = 60.0*t6;
  t8 = 20.0*xi;
  t10 = 2.0*xi;
  t12 = 2.0*t3;

  values[0] = -12.0-t2-t4+t5-t7+t8;
  values[1] = t10-3.0*t1-t12;
  values[2] = 12.0-40.0*xi+40.0*t1+80.0*t3-t5+t7;
  values[3] = -t1+t12;
  values[4] = t8-t2-t4;
  values[5] = -6.0*eta+6.0*t3+6.0*t6+1.0-t10+t1;
}

// values of the derivatives in xi-xi direction
static void N_T_P1MOD_2D_DeriveXiXi(double xi, double eta, 
                                       double *values)
{
  double t1, t5, t8;

  t1 = 40.0*eta;
  t5 = 120.0*xi;
  t8 = 2.0*eta;

  values[0] = -t1;
  values[1] = 6.0-12.0*xi-6.0*eta;
  values[2] = -60.0+t5+80.0*eta;
  values[3] = -t8;
  values[4] = 60.0-t5-t1;
  values[5] = t8;
}

// values of the derivatives in xi-eta direction
static void N_T_P1MOD_2D_DeriveXiEta(double xi, double eta, 
                                       double *values)
{
  double t3;

  t3 = -40.0*xi-40.0*eta+20.0;

  values[0] = t3;
  values[1] = 2.0-6.0*xi-2.0*eta;
  values[2] = -40.0+80.0*xi+80.0*eta;
  values[3] = -2.0*xi+2.0*eta;
  values[4] = t3;
  values[5] = 6.0*eta-2.0+2.0*xi;
}

// values of the derivatives in eta-eta direction
static void N_T_P1MOD_2D_DeriveEtaEta(double xi, double eta, 
                                       double *values)
{
  double t1, t2, t4;

  t1 = 40.0*xi;
  t2 = 120.0*eta;
  t4 = 2.0*xi;

  values[0] = -t1+60.0-t2;
  values[1] = -t4;
  values[2] = 80.0*xi-60.0+t2;
  values[3] = t4;
  values[4] = -t1;
  values[5] = -6.0+6.0*xi+12.0*eta;
}

static int N_T_P1MOD_2D_ChangeJ0[1] = { 1 };
static int N_T_P1MOD_2D_ChangeJ1[1] = { 3 };
static int N_T_P1MOD_2D_ChangeJ2[1] = { 5 };

static int *N_T_P1MOD_2D_Change1[3] = {
                N_T_P1MOD_2D_ChangeJ0, N_T_P1MOD_2D_ChangeJ1,
                N_T_P1MOD_2D_ChangeJ2 };
static int **N_T_P1MOD_2D_Change[] = { N_T_P1MOD_2D_Change1 };
