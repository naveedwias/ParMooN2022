// ***********************************************************************
// Q2 element with bubbles, conforming, 2D
// ***********************************************************************

// base function values
static void C_Q_UL2S_2D_Funct(double xi, double eta, double *values)
{
  double t1, t3, t4, t5, t7, t8, t9, t10, t14, t18;

  t1 = xi*eta;
  t3 = xi*xi;
  t4 = 1.0-t3;
  t5 = (1.0-eta)*t4;
  t7 = eta*eta;
  t8 = 1.0-t7;
  t9 = (1.0-xi)*t8;
  t10 = t4*t8;
  t14 = (1.0+xi)*t8;
  t18 = (1.0+eta)*t4;

  values[0] = 1.0/4.0-xi/4.0-eta/4.0+t1/4.0-t5/4.0-t9/4.0+t10/4.0;
  values[1] = t5/2.0-t10/2.0;
  values[2] = 1.0/4.0+xi/4.0-eta/4.0-t1/4.0-t5/4.0-t14/4.0+t10/4.0;
  values[3] = t14/2.0-t10/2.0;
  values[4] = 1.0/4.0+xi/4.0+eta/4.0+t1/4.0-t18/4.0-t14/4.0+t10/4.0;
  values[5] = t18/2.0-t10/2.0;
  values[6] = 1.0/4.0-xi/4.0+eta/4.0-t1/4.0-t18/4.0-t9/4.0+t10/4.0;
  values[7] = t9/2.0-t10/2.0;
  values[8] = t10;
}

// values of the derivatives in xi direction
static void C_Q_UL2S_2D_DeriveXi(double xi, double eta, double *values)
{
  double t1, t3, t4, t5, t6, t8, t9, t13, t16, t17;

  t1 = eta/4.0;
  t3 = (1.0-eta)*xi;
  t4 = t3/2.0;
  t5 = eta*eta;
  t6 = t5/4.0;
  t8 = xi*(1.0-t5);
  t9 = t8/2.0;
  t13 = t5/2.0;
  t16 = (1.0+eta)*xi;
  t17 = t16/2.0;

  values[0] = t1+t4-t6-t9;
  values[1] = -t3+t8;
  values[2] = -t1+t4+t6-t9;
  values[3] = 1.0/2.0-t13+t8;
  values[4] = t1+t17+t6-t9;
  values[5] = -t16+t8;
  values[6] = -t1+t17-t6-t9;
  values[7] = -1.0/2.0+t13+t8;
  values[8] = -2.0*t8;
}

// values of the derivatives in eta direction
static void C_Q_UL2S_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t5, t6, t8, t9, t11, t14, t15;

  t1 = xi/4.0;
  t2 = xi*xi;
  t3 = t2/4.0;
  t5 = (1.0-xi)*eta;
  t6 = t5/2.0;
  t8 = (1.0-t2)*eta;
  t9 = t8/2.0;
  t11 = t2/2.0;
  t14 = (1.0+xi)*eta;
  t15 = t14/2.0;

  values[0] = t1-t3+t6-t9;
  values[1] = -1.0/2.0+t11+t8;
  values[2] = -t1-t3+t15-t9;
  values[3] = -t14+t8;
  values[4] = t1+t3+t15-t9;
  values[5] = 1.0/2.0-t11+t8;
  values[6] = -t1+t3+t6-t9;
  values[7] = -t5+t8;
  values[8] = -2.0*t8;
}

// values of the derivatives in xi-xi  direction
static void C_Q_UL2S_2D_DeriveXiXi(double, double eta, double *values)
{
  double t1, t2, t3, t4;

  t1 = eta*eta;
  t2 = -eta+t1;
  t3 = 1.0-t1;
  t4 = eta+t1;

  values[0] = t2/2.0;
  values[1] = -t2;
  values[2] = t2/2.0;
  values[3] = t3;
  values[4] = t4/2.0;
  values[5] = -t4;
  values[6] = t4/2.0;
  values[7] = t3;
  values[8] = -2.0*t3;
}

// values of the derivatives in xi-eta direction
static void C_Q_UL2S_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t1, t2, t3, t5;

  t1 = xi/2.0;
  t2 = eta/2.0;
  t3 = xi*eta;
  t5 = 2.0*t3;

  values[0] = 1.0/4.0-t1-t2+t3;
  values[1] = xi-t5;
  values[2] = -1.0/4.0-t1+t2+t3;
  values[3] = -eta-t5;
  values[4] = 1.0/4.0+t1+t2+t3;
  values[5] = -xi-t5;
  values[6] = -1.0/4.0+t1-t2+t3;
  values[7] = eta-t5;
  values[8] = 4.0*t3;
}

// values of the derivatives in eta-eta direction
static void C_Q_UL2S_2D_DeriveEtaEta(double xi, double, double *values)
{
  double t1, t2, t3, t4;

  t1 = xi*xi;
  t2 = -xi+t1;
  t3 = 1.0-t1;
  t4 = xi+t1;

  values[0] = t2/2.0;
  values[1] = t3;
  values[2] = t4/2.0;
  values[3] = -t4;
  values[4] = t4/2.0;
  values[5] = t3;
  values[6] = t2/2.0;
  values[7] = -t2;
  values[8] = -2.0*t3;
}
