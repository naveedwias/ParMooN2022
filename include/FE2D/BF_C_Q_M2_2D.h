// ***********************************************************************
// Q2 element with bubbles, conforming, 2D
// ***********************************************************************

// base function values
static void C_Q_M2_2D_Funct(double xi, double eta, double *values)
{
  double t1 = xi/4.0;
  double t2 = eta/4.0;
  double t3 = xi*xi;
  double t4 = t3/2.0;
  double t5 = xi*eta;
  double t6 = t5/4.0;
  double t7 = eta*eta;
  double t8 = t7/2.0;
  double t11 = (1.0-t3)*(1.0+eta);
  double t12 = t11/4.0;
  double t15 = (1.0-t7)*(1.0+xi);
  double t16 = t15/4.0;
  double t18 = t11/2.0;
  double t21 = t15/2.0;

  values[0] = -3.0/4.0-t1-t2+t4+t6+t8+t12+t16;
  values[1] = 1.0-t3-t18;
  values[2] = -1.0/4.0+t1-t2+t4-t6+t12-t16;
  values[3] = t21;
  values[4] = 1.0/4.0+xi/4.0+eta/4.0+t5/4.0-t11/4.0-t15/4.0;
  values[5] = t18;
  values[6] = -1.0/4.0-t1+t2-t6+t8-t12+t16;
  values[7] = 1.0-t7-t21;
}

// values of the derivatives in xi direction
static void C_Q_M2_2D_DeriveXi(double xi, double eta, double *values)
{
  double t1 = eta/4.0;
  double t3 = xi*(1.0+eta);
  double t4 = t3/2.0;
  double t5 = eta*eta;
  double t6 = t5/4.0;
  double t11 = 1.0-t5;

  values[0] = xi+t1-t4-t6;
  values[1] = -2.0*xi+t3;
  values[2] = xi-t1-t4+t6;
  values[3] = t11/2.0;
  values[4] = t1+t4+t6;
  values[5] = -t3;
  values[6] = -t1+t4-t6;
  values[7] = -t11/2.0;
}

// values of the derivatives in eta direction
static void C_Q_M2_2D_DeriveEta(double xi, double eta, double *values)
{
  double t1 = xi/4.0;
  double t2 = xi*xi;
  double t3 = t2/4.0;
  double t5 = eta*(1.0+xi);
  double t6 = t5/2.0;
  double t8 = -1.0+t2;

  values[0] = t1+eta-t3-t6;
  values[1] = t8/2.0;
  values[2] = -t1-t3+t6;
  values[3] = -t5;
  values[4] = t1+t3+t6;
  values[5] = -t8/2.0;
  values[6] = -t1+eta+t3-t6;
  values[7] = -2.0*eta+t5;
}

// values of the derivatives in xi-xi  direction
static void C_Q_M2_2D_DeriveXiXi(double, double eta, double *values)
{
  double t1 = 1.0-eta;
  double t2 = 1.0+eta;

  values[0] = t1/2.0;
  values[1] = -t1;
  values[2] = t1/2.0;
  values[3] = 0.0;
  values[4] = t2/2.0;
  values[5] = -t2;
  values[6] = t2/2.0;
  values[7] = 0.0;
}

// values of the derivatives in xi-eta direction
static void C_Q_M2_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t1 = xi/2.0;
  double t2 = eta/2.0;

  values[0] = 1.0/4.0-t1-t2;
  values[1] = xi;
  values[2] = -1.0/4.0-t1+t2;
  values[3] = -eta;
  values[4] = 1.0/4.0+t1+t2;
  values[5] = -xi;
  values[6] = -1.0/4.0+t1-t2;
  values[7] = eta;
}

// values of the derivatives in eta-eta direction
static void C_Q_M2_2D_DeriveEtaEta(double xi, double, double *values)
{
  double t1 = 1.0-xi;
  double t2 = 1.0+xi;

  values[0] = t1/2.0;
  values[1] = 0.0;
  values[2] = t2/2.0;
  values[3] = -t2;
  values[4] = t2/2.0;
  values[5] = 0.0;
  values[6] = t1/2.0;
  values[7] = -t1;
}
