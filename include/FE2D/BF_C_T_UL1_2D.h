// ***********************************************************************
// P1 element with bubble, conforming, 2D
// ***********************************************************************

// base function values
static void C_T_UL1_2D_Funct(double xi, double eta, double *values)
{
  double t3, t4;
  
  t3 = xi*eta*(1.0-xi-eta);
  t4 = 20.0*t3;

  values[0] = 1.0-xi-eta-t4;
  values[1] = xi-t4;
  values[2] = eta-t4;
  values[3] = 60.0*t3;
}

// values of the derivatives in xi direction
static void C_T_UL1_2D_DeriveXi(double xi, double eta, double *values)
{
  double t2, t3, t4, t5, t8;

  t2 = eta*(1.0-xi-eta);
  t3 = 20.0*t2;
  t4 = xi*eta;
  t5 = 20.0*t4;
  t8 = -t2+t4;

  values[0] = -1.0-t3+t5;
  values[1] = 1.0-t3+t5;
  values[2] = 20.0*t8;
  values[3] = -60.0*t8;
}

// values of the derivatives in eta direction
static void C_T_UL1_2D_DeriveEta(double xi, double eta, double *values)
{
  double t2, t3, t4, t5, t7;

  t2 = xi*(1.0-xi-eta);
  t3 = 20.0*t2;
  t4 = xi*eta;
  t5 = 20.0*t4;
  t7 = -t2+t4;

  values[0] = -1.0-t3+t5;
  values[1] = 20.0*t7;
  values[2] = 1.0-t3+t5;
  values[3] = -60.0*t7;
}

// values of the derivatives in xi-xi  direction
static void C_T_UL1_2D_DeriveXiXi(double, double eta, double *values)
{
  double t1;

  t1 = 40.0*eta;

  values[0] = t1;
  values[1] = t1;
  values[2] = t1;
  values[3] = -120.0*eta;
}

// values of the derivatives in xi-eta direction
static void C_T_UL1_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t3;

  t3 = -20.0+40.0*xi+40.0*eta;

  values[0] = t3;
  values[1] = t3;
  values[2] = t3;
  values[3] = 60.0-120.0*xi-120.0*eta;
}

// values of the derivatives in eta-eta direction
static void C_T_UL1_2D_DeriveEtaEta(double xi, double, double *values)
{
  double t1;

  t1 = 40.0*xi;

  values[0] = t1;
  values[1] = t1;
  values[2] = t1;
  values[3] = -120.0*xi;
}
