// ***********************************************************************
// conforming P2 element with cell bubble
// ***********************************************************************

// base function values
static void C_T_B2_2D_Funct(double xi, double eta, double *values)
{
  double t1 = 1-xi-eta;
  double t3 = xi*eta;
  double t4 = t3*t1;

  values[0] = 2*t1*(0.5-xi-eta);
  values[1] = 4*xi*t1-20*t4;
  values[2] = 2*xi*(xi-0.5);
  values[3] = 4*eta*t1-20*t4;
  values[4] = 4*t3-20*t4;
  values[5] = 2*eta*(eta-0.5);
  values[6] = 20*t4;
}

// values of the derivatives in xi direction
static void C_T_B2_2D_DeriveXi(double xi, double eta, double *values)
{
  double t2 = eta*(1-xi-eta);
  double t3 = xi*eta;

  values[0] = -3+4*xi+4*eta;
  values[1] = 4-8*xi-4*eta-20*t2+20*t3;
  values[2] = 4*xi-1;
  values[3] = -4*eta-20*t2+20*t3;
  values[4] = 4*eta-20*t2+20*t3;
  values[5] = 0;
  values[6] = 20*t2-20*t3;
}

// values of the derivatives in eta direction
static void C_T_B2_2D_DeriveEta(double xi, double eta, double *values)
{
  double t2 = xi*(1-xi-eta);
  double t3 = xi*eta;

  values[0] = -3+4*xi+4*eta;
  values[1] = -4*xi-20*t2+20*t3;
  values[2] = 0;
  values[3] = 4-4*xi-8*eta-20*t2+20*t3;
  values[4] = 4*xi-20*t2+20*t3;
  values[5] = 4*eta-1;
  values[6] = 20*t2-20*t3;
}

// values of the derivatives in xi-xi direction
static void C_T_B2_2D_DeriveXiXi(double, double eta, double *values)
{
  values[0] = 4;
  values[1] = -8+40*eta;
  values[2] = 4;
  values[3] = 40*eta;
  values[4] = 40*eta;
  values[5] = 0;
  values[6] = -40*eta;
}

// values of the derivatives in xi-eta direction
static void C_T_B2_2D_DeriveXiEta(double xi, double eta, double *values)
{
  double t1 = -24+40*xi+40*eta;

  values[0] = 4;
  values[1] = t1;
  values[2] = 0;
  values[3] = t1;
  values[4] = -16+40*xi+40*eta;
  values[5] = 0;
  values[6] = 20-40*xi-40*eta;
}

// values of the derivatives in eta-eta direction
static void C_T_B2_2D_DeriveEtaEta(double xi, double, double *values)
{
  values[0] = 4;
  values[1] = 40*xi;
  values[2] = 0;
  values[3] = -8+40*xi;
  values[4] = 40*xi;
  values[5] = 4;
  values[6] = -40*xi;
}
