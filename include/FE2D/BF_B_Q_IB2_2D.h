// ***********************************************************************
// internal bubble of degree 2 (in the sense of Q2)
// ***********************************************************************

// base function values
static void B_Q_IB2_2D_Funct(double xi, double eta, double *values)
{
  values[0]=(1.0-xi)*(1.0-eta)*(1.0+xi)*(1.0+eta);
}

// values of the derivatives in xi direction
static void B_Q_IB2_2D_DeriveXi(double xi, double eta, double *values)
{
  values[0]=-2.*xi+2.*eta*eta*xi;
}

// values of the derivatives in eta direction
static void B_Q_IB2_2D_DeriveEta(double xi, double eta, double *values)
{
  values[0]=-2.*eta+2.*xi*xi*eta;
}
// values of the derivatives in xi-xi  direction
static void B_Q_IB2_2D_DeriveXiXi(double, double eta, double *values)
{
  values[0]=-2.+2.*eta*eta;
}
// values of the derivatives in xi-eta direction
static void B_Q_IB2_2D_DeriveXiEta(double xi, double eta, double *values)
{
  values[0]=4.*eta*xi;
}
// values of the derivatives in eta-eta direction
static void B_Q_IB2_2D_DeriveEtaEta(double xi, double, double *values)
{
  values[0]=-2.+2.*xi*xi;
}
