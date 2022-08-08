// ***********************************************************************
// P2 element, discontinous, 2D, quadrilateral
// ***********************************************************************

// base function values
static void D_Q_P2_2D_Funct(double xi, double eta, double *values)
{
  values[0] = 1;
  values[1] = 3*xi;
  values[2] = 3*eta;
  values[3] = 1.25*(3*xi*xi-1);
  values[4] = 9*xi*eta;
  values[5] = 1.25*(3*eta*eta-1);
}

// values of the derivatives in xi direction
static void D_Q_P2_2D_DeriveXi(double xi, double eta, double *values)
{
  values[0] = 0;
  values[1] = 3;
  values[2] = 0;
  values[3] = 7.5*xi;
  values[4] = 9*eta;
  values[5] = 0;
}

// values of the derivatives in eta direction
static void D_Q_P2_2D_DeriveEta(double xi, double eta, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 3;
  values[3] = 0;
  values[4] = 9*xi;
  values[5] = 7.5*eta;
}

// values of the derivatives in xi-xi direction
static void D_Q_P2_2D_DeriveXiXi(double, double, double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
  values[3]=7.5;
  values[4]=0;
  values[5]=0;
}

// values of the derivatives in eta-eta direction
static void D_Q_P2_2D_DeriveEtaEta(double, double, double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
  values[3]=0;
  values[4]=0;
  values[5]=7.5;
}

// values of the derivatives in xi-eta direction
static void D_Q_P2_2D_DeriveXiEta(double, double, double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
  values[3]=0;
  values[4]=9;
  values[5]=0;
}
