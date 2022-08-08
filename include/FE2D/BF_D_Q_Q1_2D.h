// ***********************************************************************
// Q1 element, discontinous, 2D, quadrilateral
// ***********************************************************************

// base function values
static void D_Q_Q1_2D_Funct(double xi, double eta, double *values)
{
  values[0] = 1.0;
  values[1] = xi;
  values[2] = eta;
  values[3] = xi*eta;
}

// values of the derivatives in xi direction
static void D_Q_Q1_2D_DeriveXi(double, double eta, double *values)
{
  values[0] = 0.0;
  values[1] = 1.0;
  values[2] = 0.0;
  values[3] = eta;
}

// values of the derivatives in eta direction
static void D_Q_Q1_2D_DeriveEta(double xi, double, double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 1.0;
  values[3] = xi;
}

// values of the derivatives in xi-xi direction
static void D_Q_Q1_2D_DeriveXiXi(double, double, double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
}

// values of the derivatives in xi-eta direction
static void D_Q_Q1_2D_DeriveXiEta(double, double, double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 1.0;
}

// values of the derivatives in eta-eta direction
static void D_Q_Q1_2D_DeriveEtaEta(double, double, double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
}
