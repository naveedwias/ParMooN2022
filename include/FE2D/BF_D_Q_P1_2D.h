// ***********************************************************************
// P1 element, discontinous, 2D, quadrilateral
// ***********************************************************************

// base function values
static void D_Q_P1_2D_Funct(double xi, double eta, double *values)
{
  values[0] = 1;
  values[1] = 3*xi;
  values[2] = 3*eta;
}

// values of the derivatives in xi direction
static void D_Q_P1_2D_DeriveXi(double, double, double *values)
{
  values[0] = 0;
  values[1] = 3;
  values[2] = 0;
}

// values of the derivatives in eta direction
static void D_Q_P1_2D_DeriveEta(double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 3;
}

// values of the derivatives in xi-xi direction
static void D_Q_P1_2D_DeriveXiXi(double, double, 
                                       double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
}

// values of the derivatives in eta-eta direction
static void D_Q_P1_2D_DeriveEtaEta(double, double, 
                                       double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
}

// values of the derivatives in xi-eta direction
static void D_Q_P1_2D_DeriveXiEta(double, double, 
                                       double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
}
