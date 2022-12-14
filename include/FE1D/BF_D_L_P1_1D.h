// ***********************************************************************
// P1 discontinuous element,  1D
// ***********************************************************************

// base function values
static void D_L_P1_1D_Funct(double xi, double *values)
{
  values[0]=1. + 1.5*xi;
  values[1]=3.*xi;
}

// values of the derivatives in xi direction
static void D_L_P1_1D_DeriveXi(double, double *values)
{
  values[0]=1.5;
  values[1]=3.;
}

// values of the derivatives in xi-xi  direction
static void D_L_P1_1D_DeriveXiXi(double, double *values)
{
  values[0]=0;
  values[1]=0;
}
