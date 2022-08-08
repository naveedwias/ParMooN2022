// ***********************************************************************
// P0 element, conforming, 1D
// ***********************************************************************

// base function values
static void C_L_P0_1D_Funct(double, double *values)
{
  values[0]=1;
}

// values of the derivatives in xi direction
static void C_L_P0_1D_DeriveXi(double, double *values)
{
  values[0]=0;
}

// values of the derivatives in xi-xi  direction
static void C_L_P0_1D_DeriveXiXi(double, double *values)
{
  values[0]=0;
}
