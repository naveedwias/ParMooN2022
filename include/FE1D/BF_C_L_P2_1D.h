// ***********************************************************************
// P2 element, conforming, 1D
// ***********************************************************************

// base function values
static void C_L_P2_1D_Funct(double xi, double *values)
{
  values[0]=0.5*xi*(xi-1);
  values[1]=(1+xi)*(1-xi);
  values[2]=0.5*xi*(xi+1);
}

// values of the derivatives in xi direction
static void C_L_P2_1D_DeriveXi(double xi, double *values)
{
  values[0]=xi-0.5;
  values[1]=-2*xi;
  values[2]=xi+0.5;
}

// values of the derivatives in xi-xi  direction
static void C_L_P2_1D_DeriveXiXi(double, double *values)
{
  values[0]=1;
  values[1]=-2;
  values[2]=1;
}
