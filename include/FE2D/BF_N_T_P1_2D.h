// ***********************************************************************
// P1 element, nonconforming, 2D
// ***********************************************************************

// base function values
static void N_T_P1_2D_Funct(double xi, double eta, double *values)
{
  values[0]=     -2*eta+1;
  values[1]= 2*xi+2*eta-1;
  values[2]=-2*xi      +1;
}

// values of the derivatives in xi direction
static void N_T_P1_2D_DeriveXi(double, double, double *values)
{
  values[0]= 0;
  values[1]= 2;
  values[2]=-2;
}

// values of the derivatives in eta direction
static void N_T_P1_2D_DeriveEta(double, double, double *values)
{
  values[0]=-2;
  values[1]= 2;
  values[2]= 0;
}

// values of the derivatives in xi-xi direction
static void N_T_P1_2D_DeriveXiXi(double, double, double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
}

// values of the derivatives in eta-eta direction
static void N_T_P1_2D_DeriveEtaEta(double, double, double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
}

// values of the derivatives in xi-eta direction
static void N_T_P1_2D_DeriveXiEta(double, double, double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
}
