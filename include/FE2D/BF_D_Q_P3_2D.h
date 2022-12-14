// ***********************************************************************
// P3 element, discontinous, 2D, quadrilateral
// ***********************************************************************

// base function values
static void D_Q_P3_2D_Funct(double xi, double eta, double *values)
{
  values[0] = 1;
  values[1] = 3*xi;
  values[2] = 3*eta;
  values[3] = 1.25*(3*xi*xi-1);
  values[4] = 9*xi*eta;
  values[5] = 1.25*(3*eta*eta-1);
  values[6] = 1.75*xi*(5*xi*xi-3);
  values[7] = 3.75*eta*(3*xi*xi-1);
  values[8] = 3.75*xi*(3*eta*eta-1);
  values[9] = 1.75*eta*(5*eta*eta-3);
}

// values of the derivatives in xi direction
static void D_Q_P3_2D_DeriveXi(double xi, double eta, double *values)
{
  values[0] = 0;
  values[1] = 3;
  values[2] = 0;
  values[3] = 7.5*xi;
  values[4] = 9*eta;
  values[5] = 0;
  values[6] = 26.25*xi*xi-5.25;
  values[7] = 22.5*xi*eta;
  values[8] = 11.25*eta*eta-3.75;
  values[9] = 0;
}

// values of the derivatives in eta direction
static void D_Q_P3_2D_DeriveEta(double xi, double eta, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 3;
  values[3] = 0;
  values[4] = 9*xi;
  values[5] = 7.5*eta;
  values[6] = 0;
  values[7] = 11.25*xi*xi-3.75;
  values[8] = 22.5*xi*eta;
  values[9] = 26.25*eta*eta-5.25;
}

// values of the derivatives in xi-xi direction
static void D_Q_P3_2D_DeriveXiXi(double xi, double eta, 
                                       double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
  values[3]=7.5;
  values[4]=0;
  values[5]=0;
  values[6]=52.5*xi;
  values[7]=22.5*eta;
  values[8]=0;
  values[9]=0;
}

// values of the derivatives in eta-eta direction
static void D_Q_P3_2D_DeriveEtaEta(double xi, double eta, 
                                       double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
  values[3]=0;
  values[4]=0;
  values[5]=7.5;
  values[6]=0;
  values[7]=0;
  values[8]=22.5*xi;
  values[9]=52.5*eta;
}

// values of the derivatives in xi-eta direction
static void D_Q_P3_2D_DeriveXiEta(double xi, double eta, 
                                       double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
  values[3]=0;
  values[4]=9;
  values[5]=0;
  values[6]=0;
  values[7]=22.5*xi;
  values[8]=22.5*eta;
  values[9]=0;
}
