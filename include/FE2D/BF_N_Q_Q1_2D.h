// ***********************************************************************
// Q1 element, nonconforming, 2D
// ***********************************************************************

// base function values
static void N_Q_Q1_2D_Funct(double xi, double eta, double *values)
{
  double h=0.375*(xi*xi-eta*eta);
  values[0]=-h       -0.5*eta+0.25;
  values[1]= h+0.5*xi        +0.25;
  values[2]=-h       +0.5*eta+0.25;
  values[3]= h-0.5*xi        +0.25;
}

// values of the derivatives in xi direction
static void N_Q_Q1_2D_DeriveXi(double xi, double, double *values)
{
  values[0]=-0.75*xi;
  values[1]= 0.75*xi+0.5;
  values[2]=-0.75*xi;
  values[3]= 0.75*xi-0.5;
}

// values of the derivatives in eta direction
static void N_Q_Q1_2D_DeriveEta(double, double eta, double *values)
{
  values[0]= 0.75*eta-0.5;
  values[1]=-0.75*eta;
  values[2]= 0.75*eta+0.5;
  values[3]=-0.75*eta;
}

// values of derivatives in xi-xi direction
static void N_Q_Q1_2D_DeriveXiXi(double, double, double *values)
{
  values[0]=-0.75;
  values[1]= 0.75;
  values[2]=-0.75;
  values[3]= 0.75;
}

// values of derivatives in eta-eta direction
static void N_Q_Q1_2D_DeriveEtaEta(double, double, double *values)
{
  values[0]= 0.75;
  values[1]=-0.75;
  values[2]= 0.75;
  values[3]=-0.75;
}

// values of derivatives in xi-eta direction
static void N_Q_Q1_2D_DeriveXiEta(double, double, double *values)
{
  values[0]=0;
  values[1]=0;
  values[2]=0;
  values[3]=0;
}
