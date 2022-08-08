// ***********************************************************************
// Brezzi-Douglas-Duran-Fortin element of first order on tetrahedra, 3D
// ***********************************************************************

static double N_T_BDDF1_3D_CM[144] = {
    0,-0,0,0,0,-0,0,0,-0,-1.666666667,0.3333333333,0.3333333333,
    0,0,0,-1.666666667,0.3333333333,0.3333333333,0,0,0,0,0,-0,
    -1.666666667,0.3333333333,0.3333333333,0,0,0,0,0,0,0,0,-0,
    -0.3333333333,1.666666667,-0.3333333333,-0.3333333333,-0.3333333333,1.666666667,-0.3333333333,1.666666667,-0.3333333333,1.666666667,-0.3333333333,-0.3333333333,
    0,0,0,2,0,-2,0,0,0,0,0,-0,
    2,-2,0,0,0,0,0,0,0,0,0,-0,
    -0,-0,0,0,-0,-0,0,-0,-0,2,-2,-0,
    -0.3333333333,-0.3333333333,1.666666667,1.666666667,-0.3333333333,-0.3333333333,1.666666667,-0.3333333333,-0.3333333333,-0.3333333333,1.666666667,-0.3333333333,
    2,0,-2,0,0,0,0,0,0,0,0,-0,
    0,0,0,0,0,0,0,0,0,2,0,-2,
    0,0,0,2,-2,0,0,0,0,0,0,-0,
    1.666666667,-0.3333333333,-0.3333333333,-0.3333333333,1.666666667,-0.3333333333,-0.3333333333,-0.3333333333,1.666666667,-0.3333333333,-0.3333333333,1.666666667
};

static void N_T_BDDF1_3D_Funct(double xi, double eta, double zeta,
                                                         double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={1,0,0,xi,0,0,eta,0,0,zeta,0,0};
  double mon_y[]={0,1,0,0,xi,0,0,eta,0,0,zeta,0};
  double mon_z[]={0,0,1,0,0,xi,0,0,eta,0,0,zeta};
  
  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_BDDF1_3D_DeriveXi(double, double, double, double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,1,0,0,0,0,0,0,0,0};
  double mon_y[]={0,0,0,0,1,0,0,0,0,0,0,0};
  double mon_z[]={0,0,0,0,0,1,0,0,0,0,0,0};
  
  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_BDDF1_3D_DeriveEta(double, double, double, double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,1,0,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,1,0,0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,1,0,0,0};
  
  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_BDDF1_3D_DeriveZeta(double, double, double, double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,1,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,1,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,1};
  
  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_T_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_T_BDDF1_3D_DeriveXiXi(double, double, double, double *values)
{
  int nBF = 12; // number of basis functions
  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
}

static void N_T_BDDF1_3D_DeriveXiEta(double, double, double, double *values)
{
  int nBF = 12; // number of basis functions
  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
}

static void N_T_BDDF1_3D_DeriveXiZeta(double, double, double, double *values)
{
  int nBF = 12; // number of basis functions
  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
}

static void N_T_BDDF1_3D_DeriveEtaEta(double, double, double, double *values)
{
  int nBF = 12; // number of basis functions
  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
}

static void N_T_BDDF1_3D_DeriveEtaZeta(double, double, double, double *values)
{
  int nBF = 12; // number of basis functions
  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
}

static void N_T_BDDF1_3D_DeriveZetaZeta(double, double, double, double *values)
{
  int nBF = 12; // number of basis functions
  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
}

static int N_T_BDDF1_3D_ChangeJ0[3] = { 0, 1, 2 };
static int N_T_BDDF1_3D_ChangeJ1[3] = { 3, 4, 5 };
static int N_T_BDDF1_3D_ChangeJ2[3] = { 6, 7, 8 };
static int N_T_BDDF1_3D_ChangeJ3[3] = { 9,10,11 };

static int *N_T_BDDF1_3D_Change1[4] = { N_T_BDDF1_3D_ChangeJ0, N_T_BDDF1_3D_ChangeJ1,
                                      N_T_BDDF1_3D_ChangeJ2, N_T_BDDF1_3D_ChangeJ3};
static int **N_T_BDDF1_3D_Change2 = N_T_BDDF1_3D_Change1;
static int **N_T_BDDF1_3D_Change[] = { N_T_BDDF1_3D_Change1, N_T_BDDF1_3D_Change2 };
