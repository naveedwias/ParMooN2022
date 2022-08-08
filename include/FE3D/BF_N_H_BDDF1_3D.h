// ***********************************************************************
// Brezzi-Douglas-Duran-Fortin element of first order, 3D
// ***********************************************************************

static double N_H_BDDF1_3D_CM[324] = {
//using tschebyscheff points (see NF_N_H_BDDF1_3D.h)
    -0.0441941738,0.0441941738,0,-0.0441941738,-0,0.0441941738,0,0.0625,0.0625,-0.0441941738,0,0.0441941738,0,-0.0625,-0.0625,-0.0441941738,0,0.0441941738,
    -0.0441941738,0,0.0441941738,0,-0.0625,-0.0625,-0.0441941738,0,0.0441941738,0,0.0625,0.0625,-0.0441941738,0.0441941738,0,-0.0441941738,0.0441941738,-0,
    0,-0.0625,-0.0625,-0.0441941738,0.0441941738,0,-0.0441941738,0.0441941738,0,-0.0441941738,0.0441941738,0,-0.0441941738,-0,0.0441941738,0,0.0625,0.0625,
    0,0,0,0,-0,0,0,0.0625,0.0625,0,0,0,0,0.0625,0.0625,0,0,-0,
    0,0,0,0.0883883476,-0,-0.0883883476,0,0,0,-0.0883883476,0,0.0883883476,0,0,0,0,0,-0,
    0.0883883476,-0.0883883476,0,0,-0,0,0,0,0,0,0,0,0,0,0,-0.0883883476,0,0.0883883476,
    0,0,0,0,-0,0,-0.0883883476,0,0.0883883476,0,0,0,0.0883883476,-0.0883883476,0,0,0,-0,
    0,0,0,0,0.0625,0.0625,0,0,0,0,0.0625,0.0625,0,-0,0,0,0,0,
    0.0883883476,0,-0.0883883476,0,-0,0,0,0,0,0,0,0,0,0,0,-0.0883883476,0.0883883476,-0,
    0,0,0,0,-0,0,-0.0883883476,0.0883883476,0,0,0,0,0.0883883476,0,-0.0883883476,0,0,-0,
    0,0,0,0.0883883476,-0.0883883476,0,0,0,0,-0.0883883476,0.0883883476,0,0,0,0,0,0,-0,
    0,0.0625,0.0625,0,-0,0,0,0,0,0,0,0,0,0,0,0,0.0625,0.0625,
    0,0,0,0.0441941738,-0,-0.0441941738,0,0,0,0.0441941738,0,-0.0441941738,0,0,0,0,0,-0,
    0,0,0,0,-0,0,0.0441941738,-0.0441941738,0,0,0,0,0.0441941738,0,-0.0441941738,0,0,-0,
    0.0441941738,0,-0.0441941738,0,-0,0,0,0,0,0,0,0,0,0,0,0.0441941738,-0.0441941738,-0,
    0,0,0,0,-0,0,-0.0441941738,0,0.0441941738,0,0,0,-0.0441941738,0.0441941738,0,0,0,-0,
    -0.0441941738,0.0441941738,0,0,-0,0,0,0,0,0,0,0,0,0,0,-0.0441941738,0,0.0441941738,
    0,0,0,-0.0441941738,0.0441941738,0,0,0,0,-0.0441941738,0.0441941738,0,0,0,0,0,0,0
};

static void N_H_BDDF1_3D_Funct(double xi, double eta, double zeta,
                               double *values)
{
  int nBF = 18; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={1,0,0,xi,0,0,eta,0,0,zeta,0,0,
                  xi*xi,-2*xi*zeta,0,2*xi*eta,-xi*xi,0};
  double mon_y[]={0,1,0,0,xi,0,0,eta,0,0,zeta,0,
                  -2*xi*eta,0,eta*eta,-eta*eta,0,2*eta*zeta};
  double mon_z[]={0,0,1,0,0,xi,0,0,eta,0,0,zeta,
                  0,zeta*zeta,-2*eta*zeta,0,2*xi*zeta,-zeta*zeta};
  
  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF1_3D_DeriveXi(double xi, double eta, double zeta,
                                  double *values)
{
  int nBF = 18; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,1,0,0,0,0,0,0,0,0, 2*xi,-2*zeta,0,2*eta,-2*xi,0};
  double mon_y[]={0,0,0,0,1,0,0,0,0,0,0,0, -2*eta,0,0,0,0,0};
  double mon_z[]={0,0,0,0,0,1,0,0,0,0,0,0, 0,0,0,0,2*zeta,0};
  
  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF1_3D_DeriveEta(double xi, double eta, double zeta,
                                   double *values)
{
  int nBF = 18; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,1,0,0,0,0,0, 0,0,0,2*xi,0,0};
  double mon_y[]={0,0,0,0,0,0,0,1,0,0,0,0, -2*xi,0,2*eta,-2*eta,0,2*zeta};
  double mon_z[]={0,0,0,0,0,0,0,0,1,0,0,0, 0,0,-2*zeta,0,0,0};
  
  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF1_3D_DeriveZeta(double xi, double eta, double zeta,
                                    double *values)
{
  int nBF = 18; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,1,0,0, 0,-2*xi,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,1,0, 0,0,0,0,0,2*eta};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,1, 0,2*zeta,-2*eta,0,2*xi,-2*zeta};
  
  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF1_3D_DeriveXiXi(double, double, double, double *values)
{
  int nBF = 18; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0, 2,0,0,0,-2,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0};

  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF1_3D_DeriveXiEta(double, double, double, double *values)
{
  int nBF = 18; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,2,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0, -2,0,0,0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0};

  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF1_3D_DeriveXiZeta(double, double, double, double *values)
{
  int nBF = 18; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,-2,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,2,0};

  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF1_3D_DeriveEtaEta(double, double, double, double *values)
{
  int nBF = 18; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,2,-2,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0};

  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF1_3D_DeriveEtaZeta(double, double, double, double *values)
{
  int nBF = 18; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,2};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,-2,0,0,0};

  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF1_3D_DeriveZetaZeta(double, double, double, double *values)
{
  int nBF = 18; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0, 0,2,0,0,0,-2};

  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static int N_H_BDDF1_3D_ChangeJ0[3] = { 0, 1, 2 };
static int N_H_BDDF1_3D_ChangeJ1[3] = { 3, 4, 5 };
static int N_H_BDDF1_3D_ChangeJ2[3] = { 6, 7, 8 };
static int N_H_BDDF1_3D_ChangeJ3[3] = { 9,10,11 };
static int N_H_BDDF1_3D_ChangeJ4[3] = {12,13,14 };
static int N_H_BDDF1_3D_ChangeJ5[3] = {15,16,17 };

static int *N_H_BDDF1_3D_Change1[6] = { N_H_BDDF1_3D_ChangeJ0,
                                        N_H_BDDF1_3D_ChangeJ1,
                                        N_H_BDDF1_3D_ChangeJ2,
                                        N_H_BDDF1_3D_ChangeJ3,
                                        N_H_BDDF1_3D_ChangeJ4,
                                        N_H_BDDF1_3D_ChangeJ5 };
static int **N_H_BDDF1_3D_Change2 = N_H_BDDF1_3D_Change1;
static int **N_H_BDDF1_3D_Change[] = { N_H_BDDF1_3D_Change1, N_H_BDDF1_3D_Change2 };
