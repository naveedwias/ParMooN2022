// Third order Raviart-Thomas vector element, nonconforming, 2D

static double N_T_RT3_2D_CM[576] = {
 0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00, -1.0000e+00, -3.0000e+00, -5.0000e+00, -7.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
-1.0000e+00,  3.0000e+00, -5.0000e+00,  7.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
 1.0000e+00, -4.2000e+00,  1.1000e+01, -2.3800e+01, -1.0000e+01,  6.0000e+00,  0.0000e+00,  0.0000e+00,  1.9000e+01,  4.9800e+01,  5.9000e+01,  3.2200e+01,  9.6000e+02,  2.4000e+02, -3.3600e+03, -1.2000e+03, -3.0000e+03, -8.4000e+02,  2.6880e+03,  1.1760e+03,  2.1840e+03,  6.7200e+02,  4.8720e+03,  1.8480e+03,
 0.0000e+00, -6.0000e+00,  3.0000e+01, -8.4000e+01,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
 0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  6.0000e+00,  3.0000e+01,  8.4000e+01,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
 1.9000e+01, -4.9800e+01,  5.9000e+01, -3.2200e+01, -1.0000e+01, -6.0000e+00,  0.0000e+00,  0.0000e+00,  1.0000e+00,  4.2000e+00,  1.1000e+01,  2.3800e+01,  2.4000e+02,  9.6000e+02, -8.4000e+02, -3.0000e+03, -1.2000e+03, -3.3600e+03,  6.7200e+02,  2.1840e+03,  1.1760e+03,  2.6880e+03,  1.8480e+03,  4.8720e+03,
 0.0000e+00,  7.2000e+00, -4.8000e+01,  1.8480e+02,  6.0000e+01, -7.2000e+01,  1.2000e+01,  0.0000e+00, -8.1000e+01, -1.8540e+02, -1.5900e+02, -5.4600e+01, -4.3200e+03, -1.0800e+03,  1.8144e+04,  6.9120e+03,  1.1232e+04,  3.0240e+03, -1.6128e+04, -7.5600e+03, -7.0560e+03, -2.0160e+03, -2.2680e+04, -9.0720e+03,
 0.0000e+00,  0.0000e+00, -3.0000e+01,  2.1000e+02,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
-1.8000e+01,  6.1200e+01, -1.0200e+02,  5.8800e+01,  6.0000e+01,  3.6000e+01, -2.4000e+01,  0.0000e+00,  0.0000e+00, -9.3600e+01, -3.2400e+02, -3.0240e+02, -3.2400e+03, -2.1600e+03,  9.0720e+03,  7.7760e+03,  1.6416e+04,  9.0720e+03, -6.0480e+03, -6.0480e+03, -1.5120e+04, -8.0640e+03, -2.2176e+04, -1.5120e+04,
 0.0000e+00,  9.3600e+01, -3.2400e+02,  3.0240e+02,  6.0000e+01, -3.6000e+01, -2.4000e+01,  0.0000e+00, -1.8000e+01, -6.1200e+01, -1.0200e+02, -5.8800e+01, -2.1600e+03, -3.2400e+03,  9.0720e+03,  1.6416e+04,  7.7760e+03,  9.0720e+03, -8.0640e+03, -1.5120e+04, -6.0480e+03, -6.0480e+03, -1.5120e+04, -2.2176e+04,
 0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00, -3.0000e+01, -2.1000e+02,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
-8.1000e+01,  1.8540e+02, -1.5900e+02,  5.4600e+01,  6.0000e+01,  7.2000e+01,  1.2000e+01,  0.0000e+00,  0.0000e+00, -7.2000e+00, -4.8000e+01, -1.8480e+02, -1.0800e+03, -4.3200e+03,  3.0240e+03,  1.1232e+04,  6.9120e+03,  1.8144e+04, -2.0160e+03, -7.0560e+03, -7.5600e+03, -1.6128e+04, -9.0720e+03, -2.2680e+04,
 0.0000e+00,  0.0000e+00,  4.2000e+01, -3.7800e+02, -1.0500e+02,  1.6380e+02, -6.3000e+01,  4.2000e+00,  1.1900e+02,  2.3940e+02,  1.6100e+02,  4.0600e+01,  6.0480e+03,  1.5120e+03, -2.8224e+04, -1.1088e+04, -1.3608e+04, -3.5280e+03,  2.6880e+04,  1.3104e+04,  7.5600e+03,  2.0160e+03,  3.1248e+04,  1.2600e+04,
 0.0000e+00,  0.0000e+00,  0.0000e+00, -1.4000e+02,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
 0.0000e+00, -1.0080e+02,  4.2000e+02, -3.1920e+02, -2.1000e+02,  2.5200e+01,  2.1000e+02, -2.5200e+01,  0.0000e+00,  2.7720e+02,  6.3000e+02,  3.5280e+02,  9.0720e+03,  6.0480e+03, -3.1752e+04, -3.0240e+04, -4.0320e+04, -2.1168e+04,  2.4192e+04,  2.7216e+04,  3.3264e+04,  1.6128e+04,  6.9552e+04,  5.1408e+04,
 0.0000e+00,  0.0000e+00,  2.9400e+02, -5.4600e+02, -1.0500e+02,  1.3860e+02, -2.1000e+01, -1.2600e+01,  6.3000e+01,  1.6380e+02,  1.4700e+02,  4.6200e+01,  4.5360e+03,  3.0240e+03, -2.1168e+04, -1.9656e+04, -1.2096e+04, -7.0560e+03,  2.0160e+04,  2.1168e+04,  7.5600e+03,  4.0320e+03,  2.7216e+04,  2.2680e+04,
 6.3000e+01, -1.6380e+02,  1.4700e+02, -4.6200e+01, -1.0500e+02, -1.3860e+02, -2.1000e+01,  1.2600e+01,  0.0000e+00,  0.0000e+00,  2.9400e+02,  5.4600e+02,  3.0240e+03,  4.5360e+03, -7.0560e+03, -1.2096e+04, -1.9656e+04, -2.1168e+04,  4.0320e+03,  7.5600e+03,  2.1168e+04,  2.0160e+04,  2.2680e+04,  2.7216e+04,
 0.0000e+00, -2.7720e+02,  6.3000e+02, -3.5280e+02, -2.1000e+02, -2.5200e+01,  2.1000e+02,  2.5200e+01,  0.0000e+00,  1.0080e+02,  4.2000e+02,  3.1920e+02,  6.0480e+03,  9.0720e+03, -2.1168e+04, -4.0320e+04, -3.0240e+04, -3.1752e+04,  1.6128e+04,  3.3264e+04,  2.7216e+04,  2.4192e+04,  5.1408e+04,  6.9552e+04,
 0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  1.4000e+02,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
 1.1900e+02, -2.3940e+02,  1.6100e+02, -4.0600e+01, -1.0500e+02, -1.6380e+02, -6.3000e+01, -4.2000e+00,  0.0000e+00,  0.0000e+00,  4.2000e+01,  3.7800e+02,  1.5120e+03,  6.0480e+03, -3.5280e+03, -1.3608e+04, -1.1088e+04, -2.8224e+04,  2.0160e+03,  7.5600e+03,  1.3104e+04,  2.6880e+04,  1.2600e+04,  3.1248e+04,
 0.0000e+00,  0.0000e+00,  0.0000e+00,  2.2400e+02,  5.6000e+01, -1.0080e+02,  5.6000e+01, -1.1200e+01, -5.6000e+01, -1.0080e+02, -5.6000e+01, -1.1200e+01, -2.6880e+03, -6.7200e+02,  1.3440e+04,  5.3760e+03,  5.3760e+03,  1.3440e+03, -1.3440e+04, -6.7200e+03, -2.6880e+03, -6.7200e+02, -1.3440e+04, -5.3760e+03,
 0.0000e+00,  0.0000e+00, -3.3600e+02,  3.3600e+02,  1.6800e+02, -1.0080e+02, -1.6800e+02,  1.0080e+02,  0.0000e+00, -2.0160e+02, -3.3600e+02, -1.3440e+02, -6.0480e+03, -4.0320e+03,  2.4192e+04,  2.4192e+04,  2.4192e+04,  1.2096e+04, -2.0160e+04, -2.4192e+04, -1.8144e+04, -8.0640e+03, -4.8384e+04, -3.6288e+04,
 0.0000e+00,  2.0160e+02, -3.3600e+02,  1.3440e+02,  1.6800e+02,  1.0080e+02, -1.6800e+02, -1.0080e+02,  0.0000e+00,  0.0000e+00, -3.3600e+02, -3.3600e+02, -4.0320e+03, -6.0480e+03,  1.2096e+04,  2.4192e+04,  2.4192e+04,  2.4192e+04, -8.0640e+03, -1.8144e+04, -2.4192e+04, -2.0160e+04, -3.6288e+04, -4.8384e+04,
-5.6000e+01,  1.0080e+02, -5.6000e+01,  1.1200e+01,  5.6000e+01,  1.0080e+02,  5.6000e+01,  1.1200e+01,  0.0000e+00,  0.0000e+00,  0.0000e+00, -2.2400e+02, -6.7200e+02, -2.6880e+03,  1.3440e+03,  5.3760e+03,  5.3760e+03,  1.3440e+04, -6.7200e+02, -2.6880e+03, -6.7200e+03, -1.3440e+04, -5.3760e+03, -1.3440e+04
};

static void N_T_RT3_2D_Funct(double xi, double eta, double *values)
{
  int nBF = 24; // number of basis functions
  // monomials x-component and y-component
  double mon_x[24]={1,0,  xi,0,  eta,0,  xi*xi,0,  xi*eta,0,  eta*eta,0,  xi*xi*xi,0,  xi*xi*eta,0,  xi*eta*eta,0,  eta*eta*eta,0,  xi*xi*xi*xi,  xi*xi*xi*eta,  xi*xi*eta*eta,  xi*eta*eta*eta };
  double mon_y[24]={0,1,  0,xi,  0,eta,  0,xi*xi,  0,xi*eta,  0,eta*eta,  0,xi*xi*xi,  0,xi*xi*eta,  0,xi*eta*eta,  0,eta*eta*eta,  xi*xi*xi*eta, xi*xi*eta*eta, xi*eta*eta*eta, eta*eta*eta*eta };
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT3_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT3_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of the derivatives in xi direction
static void N_T_RT3_2D_DeriveXi(double xi, double eta, double *values)
{
  int nBF = 24; // number of basis functions
  // monomials x-component and y-component
  double mon_x[24]={0,0,  1,0,  0,0,  2*xi,0,  eta,0,  0,0,  3*xi*xi,0,  2*xi*eta,0,  eta*eta,0,  0,0,  4*xi*xi*xi,  3*xi*xi*eta,  2*xi*eta*eta, eta*eta*eta};
  double mon_y[24]={0,0,  0,1,  0,0,  0,2*xi,  0,eta,  0,0,  0,3*xi*xi,  0,2*xi*eta,  0,eta*eta,  0,0,  3*xi*xi*eta, 2*xi*eta*eta, eta*eta*eta,  0 };
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT3_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT3_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of the derivatives in eta direction
static void N_T_RT3_2D_DeriveEta(double xi, double eta, double *values)
{
  int nBF = 24; // number of basis functions
  // monomials x-component and y-component
  double mon_x[24]={0,0,  0,0,  1,0,  0,0,  xi,0,  2*eta,0,  0,0,  xi*xi,0,  2*xi*eta,0,  3*eta*eta,0,  0,        xi*xi*xi,    2*xi*xi*eta,  3*xi*eta*eta};
  double mon_y[24]={0,0,  0,0,  0,1,  0,0,  0,xi,  0,2*eta,  0,0,  0,xi*xi,  0,2*xi*eta,  0,3*eta*eta,  xi*xi*xi, 2*xi*xi*eta, 3*xi*eta*eta, 4*eta*eta*eta};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT3_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT3_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in xi-xi direction
static void N_T_RT3_2D_DeriveXiXi(double xi, double eta, double *values)
{
  int nBF = 24; // number of basis functions
  // monomials x-component and y-component
  double mon_x[24]={0,0,  0,0,  0,0,  2,0,  0,0,  0,0,  6*xi,0,  2*eta,0,  0,0,  0,0, 12*xi*xi, 6*xi*eta, 2*eta*eta, 0};
  double mon_y[24]={0,0,  0,0,  0,0,  0,2,  0,0,  0,0,  0,6*xi,  0,2*eta,  0,0,  0,0, 6*xi*eta, 2*eta*eta,0, 0};
 
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT3_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT3_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in eta-eta direction
static void N_T_RT3_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  int nBF = 24; // number of basis functions
  // monomials x-component and y-component
  double mon_x[24]={0,0,  0,0,  0,0,  0,0,  0,0,  2,0,  0,0,  0,0,  2*eta,0,  6*eta,0,   0, 0, 2*xi*xi, 6*xi*eta};
  double mon_y[24]={0,0,  0,0,  0,0,  0,0,  0,0,  0,2,  0,0,  0,0,  0,2*xi,   0,6*eta,0, 2*xi*xi, 6*xi*eta, 12*eta*eta };
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT3_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT3_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in xi-eta direction
static void N_T_RT3_2D_DeriveXiEta(double xi, double eta, double *values)
{
  int nBF = 24; // number of basis functions
  // monomials x-component and y-component
  double mon_x[24]={0,0,  0,0, 0,0,  0,0,  1,0,  0,0,  0,0,  2*xi,0,   2*eta,0,  0,0,  0,3*xi*xi, 4*xi*eta, 3*eta*eta};
  double mon_y[24]={0,0,  0,0, 0,0,  0,0,  0,1,  0,0,  0,0,  0,2*xi,   0,2*eta,  0,0,  3*xi*xi, 4*xi*eta,3*eta*eta, 0};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT3_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT3_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// all dofs N(v) on the edges are integrals of (v.n p) where p is a polynomial.
// For the first and third dof on each edge p is even (p(x) = p(-x)), for the
// second and fourth it is odd (p(x)=-p(-x)). Since the direction of the 
// integration is reversed in the neighbor cell, only the first and third dof on
// each edge have to be changed with TBaseFunct2D::ChangeBF
static int N_T_RT3_2D_ChangeJ0[2] = { 0, 2 };
static int N_T_RT3_2D_ChangeJ1[2] = { 4, 6 };
static int N_T_RT3_2D_ChangeJ2[2] = { 8, 10 };

static int *N_T_RT3_2D_Change1[3] = { N_T_RT3_2D_ChangeJ0, N_T_RT3_2D_ChangeJ1,
                                     N_T_RT3_2D_ChangeJ2};
static int **N_T_RT3_2D_Change[] = { N_T_RT3_2D_Change1 };
