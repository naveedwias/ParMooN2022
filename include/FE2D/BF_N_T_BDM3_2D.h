// Third order Brezzi-Douglas-Marini vector element, nonconforming, 2D

// coefficient matrix for the degrees of freedom
static double N_T_BDM3_2D_CM[400] = {
 0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00, -1.0000e+00, -3.0000e+00, -5.0000e+00, -7.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
-1.0000e+00,  3.0000e+00, -5.0000e+00,  7.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
-2.5000e-01, -3.1250e+00,  3.1250e+00,  1.7500e+01,  6.0000e+00, -6.7500e+00,  0.0000e+00,  1.7500e+00,  1.2250e+01,  2.6875e+01,  2.6875e+01,  1.7500e+01,  1.6500e+02,  7.5000e+01, -1.3875e+02, -4.1250e+01, -1.8000e+02,  8.4000e+02, -1.0500e+03, -1.0500e+03,
 0.0000e+00, -6.0000e+00,  3.0000e+01, -8.4000e+01,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
 0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  6.0000e+00,  3.0000e+01,  8.4000e+01,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
 1.2250e+01, -2.6875e+01,  2.6875e+01, -1.7500e+01,  6.0000e+00,  6.7500e+00,  0.0000e+00, -1.7500e+00, -2.5000e-01,  3.1250e+00,  3.1250e+00, -1.7500e+01,  7.5000e+01,  1.6500e+02, -4.1250e+01, -1.3875e+02, -1.8000e+02, -8.4000e+02,  1.0500e+03,  1.0500e+03,
 0.0000e+00,  2.6250e+00, -2.0625e+01, -6.8250e+01, -1.8750e+01,  3.0000e+01, -1.1250e+01,  0.0000e+00, -3.0000e+01, -5.3625e+01, -3.5625e+01, -1.5750e+01, -4.6500e+02, -1.6500e+02,  4.3875e+02,  7.8750e+01,  4.5000e+02, -1.8900e+03,  3.1500e+03,  1.8900e+03,
 0.0000e+00,  0.0000e+00, -3.0000e+01,  2.1000e+02,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
 0.0000e+00,  1.2750e+01,  1.1250e+01, -3.1500e+01, -2.2500e+01,  0.0000e+00,  2.2500e+01,  0.0000e+00,  0.0000e+00, -5.4750e+01, -1.3875e+02, -1.3650e+02, -2.7000e+02, -2.7000e+02,  2.0250e+02,  2.0250e+02,  5.4000e+02, -3.7800e+03,  3.7800e+03,  6.3000e+03,
 0.0000e+00,  5.4750e+01, -1.3875e+02,  1.3650e+02, -2.2500e+01,  0.0000e+00,  2.2500e+01,  0.0000e+00,  0.0000e+00, -1.2750e+01,  1.1250e+01,  3.1500e+01, -2.7000e+02, -2.7000e+02,  2.0250e+02,  2.0250e+02,  5.4000e+02,  3.7800e+03, -6.3000e+03, -3.7800e+03,
 0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00, -3.0000e+01, -2.1000e+02,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
-3.0000e+01,  5.3625e+01, -3.5625e+01,  1.5750e+01, -1.8750e+01, -3.0000e+01, -1.1250e+01,  0.0000e+00,  0.0000e+00, -2.6250e+00, -2.0625e+01,  6.8250e+01, -1.6500e+02, -4.6500e+02,  7.8750e+01,  4.3875e+02,  4.5000e+02,  1.8900e+03, -1.8900e+03, -3.1500e+03,
 1.2500e+00,  3.5000e+00,  2.2500e+01,  5.7750e+01,  1.3750e+01, -2.6250e+01,  1.6250e+01, -8.7500e+00,  1.8750e+01,  2.9750e+01,  1.3750e+01,  5.2500e+00,  3.0000e+02,  9.0000e+01, -3.0000e+02, -3.7500e+01, -2.7000e+02,  1.0500e+03, -2.1000e+03, -8.4000e+02,
 0.0000e+00,  0.0000e+00,  0.0000e+00, -1.4000e+02,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
-5.0000e+00, -2.3500e+01, -1.2500e+01,  5.6000e+01,  3.0000e+01, -1.5000e+01, -3.0000e+01,  3.5000e+01,  5.0000e+00,  6.6500e+01,  1.0250e+02,  5.6000e+01,  4.2000e+02,  3.0000e+02, -3.4500e+02, -1.9500e+02, -7.2000e+02,  4.2000e+03, -5.8800e+03, -5.8800e+03,
-3.7500e+00, -1.0500e+01,  1.1250e+02, -1.7325e+02,  1.8750e+01, -1.1250e+01, -1.8750e+01,  2.6250e+01,  3.7500e+00,  7.5000e-01, -1.1250e+01, -1.5750e+01,  1.8000e+02,  9.0000e+01, -1.8000e+02, -6.7500e+01, -2.7000e+02, -3.1500e+03,  6.3000e+03,  2.5200e+03,
 3.7500e+00, -7.5000e-01, -1.1250e+01,  1.5750e+01,  1.8750e+01,  1.1250e+01, -1.8750e+01, -2.6250e+01, -3.7500e+00,  1.0500e+01,  1.1250e+02,  1.7325e+02,  9.0000e+01,  1.8000e+02, -6.7500e+01, -1.8000e+02, -2.7000e+02,  3.1500e+03, -2.5200e+03, -6.3000e+03,
 5.0000e+00, -6.6500e+01,  1.0250e+02, -5.6000e+01,  3.0000e+01,  1.5000e+01, -3.0000e+01, -3.5000e+01, -5.0000e+00,  2.3500e+01, -1.2500e+01, -5.6000e+01,  3.0000e+02,  4.2000e+02, -1.9500e+02, -3.4500e+02, -7.2000e+02, -4.2000e+03,  5.8800e+03,  5.8800e+03,
 0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  1.4000e+02,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
 1.8750e+01, -2.9750e+01,  1.3750e+01, -5.2500e+00,  1.3750e+01,  2.6250e+01,  1.6250e+01,  8.7500e+00,  1.2500e+00, -3.5000e+00,  2.2500e+01, -5.7750e+01,  9.0000e+01,  3.0000e+02, -3.7500e+01, -3.0000e+02, -2.7000e+02, -1.0500e+03,  8.4000e+02,  2.1000e+03
};


static void N_T_BDM3_2D_Funct(double xi, double eta, double *values)
{
  int nBF = 20; // number of basis functions
  // monomials x-component and y-component
  double mon_x[20]={1,0,  xi,0,  eta,0,  xi*xi,0,  xi*eta,0,  eta*eta,0,  xi*xi*xi,0,  xi*xi*eta,0,  xi*eta*eta,0,  eta*eta*eta,0 };
  double mon_y[20]={0,1,  0,xi,  0,eta,  0,xi*xi,  0,xi*eta,  0,eta*eta,  0,xi*xi*xi,  0,xi*xi*eta,  0,xi*eta*eta,  0,eta*eta*eta };
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM3_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM3_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of the derivatives in xi direction
static void N_T_BDM3_2D_DeriveXi(double xi, double eta, double *values)
{
  int nBF = 20; // number of basis functions
  // monomials x-component and y-component
  double mon_x[20]={0,0,  1,0,  0,0,  2*xi,0,  eta,0,  0,0,  3*xi*xi,0,  2*xi*eta,0,  eta*eta,0,  0,0 };
  double mon_y[20]={0,0,  0,1,  0,0,  0,2*xi,  0,eta,  0,0,  0,3*xi*xi,  0,2*xi*eta,  0,eta*eta,  0,0 };
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM3_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM3_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of the derivatives in eta direction
static void N_T_BDM3_2D_DeriveEta(double xi, double eta, double *values)
{
  int nBF = 20; // number of basis functions
  // monomials x-component and y-component
  double mon_x[20]={0,0,  0,0,  1,0,  0,0,  xi,0,  2*eta,0,  0,0,  xi*xi,0,  xi*2*eta,0,  3*eta*eta,0 };
  double mon_y[20]={0,0,  0,0,  0,1,  0,0,  0,xi,  0,2*eta,  0,0,  0,xi*xi,  0,xi*2*eta,  0,3*eta*eta };
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM3_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM3_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in xi-xi direction
static void N_T_BDM3_2D_DeriveXiXi(double xi, double eta, double *values)
{
  int nBF = 20; // number of basis functions
  // monomials x-component and y-component
  double mon_x[20]={0,0,  0,0,  0,0,  2,0,  0,0,  0,0,  6*xi,0,  2*eta,0,  0,0,  0,0 };
  double mon_y[20]={0,0,  0,0,  0,0,  0,2,  0,0,  0,0,  0,6*xi,  0,2*eta,  0,0,  0,0 };
 
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM3_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM3_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in eta-eta direction
static void N_T_BDM3_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  int nBF = 20; // number of basis functions
  // monomials x-component and y-component
  double mon_x[20]={0,0,  0,0,  0,0,  0,0,  0,0,  2,0,  0,0,  0,0,  xi*2,0,  6*eta,0 };
  double mon_y[20]={0,0,  0,0,  0,0,  0,0,  0,0,  0,2,  0,0,  0,0,  0,xi*2,  0,6*eta };
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM3_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM3_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in xi-eta direction
static void N_T_BDM3_2D_DeriveXiEta(double xi, double eta, double *values)
{
  int nBF = 20; // number of basis functions
  // monomials x-component and y-component
  double mon_x[20]={0,0,  0,0,  0,0,  0,0,  1,0,  0,0,  0,0,  2*xi,0,  2*eta,0,  0,0 };
  double mon_y[20]={0,0,  0,0,  0,0,  0,0,  0,1,  0,0,  0,0,  0,2*xi,  0,2*eta,  0,0 };
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM3_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM3_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// ***********************************************************************

// all dofs N(v) on the edges are integrals of (v.n p) where p is a polynomial.
// For the first and third dof on each edge p is even (p(x) = p(-x)), for the
// second and fourth it is odd (p(x)=-p(-x)). Since the direction of the 
// integration is reversed in the neighbor cell, only the first and third dof on
// each edge have to be changed with TBaseFunct2D::ChangeBF
static int N_T_BDM3_2D_ChangeJ0[2] = { 0, 2 };
static int N_T_BDM3_2D_ChangeJ1[2] = { 4, 6 };
static int N_T_BDM3_2D_ChangeJ2[2] = { 8, 10 };

static int *N_T_BDM3_2D_Change1[3] = {N_T_BDM3_2D_ChangeJ0, N_T_BDM3_2D_ChangeJ1,
                                     N_T_BDM3_2D_ChangeJ2};
static int **N_T_BDM3_2D_Change[] = { N_T_BDM3_2D_Change1 };
