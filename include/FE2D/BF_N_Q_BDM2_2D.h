// Second order Brezzi-Douglas-Marini vector element, nonconforming, 2D

// coefficient matrix for the degrees of freedom (this is in fact 8 times that 
// matrix)
static double N_Q_BDM2_2D_CM[196] = {
 0.,  0.,  0., -1.,  0., -5.,  0.,  0.,  0.,  1.,  0.,  5.,  3.,  0.,
 1.,  0.,  5.,  0.,  0.,  0., -1.,  0., -5.,  0.,  0.,  0.,  0.,  3.,
 0.,  0.,  5.,  2.,  0., -5.,  0.,  0.,  5.,  2.,  0., -5.,  0.,  0.,
 0., -6.,  0.,  0.,  0.,  0.,  0., -6.,  0.,  0.,  0.,  0.,  0.,  0.,
 0.,  0.,  0.,  0.,  6.,  0.,  0.,  0.,  0.,  0.,  6.,  0.,  0.,  0.,
 2.,  0., -5.,  0.,  0.,  5.,  2.,  0., -5.,  0.,  0.,  5.,  0.,  0.,
 0.,  0.,  0.,  3.,  0.,  0.,  0.,  0.,  0., -3.,  0.,  0., -3.,  0.,
 0.,  0.,-15.,  0.,  0.,  0.,  0.,  0., 15.,  0.,  0.,  0.,  0.,  0.,
 0.,  0.,  0.,  0.,  6.,  0.,  0.,  0.,  0.,  0., -6.,  0.,  0.,  0.,
 0.,  6.,  0.,  0.,  0.,  0.,  0., -6.,  0.,  0.,  0.,  0.,  0.,  0.,
 0.,  0.,  0.,  0.,  0., 15.,  0.,  0.,  0.,  0.,  0.,-15.,  0.,  0.,
-3.,  0.,  0.,  0.,  0.,  0.,  3.,  0.,  0.,  0.,  0.,  0.,  0., -3.,
 0.,  0., -5.,  0.,  0.,  0.,  0.,  0., -5.,  0.,  0.,  0.,  0.,  0.,
 0.,  0.,  0.,  0.,  0.,  5.,  0.,  0.,  0.,  0.,  0.,  5.,  0.,  0.
};

static void N_Q_BDM2_2D_Funct(double xi, double eta, double *values)
{
  int nBF = 14; // number of basis functions
  // monomials x-component and y-component
  double mon_x[14]={1,0,  xi,0,  eta,0,  xi*xi,0,  xi*eta,0,  eta*eta,0,  xi*xi*xi,3*xi*eta*eta};
  double mon_y[14]={0,1,  0,xi,  0,eta,  0,xi*xi,  0,xi*eta,  0,eta*eta,  -3*xi*xi*eta,-eta*eta*eta};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_x[j] / 8.;
      values[i+nBF] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_y[j] / 8.;
    }
  }
}

// values of the derivatives in xi direction
static void N_Q_BDM2_2D_DeriveXi(double xi, double eta, double *values)
{
  int nBF = 14; // number of basis functions
  // monomials x-component and y-component
  double mon_x[14]={0,0,  1,0,  0,0,  2*xi,0, eta,0,  0,0,   3*xi*xi,3*eta*eta};
  double mon_y[14]={0,0,  0,1,  0,0,  0,2*xi, 0,eta,  0,0,  -6*xi*eta,0      };
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_x[j] / 8.;
      values[i+nBF] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_y[j] / 8.;
    }
  }
}

// values of the derivatives in eta direction
static void N_Q_BDM2_2D_DeriveEta(double xi, double eta, double *values)
{
  int nBF = 14; // number of basis functions
  // monomials x-component and y-component
  double mon_x[14]={0,0,  0,0,  1,0,  0,0,  xi,0,  2*eta,0,  0,6*xi*eta};
  double mon_y[14]={0,0,  0,0,  0,1,  0,0,  0,xi,  0,2*eta, -3*xi*xi,-3*eta*eta};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_x[j] / 8.;
      values[i+nBF] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_y[j] / 8.;
    }
  }
}

// values of derivatives in xi-xi direction
static void N_Q_BDM2_2D_DeriveXiXi(double xi, double eta, double *values)
{
  int nBF = 14; // number of basis functions
  // monomials x-component and y-component
  double mon_x[14]={0,0, 0,0, 0,0, 2,0, 0,0, 0,0, 6*xi,0};
  double mon_y[14]={0,0, 0,0, 0,0, 0,2, 0,0, 0,0, -6*eta,0};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_x[j] / 8.;
      values[i+nBF] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_y[j] / 8.;
    }
  }
}

// values of derivatives in eta-eta direction
static void N_Q_BDM2_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  int nBF = 14; // number of basis functions
  // monomials x-component and y-component
  double mon_x[14]={0,0, 0,0, 0,0, 0,0, 0,0, 2,0, 0,6*xi};
  double mon_y[14]={0,0, 0,0, 0,0, 0,0, 0,0, 0,2, 0,-6*eta};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_x[j] / 8.;
      values[i+nBF] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_y[j] / 8.;
    }
  }
}

// values of derivatives in xi-eta direction
static void N_Q_BDM2_2D_DeriveXiEta(double xi, double eta, double *values)
{
  int nBF = 14; // number of basis functions
  // monomials x-component and y-component
  double mon_x[14]={0,0, 0,0, 0,0, 0,0, 1,0, 0,0,  0,6*eta};
  double mon_y[14]={0,0, 0,0, 0,0, 0,0, 0,1, 0,0, -6*xi,0};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_x[j] / 8.;
      values[i+nBF] += N_Q_BDM2_2D_CM[i+j*nBF]*mon_y[j] / 8.;
    }
  }
}

// all dofs N(v) on the edges are integrals of (v.n p) where p is a polynomial.
// For the first and third dof on each edge p is even (p(x) = p(-x)), for the
// second it is odd (p(x)=-p(-x)). Since the direction of the integration is
// reversed in the neighbor cell, only the first and third dof on each edge have
// to be changed with TBaseFunct2D::ChangeBF
static int N_Q_BDM2_2D_ChangeJ0[2] = { 0,2 };
static int N_Q_BDM2_2D_ChangeJ1[2] = { 3,5 };
static int N_Q_BDM2_2D_ChangeJ2[2] = { 6,8 };
static int N_Q_BDM2_2D_ChangeJ3[2] = { 9,11 };

static int *N_Q_BDM2_2D_Change1[4] = {N_Q_BDM2_2D_ChangeJ0,N_Q_BDM2_2D_ChangeJ1,
                                     N_Q_BDM2_2D_ChangeJ2,N_Q_BDM2_2D_ChangeJ3};
static int **N_Q_BDM2_2D_Change[] = { N_Q_BDM2_2D_Change1 };
