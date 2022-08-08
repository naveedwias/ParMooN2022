// Second order Raviart-Thomas vector element on quads, nonconforming, 2D

// coefficient matrix for the degrees of freedom (this is in fact 32 times that 
// matrix)
static double N_Q_RT2_2D_CM[576] = {
  0.,   0.,   0.,  -4.,   0.,  10.,   0.,   0.,   0.,   4.,   0., -10.,  27.,   0.,   0.,   0.,   0.,   0.,   0., -45.,   0.,   0.,   0.,   0.,
  4.,   0., -10.,   0.,   0.,   0.,  -4.,   0.,  10.,   0.,   0.,   0.,   0.,  27.,   0.,   0.,   0.,   0., -45.,   0.,   0.,   0.,   0.,   0.,
  0.,   0.,   0., -12.,   0.,  30.,   0.,   0.,   0., -12.,   0.,  30.,   0.,   0., 135.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,-225.,
  0.,  12.,   0.,   0.,   0.,   0.,   0.,  12.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  36.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
  0.,   0.,   0.,   0., -12.,   0.,   0.,   0.,   0.,   0., -12.,   0.,   0.,   0.,   0.,   0.,  36.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
-12.,   0.,  30.,   0.,   0.,   0., -12.,   0.,  30.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 135.,   0.,   0.,   0.,   0.,-225.,   0.,
  0.,   0.,   0.,  12.,   0., -30.,   0.,   0.,   0., -12.,   0.,  30., -27.,   0.,   0.,   0.,   0.,   0.,   0.,  45.,   0.,   0.,   0.,   0.,
  0.,   0.,  30.,   0.,   0.,   0.,   0.,   0., -30.,   0.,   0.,   0.,   0., -45.,   0.,   0.,   0.,   0., 135.,   0.,   0.,   0.,   0.,   0.,
  0.,   0.,   0.,   0., -36.,   0.,   0.,   0.,   0.,   0.,  36.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 180.,   0.,   0.,   0.,
  0., -36.,   0.,   0.,   0.,   0.,   0.,  36.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 180.,   0.,   0.,
  0.,   0.,   0.,   0.,   0., -30.,   0.,   0.,   0.,   0.,   0.,  30., -45.,   0.,   0.,   0.,   0.,   0.,   0., 135.,   0.,   0.,   0.,   0.,
-12.,   0.,  30.,   0.,   0.,   0.,  12.,   0., -30.,   0.,   0.,   0.,   0., -27.,   0.,   0.,   0.,   0.,  45.,   0.,   0.,   0.,   0.,   0.,
  0.,   0.,   0.,  20.,   0., -50.,   0.,   0.,   0.,  20.,   0., -50.,   0.,   0.,-135.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 225.,
 20.,   0., -50.,   0.,   0.,   0.,  20.,   0., -50.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,-135.,   0.,   0.,   0.,   0., 225.,   0.,
  0.,   0.,   0.,   0.,  36.,   0.,   0.,   0.,   0.,   0.,  36.,   0.,   0.,   0.,   0.,   0., -36.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
  0.,   0., -90.,   0.,   0.,   0.,   0.,   0., -90.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,-225.,   0.,   0.,   0.,   0., 675.,   0.,
  0.,   0.,   0.,   0.,   0., -90.,   0.,   0.,   0.,   0.,   0., -90.,   0.,   0.,-225.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 675.,
  0., -36.,   0.,   0.,   0.,   0.,   0., -36.,   0.,   0.,   0.,   0.,   0.,   0.,   0., -36.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
  0.,   0.,   0.,   0.,  60.,   0.,   0.,   0.,   0.,   0., -60.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,-180.,   0.,   0.,   0.,
  0.,  60.,   0.,   0.,   0.,   0.,   0., -60.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,-180.,   0.,   0.,
  0.,   0.,   0.,   0.,   0.,  90.,   0.,   0.,   0.,   0.,   0., -90.,  45.,   0.,   0.,   0.,   0.,   0.,   0.,-135.,   0.,   0.,   0.,   0.,
  0.,   0., -90.,   0.,   0.,   0.,   0.,   0.,  90.,   0.,   0.,   0.,   0.,  45.,   0.,   0.,   0.,   0.,-135.,   0.,   0.,   0.,   0.,   0.,
  0.,   0.,   0.,   0.,   0., 150.,   0.,   0.,   0.,   0.,   0., 150.,   0.,   0., 225.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,-675.,
  0.,   0., 150.,   0.,   0.,   0.,   0.,   0., 150.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 225.,   0.,   0.,   0.,   0.,-675.,   0.
};

static void N_Q_RT2_2D_Funct(double xi, double eta, double *values)
{
  int nBF = 24; // number of basis functions
  // monomials x-component and y-component
  double mon_x[24]={1,0, xi,0, eta,0, xi*xi,0, xi*eta,0, eta*eta,0, xi*xi*xi,0,    xi*xi*eta,0, xi*eta*eta,0, xi*xi*xi*eta,0,   xi*xi*eta*eta,0, xi*xi*xi*eta*eta,0};
  double mon_y[24]={0,1, 0,xi, 0,eta, 0,xi*xi, 0,xi*eta, 0,eta*eta, 0,eta*eta*eta, 0,xi*xi*eta, 0,xi*eta*eta, 0,xi*eta*eta*eta, 0,xi*xi*eta*eta, 0,xi*xi*eta*eta*eta};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_RT2_2D_CM[i+j*nBF]*mon_x[j] / 32.;
      values[i+nBF] += N_Q_RT2_2D_CM[i+j*nBF]*mon_y[j] / 32.;
    }
  }
}

// values of the derivatives in xi direction
static void N_Q_RT2_2D_DeriveXi(double xi, double eta, double *values)
{
  int nBF = 24; // number of basis functions
  // monomials x-component and y-component
  double mon_x[24]={0,0, 1,0, 0,0, 2*xi,0, eta,0, 0,0, 3*xi*xi,0, 2*xi*eta,0, eta*eta, 0,3*xi*xi*eta,0,2*xi*eta*eta,0,3*xi*xi*eta*eta,0};
  double mon_y[24]={0,0, 0,1, 0,0, 0,2*xi, 0,eta, 0,0, 0,      0, 0,2*xi*eta, 0,eta*eta,0, eta*eta*eta,0,2*xi*eta*eta,0,2*xi*eta*eta*eta};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_RT2_2D_CM[i+j*nBF]*mon_x[j] / 32.;
      values[i+nBF] += N_Q_RT2_2D_CM[i+j*nBF]*mon_y[j] / 32.;
    }
  }
}

// values of the derivatives in eta direction
static void N_Q_RT2_2D_DeriveEta(double xi, double eta, double *values)
{
  int nBF = 24; // number of basis functions
  // monomials x-component and y-component
  double mon_x[24]={0,0, 0,0, 1,0, 0,0, xi,0, 2*eta,0, 0,        0, xi*xi,0, 2*xi*eta,0, xi*xi*xi,    0, 2*xi*xi*eta,0, 2*xi*xi*xi*eta,0};
  double mon_y[24]={0,0, 0,0, 0,1, 0,0, 0,xi, 0,2*eta, 0,3*eta*eta, 0,xi*xi, 0,2*xi*eta, 0,3*xi*eta*eta, 0,2*xi*xi*eta, 0,3*xi*xi*eta*eta};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_RT2_2D_CM[i+j*nBF]*mon_x[j] / 32.;
      values[i+nBF] += N_Q_RT2_2D_CM[i+j*nBF]*mon_y[j] / 32.;
    }
  }
}

// values of derivatives in xi-xi direction
static void N_Q_RT2_2D_DeriveXiXi(double xi, double eta, double *values)
{
  int nBF = 24; // number of basis functions
  // monomials x-component and y-component
  double mon_x[24]={0,0, 0,0, 0,0, 2,0, 0,0, 0,0, 6*xi,0, 2*eta,0, 0,0, 6*xi*eta,0, 2*eta*eta,0, 6*xi*eta*eta,0};
  double mon_y[24]={0,0, 0,0, 0,0, 0,2, 0,0, 0,0, 0,   0, 0,2*eta, 0,0, 0,       0, 0,2*eta*eta, 0,2*eta*eta*eta};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_RT2_2D_CM[i+j*nBF]*mon_x[j] / 32.;
      values[i+nBF] += N_Q_RT2_2D_CM[i+j*nBF]*mon_y[j] / 32.;
    }
  }
}

// values of derivatives in eta-eta direction
static void N_Q_RT2_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  int nBF = 24; // number of basis functions
  // monomials x-component and y-component
  double mon_x[24]={0,0, 0,0, 0,0, 0,0, 2,0, 0,0, 0,    0, 0,0, 2*xi,0, 0,       0, 2*xi*xi,0, 2*xi*xi*xi,0};
  double mon_y[24]={0,0, 0,0, 0,0, 0,0, 0,2, 0,0, 0,6*eta, 0,0, 0,2*xi, 0,6*xi*eta, 0,2*xi*xi, 0,6*xi*xi*eta};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_RT2_2D_CM[i+j*nBF]*mon_x[j] / 32.;
      values[i+nBF] += N_Q_RT2_2D_CM[i+j*nBF]*mon_y[j] / 32.;
    }
  }
}

// values of derivatives in xi-eta direction
static void N_Q_RT2_2D_DeriveXiEta(double xi, double eta, double *values)
{
  int nBF = 24; // number of basis functions
  // monomials x-component and y-component
  double mon_x[24]={0,0, 0,0, 0,0, 0,0, 1,0, 0,0, 0,0, 2*xi,0, 2*eta,0, 3*xi*xi,  0, 4*xi*eta,0, 6*xi*xi*eta,0};
  double mon_y[24]={0,0, 0,0, 0,0, 0,0, 0,1, 0,0, 0,0, 0,2*xi, 0,2*eta, 0,3*eta*eta, 0,4*xi*eta, 0,6*xi*eta*eta};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_RT2_2D_CM[i+j*nBF]*mon_x[j] / 32.;
      values[i+nBF] += N_Q_RT2_2D_CM[i+j*nBF]*mon_y[j] / 32.;
    }
  }
}

// all dofs N(v) on the edges are integrals of (v.n p) where p is a polynomial.
// For the first and third dof on each edge p is even (p(x) = p(-x)), for the
// second it is odd (p(x)=-p(-x)). Since the direction of the integration is
// reversed in the neighbor cell, only the first and third dof on each edge have
// to be changed with TBaseFunct2D::ChangeBF
static int N_Q_RT2_2D_ChangeJ0[2] = { 0,2 };
static int N_Q_RT2_2D_ChangeJ1[2] = { 3,5 };
static int N_Q_RT2_2D_ChangeJ2[2] = { 6,8 };
static int N_Q_RT2_2D_ChangeJ3[2] = { 9,11 };

static int *N_Q_RT2_2D_Change1[4] = { N_Q_RT2_2D_ChangeJ0, N_Q_RT2_2D_ChangeJ1,
                                     N_Q_RT2_2D_ChangeJ2, N_Q_RT2_2D_ChangeJ3 };
static int **N_Q_RT2_2D_Change[] = { N_Q_RT2_2D_Change1 };
