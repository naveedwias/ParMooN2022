// First order Raviart-Thomas vector element, nonconforming, 2D

// coefficient matrix for the degrees of freedom
static double N_T_RT1_2D_CM[64] = {
 0.,  0.,  0.,  0., -1., -3.,  0.,  0.,
-1.,  3.,  0.,  0.,  0.,  0.,  0.,  0.,
 1., -5., -3.,  1.,  5.,  7., 16.,  8.,
 0., -6.,  0.,  0.,  0.,  0.,  0.,  0.,
 0.,  0.,  0.,  0.,  0.,  6.,  0.,  0.,
 5., -7., -3., -1.,  1.,  5.,  8., 16.,
 0.,  8.,  4., -4., -4., -4.,-16., -8.,
-4.,  4.,  4.,  4.,  0., -8., -8.,-16.
};
static void N_T_RT1_2D_Funct(double xi, double eta, double *values)
{
  int nBF = 8; // number of basis functions
  // monomials x-component and y-component
  double mon_x[]={1,0,xi,0 ,eta,0  ,xi*xi ,xi*eta };
  double mon_y[]={0,1,0 ,xi,0  ,eta,xi*eta,eta*eta};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of the derivatives in xi direction
static void N_T_RT1_2D_DeriveXi(double xi, double eta, double *values)
{
  int nBF = 8; // number of basis functions
  // monomials x-component and y-component
  double mon_x[]={0,0,1,0,0,0,2*xi,eta };
  double mon_y[]={0,0,0,1,0,0,eta ,0};
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of the derivatives in eta direction
static void N_T_RT1_2D_DeriveEta(double xi, double eta, double *values)
{
  int nBF = 8; // number of basis functions
  // monomials x-component and y-component
  double mon_x[]={0,0,0,0,1,0,0 ,xi };
  double mon_y[]={0,0,0,0,0,1,xi,2*eta};
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in xi-xi direction
static void N_T_RT1_2D_DeriveXiXi(double, double, double *values)
{
  int nBF = 8; // number of basis functions
  // monomials x-component and y-component
  double mon_x[]={0,0,0,0,0,0,2,0};
  double mon_y[]={0,0,0,0,0,0,0,0};
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in eta-eta direction
static void N_T_RT1_2D_DeriveEtaEta(double, double, double *values)
{
  int nBF = 8; // number of basis functions
  // monomials x-component and y-component
  double mon_x[]={0,0,0,0,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,2};
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in xi-eta direction
static void N_T_RT1_2D_DeriveXiEta(double, double, double *values)
{
  int nBF = 8; // number of basis functions
  // monomials x-component and y-component
  double mon_x[]={0,0,0,0,0,0,0,1};
  double mon_y[]={0,0,0,0,0,0,1,0};
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_RT1_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_RT1_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// all dofs N(v) on the edges are integrals of (v.n p) where p is a polynomial.
// For the first dof on each edge p is even (p(x) = p(-x)), for the second it is
// odd (p(x)=-p(-x)). Since the direction of the integration is reversed in the
// neighbor cell, only the first dof on each edge has to be changed with
// TBaseFunct2D::ChangeBF
static int N_T_RT1_2D_ChangeJ0[1] = { 0 };
static int N_T_RT1_2D_ChangeJ1[1] = { 2 };
static int N_T_RT1_2D_ChangeJ2[1] = { 4 };

static int *N_T_RT1_2D_Change1[3] = { N_T_RT1_2D_ChangeJ0, N_T_RT1_2D_ChangeJ1,
                                     N_T_RT1_2D_ChangeJ2};
static int **N_T_RT1_2D_Change[] = { N_T_RT1_2D_Change1 };
