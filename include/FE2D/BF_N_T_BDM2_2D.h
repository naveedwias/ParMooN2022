// Second order Brezzi-Douglas-Marini vector element, nonconforming, 2D

// coefficient matrix for the degrees of freedom
static double N_T_BDM2_2D_CM[144] = {
 0.,   0. ,   0. ,   0.,   0.,   0.,  -1.,  -3. ,  -5. ,   0.,   0.,    0.,
-1.,   3. ,  -5. ,   0.,   0.,   0.,   0.,   0. ,   0. ,   0.,   0.,    0.,
 0.,  -1.5,  -7.5,  -3.,   3.,   0.,   6.,  10.5,   7.5,  18.,   6.,   90.,
 0.,  -6. ,  30. ,   0.,   0.,   0.,   0.,   0. ,   0. ,   0.,   0.,    0.,
 0.,   0. ,   0. ,   0.,   0.,   0.,   0.,   6. ,  30. ,   0.,   0.,    0.,
 6., -10.5,   7.5,  -3.,  -3.,   0.,   0.,   1.5,  -7.5,   6.,  18.,  -90.,
 1.,   4.5,  12.5,   4.,  -6.,   5.,  -5.,  -7.5,  -2.5, -18.,  -6.,  -90.,
 0.,   0. , -30. ,   0.,   0.,   0.,   0.,   0. ,   0. ,   0.,   0.,    0.,
-2.,  -3. ,   5. ,   4.,   0., -10.,  -2., -15. , -25. , -12., -12., -180.,
-2.,  15. , -25. ,   4.,   0., -10.,  -2.,   3. ,   5. , -12., -12.,  180.,
 0.,   0. ,   0. ,   0.,   0.,   0.,   0.,   0. , -30. ,   0.,   0.,    0.,
-5.,   7.5,  -2.5,   4.,   6.,   5.,   1.,  -4.5,  12.5,  -6., -18.,   90.
};

static void N_T_BDM2_2D_Funct(double xi, double eta, double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component and y-component
  double mon_x[12]={1,0,  xi,0,  eta,0,  xi*xi,0,  xi*eta,0,  eta*eta,0};
  double mon_y[12]={0,1,  0,xi,  0,eta,  0,xi*xi,  0,xi*eta,  0,eta*eta};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of the derivatives in xi direction
static void N_T_BDM2_2D_DeriveXi(double xi, double eta, double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component and y-component
  double mon_x[12]={0,0,  1,0,  0,0,  2*xi,0,  eta,0,  0,0};
  double mon_y[12]={0,0,  0,1,  0,0,  0,2*xi,  0,eta,  0,0};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of the derivatives in eta direction
static void N_T_BDM2_2D_DeriveEta(double xi, double eta, double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component and y-component
  double mon_x[12]={0,0,  0,0,  1,0,  0,0,  xi,0,  2*eta,0};
  double mon_y[12]={0,0,  0,0,  0,1,  0,0,  0,xi,  0,2*eta};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in xi-xi direction
static void N_T_BDM2_2D_DeriveXiXi(double, double, double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component and y-component
  double mon_x[12]={0,0,  0,0,  0,0,  2,0,  0,0,  0,0};
  double mon_y[12]={0,0,  0,0,  0,0,  0,2,  0,0,  0,0};
 
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in eta-eta direction
static void N_T_BDM2_2D_DeriveEtaEta(double, double, double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component and y-component
  double mon_x[12]={0,0,  0,0,  0,0,  0,0,  0,0,  2,0};
  double mon_y[12]={0,0,  0,0,  0,0,  0,0,  0,0,  0,2};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// values of derivatives in xi-eta direction
static void N_T_BDM2_2D_DeriveXiEta(double, double, double *values)
{
  int nBF = 12; // number of basis functions
  // monomials x-component and y-component
  double mon_x[12]={0,0,  0,0,  0,0,  0,0,  1,0,  0,0};
  double mon_y[12]={0,0,  0,0,  0,0,  0,0,  0,1,  0,0};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_T_BDM2_2D_CM[i+j*nBF]*mon_x[j];
      values[i+nBF] += N_T_BDM2_2D_CM[i+j*nBF]*mon_y[j];
    }
  }
}

// ***********************************************************************

// all dofs N(v) on the edges are integrals of (v.n p) where p is a polynomial.
// For the first and third dof on each edge p is even (p(x) = p(-x)), for the
// second it is odd (p(x)=-p(-x)). Since the direction of the integration is
// reversed in the neighbor cell, only the first and third dof on each edge have
// to be changed with TBaseFunct2D::ChangeBF
static int N_T_BDM2_2D_ChangeJ0[2] = { 0,2 };
static int N_T_BDM2_2D_ChangeJ1[2] = { 3,5 };
static int N_T_BDM2_2D_ChangeJ2[2] = { 6,8 };

static int *N_T_BDM2_2D_Change1[3] = {N_T_BDM2_2D_ChangeJ0, N_T_BDM2_2D_ChangeJ1,
                                     N_T_BDM2_2D_ChangeJ2 };
static int **N_T_BDM2_2D_Change[] = { N_T_BDM2_2D_Change1 };
