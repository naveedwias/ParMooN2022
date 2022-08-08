// Third order Raviart-Thomas vector element on quads, nonconforming, 2D

// coefficient matrix for the degrees of freedom (this is in fact 768 times 
// that matrix)
static double N_Q_RT3_2D_CM[1600] = {
   0.,       0.,       0.,       0.,      72.,       0.,    -180.,       0.,       0.,       0.,       0.,       0.,     -72.,       0.,     180.,       0.,    1215.,       0.,       0.,       0.,       0.,     0.,   -2835.,       0.,   -2025.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,      -0.,       0.,    4725.,       0.,       0.,       0.,       0.,       0.,
 -72.,       0.,     180.,       0.,       0.,       0.,       0.,       0.,      72.,       0.,    -180.,       0.,       0.,       0.,       0.,       0.,       0.,    1215.,       0.,       0.,       0.,     0.,       0.,   -2025.,       0.,   -2835.,       0.,       0.,       0.,       0.,       0.,      -0.,       0.,       0.,       0.,    4725.,       0.,       0.,       0.,       0.,
   0.,       0.,       0.,       0.,    -288.,       0.,     720.,       0.,       0.,       0.,       0.,       0.,    -288.,       0.,     720.,       0.,       0.,       0.,    3240.,       0.,       0.,     0.,      -0.,       0.,      -0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -5400.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,
   0.,    -216.,       0.,     756.,       0.,       0.,       0.,       0.,       0.,    -216.,       0.,     756.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   10125.,       0.,     0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,  -14175.,       0.,       0.,       0.,  -23625.,       0.,       0.,       0.,       0.,       0.,   33075.,
   0.,       0.,       0.,       0.,       0.,     216.,       0.,    -756.,       0.,       0.,       0.,       0.,       0.,     216.,       0.,    -756.,       0.,       0.,       0.,       0.,   10125.,     0.,       0.,       0.,       0.,       0.,       0.,       0.,  -14175.,       0.,  -23625.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   33075.,       0.,
-288.,       0.,     720.,       0.,       0.,       0.,       0.,       0.,    -288.,       0.,     720.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,  3240.,       0.,      -0.,       0.,      -0.,       0.,       0.,       0.,       0.,       0.,   -5400.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,
   0.,       0.,       0.,       0.,    -720.,       0.,    1800.,       0.,       0.,       0.,       0.,       0.,     720.,       0.,   -1800.,       0.,   -4050.,       0.,      -0.,       0.,       0.,     0.,   17010.,       0.,    6750.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,  -28350.,       0.,       0.,       0.,       0.,       0.,
   0.,       0.,    -540.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     540.,       0.,       0.,       0.,       0.,       0.,       0.,   -2025.,       0.,       0.,       0.,     0.,       0.,    6075.,       0.,    4725.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,  -14175.,       0.,       0.,       0.,       0.,
   0.,       0.,       0.,       0.,       0.,    -864.,       0.,    3024.,       0.,       0.,       0.,       0.,       0.,     864.,       0.,   -3024.,       0.,       0.,       0.,       0.,       0.,     0.,       0.,       0.,       0.,       0.,   27000.,       0.,       0.,       0.,      -0.,       0.,       0.,       0.,       0.,       0.,  -37800.,       0.,       0.,       0.,
   0.,    -864.,       0.,    3024.,       0.,       0.,       0.,       0.,       0.,     864.,       0.,   -3024.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,      -0.,       0.,     0.,       0.,       0.,       0.,       0.,       0.,   27000.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,  -37800.,       0.,       0.,
   0.,       0.,       0.,       0.,       0.,       0.,     540.,       0.,       0.,       0.,       0.,       0.,      -0.,       0.,    -540.,       0.,   -2025.,       0.,      -0.,       0.,       0.,     0.,    4725.,       0.,    6075.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,  -14175.,       0.,       0.,       0.,       0.,       0.,
 720.,       0.,   -1800.,       0.,       0.,       0.,       0.,       0.,    -720.,       0.,    1800.,       0.,       0.,       0.,       0.,       0.,       0.,   -4050.,       0.,       0.,       0.,     0.,       0.,    6750.,       0.,   17010.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,  -28350.,       0.,       0.,       0.,       0.,
   0.,       0.,       0.,       0.,     480.,       0.,   -1200.,       0.,       0.,       0.,       0.,       0.,     480.,       0.,   -1200.,       0.,      -0.,       0.,   -3240.,       0.,       0.,     0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    5400.,       0.,      -0.,       0.,       0.,       0.,       0.,       0.,
   0.,       0.,       0.,   -1260.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   -1260.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,  -14175.,       0.,     0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   23625.,       0.,       0.,       0.,   33075.,       0.,       0.,       0.,       0.,       0.,  -55125.,
   0.,       0.,       0.,       0.,       0.,   -2160.,       0.,    7560.,       0.,       0.,       0.,       0.,       0.,   -2160.,       0.,    7560.,       0.,       0.,       0.,       0.,  -33750.,     0.,       0.,       0.,       0.,       0.,       0.,       0.,   47250.,       0.,  141750.,       0.,       0.,       0.,       0.,       0.,       0.,       0., -198450.,       0.,
  -0.,       0.,   -2160.,       0.,       0.,       0.,       0.,       0.,      -0.,       0.,   -2160.,       0.,       0.,       0.,       0.,       0.,       0.,      -0.,       0.,       0.,       0., -5400.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   16200.,       0.,       0.,       0.,      -0.,       0.,       0.,       0.,       0.,
   0.,       0.,       0.,       0.,      -0.,       0.,   -2160.,       0.,       0.,       0.,       0.,       0.,      -0.,       0.,   -2160.,       0.,      -0.,       0.,   -5400.,       0.,       0.,     0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   16200.,       0.,      -0.,       0.,       0.,       0.,       0.,       0.,
   0.,    2160.,       0.,   -7560.,       0.,       0.,       0.,       0.,       0.,    2160.,       0.,   -7560.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,  -33750.,       0.,     0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   47250.,       0.,       0.,       0.,  141750.,       0.,       0.,       0.,       0.,       0., -198450.,
   0.,       0.,       0.,       0.,       0.,      -0.,       0.,    1260.,       0.,       0.,       0.,       0.,       0.,      -0.,       0.,    1260.,       0.,       0.,       0.,       0.,  -14175.,     0.,       0.,       0.,       0.,       0.,       0.,       0.,   23625.,       0.,   33075.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,  -55125.,       0.,
 480.,       0.,   -1200.,       0.,       0.,       0.,       0.,       0.,     480.,       0.,   -1200.,       0.,       0.,       0.,       0.,       0.,       0.,      -0.,       0.,       0.,       0., -3240.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    5400.,       0.,       0.,       0.,      -0.,       0.,       0.,       0.,       0.,
   0.,       0.,       0.,       0.,     840.,       0.,   -2100.,       0.,       0.,       0.,       0.,       0.,    -840.,       0.,    2100.,       0.,    2835.,       0.,       0.,       0.,       0.,     0.,  -14175.,       0.,   -4725.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,      -0.,       0.,   23625.,       0.,       0.,       0.,       0.,       0.,
   0.,      -0.,       0.,   -5040.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    5040.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.,       0.,       0.,       0.,       0.,       0.,  -37800.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   63000.,       0.,       0.,
   0.,       0.,       0.,       0.,       0.,    1440.,       0.,   -5040.,       0.,       0.,       0.,       0.,       0.,   -1440.,       0.,    5040.,       0.,       0.,       0.,       0.,       0.,     0.,       0.,       0.,       0.,       0.,  -27000.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   37800.,       0.,       0.,       0.,
   0.,       0.,    5400.,       0.,       0.,       0.,       0.,       0.,      -0.,       0.,   -5400.,       0.,       0.,       0.,       0.,       0.,       0.,    6750.,       0.,       0.,       0.,     0.,       0.,  -20250.,       0.,  -28350.,       0.,       0.,       0.,       0.,       0.,      -0.,       0.,       0.,       0.,   85050.,       0.,       0.,       0.,       0.,
   0.,       0.,       0.,       0.,      -0.,       0.,   -5400.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    5400.,       0.,    6750.,       0.,       0.,       0.,       0.,     0.,  -28350.,       0.,  -20250.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,      -0.,       0.,   85050.,       0.,       0.,       0.,       0.,       0.,
   0.,    1440.,       0.,   -5040.,       0.,       0.,       0.,       0.,       0.,   -1440.,       0.,    5040.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.,       0.,       0.,       0.,       0.,       0.,  -27000.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   37800.,       0.,       0.,
   0.,       0.,       0.,       0.,       0.,      -0.,       0.,   -5040.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    5040.,       0.,       0.,       0.,       0.,       0.,     0.,       0.,       0.,       0.,       0.,  -37800.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   63000.,       0.,       0.,       0.,
-840.,       0.,    2100.,       0.,       0.,       0.,       0.,       0.,     840.,       0.,   -2100.,       0.,       0.,       0.,       0.,       0.,       0.,    2835.,       0.,       0.,       0.,     0.,       0.,   -4725.,       0.,  -14175.,       0.,       0.,       0.,       0.,       0.,      -0.,       0.,       0.,       0.,   23625.,       0.,       0.,       0.,       0.,
   0.,       0.,       0.,       0.,       0.,    2520.,       0.,   -8820.,       0.,       0.,       0.,       0.,       0.,    2520.,       0.,   -8820.,       0.,       0.,       0.,       0.,   23625.,     0.,       0.,       0.,       0.,       0.,       0.,       0.,  -33075.,       0., -118125.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,  165375.,       0.,
   0.,      -0.,       0.,   12600.,       0.,       0.,       0.,       0.,       0.,      -0.,       0.,   12600.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   47250.,       0.,     0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,  -78750.,       0.,       0.,       0., -198450.,       0.,       0.,       0.,       0.,       0.,  330750.,
   0.,       0.,       0.,       0.,       0.,       0.,    3600.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    3600.,       0.,       0.,       0.,    5400.,       0.,       0.,     0.,      -0.,       0.,      -0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,  -16200.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,
   0.,       0.,    3600.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    3600.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,  5400.,       0.,      -0.,       0.,      -0.,       0.,       0.,       0.,       0.,       0.,  -16200.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,
   0.,       0.,       0.,       0.,       0.,       0.,       0.,  -12600.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,  -12600.,       0.,       0.,       0.,       0.,   47250.,     0.,       0.,       0.,       0.,       0.,       0.,       0.,  -78750.,       0., -198450.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,  330750.,       0.,
   0.,   -2520.,       0.,    8820.,       0.,       0.,       0.,       0.,       0.,   -2520.,       0.,    8820.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   23625.,       0.,     0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,  -33075.,       0.,       0.,       0., -118125.,       0.,       0.,       0.,       0.,       0.,  165375.,
   0.,       0.,       0.,       0.,       0.,       0.,    6300.,       0.,       0.,       0.,       0.,       0.,      -0.,       0.,   -6300.,       0.,   -4725.,       0.,      -0.,       0.,       0.,     0.,   23625.,       0.,   14175.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,  -70875.,       0.,       0.,       0.,       0.,       0.,
   0.,       0.,       0.,    8400.,       0.,       0.,       0.,       0.,       0.,      -0.,       0.,   -8400.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,     0.,       0.,       0.,       0.,       0.,       0.,   37800.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,  -63000.,       0.,       0.,
   0.,       0.,       0.,       0.,       0.,       0.,       0.,    8400.,       0.,       0.,       0.,       0.,       0.,      -0.,       0.,   -8400.,       0.,       0.,       0.,       0.,       0.,     0.,       0.,       0.,       0.,       0.,   37800.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,  -63000.,       0.,       0.,       0.,
   0.,       0.,   -6300.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,    6300.,       0.,       0.,       0.,       0.,       0.,       0.,   -4725.,       0.,       0.,       0.,     0.,       0.,   14175.,       0.,   23625.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,  -70875.,       0.,       0.,       0.,       0.,
   0.,       0.,       0.,       0.,       0.,      -0.,       0.,   14700.,       0.,       0.,       0.,       0.,       0.,      -0.,       0.,   14700.,       0.,       0.,       0.,       0.,  -33075.,     0.,       0.,       0.,       0.,       0.,       0.,       0.,   55125.,       0.,  165375.,       0.,       0.,       0.,       0.,       0.,       0.,       0., -275625.,       0.,
   0.,       0.,       0.,  -14700.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,  -14700.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,  -33075.,       0.,     0.,       0.,       0.,       0.,       0.,       0.,       0.,       0.,   55125.,       0.,       0.,       0.,  165375.,       0.,       0.,       0.,       0.,       0., -275625.
};

static void N_Q_RT3_2D_Funct(double xi, double eta, double *values)
{
  int nBF = 40; // number of basis functions
  // monomials x-component and y-component
  double mon_x[40]={1,0,xi,0 ,eta,0,
  xi*xi,0, xi*eta, 0, eta*eta, 0, 
  xi*xi*xi, 0,   xi*xi*eta,    0, xi*eta*eta,    0, eta*eta*eta,   0,
  xi*xi*xi*xi,0, xi*xi*xi*eta, 0, xi*xi*eta*eta, 0, xi*eta*eta*eta,0,
  xi*xi*xi*xi*eta,0,    xi*xi*xi*eta*eta,    0,xi*xi*eta*eta*eta, 0,
  xi*xi*xi*xi*eta*eta,0,xi*xi*xi*eta*eta*eta,0,
  xi*xi*xi*xi*eta*eta*eta,0  };
  double mon_y[40]={0, 1, 0, xi, 0, eta, 
  0, xi*xi, 0, xi*eta, 0, eta*eta, 
  0, xi*xi*xi,     0, xi*xi*eta,     0, xi*eta*eta,     0,eta*eta*eta,
  0, xi*xi*xi*eta, 0, xi*xi*eta*eta, 0, xi*eta*eta*eta, 0, eta*eta*eta*eta,
  0, xi*xi*xi*eta*eta,     0, xi*xi*eta*eta*eta,    0, xi*eta*eta*eta*eta,
  0, xi*xi*xi*eta*eta*eta, 0, xi*xi*eta*eta*eta*eta,
  0, xi*xi*xi*eta*eta*eta*eta};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_RT3_2D_CM[i+j*nBF]*mon_x[j] / 768.;
      values[i+nBF] += N_Q_RT3_2D_CM[i+j*nBF]*mon_y[j] / 768.;
    }
  }
}

// values of the derivatives in xi direction
static void N_Q_RT3_2D_DeriveXi(double xi, double eta, double *values)
{
  int nBF = 40; // number of basis functions
  // monomials x-component and y-component
  double mon_x[40]={0,0,1,0 ,0,0,
  2*xi,0,   eta, 0, 0, 0, 
  3*xi*xi, 0,   2*xi*eta,    0, eta*eta,    0, 0,   0,
  4*xi*xi*xi,0, 3*xi*xi*eta, 0, 2*xi*eta*eta, 0, eta*eta*eta,0,
  4*xi*xi*xi*eta,0,    3*xi*xi*eta*eta,    0,2*xi*eta*eta*eta, 0,
  4*xi*xi*xi*eta*eta,0,3*xi*xi*eta*eta*eta,0,
  4*xi*xi*xi*eta*eta*eta,0  };
  double mon_y[40]={0, 0, 0, 1, 0, 0, 
  0, 2*xi, 0, eta, 0, 0, 
  0, 3*xi*xi,     0, 2*xi*eta,     0, eta*eta,     0,0,
  0, 3*xi*xi*eta, 0, 2*xi*eta*eta, 0, eta*eta*eta, 0, 0,
  0, 3*xi*xi*eta*eta,     0, 2*xi*eta*eta*eta,    0, eta*eta*eta*eta,
  0, 3*xi*xi*eta*eta*eta, 0, 2*xi*eta*eta*eta*eta,
  0, 3*xi*xi*eta*eta*eta*eta};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_RT3_2D_CM[i+j*nBF]*mon_x[j] / 768.;
      values[i+nBF] += N_Q_RT3_2D_CM[i+j*nBF]*mon_y[j] / 768.;
    }
  }
}

// values of the derivatives in eta direction
static void N_Q_RT3_2D_DeriveEta(double xi, double eta, double *values)
{
  int nBF = 40; // number of basis functions
  // monomials x-component and y-component
  double mon_x[40]={0,0,0,0 ,1,0,
  0,0, xi, 0, 2*eta, 0, 
  0, 0,   xi*xi,    0, 2*xi*eta,    0, 3*eta*eta,   0,
  0,0, xi*xi*xi, 0, 2*xi*xi*eta, 0, 3*xi*eta*eta,0,
  xi*xi*xi*xi,0,    2*xi*xi*xi*eta,    0,3*xi*xi*eta*eta, 0,
  2*xi*xi*xi*xi*eta,0,3*xi*xi*xi*eta*eta,0,
  3*xi*xi*xi*xi*eta*eta,0  };
  double mon_y[40]={0, 0, 0, 0, 0, 1, 
  0, 0, 0, xi, 0, 2*eta,
  0, 0,     0, xi*xi,     0, 2*xi*eta,     0,3*eta*eta,
  0, xi*xi*xi, 0, 2*xi*xi*eta, 0, 3*xi*eta*eta, 0, 4*eta*eta*eta,
  0, 2*xi*xi*xi*eta,     0, 3*xi*xi*eta*eta,    0, 4*xi*eta*eta*eta,
  0, 3*xi*xi*xi*eta*eta, 0, 4*xi*xi*eta*eta*eta,
  0, 4*xi*xi*xi*eta*eta*eta};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_RT3_2D_CM[i+j*nBF]*mon_x[j] / 768.;
      values[i+nBF] += N_Q_RT3_2D_CM[i+j*nBF]*mon_y[j] / 768.;
    }
  }
}

// values of derivatives in xi-xi direction
static void N_Q_RT3_2D_DeriveXiXi(double xi, double eta, double *values)
{
  int nBF = 40; // number of basis functions
  // monomials x-component and y-component
  double mon_x[40]={0,0,0,0 ,0,0,
  2,0,   0, 0, 0, 0, 
  6*xi, 0,   2*eta,    0, 0,    0, 0,   0,
  12*xi*xi,0, 6*xi*eta, 0, 2*eta*eta, 0, 0,0,
  12*xi*xi*eta,0,    6*xi*eta*eta,    0,2*eta*eta*eta, 0,
  12*xi*xi*eta*eta,0,6*xi*eta*eta*eta,0,
  12*xi*xi*eta*eta*eta,0  };
  double mon_y[40]={0, 0, 0, 0, 0, 0, 
  0, 2, 0, 0, 0, 0,
  0, 6*xi,     0, 2*eta,     0, 0,     0,0,
  0, 6*xi*eta, 0, 2*eta*eta, 0, 0, 0, 0,
  0, 6*xi*eta*eta,     0, 2*eta*eta*eta,    0, 0,
  0, 6*xi*eta*eta*eta, 0, 2*eta*eta*eta*eta,
  0, 6*xi*eta*eta*eta*eta};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_RT3_2D_CM[i+j*nBF]*mon_x[j] / 768.;
      values[i+nBF] += N_Q_RT3_2D_CM[i+j*nBF]*mon_y[j] / 768.;
    }
  }
}

// values of derivatives in eta-eta direction
static void N_Q_RT3_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  int nBF = 40; // number of basis functions
  // monomials x-component and y-component
  double mon_x[40]={0,0,0,0 ,0,0,
  0,0, 0, 0, 2, 0, 
  0, 0,   0,    0, 2*xi,    0, 6*eta,   0,
  0,0, 0, 0, 2*xi*xi, 0, 6*xi*eta,0,
  0,0,    2*xi*xi*xi,    0,6*xi*xi*eta, 0,
  2*xi*xi*xi*xi, 0, 6*xi*xi*xi*eta, 0,
  6*xi*xi*xi*xi*eta,0  };
  double mon_y[40]={0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 2, 
  0, 0,     0, 0,     0, 2*xi,     0,6*eta,
  0, 0, 0, 2*xi*xi, 0, 6*xi*eta, 0, 12*eta*eta,
  0, 2*xi*xi*xi,     0, 6*xi*xi*eta,    0, 12*eta*eta*eta,
  0, 6*xi*xi*xi*eta, 0, 12*xi*xi*eta*eta,
  0, 12*xi*xi*xi*eta*eta};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_RT3_2D_CM[i+j*nBF]*mon_x[j] / 768.;
      values[i+nBF] += N_Q_RT3_2D_CM[i+j*nBF]*mon_y[j] / 768.;
    }
  }
}

// values of derivatives in xi-eta direction
static void N_Q_RT3_2D_DeriveXiEta(double xi, double eta, double *values)
{
  int nBF = 40; // number of basis functions
  // monomials x-component and y-component
  double mon_x[40]={0,0,0,0 ,0,0,
  0,0,   1, 0, 0, 0, 
  0, 0,   2*xi,    0, 2*eta,    0, 0,   0,
  0,0, 3*xi*xi, 0, 4*xi*eta, 0, 3*eta*eta,0,
  4*xi*xi*xi,0,    6*xi*xi*eta,    0,6*xi*eta*eta, 0,
  8*xi*xi*xi*eta,0,9*xi*xi*eta*eta,0,
  12*xi*xi*xi*eta*eta,0  };
  double mon_y[40]={0, 0, 0, 0, 0, 0, 
  0, 0, 0, 1, 0, 0, 
  0, 0,     0, 2*xi,     0, 2*eta,     0,0,
  0, 3*xi*xi, 0, 4*xi*eta, 0, 3*eta*eta, 0, 0,
  0, 6*xi*xi*eta,     0, 6*xi*eta*eta,    0, 4*eta*eta*eta,
  0, 9*xi*xi*eta*eta, 0, 8*xi*eta*eta*eta,
  0, 12*xi*xi*eta*eta*eta};
  
  memset(values, 0.0, 2*nBF*sizeof(double)); // 2 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i    ] += N_Q_RT3_2D_CM[i+j*nBF]*mon_x[j] / 768.;
      values[i+nBF] += N_Q_RT3_2D_CM[i+j*nBF]*mon_y[j] / 768.;
    }
  }
}

// all dofs N(v) on the edges are integrals of (v.n p) where p is a polynomial.
// For the first and third dof on each edge p is even (p(x) = p(-x)), for the
// second and fourth it is odd (p(x)=-p(-x)). Since the direction of the 
// integration is reversed in the neighbor cell, only the first and third dof on
// each edge have to be changed with TBaseFunct2D::ChangeBF

static int N_Q_RT3_2D_ChangeJ0[2] = { 0, 2 };
static int N_Q_RT3_2D_ChangeJ1[2] = { 4, 6 };
static int N_Q_RT3_2D_ChangeJ2[2] = { 8, 10 };
static int N_Q_RT3_2D_ChangeJ3[2] = { 12, 14 };

static int *N_Q_RT3_2D_Change1[4] = { N_Q_RT3_2D_ChangeJ0, N_Q_RT3_2D_ChangeJ1,
                                     N_Q_RT3_2D_ChangeJ2, N_Q_RT3_2D_ChangeJ3 };
static int **N_Q_RT3_2D_Change[] = { N_Q_RT3_2D_Change1 };
