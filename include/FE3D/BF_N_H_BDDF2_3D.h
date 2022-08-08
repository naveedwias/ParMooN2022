// ***********************************************************************
// Brezzi-Douglas-Duran-Fortin element of second order on hexahedra, 3D
// ***********************************************************************

static double N_H_BDDF2_3D_CM[1521] = {
    -0,0,0,-0,-0,0,-0,0,0,0,0,-0,-0.0833333333,0.0833333333,0.0833333333,-0.0416666667,-0.0625,-0.0416666667,-0,0,0,-0,-0,-0,0.0833333333,-0.0833333333,-0.0833333333,0.0416666667,0.0625,0.0416666667,-0,0,0,-0,-0,-0,0.1875,0,0,
    -0,0,0,0,0,-0,0.0833333333,-0.0833333333,-0.0833333333,0.0416666667,0.0625,0.0416666667,-0,0,0,0,-0,-0,-0.0833333333,0.0833333333,0.0833333333,-0.0416666667,-0.0625,-0.0416666667,-0,0,0,-0,-0,0,-0,0,0,-0,0,0,0,0.1875,0,
    0.0833333333,-0.0833333333,-0.0833333333,0.0416666667,0.0625,0.0416666667,-0,0,0,-0,-0,0,-0,0,0,-0,0,0,-0,0,0,-0,-0,-0,-0,0,0,0,0,-0,-0.0833333333,0.0833333333,0.0833333333,-0.0416666667,-0.0625,-0.0416666667,0,0,0.1875,
    0.0277777778,-0.0555555556,-0,0.0277777778,0,-0,0.0277777778,-0,-0.0555555556,0,0,0.0277777778,-0,0,0,-0,0.125,-0,0.0277777778,-0,-0.0555555556,0,0,0.0277777778,-0,0,0,0,0.125,-0,0.0277777778,-0,-0.0555555556,-0,0,0.0277777778,0,0,-0,
    0,0,0,0,0,0,-0.0721687836,0.1443375673,0.1443375673,-0,-0.1443375673,-0.0721687836,-0,0,0,0,-0,0,0.0721687836,-0.1443375673,-0.1443375673,-0,0.1443375673,0.0721687836,-0,0,0,0,-0,0,0,0,0,0,0,0,0,-0,0,
    -0.0721687836,0.1443375673,0.1443375673,-0.0721687836,-0.1443375673,0,0.0833333333,-0.0833333333,-0.0833333333,0,0.0833333333,-0,0,0,0,0,0,0,0.0833333333,-0.0833333333,-0.0833333333,-0,0.0833333333,-0,0,0,0,0,0,0,0.0721687836,-0.1443375673,-0.1443375673,0,0.1443375673,0.0721687836,0,0,-0,
    0.0833333333,-0.0833333333,-0.0833333333,0,0.0833333333,0,0,0,0,0,0,0,0.0721687836,-0.1443375673,-0.1443375673,0,0.1443375673,0.0721687836,0,0,0,0,0,0,-0.0721687836,0.1443375673,0.1443375673,-0.0721687836,-0.1443375673,0,0.0833333333,-0.0833333333,-0.0833333333,-0,0.0833333333,0,0,0,0,
    0.0277777778,0,-0.0555555556,0,0,0.0277777778,-0,0,-0,0,0.125,0,0.0277777778,0,-0.0555555556,-0,0,0.0277777778,-0,0,0,-0,0.125,0,0.0277777778,-0.0555555556,0,0.0277777778,-0,0,0.0277777778,-0.0555555556,0,0.0277777778,0,0,0,-0,0,
    -0.0721687836,0.1443375673,0.1443375673,0,-0.1443375673,-0.0721687836,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0721687836,-0.1443375673,-0.1443375673,0.0721687836,0.1443375673,0,0,0,-0,
    -0,0,0,0,0,0,0,0,0,0,0,0,0.0721687836,-0.1443375673,-0.1443375673,0.0721687836,0.1443375673,0,0,0,0,0,0,0,-0.0721687836,0.1443375673,0.1443375673,0,-0.1443375673,-0.0721687836,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,-0.0721687836,0.1443375673,0.1443375673,-0.0721687836,-0.1443375673,0,0.0833333333,-0.0833333333,-0.0833333333,-0,0.0833333333,-0,0.0721687836,-0.1443375673,-0.1443375673,0.0721687836,0.1443375673,0,0.0833333333,-0.0833333333,-0.0833333333,-0,0.0833333333,0,0,0,0,0,0,0,-0,-0,0,
    0,0,0,-0,0.125,-0,0.0277777778,-0.0555555556,0,0.0277777778,0,0,0.0277777778,-0.0555555556,0,0.0277777778,0,0,0.0277777778,-0.0555555556,0,0.0277777778,0,0,0.0277777778,0,-0.0555555556,0,0,0.0277777778,0,0,0,-0,0.125,0,0,0,-0,
    -0,0,0,0,0,0,0,0,0,0,0,0,0.0833333333,-0.0833333333,-0.0833333333,0.0416666667,0.1875,0.0416666667,0,0,0,0,0,0,-0.0833333333,0.0833333333,0.0833333333,-0.0416666667,-0.1875,-0.0416666667,0,0,0,0,0,0,-0.1875,0,0,
    0,0,0,0,0,0,-0.0833333333,-0,0.1666666667,-0,0,-0.0833333333,-0,0,0,0,-0,0,0.0833333333,-0,-0.1666666667,0,0,0.0833333333,-0,0,0,0,-0,0,0,0,0,0,0,0,0,-0,0,
    -0.0833333333,0.1666666667,-0,-0.0833333333,0,0,-0,0,0,0,-0,-0,0,0,0,0,0,0,-0,0,0,0,-0,-0,0,0,0,0,0,0,0.0833333333,-0,-0.1666666667,0,0,0.0833333333,0,0,-0,
    -0,0,0,0,0,0,0,0,0,0,0,0,0.0721687836,-0.1443375673,-0.1443375673,0,0.1443375673,0.0721687836,0,0,0,0,0,0,0.0721687836,-0.1443375673,-0.1443375673,0.0721687836,0.1443375673,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0.0721687836,-0.1443375673,-0.1443375673,0,0.1443375673,0.0721687836,0,0,0,0,0,0,0.0721687836,-0.1443375673,-0.1443375673,-0,0.1443375673,0.0721687836,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,0,
    -0.1666666667,0.1666666667,0.1666666667,0,-0.1666666667,-0,0,-0,-0,0,0,-0,0,0,0,0,0,0,0,-0,-0,-0,0,-0,0,0,0,0,0,0,0.1666666667,-0.1666666667,-0.1666666667,-0,0.1666666667,0,0,0,0,
    -0,0,0,0,0,0,0,0,0,0,0,0,0.0721687836,-0.1443375673,-0.1443375673,0.0721687836,0.1443375673,0,0,0,0,0,0,0,0.0721687836,-0.1443375673,-0.1443375673,0,0.1443375673,0.0721687836,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,-0.1666666667,0.1666666667,0.1666666667,-0,-0.1666666667,-0,-0,0,0,0,-0,0,0.1666666667,-0.1666666667,-0.1666666667,-0,0.1666666667,-0,-0,0,0,0,-0,0,0,0,0,0,0,0,0,-0,0,
    0.0721687836,-0.1443375673,-0.1443375673,0.0721687836,0.1443375673,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0721687836,-0.1443375673,-0.1443375673,-0,0.1443375673,0.0721687836,0,0,-0,
    -0,0,0,0,0,0,0,0,0,0,0,0,0.0833333333,0,-0.1666666667,0,0,0.0833333333,0,0,0,0,0,0,-0.0833333333,0.1666666667,0,-0.0833333333,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,-0.0833333333,0.0833333333,0.0833333333,-0.0416666667,-0.1875,-0.0416666667,0,-0,-0,-0,0,-0,0.0833333333,-0.0833333333,-0.0833333333,0.0416666667,0.1875,0.0416666667,0,-0,-0,-0,0,0,0,0,0,0,0,0,0,-0.1875,0,
    -0.0833333333,-0,0.1666666667,0,0,-0.0833333333,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0833333333,-0.1666666667,-0,0.0833333333,0,0,0,0,0,
    -0,0,0,0,0,0,0,0,0,0,0,0,0.1666666667,-0.1666666667,-0.1666666667,0,0.1666666667,0,0,0,0,0,0,0,-0.1666666667,0.1666666667,0.1666666667,0,-0.1666666667,0,0,0,0,0,0,0,-0,0,0,
    0,0,0,0,0,0,0.0721687836,-0.1443375673,-0.1443375673,0.0721687836,0.1443375673,0,0,0,0,0,0,0,0.0721687836,-0.1443375673,-0.1443375673,0.0721687836,0.1443375673,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0.0721687836,-0.1443375673,-0.1443375673,-0,0.1443375673,0.0721687836,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0721687836,-0.1443375673,-0.1443375673,0.0721687836,0.1443375673,0,0,0,-0,
    -0,0,0,0,0,0,0,0,0,0,0,0,0.0833333333,-0.1666666667,-0,0.0833333333,0,-0,0,0,0,0,0,0,-0.0833333333,-0,0.1666666667,-0,0,-0.0833333333,0,0,0,0,0,0,0,0,0,
    -0,0,0,0,0,0,-0.0833333333,0.1666666667,0,-0.0833333333,0,0,0,0,0,0,0,0,0.0833333333,-0.1666666667,0,0.0833333333,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    -0.0833333333,0.0833333333,0.0833333333,-0.0416666667,-0.1875,-0.0416666667,0,-0,-0,0,0,0,0,0,0,0,0,0,0,-0,-0,-0,0,0,0,0,0,0,0,0,0.0833333333,-0.0833333333,-0.0833333333,0.0416666667,0.1875,0.0416666667,0,0,-0.1875,
    -0,0,0,0,0,0,-0.0277777778,0,0.0555555556,-0,0,-0.0277777778,-0,0,0,-0,-0,-0,-0.0277777778,0,0.0555555556,-0,-0,-0.0277777778,-0,0,0,0,-0,0,0,0,0,0,0,0,-0,-0,0,
    -0,0,0,0,0,0,0,0,0,0,0,0,-0.0277777778,0.0555555556,0,-0.0277777778,0,0,0,0,0,0,0,0,-0.0277777778,0,0.0555555556,0,0,-0.0277777778,0,0,0,0,0,0,0,0,0,
    -0.0277777778,0,0.0555555556,0,0,-0.0277777778,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0277777778,0.0555555556,0,-0.0277777778,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0.0833333333,-0.0833333333,-0.0833333333,-0,0.0833333333,-0,0,0,0,0,0,0,0.0833333333,-0.0833333333,-0.0833333333,0,0.0833333333,0,0,0,0,0,0,0,-0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0.0277777778,0,-0.0555555556,0,0,0.0277777778,0,0,0,0,0,0,0.0277777778,-0.0555555556,0,0.0277777778,0,0,0,0,0,0,0,0,0,0,0,
    0.0833333333,-0.0833333333,-0.0833333333,0,0.0833333333,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0833333333,-0.0833333333,-0.0833333333,-0,0.0833333333,0,0,0,0,
    0.0277777778,-0.0555555556,-0,0.0277777778,0,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0277777778,-0,-0.0555555556,-0,0,0.0277777778,0,0,-0,
    0,0,0,0,0,0,0.0833333333,-0.0833333333,-0.0833333333,0,0.0833333333,-0,0,0,0,0,0,0,0.0833333333,-0.0833333333,-0.0833333333,-0,0.0833333333,-0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0.0277777778,-0.0555555556,0,0.0277777778,0,0,0,0,0,0,0,0,0.0277777778,-0.0555555556,0,0.0277777778,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

static void N_H_BDDF2_3D_Funct(double xi, double eta, double zeta,
							 double *values)
{
  int nBF = 39; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={1,0,0,xi,0,0,eta,0,0,zeta,0,0,
  xi*xi,0,0,xi*eta,0,0,xi*zeta,0,0,eta*eta,0,0,eta*zeta,0,0,zeta*zeta,0,0,
  xi*xi*xi,-3*xi*zeta*zeta,0,2*xi*eta*zeta,3*xi*eta*eta,-xi*xi*eta,-xi*xi*xi,0,0};
  double mon_y[]={0,1,0,0,xi,0,0,eta,0,0,zeta,0,
  0,xi*xi,0,0,xi*eta,0,0,xi*zeta,0,0,eta*eta,0,0,eta*zeta,0,0,zeta*zeta,0,
  -3*xi*xi*eta,0,eta*eta*eta,-eta*eta*zeta,-eta*eta*eta,0,0,2*xi*eta*zeta,3*eta*zeta*zeta};
  double mon_z[]={0,0,1,0,0,xi,0,0,eta,0,0,zeta,
  0,0,xi*xi,0,0,xi*eta,0,0,xi*zeta,0,0,eta*eta,0,0,eta*zeta,0,0,zeta*zeta,
  0,zeta*zeta*zeta,-3*eta*eta*zeta,0,0,2*xi*eta*zeta,3*xi*xi*zeta,-xi*zeta*zeta,-zeta*zeta*zeta};
  
  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF2_3D_DeriveXi(double xi, double eta, double zeta,
                                double *values)
{
  int nBF = 39; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,1,0,0,0,0,0,0,0,0,
  2*xi,0,0,eta,0,0,zeta,0,0,0,0,0,0,0,0,0,0,0,
  3*xi*xi,-3*zeta*zeta,0,2*eta*zeta,3*eta*eta,-2*xi*eta,-3*xi*xi,0,0};
  double mon_y[]={0,0,0,0,1,0,0,0,0,0,0,0,
  0,2*xi,0,0,eta,0,0,zeta,0,0,0,0,0,0,0,0,0,0,
  -3*2*xi*eta,0,0,0,0,0,0,2*eta*zeta,0};
  double mon_z[]={0,0,0,0,0,1,0,0,0,0,0,0,
  0,0,2*xi,0,0,eta,0,0,zeta,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,2*eta*zeta,3*2*xi*zeta,-zeta*zeta,0};
  
  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF2_3D_DeriveEta(double xi, double eta, double zeta,
                                 double *values)
{
  int nBF = 39; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,1,0,0,0,0,0,
  0,0,0,xi,0,0,0,0,0,2*eta,0,0,zeta,0,0,0,0,0,
  0,0,0,2*xi*zeta,3*xi*2*eta,-xi*xi,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,1,0,0,0,0,
  0,0,0,0,xi,0,0,0,0,0,2*eta,0,0,zeta,0,0,0,0,
  -3*xi*xi,0,3*eta*eta,-2*eta*zeta,-3*eta*eta,0,0,2*xi*zeta,3*zeta*zeta};
  double mon_z[]={0,0,0,0,0,0,0,0,1,0,0,0,
  0,0,0,0,0,xi,0,0,0,0,0,2*eta,0,0,zeta,0,0,0,
  0,0,-3*2*eta*zeta,0,0,2*xi*zeta,0,0,0};
  
  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF2_3D_DeriveZeta(double xi, double eta, double zeta,
                                  double *values)
{
  int nBF = 39; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,1,0,0,
  0,0,0,0,0,0,xi,0,0,0,0,0,eta,0,0,2*zeta,0,0,
  0,-3*xi*2*zeta,0,2*xi*eta,0,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,0,xi,0,0,0,0,0,eta,0,0,2*zeta,0,
  0,0,0,-eta*eta,0,0,0,2*xi*eta,3*eta*2*zeta};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,1,
  0,0,0,0,0,0,0,0,xi,0,0,0,0,0,eta,0,0,2*zeta,
  0,3*zeta*zeta,-3*eta*eta,0,0,2*xi*eta,3*xi*xi,-xi*2*zeta,-3*zeta*zeta};
  
  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF2_3D_DeriveXiXi(double xi, double eta, double zeta,
                                  double *values)
{
  int nBF = 39; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
  2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  3*2*xi,0,0,0,0,-2*eta,-3*2*xi,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  -3*2*eta,0,0,0,0,0,0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,3*2*zeta,0,0};

  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF2_3D_DeriveXiEta(double xi, double eta, double zeta,
                                   double *values)
{
  int nBF = 39; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,2*zeta,3*2*eta,-2*xi,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
  -3*2*xi,0,0,0,0,0,0,2*zeta,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,2*zeta,0,0,0};

  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF2_3D_DeriveXiZeta(double xi, double eta, double zeta,
                                    double *values)
{
  int nBF = 39; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
  0,-3*2*zeta,0,2*eta,0,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,2*eta,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,2*eta,3*2*xi,-2*zeta,0};

  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF2_3D_DeriveEtaEta(double xi, double eta, double zeta,
                                    double *values)
{
  int nBF = 39; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,
  0,0,0,0,3*xi*2,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,
  0,0,3*2*eta,-2*zeta,-3*2*eta,0,0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,
  0,0,-3*2*zeta,0,0,0,0,0,0};

  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF2_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                                     double *values)
{
  int nBF = 39; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,
  0,0,0,2*xi,0,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,
  0,0,0,-2*eta,0,0,0,2*xi,3*2*zeta};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
  0,0,-3*2*eta,0,0,2*xi,0,0,0};

  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_BDDF2_3D_DeriveZetaZeta(double xi, double eta, double zeta,
                                      double *values)
{
  int nBF = 39; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,
  0,-3*xi*2,0,0,0,0,0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,
  0,0,0,0,0,0,0,0,3*eta*2};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,
  0,3*2*zeta,0,0,0,0,0,-xi*2,-3*2*zeta};

  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_BDDF2_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static int N_H_BDDF2_3D_ChangeJ0[6] = {  0, 1, 2, 3, 4, 5 };
static int N_H_BDDF2_3D_ChangeJ1[6] = {  6, 7, 8, 9,10,11 };
static int N_H_BDDF2_3D_ChangeJ2[6] = { 12,13,14,15,16,17 };
static int N_H_BDDF2_3D_ChangeJ3[6] = { 18,19,20,21,22,23 };
static int N_H_BDDF2_3D_ChangeJ4[6] = { 24,25,26,27,28,29 };
static int N_H_BDDF2_3D_ChangeJ5[6] = { 30,31,32,33,34,35 };

static int *N_H_BDDF2_3D_Change1[6] = { N_H_BDDF2_3D_ChangeJ0,
                                        N_H_BDDF2_3D_ChangeJ1,
                                        N_H_BDDF2_3D_ChangeJ2,
                                        N_H_BDDF2_3D_ChangeJ3,
                                        N_H_BDDF2_3D_ChangeJ4,
                                        N_H_BDDF2_3D_ChangeJ5 };
static int **N_H_BDDF2_3D_Change2 = N_H_BDDF2_3D_Change1;
static int **N_H_BDDF2_3D_Change[] = { N_H_BDDF2_3D_Change1, N_H_BDDF2_3D_Change2 };