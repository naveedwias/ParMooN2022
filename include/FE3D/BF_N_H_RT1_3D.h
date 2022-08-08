// ***********************************************************************
// Raviart-Thomas element of first order on hexahedra, 3D
// ***********************************************************************

static double N_H_RT1_3D_CM[1296] = {
//using tschebyscheff points (see NF_N_H_RT1_3D.h)
  0,0,0,0,0,0,0,0,-0.015625,-0.015625,-0.015625,-0.015625,0,0,-0,0,0.015625,0.015625,0.015625,0.015625,0,0,0,0,0.1875,-0,-0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0.015625,0.015625,0.015625,0.015625,0,0,0,0,-0.015625,-0.015625,-0.015625,-0.015625,0,0,0,0,0,0,0,0,0,0,0,0,0.1875,-0,-0,0,0,0,0,0,
  0.015625,0.015625,0.015625,0.015625,0,0,0,0,0,0,0,0,0,0,-0,0,0,0,0,0,-0.015625,-0.015625,-0.015625,-0.015625,0,0,0,0,0,0,0,0,0.1875,-0,0,0,
  0,0,0,0,0,0,0,0,0.03125,0.03125,0.03125,0.03125,0,0,0,0,0.03125,0.03125,0.03125,0.03125,0,0,0,0,-0,0,0,-0,0,0,0,0,0,0,0,0,
  0,0,0,0,-0.022097087,-0.022097087,0.022097087,0.022097087,0,0,0,0,-0.022097087,-0.022097087,0.022097087,0.022097087,0,0,0,0,0,0,0,0,0,0,0,0,-0,0.5625,-0,-0,0,0,0,0,
  -0.022097087,0.022097087,-0.022097087,0.022097087,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.022097087,0.022097087,-0.022097087,-0.022097087,0,0,0,0,0,0,0,0,-0,0.5625,-0,-0,
  0,0,0,0,0,0,0,0,0.022097087,0.022097087,-0.022097087,-0.022097087,0,0,0,0,-0.022097087,0.022097087,-0.022097087,0.022097087,0,0,0,0,0,0.5625,0,-0,0,0,0,0,0,0,0,0,
  0,0,0,0,0.03125,0.03125,0.03125,0.03125,0,0,0,0,0.03125,0.03125,0.03125,0.03125,0,0,0,0,0,0,0,0,0,0,0,0,-0,0,-0,-0,0,0,0,0,
  -0.022097087,-0.022097087,0.022097087,0.022097087,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.022097087,-0.022097087,0.022097087,-0.022097087,0,0,0,0,0,0,0,0,0,-0,0.5625,-0,
  0,0,0,0,0,0,0,0,0.022097087,-0.022097087,0.022097087,-0.022097087,0,0,0,0,-0.022097087,-0.022097087,0.022097087,0.022097087,0,0,0,0,-0,0,0.5625,-0,0,0,0,0,0,0,0,0,
  0,0,0,0,-0.022097087,0.022097087,-0.022097087,0.022097087,0,0,0,0,0.022097087,-0.022097087,0.022097087,-0.022097087,0,0,0,0,0,0,0,0,0,0,0,0,-0,0,0.5625,-0,0,0,0,0,
  0.03125,0.03125,0.03125,0.03125,0,0,0,0,0,0,0,0,0,0,-0,0,0,0,0,0,0.03125,0.03125,0.03125,0.03125,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0.046875,0.046875,0.046875,0.046875,0,0,-0,0,-0.046875,-0.046875,-0.046875,-0.046875,0,0,0,0,-0.1875,-0,-0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,-0.046875,-0.046875,-0.046875,-0.046875,0,0,0,0,0.046875,0.046875,0.046875,0.046875,0,0,0,0,0,0,0,0,0,0,0,0,-0.1875,0,0,-0,0,0,0,0,
  -0.046875,-0.046875,-0.046875,-0.046875,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.046875,0.046875,0.046875,0.046875,0,0,0,0,0,0,0,0,-0.1875,0,-0,-0,
  0,0,0,0,0,0,0,0,-0.044194174,-0.044194174,0.044194174,0.044194174,0,0,-0,0,-0.044194174,0.044194174,-0.044194174,0.044194174,0,0,0,0,0,-0,-0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,-0.044194174,-0.044194174,0.044194174,0.044194174,0,0,0,0,0.044194174,0.044194174,-0.044194174,-0.044194174,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,0,0,0,0,0,0,
  0.03125,-0.03125,-0.03125,0.03125,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.03125,0.03125,0.03125,-0.03125,0,0,0,0,0,0,0,0,0,-0,-0,1.6875,
  0,0,0,0,0,0,0,0,-0.044194174,0.044194174,-0.044194174,0.044194174,0,0,-0,0,-0.044194174,-0.044194174,0.044194174,0.044194174,0,0,0,0,0,-0,-0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0.03125,-0.03125,-0.03125,0.03125,0,0,0,0,0.03125,-0.03125,-0.03125,0.03125,0,0,0,0,0,0,0,0,0,0,0,0,-0,-0,0,1.6875,0,0,0,0,
  -0.044194174,0.044194174,-0.044194174,0.044194174,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.044194174,-0.044194174,0.044194174,0.044194174,0,0,0,0,0,0,0,0,-0,0,0,0,
  0,0,0,0,0,0,0,0,-0.03125,0.03125,0.03125,-0.03125,0,0,-0,0,0.03125,-0.03125,-0.03125,0.03125,0,0,0,0,0,-0,-0,1.6875,0,0,0,0,0,0,0,0,
  0,0,0,0,-0.044194174,0.044194174,-0.044194174,0.044194174,0,0,0,0,-0.044194174,0.044194174,-0.044194174,0.044194174,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,0,0,0,0,0,0,
  -0.044194174,-0.044194174,0.044194174,0.044194174,0,0,0,0,0,0,0,0,0,0,-0,0,0,0,0,0,-0.044194174,0.044194174,-0.044194174,0.044194174,0,0,0,0,0,0,0,0,-0,0,0,0,
  0,0,0,0,0,0,0,0,-0.066291261,-0.066291261,0.066291261,0.066291261,0,0,-0,0,0.066291261,-0.066291261,0.066291261,-0.066291261,0,0,0,0,0,-0.5625,-0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0.066291261,0.066291261,-0.066291261,-0.066291261,0,0,0,0,0.066291261,0.066291261,-0.066291261,-0.066291261,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.5625,0,0,0,0,0,0,
  0.066291261,-0.066291261,0.066291261,-0.066291261,0,0,0,0,0,0,0,0,0,0,-0,0,0,0,0,0,-0.066291261,-0.066291261,0.066291261,0.066291261,0,0,0,0,0,0,0,0,0,-0.5625,0,0,
  0,0,0,0,0,0,0,0,-0.066291261,0.066291261,-0.066291261,0.066291261,0,0,-0,0,0.066291261,0.066291261,-0.066291261,-0.066291261,0,0,0,0,-0,-0,-0.5625,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0.066291261,-0.066291261,0.066291261,-0.066291261,0,0,0,0,-0.066291261,0.066291261,-0.066291261,0.066291261,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,-0.5625,0,0,0,0,0,
  0.066291261,0.066291261,-0.066291261,-0.066291261,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.066291261,0.066291261,-0.066291261,0.066291261,0,0,0,0,0,0,0,0,-0,0,-0.5625,0,
  0,0,0,0,0,0,0,0,0.0625,-0.0625,-0.0625,0.0625,0,0,0,0,0.0625,-0.0625,-0.0625,0.0625,0,0,0,0,-0,0,0,-0,0,0,0,0,0,0,0,0,
  0,0,0,0,0.0625,-0.0625,-0.0625,0.0625,0,0,0,0,-0.0625,0.0625,0.0625,-0.0625,0,0,0,0,0,0,0,0,0,0,0,0,-0,0,-0,-0,0,0,0,0,
  0.0625,-0.0625,-0.0625,0.0625,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0625,-0.0625,-0.0625,0.0625,0,0,0,-0,0,0,0,-0,0,0,0,-0,
  0,0,0,0,0,0,0,0,0.09375,-0.09375,-0.09375,0.09375,0,0,0,0,-0.09375,0.09375,0.09375,-0.09375,0,0,0,0,0,0,0,-1.6875,0,0,0,-0,0,0,0,-0,
  0,0,0,0,-0.09375,0.09375,0.09375,-0.09375,0,0,0,0,-0.09375,0.09375,0.09375,-0.09375,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0,-1.6875,0,0,0,-0,
-0.09375,0.09375,0.09375,-0.09375,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.09375,-0.09375,-0.09375,0.09375,0,0,0,0,0,0,0,0,-0,0,0,-1.6875
};

static void N_H_RT1_3D_Funct(double xi, double eta, double zeta,
							 double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={1,0,0,xi,0,0,eta,0,0,zeta,0,0,
  xi*xi,0,0,xi*eta,0,0,xi*zeta,0,0,eta*zeta,0,0,
  xi*xi*eta,0,0,xi*xi*zeta,0,0,xi*eta*zeta,0,0,
  xi*xi*eta*zeta,0,0};
  double mon_y[]={0,1,0,0,xi,0,0,eta,0,0,zeta,0,
  0,eta*eta,0,0,xi*eta,0,0,xi*zeta,0,0,eta*zeta,0,
  0,xi*eta*eta,0,0,eta*eta*zeta,0,0,xi*eta*zeta,0,
  0,xi*eta*eta*zeta,0};
  double mon_z[]={0,0,1,0,0,xi,0,0,eta,0,0,zeta,
  0,0,zeta*zeta,0,0,xi*eta,0,0,xi*zeta,0,0,eta*zeta,
  0,0,xi*zeta*zeta,0,0,eta*zeta*zeta,0,0,xi*eta*zeta,
  0,0,xi*eta*zeta*zeta};
  
  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_RT1_3D_DeriveXi(double xi, double eta, double zeta,
                                double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,1,0,0,0,0,0,0,0,0,
  2*xi,0,0,eta,0,0,zeta,0,0,0,0,0,
  2*xi*eta,0,0,2*xi*zeta,0,0,eta*zeta,0,0,
  2*xi*eta*zeta,0,0};
  double mon_y[]={0,0,0,0,1,0,0,0,0,0,0,0,
  0,0,0,0,eta,0,0,zeta,0,0,0,0,
  0,eta*eta,0,0,0,0,0,eta*zeta,0,
  0,eta*eta*zeta,0};
  double mon_z[]={0,0,0,0,0,1,0,0,0,0,0,0,
  0,0,0,0,0,eta,0,0,zeta,0,0,0,
  0,0,zeta*zeta,0,0,0,0,0,eta*zeta,
  0,0,eta*zeta*zeta};
  
  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_RT1_3D_DeriveEta(double xi, double eta, double zeta,
                                 double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,1,0,0,0,0,0,
  0,0,0,xi,0,0,0,0,0,zeta,0,0,
  xi*xi,0,0,0,0,0,xi*zeta,0,0,
  xi*xi*zeta,0,0};
  double mon_y[]={0,0,0,0,0,0,0,1,0,0,0,0,
  0,2*eta,0,0,xi,0,0,0,0,0,zeta,0,
  0,xi*2*eta,0,0,2*eta*zeta,0,0,xi*zeta,0,
  0,xi*2*eta*zeta,0};
  double mon_z[]={0,0,0,0,0,0,0,0,1,0,0,0,
  0,0,0,0,0,xi,0,0,0,0,0,zeta,
  0,0,0,0,0,zeta*zeta,0,0,xi*zeta,
  0,0,xi*zeta*zeta};
  
  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_RT1_3D_DeriveZeta(double xi, double eta, double zeta,
                                  double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,1,0,0,
  0,0,0,0,0,0,xi,0,0,eta,0,0,
  0,0,0,xi*xi,0,0,xi*eta,0,0,
  xi*xi*eta,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,0,xi,0,0,eta,0,
  0,0,0,0,eta*eta,0,0,xi*eta,0,
  0,xi*eta*eta,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,1,
  0,0,2*zeta,0,0,0,0,0,xi,0,0,eta,
  0,0,xi*2*zeta,0,0,eta*2*zeta,0,0,xi*eta,
  0,0,xi*eta*2*zeta};
  
  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_RT1_3D_DeriveXiXi(double, double eta, double zeta,
                                  double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
  2,0,0,0,0,0,0,0,0,0,0,0,
  2*eta,0,0,2*zeta,0,0,0,0,0,
  2*eta*zeta,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,
  0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,
  0,0,0};

  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_RT1_3D_DeriveXiEta(double xi, double eta, double zeta,
                                   double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,1,0,0,0,0,0,0,0,0,
  2*xi,0,0,0,0,0,zeta,0,0,
  2*xi*zeta,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,1,0,0,0,0,0,0,0,
  0,2*eta,0,0,0,0,0,zeta,0,
  0,2*eta*zeta,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,1,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,zeta,
  0,0,zeta*zeta};

  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_RT1_3D_DeriveXiZeta(double xi, double eta, double zeta,
                                    double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,1,0,0,0,0,0,
  0,0,0,2*xi,0,0,eta,0,0,
  2*xi*eta,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,1,0,0,0,0,
  0,0,0,0,0,0,0,eta,0,
  0,eta*eta,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,1,0,0,0,
  0,0,2*zeta,0,0,0,0,0,eta,
  0,0,eta*2*zeta};

  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_RT1_3D_DeriveEtaEta(double xi, double, double zeta,
                                    double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,
  0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,2,0,0,0,0,0,0,0,0,0,0,
  0,xi*2,0,0,2*zeta,0,0,0,0,
  0,xi*2*zeta,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,
  0,0,0};

  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_RT1_3D_DeriveEtaZeta(double xi, double eta, double zeta,
                                     double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,1,0,0,
  0,0,0,0,0,0,xi,0,0,
  xi*xi,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,1,0,
  0,0,0,0,2*eta,0,0,xi,0,
  0,xi*2*eta,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,1,
  0,0,0,0,0,2*zeta,0,0,xi,
  0,0,xi*2*zeta};

  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static void N_H_RT1_3D_DeriveZetaZeta(double xi, double eta, double,
                                      double *values)
{
  int nBF = 36; // number of basis functions
  // monomials x-component, y-component and z-component
  double mon_x[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,
  0,0,0};
  double mon_y[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,
  0,0,0};
  double mon_z[]={0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,2,0,0,0,0,0,0,0,0,0,
  0,0,xi*2,0,0,eta*2,0,0,0,
  0,0,xi*eta*2};

  memset(values, 0.0, 3*nBF*sizeof(double)); // 3 is the space dimension
  for(int i=0; i<nBF; i++)
  {
    for(int j=0; j<nBF; j++)
    {
      values[i      ] += N_H_RT1_3D_CM[i+j*nBF]*mon_x[j];
      values[i+  nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_y[j];
      values[i+2*nBF] += N_H_RT1_3D_CM[i+j*nBF]*mon_z[j];
    }
  }
}

static int N_H_RT1_3D_ChangeJ0[4] = { 0, 1, 2, 3 };
static int N_H_RT1_3D_ChangeJ1[4] = { 4, 5, 6, 7 };
static int N_H_RT1_3D_ChangeJ2[4] = { 8, 9, 10, 11 };
static int N_H_RT1_3D_ChangeJ3[4] = { 12, 13, 14, 15 };
static int N_H_RT1_3D_ChangeJ4[4] = { 16, 17, 18, 19 };
static int N_H_RT1_3D_ChangeJ5[4] = { 20, 21, 22, 23 };

static int *N_H_RT1_3D_Change1[6] = { N_H_RT1_3D_ChangeJ0, N_H_RT1_3D_ChangeJ1,
                                      N_H_RT1_3D_ChangeJ2, N_H_RT1_3D_ChangeJ3,
                                      N_H_RT1_3D_ChangeJ4, N_H_RT1_3D_ChangeJ5};
static int **N_H_RT1_3D_Change2 = N_H_RT1_3D_Change1;
static int **N_H_RT1_3D_Change[] = { N_H_RT1_3D_Change1, N_H_RT1_3D_Change2 };
