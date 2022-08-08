// ***********************************************************************
// Q1 element, discontinuous, 3D
// ***********************************************************************

static void D_H_Q1_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{
  values[0] = 1;
  values[1] = xi;
  values[2] = eta;
  values[3] = zeta;
  values[4] = xi*eta;
  values[5] = xi*zeta;
  values[6] = eta*zeta;
  values[7] = xi*eta*zeta;
}

static void D_H_Q1_3D_DeriveXi(double, double eta, double zeta, double *values)
{
  values[0] = 0;
  values[1] = 1;
  values[2] = 0;
  values[3] = 0;
  values[4] = eta;
  values[5] = zeta;
  values[6] = 0;
  values[7] = eta*zeta;
}

static void D_H_Q1_3D_DeriveEta(double xi, double, double zeta, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 1;
  values[3] = 0;
  values[4] = xi;
  values[5] = 0;
  values[6] = zeta;
  values[7] = xi*zeta;
}

static void D_H_Q1_3D_DeriveZeta(double xi, double eta, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 1;
  values[4] = 0;
  values[5] = xi;
  values[6] = eta;
  values[7] = xi*eta;
}

static void D_H_Q1_3D_DeriveXiXi(double, double, double, double *values)
{
  memset(values, 0.0, 8*sizeof(double));
}

static void D_H_Q1_3D_DeriveXiEta(double, double, double zeta, double *values)
{
  memset(values, 0.0, 8*sizeof(double));
  values[4] = 1;
  values[7] = zeta;
}

static void D_H_Q1_3D_DeriveXiZeta(double, double eta, double, double *values)
{
  memset(values, 0.0, 8*sizeof(double));
  values[5] = 1;
  values[7] = eta;
}

static void D_H_Q1_3D_DeriveEtaEta(double, double, double, double *values)
{
  memset(values, 0.0, 8*sizeof(double));
}

static void D_H_Q1_3D_DeriveEtaZeta(double xi, double, double, double *values)
{
  memset(values, 0.0, 8*sizeof(double));
  values[6] = 1;
  values[7] = xi;
}

static void D_H_Q1_3D_DeriveZetaZeta(double, double, double, double *values)
{
  memset(values, 0.0, 8*sizeof(double));
}
