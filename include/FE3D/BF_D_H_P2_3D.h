// ***********************************************************************
// P2 element, discontinuous, 3D
// ***********************************************************************

static void D_H_P2_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{
  values[0] = 1;
  values[1] = xi;
  values[2] = eta;
  values[3] = zeta;
  values[4] = 3*xi*xi-1;
  values[5] = xi*eta;
  values[6] = xi*zeta;
  values[7] = 3*eta*eta-1;
  values[8] = eta*zeta;
  values[9] = 3*zeta*zeta-1;
}

static void D_H_P2_3D_DeriveXi(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 1;
  values[2] = 0;
  values[3] = 0;
  values[4] = 6*xi;
  values[5] = eta;
  values[6] = zeta;
  values[7] = 0;
  values[8] = 0;
  values[9] = 0;
}

static void D_H_P2_3D_DeriveEta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 1;
  values[3] = 0;
  values[4] = 0;
  values[5] = xi;
  values[6] = 0;
  values[7] = 6*eta;
  values[8] = zeta;
  values[9] = 0;
}

static void D_H_P2_3D_DeriveZeta(double xi, double eta, double zeta,
                             double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 1;
  values[4] = 0;
  values[5] = 0;
  values[6] = xi;
  values[7] = 0;
  values[8] = eta;
  values[9] = 6*zeta;
}

static void D_H_P2_3D_DeriveXiXi(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 6;
  values[5] = 0;
  values[6] = 0;
  values[7] = 0;
  values[8] = 0;
  values[9] = 0;
}

static void D_H_P2_3D_DeriveXiEta(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
  values[5] = 1;
  values[6] = 0;
  values[7] = 0;
  values[8] = 0;
  values[9] = 0;
}

static void D_H_P2_3D_DeriveXiZeta(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
  values[5] = 0;
  values[6] = 1;
  values[7] = 0;
  values[8] = 0;
  values[9] = 0;
}

static void D_H_P2_3D_DeriveEtaEta(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
  values[5] = 0;
  values[6] = 0;
  values[7] = 6;
  values[8] = 0;
  values[9] = 0;
}

static void D_H_P2_3D_DeriveEtaZeta(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
  values[5] = 0;
  values[6] = 0;
  values[7] = 0;
  values[8] = 1;
  values[9] = 0;
}

static void D_H_P2_3D_DeriveZetaZeta(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
  values[5] = 0;
  values[6] = 0;
  values[7] = 0;
  values[8] = 0;
  values[9] = 6;
}
