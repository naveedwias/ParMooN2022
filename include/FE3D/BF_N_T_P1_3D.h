// ***********************************************************************
// P1 element, nonconforming, 3D
// ***********************************************************************

static void N_T_P1_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{
  values[0] = 1-3*zeta;
  values[1] = 1-3*eta;
  values[2] = 3*(xi+eta+zeta)-2;
  values[3] = 1-3*xi;
}

static void N_T_P1_3D_DeriveXi(double, double, double, double *values)
{
  values[0] =  0;
  values[1] =  0;
  values[2] =  3;
  values[3] = -3;
}

static void N_T_P1_3D_DeriveEta(double, double, double, double *values)
{
  values[0] =  0;
  values[1] = -3;
  values[2] =  3;
  values[3] =  0;
}

static void N_T_P1_3D_DeriveZeta(double, double, double, double *values)
{
  values[0] = -3;
  values[1] =  0;
  values[2] =  3;
  values[3] =  0;
}

static void N_T_P1_3D_DeriveXiXi(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void N_T_P1_3D_DeriveXiEta(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void N_T_P1_3D_DeriveXiZeta(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void N_T_P1_3D_DeriveEtaEta(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void N_T_P1_3D_DeriveEtaZeta(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

static void N_T_P1_3D_DeriveZetaZeta(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}
