// ***********************************************************************
// Raviart-Thomas element of zero-th order on hexahedra, 3D
// ***********************************************************************

static void N_H_RT0_3D_Funct(double xi, double eta, double zeta,
                             double *values)
{
  // first component
  values[0] = 0;
  values[1] = 0;
  values[2] = 0.125*(xi+1);
  values[3] = 0;
  values[4] = 0.125*(xi-1);
  values[5] = 0;
  
  // second component
  values[6] = 0;
  values[7] = 0.125*(eta-1);
  values[8] = 0;
  values[9] = 0.125*(eta+1);
  values[10] = 0;
  values[11] = 0;
  
  // third component
  values[12] = 0.125*(zeta-1);
  values[13] = 0;
  values[14] = 0;
  values[15] = 0;
  values[16] = 0;
  values[17] = 0.125*(zeta+1);
}

static void N_H_RT0_3D_DeriveXi(double, double, double, double *values)
{
  // first component
  values[0] = 0;
  values[1] = 0;
  values[2] = 0.125;
  values[3] = 0;
  values[4] = 0.125;
  values[5] = 0;
  
  // second component
  values[6] = 0;
  values[7] = 0;
  values[8] = 0;
  values[9] = 0;
  values[10] = 0;
  values[11] = 0;
  
  // third component
  values[12] = 0;
  values[13] = 0;
  values[14] = 0;
  values[15] = 0;
  values[16] = 0;
  values[17] = 0;
}

static void N_H_RT0_3D_DeriveEta(double, double, double, double *values)
{
  // first component
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
  values[5] = 0;
  
  // second component
  values[6] = 0;
  values[7] = 0.125;
  values[8] = 0;
  values[9] = 0.125;
  values[10] = 0;
  values[11] = 0;
  
  // third component
  values[12] = 0;
  values[13] = 0;
  values[14] = 0;
  values[15] = 0;
  values[16] = 0;
  values[17] = 0;
}

static void N_H_RT0_3D_DeriveZeta(double, double, double, double *values)
{
  // first component
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
  values[5] = 0;
  
  // second component
  values[6] = 0;
  values[7] = 0;
  values[8] = 0;
  values[9] = 0;
  values[10] = 0;
  values[11] = 0;
  
  // third component
  values[12] = 0.125;
  values[13] = 0;
  values[14] = 0;
  values[15] = 0;
  values[16] = 0;
  values[17] = 0.125;
}

static void N_H_RT0_3D_DeriveXiXi(double, double, double, double *values)
{
  for(int i = 0; i < 18; i++)
    values[i] = 0;
}

static void N_H_RT0_3D_DeriveXiEta(double, double, double, double *values)
{
  for(int i = 0; i < 18; i++)
    values[i] = 0;
}

static void N_H_RT0_3D_DeriveXiZeta(double, double, double, double *values)
{
  for(int i = 0; i < 18; i++)
    values[i] = 0;
}

static void N_H_RT0_3D_DeriveEtaEta(double, double, double, double *values)
{
  for(int i = 0; i < 18; i++)
    values[i] = 0;
}

static void N_H_RT0_3D_DeriveEtaZeta(double, double, double, double *values)
{
  for(int i = 0; i < 18; i++)
    values[i] = 0;
}

static void N_H_RT0_3D_DeriveZetaZeta(double, double, double, double *values)
{
  for(int i = 0; i < 18; i++)
    values[i] = 0;
}

static int N_H_RT0_3D_ChangeJ0[1] = { 0 };
static int N_H_RT0_3D_ChangeJ1[1] = { 1 };
static int N_H_RT0_3D_ChangeJ2[1] = { 2 };
static int N_H_RT0_3D_ChangeJ3[1] = { 3 };
static int N_H_RT0_3D_ChangeJ4[1] = { 4 };
static int N_H_RT0_3D_ChangeJ5[1] = { 5 };

static int *N_H_RT0_3D_Change1[6] = { N_H_RT0_3D_ChangeJ0, N_H_RT0_3D_ChangeJ1,
                                      N_H_RT0_3D_ChangeJ2, N_H_RT0_3D_ChangeJ3,
                                      N_H_RT0_3D_ChangeJ4, N_H_RT0_3D_ChangeJ5};
static int **N_H_RT0_3D_Change2 = N_H_RT0_3D_Change1;
static int **N_H_RT0_3D_Change[] = { N_H_RT0_3D_Change1, N_H_RT0_3D_Change2 };
