// ***********************************************************************
// Q1 element, conforming, 3D
// ***********************************************************************

static void C_H_Q1_3D_Funct(double xi, double eta, double zeta,
                          double *values)
{
  double t1, t2, t4, t6, t9, t12, t14, t15, t17;

  t1 = 1.0/2.0-xi/2;
  t2 = 1.0/2.0-eta/2;
  t4 = 1.0/2.0-zeta/2;
  t6 = xi/2+1.0/2.0;
  t9 = eta/2+1.0/2.0;
  t12 = t6*t9;
  t14 = -1.0/2.0+xi/2;
  t15 = -1.0/2.0+eta/2;
  t17 = zeta/2+1.0/2.0;

  values[0] = t1*t2*t4;
  values[1] = t6*t2*t4;
  values[2] = t1*t9*t4;
  values[3] = t12*t4;
  values[4] = t14*t15*t17;
  values[5] = -t6*t15*t17;
  values[6] = -t14*t9*t17;
  values[7] = t12*t17;
}

static void C_H_Q1_3D_DeriveXi(double, double eta, double zeta, double *values)
{
  double t2, t3, t4, t5, t7, t8, t9;

  t2 = 1.0/2.0-zeta/2;
  t3 = (1.0/2.0-eta/2)*t2;
  t4 = eta/2+1.0/2.0;
  t5 = t4*t2;
  t7 = zeta/2+1.0/2.0;
  t8 = (-1.0/2.0+eta/2)*t7;
  t9 = t7*t4;

  values[0] = -t3/2;
  values[1] = t3/2;
  values[2] = -t5/2;
  values[3] = t5/2;
  values[4] = t8/2;
  values[5] = -t8/2;
  values[6] = -t9/2;
  values[7] = t9/2;
}

static void C_H_Q1_3D_DeriveEta(double xi, double, double zeta, double *values)
{
  double t2, t3, t4, t5, t7, t8, t9;

  t2 = 1.0/2.0-zeta/2;
  t3 = (1.0/2.0-xi/2)*t2;
  t4 = xi/2+1.0/2.0;
  t5 = t4*t2;
  t7 = zeta/2+1.0/2.0;
  t8 = (-1.0/2.0+xi/2)*t7;
  t9 = t7*t4;

  values[0] = -t3/2;
  values[1] = -t5/2;
  values[2] = t3/2;
  values[3] = t5/2;
  values[4] = t8/2;
  values[5] = -t9/2;
  values[6] = -t8/2;
  values[7] = t9/2;
}

static void C_H_Q1_3D_DeriveZeta(double xi, double eta, double, double *values)
{
  double t1, t2, t4, t6, t8, t9, t10;

  t1 = 1.0/2.0-xi/2;
  t2 = 1.0/2.0-eta/2;
  t4 = xi/2+1.0/2.0;
  t6 = eta/2+1.0/2.0;
  t8 = t4*t6;
  t9 = -1.0/2.0+xi/2;
  t10 = -1.0/2.0+eta/2;

  values[0] = -t1*t2/2;
  values[1] = -t4*t2/2;
  values[2] = -t1*t6/2;
  values[3] = -t8/2;
  values[4] = t9*t10/2;
  values[5] = -t10*t4/2;
  values[6] = -t6*t9/2;
  values[7] = t8/2;
}

static void C_H_Q1_3D_DeriveXiXi(double, double, double, double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 0.0;
  values[6] = 0.0;
  values[7] = 0.0;
}

static void C_H_Q1_3D_DeriveXiEta(double, double, double zeta, double *values)
{
  double t1, t2, t3, t4;

  t1 = 1.0/8.0-zeta/8;
  t2 = -1.0/8.0+zeta/8;
  t3 = zeta/8+1.0/8.0;
  t4 = -zeta/8-1.0/8.0;

  values[0] = t1;
  values[1] = t2;
  values[2] = t2;
  values[3] = t1;
  values[4] = t3;
  values[5] = t4;
  values[6] = t4;
  values[7] = t3;
}

static void C_H_Q1_3D_DeriveXiZeta(double, double eta, double, double *values)
{
  double t1, t2, t3, t4;

  t1 = 1.0/8.0-eta/8;
  t2 = -1.0/8.0+eta/8;
  t3 = eta/8+1.0/8.0;
  t4 = -eta/8-1.0/8.0;

  values[0] = t1;
  values[1] = t2;
  values[2] = t3;
  values[3] = t4;
  values[4] = t2;
  values[5] = t1;
  values[6] = t4;
  values[7] = t3;
}

static void C_H_Q1_3D_DeriveEtaEta(double, double, double, double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 0.0;
  values[6] = 0.0;
  values[7] = 0.0;
}

static void C_H_Q1_3D_DeriveEtaZeta(double xi, double, double, double *values)
{
  double t1, t2, t3, t4;

  t1 = 1.0/8.0-xi/8;
  t2 = xi/8+1.0/8.0;
  t3 = -1.0/8.0+xi/8;
  t4 = -xi/8-1.0/8.0;

  values[0] = t1;
  values[1] = t2;
  values[2] = t3;
  values[3] = t4;
  values[4] = t3;
  values[5] = t4;
  values[6] = t1;
  values[7] = t2;
}

static void C_H_Q1_3D_DeriveZetaZeta(double, double, double, double *values)
{
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 0.0;
  values[6] = 0.0;
  values[7] = 0.0;
}
