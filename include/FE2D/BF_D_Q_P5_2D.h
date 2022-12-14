// ***********************************************************************
// P5 element, discontinous, 2D, quadrilateral
// ***********************************************************************

// base function values
static void D_Q_P5_2D_Funct(double xi, double eta, double *values)
{
  double t1 = xi*xi;
  double t4 = xi*eta;
  double t5 = eta*eta;
  double t8 = t1*xi;
  double t13 = -1.0+3.0*t1;
  double t14 = t13*eta;
  double t17 = -1.0+3.0*t5;
  double t20 = t5*eta;
  double t24 = t1*t1;
  double t30 = xi*(5.0*t1-3.0);
  double t36 = 5.0*t5-3.0;
  double t39 = t5*t5;

  values[0] = 1.0;
  values[1] = xi;
  values[2] = eta;
  values[3] = -1.0/2.0+3.0/2.0*t1;
  values[4] = t4;
  values[5] = -1.0/2.0+3.0/2.0*t5;
  values[6] = 5.0/2.0*t8-3.0/2.0*xi;
  values[7] = t14/2.0;
  values[8] = xi*t17/2.0;
  values[9] = 5.0/2.0*t20-3.0/2.0*eta;
  values[10] = 3.0/8.0+35.0/8.0*t24-15.0/4.0*t1;
  values[11] = t30*eta/2.0;
  values[12] = t13*t17/4.0;
  values[13] = t4*t36/2.0;
  values[14] = 3.0/8.0+35.0/8.0*t39-15.0/4.0*t5;
  values[15] = 63.0/8.0*t24*xi-35.0/4.0*t8+15.0/8.0*xi;
  values[16] = (3.0+35.0*t24-30.0*t1)*eta/8.0;
  values[17] = t30*t17/4.0;
  values[18] = t14*t36/4.0;
  values[19] = xi*(3.0+35.0*t39-30.0*t5)/8.0;
  values[20] = 63.0/8.0*t39*eta-35.0/4.0*t20+15.0/8.0*eta;
}

// values of the derivatives in xi direction
static void D_Q_P5_2D_DeriveXi(double xi, double eta, double *values)
{
  double t2 = xi*xi;
  double t5 = xi*eta;
  double t7 = eta*eta;
  double t23 = 5.0*t7-3.0;
  double t26 = t2*t2;
  double t42 = t7*t7;

  values[0] = 0.0;
  values[1] = 1.0;
  values[2] = 0.0;
  values[3] = 3.0*xi;
  values[4] = eta;
  values[5] = 0.0;
  values[6] = 15.0/2.0*t2-3.0/2.0;
  values[7] = 3.0*t5;
  values[8] = -1.0/2.0+3.0/2.0*t7;
  values[9] = 0.0;
  values[10] = 35.0/2.0*t2*xi-15.0/2.0*xi;
  values[11] = 15.0/2.0*t2*eta-3.0/2.0*eta;
  values[12] = 3.0/2.0*xi*(-1.0+3.0*t7);
  values[13] = eta*t23/2.0;
  values[14] = 0.0;
  values[15] = 315.0/8.0*t26-105.0/4.0*t2+15.0/8.0;
  values[16] = 5.0/2.0*xi*(7.0*t2-3.0)*eta;
  values[17] = -15.0/4.0*t2+45.0/4.0*t2*t7+3.0/4.0-9.0/4.0*t7;
  values[18] = 3.0/2.0*t5*t23;
  values[19] = 3.0/8.0+35.0/8.0*t42-15.0/4.0*t7;
  values[20] = 0.0;
}

// values of the derivatives in eta direction
static void D_Q_P5_2D_DeriveEta(double xi, double eta, double *values)
{
  double t2 = xi*xi;
  double t5 = xi*eta;
  double t7 = eta*eta;
  double t12 = xi*(5.0*t2-3.0);
  double t26 = t2*t2;
  double t41 = t7*t7;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 1.0;
  values[3] = 0.0;
  values[4] = xi;
  values[5] = 3.0*eta;
  values[6] = 0.0;
  values[7] = -1.0/2.0+3.0/2.0*t2;
  values[8] = 3.0*t5;
  values[9] = 15.0/2.0*t7-3.0/2.0;
  values[10] = 0.0;
  values[11] = t12/2.0;
  values[12] = 3.0/2.0*(-1.0+3.0*t2)*eta;
  values[13] = 15.0/2.0*xi*t7-3.0/2.0*xi;
  values[14] = 35.0/2.0*t7*eta-15.0/2.0*eta;
  values[15] = 0.0;
  values[16] = 3.0/8.0+35.0/8.0*t26-15.0/4.0*t2;
  values[17] = 3.0/2.0*t12*eta;
  values[18] = -15.0/4.0*t7+3.0/4.0+45.0/4.0*t2*t7-9.0/4.0*t2;
  values[19] = 5.0/2.0*t5*(7.0*t7-3.0);
  values[20] = 315.0/8.0*t41-105.0/4.0*t7+15.0/8.0;
}

// values of the derivatives in xi-xi direction
static void D_Q_P5_2D_DeriveXiXi(double xi, double eta, 
                                       double *values)
{
  double t3 = xi*xi;
  double t8 = eta*eta;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 3.0;
  values[4] = 0.0;
  values[5] = 0.0;
  values[6] = 15.0*xi;
  values[7] = 3.0*eta;
  values[8] = 0.0;
  values[9] = 0.0;
  values[10] = 105.0/2.0*t3-15.0/2.0;
  values[11] = 15.0*xi*eta;
  values[12] = -3.0/2.0+9.0/2.0*t8;
  values[13] = 0.0;
  values[14] = 0.0;
  values[15] = 315.0/2.0*t3*xi-105.0/2.0*xi;
  values[16] = 105.0/2.0*t3*eta-15.0/2.0*eta;
  values[17] = -15.0/2.0*xi+45.0/2.0*xi*t8;
  values[18] = 3.0/2.0*eta*(5.0*t8-3.0);
  values[19] = 0.0;
  values[20] = 0.0;
}

// values of the derivatives in xi-eta direction
static void D_Q_P5_2D_DeriveXiEta(double xi, double eta, 
                                       double *values)
{
  double t3 = xi*xi;
  double t8 = eta*eta;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 1.0;
  values[5] = 0.0;
  values[6] = 0.0;
  values[7] = 3.0*xi;
  values[8] = 3.0*eta;
  values[9] = 0.0;
  values[10] = 0.0;
  values[11] = 15.0/2.0*t3-3.0/2.0;
  values[12] = 9.0*xi*eta;
  values[13] = -3.0/2.0+15.0/2.0*t8;
  values[14] = 0.0;
  values[15] = 0.0;
  values[16] = 5.0/2.0*xi*(7.0*t3-3.0);
  values[17] = 45.0/2.0*t3*eta-9.0/2.0*eta;
  values[18] = 45.0/2.0*xi*t8-9.0/2.0*xi;
  values[19] = 35.0/2.0*t8*eta-15.0/2.0*eta;
  values[20] = 0.0;
}

// values of the derivatives in eta-eta direction
static void D_Q_P5_2D_DeriveEtaEta(double xi, double eta, 
                                       double *values)
{
  double t3 = xi*xi;
  double t8 = eta*eta;

  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;
  values[3] = 0.0;
  values[4] = 0.0;
  values[5] = 3.0;
  values[6] = 0.0;
  values[7] = 0.0;
  values[8] = 3.0*xi;
  values[9] = 15.0*eta;
  values[10] = 0.0;
  values[11] = 0.0;
  values[12] = 9.0/2.0*t3-3.0/2.0;
  values[13] = 15.0*xi*eta;
  values[14] = 105.0/2.0*t8-15.0/2.0;
  values[15] = 0.0;
  values[16] = 0.0;
  values[17] = 3.0/2.0*xi*(5.0*t3-3.0);
  values[18] = -15.0/2.0*eta+45.0/2.0*t3*eta;
  values[19] = 105.0/2.0*xi*t8-15.0/2.0*xi;
  values[20] = 315.0/2.0*t8*eta-105.0/2.0*eta;
}
