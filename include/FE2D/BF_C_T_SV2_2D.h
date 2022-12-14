// ***********************************************************************
// SV2 element, conforming, 2D
// ***********************************************************************

// base function values
static void C_T_SV2_2D_Funct(double xi, double eta, double *values)
{
  // first triangle
  if( (xi-eta)>0 && (xi+2*eta-1)<0 )
  {
    values[0] =  2*xi*xi +8*xi*eta +8*eta*eta -3*xi -6*eta +1;
    values[1] = -4*xi*xi -4*xi*eta +8*eta*eta +4*xi -4*eta;
    values[2] =  2*xi*xi -4*xi*eta +2*eta*eta -  xi +  eta;
    values[3] = 0;
    values[4] = 0;
    values[5] = 0;
    values[6] =         -12*xi*eta-24*eta*eta       +12*eta;
    values[7] =          12*xi*eta-12*eta*eta;
    values[8] = 0;
    values[9] =                   +18*eta*eta       -3*eta;
  }
  else
    // second triangle
    if( (xi+2*eta-1)>=0 && (-2*xi-eta+1)<0 )
    {
      values[0] = 0;
      values[1] = 0;
      values[2] =  8*xi*xi  +8*xi*eta  +2*eta*eta -10*xi  -5*eta +3;
      values[3] =  8*xi*xi +20*xi*eta  +8*eta*eta -12*xi -12*eta +4;
      values[4] =  2*xi*xi  +8*xi*eta  +8*eta*eta  -5*xi -10*eta +3;
      values[5] = 0;
      values[6] = 0;
      values[7] =-24*xi*xi -36*xi*eta -12*eta*eta +36*xi +24*eta -12;
      values[8] =-12*xi*xi -36*xi*eta -24*eta*eta +24*xi +36*eta -12;
      values[9] = 18*xi*xi +36*xi*eta +18*eta*eta -33*xi -33*eta +15;
    }
    else
    {
      values[0] =   8*xi*xi  +8*xi*eta +2*eta*eta -6*xi -3*eta +1;
      values[1] = 0;
      values[2] = 0;
      values[3] = 0;
      values[4] =   2*xi*xi  -4*xi*eta +2*eta*eta +  xi -  eta;
      values[5] =   8*xi*xi  -4*xi*eta -4*eta*eta -4*xi +4*eta;
      values[6] = -24*xi*xi -12*xi*eta           +12*xi;
      values[7] = 0;
      values[8] = -12*xi*xi +12*xi*eta;
      values[9] = 18*xi*xi                        -3*xi;                   
    }
}

// values of the derivatives in xi direction
static void C_T_SV2_2D_DeriveXi(double xi, double eta, double *values)
{
  // first triangle
  if( (xi-eta)>0 && (xi+2*eta-1)<0 )
  {
    values[0] =  4*xi +8*eta -3;
    values[1] = -8*xi -4*eta +4;
    values[2] =  4*xi -4*eta -1;
    values[3] = 0;
    values[4] = 0;
    values[5] = 0;
    values[6] =      -12*eta;
    values[7] =       12*eta;
    values[8] = 0;
    values[9] = 0;
  }
  else
    // second triangle
    if( (xi+2*eta-1)>=0 && (-2*xi-eta+1)<0 )
    {
      values[0] = 0;
      values[1] = 0;
      values[2] = 16*xi  +8*eta -10;
      values[3] = 16*xi +20*eta -12;
      values[4] =  4*xi  +8*eta  -5;
      values[5] = 0;
      values[6] = 0;
      values[7] =-48*xi -36*eta +36;
      values[8] =-24*xi -36*eta +24;
      values[9] = 36*xi +36*eta -33;
    }
    else
    {
      values[0] =  16*xi  +8*eta -6;
      values[1] = 0;
      values[2] = 0;
      values[3] = 0;
      values[4] =   4*xi  -4*eta +1;
      values[5] =  16*xi  -4*eta -4;
      values[6] = -48*xi -12*eta +12;
      values[7] = 0;
      values[8] = -24*xi +12*eta;
      values[9] =  36*xi         -3;
    }
}

// values of the derivatives in eta direction
static void C_T_SV2_2D_DeriveEta(double xi, double eta, double *values)
{
  // first triangle
  if( (xi-eta)>0 && (xi+2*eta-1)<0 )
  {
    values[0] =  +8*xi +16*eta -6;
    values[1] =  -4*xi +16*eta -4;
    values[2] =  -4*xi + 4*eta +1;
    values[3] = 0;
    values[4] = 0;
    values[5] = 0;
    values[6] = -12*xi -48*eta+12;
    values[7] =  12*xi -24*eta;
    values[8] = 0;
    values[9] =        +36*eta -3;
  }
  else
    // second triangle
    if( (xi+2*eta-1)>=0 && (-2*xi-eta+1)<0 )
    {
      values[0] = 0;
      values[1] = 0;
      values[2] =  +8*xi + 4*eta - 5;
      values[3] = +20*xi +16*eta -12;
      values[4] =  +8*xi +16*eta -10;
      values[5] = 0;
      values[6] = 0;
      values[7] = -36*xi -24*eta +24;
      values[8] = -36*xi -48*eta +36;
      values[9] = +36*xi +36*eta -33;
    }
    else
    {
      values[0] =   8*xi +4*eta -3;
      values[1] = 0;
      values[2] = 0;
      values[3] = 0;
      values[4] =  -4*xi +4*eta -1;
      values[5] =  -4*xi -8*eta +4;
      values[6] = -12*xi;
      values[7] = 0;
      values[8] = +12*xi;
      values[9] = 0;
    }
}

// values of the derivatives in xi-xi  direction
static void C_T_SV2_2D_DeriveXiXi(double xi, double eta, double *values)
{
  // first triangle
  if( (xi-eta)>0 && (xi+2*eta-1)<0 )
  {
    values[0] =  4;
    values[1] = -8;
    values[2] =  4;
    values[3] =  0;
    values[4] =  0;
    values[5] =  0;
    values[6] =  0;
    values[7] =  0;
    values[8] =  0;
    values[9] =  0;
  }
  else
    // second triangle
    if( (xi+2*eta-1)>=0 && (-2*xi-eta+1)<0 )
    {
      values[0] =   0;
      values[1] =   0;
      values[2] =  16;
      values[3] =  16;
      values[4] =   4;
      values[5] =   0;
      values[6] =   0;
      values[7] = -48;
      values[8] = -24;
      values[9] =  36;
    }
    else
    {
      values[0] =  16;
      values[1] =   0;
      values[2] =   0;
      values[3] =   0;
      values[4] =   4;
      values[5] =  16;
      values[6] = -48;
      values[7] =   0;
      values[8] = -24;
      values[9] =  36;
    }
}

// values of the derivatives in xi-eta direction
static void C_T_SV2_2D_DeriveXiEta(double xi, double eta, double *values)
{
  // first triangle
  if( (xi-eta)>0 && (xi+2*eta-1)<0 )
  {
    values[0] =   8;
    values[1] =  -4;
    values[2] =  -4;
    values[3] =   0;
    values[4] =   0;
    values[5] =   0;
    values[6] = -12;
    values[7] =  12;
    values[8] =   0;
    values[9] =   0;
  }
  else
    // second triangle
    if( (xi+2*eta-1)>=0 && (-2*xi-eta+1)<0 )
    {
      values[0] =   0;
      values[1] =   0;
      values[2] =   8;
      values[3] =  20;
      values[4] =   8;
      values[5] =   0;
      values[6] =   0;
      values[7] = -36;
      values[8] = -36;
      values[9] =  36;
    }
    else
    {
      values[0] =   8;
      values[1] =   0;
      values[2] =   0;
      values[3] =   0;
      values[4] =  -4;
      values[5] =  -4;
      values[6] = -12;
      values[7] =   0;
      values[8] =  12;
      values[9] =   0;
    }
}

// values of the derivatives in eta-eta direction
static void C_T_SV2_2D_DeriveEtaEta(double xi, double eta, double *values)
{
  // first triangle
  if( (xi-eta)>0 && (xi+2*eta-1)<0 )
  {
    values[0] =  16;
    values[1] =  16;
    values[2] =   4;
    values[3] =   0;
    values[4] =   0;
    values[5] =   0;
    values[6] = -48;
    values[7] = -24;
    values[8] =   0;
    values[9] =  36;
  }
  else
    // second triangle
    if( (xi+2*eta-1)>=0 && (-2*xi-eta+1)<0 )
    {
      values[0] =   0;
      values[1] =   0;
      values[2] =   4;
      values[3] =  16;
      values[4] =  16;
      values[5] =   0;
      values[6] =   0;
      values[7] = -24;
      values[8] = -48;
      values[9] =  36;
    }
    else
    {
      values[0] =  4;
      values[1] =  0;
      values[2] =  0;
      values[3] =  0;
      values[4] =  4;
      values[5] = -8;
      values[6] =  0;
      values[7] =  0;
      values[8] =  0;
      values[9] =  0;
    }
}
