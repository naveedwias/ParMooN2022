// ***********************************************************************
// SV1 element, discont. pressure, 2D
// ***********************************************************************

// base function values
static void D_T_SV1_2D_Funct(double xi, double eta, double *values)
{
  // first triangle
  if( (xi-eta)>0 && (xi+2*eta-1)<0 )
  {
    values[0] = 1;
    values[1] = 2*xi+eta-1;
    values[2] = xi+5*eta-1;
    values[3] = 0;
    values[4] = 0;
    values[5] = 0;
    values[6] = 0;
    values[7] = 0;
    values[8] = 0;
  }
  else
    // second triangle
    if( (xi+2*eta-1)>=0 && (-2*xi-eta+1)<0 )
    {
      values[0] = 0;
      values[1] = 0;
      values[2] = 0;
      values[3] = 1;
      values[4] = -xi+eta;
      values[5] = -5*xi-4*eta+4;
      values[6] = 0;
      values[7] = 0;
      values[8] = 0;
    }
    else
    {
      values[0] = 0;
      values[1] = 0;
      values[2] = 0;
      values[3] = 0;
      values[4] = 0;
      values[5] = 0;
      values[6] = 1;
      values[7] = -xi-2*eta+1;
      values[8] = 4*xi-eta;
    }
}

// values of the derivatives in xi direction
static void D_T_SV1_2D_DeriveXi(double xi, double eta, double *values)
{
  // first triangle
  if( (xi-eta)>0 && (xi+2*eta-1)<0 )
  {
    values[0] = 0;
    values[1] = 2;
    values[2] = 1;
    values[3] = 0;
    values[4] = 0;
    values[5] = 0;
    values[6] = 0;
    values[7] = 0;
    values[8] = 0;
  }
  else
    // second triangle
    if( (xi+2*eta-1)>=0 && (-2*xi-eta+1)<0 )
    {
      values[0] = 0;
      values[1] = 0;
      values[2] = 0;
      values[3] = 0;
      values[4] = -1;
      values[5] = -5;
      values[6] = 0;
      values[7] = 0;
      values[8] = 0;
    }
    else
    {
      values[0] = 0;
      values[1] = 0;
      values[2] = 0;
      values[3] = 0;
      values[4] = 0;
      values[5] = 0;
      values[6] = 0;
      values[7] = -1;
      values[8] = 4;
    }
}

// values of the derivatives in eta direction
static void D_T_SV1_2D_DeriveEta(double xi, double eta, double *values)
{
  // first triangle
  if( (xi-eta)>0 && (xi+2*eta-1)<0 )
  {
    values[0] = 0;
    values[1] = 1;
    values[2] = 5;
    values[3] = 0;
    values[4] = 0;
    values[5] = 0;
    values[6] = 0;
    values[7] = 0;
    values[8] = 0;
  }
  else
    // second triangle
    if( (xi+2*eta-1)>=0 && (-2*xi-eta+1)<0 )
    {
      values[0] = 0;
      values[1] = 0;
      values[2] = 0;
      values[3] = 0;
      values[4] = 1;
      values[5] = -4;
      values[6] = 0;
      values[7] = 0;
      values[8] = 0;
    }
    else
    {
      values[0] = 0;
      values[1] = 0;
      values[2] = 0;
      values[3] = 0;
      values[4] = 0;
      values[5] = 0;
      values[6] = 0;
      values[7] = -2;
      values[8] = -1;
    }
}

// values of the derivatives in xi-xi  direction
static void D_T_SV1_2D_DeriveXiXi(double, double, double *values)
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
}

// values of the derivatives in xi-eta direction
static void D_T_SV1_2D_DeriveXiEta(double, double, double *values)
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
}

// values of the derivatives in eta-eta direction
static void D_T_SV1_2D_DeriveEtaEta(double, double, double *values)
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
}
