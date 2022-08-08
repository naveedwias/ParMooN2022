#include <NSE2DSUPG.h>
#include <Database.h>
#include "MooNMD_Io.h"
#include <cmath>

// ======================================================================
// Type 2, SDFEM
// ======================================================================
void NSType2SUPG(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig8;
  int i,j,N_U, N_P;
  double c0, c1, c2, c;
  double u1, u2;
  double delta, ugrad;

  double delta0 = TDatabase::ParamDB->DELTA0;
  double delta1 = TDatabase::ParamDB->DELTA1;

  MatrixA = LocMatrices[0];
  MatrixB1 = LocMatrices[1];
  MatrixB2 = LocMatrices[2];
  MatrixB1T = LocMatrices[3];
  MatrixB2T = LocMatrices[4];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y
  Orig4 = OrigValues[4];         // p_x
  Orig5 = OrigValues[5];         // p_y
  Orig6 = OrigValues[6];         // u_xx
  Orig8 = OrigValues[8];         // u_yy

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  
  // for computational comparisons of Oseen problems
  if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE == OSEEN)
  {
    u1 = coeff[3];
    u2 = coeff[4];
    c  = coeff[5];
    param[0] = u1;
    param[1] = u2;
  }

  if(c0 < hK)
    delta = delta0*hK*hK;
  else
    delta = delta1*hK*hK;

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    ugrad  = delta * (u1*test10+u2*test01);
    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz20 = Orig6[j];
      ansatz02 = Orig8[j];
      // standard terms
      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      // SD term
      val += (-c0*(ansatz20+ansatz02)
        + (u1*ansatz10+u2*ansatz01) ) * ugrad;
      // for computational comparisons of Oseen problems
      if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE == OSEEN)
      {
        ansatz00 = Orig0[j];
        val += c * ansatz00 * test00;
        val += c * ansatz00 * ugrad;
      }

      MatrixRow[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      ansatz00 = Orig1[j];

      val  = -ansatz00 * test10;
      val +=  ansatz10 * ugrad;
      MatrixRow1[j] += Mult*val;

      val  = -ansatz00 * test01;
      val +=  ansatz01 * ugrad;

      MatrixRow2[j] += Mult*val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}

// ======================================================================
// Type 2, SDFEM
// ======================================================================
void NSType2NLSUPG(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig8;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, ugrad;

  static double delta0 = TDatabase::ParamDB->DELTA0;
  static double delta1 = TDatabase::ParamDB->DELTA1;

  MatrixA = LocMatrices[0];
  MatrixB1T = LocMatrices[1];
  MatrixB2T = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y
  Orig4 = OrigValues[4];         // p_x
  Orig5 = OrigValues[5];         // p_y
  Orig6 = OrigValues[6];         // u_xx
  Orig8 = OrigValues[8];         // u_yy

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  if(c0 < hK)
    delta = delta0*hK*hK ;
  else
    delta = delta1*hK*hK ;

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    ugrad  = delta * (u1*test10+u2*test01);

    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz20 = Orig6[j];
      ansatz02 = Orig8[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += (-c0*(ansatz20+ansatz02)
                                 // SD term
	      + (u1*ansatz10+u2*ansatz01) ) * ugrad;
      MatrixRow[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      ansatz00 = Orig1[j];

      val  = -ansatz00 * test10;
      val +=  ansatz10 * ugrad;
      MatrixRow1[j] += Mult*val;

      val  = -ansatz00 * test01;
      val +=  ansatz01 * ugrad;
      MatrixRow2[j] += Mult*val;
    }
  }                              // endfor i
}

// ======================================================================
// Type 4, SDFEM, (grad u, grad v)
// ======================================================================
void NSType4SUPG(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig8;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, ugrad, r=2.0, maxu, tau;

  double delta0 = TDatabase::ParamDB->DELTA0;
  double delta1 = TDatabase::ParamDB->DELTA1;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1 = LocMatrices[4];
  MatrixB2 = LocMatrices[5];
  MatrixB1T = LocMatrices[6];
  MatrixB2T = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y
  Orig4 = OrigValues[4];         // p_x
  Orig5 = OrigValues[5];         // p_y
  Orig6 = OrigValues[6];         // u_xx
  Orig8 = OrigValues[8];         // u_yy

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // for computational comparisons of Oseen problems
  if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE == OSEEN)
  {
    u1 = coeff[3];
    u2 = coeff[4];
  }

  maxu = std::abs(u1);
  if (std::abs(u2)>maxu)
      maxu = std::abs(u2);

  // parameter from [BBJL07]
  delta =  hK*hK/(r*r*(c0+1));
  delta =  delta0*delta;

  tau =  delta1 *(1+c0);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    ugrad  = delta * (u1*test10+u2*test01);

    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz20 = Orig6[j];
      ansatz02 = Orig8[j];

      // viscous term
      val  = c0*(test10*ansatz10+test01*ansatz01);
      // convective term
      val += (u1*ansatz10+u2*ansatz01)*test00;
      // SUPG term
      val +=  (-c0*(ansatz20+ansatz02)+ (u1*ansatz10+u2*ansatz01) ) * ugrad;
      // term in both diagonal blocks
      Matrix11Row[j] += Mult * (val+tau*test10*ansatz10);
      Matrix22Row[j] += Mult * (val+tau*test01*ansatz01);
      // div-div term
      val = tau * test10*ansatz01; 
      Matrix12Row[j] += Mult * val;

      val = tau * test01*ansatz10; 
      Matrix21Row[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      ansatz00 = Orig1[j];

      val  = -ansatz00 * test10;
      val +=  ansatz10 * ugrad;
      MatrixRow1[j] += Mult*val;

      val  = -ansatz00 * test01;
      val +=  ansatz01 * ugrad;
      MatrixRow2[j] += Mult*val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}

// ======================================================================
// Type 4, SDFEM, D(u):D(v)
// ======================================================================
void NSType4SUPGDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig8;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, ugrad;

  static double delta0 = TDatabase::ParamDB->DELTA0;
  static double delta1 = TDatabase::ParamDB->DELTA1;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixB1 = LocMatrices[4];
  MatrixB2 = LocMatrices[5];
  MatrixB1T = LocMatrices[6];
  MatrixB2T = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y
  Orig4 = OrigValues[4];         // p_x
  Orig5 = OrigValues[5];         // p_y
  Orig6 = OrigValues[6];         // u_xx
  Orig8 = OrigValues[8];         // u_yy

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  if(c0 < hK)
    delta = delta0*hK*hK;
  else
    delta = delta1*hK*hK;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    ugrad  = delta * (u1*test10+u2*test01);

    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz20 = Orig6[j];
      ansatz02 = Orig8[j];

      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += (-c0*(ansatz20+ansatz02)
                                 // SD term
        + (u1*ansatz10+u2*ansatz01) ) * ugrad;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += (-c0*(ansatz20+ansatz02)
                                 // SD term
        + (u1*ansatz10+u2*ansatz01) ) * ugrad;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      ansatz00 = Orig1[j];

      val  = -ansatz00 * test10;
      val +=  ansatz10 * ugrad;
      MatrixRow1[j] += Mult*val;

      val  = -ansatz00 * test01;
      val +=  ansatz01 * ugrad;
      MatrixRow2[j] += Mult*val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}

// ======================================================================
// Type 4, SDFEM, (grad u, grad v)
// ======================================================================
void NSType4NLSUPG(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig8;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, ugrad, r=2.0, maxu, tau;

  double delta0 = TDatabase::ParamDB->DELTA0;
  double delta1 = TDatabase::ParamDB->DELTA1;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  MatrixB1T = LocMatrices[2];
  MatrixB2T = LocMatrices[3];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y
  Orig4 = OrigValues[4];         // p_x
  Orig5 = OrigValues[5];         // p_y
  Orig6 = OrigValues[6];         // u_xx
  Orig8 = OrigValues[8];         // u_yy

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // for computational comparisons of Oseen problems
  if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE == OSEEN)
  {
      u1 = coeff[3];
      u2 = coeff[4];
  }

  maxu = std::abs(u1);
  if (std::abs(u2)>maxu)
      maxu = std::abs(u2);

  // parameter from [BBJL07]
  delta =  hK*hK/(r*r*(c0+1));
  delta =  delta0*delta;

  tau =  delta1 *(1+c0);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    ugrad  = delta * (u1*test10+u2*test01);

    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz20 = Orig6[j];
      ansatz02 = Orig8[j];

      // viscous term
      val  = c0*(test10*ansatz10+test01*ansatz01);
      // convective term
      val += (u1*ansatz10+u2*ansatz01)*test00;
      // SUPG term
      val +=  (-c0*(ansatz20+ansatz02)+ (u1*ansatz10+u2*ansatz01) ) * ugrad;
      // term in both diagonal blocks
      Matrix11Row[j] += Mult * (val+tau*test10*ansatz10);
      Matrix22Row[j] += Mult * (val+tau*test01*ansatz01);
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      ansatz00 = Orig1[j];

      val  = -ansatz00 * test10;
      val +=  ansatz10 * ugrad;
      MatrixRow1[j] += Mult*val;

      val  = -ansatz00 * test01;
      val +=  ansatz01 * ugrad;
      MatrixRow2[j] += Mult*val;
    }
  }                              // endfor i
}

// ======================================================================
// Type 4, SDFEM, D(u):D(v)
// ======================================================================
void NSType4NLSUPGDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double *Orig3, *Orig4, *Orig5;
  double *Orig6, *Orig8;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, ugrad;

  static double delta0 = TDatabase::ParamDB->DELTA0;
  static double delta1 = TDatabase::ParamDB->DELTA1;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  MatrixB1T = LocMatrices[2];
  MatrixB2T = LocMatrices[3];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y
  Orig4 = OrigValues[4];         // p_x
  Orig5 = OrigValues[5];         // p_y
  Orig6 = OrigValues[6];         // u_xx
  Orig8 = OrigValues[8];         // u_yy

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  if(c0 < hK)
    delta = delta0*hK*hK;
  else
    delta = delta1*hK*hK;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    ugrad  = delta * (u1*test10+u2*test01);

    Rhs1[i] += Mult*(test00+ugrad)*c1;
    Rhs2[i] += Mult*(test00+ugrad)*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz20 = Orig6[j];
      ansatz02 = Orig8[j];

      val  = 2*c0*(test10*ansatz10+0.5*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += (-c0*(ansatz20+ansatz02)
                                 // SD term
        + (u1*ansatz10+u2*ansatz01) ) * ugrad;
      Matrix11Row[j] += Mult * val;

      val  = 2*c0*(0.5*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += (-c0*(ansatz20+ansatz02)
                                 // SD term
        + (u1*ansatz10+u2*ansatz01) ) * ugrad;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      ansatz00 = Orig1[j];

      val  = -ansatz00 * test10;
      val +=  ansatz10 * ugrad;
      MatrixRow1[j] += Mult*val;

      val  = -ansatz00 * test01;
      val +=  ansatz01 * ugrad;
      MatrixRow2[j] += Mult*val;
    }
  }                              // endfor i
}
