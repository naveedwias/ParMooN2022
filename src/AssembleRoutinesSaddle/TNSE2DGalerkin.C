#include "../../include/AssembleRoutinesSaddle/TNSE2DGalerkin.h"

#include <Hotfixglobal_AssembleNSE.h> // a temporary hotfix - check documentation!
#include <string>
#include <iostream>

#include <MooNMD_Io.h>

void TimeNSType1Galerkin(double Mult, double *coeff, double *param, double,
    double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz00 = Orig0[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz10+u2*ansatz01)*test00;
      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;
    }                            // endfor j
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

void TimeNSType2Galerkin(double Mult, double *coeff,
double *param, double,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB1T = LocMatrices[4];
  MatrixB2T = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz00 = Orig0[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz10+u2*ansatz01)*test00;
      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig1[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
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
void TimeNSType3Galerkin(double Mult, double *coeff,
double *param, double,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1  = LocMatrices[6];
  MatrixB2  = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];

    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz00 = Orig0[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j
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
void TimeNSType3GalerkinDD(double Mult, double *coeff,
double *param, double,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double **MatrixA21 = LocMatrices[2];
  double **MatrixA22 = LocMatrices[3];
  double **MatrixM11 = LocMatrices[4];
  double **MatrixM22 = LocMatrices[5];
  double **MatrixB1  = LocMatrices[6];
  double **MatrixB2  = LocMatrices[7];

  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];

  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];

  double *Orig0 = OrigValues[0];         // u
  double *Orig1 = OrigValues[1];         // p
  double *Orig2 = OrigValues[2];         // u_x
  double *Orig3 = OrigValues[3];         // u_y

  double c0 = coeff[0];                 // nu
  double c1 = coeff[1];                 // f1
  double c2 = coeff[2];                 // f2

  double u1 = param[0];                 // u1old
  double u2 = param[1];                 // u2old

  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double val;
  double val1 = 0;
  
  for(int i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(int j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz00 = Orig0[j];
      
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val1 = (u1*ansatz10+u2*ansatz01)*test00;
      val  = c0*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j
  }                              // endfor i

  for(int i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig1[i];

    for(int j=0;j<N_U;j++)
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
void TimeNSType4Galerkin(double Mult, double *coeff,
double *param, double,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22; // **MatrixA21, **MatrixA12;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
//  MatrixA12 = LocMatrices[1];
//  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[6];
  MatrixB2 = LocMatrices[7];
  MatrixB1T = LocMatrices[8];
  MatrixB2T = LocMatrices[9];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];    
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz00 = Orig0[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig1[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
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
void TimeNSType4GalerkinDD(double Mult, double *coeff,
double *param, double,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];  
  MatrixB1 = LocMatrices[6];
  MatrixB2 = LocMatrices[7];
  MatrixB1T = LocMatrices[8];
  MatrixB2T = LocMatrices[9];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // p
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

//   if(std::abs(u1)>0)
// cout << " u1 " << u1 << " u2 " << u2 <<endl;
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
   
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];
      ansatz00 = Orig0[j];

      val  = c0*(2*test10*ansatz10+test01*ansatz01);
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01);
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;      
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig1[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
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

void TimeNSType1_2NLGalerkin(double Mult, double *coeff,
double *param, double,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig2, *Orig3;
  int i,j,N_U;
  double c0;
  double u1, u2;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz10+u2*ansatz01)*test00;

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}

void TimeNSType3_4NLGalerkin(double Mult, double *coeff,
double *param, double,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig2, *Orig3;
  int i,j, N_U;
  double c0;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      //HOTFIX: Check the documentation!
      if(assemble_nse == Hotfixglobal_AssembleNSE::WITH_CONVECTION)
        val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}

void TimeNSType3_4NLGalerkinDD(double Mult, double *coeff,
double *param, double,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **)
{
  double **MatrixA11, **MatrixA22;
  double val, val1;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig2, *Orig3;
  int i,j, N_U;
  double c0;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u
  Orig2 = OrigValues[2];         // u_x
  Orig3 = OrigValues[3];         // u_y

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig2[i];
    test01 = Orig3[i];
    test00 = Orig0[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig2[j];
      ansatz01 = Orig3[j];

      val1 = (u1*ansatz10+u2*ansatz01)*test00;
      val  = c0*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// right-hand side ONLY, for NSE
// ======================================================================
void TimeNSRHS(double Mult, double *coeff,
double *, double,
double **OrigValues, int *N_BaseFuncts,
double ***, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test00;
  double *Orig0;
  int i, N_U;
  double c1, c2;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u

  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;
    //cout <<  Rhs1[i] << " " <<  Rhs2[i] << " ";
  }                              // endfor i
}
