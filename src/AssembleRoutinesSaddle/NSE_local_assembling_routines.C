#include "NSE_local_assembling_routines.h"
#include <MooNMD_Io.h>

template<int d>
void NSResistanceMassMatrixSingle(double Mult, const double *coeff,
                                  const double *, double,
                                  const double **OrigValues,
                                  const int *N_BaseFuncts,
                                  double ***LocMatrices, double **)
{
  double ** MatrixA = LocMatrices[0];
  int N_U = N_BaseFuncts[0];
  const double * u = OrigValues[0];
  double sigma = d == 2 ? coeff[4] : coeff[5];

  for(int i = 0; i < N_U; i++)
  {
    double test = u[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz = u[j];
      MatrixA[i][j] += Mult * sigma * (test * ansatz);
    }
  }
}

template <int d>
void NSResistanceMassMatrix(double Mult, const double *coeff, const double *,
                            double, const double**OrigValues,
                            const int *N_BaseFuncts, double ***LocMatrices,
                            double **)
{
  double ** MatrixA11 = LocMatrices[0];
  double ** MatrixA22 = LocMatrices[d+1];
  double ** MatrixA33 = d == 2 ? nullptr : LocMatrices[8];
  int N_U = N_BaseFuncts[0];
  const double * u = OrigValues[0];
  double sigma = d == 2 ? coeff[4] : coeff[5];

  for(int i = 0; i < N_U; i++)
  {
    double test = u[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz = u[j];
      double val = Mult * sigma * (test * ansatz);
      MatrixA11[i][j] += val;
      MatrixA22[i][j] += val;
      if(d == 3)
        MatrixA33[i][j] += val;
    }
  }
}

template<int d>
void NSLaplaceGradGradSingle(double Mult, const double *coeff, const double *,
                             double, const double **OrigValues,
                             const int *N_BaseFuncts, double ***LocMatrices,
                             double **)
{
  double ** MatrixA = LocMatrices[0];
  int N_U = N_BaseFuncts[0];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = d == 2 ? nullptr : OrigValues[4];
  double nu = coeff[0]; // = 1/reynolds_number

  for(int i = 0; i < N_U; i++)
  {
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0. : u_z[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      MatrixA[i][j] += Mult * nu * (test_x * ansatz_x + test_y * ansatz_y
                                 + test_z * ansatz_z);
    }
  }
}


template <int d>
void NSLaplaceGradGrad(double Mult, const double *coeff, const double *, double,
                       const double**OrigValues, const int *N_BaseFuncts,
                       double ***LocMatrices, double **)
{
  double ** MatrixA11 = LocMatrices[0];
  double ** MatrixA22 = LocMatrices[d+1];
  double ** MatrixA33 = d == 2 ? nullptr : LocMatrices[8];
  int N_U = N_BaseFuncts[0];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = d == 2 ? nullptr : OrigValues[4];
  double nu = coeff[0]; // = 1/reynolds_number

  for(int i = 0; i < N_U; i++)
  {
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0. : u_z[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      double val = Mult * nu* ( test_x * ansatz_x + test_y * ansatz_y
                             + test_z * ansatz_z);
      MatrixA11[i][j] += val;
      MatrixA22[i][j] += val;
      if(d == 3)
        MatrixA33[i][j] += val;
    }
  }
}

template <int d>
void NSLaplaceDeformation(double Mult, const double *coeff, const double *,
                          double, const double**OrigValues,
                          const int *N_BaseFuncts, double ***LocMatrices,
                          double **)
{
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double ** MatrixA13 = d == 2 ? nullptr : LocMatrices[2];
  double **MatrixA21 = LocMatrices[d];
  double **MatrixA22 = LocMatrices[d+1];
  double **MatrixA23 = d == 2 ? nullptr : LocMatrices[d+2];
  double ** MatrixA31 = d == 2 ? nullptr : LocMatrices[6];
  double ** MatrixA32 = d == 2 ? nullptr : LocMatrices[7];
  double ** MatrixA33 = d == 2 ? nullptr : LocMatrices[8];
  int N_U = N_BaseFuncts[0];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = d == 2 ? nullptr : OrigValues[4];
  double nu = coeff[0];
  for(int i = 0; i < N_U; i++)
  {
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0. : u_z[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      MatrixA11[i][j] += Mult * 2*nu*(test_x*ansatz_x + 0.5*test_y*ansatz_y
                                     +0.5*test_z*ansatz_z);
      MatrixA12[i][j] += Mult * nu*(test_y*ansatz_x);
      if(d == 3)
        MatrixA13[i][j] += Mult * nu*(test_z*ansatz_x);
      MatrixA21[i][j] += Mult * nu*(test_x*ansatz_y);
      MatrixA22[i][j] += Mult * 2*nu*(0.5*test_x*ansatz_x+test_y*ansatz_y
                                     +0.5*test_z*ansatz_z);
      if(d == 3)
      {
        MatrixA23[i][j] += Mult * nu*(test_z*ansatz_y);
        MatrixA31[i][j] += Mult * nu*(test_x*ansatz_z);
        MatrixA32[i][j] += Mult * nu*(test_y*ansatz_z);
        MatrixA33[i][j] += Mult * 2*nu*(0.5*test_x*ansatz_x+0.5*test_y*ansatz_y
                                       +test_z*ansatz_z);
      }
    }
  }
}

template <int d>
void NSDivergenceBlocks(double Mult, const double *, const double *, double,
                        const double **OrigValues, const int *N_BaseFuncts,
                        double ***LocMatrices, double **, int sign)
{
  double ** MatrixB1 = LocMatrices[d == 2 ? 5 : 10];
  double ** MatrixB2 = LocMatrices[d == 2 ? 6 : 11];
  double ** MatrixB3 = d == 2 ? nullptr : LocMatrices[12];
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  const double * p = OrigValues[1];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = d == 2 ? nullptr : OrigValues[4];
  for(int i = 0; i < N_P; i++)
  {
    double test = p[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      MatrixB1[i][j] -= sign * Mult * test * ansatz_x;
      MatrixB2[i][j] -= sign * Mult * test * ansatz_y;
      if(d == 3)
        MatrixB3[i][j] -=  Mult * sign * test * ansatz_z;
    }
  }
}

template <int d>
void NSGradientBlocks(double Mult, const double *, const double *, double,
                      const double**OrigValues, const int *N_BaseFuncts,
                      double ***LocMatrices, double **)
{
  double ** MatrixB1T = LocMatrices[d == 2 ? 7 : 13];
  double ** MatrixB2T = LocMatrices[d == 2 ? 8 : 14];
  double ** MatrixB3T = d == 2 ? nullptr : LocMatrices[15];
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  const double * p = OrigValues[1];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = d == 2 ? nullptr : OrigValues[4];
  for(int i = 0; i < N_U; i++)
  {
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0. : u_z[i];
    for(int j = 0; j < N_P; j++)
    {
      double ansatz = p[j];
      MatrixB1T[i][j] -= Mult * ansatz * test_x;
      MatrixB2T[i][j] -= Mult * ansatz * test_y;
      if(d == 3)
        MatrixB3T[i][j] -= Mult * ansatz * test_z;
    }
  }
}

template <int d>
void NSRightHandSide(double Mult, const double *coeff, const double *, double,
                     const double **OrigValues, const int *N_BaseFuncts,
                     double ***, double **LocRhs, int sign)
{
  double * Rhs1 = LocRhs[0];
  double * Rhs2 = LocRhs[1];
  double * Rhs3 = d == 2 ? nullptr : LocRhs[2];
  double * Rhs_div = LocRhs[d];
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  const double * u = OrigValues[0];
  const double * p = OrigValues[1];
  double f1 = coeff[1];
  double f2 = coeff[2];
  double f3 = d == 2 ? 0. : coeff[3];
  double g = coeff[d+1]; // divergence
  for(int i = 0; i < N_U; i++)
  {
    double test = u[i];
    Rhs1[i] += Mult * test * f1;
    Rhs2[i] += Mult * test * f2;
    if(d == 3)
      Rhs3[i] += Mult * test * f3;
  }
  for(int i = 0; i < N_P; i++)
  {
    double test = p[i];
    Rhs_div[i] -= Mult * sign * test * g;
  }
}

///////////////////////////////////////////////////////////////////////////////
template <int d>
void NSLaplaceGradGrad_dg(double Mult, const double *coeff, const double *,
                          double, const double **OrigValues,
                          const int *N_BaseFuncts, double ***LocMatrices,
                          double **LocRhs)
{
  double val;
  double ansatz;
  double test, test_x, test_y, test_z;
  double test_x_100, test_y_100, test_z_100;
  double test_x_010, test_y_010, test_z_010;
  double test_x_001, test_y_001, test_z_001;
  double ansatz_x_100, ansatz_y_100, ansatz_z_100;
  double ansatz_x_010, ansatz_y_010, ansatz_z_010;
  double ansatz_x_001, ansatz_y_001, ansatz_z_001;
  double test_div;

  // ( A  B1 )   ( 0 2 )
  // ( B2 C  )   ( 3 1 )

  double **MatrixA = LocMatrices[0];
//  double **MatrixC = LocMatrices[1];
  double **MatrixB1 = LocMatrices[2];
  double **MatrixB2 = LocMatrices[3];

  double *Rhs0 = LocRhs[0];
  double *Rhs1 = LocRhs[1];

  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];

  const double *Orig0 = OrigValues[0];   // u
  const double *Orig1 = OrigValues[1];   // p
  const double *Orig2 = OrigValues[2];   // u_x
  const double *Orig3 = OrigValues[3];   // u_y
  const double *Orig4 = OrigValues[4];   // u_z (unused in 2D)


  double c0 = coeff[0];  // sigma
  double f1 = coeff[1];  // f1
  double f2 = coeff[2];  // f2
  double f3 = coeff[3];  // unused 2D, f3 in 3D
  double g = d == 2 ? coeff[3] : coeff[4];  // g(x,y)

  // A, B1, B2
  for(int i=0;i<N_U;i++)
  {
    // A:
    test_x = Orig0[i];
    test_y = Orig0[N_U+i];
    test_z = d == 2 ? 0. : Orig0[2*N_U+i];

    test_x_100 = Orig2[i];
    test_y_100 = Orig2[N_U+i];
    test_z_100 = d == 2 ? 0. : Orig2[2*N_U+i];
    test_x_010 = Orig3[i];
    test_y_010 = Orig3[N_U+i];
    test_z_010 = d == 2 ? 0. : Orig3[2*N_U+i];
    test_x_001 = d == 2 ? 0. : Orig4[N_U];
    test_y_001 = d == 2 ? 0. : Orig4[N_U+i];
    test_z_001 = d == 2 ? 0. : Orig4[2*N_U+i];

    Rhs0[i] += Mult*(f1*test_x + f2*test_y + f3*test_z);

    for(int j=0;j<N_U;j++)
    {
      ansatz_x_100 = Orig2[j];
      ansatz_y_100 = Orig2[N_U+j];
      ansatz_z_100 = d == 2 ? 0. : Orig2[2*N_U+j];
      ansatz_x_010 = Orig3[j];
      ansatz_y_010 = Orig3[N_U+j];
      ansatz_z_010 = d == 2 ? 0. : Orig3[2*N_U+j];
      ansatz_x_001 = d == 2 ? 0. : Orig4[j];
      ansatz_y_001 = d == 2 ? 0. : Orig4[N_U+j];
      ansatz_z_001 = d == 2 ? 0. : Orig4[2*N_U+j];

      // A: (grad u, grad v)
      val = test_x_100*ansatz_x_100 + test_y_100*ansatz_y_100
          + test_z_100*ansatz_z_100;
      val += test_x_010*ansatz_x_010 + test_y_010*ansatz_y_010
           + test_z_010*ansatz_z_010;
      val += test_x_001*ansatz_x_001 + test_y_001*ansatz_y_001
           + test_z_001*ansatz_z_001;
      MatrixA[i][j] += Mult * c0 * val;
    }
    // B1, B2:
    test_div = test_x_100 + test_y_010 + test_z_001;
    for(int j=0;j<N_P;j++)
    {
      ansatz = Orig1[j];
      val = Mult*test_div*ansatz;
      // (p div v)
      MatrixB1[i][j] -= val;
      // (q, div u)
      MatrixB2[j][i] -= val; // slow (consider moving this into another loop)
    }
  }

  for(int i=0;i<N_P;i++)
  {
    test = Orig1[i];
    // assemble rhs: div u = g
    // rhs: -(g,q)
    Rhs1[i] -= Mult*test*g;
  }
}

///////////////////////////////////////////////////////////////////////////////
// nonlinear term
template <int d>
void NSNonlinearTerm_convective_Single(
  double Mult, const double *, const double *param, double,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **)
{
  double ** MatrixA = LocMatrices[0];
  int N_U = N_BaseFuncts[0];
  const double * u   = OrigValues[0];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = d == 2 ? nullptr : OrigValues[4];
  double u1 = param[0];
  double u2 = param[1];
  double u3 = d == 2 ? 0. : param[2];
  for(int i = 0; i < N_U; i++)
  {
    double test = u[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      MatrixA[i][j] += Mult * (u1*ansatz_x + u2*ansatz_y + u3*ansatz_z) * test;
    }
  }
}

template <int d>
void NSNonlinearTerm_convective(
  double Mult, const double *, const double *param, double,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **)
{
  double ** MatrixA11 = LocMatrices[0];
  double ** MatrixA22 = LocMatrices[d+1];
  double ** MatrixA33 = d == 2 ? nullptr : LocMatrices[8];
  int N_U = N_BaseFuncts[0];
  const double * u = OrigValues[0];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = d == 2 ? nullptr : OrigValues[4];
  double u1 = param[0];
  double u2 = param[1];
  double u3 = d == 2 ? 0. : param[2];
  for(int i = 0; i < N_U; i++)
  {
    double test = u[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      double val = Mult * (u1*ansatz_x + u2*ansatz_y + u3*ansatz_z) * test;
      MatrixA11[i][j] += val;
      MatrixA22[i][j] += val;
      if(d == 3)
        MatrixA33[i][j] += val;
    }
  }
}





///////////////////////////////////////////////////////////////////////////////
template <int d>
void NSNonlinearTerm_convective_dg(double Mult, const double *, const double *param,
                          double, const double **OrigValues,
                          const int *N_BaseFuncts, double ***LocMatrices,
                          double **)
{
  double val;
  //double ansatz;
  double test_x, test_y, test_z;
  //double test_x_100, test_y_100, test_z_100;
  //double test_x_010, test_y_010, test_z_010;
  //double test_x_001, test_y_001, test_z_001;
  double ansatz_x_100, ansatz_y_100, ansatz_z_100;
  double ansatz_x_010, ansatz_y_010, ansatz_z_010;
  double ansatz_x_001, ansatz_y_001, ansatz_z_001;
  //double test_div;
  
  // ( A  B1 )   ( 0 2 )
  // ( B2 C  )   ( 3 1 )
  
  double **MatrixA = LocMatrices[0];
//  double **MatrixC = LocMatrices[1];
//  double **MatrixB1 = LocMatrices[2];
//  double **MatrixB2 = LocMatrices[3];
  
  
  int N_U = N_BaseFuncts[0];
 // int N_P = N_BaseFuncts[1];

  const double *Orig0 = OrigValues[0];   // u
//  const double *Orig1 = OrigValues[1];   // p
  const double *Orig2 = OrigValues[2];   // u_x
  const double *Orig3 = OrigValues[3];   // u_y
  const double *Orig4 = OrigValues[4];   // u_z (unused in 2D)
  

//  double c0 = coeff[0];  // sigma
//  double f1 = coeff[1];  // f1
//  double f2 = coeff[2];  // f2
//  double f3 = coeff[3];  // unused 2D, f3 in 3D
//  double g = d == 2 ? coeff[3] : coeff[4];  // g(x,y)

  double u1 = param[0];
  double u2 = param[1];
  double u3 = d == 2 ? 0. : param[2];

  for(int i=0;i<N_U;i++)
  {

    test_x = Orig0[i];
    test_y = Orig0[N_U+i];
    test_z = d == 2 ? 0. : Orig0[2*N_U+i];
    
    
    for(int j=0;j<N_U;j++)
    {
      ansatz_x_100 = Orig2[j];
      ansatz_y_100 = Orig2[N_U+j];
      ansatz_z_100 = d == 2 ? 0. : Orig2[2*N_U+j];
      ansatz_x_010 = Orig3[j];
      ansatz_y_010 = Orig3[N_U+j];
      ansatz_z_010 = d == 2 ? 0. : Orig3[2*N_U+j];
      ansatz_x_001 = d == 2 ? 0. : Orig4[j];
      ansatz_y_001 = d == 2 ? 0. : Orig4[N_U+j];
      ansatz_z_001 = d == 2 ? 0. : Orig4[2*N_U+j];
       
      val = (u1 * ansatz_x_100 + u2 * ansatz_x_010) * test_x;
      val += (u1 * ansatz_y_100 + u2 * ansatz_y_010) * test_y;      
      
      if (d == 3)
      {
      	val += u3 * ansatz_x_001 * test_x;
      	val += u3 * ansatz_y_001 * test_y;
      	val += (u1 * ansatz_z_100 + u2 * ansatz_z_010 + u1 * ansatz_z_001) * test_z; 
      }
      
      MatrixA[i][j] += Mult * val;
    }
  } 
}

template <int d>
void NSNonlinearTerm_skew_symmetric_Single(
  double Mult, const double *, const double *param, double,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **)
{
  double ** MatrixA11 = LocMatrices[0];
  int N_U = N_BaseFuncts[0];
  const double * u = OrigValues[0];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = d == 2 ? nullptr : OrigValues[4];
  double u1 = param[0];
  double u2 = param[1];
  double u3 = d == 2 ? 0. : param[2];
  for(int i = 0; i < N_U; i++)
  {
    double test = u[i];
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0. : u_z[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz = u[j];
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      double val = (u1*ansatz_x + u2*ansatz_y + u3*ansatz_z) * test;
      val -= (u1*test_x + u2*test_y + u3*test_z) * ansatz;
      val *= 0.5 * Mult;
      MatrixA11[i][j] += val;
    }
  }
}

template <int d>
void NSNonlinearTerm_skew_symmetric(
  double Mult, const double *, const double *param, double,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **)
{
  double ** MatrixA11 = LocMatrices[0];
  double ** MatrixA22 = LocMatrices[d+1];
  double ** MatrixA33 = d == 2 ? nullptr : LocMatrices[8];
  int N_U = N_BaseFuncts[0];
  const double * u = OrigValues[0];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = d == 2 ? nullptr : OrigValues[4];
  double u1 = param[0];
  double u2 = param[1];
  double u3 = d == 2 ? 0. : param[2];
  for(int i = 0; i < N_U; i++)
  {
    double test = u[i];
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0. : u_z[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz = u[j];
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      double val = (u1*ansatz_x + u2*ansatz_y + u3*ansatz_z) * test;
      val -= (u1*test_x + u2*test_y + u3*test_z) * ansatz;
      val *= 0.5 * Mult;
      MatrixA11[i][j] += val;
      MatrixA22[i][j] += val;
      if(d == 3)
        MatrixA33[i][j] += val;
    }
  }
}

template <int d>
void NSNonlinearTerm_rotational(
  double Mult, const double *, const double *param, double,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **)
{
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double ** MatrixA13 = d == 2 ? nullptr : LocMatrices[2];
  double **MatrixA21 = LocMatrices[d];
  double **MatrixA22 = LocMatrices[d+1];
  double **MatrixA23 = d == 2 ? nullptr : LocMatrices[d+2];
  double ** MatrixA31 = d == 2 ? nullptr : LocMatrices[6];
  double ** MatrixA32 = d == 2 ? nullptr : LocMatrices[7];
  double ** MatrixA33 = d == 2 ? nullptr : LocMatrices[8];
  int N_U = N_BaseFuncts[0];
  const double * u = OrigValues[0];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = d == 2 ? nullptr : OrigValues[4];
  double u1 = param[0];
  double u2 = param[1];
  double u3 = d == 2 ? 0. : param[2];
  for(int i = 0; i < N_U; i++)
  {
    double test = u[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      MatrixA11[i][j] += Mult * (u2*ansatz_y + u3*ansatz_z) * test;
      MatrixA12[i][j] -= Mult * (u2*ansatz_x) * test;
      if(d == 3)
        MatrixA13[i][j] -= Mult * (u3*ansatz_x) * test;
      MatrixA21[i][j] -= Mult * (u1*ansatz_y) * test;
      MatrixA22[i][j] += Mult * (u1*ansatz_x + u3*ansatz_z) * test;
      if(d == 3)
      {
        MatrixA23[i][j] -= Mult * (u3*ansatz_y) * test;
        MatrixA31[i][j] -= Mult * (u1*ansatz_z) * test;
        MatrixA32[i][j] -= Mult * (u2*ansatz_z) * test;
        MatrixA33[i][j] += Mult * (u1*ansatz_x + u2*ansatz_y) * test;
      }
    }
  }
}

template <int d>
void NSNonlinearTerm_emac(
  double Mult, const double *, const double *param, double,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **)
{
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double ** MatrixA13 = d == 2 ? nullptr : LocMatrices[2];
  double **MatrixA21 = LocMatrices[d];
  double **MatrixA22 = LocMatrices[d+1];
  double **MatrixA23 = d == 2 ? nullptr : LocMatrices[d+2];
  double ** MatrixA31 = d == 2 ? nullptr : LocMatrices[6];
  double ** MatrixA32 = d == 2 ? nullptr : LocMatrices[7];
  double ** MatrixA33 = d == 2 ? nullptr : LocMatrices[8];
  int N_U = N_BaseFuncts[0];
  const double * u = OrigValues[0];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = d == 2 ? nullptr : OrigValues[4];
  double u1 = param[0];
  double u2 = param[1];
  double u3 = d == 2 ? 0. : param[2];
  for(int i = 0; i < N_U; i++)
  {
    double test = Mult * u[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      MatrixA11[i][j] += (3*u1*ansatz_x + (ansatz_y*u2 + ansatz_z*u3))*test;
      MatrixA12[i][j] += (ansatz_x*u2 + ansatz_y*u1) * test;
      if(d == 3)
        MatrixA13[i][j] += (ansatz_x*u3 + ansatz_z*u1) * test;
      MatrixA21[i][j] += (ansatz_y*u1 + ansatz_x*u2) * test;
      MatrixA22[i][j] += (3*ansatz_y*u2 + (ansatz_x*u1 + ansatz_z*u3))*test;
      if(d == 3)
      {
        MatrixA23[i][j] += (ansatz_y*u3 + ansatz_z*u2) * test;
        MatrixA31[i][j] += (ansatz_z*u1 + ansatz_x*u3) * test;
        MatrixA32[i][j] += (ansatz_z*u2 + ansatz_y*u3) * test;
        MatrixA33[i][j] += (3*ansatz_z*u3 + (ansatz_x*u1+ansatz_y*u2))*test;
      }
    }
  }
}

template <int d>
void NSNonlinearTerm_divergence(double Mult, const double *, const double *param, double,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **)
{
  double ** MatrixA11 = LocMatrices[0];
  double ** MatrixA22 = LocMatrices[d+1];
  double ** MatrixA33 = d == 2 ? nullptr : LocMatrices[8];
  int N_U = N_BaseFuncts[0];
  const double * u = OrigValues[0];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = d == 2 ? nullptr : OrigValues[4];

  const double u1 = param[0];
  const double u1x = param[1];
  const double u2 = param[d+1];
  const double u2y = param[d+3];
  const double u3 = d == 2 ? 0. : param[2*d+2];
  const double u3z = d == 2 ? 0. : param[2*d+5];

  double divu = (u1x + u2y + u3z)/2.;


  for(int i = 0; i < N_U; i++)
  {
    double test = u[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz   = u[j];
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      double val = Mult * ((u1*ansatz_x + u2*ansatz_y + u3*ansatz_z) * test + divu * ansatz * test);
      MatrixA11[i][j] += val;
      MatrixA22[i][j] += val;
      if(d == 3)
        MatrixA33[i][j] += val;
    }
  }
}

template <int d>
void NSCoriolis(double Mult, const double *coeff, const double *, double,
                const double **OrigValues, const int *N_BaseFuncts,
                double ***LocMatrices, double **)
{
  double** MatrixA12 = LocMatrices[1];
  double** MatrixA13 = LocMatrices[2];
  double** MatrixA21 = LocMatrices[3];
  double** MatrixA23 = LocMatrices[5];
  double** MatrixA31 = LocMatrices[6];
  double** MatrixA32 = LocMatrices[7];

  int N_U = N_BaseFuncts[0];

  const double* Orig3 = OrigValues[3]; // u

  double Omega1 = coeff[6];
  double Omega2 = coeff[7];
  double Omega3 = coeff[8];

  for(int i=0;i<N_U;i++)
  {
    const double test000 = Orig3[i];
    for(int j=0;j<N_U;j++)
    {
      const double ansatz000 = Orig3[j];
      const double val1 = Mult * ansatz000 * test000;
      MatrixA12[i][j] -= Omega3 * val1;
      MatrixA13[i][j] += Omega2 * val1;
      MatrixA21[i][j] += Omega3 * val1;
      MatrixA23[i][j] -= Omega1 * val1;
      MatrixA31[i][j] -= Omega2 * val1;
      MatrixA32[i][j] += Omega1 * val1;
    } // endfor j
  } // endfor i
}

template <int d>
void NSRightHandSideExplicitNL(double Mult, const double *, const double *param,
                               double, const double **OrigValues,
                               const int *N_BaseFuncts, double ***,
                               double **LocRhs)
{
  double * Rhs1 = LocRhs[0];
  double * Rhs2 = LocRhs[1];
  double * Rhs3 = d == 2 ? nullptr : LocRhs[2];
  int N_U = N_BaseFuncts[0];
  const double * u = OrigValues[0];
  double u1 = param[0];
  double u1x = param[1];
  double u1y = param[2];
  double u1z = d == 2 ? 0. : param[3];
  double u2 = param[d+1];
  double u2x = param[d+2];
  double u2y = param[d+3];
  double u2z = d == 2 ? 0. : param[d+4];
  double u3 = d == 2 ? 0. : param[2*d+2];
  double u3x = d == 2 ? 0. : param[2*d+3];
  double u3y = d == 2 ? 0. : param[2*d+4];
  double u3z = d == 2 ? 0. : param[2*d+5];

  for(int i = 0; i < N_U; i++)
  {
    double test = u[i];
    Rhs1[i] += Mult * test * (u1*u1x + u2*u1y + u3*u1z);
    Rhs2[i] += Mult * test * (u1*u2x + u2*u2y + u3*u2z);
    if(d == 3)
      Rhs3[i] += Mult * test * (u1*u3x + u2*u3y + u3*u3z);
  }
}


// general version (for Navier-Stokes and Brinkman problems)
double compute_GradDiv_delta(double stab, double nu,
                             double sigma, double characteristic_length)
{
  return stab *  (nu + sigma * characteristic_length * characteristic_length);
}

///////////////////////////////////////////////////////////////////////////////
// stabilizations
template <int d>
void NSGradDiv(double Mult, const double *coeff, const double *, double,
            const double **OrigValues, const int *N_BaseFuncts,
            double ***LocMatrices, double **, double stab,
            double characteristic_length)
{
  double ** MatrixA11 = LocMatrices[0];
  double ** MatrixA12 = LocMatrices[1];
  double ** MatrixA13 = d == 2 ? nullptr : LocMatrices[2];

  double ** MatrixA21 = LocMatrices[d == 2 ? 2 : 3];
  double ** MatrixA22 = LocMatrices[d == 2 ? 3 : 4];
  double ** MatrixA23 = d == 2 ? nullptr : LocMatrices[5];

  double ** MatrixA31 = d == 2 ? nullptr : LocMatrices[6];
  double ** MatrixA32 = d == 2 ? nullptr : LocMatrices[7];
  double ** MatrixA33 = d == 2 ? nullptr : LocMatrices[8];

  int N_U = N_BaseFuncts[0];

  // original values: u,p,u_x,u_y[,u_z]
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = d == 2 ? nullptr : OrigValues[4];

  double nu = coeff[0]; // = 1/reynolds_number
  double sigma = d == 2 ? coeff[4] : coeff[5];
  ///@todo generalize characteristic lenght. One could also give the option of
  // using L_0 = h

  double delta = compute_GradDiv_delta(stab, nu, sigma, characteristic_length);
  //(nu + sigma*characteristic_length*characteristic_length) * stab;

  for(int i = 0; i < N_U; i++)
  {
    double v_x =  u_x[i];
    double v_y =  u_y[i];
    double v_z = d == 2 ? 0. :  u_z[i];

    for(int j = 0; j < N_U; j++)
    {
      MatrixA11[i][j] += Mult * delta * v_x * u_x[j];
      MatrixA12[i][j] += Mult * delta * v_x * u_y[j]; // v_y * u_x[j];
      if (d == 3)
	       MatrixA13[i][j] += Mult * delta * v_x * u_z[j]; //v_z * u_x[j];

      MatrixA21[i][j] += Mult * delta * v_y * u_x[j]; //* v_x * u_y[j];
      MatrixA22[i][j] += Mult * delta * v_y * u_y[j];
      if (d == 3)
	       MatrixA23[i][j] += Mult * delta * v_y * u_z[j]; //* v_z * u_y[j];

      if (d == 3) {
        MatrixA31[i][j] += Mult * delta * v_z * u_x[j]; //* v_x * u_z[j];
        MatrixA32[i][j] += Mult * delta * v_z * u_y[j];//* v_y * u_z[j];
        MatrixA33[i][j] += Mult * delta * v_z * u_z[j];
      }
    }
  }
}

template <int d>
void NSGradDiv_RightHandSide(double Mult, const double *coeff, const double *,
                             double, const double **OrigValues,
                             const int *N_BaseFuncts, double ***,
                             double **LocRhs, double stab,
                             double characteristic_length)
{
  double * Rhs_u1 = LocRhs[0];
  double * Rhs_u2 = LocRhs[1];
  double * Rhs_u3 = d == 2 ? nullptr : LocRhs[2];

  int N_U = N_BaseFuncts[0];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = d == 2 ? nullptr : OrigValues[4];

  double nu = coeff[0]; // = 1/reynolds_number
  double sigma = d == 2 ? coeff[4] : coeff[5]; // = 1/permeability
  double g = d == 2 ? coeff[3] : coeff[4];
  double delta = compute_GradDiv_delta(stab, nu, sigma, characteristic_length);

  for (int i = 0; i < N_U; i++)
  {
    double v_x = u_x[i];
    double v_y = u_y[i];
    double v_z = d == 2 ? 0. : u_z[i];
    Rhs_u1[i] += Mult * delta * g * v_x;
    Rhs_u2[i] += Mult * delta * g * v_y;
    if (d == 3)
      Rhs_u3[i] += Mult * delta * g * v_z;
  }

}

double compute_PSPG_delta(double delta0, double hK, double nu)
{
  return delta0 * hK * hK / nu;
}

template <int d>
void NSPSPG(double Mult, const double *coeff, const double *, double hK,
            const double **OrigValues, const int *N_BaseFuncts,
            double ***LocMatrices, double **, double delta0)
{
  double ** MatrixB1 = LocMatrices[d == 2 ? 5 : 10];
  double ** MatrixB2 = LocMatrices[d == 2 ? 6 : 11];
  double ** MatrixB3 = d == 2 ? nullptr : LocMatrices[12];
  double ** MatrixC = LocMatrices[d == 2 ? 4 : 9];
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  const double * p_x = OrigValues[2+d];
  const double * p_y = OrigValues[3+d];
  const double * p_z = d == 2 ? nullptr : OrigValues[4+d];
  const double * u_xx = OrigValues[2+2*d];
  const double * u_yy = OrigValues[2+2*d+d];
  const double * u_zz = d == 2 ? nullptr : OrigValues[13];
  double nu = coeff[0]; // = 1/reynolds_number
  double delta = compute_PSPG_delta(delta0, hK, nu);
  for(int i = 0; i < N_P; i++)
  {
    double test_x = delta * p_x[i];
    double test_y = delta * p_y[i];
    double test_z = d == 2 ? 0. : delta * p_z[i];
    for(int j = 0; j < N_U; j++)
    {
      double laplace = nu * (u_xx[j] + u_yy[j] + (d == 2 ? 0. : u_zz[j]));
      MatrixB1[i][j] += Mult*test_x*laplace;
      MatrixB2[i][j] += Mult*test_y*laplace;
      if(d == 3)
        MatrixB3[i][j] += Mult*test_z*laplace;
    }
    for(int j = 0; j < N_P; j++)
    {
      double ansatz_x = p_x[j];
      double ansatz_y = p_y[j];
      double ansatz_z = d == 2 ? 0. : p_z[j];
      MatrixC[i][j] -= Mult * (test_x*ansatz_x + test_y*ansatz_y
                              +test_z*ansatz_z);
    }
  }
}

template <int d>
void NSPSPG_RightHandSide(double Mult, const double *coeff, const double *,
                          double hK, const double **OrigValues,
                          const int *N_BaseFuncts, double ***, double **LocRhs,
                          double delta0)
{
  double * Rhs_div = LocRhs[d];
  int N_P = N_BaseFuncts[1];
  const double * p_x = OrigValues[2+d];
  const double * p_y = OrigValues[3+d];
  const double * p_z = d == 2 ? nullptr : OrigValues[4+d];
  double nu = coeff[0]; // = 1/reynolds_number
  double f1 = coeff[1];
  double f2 = coeff[2];
  double f3 = d == 2 ? 0. : coeff[3];
  double delta = compute_PSPG_delta(delta0, hK, nu);
  for(int i = 0; i < N_P; i++)
  {
    double test_x = p_x[i];
    double test_y = p_y[i];
    double test_z = d == 2 ? 0. : p_z[i];
    Rhs_div[i] -= Mult * delta * (test_x*f1 + test_y*f2 + test_z*f3);
  }
}


// general version (for Navier-Stokes and Brinkman problems)
double compute_GLS_delta(double delta0, double hK, double nu,
                         double sigma, double characteristic_length)
{
  return delta0 * hK * hK / (nu + sigma * characteristic_length * characteristic_length);
}

// common implementation for symmetric and non-symmetric GLS:
template<int d>
void NS_GLS(double Mult, const double *coeff, const double *, double hK,
            const double **OrigValues, const int *N_BaseFuncts,
            double ***LocMatrices, double **, int sign, double delta0,
            double characteristic_length = 0.1)
{
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA22 = LocMatrices[d+1];
  double ** MatrixA33 = d == 2 ? nullptr : LocMatrices[8];
  double ** MatrixB1 = LocMatrices[d == 2 ? 5 : 10];
  double ** MatrixB2 = LocMatrices[d == 2 ? 6 : 11];
  double ** MatrixB3 = d == 2 ? nullptr : LocMatrices[12];
  double ** MatrixB1T = LocMatrices[d == 2 ? 7 : 13];
  double ** MatrixB2T = LocMatrices[d == 2 ? 8 : 14];
  double ** MatrixB3T = d == 2 ? nullptr : LocMatrices[15];
  double ** MatrixC = LocMatrices[d == 2 ? 4 : 9];
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  const double * p_x = OrigValues[2+d];
  const double * p_y = OrigValues[3+d];
  const double * p_z = d == 2 ? nullptr : OrigValues[4+d];

  const double * u = OrigValues[0];

  const double * u_xx = OrigValues[2+2*d];
  const double * u_yy = OrigValues[2+2*d+d];
  const double * u_zz = d == 2 ? nullptr : OrigValues[13];
  double nu = coeff[0]; // = 1/reynolds_number
  double sigma = d == 2 ? coeff[4] : coeff[5]; // =1/permeability

  double delta = compute_GLS_delta(delta0, hK, nu, sigma, characteristic_length);

  //cout << "NS_GLS():   delta0: " << delta0 << ", sigma: " << sigma << ", nu: " << nu << ", delta: "<< delta  << endl;


  for (int i = 0; i < N_U; i++)
  {
    double laplace_v = u_xx[i] + u_yy[i] + (d == 2 ? 0. : u_zz[i]);
    // resistance term for Brinkman-type problems
    double resistance_v = u[i];

    // multiply by -1 for non-sym GLS
    //laplace_v *= delta * sign *
    double residue_v = nu * laplace_v - sigma * resistance_v;

    for (int j = 0; j < N_U; j++)
    {
      double laplace_u = u_xx[j] + u_yy[j] + (d == 2 ? 0. : u_zz[j]);
      double resistance_u =  u[j];
      double residue_u = - nu * laplace_u + sigma * resistance_u;

      MatrixA11[i][j] += Mult * delta * sign * residue_v * residue_u;
      MatrixA22[i][j] += Mult * delta * sign * residue_v * residue_u;
      if (d == 3)
        MatrixA33[i][j] += Mult * delta * sign * residue_v * residue_u;
    }

    for (int j = 0; j < N_P; j++)
    {
      double ansatz_x = p_x[j];
      double ansatz_y = p_y[j];
      double ansatz_z = d == 2 ? 0. : p_z[j];
      MatrixB1T[i][j] += Mult * delta * sign * residue_v * ansatz_x;
      MatrixB2T[i][j] += Mult * delta * sign * residue_v * ansatz_y;
      if(d == 3)
        MatrixB3T[i][j] += Mult * delta * sign * residue_v * ansatz_z;
    }
  }

  for (int i = 0; i < N_P; i++)
  {
    double q_x = p_x[i];
    double q_y = p_y[i];
    double q_z = d == 2 ? 0. : p_z[i];

    for(int j = 0; j < N_U; j++)
    {
      double laplace_u = u_xx[j] + u_yy[j] + (d == 2 ? 0. : u_zz[j]);
      double resistance_u = u[j];
      double residue_u = - nu * laplace_u + sigma * resistance_u;

      MatrixB1[i][j] +=  Mult * delta * sign * (-q_x) * residue_u;
      MatrixB2[i][j] +=  Mult * delta * sign * (-q_y) * residue_u;
      if(d == 3)
        MatrixB3[i][j] += sign * Mult * delta * (-q_z) * residue_u;
    }
    for(int j = 0; j < N_P; j++)
    {
      double ansatz_x = p_x[j];
      double ansatz_y = p_y[j];
      double ansatz_z = d == 2 ? 0. : p_z[j];
      // positive if sign=-1 (non-sym GLS)
      MatrixC[i][j] -=  Mult * delta * sign * (q_x * ansatz_x + q_y * ansatz_y
                              + q_z * ansatz_z);
    }
  }
}

template<int d>
void NS_GLS_RightHandSide(double Mult, const double *coeff, const double *,
                          double hK, const double **OrigValues,
                          const int *N_BaseFuncts, double ***, double **LocRhs,
                          int sign, double delta0,
                          double characteristic_length = 0.1)
{
  double * Rhs1 = LocRhs[0];
  double * Rhs2 = LocRhs[1];
  double * Rhs3 = d == 2 ? nullptr : LocRhs[2];
  double * Rhs_div = LocRhs[d];
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  const double * p_x = OrigValues[2+d];
  const double * p_y = OrigValues[3+d];
  const double * p_z = d == 2 ? nullptr : OrigValues[4+d];

  const double * u = OrigValues[0];

  const double * u_xx = OrigValues[2+2*d];
  const double * u_yy = OrigValues[2+2*d+d];
  const double * u_zz = d == 2 ? nullptr : OrigValues[13];
  double nu = coeff[0]; // = 1/reynolds_number
  double f1 = coeff[1];
  double f2 = coeff[2];
  double f3 = d == 2 ? 0. : coeff[3];
  double sigma = d == 2 ? coeff[4] : coeff[5]; // = 1/permeability

  double delta = compute_GLS_delta(delta0, hK, nu, sigma, characteristic_length);

  //cout << "NS_GLS_RightHandSide():    delta0: " << delta0 << ", sigma: " << sigma << ", nu: " << nu << ", delta: "<< delta << endl;

  for (int i = 0; i < N_U; i++)
  {
    double laplace_v = u_xx[i] + u_yy[i] + (d == 2 ? 0. : u_zz[i]);
    double resistance_v = u[i];

    double residue_v = nu * laplace_v - sigma * resistance_v;
    Rhs1[i] += Mult * delta * sign * residue_v * f1;
    Rhs2[i] += Mult * delta * sign * residue_v * f2;
    if(d == 3)
      Rhs3[i] += Mult * delta * sign * residue_v * f3;
  }
  for(int i = 0; i < N_P; i++)
  {
    double q_x = p_x[i];
    double q_y = p_y[i];
    double q_z = d == 2 ? 0. : p_z[i];
    Rhs_div[i] -= Mult * delta * sign * (q_x * f1 + q_y * f2 + q_z * f3);
  }
}



template <int d>
void NSsymmGLS(double Mult, const double *coeff, const double *param, double hK,
               const double **OrigValues, const int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs, double delta0)
{
  NS_GLS<d>(Mult, coeff, param, hK, OrigValues, N_BaseFuncts, LocMatrices,
            LocRhs, 1, delta0);
}

template <int d>
void NSsymmGLS_RightHandSide(double Mult, const double *coeff,
                             const double *param, double hK,
                             const double **OrigValues, const int *N_BaseFuncts,
                             double ***LocMatrices, double **LocRhs,
                             double delta0)
{
  NS_GLS_RightHandSide<d>(Mult, coeff, param, hK, OrigValues, N_BaseFuncts,
                          LocMatrices, LocRhs, 1, delta0);
}

template <int d>
void NSnonsymmGLS(double Mult, const double *coeff, const double *param,
                  double hK, const double **OrigValues, const int *N_BaseFuncts,
                  double ***LocMatrices, double **LocRhs, double delta0,
                  double characteristic_length)
{
  NS_GLS<d>(Mult, coeff, param, hK, OrigValues, N_BaseFuncts, LocMatrices,
            LocRhs, -1, delta0, characteristic_length);
}

template <int d>
void NSnonsymmGLS_RightHandSide(double Mult, const double *coeff,
                                const double *param, double hK,
                                const double **OrigValues,
                                const int *N_BaseFuncts, double ***LocMatrices,
                                double **LocRhs, double delta0,
                                double characteristic_length)
{
  NS_GLS_RightHandSide<d>(Mult, coeff, param, hK, OrigValues, N_BaseFuncts,
                          LocMatrices, LocRhs, -1, delta0, characteristic_length);
}

template <int d>
void NS_BrezziPitkaeranta(double Mult, const double *coeff, const double *,
                          double hK, const double **OrigValues,
                          const int *N_BaseFuncts, double ***LocMatrices,
                          double **, double delta0)
{
  double ** MatrixC = LocMatrices[d == 2 ? 4 : 9];
  int N_P = N_BaseFuncts[1];
  const double * p_x = OrigValues[2+d];
  const double * p_y = OrigValues[3+d];
  const double * p_z = d == 2 ? nullptr : OrigValues[4+d];
  double nu = coeff[0]; // = 1/reynolds_number
  double delta = compute_PSPG_delta(delta0, hK, nu);
  for(int i = 0; i < N_P; i++)
  {
    double test_x = delta * p_x[i];
    double test_y = delta * p_y[i];
    double test_z = d == 2 ? 0. : delta * p_z[i];
    for(int j = 0; j < N_P; j++)
    {
      double ansatz_x = p_x[j];
      double ansatz_y = p_y[j];
      double ansatz_z = d == 2 ? 0. : p_z[j];
      MatrixC[i][j] -= Mult * (test_x*ansatz_x + test_y*ansatz_y
                              +test_z*ansatz_z);
    }
  }
}

template <int d>
void OseenSingle(double Mult, const double *coeff, const double *, double,
                       const double**OrigValues, const int *N_BaseFuncts,
                       double ***LocMatrices, double **)
{
  if(d==3)
    ErrThrow("3D is not tested so far");
  double ** MatrixA11 = LocMatrices[0];
  double ** MatrixA22 = LocMatrices[d+1];
  int N_U = N_BaseFuncts[0];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];

  const double *u = OrigValues[0];
  const double u1 = coeff[5];
  const double u2 = coeff[6];
  double sigma = coeff[7];

  for(int i = 0; i < N_U; i++)
  {
    double test   = u[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz   = u[j];
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];

      double val = Mult * (u1 * ansatz_x + u2*ansatz_y )*test;
      val += Mult * sigma * test * ansatz;
      MatrixA11[i][j] += val;
      MatrixA22[i][j] += val;
    }
  }
}

template <int d>
void OseenSingleSUPG(double Mult, const double *coeff, const double *, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **, double delta0, double /*delta1*/)
{
  if(d==3)
    ErrThrow("Only 2D case is implemented");

  double ** MatrixA11 = LocMatrices[0];
  // double ** MatrixA12 = LocMatrices[1];
  // double ** MatrixA21 = LocMatrices[2];
  double ** MatrixA22 = LocMatrices[d+1];
  int N_U = N_BaseFuncts[0];

  const double *u = OrigValues[0];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_xx = OrigValues[2+2*d];
  const double * u_yy = OrigValues[2+2*d+1];

  double nu = coeff[0]; // = 1/reynolds_number
  const double u1 = coeff[4];
  const double u2 = coeff[5];
  double sigma = coeff[6];

  double taum = delta0*hK*hK;

  for(int i = 0; i < N_U; i++)
  {
    double test   = u[i];
    double test_x = u_x[i]; // vx
    double test_y = u_y[i]; // vy

    double ugradv = taum * (u1*test_x + u2*test_y);
    for(int j = 0; j < N_U; j++)
    {
      double ansatz   = u[j];
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];

      double val = Mult * nu* ( test_x * ansatz_x + test_y * ansatz_y);// standard
      val += Mult *(u1 * ansatz_x + u2*ansatz_y )*test; // convective part
      val += Mult * sigma * test * ansatz; //reaction
      val += Mult * (-nu * (u_xx[j] + u_yy[j]) + (u1*ansatz_x+u2 *ansatz_y) + sigma *ansatz ) * ugradv;

      MatrixA11[i][j] += val;
      // MatrixA12[i][j] += Mult * delta1 * test_x * ansatz_y;
      //MatrixA21[i][j] += Mult * delta1 * test_y * ansatz_x;
      MatrixA22[i][j] += val;
    }
  }
}
template <int d>
void OseenGradBlockSUPG(double Mult, const double* coeff, const double*,
              double hK, const double ** OrigValues, const int* N_BaseFuncts,
              double *** LocMatrices, double **, double delta0)
{
  if(d==3)
    ErrThrow("3D case is not implemented");

  double ** MatrixB1T = LocMatrices[d == 2 ? 7 : 13];
  double ** MatrixB2T = LocMatrices[d == 2 ? 8 : 14];
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  const double * p = OrigValues[1];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * p_x = OrigValues[2+d];
  const double * p_y = OrigValues[3+d];

  const double u1 = coeff[4];
  const double u2 = coeff[5];
  double taum = delta0*hK*hK;

  for(int i = 0; i < N_U; i++)
  {
    double test_x = u_x[i];
    double test_y = u_y[i];

    double ugradv = taum * (u1 * test_x + u2*test_y );

    for(int j = 0; j < N_P; j++)
    {
      double ansatz = p[j];

      MatrixB1T[i][j] -= Mult * ansatz * test_x;
      MatrixB1T[i][j] += Mult * ugradv * p_x[j];
      MatrixB2T[i][j] -= Mult * ansatz * test_y;
      MatrixB2T[i][j] += Mult * ugradv * p_y[j];
    }
  }
}

template <int d>
void OseenRhsSUPG(double Mult, const double *coeff, const double *, double hK,
                  const double**OrigValues, const int *N_BaseFuncts, double ***,
                  double **LocRhs, double delta0)
{
  double * Rhs1 = LocRhs[0];
  double * Rhs2 = LocRhs[1];
  double * Rhs_div = LocRhs[2];
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  const double * u = OrigValues[0];
  const double * p = OrigValues[1];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];

  double f1 = coeff[1];
  double f2 = coeff[2];
  double g = coeff[d+1]; // divergence

  const double u1 = coeff[4];
  const double u2 = coeff[5];

  double taum = delta0*hK*hK;

  for(int i = 0; i < N_U; i++)
  {
    double test = u[i];
    double test_x = u_x[i];
    double test_y = u_y[i];
    double ugradv = taum * (u1 * test_x + u2 * test_y);
    Rhs1[i] += Mult * (test + ugradv) * f1;
    Rhs2[i] += Mult * (test + ugradv) * f2;
  }
  for(int i = 0; i < N_P; i++)
  {
    double test = p[i];
    Rhs_div[i] -= Mult * test * g;
  }
}

template <int d>
void NSLaplaceDeformationSUPG(double Mult, const double* coeff, const double* param,
  double hK, const double ** OrigValues, const int* N_BaseFuncts,
  double *** LocMatrices, double **, double delta0, double delta1)
{
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double ** MatrixA13 = d == 2 ? nullptr : LocMatrices[2];
  double **MatrixA21 = LocMatrices[d];
  double **MatrixA22 = LocMatrices[d+1];
  double **MatrixA23 = d == 2 ? nullptr : LocMatrices[d+2];
  double ** MatrixA31 = d == 2 ? nullptr : LocMatrices[6];
  double ** MatrixA32 = d == 2 ? nullptr : LocMatrices[7];
  double ** MatrixA33 = d == 2 ? nullptr : LocMatrices[8];
  int N_U = N_BaseFuncts[0];
  const double * u = OrigValues[0];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = d == 2 ? nullptr : OrigValues[4];
  const double * u_xx = OrigValues[2+2*d];
  const double * u_yy = OrigValues[2+2*d+1];
  //NOTE: please do not use 2nd derivatives for 3D, because it doesn't
  //supported and leads to Sygmentation fault
  double nu = coeff[0];

  double taum = delta0*hK*hK;
  double tauc = delta1;
  const double u1 = param[0];
  const double u2 = param[1];
  const double u3 = d==2 ? 0. : param[2];

  for(int i = 0; i < N_U; i++)
  {
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0. : u_z[i];
    double test = u[i];
    double ugradv = taum*(u1*test_x + u2*test_y + u3*test_z);
    for(int j = 0; j < N_U; j++)
    {
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      double laplac = 0.;
      if(d==2)
        laplac = -taum * nu* (u_xx[j] + u_yy[j]);
      else
        laplac = 0.;
      double val  = 2*nu*(test_x*ansatz_x + 0.5*test_y*ansatz_y
                                     +0.5*test_z*ansatz_z);
      double ugradu = (u1*ansatz_x + u2*ansatz_y + u3*ansatz_z);

      val += ugradu * test;

      // supg part
      val += laplac * ugradv;
      val += ugradu * ugradv;
      // grad-div contribution
      val += tauc * test_x * ansatz_x;

      MatrixA11[i][j] += Mult * val;

      val = nu*(test_y*ansatz_x);
      val += tauc * test_x * ansatz_y;
      MatrixA12[i][j] += Mult * val ;
      if(d == 3)
      {
        val = nu*(test_z*ansatz_x);
        val += tauc * test_x * ansatz_z;
        MatrixA13[i][j] += Mult * val;
      }
      val = nu*(test_x*ansatz_y);
      val += tauc * test_y * ansatz_x;
      MatrixA21[i][j] += Mult * val;

      val = 2*nu*(0.5*test_x*ansatz_x+test_y*ansatz_y +0.5*test_z*ansatz_z);
      // nonlinear part
      val += ugradu * test;
      // supg contribution
      val += laplac * ugradv;
      val += ugradu * ugradv;
      val += tauc * test_y * ansatz_y;

      MatrixA22[i][j] += Mult * val;

      if(d == 3)
      {
        val = nu*(test_z*ansatz_y);
        val += tauc * test_y * ansatz_z;
        MatrixA23[i][j] += Mult * val;

        val = nu*(test_x*ansatz_z);
        val += tauc * test_z * ansatz_x;
        MatrixA31[i][j] += Mult * val;

        val = nu*(test_y*ansatz_z);
        val += tauc * test_z * ansatz_y;
        MatrixA32[i][j] += Mult * val;

        val = 2*nu*(0.5*test_x*ansatz_x+0.5*test_y*ansatz_y +test_z*ansatz_z);
        // nonlinear part
        val += ugradu * test;
        // supg contribution
        val += laplac * ugradv;
        val += ugradu * ugradv;
        // grad -div
        val += tauc * test_z * ansatz_z;

        MatrixA33[i][j] += Mult * val;
      }
    }
  }
}

template <int d>
void NSGradientBlocksSUPG(double Mult, const double *, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **, double delta0)
{
  double ** MatrixB1T = LocMatrices[d == 2 ? 7 : 13];
  double ** MatrixB2T = LocMatrices[d == 2 ? 8 : 14];
  double ** MatrixB3T = d == 2 ? nullptr : LocMatrices[15];
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  const double * p = OrigValues[1];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = d == 2 ? nullptr : OrigValues[4];
  const double * p_x = OrigValues[2+d];
  const double * p_y = OrigValues[3+d];
  const double * p_z = d == 2 ? nullptr : OrigValues[4+d];

  const double u1 = param[0];
  const double  u2 = param[1];
  const double  u3 = d==2 ? 0. : param[2];

  double taum = delta0*hK*hK;

  for(int i = 0; i < N_U; i++)
  {
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0. : u_z[i];

    double ugradv = taum * (u1 * test_x + u2*test_y + u3*test_z);

    for(int j = 0; j < N_P; j++)
    {
      double ansatz = p[j];
      double ansatz_x = p_x[j];
      double ansatz_y = p_y[j];
      double ansatz_z = d==2 ? 0. : p_z[j];

      MatrixB1T[i][j] += Mult * (-ansatz * test_x + ugradv * ansatz_x) ;
      MatrixB2T[i][j] += Mult * (-ansatz * test_y + ugradv * ansatz_y) ;
      if(d == 3)
        MatrixB3T[i][j] += Mult * (-ansatz * test_z + ugradv * ansatz_z);
    }
  }
}

template <int d>
void NSRightHandSideSUPG(double Mult, const double *coeff, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***, double **LocRhs, double delta0)
{
  double * Rhs1 = LocRhs[0];
  double * Rhs2 = LocRhs[1];
  double * Rhs3 = d == 2 ? nullptr : LocRhs[2];
  double * Rhs_div = LocRhs[d];
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  const double * u = OrigValues[0];
  const double * p = OrigValues[1];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = d == 2 ? nullptr : OrigValues[4];

  double f1 = coeff[1];
  double f2 = coeff[2];
  double f3 = d == 2 ? 0. : coeff[3];
  double g = coeff[d+1]; // divergence


  const double u1 = param[0];
  const double u2 = param[1];
  const double u3 = d==2 ? 0. : param[2];
  double taum = delta0*hK*hK;

  for(int i = 0; i < N_U; i++)
  {
    double test = u[i];
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0. : u_z[i];
    double ugradv = taum * (u1 * test_x + u2 * test_y + u3 * test_z);
    Rhs1[i] += Mult * (test + ugradv) * f1;
    Rhs2[i] += Mult * (test + ugradv) * f2;
    if(d == 3)
      Rhs3[i] += Mult * (test + ugradv )*f3;
  }
  for(int i = 0; i < N_P; i++)
  {
    double test = p[i];
    Rhs_div[i] -= Mult * test * g;
  }
}

template <int d>
void NSParamsVelocity(const double *in, double *out)
{
  out[0] = in[d]; // u1old
  out[1] = in[d + 1]; // u2old

  if (d == 3)
  {
    out[2] = in[d + 2]; // u3old
  }
}

template <int d>
void NSParamsVelocityDerivatives(const double *in, double *out)
{
  out[0] = in[d]; // u1old
  out[1] = in[d + 1]; // u1old_x
  out[2] = in[d + 2]; // u1old_y

  if (d == 3)
  {
    out[3] = in[d + 3]; // u1old_z
  }

  out[d + 1] = in[2 * d + 1]; // u2old
  out[d + 2] = in[2 * d + 2]; // u2old_x
  out[d + 3] = in[2 * d + 3]; // u2old_y

  if (d == 3)
  {
    out[d + 4] = in[2 * d + 4]; // u2old_z
  }

  if (d == 3)
  {
    out[2 * d + 2] = in[3 * d + 2]; // u3old
    out[2 * d + 3] = in[3 * d + 3]; // u3old_x
    out[2 * d + 4] = in[3 * d + 4]; // u3old_y
    out[2 * d + 5] = in[3 * d + 5]; // u3old_z
  }
}

template <int d>
void NSParamsVelocityGradient(const double *in, double *out)
{
  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < d; j++)
    {
      out[i * d + j] = in[d + (d + 1) * i + 1 + j]; // d_j u_i
    }
  }
}

#ifdef __2D__
void NSParamVelocityGradients(const double *in, double *out)
{
  out[0] = in[2];   // u1
  out[1] = in[3];   // u2
  out[2] = in[4];   // u1x 
  out[3] = in[5];   // u2x
  out[4] = in[6];   // u1y 
  out[5] = in[7];   // u2y 
  out[6] = in[8];   // u1xy
  out[7] = in[9];   // u1xx
  out[8] = in[10];  // u1yy
  out[9] = in[11];  // u2xy
  out[10] = in[12]; // u2xx
  out[11] = in[13]; // u2yy
}
#endif

// explicit instatiations below:
#ifdef __2D__
template void NSResistanceMassMatrixSingle<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSResistanceMassMatrix<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSLaplaceGradGradSingle<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSLaplaceGradGrad<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSLaplaceDeformation<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSDivergenceBlocks<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, int sign);
template void NSGradientBlocks<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSRightHandSide<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, int sign);
template void NSLaplaceGradGrad_dg<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSNonlinearTerm_convective_Single<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSNonlinearTerm_convective<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSNonlinearTerm_convective_dg<2>(
double Mult, const double *, const double *param, double, const double **OrigValues,
 const int *N_BaseFuncts, double ***LocMatrices, double **);
template void NSNonlinearTerm_skew_symmetric_Single<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSNonlinearTerm_skew_symmetric<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSNonlinearTerm_rotational<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSNonlinearTerm_emac<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);

template void NSGradDiv<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0, double characteristic_length);
template void NSGradDiv_RightHandSide<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0, double characteristic_length);
template void NSPSPG<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0);
template void NSPSPG_RightHandSide<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0);
template void NSsymmGLS<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0);
template void NSsymmGLS_RightHandSide<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0);
template void NSnonsymmGLS<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0, double characteristic_length);
template void NSnonsymmGLS_RightHandSide<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0, double characteristic_length);
template void NS_BrezziPitkaeranta<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0);
template void NSNonlinearTerm_divergence<2>(
  double Mult, const double *, const double *param, double,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **);
// not yet available
// template void NSCoriolis<2>(
//   double Mult, const double *coeff, const double *param, double hK,
//   const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
//   double **LocRhs);
template void NSRightHandSideExplicitNL<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void OseenSingle<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **);
template void OseenSingleSUPG<2>(
  double Mult, const double *coeff, const double *, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **, double delta0, double delta1);
template void OseenGradBlockSUPG<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **,double delta0);
template void OseenRhsSUPG<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0);
template void NSLaplaceDeformationSUPG<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **, double delta0, double delta1);
template void NSGradientBlocksSUPG<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0);
template void NSRightHandSideSUPG<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0);

template void NSParamsVelocity<2>(const double *in, double *out);
template void NSParamsVelocityDerivatives<2>(const double *in, double *out);
template void NSParamsVelocityGradient<2>(const double *in, double *out);

#endif // 2D
#ifdef __3D__
template void NSResistanceMassMatrixSingle<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSResistanceMassMatrix<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSLaplaceGradGradSingle<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSLaplaceGradGrad<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSLaplaceDeformation<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSDivergenceBlocks<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, int sign);
template void NSGradientBlocks<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSRightHandSide<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, int sign);
template void NSLaplaceGradGrad_dg<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSNonlinearTerm_convective_Single<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSNonlinearTerm_convective<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSNonlinearTerm_convective_dg<3>(
double Mult, const double *, const double *param, double, const double **OrigValues,
 const int *N_BaseFuncts, double ***LocMatrices, double **);
template void NSNonlinearTerm_skew_symmetric_Single<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSNonlinearTerm_skew_symmetric<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSNonlinearTerm_rotational<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSNonlinearTerm_emac<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSNonlinearTerm_divergence<3>(
  double Mult, const double *, const double *param, double,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **);
template void NSCoriolis<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSRightHandSideExplicitNL<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSParamsVelocity<3>(const double *in, double *out);
template void NSParamsVelocityDerivatives<3>(const double *in, double *out);
template void NSParamsVelocityGradient<3>(const double *in, double *out);
template void NSGradDiv<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0, double characteristic_length);
template void NSGradDiv_RightHandSide<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0, double characteristic_length);
template void NSPSPG<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0);
template void NSPSPG_RightHandSide<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0);
template void NSsymmGLS<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0);
template void NSsymmGLS_RightHandSide<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0);
template void NSnonsymmGLS<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0, double characteristic_length);
template void NSnonsymmGLS_RightHandSide<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0, double characteristic_length);
template void NS_BrezziPitkaeranta<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0);
template void OseenSingle<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **);
template void OseenSingleSUPG<3>(
  double Mult, const double *coeff, const double *, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **, double delta0, double delta1);
template void OseenGradBlockSUPG<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **,double delta0);
template void OseenRhsSUPG<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0);
template void NSLaplaceDeformationSUPG<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **, double delta0, double delta1);
template void NSGradientBlocksSUPG<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0);
template void NSRightHandSideSUPG<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, double delta0);
#endif // 3D
