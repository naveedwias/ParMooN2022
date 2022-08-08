#include "DarcyMixed.h"
#include "MooNMD_Io.h"


// ======================================================================
// (DarcyType 1)
// Standard Galerkin with Raviart-Thomas (RT) or Brezzi-Douglas-Marini (BDM)
// elements
// ======================================================================
template <int d>
void BilinearAssembleDarcyGalerkin(double Mult, const double *coeff,
                                   const double *, double,
                                   const double **OrigValues,
                                   const int *N_BaseFuncts,
                                   double ***LocMatrices, double **LocRhs)
{
  double val;
  double ansatz, ansatz_x, ansatz_y, ansatz_z;
  double test, test_x, test_y, test_z;
  double test_x_100, test_y_010, test_z_001;
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
    
//     Output::print("  test function ", test_x, " ", test_y, " ", test_z);
    
    Rhs0[i] += Mult*(f1*test_x + f2*test_y + f3*test_z);
    
    for(int j=0;j<N_U;j++)
    {
      ansatz_x = Orig0[j];
      ansatz_y = Orig0[N_U+j];
      ansatz_z = d == 2 ? 0. : Orig0[2*N_U+j];

      // A: u_x v_x + u_y v_y
      val  = c0*(test_x*ansatz_x + test_y*ansatz_y + test_z*ansatz_z);
      MatrixA[i][j] += Mult * val;
    }
    // B1, B2:
    test_x_100 = Orig2[i];
    test_y_010 = Orig3[N_U+i];
    test_z_001 = d == 2 ? 0. : Orig4[2*N_U+i];
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

#ifdef __3D__
template void BilinearAssembleDarcyGalerkin<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
#else
template void BilinearAssembleDarcyGalerkin<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
#endif


