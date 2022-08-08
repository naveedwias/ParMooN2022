#include "ROM_local_assembling_routines.h"

/** ************************************************************************ */
template<int d>
void TCDDiffOnlyScale(double Mult, const double *, const double *, double hK,
                      const double **OrigValues, const int *N_BaseFuncts,
                      double ***LocMatrices, double **)
{
  // matrix for pressure/pressure term
  double **MatrixC = LocMatrices[1];
  int N_P = N_BaseFuncts[0];
  const double * p_x = OrigValues[1];
  const double * p_y = OrigValues[2];
  const double * p_z = d == 2 ? nullptr : OrigValues[3];

  for(int i=0 ; i<N_P ; i++)
  {
    double test_x = p_x[i];
    double test_y = p_y[i];
    double test_z = d == 2 ? 0. : p_z[i];
    for(int j=0 ; j<N_P ; j++)
    {
      double ansatz_x = p_x[j];
      double ansatz_y = p_y[j];
      double ansatz_z = d == 2 ? 0. : p_z[j];
      double val = Mult * (test_x*ansatz_x + test_y*ansatz_y + test_z*ansatz_z);
      MatrixC[i][j] += val * hK;
    }
  }
}


/** ************************************************************************ */
template<int d>
void NSPROM_Divergence(double Mult, const double *, const double *, double hK,
                       const double **OrigValues, const int *N_BaseFuncts,
                       double ***LocMatrices, double **)
{
  double ** MatrixB1 = LocMatrices[d == 2 ? 5 : 10];
  double ** MatrixB2 = LocMatrices[d == 2 ? 6 : 11];
  double ** MatrixB3 = d == 2 ? nullptr : LocMatrices[12];
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  const double * u = OrigValues[0];
  const double * p_x = OrigValues[2+d];
  const double * p_y = OrigValues[3+d];
  const double * p_z = d == 2 ? nullptr : OrigValues[4+d];

  for(int i=0 ; i<N_P ; i++)
  {
    double test_x = p_x[i];
    double test_y = p_y[i];
    double test_z = d == 2 ? 0. : p_z[i];
    for(int j=0 ; j<N_U ; j++)
    {
      double ansatz = u[j];
      MatrixB1[i][j] += -Mult * hK * ansatz * test_x;
      MatrixB2[i][j] += -Mult * hK * ansatz * test_y;
      if(d == 3)
      {
        MatrixB3[i][j] += -Mult * hK * ansatz * test_z;
      }
    }
  }
}


/** ************************************************************************ */
template<int d>
void NSPROM_NonlinearTerm(double Mult, const double *, const double *param,
                          double hK, const double **OrigValues,
                          const int *N_BaseFuncts, double ***LocMatrices,
                          double **)
{
  double ** MatrixB1 = LocMatrices[d == 2 ? 5 : 10];
  double ** MatrixB2 = LocMatrices[d == 2 ? 6 : 11];
  double ** MatrixB3 = d == 2 ? nullptr : LocMatrices[12];
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = d == 2 ? nullptr : OrigValues[4];
  const double * p_x = OrigValues[2+d];
  const double * p_y = OrigValues[3+d];
  const double * p_z = d == 2 ? nullptr : OrigValues[4+d];
  const double u1 = param[0];
  const double u2 = param[1];
  const double u3 = d == 2 ? 0. : param[2];

  for(int i=0;i<N_P;i++)
  {
    double test_x = p_x[i];
    double test_y  = p_y[i];
    double test_z = d == 2 ? 0. : p_z[i];
    for(int j=0;j<N_U;j++)
    {
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      double val = - Mult * hK *(u1 * ansatz_x + u2 * ansatz_y + u3 * ansatz_z);
      MatrixB1[i][j] += val * test_x;
      MatrixB2[i][j] += val * test_y;
      if(d == 3)
      {
        MatrixB3[i][j] += val * test_z;
      }
    }
  }
}


/** ************************************************************************ */
template<int d>
void NSPROM_RightHandSide(double Mult, const double *coeff, const double *,
                          double hK, const double **OrigValues,
                          const int *N_BaseFuncts, double ***,
                          double **LocRhs)
{
  double * Rhs = LocRhs[d];
  int N_P = N_BaseFuncts[1];
  const double * p_x = OrigValues[2+d];
  const double * p_y = OrigValues[3+d];
  const double * p_z = d == 2 ? nullptr : OrigValues[4+d];
  double f1 = coeff[1];
  double f2 = coeff[2];
  double f3 = d == 2 ? 0. : coeff[3];

  for(int i=0 ; i<N_P ; i++)
  {
    double test_x = p_x[i];
    double test_y = p_y[i];
    double test_z = d == 2 ? 0. : p_z[i];
    Rhs[i] += Mult * hK * (test_x*f1 + test_y*f2 + test_z*f3);
  }
}

// explicit instatiations below:
#ifdef __2D__
template void TCDDiffOnlyScale<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSPROM_Divergence<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSPROM_NonlinearTerm<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSPROM_RightHandSide<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
#endif // 2D
#ifdef __3D__
template void TCDDiffOnlyScale<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSPROM_Divergence<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSPROM_NonlinearTerm<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSPROM_RightHandSide<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
#endif // 3D
