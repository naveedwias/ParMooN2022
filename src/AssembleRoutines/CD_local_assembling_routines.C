#include <cmath>

#include "CD_local_assembling_routines.h"
#include "Database.h"
#include "ConvDiff.h"
#include <MooNMD_Io.h>

template <int d>
void TCDStiff(double Mult, const double *coeff, const double *, double ,
              const double**OrigValues, const int *N_BaseFuncts,
              double ***LocMatrices, double **)
{
  double **Matrix = LocMatrices[0];
  int N_ = N_BaseFuncts[0];
  const double * u  =OrigValues[0];
  const double * u_x=OrigValues[1];
  const double * u_y=OrigValues[2];
  const double * u_z = d == 2 ? nullptr : OrigValues[3];
  
  double eps = coeff[0];
  double b1 = coeff[1];
  double b2 = coeff[2];
  double b3 = d == 2 ? 0. : coeff[3];
  double c = coeff[d+1];
  
  for(int i=0; i<N_; i++)
  {
    double test = u[i];
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0. : u_z[i];
    for(int j=0; j<N_; j++)
    {
      double ansatz=u[j];
      double ansatz_x=u_x[j];
      double ansatz_y=u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      double val = (eps*(test_x*ansatz_x + test_y*ansatz_y + test_z*ansatz_z));
      val += (b1*ansatz_x + b2*ansatz_y+b3*ansatz_z)*test;
      val += c*test*ansatz;
      Matrix[i][j] += Mult * val;
    }
  }
}

template <int d>
void TCDStiffWithoutDiffAndLinearSource(double Mult, const double* coeff,
                                        const double*, double , const double** OrigValues,
                                        const int* N_BaseFuncts, double*** LocMatrices,
                                        double**)
{
  double **Matrix = LocMatrices[0];
  int N_ = N_BaseFuncts[0];
  const double * u   = OrigValues[0];
  const double * u_x = OrigValues[1];
  const double * u_y = OrigValues[2];
  const double * u_z = d == 2 ? nullptr : OrigValues[3];
  
  double b1 = coeff[1];
  double b2 = coeff[2];
  double b3 = d == 2 ? 0.0 : coeff[3];
  
  for(int i=0; i<N_; i++)
  {
    double test = u[i];
    
    for(int j=0; j<N_; j++)
    {
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0.0 : u_z[j];
      double val = ( b1 * ansatz_x + b2 * ansatz_y + b3 * ansatz_z ) * test;
      
      Matrix[i][j] += Mult * val;
    }
  }
}

template <int d>
void TCDDiffAndLinearSourceOnly(double Mult, const double* coeff, const double*,
                                double , const double** OrigValues,
                                const int *N_BaseFuncts, double*** LocMatrices,
                                double**, int matrixIndex)
{
  double **Matrix = LocMatrices[matrixIndex];
  int N_ = N_BaseFuncts[0];
  const double * u   = OrigValues[0];
  const double * u_x = OrigValues[1];
  const double * u_y = OrigValues[2];
  const double * u_z = d == 2 ? nullptr : OrigValues[3];
  
  double eps = coeff[0];
  double c = coeff[d+1];
  
  for(int i=0; i<N_; i++)
  {
    double test = u[i];
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0.0 : u_z[i];
    for(int j=0; j<N_; j++)
    {
      double ansatz = u[j];
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0.0 : u_z[j];
      double val = eps*(test_x*ansatz_x + test_y*ansatz_y + test_z*ansatz_z);
      val += c*test*ansatz;
      
      Matrix[i][j] += Mult * val;
    }
  }
}

template <int d>
void TCDStiff_TensorialDiffusionTerm(double Mult, const double *coeff,
                                     const double *, double,
                                     const double**OrigValues,
                                     const int *N_BaseFuncts,
                                     double ***LocMatrices, double **)
{
  double **Matrix = LocMatrices[0];
  int N_ = N_BaseFuncts[0];
  //const double * u  =OrigValues[0];
  const double * u_x=OrigValues[1];
  const double * u_y=OrigValues[2];
  const double * u_z = d == 2 ? nullptr : OrigValues[3];

  double b1 = coeff[1];
  double b2 = coeff[2];
  double b3 = d == 2 ? 0. : coeff[3];
  double a_l = coeff[d+3];

  for (int i = 0; i < N_; i++)
  {
    //double test = u[i];
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0. : u_z[i];

    for (int j = 0; j < N_; j++)
    {
      //double ansatz=u[j];
      double ansatz_x=u_x[j];
      double ansatz_y=u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];

      double val = a_l * (  ( (b1*b1*ansatz_x + b1*b2*ansatz_y + b1*b3*ansatz_z) * test_x )
                          + ( (b2*b1*ansatz_x + b2*b2*ansatz_y + b2*b3*ansatz_z) * test_y )
                          + ( (b3*b1*ansatz_x + b3*b2*ansatz_y + b3*b3*ansatz_z) * test_z )  );
      Matrix[i][j] += Mult * val;
    }
  }
}

template<int d>
void TCDMass(double Mult, const double *, const double *, double,
             const double**OrigValues,  const int *N_BaseFuncts,
             double ***LocMatrices, double **)
{
  double **Matrix = LocMatrices[1];
  int N_ = N_BaseFuncts[0];
  const double * u  =OrigValues[0];
  
  for(int i=0; i<N_; i++)
    for(int j=0; j<N_; j++)
      Matrix[i][j] += Mult * u[i] * u[j];
}

template<int d>
void TCDRhs(double Mult, const double* coeff, const double*, double,
            const double ** OrigValues, const int* N_BaseFuncts,
            double ***, double ** LocRhs)
{
  int N_ = N_BaseFuncts[0];
  const double *u = OrigValues[0];
  double f = coeff[d+2];
  double *Rhs = LocRhs[0];
  for(int i=0; i<N_; i++)
    Rhs[i] += Mult * u[i] * f;
}

template <int d>
void TCDStiffSUPG(double Mult, const double *coeff, const double *, double hK,
                  const double**OrigValues, const int *N_BaseFuncts,
                  double ***LocMatrices, double **)
{
  double **MatrixA = LocMatrices[0];
  int N_ = N_BaseFuncts[0];

  const double * u  =OrigValues[0];
  const double * u_x=OrigValues[1];
  const double * u_y=OrigValues[2];
  const double * u_z = d == 2 ? nullptr : OrigValues[3];
  const double * u_xx = d== 2 ? OrigValues[3] : OrigValues[4];
  const double * u_yy = d== 2 ? OrigValues[5] : OrigValues[7];
  const double * u_zz = d== 2 ? nullptr : OrigValues[9];

  const double eps = coeff[0];
  const double b1 = coeff[1];
  const double b2 = coeff[2];
  const double b3 = d==2 ? 0.0 : coeff[3];
  const double c = coeff[d+1];

  double theta1 = TDatabase::TimeDB->THETA1;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  int sdfem_type = TDatabase::ParamDB->SDFEM_TYPE;
  double delta;
  if(d==2)
  {
    double c00 = theta1 * time_step * eps;
    double c11 = theta1 * time_step * b1;
    double c22 = theta1 * time_step * b2;
    double c33 = 1.0 + theta1 * time_step * c;
    if(sdfem_type==8)
    {
      c33 = theta1 * time_step * c;
    }
    double bb;
    if(fabs(c11) > fabs(c22))
      bb = fabs(c11);
    else
      bb = fabs(c22);
    if(sdfem_type != 2 && sdfem_type != 12)
      delta = Compute_SDFEM_delta<2>(hK, c00, {{c11, c22}}, c33, bb);
    else
      delta = Compute_SDFEM_delta<2>(hK, eps, {{b1, b2}}, c, bb);
  }
  else
  {
    double b_norm = coeff[6];
    delta = Compute_SDFEM_delta<3>(hK, eps, {{b1, b2, b3}}, c, b_norm);
  }
  double val  = theta1 * time_step;
  if((sdfem_type!=2) && ( sdfem_type != 9)  && (sdfem_type != 10)
    && (sdfem_type != 11) && ( sdfem_type != 12))
    delta *= val;

  for(int i=0;i<N_;i++)
  {
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0. : u_z[i];

    double bgradv = delta*(b1*test_x+b2*test_y + b3*test_z);    

    for(int j=0;j<N_;j++)
    {
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      double ansatz   = u[j];
      double ansatz_xx = u_xx[j];
      double ansatz_yy = u_yy[j];
      double ansatz_zz = d==2 ? 0. : u_zz[j];

      val = (b1*ansatz_x+b2*ansatz_y + b3*ansatz_z + c*ansatz)*bgradv;
      val -= eps*(ansatz_xx + ansatz_yy + ansatz_zz) * bgradv;

      MatrixA[i][j] += Mult * val;
    }
  }
}
template <int d>
void TCDMassSUPG(double Mult, const double *coeff, const double *, double hK,
                 const double**OrigValues, const int *N_BaseFuncts,
                 double ***LocMatrices, double **)
{
  double **Matrix = LocMatrices[1];
  int N_ = N_BaseFuncts[0];
  const double * u  = OrigValues[0];
  const double * ux = OrigValues[1];
  const double * uy = OrigValues[2];
  const double * uz = d==2 ? nullptr : OrigValues[3];

  double eps = coeff[0];
  double b1 = coeff[1];
  double b2 = coeff[2];
  double b3 = d == 2 ? 0. : coeff[3];
  double c = coeff[d+1];
  double delta = 0.;
  // stabilization parameter
  int sdfem_type = TDatabase::ParamDB->SDFEM_TYPE;
  double theta1 = TDatabase::TimeDB->THETA1;
  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  if(d==2)
  {
    double c00 = theta1 * tau * eps;
    double c11 = theta1 * tau * b1;
    double c22 = theta1 * tau * b2;
    double c33 = 1.0 + theta1 * tau * c;
    if(sdfem_type==8)
    {
      c33 = theta1 * tau * c;
    }
    double bb;
    if(fabs(c11) > fabs(c22))
      bb = fabs(c11);
    else
      bb = fabs(c22);
    if(sdfem_type != 2 && sdfem_type != 12)
      delta = Compute_SDFEM_delta<2>(hK, c00, {{c11, c22}}, c33, bb);
    else
      delta = Compute_SDFEM_delta<2>(hK, eps, {{b1, b2}}, c, bb);
  }
  else
  {
    double b_norm = coeff[6];
    delta = Compute_SDFEM_delta<3>(hK, eps, {{b1, b2, b3}}, c, b_norm);
  }
  double val  = theta1 * tau;
  if((sdfem_type!=2) && ( sdfem_type != 9)  && (sdfem_type != 10)
    && (sdfem_type != 11) && ( sdfem_type != 12))
    delta *= val;

  for(int i=0; i<N_; i++)
  {
    double test_x = ux[i];
    double test_y = uy[i];
    double test_z = d==2 ? 0.0 : uz[i];
    double bgradv = delta*(b1 * test_x + b2 * test_y + b3 * test_z);
    
    for(int j=0; j<N_; j++)
      Matrix[i][j] += Mult * u[j] * bgradv;
  }
}

template <int d>
void TCDRhsSUPG(double Mult, const double *coeff, const double *, double hK,
                const double**OrigValues, const int *N_BaseFuncts, double ***,
                double **LocRhs)
{
  int N_ = N_BaseFuncts[0];
  
  const double *u_x = OrigValues[1];
  const double *u_y = OrigValues[2];
  const double *u_z = d==2 ? nullptr : OrigValues[3];
  double f = coeff[d+2];
  double *Rhs = LocRhs[0];
  double eps = coeff[0];
  double b1 = coeff[1];
  double b2 = coeff[2];
  double b3 = d==2 ? 0.0 : coeff[3];
  double c = coeff[d+1];
  double delta;
  int sdfem_type = TDatabase::ParamDB->SDFEM_TYPE;
  double theta1 = TDatabase::TimeDB->THETA1;
  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  if(d==2)
  {
    double c00 = theta1 * tau * eps;
    double c11 = theta1 * tau * b1;
    double c22 = theta1 * tau * b2;
    double c33 = 1.0 + theta1 * tau * c;
    if(sdfem_type==8)
    {
      c33 = theta1 * tau * c;
    }
    double bb;
    if(fabs(c11) > fabs(c22))
      bb = fabs(c11);
    else
      bb = fabs(c22);
    if(sdfem_type != 2 && sdfem_type != 12)
      delta = Compute_SDFEM_delta<2>(hK, c00, {{c11, c22}}, c33, bb);
    else
      delta = Compute_SDFEM_delta<2>(hK, eps, {{b1, b2}}, c, bb);
  }
  else
  {
    double b_norm = coeff[6];
    delta = Compute_SDFEM_delta<3>(hK, eps, {{b1, b2, b3}}, c, b_norm);
  }
  double val  = theta1 * tau;
  if((sdfem_type!=2) && ( sdfem_type != 9)  && (sdfem_type != 10)
    && (sdfem_type != 11) && ( sdfem_type != 12))
    delta *= val;

  for(int i=0; i<N_; i++)
  {
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d==2 ? 0.0 : u_z[i];
    double bgradv = delta * (b1 * test_x + b2 * test_y + b3 * test_z);
    Rhs[i] += Mult * f * bgradv;
  }
}

template<int d>
void TCDGradGrad(double Mult, const double *, const double *, double,
                 const double**OrigValues, const int *N_BaseFuncts,
                 double ***LocMatrices, double **)
{
  double **Matrix = LocMatrices[0];
  int N_ = N_BaseFuncts[0];
  const double * u_x=OrigValues[1];
  const double * u_y=OrigValues[2];
  const double * u_z = d == 2 ? nullptr : OrigValues[3];

  for(int i=0; i<N_; i++)
  {
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0. : u_z[i];
    for(int j=0; j<N_; j++)
    {
      double ansatz_x=u_x[j];
      double ansatz_y=u_y[j];
      double ansatz_z = d == 2 ? 0. : u_z[j];
      double val = (test_x*ansatz_x + test_y*ansatz_y + test_z*ansatz_z);
      Matrix[i][j] += Mult * val;
    }
  }
}

#ifdef __2D__
template void TCDStiff<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void TCDStiffWithoutDiffAndLinearSource<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void TCDDiffAndLinearSourceOnly<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, int matrixIndex);
template void TCDStiff_TensorialDiffusionTerm<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void TCDMass<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void TCDRhs<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void TCDStiffSUPG<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void TCDMassSUPG<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void TCDRhsSUPG<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void TCDGradGrad<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
#endif // 2D
#ifdef __3D__
template void TCDStiff<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void TCDStiffWithoutDiffAndLinearSource<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void TCDDiffAndLinearSourceOnly<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs, int matrixIndex);
template void TCDStiff_TensorialDiffusionTerm<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void TCDMass<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void TCDRhs<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void TCDStiffSUPG<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void TCDMassSUPG<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void TCDRhsSUPG<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void TCDGradGrad<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double**OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
#endif // 3D
