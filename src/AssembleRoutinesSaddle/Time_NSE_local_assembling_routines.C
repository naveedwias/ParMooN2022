#include "Time_NSE_local_assembling_routines.h"
#include "Database.h"
#include "CommonRoutineTNSE3D.h"
#include "PointwiseAssemblyData.h"
#include <MooNMD_Io.h>
#include <cmath>

template <int d>
void NSMassMatrixSingle(double Mult, const double *, const double *, double,
                        const double **OrigValues, const int *N_BaseFuncts,
                        double ***LocMatrices, double **)
{
  double ** MatrixM = LocMatrices[0];
  int N_U = N_BaseFuncts[0];
  const double * u = OrigValues[0];  
  
  for(int i = 0; i < N_U; i++)
  {
    double test = u[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz = u[j];
      MatrixM[i][j] += Mult * (test * ansatz);
    }
  }
}

template<int d>
void NSMassMatrix(double Mult, const double *, const double *, double,
                  const double **OrigValues, const int *N_BaseFuncts,
                  double ***LocMatrices, double **)
{
  double ** MatrixM11 = LocMatrices[0];
  double ** MatrixM22 = LocMatrices[1];
  double ** MatrixM33 = d == 2 ? nullptr : LocMatrices[2];
  int N_U = N_BaseFuncts[0];
  const double * u = OrigValues[0];
  
  for(int i = 0; i < N_U; i++)
  {
    double test = u[i];
    for(int j = 0; j < N_U; j++)
    {
      double ansatz = u[j];
      double val = Mult * (test * ansatz);
      MatrixM11[i][j] += val;
      MatrixM22[i][j] += val;
      if(d == 3)
        MatrixM33[i][j] += val;
    }
  }
}

template<int d> 
void NSLaplaceGradGradSingleSmagorinsky(double Mult, const double *,
                                        const double *param, double hK,
                                        const double **OrigValues,
                                        const int *N_BaseFuncts,
                                        double ***LocMatrices, double **)
{
  if(d==2)
    ErrThrow("smagorinsky is not supported in 2D yet");
  double ** MatrixA = LocMatrices[0];
  int N_U = N_BaseFuncts[0];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = d == 2 ? nullptr : OrigValues[4];
  const double *x = &param[12];
  const double *y = &param[13];
  const double *z = &param[14];
  const double *u = &param[0];
  const double *gradu = &param[3];
  const double *uConv = &param[0];
  
  double mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711);
  
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
      MatrixA[i][j] += Mult * mu * (test_x * ansatz_x + test_y * ansatz_y
                                 + test_z * ansatz_z);
    }
  }
}

template <int d>
void NSLaplaceGradGradSmagorinsky(double Mult, const double *,
                                  const double *param, double hK,
                                  const double **OrigValues,
                                  const int *N_BaseFuncts,
                                  double ***LocMatrices, double **)
{
  if(d==2)
    ErrThrow("smagorinsky is not supported in 2D yet");
  double ** MatrixA11 = LocMatrices[0];
  double ** MatrixA22 = LocMatrices[d+1];
  double ** MatrixA33 = d == 2 ? nullptr : LocMatrices[8];
  int N_U = N_BaseFuncts[0];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = d == 2 ? nullptr : OrigValues[4];
  
  const double *x = &param[12];
  const double *y = &param[13];
  const double *z = &param[14];
  const double *u = &param[0];
  const double *gradu = &param[3];
  const double *uConv = &param[0];
  double mu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711);
  
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
      double val = Mult * mu* ( test_x * ansatz_x + test_y * ansatz_y
                             + test_z * ansatz_z);
      MatrixA11[i][j] += val;
      MatrixA22[i][j] += val;
      if(d == 3)
        MatrixA33[i][j] += val;
    }
  }
}

template <int d>
void NSLaplaceDeformationSmagorinsky(double Mult, const double *,
                                     const double *param, double hK,
                                     const double **OrigValues,
                                     const int *N_BaseFuncts,
                                     double ***LocMatrices, double **)
{
  if(d==2)
    ErrThrow("smagorinsky is not supported in 2D yet");
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
  
  const double *x = &param[12];
  const double *y = &param[13];
  const double *z = &param[14];
  const double *u = &param[0];
  const double *gradu = &param[3];
  const double *uConv = &param[0];
  double nu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, -4711)/2.;
  
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
void NSLaplaceDeformationVariationalMS(double Mult, const double *,
                                       const double *param, double hK,
                                       const double**OrigValues,
                                       const int *N_BaseFuncts,
                                       double ***LocMatrices, double **)
{
  if(d==2)
    ErrThrow("Variational MultiScale method is not supported in 2D yet");
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double **MatrixA13 = d == 2 ? nullptr : LocMatrices[2];
  double **MatrixA21 = LocMatrices[d];
  double **MatrixA22 = LocMatrices[d+1];
  double **MatrixA23 = d == 2 ? nullptr : LocMatrices[d+2];
  double **MatrixA31 = d == 2 ? nullptr : LocMatrices[6];
  double **MatrixA32 = d == 2 ? nullptr : LocMatrices[7];
  double **MatrixA33 = d == 2 ? nullptr : LocMatrices[8];
  int N_U = N_BaseFuncts[0];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = d == 2 ? nullptr : OrigValues[4];
  
  const double *x = &param[12];
  const double *y = &param[13];
  const double *z = &param[14];
  const double *u = &param[0];
  const double *gradu = &param[3];
  const double *uConv = &param[0];
  const double *projection_space_label = &param[15];
  //VMS0: at the moment this is same as the smagorinsky model
  double nu_t = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, projection_space_label[0]);
  nu_t = nu_t/2.;
 
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
      MatrixA11[i][j] += Mult * 2*nu_t*(test_x*ansatz_x + 0.5*test_y*ansatz_y 
                                     +0.5*test_z*ansatz_z);
      MatrixA12[i][j] += Mult * nu_t*(test_y*ansatz_x);
      if(d == 3)
        MatrixA13[i][j] += Mult * nu_t*(test_z*ansatz_x);
      MatrixA21[i][j] += Mult * nu_t*(test_x*ansatz_y);
      MatrixA22[i][j] += Mult * 2*nu_t*(0.5*test_x*ansatz_x+test_y*ansatz_y
                                     +0.5*test_z*ansatz_z);
      if(d == 3)
      {
        MatrixA23[i][j] += Mult * nu_t*(test_z*ansatz_y);
        MatrixA31[i][j] += Mult * nu_t*(test_x*ansatz_z);
        MatrixA32[i][j] += Mult * nu_t*(test_y*ansatz_z);
        MatrixA33[i][j] += Mult * 2*nu_t*(0.5*test_x*ansatz_x+0.5*test_y*ansatz_y
                                       +test_z*ansatz_z);
      }
    }
  }
}


template<>
void NSParamsVariationalMSLargeScale<2>(const double *, double *)
{
  ErrThrow("VMS Model in 2-dimensional case is not supported so far!!");
}

template<int d>
void NSLumpMassMatrix(double Mult, const double *, const double *, double ,
                      const double**OrigValues, const int *N_BaseFuncts,
                      double ***LocMatrices, double **)
{
  if(d==2)
    ErrThrow("Variational MultiScale method is not supported in 2D yet");
  double **MatrixL = LocMatrices[9];
  const double *l = OrigValues[5];
  int N_ = N_BaseFuncts[2];
  for(int i=0;i<N_;i++)
  {
     double test = l[i];     
     for(int j=0;j<N_;j++)
     {
        double ansatz = l[j];
        MatrixL[i][j] += Mult * test * ansatz;
     }
  }

}
template<int d>
void NSVariationlMS_GMatrices(double Mult, const double *, const double *,
                              double , const double**OrigValues,
                              const int *N_BaseFuncts, double ***LocMatrices,
                              double **)
{
  if(d==2)
    ErrThrow("Variational MultiScale method is not supported in 2D yet");
  double **MatrixG11 = LocMatrices[19];
  double **MatrixG22 = LocMatrices[20];
  double **MatrixG33 = LocMatrices[21];
  
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = OrigValues[4];
  const double *l = OrigValues[5];
  
  int N_U = N_BaseFuncts[0];
  int N_L = N_BaseFuncts[2];
  
  for(int i=0;i<N_L;i++)
  {
    double test = l[i];     
    double val =  Mult * test;
    for(int j=0;j<N_U;j++)
    {
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = u_z[j];
      
      MatrixG11[i][j] -= val * ansatz_x;
      MatrixG22[i][j] -= val * ansatz_y;
      MatrixG33[i][j] -= val * ansatz_z;
    }
  }
}

template<int d>
void NSVariationlMS_GTildeMatrices(double Mult, const double *,
                                   const double *param, double hK,
                                   const double**OrigValues,
                                   const int *N_BaseFuncts,
                                   double ***LocMatrices, double **)
{
  if(d==2)
    ErrThrow("Variational MultiScale method is not supported in 2D yet");
  double **MatrixGT11 = LocMatrices[16];
  double **MatrixGT22 = LocMatrices[17];
  double **MatrixGT33 = LocMatrices[18];
  
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = OrigValues[4];
  const double *l = OrigValues[5];
  
  int N_U = N_BaseFuncts[0];
  int N_L = N_BaseFuncts[2];
  
  const double *x = &param[12];
  const double *y = &param[13];
  const double *z = &param[14];
  const double *u = &param[0];
  const double *gradu = &param[3];
  const double *uConv = &param[0];
  const double *projection_space_label = &param[15];
  //VMS0: at the moment this is same as the smagorinsky model
  double nu = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, projection_space_label[0])/2.;
  
  for(int i=0;i<N_U;i++)
  {
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = u_z[i];
    
    for(int j=0;j<N_L;j++)
    {
      double ansatz = l[j];
      double val = Mult * 2. * nu * ansatz;
      MatrixGT11[i][j] -= val * test_x;
      MatrixGT22[i][j] -= val * test_y;
      MatrixGT33[i][j] -= val * test_z;
    }
  }
}


template<>
void NSParamVelGradSmagorinsky<2>(const double * /*in*/, double * /*out*/)
{
  ErrThrow("Smagorinsky Model in 2-dimensional case is not supported so far!!");
}

template<>
void NSParamVelGradSmagorinsky<3>(const double *in, double *out)
{
  out[0] = in[3]; // u1old
  out[1] = in[4]; // u2old
  out[2] = in[5]; // u3old
  
  out[3] = in[6]; // D1u1
  out[4] = in[7]; // D1u2
  out[5] = in[8]; // D1u3
  out[6] = in[9]; // D2u1
  out[7] = in[10]; // D2u2
  out[8] = in[11]; // D2u3
  out[9] = in[12]; // D3u1
  out[10] = in[13]; // D3u2
  out[11] = in[14]; // D3u3
  
  out[12] = in[0]; // x - coordinate for van Driest damping
  out[13] = in[1]; // y - coordinate for van Driest damping
  out[14] = in[2]; // z - coordinate for van Driest damping
}

template<>
void NSParamsVariationalMSLargeScale<3>(const double *in, double *out)
{
  out[0] = in[3]; // u1old
  out[1] = in[4]; // u2old
  out[2] = in[5]; // u3old

  out[3] = in[6]; // D1u1
  out[4] = in[7]; // D1u2
  out[5] = in[8]; // D1u3
  out[6] = in[9]; // D2u1
  out[7] = in[10]; // D2u2
  out[8] = in[11]; // D2u3
  out[9] = in[12]; // D3u1
  out[10] = in[13]; // D3u2
  out[11] = in[14]; // D3u3

  out[12] = in[0]; // x - coordinate for van Driest damping
  out[13] = in[1]; // y - coordinate for van Driest damping
  out[14] = in[2]; // z - coordinate for van Driest damping
  
  out[15] = in[15]; // projection space label
}


template <int d>
void NS_SUPG(double Mult, const double *coeff, const double *param, 
                    double hK, const double **OrigValues, const int *N_BaseFuncts,
                    double ***LocMatrices, double **, double delta0, double delta1)
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
  // original values: u,p,u_x,u_y[,u_z]
  const double *u = OrigValues[0]; // u
  const double *u_x = OrigValues[2]; // u_x
  const double *u_y = OrigValues[3]; // u_y
  const double * u_z = d == 2 ? nullptr : OrigValues[4];
    
  int N_U = N_BaseFuncts[0];
  
  double c0=coeff[0];
  
  const double u1=param[0];
  const double u2=param[1];
  const double u3= d == 2 ? 0 : param[2];

  double test_x, test_y, test_z, test;
  double ansatz_x, ansatz_y, ansatz_z;
  
  // stabilization parameters
  double tau_m, tau_c;
  if(TDatabase::ParamDB->NSTYPE==4)
  {
    tau_m = delta0*hK*hK;
    tau_c = delta1;
  }
  else
  {
    double stab_params[2];// tau_m = stab_params[0]; tau_c=stab_params[1];
    compute_stabilization_parameters<d>(param, coeff, stab_params);
    tau_m = stab_params[0];
    tau_c = stab_params[1];
  }
  
  for(int i=0; i<N_U; ++i)
  {
    test = u[i];
    test_x = u_x[i];
    test_y = u_y[i];
    test_z = d == 2 ? 0. : u_z[i];
    double ugradv = tau_m * (u1*test_x+u2*test_y+u3*test_z);
    
    // velocity-velocity blocks
    for(int j=0; j<N_U; ++j)
    {
      ansatz_x = u_x[j];
      ansatz_y = u_y[j];
      ansatz_z = d == 2 ? 0 : u_z[j];
      // Galerkin part
      double val  = c0*(2*test_x*ansatz_x +test_y*ansatz_y +test_z*ansatz_z); // diffusion term
      double conv = (u1*ansatz_x + u2*ansatz_y + u3*ansatz_z)*test; // convective term
      val += conv;
      // supg contribution
      val += (u1*ansatz_x+u2*ansatz_y+u3*ansatz_z)*ugradv; 
      // grad div contribution
      val += tau_c*test_x*ansatz_x;
      MatrixA11[i][j] += Mult * val;
            
      val  = c0*(test_y*ansatz_x);
      // grad div contribution
      val += tau_c * test_x * ansatz_y;
      MatrixA12[i][j] += Mult * val;
      if(d==3)
      {
        val  = c0*(test_z*ansatz_x);
        // grad div contribution
        val += tau_c * test_x * ansatz_z;
        MatrixA13[i][j] += Mult * val;
      }

      val  = c0*(test_x*ansatz_y);
      // grad div contribution
      val += tau_c * test_y * ansatz_x;
      MatrixA21[i][j] += Mult * val;

      val  = 2*c0*(0.5*test_x*ansatz_x+test_y*ansatz_y
                   +0.5*test_z*ansatz_z);
      // convective
      val += conv;
      // supg contribution
      val += (u1*ansatz_x+u2*ansatz_y+u3*ansatz_z)*ugradv; 
      // grad div contribution
      val += tau_c * test_y * ansatz_y;
      MatrixA22[i][j] += Mult * val;

      if(d==3)
      {
        val  = c0*(test_z*ansatz_y) + tau_c * test_y * ansatz_z;
        MatrixA23[i][j] += Mult * val;
      
        val  = c0*(test_x*ansatz_z);
        // grad div contribution
        val += tau_c * test_z * ansatz_x;
        MatrixA31[i][j] += Mult * val;
        
        val  = c0*(test_y*ansatz_z);
        // grad div contribution
        val += tau_c * test_z * ansatz_y;
        MatrixA32[i][j] += Mult * val;

        val  = 2*c0*(0.5*test_x*ansatz_x +0.5*test_y*ansatz_y +test_z*ansatz_z);
        // convective term
        val += (u1*ansatz_x+u2*ansatz_y+u3*ansatz_z)*test;
        // supg contribution
        val += (u1*ansatz_x+u2*ansatz_y+u3*ansatz_z)*ugradv;
        // grad div contribution
        val += tau_c * test_z * ansatz_z;
        MatrixA33[i][j] += Mult * val;
      }
    }
  }
}


template<int d>
void NS_SUPG_GradientBlocks(double Mult, const double *, const double *param, double hK, 
                             const double **OrigValues, const int *N_BaseFuncts, 
                             double ***LocMatrices, double **, double delta0)
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
  
  const double u1=param[0];
  const double u2=param[1];
  const double u3= d == 2 ? 0 : param[2];
  
  double tau_m = delta0*hK*hK;
  
  for(int i = 0; i < N_U; i++)
  {
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0. : u_z[i];
    
    double ugradv = tau_m * (u1*test_x+u2*test_y+u3*test_z);
    
    for(int j = 0; j < N_P; j++)
    {
      double ansatz = p[j];
      double ansatz_x = p_x[j];
      double ansatz_y = p_y[j];
      double ansatz_z = d==2 ? 0.0 : p_z[j];
      MatrixB1T[i][j] -= Mult * (ansatz * test_x + ansatz_x * ugradv);
      MatrixB2T[i][j] -= Mult * (ansatz * test_y + ansatz_y * ugradv);;
      if(d == 3)
        MatrixB3T[i][j] -= Mult * (ansatz * test_z + ansatz_z * ugradv);;
    }
  }
}


template<int d> 
void NS_SUPG_RightHandSide_InfSup(double Mult, const double *coeff, const double *param, 
                           double hK, const double **OrigValues, 
                           const int *N_BaseFuncts, double ***, 
                           double **LocRhs, double delta1)
{
  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];
  double *Rhs3 = d==2 ? nullptr : LocRhs[2];
  
  const double *u = OrigValues[0]; // u
  const double *u_x = OrigValues[2]; // u_x
  const double *u_y = OrigValues[3]; // u_y
  const double *u_z = d==2 ? nullptr : OrigValues[4]; // u_z
  
  double f1 = coeff[1];
  double f2 = coeff[2];
  double f3 = d==2 ? 0. : coeff[3];

  const double u1=param[0];
  const double u2=param[1];
  const double u3= d==2 ? 0. : param[2];

  int N_U = N_BaseFuncts[0];
  double tau_m = delta1*hK*hK;
  for(int i=0; i<N_U; ++i)
  {
    double test = u[i];
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0. : u_z[i];
    
    double ugradv = tau_m * (u1*test_x+u2*test_y+u3*test_z);
    
    Rhs1[i] += Mult*f1*(test + ugradv);
    Rhs2[i] += Mult*f2*(test + ugradv);
    if(d==3)
      Rhs3[i] += Mult*f3*(test + ugradv);
  }
}

template <int d> 
void NS_SUPG_MassMatrix(double Mult, const double *coeff, const double *param, 
                        double hK, const double **OrigValues, 
                        const int *N_BaseFuncts, double ***LocMatrices, 
                        double **, double delta0)
{
  double ** MatrixM11 = LocMatrices[0];
  double ** MatrixM22 = LocMatrices[1];
  double ** MatrixM33 = d == 2 ? nullptr : LocMatrices[2];
  int N_U = N_BaseFuncts[0];
  
  const double * u = OrigValues[0];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = d==2 ? nullptr : OrigValues[4];
  const double u1=param[0];
  const double u2=param[1];
  const double u3= d==2 ? 0. : param[2];
  
  double stab_params[2];
  // find other solution for that
  if(TDatabase::ParamDB->NSTYPE==4)
    stab_params[0] = delta0*hK*hK;
  else // tau_m = stab_params[0]; tau_c=stab_params[1];
    compute_stabilization_parameters<d>(param, coeff, stab_params);
      
  for(int i = 0; i < N_U; i++)
  {
    double test = u[i];
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0. : u_z[i];
    
    double ugradv = stab_params[0] * (u1*test_x+u2*test_y+u3*test_z);
    for(int j = 0; j < N_U; j++)
    {
      double ansatz = u[j];
      double val = Mult * (test + ugradv )*ansatz;
      MatrixM11[i][j] += val;
      MatrixM22[i][j] += val;
      if(d == 3)
        MatrixM33[i][j] += val;
    }
  }
}

template<int d>
void NS_SUPG_EquOrder_Gradient_DivergenceBlocks(
  double Mult, const double *coeff, const double *param, double , 
  const double **OrigValues, const int *N_BaseFuncts, 
  double ***LocMatrices, double **, double factor)
{
  double ** MatrixB1 = LocMatrices[d == 2 ? 5 : 10];
  double ** MatrixB2 = LocMatrices[d == 2 ? 6 : 11];
  double ** MatrixB3 = d == 2 ? nullptr : LocMatrices[12];  
  double ** MatrixB1T = LocMatrices[d == 2 ? 7 : 13];
  double ** MatrixB2T = LocMatrices[d == 2 ? 8 : 14];
  double ** MatrixB3T = d == 2 ? nullptr : LocMatrices[15];
  double ** MatrixC = LocMatrices[d == 2 ? 4 : 9];
  
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];
  
  const double * u = OrigValues[0];
  const double * p = OrigValues[1];
  const double * u_x = OrigValues[2];
  const double * u_y = OrigValues[3];
  const double * u_z = d == 2 ? nullptr : OrigValues[4];
  const double * p_x = OrigValues[2+d];
  const double * p_y = OrigValues[3+d];
  const double * p_z = d == 2 ? nullptr : OrigValues[4+d];
  
  const double u1=param[0];
  const double u2=param[1];
  const double u3= d == 2 ? 0 : param[2];
  
  double stab_params[2]; // tau_m = stab_params[0]; tau_c=stab_params[1];
  compute_stabilization_parameters<d>(param, coeff, stab_params);
  double tau_m = stab_params[0];
  
  double ansatz, ansatz_x, ansatz_y, ansatz_z;
  double test, test_x, test_y, test_z;
  for(int i = 0; i < N_U; i++)
  {
    test_x = u_x[i];
    test_y = u_y[i];
    test_z = d == 2 ? 0. : u_z[i];
    
    double ugradv = stab_params[0] * (u1*test_x+u2*test_y+u3*test_z);
    
    for(int j = 0; j < N_P; j++)
    {
      ansatz = p[j];
      ansatz_x = p_x[j];
      ansatz_y = p_y[j];
      ansatz_z = d==2 ? 0.0 : p_z[j];
      MatrixB1T[i][j] -= Mult * (ansatz * test_x + ansatz_x * ugradv);
      MatrixB2T[i][j] -= Mult * (ansatz * test_y + ansatz_y * ugradv);;
      if(d == 3)
        MatrixB3T[i][j] -= Mult * (ansatz * test_z + ansatz_z * ugradv);;
    }
  }
  
  for(int i=0; i<N_P; ++i)
  { 
    test   = p[i];
    test_x = p_x[i];
    test_y = p_y[i];
    test_z = d == 2 ? 0. : p_z[i];

    // velocity-pressure block
    for(int j=0;j<N_U;j++)
    {
      // velocity ansatz functions
      ansatz   = u[j];
      ansatz_x = u_x[j];
      ansatz_y = u_y[j];
      ansatz_z = d==2 ? 0. : u_z[j];

      double val = -test*ansatz_x;
      // supg contribution
      val -= tau_m*(ansatz * factor + u1*ansatz_x + u2 * ansatz_y + u3*ansatz_z)*test_x;
      MatrixB1[i][j] -= Mult*val;

      val = -test*ansatz_y;
      // supg contribution
      val -= tau_m*(ansatz * factor + u1*ansatz_x + u2 * ansatz_y + u3*ansatz_z)*test_y;
      MatrixB2[i][j] -= Mult*val;
      if(d==3)
      {
        val = -test*ansatz_z;
        val -= tau_m*(ansatz * factor + u1*ansatz_x + u2 * ansatz_y + u3*ansatz_z)*test_z;
        MatrixB3[i][j] -= Mult*val;
      }
    }
    // pressure-pressure block
    for(int j=0; j<N_P; j++)
    {
      // pressure ansatz
      ansatz_x = p_x[j];
      ansatz_y = p_y[j];
      ansatz_z = d == 2 ? 0. : p_z[j];

      double val = tau_m * (ansatz_x * test_x + ansatz_y * test_y + ansatz_z * test_z);
      MatrixC[i][j] += Mult*val;
    }
  }
}

template<int d> 
void NS_SUPG_RightHandSide_EquOrder(double Mult, const double *coeff, 
         const double *param, double /*hK*/, const double **OrigValues, 
         const int *N_BaseFuncts, double ***, double **LocRhs, double factor)
{
  double *Rhs1 = LocRhs[0];
  double *Rhs2 = LocRhs[1];
  double *Rhs3 = d==2 ? nullptr : LocRhs[2];
  double *Rhs_div = LocRhs[d];
  
  const double *u = OrigValues[0]; // u
  const double *u_x = OrigValues[2]; // u_x
  const double *u_y = OrigValues[3]; // u_y
  const double * u_z = d == 2 ? nullptr : OrigValues[4];
  // const double * p   = OrigValues[1];
  const double * p_x = OrigValues[2+d];
  const double * p_y = OrigValues[3+d];
  const double * p_z = d == 2 ? nullptr : OrigValues[4+d];
  
  double f1 = coeff[1];
  double f2 = coeff[2];
  double f3 = d==2 ? 0. : coeff[3];

  const double u1=param[0];
  const double u2=param[1];
  const double u3= d==2 ? 0. : param[2];
  
  const double u1_old_time = param[d];
  const double u2_old_time = param[d+1];
  const double u3_old_time = d==2 ? 0. : param[5];

  int N_U = N_BaseFuncts[0];
  // stabilization parameters
  double stab_params[2];
  compute_stabilization_parameters<d>(param, coeff, stab_params);
  double tau_m = stab_params[0];
  
  for(int i=0; i<N_U; ++i)
  {
    double test = u[i];
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0. : u_z[i];
    
    double ugradv = tau_m * (u1*test_x+u2*test_y+u3*test_z);
    
    Rhs1[i] += Mult*f1*(test + ugradv);
    Rhs2[i] += Mult*f2*(test + ugradv);
    if(d==3)
      Rhs3[i] += Mult*f3*(test + ugradv);
  }
  
  // double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH; 
  int N_P = N_BaseFuncts[1];
  // pressure test functions
  for(int i=0; i<N_P; ++i)
  {
    // pressure test
    double test_x = p_x[i];
    double test_y = p_y[i];
    double test_z = d==2 ? 0. : p_z[i];
    
    double val  = (u1_old_time/factor + f1)*test_x;
    val += (u2_old_time/factor +f2)*test_y;
    if(d==3)
      val += (u3_old_time/factor +f3)*test_z;
    Rhs_div[i] += Mult * stab_params[0] * val;
  }
}

template <int d>
void NS_SUPG_skew_symmetric(double Mult, const double *coeff, 
        const double *param, double hK, const double **OrigValues, 
        const int *N_BaseFuncts, double ***LocMatrices, double **, 
        double delta0, double delta1)
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
  // original values: u,p,u_x,u_y[,u_z]
  const double *u = OrigValues[0]; // u
  const double *u_x = OrigValues[2]; // u_x
  const double *u_y = OrigValues[3]; // u_y
  const double * u_z = d == 2 ? nullptr : OrigValues[4];
    
  int N_U = N_BaseFuncts[0];
  
  double c0=coeff[0];
  
  const double u1=param[0];
  const double u2=param[1];
  const double u3= d == 2 ? 0 : param[2];

  double test_x, test_y, test_z, test;
  double ansatz_x, ansatz_y, ansatz_z, ansatz;
  
  
  // stabilization parameters
  double tau_m, tau_c;
  if(TDatabase::ParamDB->NSTYPE==4)
  {
    tau_m = delta0*hK*hK;
    tau_c = delta1;
  }
  else
  {
    double stab_params[2];// tau_m = stab_params[0]; tau_c=stab_params[1];
    compute_stabilization_parameters<d>(param, coeff, stab_params);
    tau_m = stab_params[0];
    tau_c = stab_params[1];
  }
  
  for(int i=0; i<N_U; ++i)
  {
    test = u[i];
    test_x = u_x[i];
    test_y = u_y[i];
    test_z = d == 2 ? 0. : u_z[i];
    double ugradv = tau_m * (u1*test_x+u2*test_y+u3*test_z);
    
    // velocity-velocity blocks
    for(int j=0; j<N_U; ++j)
    {
      ansatz = u[j];
      ansatz_x = u_x[j];
      ansatz_y = u_y[j];
      ansatz_z = d == 2 ? 0 : u_z[j];
      // Galerkin part
      double val  = c0*(2*test_x*ansatz_x +test_y*ansatz_y +test_z*ansatz_z); // diffusion term
      double conv = (u1*ansatz_x + u2*ansatz_y + u3*ansatz_z)*test; // convective term
      conv -= (u1*test_x + u2*test_y + u3*test_z) * ansatz;
      conv *= 0.5;
      val += conv;
      // supg contribution
      val += (u1*ansatz_x+u2*ansatz_y+u3*ansatz_z)*ugradv; 
      // grad div contribution
      val += tau_c*test_x*ansatz_x;
      MatrixA11[i][j] += Mult * val;
            
      val  = c0*(test_y*ansatz_x);
      // grad div contribution
      val += tau_c * test_x * ansatz_y;
      MatrixA12[i][j] += Mult * val;
      if(d==3)
      {
        val  = c0*(test_z*ansatz_x);
        // grad div contribution
        val += tau_c * test_x * ansatz_z;
        MatrixA13[i][j] += Mult * val;
      }

      val  = c0*(test_x*ansatz_y);
      // grad div contribution
      val += tau_c * test_y * ansatz_x;
      MatrixA21[i][j] += Mult * val;

      val  = 2*c0*(0.5*test_x*ansatz_x+test_y*ansatz_y
                   +0.5*test_z*ansatz_z);
      // convective
      val += conv;
      // supg contribution
      val += (u1*ansatz_x+u2*ansatz_y+u3*ansatz_z)*ugradv; 
      // grad div contribution
      val += tau_c * test_y * ansatz_y;
      MatrixA22[i][j] += Mult * val;

      if(d==3)
      {
        val  = c0*(test_z*ansatz_y) + tau_c * test_y * ansatz_z;
        MatrixA23[i][j] += Mult * val;
      
        val  = c0*(test_x*ansatz_z);
        // grad div contribution
        val += tau_c * test_z * ansatz_x;
        MatrixA31[i][j] += Mult * val;
        
        val  = c0*(test_y*ansatz_z);
        // grad div contribution
        val += tau_c * test_z * ansatz_y;
        MatrixA32[i][j] += Mult * val;

        val  = 2*c0*(0.5*test_x*ansatz_x +0.5*test_y*ansatz_y +test_z*ansatz_z);
        // convective term
        val += conv;
        // supg contribution
        val += (u1*ansatz_x+u2*ansatz_y+u3*ansatz_z)*ugradv;
        // grad div contribution
        val += tau_c * test_z * ansatz_z;
        MatrixA33[i][j] += Mult * val;
      }
    }
  }
}


template <int d>
void NS_SUPG_rotational(double Mult, const double *coeff, 
        const double *param, double hK, const double **OrigValues, 
        const int *N_BaseFuncts, double ***LocMatrices, double **, 
        double delta0, double delta1)
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
  // original values: u,p,u_x,u_y[,u_z]
  const double *u = OrigValues[0]; // u
  const double *u_x = OrigValues[2]; // u_x
  const double *u_y = OrigValues[3]; // u_y
  const double * u_z = d == 2 ? nullptr : OrigValues[4];
    
  int N_U = N_BaseFuncts[0];
  
  double c0=coeff[0];
  
  const double u1=param[0];
  const double u2=param[1];
  const double u3= d == 2 ? 0 : param[2];

  double test_x, test_y, test_z, test;
  double ansatz_x, ansatz_y, ansatz_z;
  
  // stabilization parameters
  double tau_m, tau_c;
  if(TDatabase::ParamDB->NSTYPE==4)
  {
    tau_m = delta0*hK*hK;
    tau_c = delta1;
  }
  else
  {
    double stab_params[2];// tau_m = stab_params[0]; tau_c=stab_params[1];
    compute_stabilization_parameters<d>(param, coeff, stab_params);
    tau_m = stab_params[0];
    tau_c = stab_params[1];
  }
  
  for(int i=0; i<N_U; ++i)
  {
    test = u[i];
    test_x = u_x[i];
    test_y = u_y[i];
    test_z = d == 2 ? 0. : u_z[i];
    double ugradv = tau_m * (u1*test_x+u2*test_y+u3*test_z);
    
    // velocity-velocity blocks
    for(int j=0; j<N_U; ++j)
    {
      ansatz_x = u_x[j];
      ansatz_y = u_y[j];
      ansatz_z = d == 2 ? 0 : u_z[j];
      // Galerkin part
      double val  = c0*(2*test_x*ansatz_x +test_y*ansatz_y +test_z*ansatz_z); // diffusion term
      val += (u2*ansatz_y + u3*ansatz_z) * test; // convective term
      // supg contribution
      val += (u1*ansatz_x+u2*ansatz_y+u3*ansatz_z)*ugradv; 
      // grad div contribution
      val += tau_c*test_x*ansatz_x;
      MatrixA11[i][j] += Mult * val;
            
      val  = c0*(test_y*ansatz_x);
      // grad div contribution
      val += tau_c * test_x * ansatz_y;
      val -= u2*ansatz_x * test;
      MatrixA12[i][j] += Mult * val;
      if(d==3)
      {
        val  = c0*(test_z*ansatz_x);
        // grad div contribution
        val += tau_c * test_x * ansatz_z;
        val -= u3*ansatz_x * test;
        MatrixA13[i][j] += Mult * val;
      }

      val  = c0*(test_x*ansatz_y);
      // convective term 
      val -= u1 * ansatz_y * test;
      // grad div contribution
      val += tau_c * test_y * ansatz_x;
      MatrixA21[i][j] += Mult * val;

      val  = 2*c0*(0.5*test_x*ansatz_x+test_y*ansatz_y
                   +0.5*test_z*ansatz_z);
      // convective term
      val += (u1*ansatz_x + u3*ansatz_z) * test;
      // supg contribution
      val += (u1*ansatz_x+u2*ansatz_y+u3*ansatz_z)*ugradv; 
      // grad div contribution
      val += tau_c * test_y * ansatz_y;
      MatrixA22[i][j] += Mult * val;

      if(d==3)
      {
        val  = c0*(test_z*ansatz_y) + tau_c * test_y * ansatz_z;
        // convective term
        val -= u3*ansatz_y * test;
        MatrixA23[i][j] += Mult * val;
        
        val  = c0*(test_x*ansatz_z);
        // convective term
        val -= u1*ansatz_z * test;
        // grad div contribution
        val += tau_c * test_z * ansatz_x;
        MatrixA31[i][j] += Mult * val;
        
        val  = c0*(test_y*ansatz_z);
        // convective term
        val -= u2 * ansatz_z * test;
        // grad div contribution
        val += tau_c * test_z * ansatz_y;
        MatrixA32[i][j] += Mult * val;

        val  = 2*c0*(0.5*test_x*ansatz_x +0.5*test_y*ansatz_y +test_z*ansatz_z);
        // convective term
        val += (u1*ansatz_x + u2*ansatz_y) * test;;
        // supg contribution
        val += (u1*ansatz_x+u2*ansatz_y+u3*ansatz_z)*ugradv;
        // grad div contribution
        val += tau_c * test_z * ansatz_z;
        MatrixA33[i][j] += Mult * val;
      }
    }
  }
}


template <int d>
void NS_SUPG_emac(double Mult, const double *coeff, 
        const double *param, double hK, const double **OrigValues, 
        const int *N_BaseFuncts, double ***LocMatrices, double **, 
        double delta0, double delta1)
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
  // original values: u,p,u_x,u_y[,u_z]
  const double *u = OrigValues[0]; // u
  const double *u_x = OrigValues[2]; // u_x
  const double *u_y = OrigValues[3]; // u_y
  const double * u_z = d == 2 ? nullptr : OrigValues[4];
    
  int N_U = N_BaseFuncts[0];
  
  double c0=coeff[0];
  
  const double u1=param[0];
  const double u2=param[1];
  const double u3= d == 2 ? 0 : param[2];

  double test_x, test_y, test_z, test;
  double ansatz_x, ansatz_y, ansatz_z;
  
  // stabilization parameters
  double tau_m, tau_c;
  if(TDatabase::ParamDB->NSTYPE==4)
  {
    tau_m = delta0*hK*hK;
    tau_c = delta1;
  }
  else
  {
    double stab_params[2];// tau_m = stab_params[0]; tau_c=stab_params[1];
    compute_stabilization_parameters<d>(param, coeff, stab_params);
    tau_m = stab_params[0];
    tau_c = stab_params[1];
  }
  
  for(int i=0; i<N_U; ++i)
  {
    test = u[i];
    test_x = u_x[i];
    test_y = u_y[i];
    test_z = d == 2 ? 0. : u_z[i];
    double ugradv = tau_m * (u1*test_x+u2*test_y+u3*test_z);
    
    // velocity-velocity blocks
    for(int j=0; j<N_U; ++j)
    {
      ansatz_x = u_x[j];
      ansatz_y = u_y[j];
      ansatz_z = d == 2 ? 0 : u_z[j];
      // Galerkin part
      double val  = c0*(2*test_x*ansatz_x +test_y*ansatz_y +test_z*ansatz_z); // diffusion term
      val += (3.*u1*ansatz_x + ansatz_y*u2 + ansatz_z*u3)*test; // convective term
      // supg contribution
      val += (u1*ansatz_x+u2*ansatz_y+u3*ansatz_z)*ugradv; 
      // grad div contribution
      val += tau_c*test_x*ansatz_x;
      MatrixA11[i][j] += Mult * val;
            
      val  = c0*(test_y*ansatz_x);
      // grad div contribution
      val += tau_c * test_x * ansatz_y;
      val += (ansatz_x*u2 + ansatz_y*u1) * test;
      MatrixA12[i][j] += Mult * val;
      if(d==3)
      {
        val  = c0*(test_z*ansatz_x);
        // grad div contribution
        val += tau_c * test_x * ansatz_z;
        val += (ansatz_x*u3 + ansatz_z*u1) * test;
        MatrixA13[i][j] += Mult * val;
      }

      val  = c0*(test_x*ansatz_y);
      // convective term 
      val += (ansatz_y*u1 + ansatz_x*u2) * test;
      // grad div contribution
      val += tau_c * test_y * ansatz_x;
      MatrixA21[i][j] += Mult * val;

      val  = 2*c0*(0.5*test_x*ansatz_x+test_y*ansatz_y
                   +0.5*test_z*ansatz_z);
      // convective term
      val += (ansatz_x*u1 + 3.*ansatz_y*u2 + ansatz_z*u3) * test;
      // supg contribution
      val += (u1*ansatz_x+u2*ansatz_y+u3*ansatz_z)*ugradv; 
      // grad div contribution
      val += tau_c * test_y * ansatz_y;
      MatrixA22[i][j] += Mult * val;

      if(d==3)
      {
        val  = c0*(test_z*ansatz_y) + tau_c * test_y * ansatz_z;
        // convective term
        val += (ansatz_y*u3 + ansatz_z*u2) * test;
        MatrixA23[i][j] += Mult * val;
        
        val  = c0*(test_x*ansatz_z);
        // convective term
        val += (ansatz_z*u1 + ansatz_x*u3) * test;
        // grad div contribution
        val += tau_c * test_z * ansatz_x;
        MatrixA31[i][j] += Mult * val;
        
        val  = c0*(test_y*ansatz_z);
        // convective term
        val += (ansatz_z*u2 + ansatz_y*u3) * test;
        // grad div contribution
        val += tau_c * test_z * ansatz_y;
        MatrixA32[i][j] += Mult * val;

        val  = 2*c0*(0.5*test_x*ansatz_x +0.5*test_y*ansatz_y +test_z*ansatz_z);
        // convective term
        val += ( ansatz_x*u1 + ansatz_y*u2 + 3.*ansatz_z*u3) * test;;
        // supg contribution
        val += (u1*ansatz_x+u2*ansatz_y+u3*ansatz_z)*ugradv;
        // grad div contribution
        val += tau_c * test_z * ansatz_z;
        MatrixA33[i][j] += Mult * val;
      }
    }
  }
}

template <int d>
void NSParamsVelocityDerivatives_SUPG_inf_sup(const double *in, double *out)
{
  out[0] = in[d];
  out[1] = in[d+1];
  if(d==3)
    out[2] = in[d+2];  
}
template <int d>
void NSParamsVelocityDerivatives_SUPG_equal_order(const double *in, double *out)
{
  out[0] = in[d];
  out[1] = in[d+1];
  out[2] = in[d+2];
  out[3] = in[d+3];
  if(d==3)
  {
    out[4] = in[d+4]; 
    out[5] = in[d+5];
  }
}

template <int d>
void NSParamsVMSResidualsLaplacian(const double *in, double *out)
{
  // expected input layout:
  // u, Du, Lu, Dp, u_old

  // + d indices for location
  // + d indices for u
  // + d * d indices for Du
  const int Lu1_ofs = (d + 2) * d;
  const int Lu2_ofs = (d + 3) * d;
  const int Lu3_ofs = (d + 4) * d;

  double u1_xx = in[Lu1_ofs + 0];
  double u1_yy = in[Lu1_ofs + 1];
  double u1_zz = d == 3 ? in[Lu1_ofs + 2] : 0.0;

  double u2_xx = in[Lu2_ofs + 0];
  double u2_yy = in[Lu2_ofs + 1];
  double u2_zz = d == 3 ? in[Lu2_ofs + 2] : 0.0;

  double u3_xx = d == 3 ? in[Lu3_ofs + 0] : 0.0;
  double u3_yy = d == 3 ? in[Lu3_ofs + 1] : 0.0;
  double u3_zz = d == 3 ? in[Lu3_ofs + 2] : 0.0;

  out[0] = u1_xx + u1_yy + u1_zz;
  out[1] = u2_xx + u2_yy + u2_zz;

  if (d == 3)
  {
    out[2] = u3_xx + u3_yy + u3_zz;
  }
}

template <int d>
void NSParamsVMSResidualsWithoutLaplacianAndRHS(const double *in, double *out, bool extend_advection)
{
  // expected input layout:
  // u, Du, Lu, Dp, u_old

  // + d indices for location
  const int u_ofs = d;

  // + d indices for u
  const int Du1_ofs = 2 * d;
  const int Du2_ofs = 3 * d;
  const int Du3_ofs = 4 * d;

  // + d * d indices for Du
  // + d * d indices for Lu
  const int Dp_ofs = (2 * d + 2) * d;

  // + d indices for Dp
  const int uo_ofs = (2 * d + 3) * d;

  double u1 = in[u_ofs + 0];
  double u2 = in[u_ofs + 1];
  double u3 = in[u_ofs + 2];

  double u1_x = in[Du1_ofs + 0];
  double u1_y = in[Du1_ofs + 1];
  double u1_z = d == 3 ? in[Du1_ofs + 2] : 0.0;

  double u2_x = in[Du2_ofs + 0];
  double u2_y = in[Du2_ofs + 1];
  double u2_z = d == 3 ? in[Du2_ofs + 2] : 0.0;

  double u3_x = d == 3 ? in[Du3_ofs + 0] : 0.0;
  double u3_y = d == 3 ? in[Du3_ofs + 1] : 0.0;
  double u3_z = d == 3 ? in[Du3_ofs + 2] : 0.0;

  double p_x = in[Dp_ofs + 0];
  double p_y = in[Dp_ofs + 1];
  double p_z = d == 3 ? in[Dp_ofs + 2] : 0.0;

  double u1_old = in[uo_ofs + 0];
  double u2_old = in[uo_ofs + 1];
  double u3_old = d == 3 ? in[uo_ofs + 2] : 0.0;

  double u1_t = (u1 - u1_old) / TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double u2_t = (u2 - u2_old) / TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double u3_t = (u3 - u3_old) / TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  if (extend_advection)
  {
    const double* u_prime = PointwiseAssemblyData::GetOldData();

    if (u_prime == nullptr)
    {
      ErrThrow("u' is not available!");
    }

    u1 += u_prime[0];
    u2 += u_prime[1];
    if (d == 3)
    {
      u3 += u_prime[2];
    }
  }

  out[0] = u1_t
    + u1 * u1_x + u2 * u1_y + u3 * u1_z
    + p_x;

  out[1] = u2_t
    + u1 * u2_x + u2 * u2_y + u3 * u2_z
    + p_y;

  if (d == 3)
  {
    out[2] = u3_t
      + u1 * u3_x + u2 * u3_y + u3 * u3_z
      + p_z;
  }

  out[d] = u1_x + u2_y + u3_z;
}

template <int d>
void NSParamsOldVelocity(const double *in, double *out)
{
  // expected input layout:
  // u, Du, Lu, Dp, u_old

  // + d indices for location
  // + d indices for u
  // + d * d indices for Du
  // + d * d indices for Lu
  // + d indices for Dp
  const int uo_ofs = (2 * d + 3) * d;

  out[0] = in[uo_ofs + 0];
  out[1] = in[uo_ofs + 1];

  if (d == 3)
  {
    out[2] = in[uo_ofs + 2];
  }
}

#define RBVMS_SGS_ALTERNATE_FORM

template<int d>
void
NS_GetVMSResiduals(const double *coeff, const double *param, double hK,
  int inv_jacobian_offset, double *res_m, double &res_c,
  double &tau_m, double &tau_c, RBVMS_Settings settings,
  bool premultiply)
{
  const double *Jinv = coeff + inv_jacobian_offset;
  const double *w = param + 0; // velocity starts at 0
  const double *res_m_lap = param + d; // laplacian
  const double *res_m_base = param + 2 * d; // res_m without laplacian and rhs

  double nu = coeff[0];

  switch (settings.mode)
  {
    case RBVMS_ParamMode::G:
    {
      double G[d * d];
      double g[d];

      for (int i = 0; i < d; i++)
      {
        g[i] = 0.0;

        for (int j = 0; j < d; j++)
        {
          g[i] += Jinv[d * i + j]; // d x^_j / d x_i

          int m = i + d * j; // G_ij
          G[m] = 0.0;
          for (int k = 0; k < d; k++)
          {
            // d x^_k / d x_i * d x^_k / d x_j
            G[m] += Jinv[d * i + k] * Jinv[d * j + k];
          }
        }
      }

      auto parameters = RBVMS_Param_G<d>(w, G, g,
        TDatabase::TimeDB->CURRENTTIMESTEPLENGTH, nu, settings);

      tau_m = parameters.first;
      tau_c = parameters.second;
      break;
    }

    case RBVMS_ParamMode::H:
    {
      auto parameters = RBVMS_Param_H<d>(hK, settings);

      tau_m = parameters.first;
      tau_c = parameters.second;
      break;
    }

    case RBVMS_ParamMode::Codina:
    {
      const double* u_prime = PointwiseAssemblyData::GetOldData();

      auto parameters = RBVMS_Param_Codina<d>(hK, nu, w, u_prime, settings);

      tau_m = parameters.first;
      tau_c = parameters.second;
      break;
    }

    default:
    {
      ErrThrow("Unsupported parameter mode!");

      tau_m = 0.0;
      tau_c = 0.0;
      break;
    }
  }

  if (premultiply)
  {
    // coeff[d + 1] is the divergence rhs
    res_c = tau_c * (param[3 * d] - coeff[d + 1]);
  }
  else
  {
    // coeff[d + 1] is the divergence rhs
    res_c = param[3 * d] - coeff[d + 1];
  }

  double mul = premultiply ? tau_m : 1.0;

  for (int i = 0; i < d; i++)
  {
    // rhs is in coeff[1] ff.
    res_m[i] = mul * (res_m_base[i] - nu * res_m_lap[i] - coeff[1 + i]);
  }
}

template<int d>
void NSVMSResiduals_MassMatrix(double Mult, const double *coeff, const double *param,
  double hK, const double **OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **, int inv_jacobian_offset,
  RBVMS_Settings settings)
{
  const bool residual_terms = !settings.explicit_time_derivative;
  const bool has_B_blocks = settings.momentum_pressure_coupling_B;

  double **MatrixM11 = LocMatrices[0];
  double **MatrixM12 = LocMatrices[1];
  double **MatrixM13 = d == 2 ? nullptr : LocMatrices[2];

  double **MatrixM21 = LocMatrices[d];
  double **MatrixM22 = LocMatrices[d + 1];
  double **MatrixM23 = d == 2 ? nullptr : LocMatrices[d + 2];

  double **MatrixM31 = d == 2 ? nullptr : LocMatrices[2 * d];
  double **MatrixM32 = d == 2 ? nullptr : LocMatrices[2 * d + 1];
  double **MatrixM33 = d == 2 ? nullptr : LocMatrices[2 * d + 2];

  double **MatrixQM1 = LocMatrices[d * d];
  double **MatrixQM2 = LocMatrices[d * d + 1];
  double **MatrixQM3 = d == 2 ? nullptr : LocMatrices[d * d + 2];

  int N_V = N_BaseFuncts[0];
  int N_Q = N_BaseFuncts[1];

  const double *v = OrigValues[0];

  const double *v_x = OrigValues[2];
  const double *v_y = OrigValues[3];
  const double *v_z = d == 2 ? nullptr : OrigValues[4];

  const double *q_x = OrigValues[d + 2];
  const double *q_y = OrigValues[d + 3];
  const double *q_z = d == 2 ? nullptr : OrigValues[d + 4];

  const double *u = param + 0;

  double res_c;
  double res_m[d];

  double tau_m, tau_c;

  NS_GetVMSResiduals<d>(coeff, param, hK, inv_jacobian_offset, res_m, res_c,
    tau_m, tau_c, settings, true);

  // velocity rows
  for (int i = 0; i < N_V; i++)
  {
    // test function: v e_m
    //
    // ->
    //
    // v e_m
    // + \tau_m (u \cdot \nabla) v e_m
    // + \tau_m u_m \nabla v

    double test = v[i];
    double test_x = v_x[i];
    double test_y = v_y[i];
    double test_z = d == 2 ? 0.0 : v_z[i];

    if (residual_terms)
    {
      double conv = u[0] * test_x + u[1] * test_y;
      if (d == 3)
      {
        conv += u[2] * test_z;
      }

      test += conv * tau_m;
    }

    for (int j = 0; j < N_V; j++)
    {
      // ansatz function: v' e_k

      double ansatz = Mult * v[j];

      // (v e_m + \tau_m (u \cdot \nabla) v e_m) \cdot (v' e_k)
      MatrixM11[i][j] += test * ansatz;
      MatrixM22[i][j] += test * ansatz;
      if (d == 3)
      {
        MatrixM33[i][j] += test * ansatz;
      }

      if (residual_terms)
      {
        // \tau_m (u_m \nabla v) \cdot (v' e_k)
        if (d == 2)
        {
          MatrixM11[i][j] += tau_m * u[0] * test_x * ansatz;
          MatrixM12[i][j] += tau_m * u[0] * test_y * ansatz;

          MatrixM21[i][j] += tau_m * u[1] * test_x * ansatz;
          MatrixM22[i][j] += tau_m * u[1] * test_y * ansatz;
        }
        else
        {
          MatrixM11[i][j] += tau_m * u[0] * test_x * ansatz;
          MatrixM12[i][j] += tau_m * u[0] * test_y * ansatz;
          MatrixM13[i][j] += tau_m * u[0] * test_z * ansatz;

          MatrixM21[i][j] += tau_m * u[1] * test_x * ansatz;
          MatrixM22[i][j] += tau_m * u[1] * test_y * ansatz;
          MatrixM23[i][j] += tau_m * u[1] * test_z * ansatz;

          MatrixM31[i][j] += tau_m * u[2] * test_x * ansatz;
          MatrixM32[i][j] += tau_m * u[2] * test_y * ansatz;
          MatrixM33[i][j] += tau_m * u[2] * test_z * ansatz;
        }
      }

      if (residual_terms)
      {
#ifdef RBVMS_SGS_ALTERNATE_FORM
        double prod = res_m[0] * test_x + res_m[1] * test_y;
        if (d == 3)
        {
          prod += res_m[2] * test_z;
        }
        prod *= tau_m;

        // -\tau_m ((v' e_k) \otimes res_m') \cdot (\nabla (v e_m))
        if (d == 2)
        {
          MatrixM11[i][j] -= ansatz * prod;
          MatrixM22[i][j] -= ansatz * prod;
        }
        else
        {
          MatrixM11[i][j] -= ansatz * prod;
          MatrixM22[i][j] -= ansatz * prod;
          MatrixM33[i][j] -= ansatz * prod;
        }
#else
        double ansatz_d = tau_m * ansatz;

        // -\tau_m ((res_m' \otimes (v' e_k)) \cdot (\nabla (v e_m))
        if (d == 2)
        {
          MatrixM11[i][j] -= res_m[0] * ansatz_d * test_x;
          MatrixM12[i][j] -= res_m[0] * ansatz_d * test_y;

          MatrixM21[i][j] -= res_m[1] * ansatz_d * test_x;
          MatrixM22[i][j] -= res_m[1] * ansatz_d * test_y;
        }
        else
        {
          MatrixM11[i][j] -= res_m[0] * ansatz_d * test_x;
          MatrixM12[i][j] -= res_m[0] * ansatz_d * test_y;
          MatrixM13[i][j] -= res_m[0] * ansatz_d * test_z;

          MatrixM21[i][j] -= res_m[1] * ansatz_d * test_x;
          MatrixM22[i][j] -= res_m[1] * ansatz_d * test_y;
          MatrixM23[i][j] -= res_m[1] * ansatz_d * test_z;

          MatrixM31[i][j] -= res_m[2] * ansatz_d * test_x;
          MatrixM32[i][j] -= res_m[2] * ansatz_d * test_y;
          MatrixM33[i][j] -= res_m[2] * ansatz_d * test_z;
        }
#endif
      }
    }
  }

  if (residual_terms && has_B_blocks)
  {
    // pressure rows
    for (int i = 0; i < N_Q; i++)
    {
      // test function: \nabla q

      double test_x = q_x[i];
      double test_y = q_y[i];
      double test_z = d == 2 ? 0.0 : q_z[i];

      for (int j = 0; j < N_V; j++)
      {
        // ansatz function: v e_k

        double ansatz = Mult * v[j];

        // \tau_m (v e_k) \dot (\nabla q)
        MatrixQM1[i][j] += tau_m * ansatz * test_x;
        MatrixQM2[i][j] += tau_m * ansatz * test_y;

        if (d == 3)
        {
          MatrixQM3[i][j] += tau_m * ansatz * test_z;
        }
      }
    }
  }
}

template<int d>
void NS_VMSResiduals_ContinuityTerms_RHS(double Mult, const double *coeff,
  const double *param, const double **OrigValues, const int *N_BaseFuncts,
  double **LocRhs, double tau_m, double tau_c, RBVMS_Settings settings)
{
  const bool has_ut_terms = settings.explicit_time_derivative;
  const bool has_B_blocks = settings.momentum_pressure_coupling_B;

  double *Rhs_u1 = LocRhs[0];
  double *Rhs_u2 = LocRhs[1];
  double *Rhs_u3 = d == 2 ? nullptr : LocRhs[2];
  double *Rhs_p = LocRhs[d];

  int N_V = N_BaseFuncts[0];
  int N_Q = N_BaseFuncts[1];

  const double *q_x = OrigValues[d + 2];
  const double *q_y = OrigValues[d + 3];
  const double *q_z = d == 2 ? nullptr : OrigValues[d + 4];

  const double *v_x = OrigValues[2];
  const double *v_y = OrigValues[3];
  const double *v_z = d == 2 ? nullptr : OrigValues[4];

  if (has_B_blocks)
  {
    const double *f = coeff + 1;

    double rhs_1 = f[0];
    double rhs_2 = f[1];
    double rhs_3 = d == 2 ? 0.0 : f[2];

    if (has_ut_terms)
    {
      const double *u_old = param + (3 * d + 1);
      double inv_tau = 1.0 / TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
      rhs_1 += inv_tau * u_old[0];
      rhs_2 += inv_tau * u_old[1];
      if (d == 3)
      {
        rhs_3 += inv_tau * u_old[2];
      }
    }

    double m = Mult * tau_m;

    // \tau_m (f, \nabla q)
    //
    // - q is a pressure test function
    for (int i = 0; i < N_Q; i++)
    {
      // test term:
      //
      // \nabla q

      double test_x = m * q_x[i];
      double test_y = m * q_y[i];
      double test_z = d == 2 ? 0.0 : (m * q_z[i]);

      Rhs_p[i] += test_x * rhs_1 + test_y * rhs_2 + test_z * rhs_3;
    }
  }

  double gm = Mult * tau_c * coeff[d + 1];

  // velocity rows:
  //
  // \tau_c g \nabla \cdot (v e_m)
  for (int i = 0; i < N_V; i++)
  {
    Rhs_u1[i] += gm * v_x[i];
    Rhs_u2[i] += gm * v_y[i];
    if (d == 3)
    {
      Rhs_u3[i] += gm * v_z[i];
    }
  }
}

template<int d>
void NS_VMSResiduals_CrossTerms_RHS(double Mult, const double *coeff,
  const double *param, const double **OrigValues, const int *N_BaseFuncts,
  double **LocRhs, double tau_m, RBVMS_Settings settings)
{
  const bool has_ut_terms = settings.explicit_time_derivative;

  double *Rhs_u1 = LocRhs[0];
  double *Rhs_u2 = LocRhs[1];
  double *Rhs_u3 = d == 2 ? nullptr : LocRhs[2];

  int N_V = N_BaseFuncts[0];

  const double *v_x = OrigValues[2];
  const double *v_y = OrigValues[3];
  const double *v_z = d == 2 ? nullptr : OrigValues[4];

  const double *u = param + 0;
  const double *f = coeff + 1;

  double m = Mult * tau_m;

  double rhs_1 = f[0];
  double rhs_2 = f[1];
  double rhs_3 = d == 2 ? 0.0 : f[2];

  if (has_ut_terms)
  {
    const double *u_old = param + (3 * d + 1);
    double inv_tau = 1.0 / TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    rhs_1 += inv_tau * u_old[0];
    rhs_2 += inv_tau * u_old[1];
    if (d == 3)
    {
      rhs_3 += inv_tau * u_old[2];
    }
  }

  // \tau_m (f, (u \cdot \nabla) v)
  //
  // - u is the previous velocity
  // - v is a velocity test function
  for (int i = 0; i < N_V; i++)
  {
    // test term:
    // (u \cdot \nabla) v e_m

    double test_conv = u[0] * v_x[i] + u[1] * v_y[i];
    if (d == 3)
    {
      test_conv += u[2] * v_z[i];
    }
    test_conv *= m;

    Rhs_u1[i] += test_conv * rhs_1;
    Rhs_u2[i] += test_conv * rhs_2;

    if (d == 3)
    {
      Rhs_u3[i] += test_conv * rhs_3;
    }
  }

  // \tau_m (f, (\nabla v)^T u)
  //
  // - u is the previous velocity
  // - v is a velocity test function
  for (int i = 0; i < N_V; i++)
  {
    // test term:
    //
    // (\nabla (v e_m))^T u = (St d/dx_t v E_tm) u
    // = u_m \nabla v

    double test_x = m * v_x[i];
    double test_y = m * v_y[i];
    double test_z = d == 2 ? 0.0 : (m * v_z[i]);

    double prod = rhs_1 * test_x + rhs_2 * test_y + rhs_3 * test_z;

    Rhs_u1[i] += u[0] * prod;
    Rhs_u2[i] += u[1] * prod;

    if (d == 3)
    {
      Rhs_u3[i] += u[2] * prod;
    }
  }
}

template<int d>
void NS_VMSResiduals_SubgridScaleTerm_RHS(double Mult, const double *coeff,
  const double *param, const double **OrigValues, const int *N_BaseFuncts,
  double **LocRhs, const double *res_m, double tau_m, RBVMS_Settings settings)
{
  const bool has_ut_terms = settings.explicit_time_derivative;

  double *Rhs_u1 = LocRhs[0];
  double *Rhs_u2 = LocRhs[1];
  double *Rhs_u3 = d == 2 ? nullptr : LocRhs[2];

  int N_V = N_BaseFuncts[0];

  const double *v_x = OrigValues[2];
  const double *v_y = OrigValues[3];
  const double *v_z = d == 2 ? nullptr : OrigValues[4];

  const double *f = coeff + 1;

  double m = Mult * tau_m;

  double rhs_1 = f[0];
  double rhs_2 = f[1];
  double rhs_3 = d == 2 ? 0.0 : f[2];

  if (has_ut_terms)
  {
    const double *u_old = param + (3 * d + 1);
    double inv_tau = 1.0 / TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    rhs_1 += inv_tau * u_old[0];
    rhs_2 += inv_tau * u_old[1];
    if (d == 3)
    {
      rhs_3 += inv_tau * u_old[2];
    }
  }

#ifdef RBVMS_SGS_ALTERNATE_FORM
  // -\tau_m (f \otimes res_m', \nabla v)
  //
  // - res_m' is fully extrapolated
  // - v is a velocity test function
  for (int i = 0; i < N_V; i++)
  {
    // test term:
    //
    // \nabla (v e_m) = St d/dx_t v E_mt

    double prod = res_m[0] * v_x[i] + res_m[1] * v_y[i];
    if (d == 3)
    {
      prod += res_m[2] * v_z[i];
    }
    prod *= m;

    Rhs_u1[i] -= rhs_1 * prod;
    Rhs_u2[i] -= rhs_2 * prod;

    if (d == 3)
    {
      Rhs_u3[i] -= rhs_3 * prod;
    }
  }
#else
  // -\tau_m (res_m' \otimes f, \nabla v)
  //
  // - res_m' is fully extrapolated
  // - v is a velocity test function
  for (int i = 0; i < N_V; i++)
  {
    // test term:
    //
    // \nabla (v e_m) = St d/dx_t v E_mt

    double prod = rhs_1 * v_x[i] + rhs_2 * v_y[i];
    if (d == 3)
    {
      prod += rhs_3 * v_z[i];
    }
    prod *= m;

    Rhs_u1[i] -= res_m[0] * prod;
    Rhs_u2[i] -= res_m[1] * prod;

    if (d == 3)
    {
      Rhs_u3[i] -= res_m[2] * prod;
    }
  }

#endif
}

template<int d>
void NS_VMSResiduals_ContinuityTerms(double Mult, const double *coeff,
  const double *param, const double **OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double tau_m, double tau_c, RBVMS_Settings settings)
{
  const bool has_ut_terms = settings.explicit_time_derivative;
  const bool has_B_blocks = settings.momentum_pressure_coupling_B;
  const bool has_C_block = settings.momentum_pressure_coupling_C;

  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double **MatrixA13 = d == 2 ? nullptr : LocMatrices[2];
  double **MatrixA21 = LocMatrices[d];
  double **MatrixA22 = LocMatrices[d + 1];
  double **MatrixA23 = d == 2 ? nullptr : LocMatrices[d + 2];
  double **MatrixA31 = d == 2 ? nullptr : LocMatrices[2 * d];
  double **MatrixA32 = d == 2 ? nullptr : LocMatrices[2 * d + 1];
  double **MatrixA33 = d == 2 ? nullptr : LocMatrices[2 * d + 2];

  double **MatrixB1 = LocMatrices[d * d + 1];
  double **MatrixB2 = LocMatrices[d * d + 2];
  double **MatrixB3 = d == 2 ? nullptr : LocMatrices[d * d + 3];

  double **MatrixC = LocMatrices[d * d];

  int N_V = N_BaseFuncts[0];
  int N_Q = N_BaseFuncts[1];

  const double *v = OrigValues[0];

  const double *v_x = OrigValues[2];
  const double *v_y = OrigValues[3];
  const double *v_z = d == 2 ? nullptr : OrigValues[4];

  const double *v_xx = OrigValues[2 * d + 2];
  const double *v_yy = OrigValues[2 * d + 3];
  const double *v_zz = d == 2 ? nullptr : OrigValues[2 * d + 4];

  const double *q_x = OrigValues[d + 2];
  const double *q_y = OrigValues[d + 3];
  const double *q_z = d == 2 ? nullptr : OrigValues[d + 4];

  const double *u = param + 0;

  double nu = coeff[0];

  double m = Mult * tau_m;

  double inv_tau = 1.0 / TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  if (has_B_blocks || has_C_block)
  {
    // pressure rows:
    //
    // (res_m, \nabla q)
    //
    // - res_m = \tau_m (\partial_t u - \nu \Delta u + (u \cdot \nabla) u + \nabla p - f)
    //   is an ansatz term (with linearized convection and extrapolated time derivative)
    // - q is a pressure test function
    // - the \partial_t u and f terms are moved to the mass matrix and rhs
    for (int i = 0; i < N_Q; i++)
    {
      // test term:
      //
      // \nabla q

      double test_x = m * q_x[i];
      double test_y = m * q_y[i];
      double test_z = d == 2 ? 0.0 : (m * q_z[i]);

      if (has_B_blocks)
      {
        // velocity columns
        for (int j = 0; j < N_V; j++)
        {
          // ansatz term:
          // (-\nu \Delta + (u \cdot \nabla)) v' e_k

          double ansatz_lap = v_xx[j] + v_yy[j];
          double ansatz_conv = u[0] * v_x[j] + u[1] * v_y[j];
          if (d == 3)
          {
            ansatz_lap += v_zz[j];
            ansatz_conv += u[2] * v_z[j];
          }

          double ansatz_term = -nu * ansatz_lap + ansatz_conv;

          if (has_ut_terms)
          {
            ansatz_term += inv_tau * v[j];
          }

          MatrixB1[i][j] += ansatz_term * test_x;
          MatrixB2[i][j] += ansatz_term * test_y;

          if (d == 3)
          {
            MatrixB3[i][j] += ansatz_term * test_z;
          }
        }
      }

      if (has_C_block)
      {
        // pressure columns
        for (int j = 0; j < N_Q; j++)
        {
          // ansatz term:
          // \nabla q'

          if (d == 2)
          {
            MatrixC[i][j] += q_x[j] * test_x + q_y[j] * test_y;
          }
          else
          {
            MatrixC[i][j] += q_x[j] * test_x + q_y[j] * test_y  + q_z[j] * test_z;
          }
        }
      }
    }
  }

  m = Mult * tau_c;

  // velocity rows:
  //
  // (res_c, \nabla \cdot v)
  //
  // - res_c = tau_c (\nabla \cdot u - g)
  //   is an ansatz term
  // - v is a velocity test function
  // - the g term is moved to the rhs
  for (int i = 0; i < N_V; i++)
  {
    // test term: div(v e_m) = d/dx_m v

    double test_x = m * v_x[i];
    double test_y = m * v_y[i];
    double test_z = d == 2 ? 0.0 : (m * v_z[i]);

    for (int j = 0; j < N_V; j++)
    {
      // ansatz term: div(v' e_k) = d/dx_k v'

      double ansatz_x = v_x[j];
      double ansatz_y = v_y[j];
      double ansatz_z = d == 2 ? 0.0 : v_z[j];

      if (d == 2)
      {
        MatrixA11[i][j] += test_x * ansatz_x;
        MatrixA12[i][j] += test_x * ansatz_y;

        MatrixA21[i][j] += test_y * ansatz_x;
        MatrixA22[i][j] += test_y * ansatz_y;
      }
      else
      {
        MatrixA11[i][j] += test_x * ansatz_x;
        MatrixA12[i][j] += test_x * ansatz_y;
        MatrixA13[i][j] += test_x * ansatz_z;

        MatrixA21[i][j] += test_y * ansatz_x;
        MatrixA22[i][j] += test_y * ansatz_y;
        MatrixA23[i][j] += test_y * ansatz_z;

        MatrixA31[i][j] += test_z * ansatz_x;
        MatrixA32[i][j] += test_z * ansatz_y;
        MatrixA33[i][j] += test_z * ansatz_z;
      }
    }
  }
}

template<int d>
void NS_VMSResiduals_CrossTerms(double Mult, const double *coeff,
  const double *param, const double **OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double tau_m, RBVMS_Settings settings)
{
  const bool has_ut_terms = settings.explicit_time_derivative;

  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double **MatrixA13 = d == 2 ? nullptr : LocMatrices[2];
  double **MatrixA21 = LocMatrices[d];
  double **MatrixA22 = LocMatrices[d + 1];
  double **MatrixA23 = d == 2 ? nullptr : LocMatrices[d + 2];
  double **MatrixA31 = d == 2 ? nullptr : LocMatrices[2 * d];
  double **MatrixA32 = d == 2 ? nullptr : LocMatrices[2 * d + 1];
  double **MatrixA33 = d == 2 ? nullptr : LocMatrices[2 * d + 2];

  double **MatrixB1T = LocMatrices[d * d + d + 1];
  double **MatrixB2T = LocMatrices[d * d + d + 2];
  double **MatrixB3T = d == 2 ? nullptr : LocMatrices[d * d + d + 3];

  int N_V = N_BaseFuncts[0];
  int N_Q = N_BaseFuncts[1];

  const double *v = OrigValues[0];

  const double *v_x = OrigValues[2];
  const double *v_y = OrigValues[3];
  const double *v_z = d == 2 ? nullptr : OrigValues[4];

  const double *v_xx = OrigValues[2 * d + 2];
  const double *v_yy = OrigValues[2 * d + 3];
  const double *v_zz = d == 2 ? nullptr : OrigValues[2 * d + 4];

  const double *q_x = OrigValues[d + 2];
  const double *q_y = OrigValues[d + 3];
  const double *q_z = d == 2 ? nullptr : OrigValues[d + 4];

  const double *u = param + 0;

  double nu = coeff[0];
  double m = Mult * tau_m;

  double inv_tau = 1.0 / TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  // -n(res_m; u, v) = (res_m, (u \cdot \nabla) v)
  //
  // - u is the previous velocity
  // - res_m is ansatz (with linearized convection)
  // - v is a velocity test function
  // - the \partial_t u and f terms are moved to the mass matrix and rhs
  for (int i = 0; i < N_V; i++)
  {
    // test term:
    // (u \cdot \nabla) v e_m

    double test_conv = u[0] * v_x[i] + u[1] * v_y[i];
    if (d == 3)
    {
      test_conv += u[2] * v_z[i];
    }
    test_conv *= m;

    // velocity columns
    for (int j = 0; j < N_V; j++)
    {
      // ansatz term:
      //
      // (-\nu \Delta + (u \cdot \nabla)) v' e_k

      double ansatz_lap = v_xx[j] + v_yy[j];
      double ansatz_conv = u[0] * v_x[j] + u[1] * v_y[j];
      if (d == 3)
      {
        ansatz_lap += v_zz[j];
        ansatz_conv += u[2] * v_z[j];
      }

      double ansatz_term = -nu * ansatz_lap + ansatz_conv;

      if (has_ut_terms)
      {
        ansatz_term += inv_tau * v[j];
      }

      if (d == 2)
      {
        MatrixA11[i][j] += ansatz_term * test_conv;
        MatrixA22[i][j] += ansatz_term * test_conv;
      }
      else
      {
        MatrixA11[i][j] += ansatz_term * test_conv;
        MatrixA22[i][j] += ansatz_term * test_conv;
        MatrixA33[i][j] += ansatz_term * test_conv;
      }
    }

    // pressure columns
    for (int j = 0; j < N_Q; j++)
    {
      // ansatz term:
      //
      // \nabla q'

      MatrixB1T[i][j] += q_x[j] * test_conv;
      MatrixB2T[i][j] += q_y[j] * test_conv;

      if (d == 3)
      {
        MatrixB3T[i][j] += q_z[j] * test_conv;
      }
    }
  }

  // -n(u; res_m, v) = (res_m, (\nabla v)^T u)
  //
  // - u is the previous velocity
  // - res_m is ansatz (with linearized convection)
  // - v is a velocity test function
  // - the \partial_t u and f terms are moved to the mass matrix and rhs
  for (int i = 0; i < N_V; i++)
  {
    // test term:
    //
    // (\nabla (v e_m))^T u = (St d/dx_t v E_tm) u
    // = u_m \nabla v

    double test_x = m * v_x[i];
    double test_y = m * v_y[i];
    double test_z = d == 2 ? 0.0 : (m * v_z[i]);

    // velocity columns
    for (int j = 0; j < N_V; j++)
    {
      // ansatz term:
      //
      // (-\nu \Delta + (u \cdot \nabla)) v' e_k

      double ansatz_lap = v_xx[j] + v_yy[j];
      double ansatz_conv = u[0] * v_x[j] + u[1] * v_y[j];
      if (d == 3)
      {
        ansatz_lap += v_zz[j];
        ansatz_conv += u[2] * v_z[j];
      }

      double ansatz_term = -nu * ansatz_lap + ansatz_conv;

      if (has_ut_terms)
      {
        ansatz_term += inv_tau * v[j];
      }

      if (d == 2)
      {
        MatrixA11[i][j] += ansatz_term * u[0] * test_x;
        MatrixA12[i][j] += ansatz_term * u[0] * test_y;

        MatrixA21[i][j] += ansatz_term * u[1] * test_x;
        MatrixA22[i][j] += ansatz_term * u[1] * test_y;
      }
      else
      {
        MatrixA11[i][j] += ansatz_term * u[0] * test_x;
        MatrixA12[i][j] += ansatz_term * u[0] * test_y;
        MatrixA13[i][j] += ansatz_term * u[0] * test_z;

        MatrixA21[i][j] += ansatz_term * u[1] * test_x;
        MatrixA22[i][j] += ansatz_term * u[1] * test_y;
        MatrixA23[i][j] += ansatz_term * u[1] * test_z;

        MatrixA31[i][j] += ansatz_term * u[2] * test_x;
        MatrixA32[i][j] += ansatz_term * u[2] * test_y;
        MatrixA33[i][j] += ansatz_term * u[2] * test_z;
      }
    }

    // pressure columns
    for (int j = 0; j < N_Q; j++)
    {
      // ansatz term:
      //
      // \nabla q'

      double prod = q_x[j] * test_x + q_y[j] * test_y;
      if (d == 3)
      {
        prod += q_z[j] * test_z;
      }

      MatrixB1T[i][j] += u[0] * prod;
      MatrixB2T[i][j] += u[1] * prod;

      if (d == 3)
      {
        MatrixB3T[i][j] += u[2] * prod;
      }
    }
  }
}

template<int d>
void NS_VMSResiduals_SubgridScaleTerm(double Mult, const double *coeff,
  const double *param, const double **OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, const double *res_m, double tau_m,
  RBVMS_Settings settings)
{
  const bool has_ut_terms = settings.explicit_time_derivative;

  double **MatrixA11 = LocMatrices[0];
  double **MatrixA22 = LocMatrices[d + 1];
  double **MatrixA33 = d == 2 ? nullptr : LocMatrices[2 * d + 2];

#ifndef RBVMS_SGS_ALTERNATE_FORM
  double **MatrixA12 = LocMatrices[1];
  double **MatrixA13 = d == 2 ? nullptr : LocMatrices[2];

  double **MatrixA21 = LocMatrices[d];
  double **MatrixA23 = d == 2 ? nullptr : LocMatrices[d + 2];

  double **MatrixA31 = d == 2 ? nullptr : LocMatrices[2 * d];
  double **MatrixA32 = d == 2 ? nullptr : LocMatrices[2 * d + 1];
#endif

  double **MatrixB1T = LocMatrices[d * d + d + 1];
  double **MatrixB2T = LocMatrices[d * d + d + 2];
  double **MatrixB3T = d == 2 ? nullptr : LocMatrices[d * d + d + 3];

  int N_V = N_BaseFuncts[0];
  int N_Q = N_BaseFuncts[1];

  const double *v = OrigValues[0];

  const double *v_x = OrigValues[2];
  const double *v_y = OrigValues[3];
  const double *v_z = d == 2 ? nullptr : OrigValues[4];

  const double *v_xx = OrigValues[2 * d + 2];
  const double *v_yy = OrigValues[2 * d + 3];
  const double *v_zz = d == 2 ? nullptr : OrigValues[2 * d + 4];

  const double *q_x = OrigValues[d + 2];
  const double *q_y = OrigValues[d + 3];
  const double *q_z = d == 2 ? nullptr : OrigValues[d + 4];

  const double *u = param + 0;

  double nu = coeff[0];

  double m = Mult * tau_m;

  double inv_tau = 1.0 / TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

#ifdef RBVMS_SGS_ALTERNATE_FORM
  // n(res_m'; res_m, v) = -(res_m \otimes res_m', \nabla v)
  //
  // - res_m' is fully extrapolated
  // - res_m is ansatz (with linearized convection)
  // - v is a velocity test function
  // - the \partial_t u and f terms are moved to the mass matrix and rhs
  for (int i = 0; i < N_V; i++)
  {
    // test term:
    //
    // \nabla (v e_m) = St d/dx_t v E_mt

    double test_x = v_x[i];
    double test_y = v_y[i];
    double test_z = d == 2 ? 0.0 : v_z[i];

    double prod = res_m[0] * test_x + res_m[1] * test_y;
    if (d == 3)
    {
      prod += res_m[2] * test_z;
    }
    prod *= m;

    // velocity columns
    for (int j = 0; j < N_V; j++)
    {
      // ansatz term:
      //
      // ( (-\nu \Delta + (u \cdot \nabla)) v' e_k ) \otimes res_m'

      double ansatz_lap = v_xx[j] + v_yy[j];
      double ansatz_conv = u[0] * v_x[j] + u[1] * v_y[j];
      if (d == 3)
      {
        ansatz_lap += v_zz[j];
        ansatz_conv += u[2] * v_z[j];
      }

      double ansatz_term = -nu * ansatz_lap + ansatz_conv;

      if (has_ut_terms)
      {
        ansatz_term += inv_tau * v[j];
      }

      if (d == 2)
      {
        MatrixA11[i][j] -= prod * ansatz_term;
        MatrixA22[i][j] -= prod * ansatz_term;
      }
      else
      {
        MatrixA11[i][j] -= prod * ansatz_term;
        MatrixA22[i][j] -= prod * ansatz_term;
        MatrixA33[i][j] -= prod * ansatz_term;
      }
    }

    // pressure columns
    for (int j = 0; j < N_Q; j++)
    {
      // ansatz term:
      //
      // ( \nabla q' ) \otimes res_m'

      MatrixB1T[i][j] -= prod * q_x[j];
      MatrixB2T[i][j] -= prod * q_y[j];

      if (d == 3)
      {
        MatrixB3T[i][j] -= prod * q_z[j];
      }
    }
  }
#else
  // n(res_m'; res_m, v) = -(res_m' \otimes res_m, \nabla v)
  //
  // - res_m' is fully extrapolated
  // - res_m is ansatz (with linearized convection)
  // - v is a velocity test function
  // - the \partial_t u and f terms are moved to the mass matrix and rhs
  for (int i = 0; i < N_V; i++)
  {
    // test term:
    //
    // \nabla (v e_m) = St d/dx_t v E_mt

    double test_x = m * v_x[i];
    double test_y = m * v_y[i];
    double test_z = d == 2 ? 0.0 : (m * v_z[i]);

    // velocity columns
    for (int j = 0; j < N_V; j++)
    {
      // ansatz term:
      //
      // res_m' \otimes ( (-\nu \Delta + (u \cdot \nabla)) v' e_k )

      double ansatz_lap = v_xx[j] + v_yy[j];
      double ansatz_conv = u[0] * v_x[j] + u[1] * v_y[j];
      if (d == 3)
      {
        ansatz_lap += v_zz[j];
        ansatz_conv += u[2] * v_z[j];
      }

      double ansatz_term = -nu * ansatz_lap + ansatz_conv;

      if (has_ut_terms)
      {
        ansatz_term += inv_tau * v[j];
      }

      if (d == 2)
      {
        MatrixA11[i][j] -= res_m[0] * ansatz_term * test_x;
        MatrixA12[i][j] -= res_m[0] * ansatz_term * test_y;

        MatrixA21[i][j] -= res_m[1] * ansatz_term * test_x;
        MatrixA22[i][j] -= res_m[1] * ansatz_term * test_y;
      }
      else
      {
        MatrixA11[i][j] -= res_m[0] * ansatz_term * test_x;
        MatrixA12[i][j] -= res_m[0] * ansatz_term * test_y;
        MatrixA13[i][j] -= res_m[0] * ansatz_term * test_z;

        MatrixA21[i][j] -= res_m[1] * ansatz_term * test_x;
        MatrixA22[i][j] -= res_m[1] * ansatz_term * test_y;
        MatrixA23[i][j] -= res_m[1] * ansatz_term * test_z;

        MatrixA31[i][j] -= res_m[2] * ansatz_term * test_x;
        MatrixA32[i][j] -= res_m[2] * ansatz_term * test_y;
        MatrixA33[i][j] -= res_m[2] * ansatz_term * test_z;
      }
    }

    // pressure columns
    for (int j = 0; j < N_Q; j++)
    {
      // ansatz term:
      //
      // res_m' \otimes ( \nabla q' )

      double prod = q_x[j] * test_x + q_y[j] * test_y;
      if (d == 3)
      {
        prod += q_z[j] * test_z;
      }

      MatrixB1T[i][j] -= res_m[0] * prod;
      MatrixB2T[i][j] -= res_m[1] * prod;

      if (d == 3)
      {
        MatrixB3T[i][j] -= res_m[2] * prod;
      }
    }
  }
#endif
}

template<int d>
void NS_VMSResiduals_GalerkinContinuityBlocks(double Mult,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices)
{
  double **MatrixB1 = LocMatrices[d * d + 1];
  double **MatrixB2 = LocMatrices[d * d + 2];
  double **MatrixB3 = d == 2 ? nullptr : LocMatrices[d * d + 3];

  double **MatrixB1T = LocMatrices[d * d + d + 1];
  double **MatrixB2T = LocMatrices[d * d + d + 2];
  double **MatrixB3T = d == 2 ? nullptr : LocMatrices[d * d + d + 3];

  int N_V = N_BaseFuncts[0];
  int N_Q = N_BaseFuncts[1];

  const double *q = OrigValues[1];

  const double *v_x = OrigValues[2];
  const double *v_y = OrigValues[3];
  const double *v_z = d == 2 ? nullptr : OrigValues[4];

  // velocity rows:
  // -(\nabla \cdot v, p)
  for (int i = 0; i < N_V; i++)
  {
    double test_x = v_x[i];
    double test_y = v_y[i];
    double test_z = d == 2 ? 0.0 : v_z[i];

    for (int j = 0; j < N_Q; j++)
    {
      double ansatz = Mult * q[j];

      MatrixB1T[i][j] -= ansatz * test_x;
      MatrixB2T[i][j] -= ansatz * test_y;

      if (d == 3)
      {
        MatrixB3T[i][j] -= ansatz * test_z;
      }
    }
  }

  // pressure rows:
  // (\nabla \cdot u, q)
  for (int i = 0; i < N_Q; i++)
  {
    double test = Mult * q[i];

    for (int j = 0; j < N_V; j++)
    {
      MatrixB1[i][j] += test * v_x[j];
      MatrixB2[i][j] += test * v_y[j];

      if (d == 3)
      {
        MatrixB3[i][j] += test * v_z[j];
      }
    }
  }
}

template<int d>
void NSVMSResiduals_RHS(double Mult, const double *coeff, const double *param,
  double hK, const double **OrigValues, const int *N_BaseFuncts,
  double ***, double **LocRhs, int inv_jacobian_offset,
  RBVMS_Settings settings)
{
  double res_c;
  double res_m[d];

  double tau_m, tau_c;

  NS_GetVMSResiduals<d>(coeff, param, hK, inv_jacobian_offset, res_m, res_c,
    tau_m, tau_c, settings, true);

  NS_VMSResiduals_ContinuityTerms_RHS<d>(Mult,
    coeff, param,
    OrigValues, N_BaseFuncts,
    LocRhs, tau_m, tau_c,
    settings);

  NS_VMSResiduals_CrossTerms_RHS<d>(Mult,
    coeff, param,
    OrigValues, N_BaseFuncts,
    LocRhs, tau_m,
    settings);

  NS_VMSResiduals_SubgridScaleTerm_RHS<d>(Mult,
    coeff, param,
    OrigValues, N_BaseFuncts,
    LocRhs, res_m, tau_m,
    settings);
}

template<int d>
void NSVMSResiduals(double Mult, const double *coeff, const double *param,
  double hK, const double **OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **, int inv_jacobian_offset,
  RBVMS_Settings settings)
{
  double res_c;
  double res_m[d];

  double tau_m, tau_c;

  NS_GetVMSResiduals<d>(coeff, param, hK, inv_jacobian_offset, res_m, res_c,
    tau_m, tau_c, settings, true);

  NS_VMSResiduals_ContinuityTerms<d>(Mult,
    coeff, param,
    OrigValues, N_BaseFuncts,
    LocMatrices, tau_m, tau_c,
    settings);

  NS_VMSResiduals_CrossTerms<d>(Mult,
    coeff, param,
    OrigValues, N_BaseFuncts,
    LocMatrices, tau_m,
    settings);

  NS_VMSResiduals_SubgridScaleTerm<d>(Mult,
    coeff, param,
    OrigValues, N_BaseFuncts,
    LocMatrices, res_m, tau_m,
    settings);

  NS_VMSResiduals_GalerkinContinuityBlocks<d>(Mult, OrigValues, N_BaseFuncts,
    LocMatrices);
}

template<int d>
void TNSE_RBVMS_Time_Convection(const double *u,
  const double *u_prime_old, const double* u_prime,
  const double *f, double tau_m, double delta_time,
  double* b, double* c, double* c_prime,
  double& mU, double& mR, double& mR_prime)
{
  mU = tau_m / (tau_m + delta_time);
  mR = delta_time * mU;

  double mU_prime = -1.0 / (tau_m + delta_time); // -mU / tau_m
  mR_prime = 1.0 - delta_time / (tau_m + delta_time); // 1 - mR / tau_m

  for (int i = 0; i < d; i++)
  {
    b[i] = u[i] + u_prime[i];
    c[i] = mU * u_prime_old[i] + mR * f[i];
    c_prime[i] = mU_prime * u_prime_old[i] + mR_prime * f[i];
  }
}

template<int d>
void TNSE_RBVMS_Time_MassMatrix(double Mult, const double *coeff, const double *param,
  double hK, const double **OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **, int inv_jacobian_offset,
  RBVMS_Settings settings)
{
  double **MatrixM11 = LocMatrices[0];
  double **MatrixM12 = LocMatrices[1];
  double **MatrixM13 = d == 2 ? nullptr : LocMatrices[2];

  double **MatrixM21 = LocMatrices[d];
  double **MatrixM22 = LocMatrices[d + 1];
  double **MatrixM23 = d == 2 ? nullptr : LocMatrices[d + 2];

  double **MatrixM31 = d == 2 ? nullptr : LocMatrices[2 * d];
  double **MatrixM32 = d == 2 ? nullptr : LocMatrices[2 * d + 1];
  double **MatrixM33 = d == 2 ? nullptr : LocMatrices[2 * d + 2];

  double **MatrixQM1 = LocMatrices[d * d];
  double **MatrixQM2 = LocMatrices[d * d + 1];
  double **MatrixQM3 = d == 2 ? nullptr : LocMatrices[d * d + 2];

  int N_V = N_BaseFuncts[0];
  int N_Q = N_BaseFuncts[1];

  const double *v = OrigValues[0];

  const double *v_x = OrigValues[2];
  const double *v_y = OrigValues[3];
  const double *v_z = d == 2 ? nullptr : OrigValues[4];

  const double *v_xx = OrigValues[2 * d + 2];
  const double *v_yy = OrigValues[2 * d + 3];
  const double *v_zz = d == 2 ? nullptr : OrigValues[2 * d + 4];

  double nu = coeff[0];

  const double *q_x = OrigValues[d + 2];
  const double *q_y = OrigValues[d + 3];
  const double *q_z = d == 2 ? nullptr : OrigValues[d + 4];

  double res_c;
  double res_m[d];

  double tau_m, tau_c;

  NS_GetVMSResiduals<d>(coeff, param, hK, inv_jacobian_offset, res_m, res_c,
    tau_m, tau_c, settings, true);

  for (int i = 0; i < N_V; i++)
  {
    double test = Mult * v[i];

    // test function: v e_m
    for (int j = 0; j < N_V; j++)
    {
      // ansatz function: v' e_k
      double ansatz = v[j];

      // (v e_m) \cdot (v' e_k)
      MatrixM11[i][j] += test * ansatz;
      MatrixM22[i][j] += test * ansatz;
      if (d == 3)
      {
        MatrixM33[i][j] += test * ansatz;
      }
    }
  }

  if (!settings.explicit_time_derivative)
  {
    const double* old_data = PointwiseAssemblyData::GetOldData();
    double* current_data = PointwiseAssemblyData::GetCurrentData();

    const double* old_subgrid_velocity = old_data;
    const double* old_residual = old_data + d;

    double* current_subgrid_velocity = current_data;
    double* current_residual = current_data + d;

    for (int i = 0; i < d; i++)
    {
      current_residual[i] = res_m[i];
    }

    RBMVS_Time_EvolveSubscale<d>(old_subgrid_velocity, current_subgrid_velocity,
      old_residual, current_residual,
      tau_m, TDatabase::TimeDB->CURRENTTIMESTEPLENGTH,
      settings);

    double mU, mR, mR_prime;
    double b[d];
    double c[d];
    double c_prime[d];

    TNSE_RBVMS_Time_Convection<d>(param + 0,
      old_subgrid_velocity, current_subgrid_velocity, coeff + 1,
      tau_m, TDatabase::TimeDB->CURRENTTIMESTEPLENGTH,
      b, c, c_prime, mU, mR, mR_prime);

    const double *Du = param + 4 * d + 1;

    for (int i = 0; i < N_V; i++)
    {
      double test = Mult * v[i];
      double test_L;

      if (d == 2)
      {
        test_L = Mult * (b[0] * v_x[i] + b[1] * v_y[i]
          - nu * (v_xx[i] + v_yy[i]));
      }
      else
      {
        test_L = Mult * (b[0] * v_x[i] + b[1] * v_y[i] + b[2] * v_z[i]
          - nu * (v_xx[i] + v_yy[i] + v_zz[i]));
      }

      for (int j = 0; j < N_V; j++)
      {
        // -mR' (v, w)
        //
        // M_mk_ij
        // = -mR' (vj e_k, vi e_m)
        // = -mR' \delta_mk vj vi

        if (d == 2)
        {
          double val = mR_prime * test * v[j];
          MatrixM11[i][j] -= val;
          MatrixM22[i][j] -= val;
        }
        else
        {
          double val = mR_prime * test * v[j];
          MatrixM11[i][j] -= val;
          MatrixM22[i][j] -= val;
          MatrixM33[i][j] -= val;
        }

        // -mR ((w \cdot \nabla) u, v)
        //
        // M_mk_ij
        // = -mR (((vj e_k) \cdot \nabla) u , vi e_m)
        // = -mR (vj d_k u, vi e_m)
        // = -mR vj vi d_k um

        double test_ansatz = mR * test * v[j];

        if (d == 2)
        {
          MatrixM11[i][j] -= test_ansatz * Du[0 * d + 0];
          MatrixM12[i][j] -= test_ansatz * Du[0 * d + 1];

          MatrixM21[i][j] -= test_ansatz * Du[1 * d + 0];
          MatrixM22[i][j] -= test_ansatz * Du[1 * d + 1];
        }
        else
        {
          MatrixM11[i][j] -= test_ansatz * Du[0 * d + 0];
          MatrixM12[i][j] -= test_ansatz * Du[0 * d + 1];
          MatrixM13[i][j] -= test_ansatz * Du[0 * d + 2];

          MatrixM21[i][j] -= test_ansatz * Du[1 * d + 0];
          MatrixM22[i][j] -= test_ansatz * Du[1 * d + 1];
          MatrixM23[i][j] -= test_ansatz * Du[1 * d + 2];

          MatrixM31[i][j] -= test_ansatz * Du[2 * d + 0];
          MatrixM32[i][j] -= test_ansatz * Du[2 * d + 1];
          MatrixM33[i][j] -= test_ansatz * Du[2 * d + 2];
        }

        // mR (v, Lw)
        //
        // M_mk_ij
        // = mR (vj e_k, L (vi e_m))
        // = mR (vj e_k, (L vi) e_m)
        // = mR \delta_mk vj (L vi)

        double test_L_ansatz = mR * test_L * v[j];

        if (d == 2)
        {
          MatrixM11[i][j] += test_L_ansatz;
          MatrixM22[i][j] += test_L_ansatz;
        }
        else
        {
          MatrixM11[i][j] += test_L_ansatz;
          MatrixM22[i][j] += test_L_ansatz;
          MatrixM33[i][j] += test_L_ansatz;
        }
      }
    }

    double m = Mult * mR;

    // m_R (v, \nabla q)
    for (int i = 0; i < N_Q; i++)
    {
      // test function: \nabla q

      double test_x = m * q_x[i];
      double test_y = m * q_y[i];
      double test_z = d == 2 ? 0.0 : (m * q_z[i]);

      for (int j = 0; j < N_V; j++)
      {
        // ansatz function: v e_k

        double ansatz = v[j];

        // m_R (v e_k, \nabla q)
        // = m_R v d_k q
        MatrixQM1[i][j] += ansatz * test_x;
        MatrixQM2[i][j] += ansatz * test_y;

        if (d == 3)
        {
          MatrixQM3[i][j] += ansatz * test_z;
        }
      }
    }
  }
}

template<int d>
void TNSE_RBVMS_Time_Pressure_RHS(double Mult, const double* u_prime,
  const double **OrigValues, const int *N_BaseFuncts, double **LocRhs)
{
  double *Rhs_p = LocRhs[d];
  int N_Q = N_BaseFuncts[1];

  const double *q_x = OrigValues[d + 2];
  const double *q_y = OrigValues[d + 3];
  const double *q_z = d == 2 ? nullptr : OrigValues[d + 4];

  double rhs_1 = u_prime[0];
  double rhs_2 = u_prime[1];
  double rhs_3 = d == 2 ? 0.0 : u_prime[2];

  // (u', \nabla q)
  //
  // - q is a pressure test function
  for (int i = 0; i < N_Q; i++)
  {
    // test term:
    //
    // \nabla q

    double test_x = q_x[i];
    double test_y = q_y[i];
    double test_z = d == 2 ? 0.0 : q_z[i];

    Rhs_p[i] += Mult * (test_x * rhs_1 + test_y * rhs_2 + test_z * rhs_3);
  }
}

template<int d>
void TNSE_RBVMS_Time_SGS_RHS(double Mult,
  const double* u_prime, const double* u_prime_old, const double* coeff,
  const double **OrigValues, const int *N_BaseFuncts, double **LocRhs)
{
  double *Rhs_u1 = LocRhs[0];
  double *Rhs_u2 = LocRhs[1];
  double *Rhs_u3 = d == 2 ? nullptr : LocRhs[2];

  int N_V = N_BaseFuncts[0];

  const double *v = OrigValues[0];

  const double *v_x = OrigValues[2];
  const double *v_y = OrigValues[3];
  const double *v_z = d == 2 ? nullptr : OrigValues[4];

  const double *v_xx = OrigValues[2 * d + 2];
  const double *v_yy = OrigValues[2 * d + 3];
  const double *v_zz = d == 2 ? nullptr : OrigValues[2 * d + 4];

  double nu = coeff[0];

  double inv_tau = 1.0 / TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double dt_u_prime_1 = inv_tau * (u_prime[0] - u_prime_old[0]);
  double dt_u_prime_2 = inv_tau * (u_prime[1] - u_prime_old[1]);
  double dt_u_prime_3 = d == 2 ? 0.0 : (inv_tau * (u_prime[2] - u_prime_old[2]));

  // -(\partial_t u', v)
  for (int i = 0; i < N_V; i++)
  {
    double test = Mult * v[i];

    Rhs_u1[i] -= dt_u_prime_1 * test;
    Rhs_u2[i] -= dt_u_prime_2 * test;

    if (d == 3)
    {
      Rhs_u3[i] -= dt_u_prime_3 * test;
    }
  }

  // \nu (u', \Delta_h v)
  for (int i = 0; i < N_V; i++)
  {
    // test term:
    //
    // \Delta (v e_m)

    double val = Mult * nu * (v_xx[i] + v_yy[i] + v_zz[i]);

    Rhs_u1[i] += u_prime[0] * val;
    Rhs_u2[i] += u_prime[1] * val;

    if (d == 3)
    {
      Rhs_u3[i] += u_prime[2] * val;
    }
  }

  // n(u'; v, u') = (u', (u' \cdot \nabla) v)
  for (int i = 0; i < N_V; i++)
  {
    // test term:
    //
    // \nabla (v e_m)

    double prod = u_prime[0] * v_x[i] + u_prime[1] * v_y[i];
    if (d == 3)
    {
      prod += u_prime[2] * v_z[i];
    }
    prod *= Mult;

    Rhs_u1[i] += prod * u_prime[0];
    Rhs_u2[i] += prod * u_prime[1];
    if (d == 3)
    {
      Rhs_u3[i] += prod * u_prime[2];
    }
  }
}

template<int d>
void TNSE_RBVMS_Time_GradDiv_RHS(double Mult, double tau_c, const double* coeff,
  const double **OrigValues, const int *N_BaseFuncts, double **LocRhs)
{
  double *Rhs_u1 = LocRhs[0];
  double *Rhs_u2 = LocRhs[1];
  double *Rhs_u3 = d == 2 ? nullptr : LocRhs[2];

  int N_V = N_BaseFuncts[0];

  const double *v_x = OrigValues[2];
  const double *v_y = OrigValues[3];
  const double *v_z = d == 2 ? nullptr : OrigValues[4];

  double gm = Mult * tau_c * coeff[d + 1];

  // \tau_c g \nabla \cdot (v e_m)
  for (int i = 0; i < N_V; i++)
  {
    Rhs_u1[i] += gm * v_x[i];
    Rhs_u2[i] += gm * v_y[i];
    if (d == 3)
    {
      Rhs_u3[i] += gm * v_z[i];
    }
  }
}

template<int d>
void TNSE_RBVMS_Time_ImplicitResiduals_RHS(double Mult,
  const double* b, const double* c, const double* c_prime,
  const double* coeff, const double* param,
  const double **OrigValues, const int *N_BaseFuncts, double **LocRhs)
{
  double *Rhs_u1 = LocRhs[0];
  double *Rhs_u2 = LocRhs[1];
  double *Rhs_u3 = d == 2 ? nullptr : LocRhs[2];
  double *Rhs_p = LocRhs[d];

  int N_V = N_BaseFuncts[0];
  int N_Q = N_BaseFuncts[1];

  const double *v = OrigValues[0];

  const double *v_x = OrigValues[2];
  const double *v_y = OrigValues[3];
  const double *v_z = d == 2 ? nullptr : OrigValues[4];

  const double *v_xx = OrigValues[2 * d + 2];
  const double *v_yy = OrigValues[2 * d + 3];
  const double *v_zz = d == 2 ? nullptr : OrigValues[2 * d + 4];

  const double *q_x = OrigValues[d + 2];
  const double *q_y = OrigValues[d + 3];
  const double *q_z = d == 2 ? nullptr : OrigValues[d + 4];

  const double *Du = param + 4 * d + 1;

  double nu = coeff[0];

  // -(c', v)
  for (int i = 0; i < N_V; i++)
  {
    double test = Mult * v[i];

    if (d == 2)
    {
      Rhs_u1[i] -= test * c_prime[0];
      Rhs_u2[i] -= test * c_prime[1];
    }
    else
    {
      Rhs_u1[i] -= test * c_prime[0];
      Rhs_u2[i] -= test * c_prime[1];
      Rhs_u3[i] -= test * c_prime[2];
    }
  }

  for (int i = 0; i < N_V; i++)
  {
    // -((c \cdot \nabla) u, v)
    //
    // rhs_m_i
    // = -((c \cdot \nabla) u, vi e_m )
    // = -vi (c \cdot \nabla) um
    if (d == 2)
    {
      double test = Mult * v[i];
      Rhs_u1[i] -= test * (c[0] * Du[0 * d + 0] + c[1] * Du[0 * d + 1]);
      Rhs_u2[i] -= test * (c[0] * Du[1 * d + 0] + c[1] * Du[1 * d + 1]);
    }
    else
    {
      double test = Mult * v[i];
      Rhs_u1[i] -= test
        * (c[0] * Du[0 * d + 0] + c[1] * Du[0 * d + 1] + c[2] * Du[0 * d + 2]);
      Rhs_u2[i] -= test
        * (c[0] * Du[1 * d + 0] + c[1] * Du[1 * d + 1] + c[2] * Du[1 * d + 2]);
      Rhs_u3[i] -= test
        * (c[0] * Du[2 * d + 0] + c[1] * Du[2 * d + 1] + c[2] * Du[2 * d + 2]);
    }

    // (c, Lv)
    //
    // rhs_m_i
    // = (c, L (vi e_m))
    // = cm (L vi)
    if (d == 2)
    {
      double test_L = Mult * (b[0] * v_x[i] + b[1] * v_y[i]
        - nu * (v_xx[i] + v_yy[i]));
      Rhs_u1[i] += c[0] * test_L;
      Rhs_u2[i] += c[1] * test_L;
    }
    else
    {
      double test_L = Mult
        * (b[0] * v_x[i] + b[1] * v_y[i] + b[2] * v_z[i]
          - nu * (v_xx[i] + v_yy[i] + v_zz[i]));
      Rhs_u1[i] += c[0] * test_L;
      Rhs_u2[i] += c[1] * test_L;
      Rhs_u3[i] += c[2] * test_L;
    }
  }

  for (int i = 0; i < N_Q; i++)
  {
    // (c, \nabla q)

    if (d == 2)
    {
      Rhs_p[i] += Mult * (q_x[i] * c[0] + q_y[i] * c[1]);
    }
    else
    {
      Rhs_p[i] += Mult * (q_x[i] * c[0] + q_y[i] * c[1] + q_z[i] * c[2]);
    }
  }
}

template<int d>
void TNSE_RBVMS_Time_SubgridTerms_RHS(double Mult, const double *coeff,
  const double *param, double hK, const double **OrigValues,
  const int *N_BaseFuncts, double ***, double **LocRhs,
  int inv_jacobian_offset, RBVMS_Settings settings)
{
  const double* old_data = PointwiseAssemblyData::GetOldData();
  double* current_data = PointwiseAssemblyData::GetCurrentData();

  const double* old_subgrid_velocity = old_data;
  const double* old_residual = old_data + d;

  double* current_subgrid_velocity = current_data;
  double* current_residual = current_data + d;

  double res_c;
  double res_m[d];

  double tau_m, tau_c;

  NS_GetVMSResiduals<d>(coeff, param, hK, inv_jacobian_offset, res_m, res_c,
    tau_m, tau_c, settings, false);

  for (int i = 0; i < d; i++)
  {
    current_residual[i] = res_m[i];
  }

  RBMVS_Time_EvolveSubscale<d>(old_subgrid_velocity, current_subgrid_velocity,
    old_residual, current_residual,
    tau_m, TDatabase::TimeDB->CURRENTTIMESTEPLENGTH,
    settings);

  if (!settings.explicit_time_derivative)
  {
    double mU, mR, mR_prime;
    double b[d];
    double c[d];
    double c_prime[d];

    TNSE_RBVMS_Time_Convection<d>(param + 0,
      old_subgrid_velocity, current_subgrid_velocity, coeff + 1,
      tau_m, TDatabase::TimeDB->CURRENTTIMESTEPLENGTH,
      b, c, c_prime, mU, mR, mR_prime);

    TNSE_RBVMS_Time_ImplicitResiduals_RHS<d>(Mult, b, c, c_prime,
      coeff, param, OrigValues, N_BaseFuncts, LocRhs);
  }
  else
  {
    // (u', \nabla q)
    TNSE_RBVMS_Time_Pressure_RHS<d>(Mult, current_subgrid_velocity,
      OrigValues, N_BaseFuncts, LocRhs);

    // -(\partial_t u', v) + \nu (u', \Delta_h v) + (u', (u' \cdot \nabla) v)
    TNSE_RBVMS_Time_SGS_RHS<d>(Mult,
      current_subgrid_velocity, old_subgrid_velocity, coeff,
      OrigValues, N_BaseFuncts, LocRhs);

    // \tau_c (g, \nabla \cdot v)
    TNSE_RBVMS_Time_GradDiv_RHS<d>(Mult, tau_c, coeff,
      OrigValues, N_BaseFuncts, LocRhs);
  }
}

template<int d>
void TNSE_RBVMS_Time_CrossTerms(double Mult, const double* u_prime,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices)
{
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double **MatrixA13 = d == 2 ? nullptr : LocMatrices[2];
  double **MatrixA21 = LocMatrices[d];
  double **MatrixA22 = LocMatrices[d + 1];
  double **MatrixA23 = d == 2 ? nullptr : LocMatrices[d + 2];
  double **MatrixA31 = d == 2 ? nullptr : LocMatrices[2 * d];
  double **MatrixA32 = d == 2 ? nullptr : LocMatrices[2 * d + 1];
  double **MatrixA33 = d == 2 ? nullptr : LocMatrices[2 * d + 2];

  int N_V = N_BaseFuncts[0];

  const double *v = OrigValues[0];

  const double *v_x = OrigValues[2];
  const double *v_y = OrigValues[3];
  const double *v_z = d == 2 ? nullptr : OrigValues[4];

  // ((u' \cdot \nabla) u, v)
  for (int i = 0; i < N_V; i++)
  {
    // test function: v e_m

    double test = Mult * v[i];

    for (int j = 0; j < N_V; j++)
    {
      // ansatz function: v' e_k
      //
      // product: A_mk_ij
      // = ((u' \cdot \nabla) (v' e_k), v e_m)
      // = (((u' \cdot \nabla) v') e_k, v e_m)
      // = \delta_mk v (u' \cdot \nabla) v'

      double val = u_prime[0] * v_x[j] + u_prime[1] * v_y[j];
      if (d == 3)
      {
        val += u_prime[2] * v_z[j];
      }
      val *= test;

      MatrixA11[i][j] += val;
      MatrixA22[i][j] += val;
      if (d == 3)
      {
        MatrixA33[i][j] += val;
      }
    }
  }

  // -(u', (u \cdot \nabla) v)
  for (int i = 0; i < N_V; i++)
  {
    // test function: v e_m

    double test_x = Mult * v_x[i];
    double test_y = Mult * v_y[i];
    double test_z = d == 2 ? 0.0 : (Mult * v_z[i]);

    for (int j = 0; j < N_V; j++)
    {
      // ansatz function: v' e_k
      //
      // product: A_mk_ij
      // = (-u', ((v' e_k) \cdot \nabla) (v e_m) )
      // = (-u', (v' d_k v) e_m)
      // = -u'_m v' d_k v

      double ansatz = v[j];

      if (d == 2)
      {
        MatrixA11[i][j] -= u_prime[0] * ansatz * test_x;
        MatrixA12[i][j] -= u_prime[0] * ansatz * test_y;

        MatrixA21[i][j] -= u_prime[1] * ansatz * test_x;
        MatrixA22[i][j] -= u_prime[1] * ansatz * test_y;
      }
      else
      {
        MatrixA11[i][j] -= u_prime[0] * ansatz * test_x;
        MatrixA12[i][j] -= u_prime[0] * ansatz * test_y;
        MatrixA13[i][j] -= u_prime[0] * ansatz * test_z;

        MatrixA21[i][j] -= u_prime[1] * ansatz * test_x;
        MatrixA22[i][j] -= u_prime[1] * ansatz * test_y;
        MatrixA23[i][j] -= u_prime[1] * ansatz * test_z;

        MatrixA31[i][j] -= u_prime[2] * ansatz * test_x;
        MatrixA32[i][j] -= u_prime[2] * ansatz * test_y;
        MatrixA33[i][j] -= u_prime[2] * ansatz * test_z;
      }
    }
  }
}

template<int d>
void TNSE_RBVMS_Time_GradDiv(double Mult, double tau_c,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices)
{
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double **MatrixA13 = d == 2 ? nullptr : LocMatrices[2];
  double **MatrixA21 = LocMatrices[d];
  double **MatrixA22 = LocMatrices[d + 1];
  double **MatrixA23 = d == 2 ? nullptr : LocMatrices[d + 2];
  double **MatrixA31 = d == 2 ? nullptr : LocMatrices[2 * d];
  double **MatrixA32 = d == 2 ? nullptr : LocMatrices[2 * d + 1];
  double **MatrixA33 = d == 2 ? nullptr : LocMatrices[2 * d + 2];

  int N_V = N_BaseFuncts[0];

  const double *v_x = OrigValues[2];
  const double *v_y = OrigValues[3];
  const double *v_z = d == 2 ? nullptr : OrigValues[4];

  double m = tau_c * Mult;

  // \tau_c (\nabla \cdot u, \nabla \cdot v)
  for (int i = 0; i < N_V; i++)
  {
    // test function: v e_m

    double test_x = m * v_x[i];
    double test_y = m * v_y[i];
    double test_z = d == 2 ? 0.0 : (m * v_z[i]);

    for (int j = 0; j < N_V; j++)
    {
      // ansatz function: v' e_k
      //
      // product: A_mk_ij
      // = -\tau_c (\nabla \cdot (v' e_k), \nabla \cdot (v e_m))
      // = -\tau_c d_k v' d_m v

      double ansatz_x = v_x[j];
      double ansatz_y = v_y[j];
      double ansatz_z = d == 2 ? 0.0 : v_z[j];

      if (d == 2)
      {
        MatrixA11[i][j] += test_x * ansatz_x;
        MatrixA12[i][j] += test_x * ansatz_y;

        MatrixA21[i][j] += test_y * ansatz_x;
        MatrixA22[i][j] += test_y * ansatz_y;
      }
      else
      {
        MatrixA11[i][j] += test_x * ansatz_x;
        MatrixA12[i][j] += test_x * ansatz_y;
        MatrixA13[i][j] += test_x * ansatz_z;

        MatrixA21[i][j] += test_y * ansatz_x;
        MatrixA22[i][j] += test_y * ansatz_y;
        MatrixA23[i][j] += test_y * ansatz_z;

        MatrixA31[i][j] += test_z * ansatz_x;
        MatrixA32[i][j] += test_z * ansatz_y;
        MatrixA33[i][j] += test_z * ansatz_z;
      }
    }
  }
}

template<int d>
void TNSE_RBVMS_Time_ImplicitResiduals(double Mult, double mR, double mR_prime,
  const double *b, const double *coeff, const double *param,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices)
{
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA12 = LocMatrices[1];
  double **MatrixA13 = d == 2 ? nullptr : LocMatrices[2];
  double **MatrixA21 = LocMatrices[d];
  double **MatrixA22 = LocMatrices[d + 1];
  double **MatrixA23 = d == 2 ? nullptr : LocMatrices[d + 2];
  double **MatrixA31 = d == 2 ? nullptr : LocMatrices[2 * d];
  double **MatrixA32 = d == 2 ? nullptr : LocMatrices[2 * d + 1];
  double **MatrixA33 = d == 2 ? nullptr : LocMatrices[2 * d + 2];

  double **MatrixB1 = LocMatrices[d * d + 1];
  double **MatrixB2 = LocMatrices[d * d + 2];
  double **MatrixB3 = d == 2 ? nullptr : LocMatrices[d * d + 3];

  double **MatrixB1T = LocMatrices[d * d + d + 1];
  double **MatrixB2T = LocMatrices[d * d + d + 2];
  double **MatrixB3T = d == 2 ? nullptr : LocMatrices[d * d + d + 3];

  double **MatrixC = LocMatrices[d * d];

  int N_V = N_BaseFuncts[0];
  int N_Q = N_BaseFuncts[1];

  const double *Du = param + 4 * d + 1;

  const double *v = OrigValues[0];

  const double *v_x = OrigValues[2];
  const double *v_y = OrigValues[3];
  const double *v_z = d == 2 ? nullptr : OrigValues[4];

  const double *v_xx = OrigValues[2 * d + 2];
  const double *v_yy = OrigValues[2 * d + 3];
  const double *v_zz = d == 2 ? nullptr : OrigValues[2 * d + 4];

  const double *q_x = OrigValues[d + 2];
  const double *q_y = OrigValues[d + 3];
  const double *q_z = d == 2 ? nullptr : OrigValues[d + 4];

  double nu = coeff[0];

  // velocity columns
  for (int j = 0; j < N_V; j++)
  {
    double ansatz_L;

    if (d == 2)
    {
      ansatz_L = -nu * (v_xx[j] + v_yy[j]) + b[0] * v_x[j] + b[1] * v_y[j];
    }
    else
    {
      ansatz_L = -nu * (v_xx[j] + v_yy[j] + v_zz[j])
        + b[0] * v_x[j] + b[1] * v_y[j] + b[2] * v_z[j];
    }

    ansatz_L *= Mult;

    // velocity rows
    for (int i = 0; i < N_V; i++)
    {
      double test = v[i];
      double test_x = v_x[i];
      double test_y = v_y[i];
      double test_z = d == 2 ? 0.0 : v_z[i];
      double test_L;

      if (d == 2)
      {
        test_L = -nu * (v_xx[i] + v_yy[i]) + b[0] * test_x + b[1] * test_y;
      }
      else
      {
        test_L = -nu * (v_xx[i] + v_yy[i] + v_zz[i])
         + b[0] * test_x + b[1] * test_y + b[2] * test_z;
      }

      // mR ((-Lu \cdot \nabla) u, v)
      //
      // A_mk_ij
      // = mR ((-L(vj e_k) \cdot \nabla) u, vi e_m)
      // = -mR ( (L vj) d_k u, vi e_m )
      // = -mR (L vj) vi d_k um

      if (d == 2)
      {
        double val = mR * ansatz_L * test;

        MatrixA11[i][j] -= val * Du[0 * d + 0];
        MatrixA12[i][j] -= val * Du[0 * d + 1];

        MatrixA21[i][j] -= val * Du[1 * d + 0];
        MatrixA22[i][j] -= val * Du[1 * d + 1];
      }
      else
      {
        double val = mR * ansatz_L * test;

        MatrixA11[i][j] -= val * Du[0 * d + 0];
        MatrixA12[i][j] -= val * Du[0 * d + 1];
        MatrixA13[i][j] -= val * Du[0 * d + 2];

        MatrixA21[i][j] -= val * Du[1 * d + 0];
        MatrixA22[i][j] -= val * Du[1 * d + 1];
        MatrixA23[i][j] -= val * Du[1 * d + 2];

        MatrixA31[i][j] -= val * Du[2 * d + 0];
        MatrixA32[i][j] -= val * Du[2 * d + 1];
        MatrixA33[i][j] -= val * Du[2 * d + 2];
      }

      // mR (Lu, Lv)
      //
      // A_mk_ij
      // = mR (L (vj e_k), L (vi e_m))
      // = mR ((L vj) e_k, (L vi) e_m)
      // = mR \delta_mk (L vj) (L vi)

      if (d == 2)
      {
        double val = mR * ansatz_L * test_L;

        MatrixA11[i][j] += val;
        MatrixA22[i][j] += val;
      }
      else
      {
        double val = mR * ansatz_L * test_L;

        MatrixA11[i][j] += val;
        MatrixA22[i][j] += val;
        MatrixA33[i][j] += val;
      }

      // -mR' (Lu, v)
      //
      // A_mk_ij
      // = -mR' (L (vj e_k), vi e_m)
      // = -mR' \delta_mk vi (L vj)
      if (d == 2)
      {
        double val = mR_prime * ansatz_L * test;

        MatrixA11[i][j] -= val;
        MatrixA22[i][j] -= val;
      }
      else
      {
        double val = mR_prime * ansatz_L * test;

        MatrixA11[i][j] -= val;
        MatrixA22[i][j] -= val;
        MatrixA33[i][j] -= val;
      }
    }

    // pressure rows
    for (int i = 0; i < N_Q; i++)
    {
      // mR (Lu, \nabla q)
      MatrixB1[i][j] += mR * ansatz_L * q_x[j];
      MatrixB2[i][j] += mR * ansatz_L * q_y[j];
      if (d == 3)
      {
        MatrixB3[i][j] += mR * ansatz_L * q_z[j];
      }
    }
  }

  // pressure columns
  for (int j = 0; j < N_Q; j++)
  {
    double ansatz_x = Mult * q_x[j];
    double ansatz_y = Mult * q_y[j];
    double ansatz_z = d == 2 ? 0.0 : (Mult * q_z[j]);

    // velocity rows
    for (int i = 0; i < N_V; i++)
    {
      // mR (-(\nabla p \cdot \nabla) u, v)
      //
      // BT_m_ij
      // = mR (-(\nabla qj \cdot \nabla) u, vi e_m)
      // = -mR vi (\nabla qj \cdot \nabla) um

      if (d == 2)
      {
        double test = mR * v[i];
        double conv_1 = ansatz_x * Du[0 * d + 0] + ansatz_y * Du[0 * d + 1];
        double conv_2 = ansatz_x * Du[1 * d + 0] + ansatz_y * Du[1 * d + 1];

        MatrixB1T[i][j] -= conv_1 * test;
        MatrixB2T[i][j] -= conv_2 * test;
      }
      else
      {
        double test = mR * v[i];
        double conv_1 = ansatz_x * Du[0 * d + 0]
          + ansatz_y * Du[0 * d + 1]
          + ansatz_z * Du[0 * d + 2];

        double conv_2 = ansatz_x * Du[1 * d + 0]
          + ansatz_y * Du[1 * d + 1]
          + ansatz_z * Du[1 * d + 2];

        double conv_3 = ansatz_x * Du[2 * d + 0]
          + ansatz_y * Du[2 * d + 1]
          + ansatz_z * Du[2 * d + 2];

        MatrixB1T[i][j] -= conv_1 * test;
        MatrixB2T[i][j] -= conv_2 * test;
        MatrixB3T[i][j] -= conv_3 * test;
      }

      // mR (\nabla p, L v)
      //
      // BT_m_ij
      // = mR (\nabla qj, L (vi e_m))
      // = mR (\nabla qj, (L v_i) e_m)
      // = mR (d_m qj) (L v_i)
      if (d == 2)
      {
        double test_L = mR * (b[0] * v_x[i] + b[1] * v_y[i]
          - nu * (v_xx[i] + v_yy[i]));

        MatrixB1T[i][j] += ansatz_x * test_L;
        MatrixB2T[i][j] += ansatz_y * test_L;
      }
      else
      {
        double test_L = mR
          * (b[0] * v_x[i] + b[1] * v_y[i] + b[2] * v_z[i]
            - nu * (v_xx[i] + v_yy[i] + v_zz[i]));

        MatrixB1T[i][j] += ansatz_x * test_L;
        MatrixB2T[i][j] += ansatz_y * test_L;
        MatrixB3T[i][j] += ansatz_z * test_L;
      }

      // -mR' (\nabla p, v)
      //
      // BT_m_ij
      // = -mR' (\nabla qj, vi e_m)
      // = vi (d_m qj)
      if (d == 2)
      {
        double test = mR_prime * v[i];
        MatrixB1T[i][j] -= test * ansatz_x;
        MatrixB2T[i][j] -= test * ansatz_y;
      }
      else
      {
        double test = mR_prime * v[i];
        MatrixB1T[i][j] -= test * ansatz_x;
        MatrixB2T[i][j] -= test * ansatz_y;
        MatrixB3T[i][j] -= test * ansatz_z;
      }
    }

    // pressure rows
    for (int i = 0; i < N_Q; i++)
    {
      // mR (\nabla p, \nabla q)
      if (d == 2)
      {
        MatrixC[i][j] += mR * (ansatz_x * q_x[i] + ansatz_y * q_y[i]);
      }
      else
      {
        MatrixC[i][j] += mR
          * (ansatz_x * q_x[i] + ansatz_y * q_y[i] + ansatz_z * q_z[i]);
      }
    }
  }
}

template<int d>
void TNSE_RBVMS_Time_SubgridTerms_Matrices(double Mult, const double *coeff,
  const double *param, double hK, const double **OrigValues,
  const int *N_BaseFuncts, double ***LocMatrices, double **,
  int inv_jacobian_offset, RBVMS_Settings settings)
{
  const double* old_data = PointwiseAssemblyData::GetOldData();
  double* current_data = PointwiseAssemblyData::GetCurrentData();

  const double* old_subgrid_velocity = old_data;
  const double* old_residual = old_data + d;

  double* current_subgrid_velocity = current_data;
  double* current_residual = current_data + d;

  double res_c;
  double res_m[d];

  double tau_m, tau_c;

  NS_GetVMSResiduals<d>(coeff, param, hK, inv_jacobian_offset, res_m, res_c,
    tau_m, tau_c, settings, false);

  for (int i = 0; i < d; i++)
  {
    current_residual[i] = res_m[i];
  }

  RBMVS_Time_EvolveSubscale<d>(old_subgrid_velocity, current_subgrid_velocity,
    old_residual, current_residual,
    tau_m, TDatabase::TimeDB->CURRENTTIMESTEPLENGTH,
    settings);

  if (!settings.explicit_time_derivative)
  {
    double mU, mR, mR_prime;
    double b[d];
    double c[d];
    double c_prime[d];

    TNSE_RBVMS_Time_Convection<d>(param + 0,
      old_subgrid_velocity, current_subgrid_velocity, coeff + 1,
      tau_m, TDatabase::TimeDB->CURRENTTIMESTEPLENGTH,
      b, c, c_prime, mU, mR, mR_prime);

    TNSE_RBVMS_Time_ImplicitResiduals<d>(Mult,
      mR, mR_prime, b, coeff, param,
      OrigValues, N_BaseFuncts, LocMatrices);
  }
  else
  {
    // n(u'; u, v) + n(u; u', v)
    TNSE_RBVMS_Time_CrossTerms<d>(Mult, current_subgrid_velocity,
      OrigValues, N_BaseFuncts, LocMatrices);

    // -\tau_c (\nabla \cdot u, \nabla \cdot v)
    TNSE_RBVMS_Time_GradDiv<d>(Mult, tau_c,
      OrigValues, N_BaseFuncts, LocMatrices);
  }

  NS_VMSResiduals_GalerkinContinuityBlocks<d>(Mult, OrigValues, N_BaseFuncts,
    LocMatrices);
}

template <>
void compute_stabilization_parameters<2>(const double* u, const double* coeff, double* params)
{
  double x0 = TDatabase::ParamDB->INTERNAL_VERTEX_X[0];
  double y0 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[0];
  double x1 = TDatabase::ParamDB->INTERNAL_VERTEX_X[1];
  double y1 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[1];
  double x2 = TDatabase::ParamDB->INTERNAL_VERTEX_X[2];
  double y2 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[2];
 
  double d11, d12, d21, d22;
  // double rec_detjk = 1./coeff[19];
  
  // triangle
  if (TDatabase::ParamDB->INTERNAL_VERTEX_X[3]== -4711)
  {
    double xc1=x1-x0;
    double xc2=x2-x0;
    double yc1=y1-y0;
    double yc2=y2-y0;
    double rec_detjk = 1./(xc1*yc2-xc2*yc1);

    d11 = (y2-y0) * rec_detjk;  //dxi/dx
    d12 = (x0-x2) * rec_detjk;  //dxi/dy
    d21 = (y0-y1) * rec_detjk;  //deta/dx
    d22 = (x1-x0) * rec_detjk;  //deta/dy
  }
  else
  {
    double x3 = TDatabase::ParamDB->INTERNAL_VERTEX_X[3];
    double y3 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[3];
    double xc1=(-x0 + x1 + x2 - x3) * 0.25;
    double xc2=(-x0 - x1 + x2 + x3) * 0.25;
    double yc1=(-y0 + y1 + y2 - y3) * 0.25;
    double yc2=(-y0 - y1 + y2 + y3) * 0.25;
    
    double rec_detjk = 1./(xc1*yc2 - xc2*yc1);
    // quadrilateral
    d11 = yc2 * rec_detjk;  //dxi/dx
    d12 = -xc2 * rec_detjk;  //dxi/dy
    d21 = -yc1 * rec_detjk;  //deta/dx
    d22 = xc1 * rec_detjk;  //deta/dy
  }

  double time_step_length = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double cinv = TDatabase::ParamDB->DELTA0;
  double tau_m, tau_c;

  double g11 = d11*d11 + d21*d21;
  double g12 = d11*d12 + d21*d22;
  double g22 = d12*d12 + d22*d22;

  tau_m = g11*g11 + 2.*g12*g12 + g22*g22; // G:G
  tau_m *= cinv*coeff[0]*coeff[0];
  tau_m +=  4./(time_step_length*time_step_length);
  tau_m += u[0] * (g11*u[0]+g12*u[1]) + u[1]*(g12*u[0]+g22*u[1]);
  tau_m = 1./std::sqrt(tau_m);

  tau_c = (d11+d21)*(d11+d21)+(d12+d22)*(d12+d22);
  tau_c *= tau_m;
  tau_c = 1./tau_c;

  params[0] = tau_m;
  params[1] = tau_c;
  
  // for output 
  // TDatabase::ParamDB->P14 = tau_m;
  // TDatabase::ParamDB->P15 = tau_c;
}

template <>
void compute_stabilization_parameters<3>(const double* u, const double* coeff, double* params)
{
  double eps  = 1e-12;
  double time_step_length = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double cinv = TDatabase::ParamDB->DELTA0;
  
  double x0 = TDatabase::ParamDB->INTERNAL_VERTEX_X[0];
  double y0 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[0];
  double z0 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[0];
  
  double x1 = TDatabase::ParamDB->INTERNAL_VERTEX_X[1];
  double y1 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[1];
  double z1 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[1];
  
  double x2 = TDatabase::ParamDB->INTERNAL_VERTEX_X[2];
  double y2 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[2];
  double z2 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[2];
  
  
  double nu = coeff[0];
  double rec_detjk = coeff[19];
  rec_detjk = 1/rec_detjk;
  double d11, d12, d13, d21, d22, d23, d31, d32, d33;
  // tetrahedron
  if (TDatabase::ParamDB->INTERNAL_VERTEX_X[4] == -4711)
  {
    double x3 = TDatabase::ParamDB->INTERNAL_VERTEX_X[3];
    double y3 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[3];
    double z3 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[3];
    
    d11 = ((y2-y0)*(z3-z0)+(y3-y0)*(z0-z2)) * rec_detjk;  //dxi/dx
    d12 = ((x3-x0)*(z2-z0)+(x2-x0)*(z0-z3)) * rec_detjk;  //dxi/dy
    d13 = ((x2-x0)*(y3-y0)+(x3-x0)*(y0-y2)) * rec_detjk;  //dxi/dz
    
    d21 = ((y3-y0)*(z1-z0)+(y1-y0)*(z0-z3)) * rec_detjk;  //deta/dx
    d22 = ((x1-x0)*(z3-z0)+(x3-x0)*(z0-z1)) * rec_detjk;  //deta/dy
    d23 = ((x3-x0)*(y1-y0)+(x1-x0)*(y0-y3)) * rec_detjk;  //deta/dz
    
    d31 = ((y1-y0)*(z2-z0)+(y2-y0)*(z0-z1)) * rec_detjk;  //dzeta/dx
    d32 = ((x2-x0)*(z1-z0)+(x1-x0)*(z0-z2)) * rec_detjk;  //dzeta/dy
    d33 = ((x1-x0)*(y2-y0)+(x2-x0)*(y0-y1)) * rec_detjk;  //dzeta/dz	
  }
  else // hexahedron
  {
    double x4 = TDatabase::ParamDB->INTERNAL_VERTEX_X[4];
    double y4 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[4];
    double z4 = TDatabase::ParamDB->INTERNAL_VERTEX_Z[4];
    
    d11 = ((y1-y0)*(z4-z0)+(y0-y4)*(z1-z0)) * 0.5 * rec_detjk;  //dxi/dx
    d12 = ((x4-x0)*(z1-z0)+(x1-x0)*(z0-z4)) * 0.5 * rec_detjk;  //dxi/dy
    d13 = ((x1-x0)*(y4-y0)+(x0-x4)*(y1-y0)) * 0.5 * rec_detjk;  //dxi/dz
    
    d21 = ((y4-y0)*(z1-z2)+(y1-y2)*(z0-z4)) * 0.5 * rec_detjk;  //deta/dx
    d22 = ((x1-x2)*(z4-z0)+(x0-x4)*(z1-z2)) * 0.5 * rec_detjk;  //deta/dy
    d23 = ((x4-x0)*(y1-y2)+(x1-x2)*(y0-y4)) * 0.5 * rec_detjk;  //deta/dz
    
    d31 = ((y1-y2)*(z1-z0)+(y1-y0)*(z2-z1)) * 0.5 * rec_detjk;  //dzeta/dx
    d32 = ((x1-x0)*(z1-z2)+(x1-x2)*(z0-z1)) * 0.5 * rec_detjk;  //dzeta/dy
    d33 = ((x1-x2)*(y1-y0)+(x1-x0)*(y2-y1)) * 0.5 * rec_detjk;  //dzeta/dz
  }
  
  double g11 = d11*d11 + d21*d21 + d31*d31;
  double g12 = d11*d12 + d21*d22 + d31*d32;
  double g13 = d11*d13 + d21*d23 + d31*d33;
  double g22 = d12*d12 + d22*d22 + d32*d32;
  double g23 = d12*d13 + d22*d23 + d32*d33;
  double g33 = d13*d13 + d23*d23 + d33*d33;
  
  // G : G
  double tau_m = g11*g11 + 2*g12*g12 + 2*g13*g13 + g22*g22 + 2*g23*g23 + g33*g33; 
  
  tau_m *= cinv*nu*nu;
  tau_m +=  4/(time_step_length*time_step_length); 
  tau_m += u[0] * (g11*u[0]+g12*u[1]+g13*u[2]) + u[1]*(g12*u[0]+g22*u[1]+g23*u[2])
           + u[2]*(g13*u[0]+g23*u[1]+g33*u[2]);
  if (tau_m < eps)
  {
    params[0] = 0;
    params[1] = 0;
    return;
  }
  
  // parameter for the mementum equation
  tau_m = 1./std::sqrt(tau_m);
  double tau_c;
  tau_c  = (d11+d21+d31)*(d11+d21+d31);
  tau_c += (d12+d22+d32)*(d12+d22+d32);
  tau_c += (d13+d23+d33)*(d13+d23+d33);
  
  tau_c *= tau_m;
  tau_c = 1./tau_c;
  
  params[0] = tau_m;
  params[1] = tau_c;
}
#ifdef __2D__
template void NSMassMatrixSingle<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSMassMatrix<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSLaplaceDeformationVariationalMS<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSLumpMassMatrix<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSVariationlMS_GMatrices<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSVariationlMS_GTildeMatrices<2>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NS_SUPG<2>(double Mult, const double *coeff, const double *param, 
  double, const double **OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **, double delta0, double delta1);
template void NS_SUPG_GradientBlocks<2>(
  double Mult, const double *, const double *param, 
  double hK, const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices, 
  double **, double delta0);
template void NS_SUPG_RightHandSide_InfSup<2>(double Mult, const double *coeff, 
  const double *param, double hK, const double **OrigValues, const int *N_BaseFuncts, double ***, 
  double **LocRhs, double delta0);
template void NS_SUPG_MassMatrix<2>(double Mult, const double *coeff, const double *param, 
  double hK, const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices, 
  double **, double delta0);
template void NS_SUPG_EquOrder_Gradient_DivergenceBlocks<2>(double Mult, const double *coeff, 
  const double *param, double, const double **OrigValues, const int *N_BaseFuncts, 
  double ***LocMatrices, double **, double factor);
template void NS_SUPG_RightHandSide_EquOrder<2>(double Mult, const double *coeff, 
  const double *param, double hK, const double **OrigValues, const int *N_BaseFuncts, double ***, 
  double **LocRhs, double factor);
template void NS_SUPG_skew_symmetric<2>(double Mult, const double *coeff, const double *param, 
  double, const double **OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **, double delta0, double delta1);
template void NS_SUPG_rotational<2>(double Mult, const double *coeff, const double *param, 
  double, const double **OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **, double delta0, double delta1);
template void NS_SUPG_emac<2>(double Mult, const double *coeff, const double *param, 
  double, const double **OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **, double delta0, double delta1);
template void NSParamsVelocityDerivatives_SUPG_inf_sup<2>(const double *in, double *out);
template void NSParamsVelocityDerivatives_SUPG_equal_order<2>(const double *in, double *out);

template void NSParamsVMSResidualsWithoutLaplacianAndRHS<2>(const double *in, double *out, bool extend_advection);
template void NSParamsVMSResidualsLaplacian<2>(const double *in, double *out);
template void NSParamsOldVelocity<2>(const double *in, double *out);

template void NSVMSResiduals_MassMatrix<2>(double Mult, const double *coeff, const double *param,
  double hK, const double**OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **LocRhs, int inv_jacobian_offset,
  RBVMS_Settings settings);
template void NSVMSResiduals<2>(double Mult, const double *coeff, const double *param,
  double hK, const double**OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **LocRhs, int inv_jacobian_offset,
  RBVMS_Settings settings);
template void NSVMSResiduals_RHS<2>(double Mult, const double *coeff, const double *param,
  double hK, const double**OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **LocRhs, int inv_jacobian_offset,
  RBVMS_Settings settings);

template void TNSE_RBVMS_Time_SubgridTerms_Matrices<2>(double Mult, const double *coeff,
  const double *param, double hK, const double **OrigValues,
  const int *N_BaseFuncts, double ***LocMatrices, double **,
  int inv_jacobian_offset, RBVMS_Settings settings);
template void TNSE_RBVMS_Time_SubgridTerms_RHS<2>(double Mult, const double *coeff,
  const double *param, double hK, const double **OrigValues,
  const int *N_BaseFuncts, double ***LocMatrices, double **,
  int inv_jacobian_offset, RBVMS_Settings settings);
template void TNSE_RBVMS_Time_MassMatrix<2>(double Mult, const double *coeff, const double *param,
  double hK, const double **OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **, int inv_jacobian_offset,
  RBVMS_Settings settings);

#endif
#ifdef __3D__
template void NSMassMatrixSingle<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSMassMatrix<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSLaplaceGradGradSingleSmagorinsky<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSLaplaceGradGradSmagorinsky<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSLaplaceDeformationSmagorinsky<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSLaplaceDeformationVariationalMS<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSLumpMassMatrix<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSVariationlMS_GMatrices<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NSVariationlMS_GTildeMatrices<3>(
  double Mult, const double *coeff, const double *param, double hK,
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices,
  double **LocRhs);
template void NS_SUPG<3>(double Mult, const double *coeff, const double *param, 
  double, const double **OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **, 
  double delta0, double tau_c);
template void NS_SUPG_GradientBlocks<3>(double Mult, const double *, const double *param, 
  double hK, const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices, 
  double **, double delta0);
template void NS_SUPG_RightHandSide_InfSup<3>(double Mult, const double *coeff, 
  const double *param, double hK, const double **OrigValues, const int *N_BaseFuncts, double ***, 
  double **LocRhs, double delta0);
template void NS_SUPG_MassMatrix<3>(double Mult, const double *coeff, const double *param, 
  double hK, const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices, 
  double **, double delta0);
template void NS_SUPG_EquOrder_Gradient_DivergenceBlocks<3>(double Mult, const double *coeff, 
  const double *param, double, const double **OrigValues, const int *N_BaseFuncts, 
  double ***LocMatrices, double **, double factor);
template void NS_SUPG_RightHandSide_EquOrder<3>(double Mult, const double *coeff, 
  const double *param, double hK, const double **OrigValues, const int *N_BaseFuncts, double ***, 
  double **LocRhs, double factor);
template void NS_SUPG_skew_symmetric<3>(double Mult, const double *coeff, const double *param, 
  double, const double **OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **, double delta0, double delta1);
template void NS_SUPG_rotational<3>(double Mult, const double *coeff, const double *param, 
  double, const double **OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **, double delta0, double delta1);
template void NS_SUPG_emac<3>(double Mult, const double *coeff, const double *param, 
  double, const double **OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **, double delta0, double delta1);
template void NSParamsVelocityDerivatives_SUPG_inf_sup<3>(const double *in, double *out);
template void NSParamsVelocityDerivatives_SUPG_equal_order<3>(const double *in, double *out);

template void NSParamsVMSResidualsWithoutLaplacianAndRHS<3>(const double *in, double *out, bool extend_advection);
template void NSParamsVMSResidualsLaplacian<3>(const double *in, double *out);
template void NSParamsOldVelocity<3>(const double *in, double *out);

template void NSVMSResiduals_MassMatrix<3>(double Mult, const double *coeff, const double *param,
  double hK, const double**OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **LocRhs, int inv_jacobian_offset,
  RBVMS_Settings settings);
template void NSVMSResiduals<3>(double Mult, const double *coeff, const double *param,
  double hK, const double**OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **LocRhs, int inv_jacobian_offset,
  RBVMS_Settings settings);
template void NSVMSResiduals_RHS<3>(double Mult, const double *coeff, const double *param,
  double hK, const double**OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **LocRhs, int inv_jacobian_offset,
  RBVMS_Settings settings);

template void TNSE_RBVMS_Time_SubgridTerms_Matrices<3>(double Mult, const double *coeff,
  const double *param, double hK, const double **OrigValues,
  const int *N_BaseFuncts, double ***LocMatrices, double **,
  int inv_jacobian_offset, RBVMS_Settings settings);
template void TNSE_RBVMS_Time_SubgridTerms_RHS<3>(double Mult, const double *coeff,
  const double *param, double hK, const double **OrigValues,
  const int *N_BaseFuncts, double ***LocMatrices, double **,
  int inv_jacobian_offset, RBVMS_Settings settings);
template void TNSE_RBVMS_Time_MassMatrix<3>(double Mult, const double *coeff, const double *param,
  double hK, const double **OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **, int inv_jacobian_offset,
  RBVMS_Settings settings);

#endif
