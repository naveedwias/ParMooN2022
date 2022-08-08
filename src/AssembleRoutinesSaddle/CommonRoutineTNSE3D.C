#include <CommonRoutineTNSE3D.h>

#include <Database.h>
#include <Convolution.h>
#include <ConvDiff.h>
#include "MooNMD_Io.h"
#include <cmath>
#include <string.h>

extern "C"
{
  void dgesvd_(char *jobu, char *jobvt, int *m, int *n, double *a, int *lda, double *s, double *u, int *ldu, double *vt, int *ldvt, double *work, int *lwork, int *info);
}

template<int r>
bool svdvals_(double* A, double* s)
{
  int rk = r;
  int lwork = 5 * r;
  int info = 0;

  double work[5 * r];

  dgesvd_("N", "N", &rk, &rk, A, &rk, s, nullptr, &rk, nullptr, &rk, work, &lwork, &info);

  if (info != 0)
  {
    Output::warn("LAPACK", "dgesvd returned info = ", info);
  }

  return info == 0;
}

template<int r>
bool svdvals(const double* A, double* s)
{
  double AA[r * r];

  memcpy(AA, A, r * r * sizeof(double));

  return svdvals_<r>(AA, s);
}

double standard(double delta, double frobeniusNorm)
{
  double viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;

  double viscosity = viscosity_constant * delta * delta
    * std::sqrt(frobeniusNorm);

  return viscosity;
}

double laplacian(double delta, double frobeniusNorm)
{
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  double  viscosity_power = TDatabase::ParamDB->TURBULENT_VISCOSITY_POWER;

  double viscosity = viscosity_constant * delta * delta * 
                      std::pow(frobeniusNorm, viscosity_power/2.0);
  return viscosity;
}

// Layton SIAM J. Sci. Comput. 1996 
double getViscosityLayton96(double delta, double frobeniusNorm)
{
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  double  viscosity_power = TDatabase::ParamDB->TURBULENT_VISCOSITY_POWER;
  double  viscosity_sigma = TDatabase::ParamDB->TURBULENT_VISCOSITY_SIGMA;
  
  double viscosity = viscosity_constant * std::pow(delta,viscosity_sigma) * 
                 std::pow(frobeniusNorm,(viscosity_power-2)/2.0)
                 /std::pow(std::abs(std::log(delta)),2*(viscosity_power-1)/3);
  return viscosity;
}
// \|u-g_\delta \ast u\|_2 
// give better name
double normof_velocityminus_gdeltaTimesuStar(double delta, const double* u,
                                             const double *uConv)
{
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  double viscosity = viscosity_constant * delta * 
              std::sqrt((u[0]-uConv[0])*(u[0]-uConv[0]) +
                        (u[1]-uConv[1])*(u[1]-uConv[1]) +
                        (u[2]-uConv[2])*(u[2]-uConv[2]));
  return viscosity;
}
// \|Du^h-G^H\|_2
// the values of the off diagonal entries of G^H has to be divided by 2
// since the basis functions are l^H/2 !!!
double normof_Defvelocityminus_gdeltaTimesuStar(double delta,
                                                const double* gradu,
                                                const double *uConv)
{
  double a11, a12, a13, a22, a23, a33;
  a11 = gradu[0] - uConv[0];
  a12 = (gradu[1]+gradu[3]-uConv[1])/2.0;
  a13 = (gradu[2]+gradu[6]-uConv[2])/2.0;
  a22 = gradu[4]-uConv[3];
  a23 = (gradu[5]+gradu[7]-uConv[4])/2.0;
  a33 = gradu[8]-uConv[5];
  
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  double viscosity = viscosity_constant * delta * delta * 
              std::sqrt(a11*a11+2*a12*a12+2*a13*a13 
                       +a22*a22+2*a23*a23+a33*a33);
  return viscosity;
}

// all parameters free
double parametersFree(double delta, double frobeniusNorm)
{
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  double  viscosity_power = TDatabase::ParamDB->TURBULENT_VISCOSITY_POWER;
  double viscosity = viscosity_constant * std::pow(delta,viscosity_power) * 
      std::pow(frobeniusNorm,viscosity_power/2.0);
  return viscosity;
}
// walls at z=0, and z = 2

double vanDriestDampingChannel(double delta, double frobeniusNorm,
                               const double* z)
{
  // walls at z=0 and z=2
  double reynolds_number = TDatabase::ParamDB->RE_NR;
  double viscosity;
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  double van_driest_damping = 26.0;
  double zplus = reynolds_number*(1-std::abs(1-z[0]));
  if(zplus < 5)
  {
    viscosity = viscosity_constant * delta * delta * (1-std::exp(-zplus/van_driest_damping))
                * (1-std::exp(-zplus/van_driest_damping)) * std::sqrt(frobeniusNorm);
  }
  else
  {
    viscosity = viscosity_constant * delta * delta * std::sqrt(frobeniusNorm);
  }
  return viscosity;
}

// van Driest damping for cylinder with squared cross--section
// left and right wall at the cylinder
double vanDriestDampingCylinder(double delta, double frobeniusNorm,
                                const double* x, const double* y)
{
  double x0, y0;
  double van_driest_damping = 26.0;
  double eps = 1e-06;
  double zplus = 1000;
  
  if ((x[0] > 0.45 - eps) && (x[0] < 0.55 + eps))
  {
    // distance to the wall
    if (y[0] > 0.7)
      y0 = y[0] - 0.75;
    else
      y0 = 0.65 - y[0];
    // wall units 
    zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_SIDES * y0;
  }
  if ((y[0] > 0.65 - eps) && (y[0] < 0.75 + eps))
  {
    // distance to the wall
    if (x[0] < 0.5)
    {
      x0 = 0.45 - x[0];
      // wall units 
      zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_FRONT * x0;
    }
    else
    {
      x0 = x[0] - 0.55;
      // wall units 
      zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_BACK * x0;
    }
  }
  double viscosity;
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  if(zplus < 5)
  {
    viscosity = viscosity_constant * delta * delta * (1-std::exp(-zplus/van_driest_damping)) *
                (1-std::exp(-zplus/van_driest_damping)) * std::sqrt(frobeniusNorm);
  }
  else
  {
    viscosity = viscosity_constant * delta * delta * std::sqrt(frobeniusNorm);
  }
  return viscosity;
}

// walls at z=0 and z=2
double vanDriestDampingChannelContinuous(double delta, double frobeniusNorm,
                                         const double* z)
{
  double reynolds_number = TDatabase::ParamDB->RE_NR;
  double zplus = reynolds_number*(1-std::abs(1-z[0]));
  double van_driest_damping = 26.0;
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  double viscosity = viscosity_constant * delta * delta * (1-std::exp(-zplus/van_driest_damping)) *
              (1-std::exp(-zplus/van_driest_damping)) * std::sqrt(frobeniusNorm);

  return viscosity;
}

// van Driest damping for channel flow (paper: Rudman, Blackburn'99)
  // walls at z=0 and z=2
double vanDriestDampingChannelRB99(double delta, double frobeniusNorm,
                                   const double* z)
{
  double reynolds_number = TDatabase::ParamDB->RE_NR;
  double zplus = reynolds_number*(1-std::abs(1-z[0]));
  double van_driest_damping = 26.0;
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  
  double viscosity = viscosity_constant * delta * delta * 
              (1-std::exp(-(zplus/van_driest_damping)*(zplus/van_driest_damping) *
              (zplus/van_driest_damping))) * std::sqrt(frobeniusNorm);

  return viscosity;
}

// van Driest damping for channel flow (paper: Rudman, Blackburn'99) with diff A+
double vanDriestDampingChannelRB99APlus(double delta, double frobeniusNorm,
                                        const double* z)
{
  double reynolds_number = TDatabase::ParamDB->RE_NR;
  double van_driest_damping = 17.0;
  // walls at z=0 and z=2
  double zplus = reynolds_number*(1-std::abs(1-z[0]));
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  double viscosity = viscosity_constant * delta * delta * (1-std::exp(-(zplus/van_driest_damping) *
              (zplus/van_driest_damping)*(zplus/van_driest_damping))) * 
              std::sqrt(frobeniusNorm);
  return viscosity;
}

// van Driest damping (continuous, classical) for cylinder with squared cross--section
double vanDriestDampingCylinderContinuous(double delta, double frobeniusNorm,
                                          const double* x, const double* y)
{
  // left and right wall at the cylinder
  double x0, y0;
  double zplus = 1000;
  double eps = 1e-06;
  
  if ((x[0] > 0.45 - eps) && (x[0] < 0.55 + eps))
  {
    // distance to the wall
    if (y[0] > 0.7)
      y0 = y[0] - 0.75;
    else
      y0 = 0.65 - y[0];
    // wall units
    zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_SIDES * y0;
  }
  if ((y[0] > 0.65 - eps) && (y[0] < 0.75 + eps))
  {
    // distance to the wall
    if (x[0] < 0.5)
    {
      x0 = 0.45 - x[0];
      // wall units
      zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_FRONT * x0;
    }
    else
    {
      x0 = x[0] - 0.55;
      // wall units
      zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_BACK * x0;
    }
  }
  double van_driest_damping = 26.0;
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  double viscosity = viscosity_constant * delta * delta * (1-std::exp(-zplus/van_driest_damping)) *
              (1-std::exp(-zplus/van_driest_damping)) * std::sqrt(frobeniusNorm);
  
  return viscosity;
}

// van Driest damping (paper: Rudman, Blackburn'99) for cylinder with 
// squared cross--section
double vanDriestDampingCylinderRB99(double delta, double frobeniusNorm,
                                    const double* x, const double *y)
{
  // left and right wall at the cylinder
  double x0, y0, zplus = 1000;
  double eps = 1e-06;
  
  if ((x[0] > 0.45 - eps) && (x[0] < 0.55 + eps))
  {
    // distance to the wall
    if (y[0] > 0.7)
      y0 = y[0] - 0.75;
    else
      y0 = 0.65 - y[0];
    // wall units
    zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_SIDES * y0;
  }
  if ((y[0] > 0.65 - eps) && (y[0] < 0.75 + eps))
  {
    // distance to the wall
    if (x[0] < 0.5)
    {
      x0 = 0.45 - x[0];
      // wall units
      zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_FRONT * x0;
    }
    else
    {
      x0 = x[0] - 0.55;
      // wall units
      zplus = TDatabase::ParamDB->CYLINDER_22000_YPLUS_BACK * x0;
    }
  }
  double  viscosity_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  double van_driest_damping = 26.0;
  
  double viscosity = viscosity_constant* delta * delta *
             (1-std::exp(-(zplus/van_driest_damping)*(zplus/van_driest_damping) *
             (zplus/van_driest_damping))) *  std::sqrt(frobeniusNorm);
  return viscosity;
}

/// Verstappen model (J Sci Comput'11) 
/// C = 1/mu_max as on page 107
double VerstappenViscositymodelCmuMax(double delta, const double* gradu,
                                      double frobeniusNorm)
{
  double a11, a12, a13, a22, a23, a33;
  // compute filter width in coordinate directions
  // hk is just the value if the filter width would be zero (this should not happen)
  double hk = delta/2.0;
  double delta_x = Mesh_size_in_convection_direction<3>(hk, {{1,0,0}});
  delta_x *= delta_x;
  double delta_y = Mesh_size_in_convection_direction<3>(hk, {{0,1,0}});
  delta_y *= delta_y;
  double delta_z = Mesh_size_in_convection_direction<3>(hk, {{0,0,1}});
  delta_z *= delta_z;
  
  double mu_max = 4 * ( 1./delta_x + 1./delta_y + 1./delta_z );
  if(TDatabase::ParamDB->TURBULENT_VISCOSITY_TENSOR == 0)
  {
    // 0.5*(grad(u) + grad(u)^T)
    a11 = gradu[0]+gradu[0];
    a11 /= 2.;
    a12 = gradu[1]+gradu[3];
    a12 /= 2.;
    a13 = gradu[2]+gradu[6];
    a13 /= 2.;
    a22 = gradu[4]+gradu[4];
    a22 /= 2.;
    a23 = gradu[5]+gradu[7];
    a23 /= 2.;
    a33 = gradu[8]+gradu[8];
    a33 /= 2.;
  }
  else
    ErrThrow("ERROR: Verstappen model needs a symmetric stress tensor!");
  
  double invariant_3 = - (a11*a22*a33 + 2.*a12*a23*a13 - a12*a12*a33 - 
                          a23*a23*a11 - a13*a13*a22);
  double invariant_2 = 0.5 * frobeniusNorm;
      
  double viscosity = (1.5 * std::abs(invariant_3) ) / (mu_max * invariant_2);
  return viscosity;
}

// Verstappen model (J Sci Comput'11) 
// C = (h/pi)^2 as on page as on page 97, where Delta = (h_x*h_y*h_z)^(1/3) 
double VerstappenViscositymodelChPhi(double delta, const double* gradu,
                                     double frobeniusNorm)
{
  double a11, a12, a13, a22, a23, a33;
  // compute filter width in coordinate directions
  // hk is just the value if the filter width would be zero (this should not happen)
  double hk = delta/2.0;
  double delta_x = Mesh_size_in_convection_direction<3>(hk, {{1,0,0}});
  double delta_y = Mesh_size_in_convection_direction<3>(hk, {{0,1,0}});
  double delta_z = Mesh_size_in_convection_direction<3>(hk, {{0,0,1}});
  
  // TODO: change cell width hk using CELL_MEASURE (more elegant), now too slow! 
  
  hk = delta_x*delta_y*delta_z;
  hk = std::pow(hk,1.0/3.0);
  
  if(TDatabase::ParamDB->TURBULENT_VISCOSITY_TENSOR==0)
  {
    // 0.5*(grad(u) + grad(u)^T)
    a11 = gradu[0]+gradu[0];
    a11 /= 2.;
    a12 = gradu[1]+gradu[3];
    a12 /= 2.;
    a13 = gradu[2]+gradu[6];
    a13 /= 2.;
    a22 = gradu[4]+gradu[4];
    a22 /= 2.;
    a23 = gradu[5]+gradu[7];
    a23 /= 2.;
    a33 = gradu[8]+gradu[8];
    a33 /= 2.;
  }
  else
    ErrThrow("ERROR: Verstappen model needs a symmetric stress tensor!");
  
  double invariant_3 = - (a11*a22*a33 + 2.*a12*a23*a13 - a12*a12*a33 -
                          a23*a23*a11 - a13*a13*a22);
  double invariant_2 = 0.5 * frobeniusNorm;
      
  double viscosity = ( 1.5 * hk * hk * std::abs(invariant_3) )  / ( M_PI * M_PI * invariant_2 );  
  return viscosity;
}

double VerstappenModelSimple(double delta, const double* gradu,
                             double frobeniusNorm)
{
  double a11, a12, a13, a22, a23, a33;
  if(TDatabase::ParamDB->TURBULENT_VISCOSITY_TENSOR==0)
  {
    a11 = gradu[0]+gradu[0];
    a12 = gradu[1]+gradu[3];
    a13 = gradu[2]+gradu[6];
    a22 = gradu[4]+gradu[4];
    a23 = gradu[5]+gradu[7];
    a33 = gradu[8]+gradu[8];
  }
  else
    ErrThrow("ERROR: Verstappen model needs a symmetric stress tensor!");
  double invariant_3 = - (a11 * a22 * a33 + 2.* a12 * a23 * a13 -
                          a12 * a12 * a33 - a23 * a23 * a11 - a13 * a13 * a22);
  invariant_3 /= 8.0;
  double invariant_2 = frobeniusNorm;
  double viscosity = (6.0 * delta * delta * std::abs(invariant_3) );
  if (invariant_2>0)
    viscosity /= (M_PI * M_PI * invariant_2);
  else
    viscosity /= (M_PI * M_PI);   
  return viscosity;
}

// eddy viscosity model: Vreman, Phys. Fluids 16 (10), 3670 -3681, 2004
// frobenius norm of gradient of velocity
// use same notations as in paper
double vremanViscosityModel(const double *gradu)
{
  // Alpha <- |Du|^2_F
  double Alpha = 0.0;

  for (int i = 0; i < 9; i++)
  {
    Alpha += gradu[i] * gradu[i];
  }

  if (std::abs(Alpha) < 1.e-12)
  {
    // velocity is locally near-constant - return zero eddy viscosity
    return 0.0;
  }
  else
  {
    // compute filter width in coordinate directions

    double min_x = TDatabase::ParamDB->INTERNAL_VERTEX_X[0];
    double max_x = min_x;
    double min_y = TDatabase::ParamDB->INTERNAL_VERTEX_Y[0];
    double max_y = min_y;
    double min_z = TDatabase::ParamDB->INTERNAL_VERTEX_Z[0];
    double max_z = min_z;

    for (int i = 1; i < 4; i++)
    {
      double x = TDatabase::ParamDB->INTERNAL_VERTEX_X[i];
      double y = TDatabase::ParamDB->INTERNAL_VERTEX_Y[i];
      double z = TDatabase::ParamDB->INTERNAL_VERTEX_Z[i];

      min_x = std::min(x, min_x);
      max_x = std::max(x, max_x);

      min_y = std::min(y, min_y);
      max_y = std::max(y, max_y);

      min_z = std::min(z, min_z);
      max_z = std::max(z, max_z);
    }

    // \delta <- (d_x^2, d_y^2, d_z^2)

    double delta_x = max_x - min_x;
    delta_x *= delta_x;

    double delta_y = max_y - min_y;
    delta_y *= delta_y;

    double delta_z = max_z - min_z;
    delta_z *= delta_z;

    // compute second invariant of gradient of velocity, scaled with
    // filter width in coordinate directions

    // b_ij <- \sum_k \delta_k \partial_i u_k \partial_j u_k
    double b11 = delta_x * gradu[0] * gradu[0] + delta_y * gradu[3] * gradu[3] + delta_z * gradu[6] * gradu[6];
    double b22 = delta_x * gradu[1] * gradu[1] + delta_y * gradu[4] * gradu[4] + delta_z * gradu[7] * gradu[7];
    double b33 = delta_x * gradu[2] * gradu[2] + delta_y * gradu[5] * gradu[5] + delta_z * gradu[8] * gradu[8];
    double b12 = delta_x * gradu[0] * gradu[1] + delta_y * gradu[3] * gradu[4] + delta_z * gradu[6] * gradu[7];
    double b13 = delta_x * gradu[0] * gradu[2] + delta_y * gradu[3] * gradu[5] + delta_z * gradu[6] * gradu[8];
    double b23 = delta_x * gradu[1] * gradu[2] + delta_y * gradu[4] * gradu[5] + delta_z * gradu[7] * gradu[8];

    // Bbeta <- \sum_{i > j} (b_ii b_jj - b_ij^2)
    double Bbeta = b11 * b22 - b12 * b12 + b11 * b33 - b13 * b13 + b22 * b33 - b23 * b23;

    // check for round-off errors
    if (Bbeta < 0.0)
    {
      Bbeta = 0.0;
    }

    return TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT * std::sqrt(Bbeta / Alpha);
  }
}

double sigmaViscosityModel(double delta, const double* gradu)
{
  double sigma[3];
  svdvals<3>(gradu, sigma);

  if (sigma[0] < 1.e-16)
  {
    return 0.0;
  }

  double c = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT * delta / sigma[0];

  return c * c * sigma[2] * (sigma[0] - sigma[1]) * (sigma[1] - sigma[2]);
}

double frobeniusNormTensor(const double*, const double* gradu, const double* uConv, int)
{
  if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE == 0)
  {
    // no turbulent viscosity
    return 0.0;
  }

  int viscosityTensor = TDatabase::ParamDB->TURBULENT_VISCOSITY_TENSOR;

  // compute square of the Frobenius norm of the tensor using 
  // deformation tensor
  switch (viscosityTensor)
  {
    case 0: // 0.5*(grad(u) + grad(u)^T)
    {
      double a11 = gradu[0] + gradu[0];
      double a12 = gradu[1] + gradu[3];
      double a13 = gradu[2] + gradu[6];
      double a22 = gradu[4] + gradu[4];
      double a23 = gradu[5] + gradu[7];
      double a33 = gradu[8] + gradu[8];

      return 0.5 * (a12 * a12 + a13 * a13 + a23 * a23)
        + 0.25 * (a11 * a11 + a22 * a22 + a33 * a33);
    }

    case 1: // grad u form
    {
      double frobenius_norm_tensor = 0.0;

      for (int i = 0; i < 9; i++)
      {
        frobenius_norm_tensor += gradu[i] * gradu[i];
      }

      return frobenius_norm_tensor;
    }

    case 2:
    {
      // deformation tensor of small scales
      // compute (grad(u)+grad(u)^T)/2 - G^H works only with VMS methods
      if (uConv == nullptr)
      {
        ErrThrow("TURBULENT_VISCOSITY_TENSOR 2 works only with VMS methods !!!");
      }

      double a11 = 0.5 * (gradu[0] + gradu[0]) - uConv[0];
      double a12 = 0.5 * (gradu[1] + gradu[3]) - uConv[1];
      double a13 = 0.5 * (gradu[2] + gradu[6]) - uConv[2];
      double a22 = 0.5 * (gradu[4] + gradu[4]) - uConv[3];
      double a23 = 0.5 * (gradu[5] + gradu[7]) - uConv[4];
      double a33 = 0.5 * (gradu[8] + gradu[8]) - uConv[5];

      return 2.0 * (a12 * a12 + a13 * a13 + a23 * a23)
        + a11 * a11 + a22 * a22 + a33 * a33;
    }

    default:
      ErrThrow("Turbulent viscosity tensor: ", viscosityTensor,
        " is not implemented so far");
  }

  return 0.0;
}

double turbulentViscosityDelta(double hK)
{
  double filter_constant = TDatabase::ParamDB->FILTER_WIDTH_CONSTANT;
  double filter_power = TDatabase::ParamDB->FILTER_WIDTH_POWER;

  if (filter_power == 1)
  {
    return filter_constant * hK;
  }
  else
  {
    return filter_constant*std::pow(hK,filter_power);
  }
}

double sgsKineticEnergy3D(double hK, const double* u, const double* gradu,
                            const double* uConv, const double* x,
                            const double* y, const double* z, double proj,
                            double scale)
{
  double delta = turbulentViscosityDelta(hK);
  double nu_T = turbulentViscosity3D(hK, u, gradu, uConv, x, y, z, proj);
  double sqrt_k = 0.5 * nu_T / (scale * delta);

  return sqrt_k * sqrt_k;
}

double turbulentViscosity3D(double hK, const double* u, const double* gradu,
                            const double* uConv, const double* x,
                            const double* y, const double* z, double)
{
  double delta = turbulentViscosityDelta(hK);
  double viscosity;

  // compute the frobenius norm of the relevant tensor
  double frobenius_norm_tensor = frobeniusNormTensor(u, gradu, uConv);

  int viscosityType = TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE;
  switch(viscosityType)
  {
    case 1:// Smagorinsky
      viscosity = standard(delta, frobenius_norm_tensor);
      break;
    case 2: // p laplacian
      viscosity = laplacian(delta, frobenius_norm_tensor);
      break;
    case 3: // Layton SIAM J. Sci. Comput. 1996 
      viscosity = getViscosityLayton96(delta, frobenius_norm_tensor);
      break;
    case 4: // \|u-g_\delta \ast u\|_2
      viscosity = normof_velocityminus_gdeltaTimesuStar(delta, u, uConv);
      break;
    case 5: // \|Du^h-G^H\|_2
      viscosity = normof_Defvelocityminus_gdeltaTimesuStar(delta, gradu, uConv);
      break;
    case 6: // all parameters free
      viscosity = parametersFree(delta, frobenius_norm_tensor);
      break;
    case 7: // van Driest damping for channel
      viscosity = vanDriestDampingChannel(delta, frobenius_norm_tensor, z);
      break;
    case 8: // van Driest damping for cylinder with squared cross--section
            // left and right wall at the cylinder
      viscosity = vanDriestDampingCylinder(delta, frobenius_norm_tensor, x, y);
      break;
    case 9:  // van Driest damping for channel flow (continuous)
      viscosity = vanDriestDampingChannelContinuous(delta, frobenius_norm_tensor, z);
      break;
    case 10: // van Driest damping for channel flow (paper: Rudman, Blackburn'99)
      viscosity = vanDriestDampingChannelRB99(delta, frobenius_norm_tensor, z);
      break;
    case 11:// van Driest damping for channel flow 
            // (paper: Rudman, Blackburn'99) with diff A+
      viscosity = vanDriestDampingChannelRB99APlus(delta, frobenius_norm_tensor, z);
      break;
    case 12: // van Driest damping (continuous, classical) for cylinder 
             // with squared cross--section
      viscosity = vanDriestDampingCylinderContinuous(delta, frobenius_norm_tensor, x, y);
      break;
    case 13:  // van Driest damping (paper: Rudman, Blackburn'99) 
              // for cylinder with squared cross--section
      viscosity = vanDriestDampingCylinderRB99(delta, frobenius_norm_tensor, x, y);
      break;
    case 14: // Verstappen model (J Sci Comput'11)
      viscosity = VerstappenViscositymodelCmuMax(delta, gradu, frobenius_norm_tensor);
      break;
    case 15: // Verstappen model (J Sci Comput'11)
      viscosity = VerstappenViscositymodelChPhi(delta, gradu, frobenius_norm_tensor);
      break;
    case 16: // eddy viscosity model: Vreman, Phys. Fluids 16 (10), 3670 -3681, 2004
             // frobenius norm of gradient of velocity
             // use same notations as in paper
      viscosity = vremanViscosityModel(gradu);
      break;
    case 17: // simple Verstappen model (J Sci Comput'11, p. 97)
      viscosity = VerstappenModelSimple(delta, gradu, frobenius_norm_tensor);
      break;
    case 18: // sigma model (Nicoud et al., Using singular values to build a subgrid-scale model for large eddy simulations, 2011)
      viscosity = sigmaViscosityModel(delta, gradu);
      break;
    case 110: // local Smagorinsky for husbandry
      viscosity = standard(delta, frobenius_norm_tensor);
      if (x[0] < TDatabase::ParamDB->P8)
      {
         viscosity *= TDatabase::ParamDB->P7;
      }
      break;
    default:
      ErrThrow("Turbulent viscosity for viscosity type : ", viscosityType, " is not implemented");
  }

  return (viscosity);
}

double DivDivStab3D(double u1, double u2, double u3, double hK, double eps) 
{
  int divdiv_type = TDatabase::ParamDB->DIV_DIV_STAB_TYPE;
  double c_1 = TDatabase::ParamDB->DIV_DIV_STAB_C1;
  double c_2 = TDatabase::ParamDB->DIV_DIV_STAB_C2;
  double tau;
    
  switch(divdiv_type)
  {
    // constant
    case 0:
      tau = c_2;
      break;
      // for non inf-sup stable fe, Codina, Gravemeier
    case 1:
      tau = c_2*std::sqrt(u1*u1+u2*u2+u3*u3)*hK/c_1;
      tau = tau*tau + eps*eps ;
      tau = std::sqrt(tau);
      break;
    case 2:
      tau = std::sqrt(u1*u1+u2*u2+u3*u3)*hK/c_1;
      tau = tau*tau + eps*eps ;
      tau = c_2*std::sqrt(tau);
      break;
    default:
      ErrThrow("div-div stabilization ", divdiv_type, " not implemented !!!");
      break;
  }
  return(tau);
}
