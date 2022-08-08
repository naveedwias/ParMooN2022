#include <MooNMD_Io.h>
#include <NonNewtonianViscosity.h>
#include <cmath>
#include <string.h>

extern "C" 
{
  //  computes the eigenvalues and eigenvectors
  void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w,
              double* work, int* lwork, int* info);
}

// Carreau-Yasuda fluid:
// nu = nu_infty + (nu_0 - nu_infty) (1 + (lambda g)^a )^((n - 1)/ a),
// where g = \dot \gamma is the shear rate
double CarreauYasudaViscosity(double shear_rate, const ViscositySettings& settings)
{
  if (settings.lambda * shear_rate > 0.0)
  {
    return settings.nu_infty
    + (settings.nu_0 - settings.nu_infty)
    * std::pow(
      1.0 + std::pow(settings.lambda * shear_rate, settings.a),
      (settings.n - 1.0) / settings.a);
  }
  else
  {
    return settings.nu_0;
  }
}

// Clamped power-law fluid: nu = clamp(k g^(n - 1), nu_0, nu_infty),
// where g = \dot \gamma is the shear rate
double PowerLawViscosity(double shear_rate, const ViscositySettings& settings)
{
  double nu = settings.k * std::pow(std::max(shear_rate, 1.e-16), settings.n - 1.0);

  return std::max(settings.nu_0, std::min(settings.nu_infty, nu));
}

// Cross fluid: nu = nu_infty + (nu_0 - nu_infty) / (1 + (kg)^n),
// where g = \dot \gamma is the shear rate
double CrossViscosity(double shear_rate, const ViscositySettings& settings)
{
  if (settings.k * shear_rate > 0.0)
  {
    return settings.nu_infty
    + (settings.nu_0 - settings.nu_infty)
    / (1.0 + std::pow(settings.k * shear_rate, settings.n));
  }
  else
  {
    return settings.nu_0;
  }
}

// regularized Herschel-Bulkley fluid:
// nu = tau_0 / g + k g^(n - 1),
// where g = max(g_0, \dot \gamma) is the regularized shear rate.
// note that g_0 is precomputed from nu_0.
double HerschelBulkleyViscosity(double shear_rate, const ViscositySettings& settings)
{
  shear_rate = std::max(shear_rate, settings.g_0);

  return settings.tau_0 / shear_rate
  + settings.k * std::pow(shear_rate, settings.n - 1.0);
}

ViscositySettings::ViscositySettings()
  : mode(ViscosityMode::Newtonian),
  nu_0(1.0), nu_infty(1.0),
  tau_0(0.0), k(1.0), n(1.0),
  lambda(1.0), a(2.0), g_0(1.0)
{
}

ViscositySettings::ViscositySettings(const ParameterDatabase &db)
{
  mode = ViscosityMode::Newtonian;

  if (db["viscosity_mode"].is("newtonian"))
  {
    mode = ViscosityMode::Newtonian;
  }
  else if (db["viscosity_mode"].is("carreau_yasuda"))
  {
    mode = ViscosityMode::CarreauYasuda;
  }
  else if (db["viscosity_mode"].is("power_law"))
  {
    mode = ViscosityMode::PowerLaw;
  }
  else if (db["viscosity_mode"].is("cross"))
  {
    mode = ViscosityMode::Cross;
  }
  else if (db["viscosity_mode"].is("herschel_bulkley"))
  {
    mode = ViscosityMode::HerschelBulkley;
  }
  else
  {
    ErrThrow("Unknown viscosity_mode: '", db["viscosity_mode"], "'");
  }

  nu_0 = db["viscosity_nu_0"];
  nu_infty = db["viscosity_nu_infty"];
  tau_0 = db["viscosity_tau_0"];
  k = db["viscosity_k"];
  n = db["viscosity_n"];
  lambda = db["viscosity_lambda"];
  a = db["viscosity_a"];

  g_0 = db["viscosity_g_0"];

  if (mode == ViscosityMode::HerschelBulkley)
  {
    nu_0 = HerschelBulkleyViscosity(g_0, *this);

    Output::root_info("Herschel-Bulkley fluid", "Adjusted limiting viscosity "
        "to match regularization threshold:"
        "\n\\dot \\gamma_0 = ", g_0,
        "\n\\nu_0 = ", nu_0);
  }
}

typedef double ViscosityFunc(double, const ViscositySettings&);

template<int d>
void NonNewtonianViscosityAnisotropic(const double* grad_u, double* nu,
  const ViscositySettings& settings)
{
  ViscosityFunc* visc = nullptr;

  switch (settings.mode)
  {
    case ViscosityMode::Newtonian:
      for (int i = 0; i < d * d; i++)
      {
        nu[i] = settings.nu_0;
      }
      return;

    case ViscosityMode::CarreauYasuda:
      visc = &CarreauYasudaViscosity;
      break;

    case ViscosityMode::PowerLaw:
      visc = &PowerLawViscosity;
      break;

    case ViscosityMode::Cross:
      visc = &CrossViscosity;
      break;

    case ViscosityMode::HerschelBulkley:
      visc = &HerschelBulkleyViscosity;
      break;

    default:
      ErrThrow("Unsupported viscosity mode!");
      break;
  }

  for (int i = 0; i < d; i++)
  {
    nu[d * i + i] = visc(2.0 * std::abs(grad_u[d * i + i]), settings);

    for (int j = i + 1; j < d; j++)
    {
      double nu_ij = visc(std::abs(grad_u[d * i + j] + grad_u[d * j + i]), settings);
      nu[d * i + j] = nu_ij;
      nu[d * j + i] = nu_ij;
    }
  }
}

template<int d>
double NonNewtonianViscosity(const double* grad_u,
  const ViscositySettings& settings)
{
  int r = d;
  int info = 0;
  double D[d * d];
  double ev[d];
  double work[3 * d];
  int lwork = 3 * d;
  char N = 'N';
  char U = 'U';

  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < d; j++)
    {
      D[d * i + j] = grad_u[d * i + j] + grad_u[d * j + i];
    }
  }

  dsyev_(&N, &U, &r, D, &r, ev, work, &lwork, &info);

  double shear_rate = 0.5 * std::max(std::abs(ev[d - 1]), std::abs(ev[0]));

  switch (settings.mode)
  {
    case ViscosityMode::Newtonian:
      return settings.nu_0;

    case ViscosityMode::CarreauYasuda:
      return CarreauYasudaViscosity(shear_rate, settings);

    case ViscosityMode::PowerLaw:
      return PowerLawViscosity(shear_rate, settings);

    case ViscosityMode::Cross:
      return CrossViscosity(shear_rate, settings);

    case ViscosityMode::HerschelBulkley:
      return HerschelBulkleyViscosity(shear_rate, settings);

    default:
      ErrThrow("Unsupported viscosity mode!");
      return 0.0;
  }
}

template<int d>
void NSLaplaceGradGradSingleNonNewtonian(double Mult, const double *,
                                        const double *param, double,
                                        const double **OrigValues,
                                        const int *N_BaseFuncts,
                                        double ***LocMatrices, double **,
                                        const ViscositySettings& settings)
{
  double **MatrixA = LocMatrices[0];

  int N_U = N_BaseFuncts[0];

  const double *u_x = OrigValues[2];
  const double *u_y = OrigValues[3];
  const double *u_z = d == 2 ? nullptr : OrigValues[4];

  const double *grad_u = param + d;

  const double m_nu = Mult * NonNewtonianViscosity<d>(grad_u, settings);

  for(int i = 0; i < N_U; i++)
  {
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0.0 : u_z[i];

    for(int j = 0; j < N_U; j++)
    {
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0.0 : u_z[j];

      MatrixA[i][j] += m_nu * (test_x * ansatz_x + test_y * ansatz_y
                               + test_z * ansatz_z);
    }
  }
}

template <int d>
void NSLaplaceGradGradNonNewtonian(double Mult, const double *,
                                  const double *param, double,
                                  const double **OrigValues,
                                  const int *N_BaseFuncts,
                                  double ***LocMatrices, double **,
                                  const ViscositySettings& settings)
{
  double **MatrixA11 = LocMatrices[0];
  double **MatrixA22 = LocMatrices[d + 1];
  double **MatrixA33 = d == 2 ? nullptr : LocMatrices[2 * d + 2];

  int N_U = N_BaseFuncts[0];

  const double *u_x = OrigValues[2];
  const double *u_y = OrigValues[3];
  const double *u_z = d == 2 ? nullptr : OrigValues[4];

  const double *grad_u = param + d;

  const double m_nu = Mult * NonNewtonianViscosity<d>(grad_u, settings);

  for (int i = 0; i < N_U; i++)
  {
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0.0 : u_z[i];

    for (int j = 0; j < N_U; j++)
    {
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0.0 : u_z[j];

      double val = m_nu * (test_x * ansatz_x + test_y * ansatz_y
                         + test_z * ansatz_z);

      MatrixA11[i][j] += val;
      MatrixA22[i][j] += val;

      if (d == 3)
      {
        MatrixA33[i][j] += val;
      }
    }
  }
}

template <int d>
void NSLaplaceDeformationNonNewtonian(double Mult, const double *,
                                     const double *param, double,
                                     const double **OrigValues,
                                     const int *N_BaseFuncts,
                                     double ***LocMatrices, double **,
                                     const ViscositySettings& settings)
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

  int N_U = N_BaseFuncts[0];

  const double *u_x = OrigValues[2];
  const double *u_y = OrigValues[3];
  const double *u_z = d == 2 ? nullptr : OrigValues[4];

  const double *grad_u = param + d;

  const double m_nu = Mult * NonNewtonianViscosity<d>(grad_u, settings);

  for (int i = 0; i < N_U; i++)
  {
    double test_x = u_x[i];
    double test_y = u_y[i];
    double test_z = d == 2 ? 0.0 : u_z[i];

    for (int j = 0; j < N_U; j++)
    {
      double ansatz_x = u_x[j];
      double ansatz_y = u_y[j];
      double ansatz_z = d == 2 ? 0.0 : u_z[j];

      if (d == 2)
      {
        MatrixA11[i][j] += m_nu * (2.0 * test_x * ansatz_x + test_y * ansatz_y);
        MatrixA12[i][j] += m_nu * test_y * ansatz_x;

        MatrixA21[i][j] += m_nu * test_x * ansatz_y;
        MatrixA22[i][j] += m_nu * (test_x * ansatz_x + 2.0 * test_y * ansatz_y);
      }
      else
      {
        MatrixA11[i][j] += m_nu * (2.0 * test_x * ansatz_x + test_y * ansatz_y + test_z * ansatz_z);
        MatrixA12[i][j] += m_nu * test_y * ansatz_x;
        MatrixA13[i][j] += m_nu * test_z * ansatz_x;

        MatrixA21[i][j] += m_nu * test_x * ansatz_y;
        MatrixA22[i][j] += m_nu * (test_x * ansatz_x + 2.0 * test_y * ansatz_y + test_z * ansatz_z);
        MatrixA23[i][j] += m_nu * test_z * ansatz_y;

        MatrixA31[i][j] += m_nu * test_x * ansatz_z;
        MatrixA32[i][j] += m_nu * test_y * ansatz_z;
        MatrixA33[i][j] += m_nu * (test_x * ansatz_x + test_y * ansatz_y + 2.0 * test_z * ansatz_z);
      }
    }
  }
}

template double NonNewtonianViscosity<2>(const double*, const ViscositySettings&);
template double NonNewtonianViscosity<3>(const double*, const ViscositySettings&);

template void NSLaplaceGradGradSingleNonNewtonian<2>(double, const double *,
  const double *, double, const double **, const int *, double ***, double **,
  const ViscositySettings&);
template void NSLaplaceGradGradSingleNonNewtonian<3>(double, const double *,
  const double *, double, const double **, const int *, double ***, double **,
  const ViscositySettings&);

template void NSLaplaceGradGradNonNewtonian<2>(double, const double *,
  const double *, double, const double **, const int *, double ***, double **,
  const ViscositySettings&);
template void NSLaplaceGradGradNonNewtonian<3>(double, const double *,
  const double *, double, const double **, const int *, double ***, double **,
  const ViscositySettings&);

template void NSLaplaceDeformationNonNewtonian<2>(double, const double *,
  const double *, double, const double **, const int *, double ***, double **,
  const ViscositySettings&);
template void NSLaplaceDeformationNonNewtonian<3>(double, const double *,
  const double *, double, const double **, const int *, double ***, double **,
  const ViscositySettings&);