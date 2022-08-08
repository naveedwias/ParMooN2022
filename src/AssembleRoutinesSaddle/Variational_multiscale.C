#include <cmath>

#include "Database.h"
#include "Variational_multiscale.h"
#include "MooNMD_Io.h"

void matrices_a_tilde(std::vector<std::shared_ptr<FEMatrix>> &a,
                        std::shared_ptr<const FEMatrix> g, std::shared_ptr<const FEMatrix> gt,
                        std::shared_ptr<const FEMatrix> m, std::string flag)
{
  int nActive=0, nDof=0;

#ifdef __2D__
   nActive = a[0]->GetAnsatzSpace2D()->get_n_active();
   nDof = a[0]->GetAnsatzSpace2D()->get_n_dof();
#else
   nActive = a[0]->GetAnsatzSpace3D()->get_n_active();
   nDof = a[0]->GetAnsatzSpace3D()->get_n_dof();
#endif

  std::vector<double> value(nDof,0);

  int begin, end;
  double* entries_a11 = nullptr;
  double* entries_a22 = nullptr;
  double* entries_a33 = nullptr;
  double* entries_a = nullptr;

  // matrix entries
  if (a.size() == 3)
  {
    entries_a11 = a[0]->GetEntries(); // entries of matrix a11
    entries_a22 = a[1]->GetEntries(); // entries of matrix a22
    entries_a33 = a[2]->GetEntries(); // entries of matrix a33
  }
  else
  {
    entries_a = a[0]->GetEntries(); // off diagonal matrix a
  }

  const double* entries_gt = gt->GetEntries(); // entries of matrix G tilde
  const double* entries_g = g->GetEntries(); // entries of matrix G
  const double* entries_m_vms = m->GetEntries(); // mass matrix from vms

  for (int i = 0; i < nActive; ++i)
  {
    // multiplication of matrices
    {
      // i-th row of  Matrix tildeG
      begin = gt->get_row_ptr()[i];
      end = gt->get_row_ptr()[i + 1];

      for (int j = begin; j < end; ++j)
      {
        // row of the matrix G11

        int index = gt->get_vector_columns()[j];
        int begin_1 = g->get_row_ptr()[index];
        int end_1 = g->get_row_ptr()[index + 1];

        double val = entries_gt[j] / entries_m_vms[index];

        for (int k = begin_1; k < end_1; ++k)
        {
          int index1 = g->get_vector_columns()[k];
          value.at(index1) += val * entries_g[k];
        }
      }
    }

    // add matrix entries to A
    begin = a[0]->get_row_ptr()[i];
    end = a[0]->get_row_ptr()[i + 1];

    for (int j = begin; j < end; ++j)
    {
      int index = a[0]->get_vector_columns()[j];

      if (flag.compare("a11_tilde") == 0)
      {
        entries_a11[j] -= value.at(index);
        entries_a22[j] -= 0.5 * value.at(index);
        entries_a33[j] -= 0.5 * value.at(index);
      }
      else if (flag.compare("a22_tilde") == 0)
      {
        entries_a11[j] -= 0.5 * value.at(index);
        entries_a22[j] -= value.at(index);
        entries_a33[j] -= 0.5 * value.at(index);
      }
      else if (flag.compare("a33_tilde") == 0)
      {
        entries_a11[j] -= 0.5 * value.at(index);
        entries_a22[j] -= 0.5 * value.at(index);
        entries_a33[j] -= value.at(index);
      }
      else
      {
        entries_a[j] -= 0.5 * value.at(index);
      }

      value.at(index) = 0.0;
    }
  }
}
//====================================================================================

template<int d>
void LumpMassMatrixToDiagonalMatrix(std::shared_ptr<FEMatrix> & matrix)
{
  if (d == 2)
  {
    ErrThrow("Variational MultiScale method is not supported in 2D yet");
  }

  double* entries = matrix->GetEntries();

  int* rowptr = matrix->get_row_ptr();
  int* kcol = matrix->get_vector_columns();

  int rows = matrix->get_n_rows();

  for (int i = 0; i < rows; ++i)
  {
    int begin = rowptr[i];
    int end   = rowptr[i + 1];

    for (int j = begin; j < end; ++j)
    {
      // diagonal entry
      if (kcol[j] == i)
      {
        rowptr[i] = i;
        entries[i] = entries[j];
        kcol[i] = i;
        break;
      }
    }

    if (std::abs(entries[i]) < 1e-10)
    {
      entries[i] = 1.0;
    }
  }

  rowptr[rows] = rows;
}

template<int d>
void VMS_ProjectionUpdateMatrices(std::vector<std::shared_ptr< FEMatrix>>& blocks,
                                  std::vector<std::shared_ptr<FEMatrix>> matrices_vms)
{
  if (d == 2)
  {
    ErrThrow("Variational MultiScale method is not supported in 2D yet");
  }

  std::vector<std::shared_ptr<FEMatrix>> A(3);
  A.at(0) = blocks.at(0);  // A11
  A.at(1) = blocks.at(5);  // A22
  A.at(2) = blocks.at(10); // A33

  std::shared_ptr<FEMatrix> M = matrices_vms.at(6);

  std::shared_ptr<FEMatrix> GT = matrices_vms.at(0);
  std::shared_ptr<FEMatrix> G = matrices_vms.at(3);

  // diagonal A tilde matrices
  // (tilde_A11, tilde_A22, tilde_A33)

  matrices_a_tilde(A, G, GT, M, "a11_tilde");

  GT = matrices_vms.at(1);
  G = matrices_vms.at(4);
  matrices_a_tilde(A, G, GT, M, "a22_tilde");

  GT = matrices_vms.at(1);
  G = matrices_vms.at(4);
  matrices_a_tilde(A, G, GT, M, "a33_tilde");

  // off diagonal A tilde matrices
  // tilde_A12
  A.resize(1);
  A.at(0) = blocks.at(1);  // A12
  GT = matrices_vms.at(1); // GT22
  G = matrices_vms.at(3);  // G11
  matrices_a_tilde(A, G, GT, M, "a12_tilde");

  // tilde_A13
  A.resize(1);
  A.at(0) = blocks.at(2);
  GT = matrices_vms.at(2); // GT33
  G = matrices_vms.at(3);  // G11
  matrices_a_tilde(A, G, GT, M, "a13_tilde");

  // tilde_A21
  A.resize(1);
  A.at(0) = blocks.at(4);
  GT = matrices_vms.at(0); // GT11
  G = matrices_vms.at(4);  // G22
  matrices_a_tilde(A, G, GT, M, "a21_tilde");

  //tilde_A23
  A.resize(1);
  A.at(0) = blocks.at(6);
  GT = matrices_vms.at(2); // GT33
  G = matrices_vms.at(4);  // G22
  matrices_a_tilde(A, G, GT, M, "a23_tilde");

  // tilde_A31
  A.resize(1);
  A.at(0) = blocks.at(8);
  GT = matrices_vms.at(0); // GT11
  G = matrices_vms.at(5);  // G22
  matrices_a_tilde(A, G, GT, M, "a31_tilde");

  // tilde_A32
  A.resize(1);
  A.at(0) = blocks.at(9);
  GT = matrices_vms.at(1); // GT22
  G = matrices_vms.at(5);  // G33
  matrices_a_tilde(A, G, GT, M, "a32_tilde");
}

RBVMS_Settings::RBVMS_Settings()
  : mode(RBVMS_ParamMode::G),
  time_discretization(RBVMS_TimeDiscretization::BackwardEuler),
  delta_0(0.25), delta_1(0.1),
  C_inv(1.0), tau_mul(0.25),
  tau_m_time(false), square_hk(true),
  explicit_time_derivative(false),
  momentum_pressure_coupling_B(false), momentum_pressure_coupling_C(false)
{
}

RBVMS_Settings::RBVMS_Settings(const ParameterDatabase& db)
{
  if (db["rbvms_param_mode"].is("g"))
  {
    mode = RBVMS_ParamMode::G;
  }
  else if (db["rbvms_param_mode"].is("h"))
  {
    mode = RBVMS_ParamMode::H;
  }
  else if (db["rbvms_param_mode"].is("codina"))
  {
    mode = RBVMS_ParamMode::Codina;
  }
  else
  {
    mode = RBVMS_ParamMode::G;
    ErrThrow("Unknown RBMVS parameter mode ", db["rbvms_param_mode"]);
  }

  if (db["rbvms_time_discretization"].is("backward_euler"))
  {
    time_discretization = RBVMS_TimeDiscretization::BackwardEuler;
  }
  else if (db["rbvms_time_discretization"].is("crank_nicolson"))
  {
    time_discretization = RBVMS_TimeDiscretization::CrankNicolson;
  }
  else
  {
    time_discretization = RBVMS_TimeDiscretization::BackwardEuler;
    ErrThrow("Unknown RBVMS time discretization ", db["rbvms_time_discretization"]);
  }

  delta_0 = db["rbvms_delta_0"];
  delta_1 = db["rbvms_delta_1"];

  tau_m_time = db["rbvms_tau_m_time"];
  square_hk = db["rbvms_square_hk"];

  C_inv = db["rbvms_c_inv"];
  tau_mul = db["rbvms_tau_mul"];

  explicit_time_derivative = db["rbvms_explicit_time_derivative"];
  momentum_pressure_coupling_B = db["rbvms_momentum_pressure_coupling_B"];
  momentum_pressure_coupling_C = db["rbvms_momentum_pressure_coupling_C"];
}

template<int d>
std::pair<double, double> RBVMS_Param_G(const double* w, const double* G,
  const double* g, double tau, double nu, const RBVMS_Settings& settings)
{
  double inv_tau_m_squared = 0.0;

  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < d; j++)
    {
      inv_tau_m_squared += w[i] * w[j] * G[i * d + j];
    }
  }

  inv_tau_m_squared = std::abs(inv_tau_m_squared);

  double GF2 = 0.0;

  for (int i = 0; i < d * d; i++)
  {
    GF2 += G[i] * G[i];
  }

  inv_tau_m_squared += settings.C_inv * nu * nu * GF2;
  inv_tau_m_squared += settings.tau_mul / (tau * tau);

  double tau_m = 1.0 / std::sqrt(inv_tau_m_squared);

  double tau_c = 0.0;
  for (int i = 0; i < d; i++)
  {
    tau_c += g[i] * g[i];
  }

  tau_c = 1.0 / (tau_m * tau_c);

  return std::make_pair(tau_m, tau_c);
}

template<int d>
std::pair<double, double> RBVMS_Param_H(double hK, const RBVMS_Settings& settings)
{
  double tau_m = settings.delta_0 * hK;

  if (settings.square_hk)
  {
    tau_m *= hK;
  }

  double tau_c = settings.delta_1;

  if (settings.tau_mul > 0.0 && settings.tau_m_time)
  {
    tau_m = std::max(tau_m, std::sqrt(1.0 / settings.tau_mul) * TDatabase::TimeDB->CURRENTTIMESTEPLENGTH);
  }

  return std::make_pair(tau_m, tau_c);
}

template<int d>
std::pair<double, double> RBVMS_Param_Codina(double hK, double nu,
  const double* u, const double* u_prime, const RBVMS_Settings& settings)
{
  if (u_prime == nullptr)
  {
    ErrThrow("Cannot use Codina parameters without available pointwise data!");
  }

  double u1 = u[0] + u_prime[0];
  double u2 = u[1] + u_prime[1];
  double u3 = d == 2 ? 0.0 : (u[2] + u_prime[2]);
  double mag_u = std::sqrt(u1 * u1 + u2 * u2 + u3 * u3);

  double tau_m = hK * hK / (settings.delta_1 * nu + hK * settings.delta_0 * mag_u);
  double tau_c = hK * hK / (settings.delta_1 * tau_m);

  return std::make_pair(tau_m, tau_c);
}

template<int d>
void RBMVS_Time_EvolveSubscale(const double* old_data, double* new_data,
  const double* old_res, const double* new_res,
  double tau_m, double delta_t, const RBVMS_Settings& settings)
{
  switch (settings.time_discretization)
  {
    case RBVMS_TimeDiscretization::CrankNicolson:
    {
      double u0_mul = (2 * tau_m - delta_t) / (2 * tau_m + delta_t);
      double r_mul = -tau_m * delta_t / (2 * tau_m + delta_t);

      for (int i = 0; i < d; i++)
      {
        new_data[i] = u0_mul * old_data[i] + r_mul * (new_res[i] + old_res[i]);
      }

      break;
    }

    case RBVMS_TimeDiscretization::BackwardEuler:
    {
      double u0_mul = tau_m / (tau_m + delta_t);
      double r1_mul = -delta_t * u0_mul;

      for (int i = 0; i < d; i++)
      {
        new_data[i] = u0_mul * old_data[i] + r1_mul * new_res[i];
      }

      break;
    }

    default:
    {
      ErrThrow("Unsupported time discretization for rbvms_time");
      break;
    }
  }
}

#ifdef __2D__

template void LumpMassMatrixToDiagonalMatrix<2>(std::shared_ptr<FEMatrix> & matrix);
template void VMS_ProjectionUpdateMatrices<2>(std::vector<std::shared_ptr<FEMatrix>>& blocks,
                     std::vector<std::shared_ptr<FEMatrix>> matrices_vms);

template std::pair<double, double> RBVMS_Param_G<2>(const double* w, const double* G,
  const double* g, double tau, double nu, const RBVMS_Settings& settings);
template std::pair<double, double> RBVMS_Param_H<2>(double hK,
  const RBVMS_Settings& settings);
template std::pair<double, double> RBVMS_Param_Codina<2>(double hK, double nu,
  const double* u, const double* u_prime, const RBVMS_Settings& settings);

template void RBMVS_Time_EvolveSubscale<2>(const double* old_data, double* new_data,
  const double* old_res, const double* new_res,
  double tau_m, double delta_t, const RBVMS_Settings& settings);

#endif
#ifdef __3D__

template void LumpMassMatrixToDiagonalMatrix<3>(std::shared_ptr<FEMatrix> & matrix);
template void VMS_ProjectionUpdateMatrices<3>(std::vector<std::shared_ptr<FEMatrix>>& blocks,
                     std::vector<std::shared_ptr<FEMatrix>> matrices_vms);

template std::pair<double, double> RBVMS_Param_G<3>(const double* w, const double* G,
  const double* g, double tau, double nu, const RBVMS_Settings& settings);
template std::pair<double, double> RBVMS_Param_H<3>(double hK,
  const RBVMS_Settings& settings);
template std::pair<double, double> RBVMS_Param_Codina<3>(double hK, double nu,
  const double* u, const double* u_prime, const RBVMS_Settings& settings);

template void RBMVS_Time_EvolveSubscale<3>(const double* old_data, double* new_data,
  const double* old_res, const double* new_res,
  double tau_m, double delta_t, const RBVMS_Settings& settings);

#endif