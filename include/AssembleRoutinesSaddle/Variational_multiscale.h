// =======================================================================
// Variational_multiScale.h
//
// Purpose:     routines for projection-based VMS
// =======================================================================

#ifndef __Variational_multiScale__
#define __Variational_multiScale__

#include <memory>
#include<FEMatrix.h>
#include<BlockFEMatrix.h>
#include "ParameterDatabase.h"

#ifdef __2D__
#include <FEVectFunct2D.h>
#endif

#ifdef __3D__
#include<FEVectFunct3D.h>
#endif

#include "templateNames.h"
#include "vector"

/** creates a diagonal mass matrix */
template<int d>
void LumpMassMatrixToDiagonalMatrix(std::shared_ptr<FEMatrix> & matrix);

/**
 *  block=
 *  [ A11  A12  A13  B1T ]
 *  [ A21  A22  A23  B2T ]
 *  [ A31  A32  A33  B3T ]
 *  [ B1   B2   B3   C   ]
 */

/**
 * matrices_vms[] ={Gt11, GT22, GT33, G11, G22, G33, M};
 */
template<int d>
void VMS_ProjectionUpdateMatrices(std::vector<std::shared_ptr<FEMatrix>>& blocks,
                                  std::vector<std::shared_ptr<FEMatrix>> matrices_vms);

enum class RBVMS_ParamMode { G, H, Codina };
enum class RBVMS_TimeDiscretization { BackwardEuler, CrankNicolson };

struct RBVMS_Settings
{
  public:
    RBVMS_Settings();
    RBVMS_Settings(const ParameterDatabase &db);

    RBVMS_ParamMode mode;
    RBVMS_TimeDiscretization time_discretization;

    double delta_0;
    double delta_1;

    double C_inv;
    double tau_mul;

    bool tau_m_time;
    bool square_hk;

    bool explicit_time_derivative;
    bool momentum_pressure_coupling_B;
    bool momentum_pressure_coupling_C;
};

/// @brief Computes RBVMS stabilization parameters \tau_m and \tau_c,
/// where
///
/// \tau_m^-2 = 4 / \tau^2 + (w, Gw) + Cinv nu^2 |G|_F, and
/// \tau_c = 1 / (\tau_m |g|^2)
///
/// where \tau is the time step length, nu is the viscosity, Cinv is a constant
/// from an inverse estimate, w is the resolved velocity estimate, and G and g
/// are a matrix and vector derived from the inverse reference transformation's
/// Jacobian. Specifically, if X_i are reference coordinates and x_i are world
/// coordinates, let J_ij = \partial X_i / \partial x_j and
///
/// g_i =  \sum_k J_ki = \mathbb{1}^T J,
/// G_ij = \sum_k J_ki J_kj = J^T J.
///
/// See: Bazilevs et al., Variational multiscale residual-based turbulence
///      modeling for large eddy simulation of incompressible flows (2007)
template<int d>
std::pair<double, double> RBVMS_Param_G(const double* w, const double* G,
  const double* g, double tau, double nu,
  const RBVMS_Settings& settings);

/// @brief Computes RBVMS stabilization parameters \tau_m and \tau_c,
/// where
///
/// \tau_m = \delta_0 h_K^2
/// \tau_c = \delta_1
///
/// with nonnegative \delta_0, \delta_1. If tau_m_time is true:
///
/// \tau_m = \max \{ \delta_0 h_K^2, 1/2 \tau \}.
template<int d>
std::pair<double, double> RBVMS_Param_H(double hK,
  const RBVMS_Settings& settings);

template<int d>
std::pair<double, double> RBVMS_Param_Codina(double hK, double nu,
  const double* u, const double* u_prime, const RBVMS_Settings& settings);

template<int d>
void RBMVS_Time_EvolveSubscale(const double* old_data, double* new_data,
  const double* old_res, const double* new_res,
  double tau_m, double delta_t, const RBVMS_Settings& settings);

#endif
