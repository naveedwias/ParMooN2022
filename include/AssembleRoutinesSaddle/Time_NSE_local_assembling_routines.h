#ifndef INCLUDE_ASSEMBLEROUTINESSADDLE_TIME_NSE_LOCAL_ASSEMBLING_ROUTINES_H
#define INCLUDE_ASSEMBLEROUTINESSADDLE_TIME_NSE_LOCAL_ASSEMBLING_ROUTINES_H

#include <Variational_multiscale.h>

/** Store only the local assembling routines that are additional from the 
 * time stepping schemes, e.g., the mass matrices, additional matrices in the 
 * case of stabilization schemes. Note that: we can still use the local assembly
 * of velocity-velocity coupling from the steady case but due 
 * to the ordering of matrices in the Assemble function, rectangular blocks 
 * will have different number's because of the additional mass matrices. Keeping 
 * this in mind, velocity-pressure coupling blocks are copied here with correct numbering. 
 */

template<int d>
void NSMassMatrixSingle(double Mult, const double *coeff, const double *,
                        double, const double **OrigValues,
                        const int *N_BaseFuncts, double ***LocMatrices,
                        double **);

template<int d>
void NSMassMatrix(double Mult, const double *, const double *, double,
                  const double **OrigValues, const int *N_BaseFuncts,
                  double ***LocMatrices, double **);

template<int d> 
void NSLaplaceGradGradSingleSmagorinsky(double Mult, const double *coeff,
                                        const double *param, double hK,
                                        const double **OrigValues,
                                        const int *N_BaseFuncts,
                                        double ***LocMatrices, double **);

template<int d> 
void NSLaplaceGradGradSmagorinsky(double Mult, const double *coeff,
                                  const double *param, double hK,
                                  const double**OrigValues,
                                  const int *N_BaseFuncts,
                                  double ***LocMatrices, double **);

template<int d> 
void NSLaplaceDeformationSmagorinsky(double Mult, const double *coeff,
                                     const double *param, double hK,
                                     const double**OrigValues,
                                     const int *N_BaseFuncts,
                                     double ***LocMatrices, double **);

template<int d> 
void NSLaplaceDeformationVariationalMS(double Mult, const double *,
                                       const double *param, double hK,
                                       const double**OrigValues,
                                       const int *N_BaseFuncts,
                                       double ***LocMatrices, double **);

template<int d>
void NSLumpMassMatrix(double Mult, const double *, const double *, double ,
                      const double**OrigValues, const int *N_BaseFuncts,
                      double ***LocMatrices, double **);

template<int d>
void NSVariationlMS_GMatrices(double Mult, const double *, const double *,
                              double , const double**OrigValues,
                              const int *N_BaseFuncts, double ***LocMatrices,
                              double **);

template<int d>
void NSVariationlMS_GTildeMatrices(double Mult, const double *,
                                   const double *param, double hK,
                                   const double**OrigValues,
                                   const int *N_BaseFuncts,
                                   double ***LocMatrices, double **);

template<int d>
void NSParamVelGradSmagorinsky(const double *in, double *out);

template<int d>
void NSParamsVariationalMSLargeScale(const double *in, double *out);

///////////////////////////////////////////////////////////////////////////////
// Residual-based VMS terms (residual_based_vms)

// the residuals are separated out like this because the viscosity \nu and the
// RHS f are retrieved from coeffs later during assembly, which 

// residuals:
//
// first d components:
// out[0..d] = res_m + \nu \Delta u + f = du/dt + (u \cdot \nabla) u + \nabla p,
//
// d-th component:
// out[d] = res_c = \nabla \cdot u
template <int d>
void NSParamsVMSResidualsWithoutLaplacianAndRHS(const double *in, double *out,
  bool extend_advection = false);

// laplacian term:
//
// d components:
//
// out = \Delta u
template <int d>
void NSParamsVMSResidualsLaplacian(const double *in, double *out);

// used for ansatz time derivative
//
// d components:
//
// out = u_old
template<int d>
void NSParamsOldVelocity(const double *in, double *out);

// completes the residuals computed in the parameter stage by retrieving the
// viscosity, RHS, and reference transform data and computes the LHS residual
// terms:
//
// (res_m, \nabla q) + (res_c, \nabla \cdot v)
// - n(res_m; u, v) - n(u; res_m, v)
// + n(res_m; res_m, v),
//
// where the terms in the first line represent the subgrid continuity
// constraint,
// the terms in the second line represent the coupling between subgrid
// and coarse components,
// and the third line term represents subgrid convection.
template<int d>
void NSVMSResiduals(double Mult, const double *coeff, const double *param,
  double hK, const double **OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **LocRhs, int inv_jacobian_offset,
  RBVMS_Settings settings);

template<int d>
void NSVMSResiduals_MassMatrix(double Mult, const double *coeff, const double *param,
  double hK, const double **OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **LocRhs, int inv_jacobian_offset,
  RBVMS_Settings settings);

template<int d>
void NSVMSResiduals_RHS(double Mult, const double *coeff, const double *param,
  double hK, const double **OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **LocRhs, int inv_jacobian_offset,
  RBVMS_Settings settings);

template<int d>
void TNSE_RBVMS_Time_MassMatrix(double Mult, const double *coeff, const double *param,
  double hK, const double **OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **, int inv_jacobian_offset,
  RBVMS_Settings settings);

template<int d>
void TNSE_RBVMS_Time_SubgridTerms_RHS(double Mult, const double *coeff,
  const double *param, double hK, const double **OrigValues,
  const int *N_BaseFuncts, double ***, double **LocRhs,
  int inv_jacobian_offset, RBVMS_Settings settings);

template<int d>
void TNSE_RBVMS_Time_SubgridTerms_Matrices(double Mult, const double *coeff,
  const double *param, double hK, const double **OrigValues,
  const int *N_BaseFuncts, double ***LocMatrices, double **,
  int inv_jacobian_offset, RBVMS_Settings settings);

/** everything below is for the Residual based methods
 * SUPG
 * Rb-VMS
 * NOTE: these methods are only implemented for 
 * NSTYPE 4   inf-sup stable pair
 * NSTYPE 14  equal-order pairs with pressure stabilization
 */
// SUPG (Streamline-Upwind Petrov-Galerkin
template<int d>
void NS_SUPG(
  double Mult, const double *coeff, const double *param, double hK, 
  const double **OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **, 
double delta0, double delta1);

template<int d>
void NS_SUPG_GradientBlocks(
  double Mult, const double *, const double *param, double hK, 
  const double **OrigValues, const int *N_BaseFuncts, 
  double ***LocMatrices, double **, double delta0);

template<int d> 
void NS_SUPG_RightHandSide_InfSup(
  double Mult, const double *coeff, const double *param, double hK, 
  const double **OrigValues, const int *N_BaseFuncts, 
  double ***, double **LocRhs, double delta1);

template <int d> 
void NS_SUPG_MassMatrix(
  double Mult, const double *coeff, const double *param, double hK, 
  const double **OrigValues, const int *N_BaseFuncts, 
  double ***LocMatrices, double **, double delta0);

template<int d>
void NS_SUPG_EquOrder_Gradient_DivergenceBlocks(
  double Mult, const double *coeff, const double *param, double , 
  const double **OrigValues, const int *N_BaseFuncts, double ***LocMatrices, 
  double **, double factor);

template<int d> 
void NS_SUPG_RightHandSide_EquOrder(
  double Mult, const double *coeff, const double *param, double hK, 
  const double **OrigValues, const int *N_BaseFuncts, double ***, 
  double **LocRhs, double factor);
////////////////////////////////////////////////////////////////////////////////
// local assembling different nonlinear forms for supg
template <int d>
void NS_SUPG_skew_symmetric(
  double Mult, const double *coeff, const double *param, double hK, 
  const double **OrigValues, const int *N_BaseFuncts, 
  double ***LocMatrices, double **, double delta0, double delta1);

template <int d>
void NS_SUPG_rotational(
  double Mult, const double *coeff, const double *param, double hK, 
  const double **OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **, double delta0, double delta1);

template <int d>
void NS_SUPG_emac(
  double Mult, const double *coeff, const double *param, double hK, 
  const double **OrigValues, const int *N_BaseFuncts,
  double ***LocMatrices, double **, double delta0, double delta1);
////////////////////////////////////////////////////////////////////////////////

template<int d>
void NSParamsVelocityDerivatives_SUPG_inf_sup(const double *in, double *out);
template<int d>
void NSParamsVelocityDerivatives_SUPG_equal_order(const double *in, double *out);
//===============================================================================
// stabilization parameters for equal order case
//===============================================================================
template<int d>
void compute_stabilization_parameters(const double* u, const double* coeff, 
                                            double* params);
#endif
