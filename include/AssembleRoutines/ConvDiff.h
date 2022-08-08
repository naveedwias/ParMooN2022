// ======================================================================
// @(#)ConvDiff.h        12/06/26
//
// common declaration for all convection diffusion problems
// ======================================================================

#ifndef __CONVDIFF__
#define __CONVDIFF__

#include <Collection.h>
#include "Constants.h"

#define HMM86            0 
#define TP86_1           1
#define TP86_2           2 
#define JSW87            3 
#define GdC88            4 
#define dCG91            5
#define dCA03            6 
#define AS97             7
#define C93              8 
#define KLR02_1          9
#define KLR02_2          10
#define KLR02_3          11
#define KLR02_4          12
#define J90              13
#define BE02_1           14 
#define BE02_2           15 
#define BH04             16 
#define BE05_1           17 
#define BE05_2           18 
#define LP96             19
#define CS99             20
#define MH_Kno06         21
#define BE02_3           22 
#define Y_Z_BETA         24 
#define JSW87_1          25 
#define FEM_TVD          52
#define GENERAL_SOLD    200

template <int d>
double Mesh_size_in_convection_direction_without_storing(
  double hK, std::array<double, d> b);

template <int d>
double Mesh_size_in_convection_direction(double hK, std::array<double, d> b);

template <int d>
double Compute_SDFEM_delta(double hK, double eps, std::array<double, d> b,
                           double react, double linfb);

template <int d>
double Compute_SOLD_sigma(double hK, double eps, std::array<double, d> b,
                          double c, double f, double linfb, double deltaK,
                          double *param, double residual, int residual_computed,
                          int time_dependent_problem);

/** coercivity constant of problem 
 * used eg for residual based estimator of Verf"uhrt 2005 */
double EstimateCoercivityConstant(TCollection *Coll, 
#ifdef __2D__
                                  const CoeffFct2D& Coeff
#else // 3D
                                  const CoeffFct3D& Coeff
#endif
                                 );

/** it should be i>=0 and i<= 31, otherwise error*/
void SetSoldParameters(int i);


#ifdef __2D__
void EdgeStabilization(TFESpace2D *fespace,  TFEFunction2D *u, 
                       const CoeffFct2D& Coeffs, double *rhs,
                       int time_dependent,
                       double *time_step, TFEFunction2D *old_u);
#endif

double ComputeAlpha(double hK);

//=============================================================================
// local assembling routines:

template<int d>
void BilinearAssembleGalerkin(double Mult, const double *coeff,
                              const double* param, double hK,
                              const double **OrigValues,
                              const int *N_BaseFuncts, double ***LocMatrices,
                              double **LocRhs);
// SUPG/SDFEM
// SDFEM - Streamline Diffusion Finite Element Method,
// SUPG - Streamline Upwind Petrov Galerkin
template<int d>
void BilinearAssemble_SD(double Mult, const double *coeff, const double* param,
                         double hK, const double **OrigValues,
                         const int *N_BaseFuncts, double ***LocMatrices,
                         double **LocRhs);
// Galerking least squares
template<int d>
void BilinearAssemble_GLS(double Mult, const double *coeff, const double* param,
                         double hK, const double **OrigValues,
                          const int *N_BaseFuncts, double ***LocMatrices,
                          double **LocRhs);

//=============================================================================
// local error computing routines:
template<int d>
void conv_diff_l2_h1_linf_error(int N_Points, std::array<const double*, d> xyz,
                                const double *AbsDetjk, const double *Weights,
                                double hK, const double *const* Der,
                                const double *const* Exact,
                                const double *const* coeffs, double *LocError);

template<int d>
void PoissonAssemble(double Mult, const double *coeff,
                                 const double*, double,
                                 const double **OrigValues,
                                 const int *N_BaseFuncts, double ***LocMatrices,
                                 double **LocRhs);
template<int d>
void AssembleCDRHomogeneous(double Mult, const double *coeff,
                                 const double*, double,
                                 const double **OrigValues,
                                 const int *N_BaseFuncts, double ***LocMatrices,
                                 double **LocRhs);

#endif // __CONVDIFF__
