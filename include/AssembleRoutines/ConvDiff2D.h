// ======================================================================
// @(#)ConvDiff2D.h        1.13 10/19/99
//
// common declaration for all convection diffusion problems
// ======================================================================

#ifndef __CONVDIFF2D__
#define __CONVDIFF2D__

/** the local assembling routines. Each of them corresponds to one 
 * LocalAssembling2D_type */

/** ========================================================================= */
/** ========================================================================= */
// CD2D: stationary convection diffusion problems

void BilinearAssemble_Axial3D(double Mult, double *coeff, double* param,
                              double hK, double **OrigValues, int *N_BaseFuncts,
                              double ***LocMatrices, double **LocRhs);


/** ========================================================================= */
/** ========================================================================= */
// TCD2D: time dependent convection diffusion problems
// Galerkin: local assembling routine for mass matrix 
void LocalMatrixM(double Mult, double *coeff, double *param, double hK, 
                        double **OrigValues, int *N_BaseFuncts, 
                        double ***LocMatrices, double **LocRhs);
//Galerkin: local assembling routine for stiffness matrix and right hand side
void LocalMatrixARhs(double Mult, double *coeff, double *param, double hK, 
                        double **OrigValues, int *N_BaseFuncts,
                        double ***LocMatrices, double **LocRhs);

// SUPG: local assembling routine for stiffness matrix and right hand side
void LocalMatrixARhs_SUPG(double Mult, double *coeff, double *param,
                                double hK, double **OrigValues,
                                int *N_BaseFuncts,double ***LocMatrices,
                                double **LocRhs);
// SUPG: local routine which assembles the mass matrix and a 
// matrix which comes from the SUPG part of the time derivative
void LocalMatrixM_SUPG(double Mult, double *coeff, double *param,
                             double hK, double **OrigValues,int *N_BaseFuncts,
                             double ***LocMatrices, double **LocRhs);

/** ========================================================================= */
/** ========================================================================= */
/** ========================================================================= */
/** ========================================================================= */


void BilinearAssemble_UPW1(double Mult, double *coeff, double* param, double hK,
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs);

void BilinearAssemble_UPW2(double Mult, double *coeff, double* param, double hK,
                           double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs);

void BilinearAssemble_SOLD(double Mult, double *coeff, double *param,
                           double hK, double **OrigValues, int *N_BaseFuncts,
                           double ***LocMatrices, double **LocRhs);

void BilinearAssemble_SOLD_Orthogonal(double Mult, double *coeff,
                                      double *param, double hK,
                                      double **OrigValues, int *N_BaseFuncts,
                                      double ***LocMatrices, double **LocRhs);

void RhsAssemble_LP96(double Mult, double *coeff, double *param, double hK,
                      double **OrigValues, int *N_BaseFuncts,
                      double ***LocMatrices, double **LocRhs);

void BilinearAssemble_MH_Kno06(double Mult, double *coeff, double *param,
                               double hK,
                               double **OrigValues, int *N_BaseFuncts,
                               double ***LocMatrices, double **LocRhs);

void BilinearAssemble_SD_SOLD(double Mult, double *coeff, double *param, 
                              double hK, double **OrigValues, int *N_BaseFuncts,
                              double ***LocMatrices, double **LocRhs);
            
void BilinearAssemble2LevelLPS_Q0(double Mult, double *coeff, double *param,
                                  double hK, double **OrigValues,
                                  int *N_BaseFuncts, double ***LocMatrices,
                                  double **LocRhs);

void RhsAssemble_RhsAdjointEnergyEstimate(double Mult, double *coeff, 
                                          double *param, double hK,
                                          double **OrigValues, 
                                          int *N_BaseFuncts,
                                          double ***LocMatrices,
                                          double **LocRhs);
void RhsAssemble_RhsAdjointTV(double Mult, double *coeff, double *param,
                              double hK, double **OrigValues, int *N_BaseFuncts,
                              double ***LocMatrices, double **LocRhs);

void RhsAssemble_RhsAdjointTV2(double Mult, double *coeff, double *param,
                               double hK, double **OrigValues,
                               int *N_BaseFuncts, double ***LocMatrices,
                               double **LocRhs);

void RhsAssemble_RhsAdjointNormBL1_NormBorthL1(double Mult, double *coeff,
                                               double *param, double hK,
                                               double **OrigValues,
                                               int *N_BaseFuncts,
                                               double ***LocMatrices,
                                               double **LocRhs);

void RhsAssemble_RhsAdjointNormResidualL1_NormBorthL1(double Mult,
                                                      double *coeff,
                                                      double *param, double hK,
                                                      double **OrigValues,
                                                      int *N_BaseFuncts,
                                                      double ***LocMatrices,
                                                      double **LocRhs);

void RhsAssemble_RhsAdjointL2Error(double Mult, double *coeff, double *param,
                                   double hK, double **OrigValues,
                                   int *N_BaseFuncts, double ***LocMatrices,
                                   double **LocRhs);

void RhsAssemble_RhsAdjointH1Error(double Mult, double *coeff, double *param,
                                   double hK, double **OrigValues, 
                                   int *N_BaseFuncts, double ***LocMatrices,
                                   double **LocRhs);


// ========================================================================
// parameter functions
// ========================================================================



#endif // __CONVDIFF2D__
