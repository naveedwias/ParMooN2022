/** ****************************************************************************
*
* @name       ROM_local_assembling_routines
* @brief      specific local assembling routines for TNSE ROM problems
*
* Store only the local assembling routines that are additional from the
* Time Navier Stokes ROM, the terms used for the assembling of the pressure.
*
* @note all the following functions but one are expecting the use of the type
*       BlockFEMatrix::NSE*D_Type*. Only the function TCDDiffOnlyScale, which is
*      comparable to TCDDiffOnly, is expecting type BlockFEMatrix::CD*D
*
*******************************************************************************/

#ifndef ROM_LOCAL_ASSEMBLING_ROUTINES_H
#define ROM_LOCAL_ASSEMBLING_ROUTINES_H


// Assembling the pressure matrix: hK * (grad p, grad q)
// similar to TCDDiffOnlyScale (TO BE USED with BlockFEMatrix::CD*D)
template<int d>
void TCDDiffOnlyScale(double Mult, const double *coeff, const double *param,
                      double hK, const double **OrigValues,
                      const int *N_BaseFuncts, double ***LocMatrices,
                      double **LocRhs);

// Assembling the pressure-velocity matrix: hK * (u, grad q)
template<int d>
void NSPROM_Divergence(double Mult, const double *coeff, const double *param,
                       double hK, const double **OrigValues,
                       const int *N_BaseFuncts, double ***LocMatrices,
                       double **LocRhs);

// Assembling the pressure-velocity matrix: hK * ( (u.grad)u, grad q )
template<int d>
void NSPROM_NonlinearTerm(double Mult, const double *coeff, const double *param,
                          double hK, const double **OrigValues,
                          const int *N_BaseFuncts, double ***LocMatrices,
                          double **LocRhs);

// Assembling the pressure rhs: hK * (f, grad q)
// similar to NSPSPG_RightHandSide()
template<int d>
void NSPROM_RightHandSide(double Mult, const double *coeff, const double *param,
                          double hK, const double **OrigValues,
                          const int *N_BaseFuncts, double ***LocMatrices,
                          double **LocRhs);

#endif // ROM_LOCAL_ASSEMBLING_ROUTINES_H
