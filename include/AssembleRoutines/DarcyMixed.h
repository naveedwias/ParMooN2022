#ifndef __DARCY_MIXED_2D__ 
#define __DARCY_MIXED_2D__

/** the local assembling routines. Each of them corresponds to one 
 * LocalAssembling2D_type */

// ======================================================================
// Standard Galerkin with Raviart-Thomas (RT) or Brezzi-Douglas-Marini (BDM)
// elements
template <int d>
void BilinearAssembleDarcyGalerkin(double Mult, const double *coeff,
                                   const double *param, double hK,
                                   const double **OrigValues,
                                   const int *N_BaseFuncts,
                                   double ***LocMatrices, double **LocRhs);

#endif // __DARCY_MIXED_2D__
