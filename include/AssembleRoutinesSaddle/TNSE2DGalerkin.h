#ifndef TNSE2DGALERKIN_H
#define TNSE2DGALERKIN_H

// ===============================================================================
// Type 1, 2, 3 and 4, Standard Galerkin, grad(u), grad(v)
// ===============================================================================
void TimeNSType1Galerkin(double Mult, double *coeff, double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                double **LocRhs);

void TimeNSType2Galerkin(double Mult, double *coeff, double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                double **LocRhs);

void TimeNSType1_2NLGalerkin(double Mult, double *coeff, double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                double **LocRhs);

void TimeNSType3Galerkin(double Mult, double *coeff, double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                double **LocRhs);

void TimeNSType4Galerkin(double Mult, double *coeff, double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                double **LocRhs);

void TimeNSType3_4NLGalerkin(double Mult, double *coeff, double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                double **LocRhs);
// ===============================================================================
// Type 3 and 4, Standard Galerkin, D(u):D(v)
// ===============================================================================
void TimeNSType3GalerkinDD(double Mult, double *coeff, double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                double **LocRhs);
void TimeNSType4GalerkinDD(double Mult, double *coeff, double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                double **LocRhs);
void TimeNSType3_4NLGalerkinDD(double Mult, double *coeff, double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, 
                double **LocRhs);

// ======================================================================
// right-hand side ONLY
// ======================================================================
void TimeNSRHS(double Mult, double *coeff, 
               double *param, double hK, 
               double **OrigValues, int *N_BaseFuncts,
               double ***LocMatrices, double **LocRhs);

#endif // TNSE2DGALERKIN_H
