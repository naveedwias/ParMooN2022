#ifndef INCLUDE_ASSEMBLEROUTINES_CD_LOCAL_ASSEMBLING_ROUTINES_H
#define INCLUDE_ASSEMBLEROUTINES_CD_LOCAL_ASSEMBLING_ROUTINES_H

///////////////////////////////////////////////////////////////////////////////
// standard terms (needed for Galerkin)

template <int d>
void TCDStiff(double Mult, const double *coeff, const double *param, double hK,
              const double**OrigValues, const int *N_BaseFuncts,
              double ***LocMatrices, double **LocRhs);
// add a term of the form: div (a_l * b*b^T grad u), for a vector b, and a scalar a_l
template <int d>
void TCDStiff_TensorialDiffusionTerm(double Mult, const double *coeff,
                                     const double *param, double hK,
                                     const double**OrigValues,
                                     const int *N_BaseFuncts,
                                     double ***LocMatrices, double **LocRhs);
template <int d>
void TCDMass(double Mult, const double *coeff, const double *param, double hK,
             const double**OrigValues, const int *N_BaseFuncts,
             double ***LocMatrices, double **LocRhs);

template <int d>
void TCDRhs(double Mult, const double *coeff, const double *param, double hK,
            const double**OrigValues, const int *N_BaseFuncts,
            double ***LocMatrices, double **LocRhs);

///////////////////////////////////////////////////////////////////////////////
// stabilization terms
template <int d>
void TCDStiffSUPG(double Mult, const double *coeff, const double *param,
                  double hK, const double**OrigValues, const int *N_BaseFuncts,
                  double ***LocMatrices, double **LocRhs);
template <int d>
void TCDMassSUPG(double Mult, const double *coeff, const double *param,
                 double hK, const double**OrigValues, const int *N_BaseFuncts,
                 double ***LocMatrices, double **LocRhs);
template <int d>
void TCDRhsSUPG(double Mult, const double *coeff, const double *param,
                double hK, const double**OrigValues, const int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs);

///////////////////////////////////////////////////////////////////////////////
// more special terms

// for ROM regularized initial condition (as solution of Helmoltz equation)
template <int d>
void TCDGradGrad(double Mult, const double *coeff, const double *param,
                 double hK, const double**OrigValues, const int *N_BaseFuncts,
                 double ***LocMatrices, double **LocRhs);
template <int d>
void TCDStiffWithoutDiffAndLinearSource(double Mult, const double* coeff,
                                        const double *param, double hK, const double** OrigValues,
                                        const int* N_BaseFuncts, double*** LocMatrices,
                                        double **LocRhs);
template <int d>
void TCDDiffAndLinearSourceOnly(double Mult, const double* coeff, const double *param,
                                double hK, const double** OrigValues,
                                const int *N_BaseFuncts, double*** LocMatrices,
                                double **LocRhs, int matrixIndex);

#endif // INCLUDE_ASSEMBLEROUTINES_CD_LOCAL_ASSEMBLING_ROUTINES_H
