// =======================================================================
// LocalProjection.h
//
// Purpose:   routines for local projection stabilization
//
// Author:    Gunar Matthies  2007/03/06
//
// =======================================================================

#ifndef __LOCAL_PROJECTION__
#define __LOCAL_PROJECTION__

#include <LinAlg.h>

#ifdef __2D__
void CoupledDefect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        TMatrix *B1T, TMatrix *B2T, TMatrix *C,
        double *x, double *b, double *r);

//void Defect_NSE2C(TSquareMatrix **A, TMatrix **B, double *x,
//                  double *b, double *r);

void CoupledMatVect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        TMatrix *B1T, TMatrix *B2T, TMatrix *C,
        double *x, double *y);

void MatVect_NSE2C(TSquareMatrix **A, TMatrix **B, double *x, double *y);

void UltraLocalProjection(void* A, bool ForPressure, CoeffFct2D *Coeff);
void UltraLocalProjection(void* A, bool ForPressure);
void UltraLocalProjectionSD(void* A, bool ForPressure);

double UltraLocalError(TFEFunction2D *uh, DoubleFunct2D *ExactU,
        double lpcoeff, double lpexponent, int OrderDiff);

double UltraLocalErrorDivergence(TFEFunction2D *uh1, TFEFunction2D *uh2,
                       DoubleFunct2D *ExactU1, DoubleFunct2D *ExactU2,
                       double lpcoeff, double lpexponent, int OrderDiff);

double UltraLocalErrorStreamline(TFEFunction2D *uh, DoubleFunct2D *ExactU,
                       TFEFunction2D *b1, TFEFunction2D *b2,
                       double lpcoeff, double lpexponent, int OrderDiff);

double UltraLocalErrorStreamlinePWConst(TFEFunction2D *uh, DoubleFunct2D *ExactU,
                       TFEFunction2D *b1, TFEFunction2D *b2,
                       double lpcoeff, double lpexponent, int OrderDiff);

void AddStreamlineTerm(TSquareMatrix2D* A, TFEFunction2D *uh1,
                       TFEFunction2D *uh2,
                       double lpcoeff, double lpexponent, int OrderDiff);

void AddStreamlineTermPWConst(TSquareMatrix2D* A, TFEFunction2D *uh1,
                              TFEFunction2D *uh2,
                              double lpcoeff, double lpexponent, int OrderDiff);

void AddDivergenceTerm(TSquareMatrix2D *A11,TSquareMatrix2D *A12,
                       TSquareMatrix2D *A21,TSquareMatrix2D *A22,
                       double lpcoeff, double lpexponent, int OrderDiff);

void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A21,
                    TSquareMatrix *A22, TMatrix *B1, TMatrix *B2,
                    TMatrix *B1T, TMatrix *B2T,
                    TMatrix *C,
                    double *x, double *y);

void MatVect_NSE4C(TSquareMatrix **A, TMatrix **B, double *x, double *y);

void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A21,
                   TSquareMatrix *A22, TMatrix *B1, TMatrix *B2,
                   TMatrix *B1T, TMatrix *B2T,
                   TMatrix *C,
                   double *x, double *b, double *r);

void Defect_NSE4C(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r);

double UltraLocalErrorSmooth(TFEFunction2D *uh, DoubleFunct2D *ExactU,
                             double lpcoeff, double lpexponent, int OrderDiff);
bool TestCell(TBaseCell *cell);

void UltraLocalProjectionFunction(void* A, bool ForPressure);

void UltraLocalProjectionStreamlinePLaplacian(TSquareMatrix2D* A, 
                                              TFEFunction2D *uh, 
                                              const CoeffFct2D& Coeff);

void LocalProjectionCoarseGridQ0(TFEFunction2D *uh, TFEFunction2D *uh_proj,
                                 const CoeffFct2D& Coeff, int convection_flag);  
         
void LocalProjectionCrossWindCoarseGridQ0(TDomain *Domain, int mg_level,
                                          TFEFunction2D *uh,
                                          TFEFunction2D *uh_proj,
                                          const CoeffFct2D& Coeff,
                                          double *rhs, int convection_flag); 


void AdaptivePostProcess(TFEFunction2D *FeFunction, double *PostSol, bool DirichletBC);


void AddALEStreamlineLPS(TSquareMatrix2D* A, int N_FeFunct, TFEFunction2D **FeFunct,
                         double lpcoeff, double lpexponent, int OrderDiff);
#else // __3D__ 

void AddStreamlineTerm(TSquareMatrix3D* A, TFEFunction3D *uh1,
                       TFEFunction3D *uh2, TFEFunction3D *uh3,
                       double lpcoeff, double lpexponent, int OrderDiff); 

void UltraLocalProjection(TSquareMatrix3D* A, 
                          double lpcoeff, double lpexponent, int OrderDiff);

FE_type GetElement3D(TBaseCell *cell, int CoarseOrder);

void UltraLocalProjection3D(void* A, bool ForPressure);

double UltraLocalError3D(TFEFunction3D *uh, DoubleFunct3D *ExactU,
        double lpcoeff, double lpexponent, int OrderDiff);


#endif // __3D__

#endif // __LOCAL_PROJECTION__
