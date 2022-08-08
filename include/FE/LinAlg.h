#ifndef __LINALG__
#define __LINALG__

// =============================================================================
// BLAS 1
// =============================================================================

/** return inner product (x,y) */
double Ddot(int n, const double *x, const double *y);

/** y := alpha*x + y */
void Daxpy(int n, double alpha, const double* const x, double* y);

/** x := alpha*x */
void Dscal(int n, double alpha, double *x);

/** return Euclidian norm of x */
double Dnorm(int n, const double *x);

// =============================================================================
// BLAS 2
// =============================================================================

/** solve a system of linear equations */
void SolveLinearSystemLapack(double *a, double *b, int N_Eqn, int LDA);

/* subroutine for solving a multiple systems of linear equations */
/* solution is transposed */
void SolveMultipleSystems(double *a, double *b, int N_Eqn,
                          int LDA, int LDB, int N_Rhs);

/* subroutine for solving a multiple systems of linear equations */
void SolveMultipleSystemsNew(double *a, double *b, int N_Eqn,
                             int LDA, int LDB, int N_Rhs);

// =============================================================================
// Multi grid components
// =============================================================================


#endif
