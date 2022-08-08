// =============================================================================
// Blas1.C
//
// Purpose:     basic routines for linear algebra
//
// Author:      Gunar Matthies          17.01.2003
//
// =============================================================================

#include <string.h>
//#include <math.h>
#include <cmath>

extern "C" {
// vector-vector sum: y = alpha * x + y
void daxpy_(int* n, double* alpha, const double* const x, int* incx,
            const double* y, int* incy);
}

/** return inner product (x,y) */
double Ddot(int n, const double *x, const double *y)
{
  double r;
  int i;
  const double *a, *b;

  r = 0.0;
  a = x;
  b = y;
  for(i=0;i<n;i++)
  {
    r += *a * *b;
    a++;
    b++;
  }

  return r;
}

/** y := alpha*x + y */
void Daxpy(int n, double alpha, const double* const x, double* y)
{
//   int i;
//   const double *a;
//   double *b;
//   double scal;
// 
//   a = x;
//   b = y;
//   scal = alpha;
//   for(i=0;i<n;i++)
//   {
//     *b += scal * *a;
//     a++;
//     b++;
//   }

  int incx = 1;
  int incy = 1;
  daxpy_(&n, &alpha, x, &incx, y, &incy);

}

/** x := alpha*x */
void Dscal(int n, double alpha, double *x)
{
  int i;
  double scal;
  double *a;

  scal = alpha;
  a = x;
  for(i=0; i<n; i++)
  {
    *a *= scal;
    a++;
  }
}

/** return Euclidian norm of x */
double Dnorm(int n, const double *x)
{
  double r;
  int i;
  const double *a;
  double z;

  a = x;
  r = 0.0;
  for(i=0; i<n; i++)
  {
    z = *a;
    r += z*z;
    a++;
  }

  return std::sqrt(r);
}