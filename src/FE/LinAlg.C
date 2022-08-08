// =============================================================================
// @(#)LinAlg.C        1.18 07/03/00
//
// Purpose:     basic routines for linear algebra
//
// Author:      Gunar Matthies          27.01.1999
//              Volker John             27.10.1999
//              Sashikumaar Ganesan     08.10.2009 (eigen values)
// =============================================================================

#include <string.h>
//#include <math.h>
#include <cmath>
#include <stdlib.h>

#include <MooNMD_Io.h>
#include "BaseCell.h"
#include <LinAlg.h>


extern "C" {

void dgetrf_(int *m, int *n, double *A, int *lda, int *ipiv, int *info);

void dgetrs_(char *trans, int *n, int *nrhs, double *A, int *lda, int *ipiv,
             double *B, int *ldb, int *info);

}


// =============================================================================
// general routines independent of dimension
// =============================================================================
#define AT(i,j) (a[j*LDA+i])
#define A(i,j) (a[i*LDA+j])

void SolveLinearSystemLapack(double *a, double *b, int N_Eqn, int)
{
// Arguments:
//    a         double array which contains the matrix columnwise
//              a[i,j] = a[i*LDA+j]
//    b         on input: rhs
//              on output: solution
//    N_Eqn     number of equations
//    LDA       leading dimension of matrix a
  int m, n, nrhs, lda, ldb;
  int *ipivot, info;
  char t='n';

  m = N_Eqn;
  n = N_Eqn;
  lda = N_Eqn;
  ldb = N_Eqn;
  nrhs = 1;

  ipivot = new int[n];
  
  dgetrf_(&m, &n, a, &lda, ipivot, &info);
  dgetrs_(&t, &n, &nrhs, a, &lda, ipivot, b, &ldb, &info);

  delete ipivot;
}

void SolveMultipleSystems(double *a, double *b, int N_Eqn,
                          int LDA, int LDB, int N_Rhs)
// Arguments:
//    a         double array which contains the matrix row wise
//              a[i,j] = a[i*LDA+j]
//    b         on input: rhs
//              on output: solution
//    N_Eqn     number of equations
//    LDA       leading dimension of matrix a
//    LDB       leading dimension of vector b
//    N_Rhs     number of right hand sides
//
// This works for a -- stored row wise
//                b -- stored row wise, b = (b1,b2,b3, ...)
//                LDA = LDB
// The result will be stored row wise
{
  int i,j,l,m, row;
  double pivot, tmp, *f;

//  int ii, jj;

/*
  for(ii=0;ii<N_Eqn;ii++)
  {
    cout << ii;
    for(jj=0;jj<N_Eqn;jj++)
      cout << setw(8) << A(ii,jj);
    cout << endl;
  }
  cout << endl;

  cout << "number of Rhs: " << N_Rhs << endl;
*/

  // for all columns
  for(i=0;i<N_Eqn-1;i++)
  {
    // compute pivot element
    pivot = 0;
    row = i;
    // check current column from diagonal element
    for(l=i;l<N_Eqn;l++)
    {
      // a_li
      tmp = std::abs(A(l,i));
      if(tmp > pivot)
      {
        pivot = tmp;
        row = l;
      } // endif
    } // endfor l
    if(pivot == 0.0)
    {
      ErrThrow("Error in solving multiple Systems");
    }
    // change rows i and 'row' if necessary
    if(i<row)
    {
      for(l=0;l<N_Eqn;l++)
      {
        // N_Eqn subsequent entries since a is stored row wise
        pivot = A(i,l);
        A(i,l) = A(row,l);
        A(row, l) = pivot;
      }
      // same for rhs
      for(j=0;j<N_Rhs;j++)
      {
        // since rhs is stored row wise, find corresponding index in each rhs
        f = b+j*LDB;
        tmp = f[i];
        f[i] = f[row];
        f[row] = tmp;
      }
    } // endif

    // apply pivoting
    tmp = A(i,i);
    // current column
    for(l=i+1;l<N_Eqn;l++)
    {
      A(l,i) /= tmp;
    }
    // remainder of the matrix
    for(l=i+1;l<N_Eqn;l++)
    {
      tmp = A(i,l);
      for(m=i+1;m<N_Eqn;m++)
      {
        A(m,l) -= A(m,i) * tmp;
      }
    }
  } // endfor i

  for(i=0;i<N_Eqn;i++)
  {
    for(j=0;j<N_Rhs;j++)
    {
      f = b+j*LDB;
      tmp = f[i];
      for(l=i+1;l<N_Eqn;l++)
      {
        f[l] -= A(l,i)*tmp;
      }
    }
  }

  for(i=N_Eqn-1;i>=0;i--)
  {
    for(j=0;j<N_Rhs;j++)
    {
      f = b+j*LDB;
      f[i] /= A(i,i);
      tmp = f[i];
      for(l=0;l<i;l++)
      {
        f[l] -= A(l,i)*tmp;
      }
    }
  }
}

/* subroutine for solving a multiple systems of linear equations */
void SolveMultipleSystemsNew(double *a, double *b, int N_Eqn,
                             int LDA, int LDB, int N_Rhs)
// Arguments:
//    a         double array which contains the matrix row wise
//    b         on input: rhs
//              on output: solution
//    N_Eqn     number of equations
//    LDA       leading dimension of matrix a
//    LDB       leading dimension of vector b
//    N_Rhs     number of right hand sides
// This works for a -- stored row wise
//                b -- stored column wise, b = (b1_1,..., b_NRhs_1, b1_2, ...)
//                LDB = N_Rhs
// The result will be stored column wise
{
  int i,j,l,m, row;  // k, info, *ipiv;
  double pivot, tmp, *f, *frow,eps=0.0;

//  int ii, jj;
/*
  for(ii=0;ii<N_Eqn;ii++)
  {
    for(jj=0;jj<N_Eqn;jj++)
    {
      cout << "a(" << ii+1 << "," << jj+1 << ") = " << setw(8);
      cout << A(ii,jj) << ";" << endl;
    }
  }
  cout << endl;
  for(jj=0;jj<N_Rhs;jj++)
  {
    cout << jj+1 << ": ";
    for(ii=0;ii<N_Eqn;ii++)
      cout << setw(15) << b[ii*LDB+jj];
    cout << endl;
  }
*/

  // LU decomposition of matrix A with pivot search
  for(i=0;i<N_Eqn-1;i++)
  {
    pivot = 0;
    row = i;
    // find pivot
    for(l=i;l<N_Eqn;l++)
    {
      tmp = std::abs(A(l,i));
      if(tmp > pivot)
      {
        pivot = tmp;
        row = l;
      } // endif
    } // endfor l
    if(pivot <= eps)
    {
      ErrThrow("Error in solving multiple Systems ", "equation: ", i);
    }
    // change rows if necessary
    if(i<row)
    {
      for(l=0;l<N_Eqn;l++)
      {
        pivot = A(i,l);
        A(i,l) = A(row,l);
        A(row, l) = pivot;
      }
      // row-th row of rhs
      frow = b+row*LDB;
      // i-th row of rhs
      f = b+i*LDB;
      for(j=0;j<N_Rhs;j++)
      {
        tmp = f[j];
        f[j] = frow[j];
        frow[j] = tmp;
      }
    } // endif

    tmp = A(i,i);
    for(l=i+1;l<N_Eqn;l++)
    {
      A(l,i) /= tmp;
    }

    for(l=i+1;l<N_Eqn;l++)
    {
      tmp = A(i,l);
      for(m=i+1;m<N_Eqn;m++)
      {
        A(m,l) -= A(m,i) * tmp;
      }
    }

  } // endfor i

  frow = b;
  // solve left lower system
  for(i=0;i<N_Eqn;i++)
  {
    // i-th row of rhs
    f = b+i*LDB;
    // for all rhs
    for(j=0;j<N_Rhs;j++)
    {
      tmp = f[j];
      // for remaining rows
      for(l=i+1;l<N_Eqn;l++)
      {
        frow[l*LDB+j] -= A(l,i)*tmp;
      }
    }
  }
  // solve right upper system
  for(i=N_Eqn-1;i>=0;i--)
  {
    // i-th row of rhs
    f = b+i*LDB;
    // for all rhs
    for(j=0;j<N_Rhs;j++)
    {
      f[j] /= A(i,i);
      tmp = f[j];
      for(l=0;l<i;l++)
      {
        frow[l*LDB+j] -= A(l,i)*tmp;
      }
    }
  }
}
