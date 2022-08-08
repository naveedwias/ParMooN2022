/** **************************************************************************** 
*
* @name       DenseMatrix
* @brief      Implements class DenseMatrix declared in DenseMatrix.h.
*
* @author     Clemens Bartsch
* @date       2015/06/03
*
*******************************************************************************/

#include <MooNMD_Io.h>
#include <cmath>
#include <cassert> //not the most sophisticated solution, but fit four our case
#include <cstring>
#include <iterator>
#include <stdexcept>

#include <DenseMatrix.h>

// Extern declaration of Lapack methods.
extern "C" {
// LU factorize...
void dgetrf_(int *m, int *n, double *A, int *lda, int *ipiv, int *info);
// ... and solve!
void dgetrs_(char* transa, int* n, int* nrhs, double *a, int* lda,
             int *ipiv, double *b, int* ldb, int *info);

// matrix-matrix product C = alpha * op(A) * op(A) + beta * C
void dsyrk_(char* uplo, char* trans, int* n, int* k, double* alpha,
            const double* const a, int* lda,
            double* beta, double* c, int* ldc);

// matrix-matrix product C = alpha * op(A) * op(B) + beta * C
void dgemm_(char* transa, char* transb, int* m, int* n, int* k, double* alpha,
            const double* const a, int* lda, const double* const b, int* ldb,
            double* beta, double* c, int* ldc);

// matrix-vector product y = alpha * op(A) * x + beta * y
void dgemv_(char* trans, int* m,  int* n, double* alpha, const double* const a,
            int* lda, const double* const x,  int* incx, double* beta,
            double* y, int* incy);

// computes the selected eigenvalues and optionally eigenvectors
// using a Relatively Robust Representation
// (see: http://www.netlib.org/lapack/lug/node30.html)
void dsyevr_(char* job_type, char* range, char* uplo, int* n,
             const double* const a, int* lda, double* vl, double* vu,
             int* il, int* iu, double* tolerance, int* m, double* w,
             const double* const z, int* ldz, int* isuppz,
             double* work, int* lwork, int* iwork, int* liwork, int* info);
}

/** ***************************************************************************/
DenseMatrix::DenseMatrix(const size_t nRows, const size_t nColumns) :
  nRows_(nRows),
  nColumns_(nColumns)
{
  // The way we intend to used the matrix, leading dimension is number of rows.
  leadingDimension_ = nRows_;

  // Allocate the entries_ array and fill it with zeroes. "()".
  entries_ = new double[leadingDimension_*nColumns_]();

  // The memory needed for LU factorization is not allocated yet.
  pivotsLU_= nullptr;
  entriesLU_ = nullptr;
}

/** ***************************************************************************/
DenseMatrix::DenseMatrix(size_t nRows,
                         size_t nColumns,
                         size_t leadingDimension) :
  nRows_(nRows),
  nColumns_(nColumns),
  leadingDimension_(leadingDimension)
{
  if(leadingDimension < nRows)
  {
    // Leading dimension may not be lesser than number of rows.
    ErrThrow("Error in DenseMatrix leading dimension must be greater"
             " or equal number of rows.");
  }
  // Allocate the entries_ array and fill it with zeroes. The latter is
  // achieved due to "()".
  entries_ = new double[leadingDimension_*nColumns_]();

  // The memory needed for LU factorization is not allocated yet.
  pivotsLU_= nullptr;
  entriesLU_ = nullptr;
}

// TODO
/** ***************************************************************************/
DenseMatrix::DenseMatrix(const DenseMatrix &obj) :
  nRows_(obj.getNRows()), nColumns_(obj.getNColumns()),
  leadingDimension_(obj.getLeadingDimension())
{
  // Allocate the entries_ array.
  entries_ = new double[leadingDimension_*nColumns_];
  // And fill it with the entries from obj.
  memcpy(entries_,obj.entries_,leadingDimension_*nColumns_*sizeof(double));

  pivotsLU_= nullptr;
  entriesLU_ = nullptr;

  if(obj.pivotsLU_)
  {
    pivotsLU_ = new int[nColumns_];
    memcpy(pivotsLU_,obj.getPivotsLU(),nColumns_*sizeof(int));
  }

  if(obj.entriesLU_)
  {
    entriesLU_ = new double[leadingDimension_*nColumns_];
    memcpy(entriesLU_,obj.entries_,leadingDimension_*nColumns_*sizeof(double));
  }
}

/** ***************************************************************************/
DenseMatrix::DenseMatrix(DenseMatrix&& other)
{
  *this = std::move(other);
}

/** ***************************************************************************/
DenseMatrix& DenseMatrix::operator=(DenseMatrix other)
{
  //do a swap with the copy constructed object "other"
  swap(*this, other);

  return *this;
}

/** ***************************************************************************/
DenseMatrix::~DenseMatrix()
{
  delete[] entries_;
  delete[] pivotsLU_;
  delete[] entriesLU_;
}

/** ***************************************************************************/
DenseMatrix& DenseMatrix::operator+=(const DenseMatrix& B)
{
  this->add( B, 1.);
  return *this;
}

/** ***************************************************************************/
DenseMatrix& DenseMatrix::operator-=(const DenseMatrix& B)
{
  this->add( B, -1.);
  return *this;
}

/** ***************************************************************************/
void DenseMatrix::add(const std::vector<double>* const v,
                      double a, bool transpose)
{
  if(transpose == false)
  {
    if(v->size() != nRows_)
    {
      ErrThrow("Dimension mismatch during matrix-vector addition: ",
               nRows_, " rows in this->matrix and ",
               v->size(), " elements in the vector.");
    }
    for(int i =0; i < (int)nColumns_; i++)
    {
      for(int j = 0; j < (int)nRows_; j++)
      {
        entries_[j + i*leadingDimension_] += a * (*v)[j];
      }
    }
  }
  else
  {
    if(v->size() != nColumns_)
    {
      ErrThrow("Dimension mismatch during matrix-vector addition: ",
                nColumns_, " columns in this->matrix and ",
                v->size(), " elements in the vector.");
    }
    for(int i =0; i < (int)nColumns_; i++)
    {
      for(int j = 0; j < (int)nRows_; j++)
      {
        entries_[j + i*leadingDimension_] += a * (*v)[i];
      }
    }
  }
  // reset old LU values
  this->resetLU();
}

void DenseMatrix::addToColumn(const std::vector<double>* const v,
  int col_index, double a)
{
  if (v->size() != nRows_)
  {
    ErrThrow("Dimension mismatch during matrix-vector addition: ",
             nRows_, " rows in this->matrix and ",
             v->size(), " elements in the vector.");
  }

  for (int j = 0; j < (int)nRows_; j++)
  {
    entries_[j + col_index * leadingDimension_] += a * (*v)[j];
  }
}

void DenseMatrix::addToRow(const std::vector<double>* const v,
  int row_index, double a)
{
  if (v->size() != nColumns_)
  {
    ErrThrow("Dimension mismatch during matrix-vector addition: ",
             nColumns_, " columns in this->matrix and ",
             v->size(), " elements in the vector.");
  }

  for (int i = 0; i < (int)nColumns_; i++)
  {
    entries_[row_index + i * leadingDimension_] += a * (*v)[i];
  }
}

void DenseMatrix::scaleColumn(int col_index, double a)
{
  for (int j = 0; j < (int)nRows_; j++)
  {
    entries_[j + col_index * leadingDimension_] *= a;
  }
}

void DenseMatrix::scaleRow(int row_index, double a)
{
  for (int i = 0; i < (int)nColumns_; i++)
  {
    entries_[row_index + i * leadingDimension_] *= a;
  }
}

/** ***************************************************************************/
void DenseMatrix::add(const DenseMatrix& B, double a)
{
  if(  (int)this->nRows_ != B.getNRows()
    || (int)this->nColumns_ != B.getNColumns()) // compare sizes
  {
    ErrThrow("DenseMatrix::add: the two matrices do not match.");
  }

  int n_entries = nRows_ * nColumns_;
  const double *BEntries = B.get_entries();
  for(int i = 0; i < n_entries; ++i)
  {
    this->entries_[i] += a * BEntries[i];
  }
  // reset old LU values
  this->resetLU();
}

/** ***************************************************************************/
void DenseMatrix::multiply(double alpha)
{
  int n_entries = nRows_ * nColumns_;

  for(int i=0 ; i<n_entries; ++i)
  {
    this->entries_[i] *= alpha;
  }
}

/** ***************************************************************************/
std::vector<double> DenseMatrix::multiply(
                                      std::vector<double>::const_iterator first,
                                      std::vector<double>::const_iterator last,
                                      double alpha,
                                      bool transpose) const
{
  char trans[1] = {'n'};

  double beta  = 0.;

  int nRows_a = nRows_;
  int nCols_a = nColumns_;
  int LDA     = leadingDimension_;
  int size_x  = distance(first, last);
  int size_y  = nRows_a;
  int size_a  = nCols_a;
  int incx    = 1;
  int incy    = 1;

  if(transpose)
  {
     trans[0] = 't';
     size_y   = nCols_a;
     size_a   = nRows_a;
  }

  if(size_x != size_a)
  {
    ErrThrow("Dimension mismatch during matrix-vector multiplication: ",
             size_a, " columns in this->matrix and ",
             size_x, " elements in the vector.");
  }

  std::vector<double> y(size_y);

  dgemv_(trans,
         &nRows_a,
         &nCols_a,
         &alpha,
         entries_,
         &LDA,
         &(*first),
         &incx,
         &beta,
         &y[0],
         &incy);

   return y;
}

/** ***************************************************************************/
std::vector<double> DenseMatrix::multiply(const std::vector<double>* const x,
                                          double alpha,
                                          bool transpose) const
{
  return this->multiply(x->begin(), x->end(), alpha, transpose);
}

/** ***************************************************************************/
std::shared_ptr<DenseMatrix> DenseMatrix::multiply(const DenseMatrix* const B,
                                                   bool  transp_A,
                                                   bool  transp_B) const
{
  char transa[1] = {'n'};
  char transb[1] = {'n'};

  double alpha = 1.;
  double beta  = 0.;

  int nRows_a = nRows_;
  int nCols_a = nColumns_;
  int LDA     = leadingDimension_;
  int nRows_b = B->getNRows();
  int nCols_b = B->getNColumns();
  int LDB     = B->getLeadingDimension();
  int nRows_c = nRows_a;
  int nCols_c = nCols_b;

  if(transp_A)
  {
     transa[0] = 't';
     nRows_a = nCols_a;
     nCols_a = nRows_c;
     nRows_c = nRows_a;
  }

  if(transp_B)
  {
     transb[0] = 't';
     nCols_b = nRows_b;
     nRows_b = nCols_c;
     nCols_c = nCols_b;
  }

  int LDC = nRows_c;

  if(nCols_a != nRows_b)
  {
    ErrThrow("Dimension mismatch during matrix-matrix multiplication: ",
             nCols_a, " columns in this->matrix and ",
             nRows_b, " rows in the other one.");
  }

  auto c = std::make_shared<DenseMatrix>(nRows_c, nCols_c);

  const double* const b_entries = B->get_entries();
  double* c_entries = c->get_entries();

   dgemm_(transa,
          transb,
          &nRows_a,
          &nCols_b,
          &nCols_a,
          &alpha,
          entries_,
          &LDA,
          b_entries,
          &LDB,
          &beta,
          c_entries,
          &LDC);

   return c;
}

/** ***************************************************************************/
void DenseMatrix::decomposeLU()
{
  // get memory for LU fact if necessary and copy entries into it
  if(!entriesLU_)
  {
    entriesLU_ = new double[leadingDimension_*nColumns_]();
  }
  if(!pivotsLU_)
  {
    pivotsLU_= new int[nColumns_]();
  }

  memcpy(entriesLU_, entries_, leadingDimension_*nColumns_*sizeof(double));

  //Do the LU decomp and store in entries_ and pivotsLU_
  int info = 0;
  int nR = nRows_;
  int nC = nColumns_;
  int LDA = leadingDimension_;
  dgetrf_(&nR, &nC, entriesLU_, &LDA, pivotsLU_, &info);

  if(info != 0)
  {
    ErrThrow("LAPACK dgetrf (factorize) failed with info = ", info);
  }
}

/** ***************************************************************************/
void DenseMatrix::solve(double* rhsToSolution) const
{
  // Throw if the matrix is not quadratic.
  if(!(nRows_ == nColumns_ && nColumns_ == leadingDimension_))
  {
    throw std::runtime_error("DenseMatrix::solve: Solver is only operating when"
                           " nRows_ == nColumns_ == leadingDimension_ so far.");
  }

  //Solve, assuming the LU decomposition was performed before.
  if(!entriesLU_ || !pivotsLU_)
  {
    ErrThrow("DenseMatrix::decomposeLU must be called"
             " before DenseMatrix::solve");
  }

  // Call the solve routine.
  int info;
  char control[1] = {'n'};
  int nC = nColumns_;
  int LDA = leadingDimension_;
  int nrhs = 1;
  dgetrs_(control , &nC, &nrhs, entriesLU_, &LDA,
          pivotsLU_, rhsToSolution, &nC, &info);

  if(info != 0)
  {
    ErrThrow("LAPACK dgetrs (solve) failed with info = ", info);
  }
}

/** ***************************************************************************/
void DenseMatrix::reset()
{
  memset(entries_, 0., leadingDimension_*nColumns_* sizeof(double));
  this->resetLU();
}

/** ***************************************************************************/
void DenseMatrix::resetLU()
{
  if(pivotsLU_ != nullptr)
  {
    memset(pivotsLU_, 0, nColumns_*sizeof(int));
  }

  if(entriesLU_ != nullptr)
  {
    memset(entriesLU_, 0., leadingDimension_*nColumns_*sizeof(double));
  }
}

/** ***************************************************************************/
void DenseMatrix::clearLU()
{
  if(pivotsLU_ != nullptr)
  {
    delete[] pivotsLU_;
    pivotsLU_ = nullptr;
  }

  if(entriesLU_ != nullptr)
  {
    delete[] entriesLU_;
    entriesLU_ = nullptr;
  }
}

/** ***************************************************************************/
void DenseMatrix::setEntry(size_t r, size_t c, double value)
{
  // Throw if wrong line or column number is entered.
  assert(r < nRows_ && c < nColumns_);

  //Write the value to the right place.
  entries_[c*leadingDimension_ + r] = value;
}

/** ***************************************************************************/
double* DenseMatrix::getEntry_ptr(size_t r, size_t c)
{
  // Throw if wrong line or column number is entered.
  assert (r < nRows_ && c < nColumns_);

  //Get the pointer to the specified entry.
  return entries_+(c*leadingDimension_ + r);
}

/** ***************************************************************************/
double DenseMatrix::getEntry( size_t r, size_t c) const
{
  // Throw if wrong line or column number is entered.
  assert (r < nRows_ && c < nColumns_);

  //Get the entry from the right place.
  return entries_[c*leadingDimension_ + r];
}

void DenseMatrix::getRow(size_t row_index, double* buffer) const
{
  for (unsigned int i = 0; i < nColumns_; i++)
  {
    buffer[i] = getEntry(row_index, i);
  }
}

void DenseMatrix::getColumn(size_t col_index, double* buffer) const
{
  for (unsigned int i = 0; i < nRows_; i++)
  {
    buffer[i] = getEntry(i, col_index);
  }
}

void DenseMatrix::setRow(size_t row_index, const double* buffer)
{
  for (unsigned int i = 0; i < nColumns_; i++)
  {
    setEntry(row_index, i, buffer[i]);
  }
}

void DenseMatrix::setColumn(size_t col_index, const double* buffer)
{
  for (unsigned int i = 0; i < nRows_; i++)
  {
    setEntry(i, col_index, buffer[i]);
  }
}

/** ***************************************************************************/
double DenseMatrix::getLUEntry(const size_t r, size_t c) const
{
  if(!entriesLU_)
  {
    ErrThrow("You have to factorize before you can view the factorization!");
  }

  // Throw if wrong line or column number is entered.
  assert (r < nRows_ && c < nColumns_);

  //Get the entry from the right place.
  return entriesLU_[c*leadingDimension_ + r];
}

/** ***************************************************************************/
std::vector<double> DenseMatrix::get_matrix_row(size_t i) const
{
  if((int) i > (int)nRows_)
  {
    ErrThrow("Requested row number greater than number of rows.");
  }

  std::vector<double> row(nColumns_ , 0.0);
  for(int j = 0; j < (int)nColumns_ ; j++)
  {
    row.at(j) = entries_[j*leadingDimension_ + i];
  }
  return row;
}

/** ***************************************************************************/
std::vector<double> DenseMatrix::get_matrix_column(size_t j) const
{
  if((int) j > (int)nColumns_)
  {
    ErrThrow("Requested column number greater than number of columns.");
  }

  std::vector<double> col(nRows_, 0.0);
  for(int i =0; i < (int)nRows_; i++)
  {
    col.at(i) = entries_[j*leadingDimension_ + i];
  }
  return col;
}

/** ***************************************************************************/
void DenseMatrix::print(const std::string& name) const
{
  Output::info("DenseMatrix: ", name);
  for (size_t i = 0; i < nRows_ ;++i)
  {
    std::string row = "Row " + std::to_string(i) + " ";
    for(size_t j =0; j < nColumns_; ++j)
    {
      row += std::to_string(getEntry(i,j)) + " ";
    }
    Output::print(row);
  }
}

/** ***************************************************************************/
void DenseMatrix::printLU(const std::string& name) const
{
  assert(entriesLU_ && pivotsLU_);

  Output::info("DenseMatrix LU decomp: ", name);
  for (size_t i = 0; i < nRows_ ;++i)
  {
    std::string row = "Row " + std::to_string(i) + " ";
    for(size_t j =0; j < nColumns_; ++j)
    {
      row += std::to_string(getLUEntry(i,j)) + " ";
    }
    Output::print(row);
  }

  Output::info("DenseMatrix LU pivots: ", name);
  std::string lu;
  for(size_t j =0; j < nColumns_; ++j)
    lu += std::to_string(pivotsLU_[j]) + " ";
  Output::print(lu);
}

/** ***************************************************************************/
double DenseMatrix::norm() const
{
  double sum = 0;
  for(size_t k = 0; k < leadingDimension_*nColumns_ ;++k)
  {
    sum += entries_[k]*entries_[k];
  }
  return std::sqrt(sum);
}

/** ***************************************************************************/
double DenseMatrix::norm2(double tolerance) const
{
#ifdef _MPI
  ErrThrow("Spectral norm of dense matrix is not tested with MPI yet.");
#endif

  // compute A^T * A
  int n = nColumns_;
  char UPLO = 'u';
  char trans = 't';
  double alpha = 1.;
  double beta  = 0.;
  int nCols_a = nColumns_;
  int nRows_a = nRows_;
  int LDA     = leadingDimension_;
  auto AtA_entries = new double[n*n]();
  int LDAtA     = n;

  dsyrk_(&UPLO, &trans, &nCols_a, &nRows_a, &alpha, entries_, &LDA,
         &beta, AtA_entries, &LDAtA);

  // compute largest eigenvalue max_eig of A^T * A
  char job_type = 'N'; // Compute eigenvalues only
  char range = 'I';    // the Imin-th through Imax-th eigenvalues
  double Vmin = 0.;    // dummy value, not used with range = 'I'
  double Vmax = 0.;    // dummy value, not used with range = 'I'
  int Imin = n;        // in order to only compute the largest eigenvalue
  int Imax = n;        // since they are sorted in ascending order
  int m = 1;           // nb eigenvalues found
  double max_eig = 0.; // the largest eigenvalue
  double* Z = nullptr; // dummy value, not used with range = 'I'
  int LDZ = 1;
  int isuppz[2];
  int lwork;
  int liwork;
  int info = 0;

  /// workspace query: the routine only calculates the optimal size of work
  lwork = -1;
  liwork = -1;
  double wkopt;
  int iwkopt;

  dsyevr_(&job_type, &range, &UPLO, &n, AtA_entries, &LDAtA, &Vmin, &Vmax,
          &Imin, &Imax, &tolerance, &m, &max_eig, Z, &LDZ, isuppz,
          &wkopt, &lwork, &iwkopt, &liwork, &info);
  lwork = (int)wkopt;
  auto work = new double[lwork]();
  liwork = (int)iwkopt;
  auto iwork = new int[liwork]();

  /// solve eigenproblem, the eigenvalues are stored in ascending order
  dsyevr_(&job_type, &range, &UPLO, &n, AtA_entries, &LDAtA, &Vmin, &Vmax,
          &Imin, &Imax, &tolerance, &m, &max_eig, Z, &LDZ, isuppz,
          work, &lwork, iwork, &liwork, &info);

  /// check for convergence
  if( info > 0 )
  {
    ErrThrow("The algorithm failed to compute eigenvalues.");
  }

  delete[] work;
  delete[] iwork;
  delete[] AtA_entries;

  return std::sqrt(max_eig);
}

/** ***************************************************************************/
void DenseMatrix::norm2_bounds(double& min, double& max) const
{
  // lower bound: square root of the maximum column sum of squares
  double sum;
  min = 0.;
  for(size_t i=0 ; i<nColumns_ ; ++i)
  {
    sum = 0.;
    for(size_t j=0 ; j<leadingDimension_ ; ++j)
    {
      sum += entries_[i*leadingDimension_+j]*entries_[i*leadingDimension_+j];
    }
    if(sum>min)
    {
      min = sum;
    }
  }
  min = std::sqrt(min);

  // upper bound: Frobenius norm
  max = this->norm();
}

/** ***************************************************************************/
void swap(DenseMatrix& first, DenseMatrix& second)
{
  //std::swap all members
  std::swap(first.nRows_, second.nRows_);
  std::swap(first.nColumns_, second.nColumns_);
  std::swap(first.leadingDimension_, second.leadingDimension_);

  std::swap(first.entries_, second.entries_);
  std::swap(first.entriesLU_, second.entriesLU_);
  std::swap(first.pivotsLU_, second.pivotsLU_);
}
