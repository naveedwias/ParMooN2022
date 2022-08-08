/** **************************************************************************** 
*
* @name       DenseMatrix
* @brief      Provides a handy wrapper for general dense band matrices stored
*             as packed matrix in column-major order.
*
* Note that row and column indexing starts at 0 (C++ style), as in the rest of
* ParMooN.
*
* It is not suited to do anything but being set, decomposed and solved.
* Should you require more functionality, you'll have to take the trouble and
*
* @author     Clemens Bartsch
* @date       2015/06/03, import to ParMooN: 2016/06/21
*
*******************************************************************************/

#ifndef DENSEMATRIX_H_
#define DENSEMATRIX_H_

#include <cstddef>
#include <memory> // std::shared_ptr
#include <string>
#include <vector>

class DenseMatrix
{
public:
  /*!
   * @brief Construct a nRows x nColumns matrix filled with 0.
   * @param[in] nRows The number of rows.
   * @param[in] nColumns The number of columns.
   */
  DenseMatrix(size_t nRows, size_t nColumns);

  /*!
   * @brief Construct a nRows x nColumns matrix filled with 0
   * and manually set leading dimension.
   * @param[in] nRows The number of rows.
   * @param[in] nColumns The number of columns.
   */
  DenseMatrix(size_t nRows, size_t nColumns, size_t leadingDimension);

  //! Copy constructor. Performs deep copy.
  DenseMatrix(const DenseMatrix &obj);

  //! Move constructor.
  DenseMatrix(DenseMatrix&&);

  //! Assignment operator, using copy-and-swap idiom for deep copy.
  DenseMatrix& operator=(DenseMatrix other);

  //! Standard destructor.
  ~DenseMatrix();

  /** @brief add another matrix to this one
   *
   * This is of course only possible if the sizes are the same.
   */
  DenseMatrix& operator+=(const DenseMatrix& B);

  friend DenseMatrix operator+(DenseMatrix A, const DenseMatrix& B)
  {
    A += B;
    return A;
  }

  /** @brief substract another matrix to this one
   *
   * This is of course only possible if the sizes are the same.
   */
  DenseMatrix& operator-=(const DenseMatrix& B);

  friend DenseMatrix operator-(DenseMatrix A, const DenseMatrix& B)
  {
    A -= B;
    return A;
  }

  /**
   * @brief compute matrix-vector addition A += a * op(v)
   *
   * if transpose is to false, the addition is applied column-wise,
   * else op(v) is the transpose of the vector v and the addition is applied
   * row-wise
   *
   * @param[in] v         vector to be added to this matrix
   * @param[in] a         scaling factor, default is 1.0
   * @param[in] transpose boolean which specifies the form of op(v),
   *                      default false
   */
  void add(const std::vector<double>* const v,
           double a = 1., bool transpose = false);

  /**
   * @brief compute matrix-matrix addition A += a * B
   *
   * @param[in] B         matrix to be added to this matrix
   * @param[in] a         scaling factor, default is 1.0
   */
  void add(const DenseMatrix& B, double a = 1.);

  void addToColumn(const std::vector<double>* const v, int col_index,
           double a = 1);

  void addToRow(const std::vector<double>* const v, int row_index,
           double a = 1);

  void scaleColumn(int col_index, double a);

  void scaleRow(int row_index, double a);

  /**
   * @brief compute matrix-scalar product A = alpha * A
   *
   * @param[in] alpha scaling factor
   */
  void multiply(double alpha);

  /**
   * @brief compute matrix-vector product y = alpha * op(A) * x
   *
   * op(A) is this matrix if transpose is false, or else op(A) = A^T.
   *
   * iterator pointing to the first element
   * @param[in] first     iterator pointing to the fisrt element of the vector x
   *                      (to be multiplied (from right) to this matrix)
   * @param[in] last      iterator pointing to the last element of the vector x
   * @param[in] alpha     scale coefficient
   * @param[in] transpose boolean which specifies the form of op(A),
   *                      default false
   * @return A pointer to the product vector alpha * op(A) * x.
   */
  std::vector<double> multiply(std::vector<double>::const_iterator first,
                               std::vector<double>::const_iterator last,
                               double alpha = 1.,
                               bool  transpose = false) const;

  /**
   * @brief compute matrix-vector product y = alpha * op(A) * x
   *
   * op(A) is this matrix if transpose is false, or else op(A) = A^T.
   *
   * @param[in] x         vector to be multiplied (from right) to this matrix
   * @param[in] alpha     scale coefficient   * 
   * @param[in] transpose boolean which specifies the form of op(A),
   *                      default false
   * @return A pointer to the product vector alpha * op(A) * x.
   */
  std::vector<double> multiply(const std::vector<double>* const x,
                               double alpha = 1.,
                               bool  transpose = false) const;

  /**
   * @brief compute matrix-matrix product C = op(A) * op(B)
   *
   * op(A) is this matrix if transpose is false, or else op(A) = A^T.
   * op(B) is this matrix if transpose is false, or else op(B) = B^T.
   *
   * The DenseMatrix 'C' is created in this function.
   *
   * @param[in] B        matrix to be multiplied (from right) to this matrix
   * @param[in] transp_A boolean which specifies the form of op(A),
   *                     default false
   * @param[in] transp_B boolean which specifies the form of op(B),
   *                     default false
   * @return A pointer to the product matrix op(A) * B.
   */
  std::shared_ptr<DenseMatrix> multiply(const DenseMatrix* const B,
                                        bool  transp_A = false,
                                        bool  transp_B = false) const;

  /**
   * @brief Do and store the LU decomposition with LAPACK's dgetrf routine.
   */
  void decomposeLU();

  /*!
   * @brief Solve a linear system with LAPACK's dgetrs routine.
   *
   * Must make sure that decomposeLU has been called before - the method
   * will complain if not so.
   * The memory space for the in/out array and the correct array length must be
   * taken care of outside the class.
   *
   * @param[in,out] rhsToSolution An array which contains the rhs on input and
   *                              the solution calculated by FortranStyle on
   *                              output.
   */
  void solve(double* rhsToSolution) const;

  /// @brief reset all matrix entries to zero
  void reset();

  /// @brief reset all LU entries and pivots to zero
  void resetLU();

  /// @brief reset all LU pointers to nullptr
  void clearLU();

  //! @brief return the number of rows of this matrix
  int getNRows() const
  { return nRows_; }

  //! @brief return the number of columns of this matrix
  int getNColumns() const
  { return nColumns_; }

  //! @brief return the leading dimension of this matrix
  int getLeadingDimension() const
  { return leadingDimension_; }

  //! @brief Write a value into the matrix.
  void setEntry(size_t rowNumber, size_t columnNumber, double valueToWrite);

  //! @brief returns a pointer to the entries_ array (const version).
  const double* get_entries() const
  { return entries_; }

  //! @brief returns a pointer to the entries_ array.
  double* get_entries()
  { return entries_; }

  //! @brief returns a pointer to the value at a certain position.
  double* getEntry_ptr(size_t rowNumber, size_t columnNumber);

  /// Get a value at a certain position-
  /// this is merely used for debugging and testing.
  double getEntry(size_t rowNumber, size_t columnNumber) const;

  void getRow(size_t row_index, double* buffer) const;
  void getColumn(size_t col_index, double* buffer) const;

  void setRow(size_t row_index, const double* buffer);
  void setColumn(size_t col_index, const double* buffer);

  //! @brief returns a pointer to the entriesLU_ array (const version).
  const double* getEntriesLU() const
  { return entriesLU_; }

  /// Get a value at a certain position from the LU matrix -
  /// this is merely used for debugging and testing.
  double getLUEntry(const size_t r, size_t c) const;

  //! @brief returns a pointer to the pivotsLU_ array (const version).
  const int* getPivotsLU() const
  { return pivotsLU_; }

  /**
   * @brief This method will return a row of this matrix as a vector.
   *
   * @param[in] i the row index.
   * @return A matrix row as an std::vector<double>.
   */
  std::vector<double> get_matrix_row(size_t i) const;

  /**
   * @brief This method will return a column of this matrix as a vector.
   *
   * @param[in] j the column index.
   * @return A matrix column as an std::vector<double>.
   */
  std::vector<double> get_matrix_column(size_t j) const;

  /// Print the matrix the console.
  /// The method is not perfect and is intended for testing and debugging only.
  void print(const std::string& name) const;

  void printLU(const std::string& name) const;

  /// Calculates the matrix' Frobenius norm. We implement only Frobenius,
  /// because this one is independent on d.o.f. numbering and therefore
  /// handy for debugging in MPI.
  double norm() const;

  /// Calculates the spectral norm of the matrix A: sqrt(lambda_max)
  /// with lambda_max the largest eigenvalue of A^T * A
  double norm2(double tolerance=1e-12) const;

  /// Calculates lower and upper bounds for the spectral norm of the matrix A:
  ///  - lower bound is the square root of the maximum column sum of squares
  ///  - upper bound is the Frobenius norm computed with norm()
  void norm2_bounds(double& min, double& max) const;


private:
  friend void swap(DenseMatrix& one, DenseMatrix& other);

  //! The number of rows.
  size_t nRows_;

  //! The number of columns.
  size_t nColumns_;

  /*! The leading dimension, for a detailed explanation see e.g.:
  * http://www-01.ibm.com/support/knowledgecenter/SSFHY8_5.3.0/com.ibm.cluster.essl.v5r3.essl100.doc/am5gr_leaddi.htm
  */
  size_t leadingDimension_;

  /*!
   * An array which contains the entries of the matrix columnwise (FORTRAN!).
   * Each column consists of leadingDimension_ subsequent locations in the
   * array, but only the first nRows_ locations of these are used.
   * That is why we require leadingDimension >= nRows_ (usually =).
   */
  double* entries_;

  /// This should store the entries of the LU factorization, if computed.
  double* entriesLU_;

  //! The LAPACK LU decomposition calculates and stores pivot entries,
  //! which have to be passed to the solver.
  int* pivotsLU_;

};

#endif /* DENSEMATRIX_H_ */
