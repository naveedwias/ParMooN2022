/** ************************************************************************
*
* @class     TMatrix
* @brief     store a sparse matrix
*
* Sparse matrices in ParMooN are stored in <em>compressed row storage</em>.
* Basically the matrix consists of a structure and entries and can perform
* algebraic operations.
*
* To create a TMatrix either copy or move an existing one or use a TStructure
* to construct a new one.
*
* Objects ob type TMatrix are not supposed to be copied. If you want that
* copy/move construct a new TMatrix.
*
* \todo apply naming convention to all methods
*
* @ruleof0
*
 ************************************************************************  */

#ifndef __MATRIX__
#define __MATRIX__

#include <DenseMatrix.h>
#include <Structure.h>
#include <string>
#include <map>
#include <vector>

#ifdef _MPI
class TParFECommunicator3D;
/**
 * In MPI case, different sparsity patterns may need adaptations in
 * other program parts (MumpsWrapper!)
 *
 * STANDARD: Is the usual case - non-zero entries M_{ij} occur only
 * when i and j are d.o.f. sharing the same cell.
 *
 * B_TIMES_BT: This sparsity pattern occurs, when a STANDARD matrix
 * B gets multiplied with its own transpose. The product matrix P
 * will have additional entries P_{ij} with i and j in adjacent cells,
 * which necessitates changes in the parallel structure of the code.
 */
enum class SparsityType{STANDARD, B_TIMES_BT};
#endif

class TMatrix
{
  protected:
    /** @brief structure of the matrix
     *
     * This stores the information where this matrix has entries. Many objects
     * of type TMatrix can have the same TStructure, therefore we use a
     * `std::shared_ptr` here.
     */
    std::shared_ptr<TStructure> structure;

    /** @brief matrix entries
     *
     * Its size is determined by the TMatrix::structure.
     */
    std::vector<double> entries;

#ifdef _MPI
    /// The sparsity type of the matrix, holding information for
    /// other program parts, esp. the solver.
    SparsityType sparse_type;
#endif

    /**
     * @brief replace the structure by a copy
     *
     * If you plan on calling a method which changes the structure, such as
     * TMatrix::remove_zeros, TMatrix::changeRows or TMatrix::reorderMatrix,
     * this method has to be called before. Otherwise other matrices sharing
     * the same TStructure will become corrupted.
     *
     * This method makes a deep copy of the structure. Then this matrix only
     * knows its copy and is the only one knowing this copy. Then this matrix
     * can safely call any structure changing methods.
     *
     * Note that this method does nothing, if this matrix is the only one
     * sharing its structure.
     */
    void copyOwnStructure();

  public:
    /** @brief generate the matrix, initialize entries with zeros */
    explicit TMatrix(std::shared_ptr<TStructure> structure
#ifdef _MPI
            , const SparsityType& sparse_type = SparsityType::STANDARD
#endif
    );

    /**
     * @brief Generates an empty `nRows`*`nCols` Matrix with no entries
     */
    TMatrix(int nRows, int nCols);

    /**
     * @brief Generates a diagonal matrix with entries from input std::vector<double> diag
     */
    TMatrix(const std::vector<double>& diag,
            std::shared_ptr<TStructure> Structure);

    /// @brief Default copy constructor
    ///
    /// Performs a deep copy of the entries, but not of the structure.
    TMatrix(const TMatrix&) = default;

    /// @brief Default move constructor.
    TMatrix(TMatrix&&) = default;

    /// @brief copy assignment, deep copy of the entries, shallow copy of the
    /// structure
    TMatrix & operator=(const TMatrix& A) = default;

    /// @brief Default move assignment
    TMatrix& operator=(TMatrix&&) = default;

    /// @brief Default destructor.
    virtual ~TMatrix() = default;


    /// @brief reset all matrix entries to zero
    void reset();

    int clean_denormals();

    /// @brief return number of rows
    int get_n_rows() const
    { return structure->get_n_rows(); }

    /// @brief return number of columns
    int get_n_columns() const
    { return structure->get_n_columns(); }

    /// @brief return number of matrix entries
    int get_n_entries() const
    { return structure->get_n_entries(); }

    /// @brief return number of entries in a specified row
    size_t get_n_entries_in_row(size_t row_index) const
    { return structure->get_n_entries_in_row(row_index); }

    /// @brief return the column pointer in the TStructure of this matrix
    const int *get_vector_columns() const
    { return structure->get_vector_columns(); }
    /** @brief return the column pointer in the TStructure of this matrix
     *
     * This version should never be used. It only exists because some other
     * parts of the software do not respect the const keyword (like AMG) or are
     * not well implemented (changing the structure of a matrix).
     */
    int *get_vector_columns()
    { return structure->get_vector_columns(); }

    /// @brief return the row pointer in the TStructure of this matrix
    const int *get_row_ptr() const
    { return structure->get_row_ptr(); }

    /** @brief return the row pointer in the TStructure of this matrix
     *
     * This version should never be used. It only exists because some other
     * parts of the software do not respect the const keyword (like AMG) or are
     * not well implemented (changing the structure of a matrix).
     */
    int *get_row_ptr()
    { return structure->get_row_ptr(); }

    const std::vector<int>& get_row_array() const
    { return structure->get_row_array(); }

    /// @brief return structure
    const TStructure& GetStructure() const
    { return *structure; }

    /// @brief return matrix entries as a pointer to const double
    /// @note try to avoid this method and use TMatrix::get_entries() instead
    const double *GetEntries() const
    { return &entries[0]; }

    /// @brief return matrix entries as a pointer to double
    double *GetEntries()
    { return &entries[0]; }

    /// @brief return the entries as a std::vector
    const std::vector<double>& get_entries() const
    { return entries; }

    /// @brief get a vector of the diagonal entries of this matrix.
    /// This will work even if the structure does not have entries on the
    /// diagonal. In that case it will have a zero there. The length of this
    /// vector is the minimum of the number of rows and columns.
    std::vector<double> get_diagonal() const;

    std::vector<double> get_row_sums() const;
    std::vector<double> get_col_sums() const;

    /** @brief return the norm of the matrix
     *
     * The parameter \p p determines which norm to compute. Choose \p as
     * -2 for Frobenius norm
     * -1 for maximum absolute row sum
     *  0 for maximum entry
     *  1 for maximum absolute column sum (not yet implemented)
     *  2 for euclidean norm, (not yet implemented)
     */
    double GetNorm(int p=-1) const;

    /** @brief write matrix into file in MatrixMarket format
     *
     * Writes the matrix into the file whose path and name is given as input
     * parameter.
     * So far there is only one format: the MatrixMarket coordinate format.
     * It is human readable and explained at
     *  http://math.nist.gov/MatrixMarket/formats.html.
     * There is a nice MATLAB read-in function (mmread.m) for that format on
     * the same homepage. You can use "write" in combination with that function
     * to examine your TMatrix in MATLAB, which is great for debugging.
     *
     * @param[in] filename The filename where to write the matrix.
     */
    void write(const std::string& filename) const;

    /// @brief Print matrix into the shell
    void Print(const char *name = "a") const;

    /** @brief print the full matrix, including all zeros
     *
     * This is only meaningful for very small matrices.
     */
    void PrintFull(const std::string& name="", int fieldWidth=4) const;

    /// @brief add a value at selected entry
    void add(int i,int j, double val);
    /** @brief add values in row 'i' given by the map 'vals', multiplied by
     * 'factor'
     *
     * This should be faster than adding all values in 'vals' individually
     */
    void add(int i, const std::map<int,double>& vals, double factor = 1.0);
    /** @brief add values `vals[i][j]` to this matrix at the positions `(i,j)`
     * for all `i,j` defined in the map `vals`
     */
    void add(const std::map<int, std::map<int,double>>& vals, double factor = 1.0);
    /// @brief set a value at selected entry
    void set(int i, int j, double val);
    /// @brief get a value at selected entry
    const double& get(int i,int j) const;
    /// @brief  get a value at selected entry (you may change that value)
    double& get(int i,int j);

    /** This method will return a column of this CSR matrix as a vector.
     *@param[in] j The column index.
     *@return A matrix column as an std::vector<double>.
     */
    std::vector<double> get_matrix_column(size_t j) const;

    /** This method will return a row of a CSR matrix as a vector.
     *@param[in] i The row index.
     *@return A matrix row as an std::vector.
     */
    std::vector<double> get_matrix_row(size_t i) const;

    /**
     *
     * This method will replace a column of this matrix by the given vector.
     * The program will quit if the vector dfoes not fit the sparsity structure
     * of the matrix or if the given column index is out of range.
     *
     *@param[in] j The column index.
     *@return A matrix column as an std::vector<double>.
     */
    void set_matrix_column(size_t j, const std::vector<double>& col);

    /// @brief reset the entries, the given vector must be of the same size
    void setEntries(const std::vector<double>& entries);

    /** @brief reorders the Matrix to comply with direct solvers.
     *
     * @warning This changes the structure of the matrix
     */
    void reorderMatrix();

    /////////////// Routines for periodic boundary conditions /////////////////
    /**
    * @brief replace several rows in the matrix with new entries.
    *
    * Replace rows by new ones. This creates a new structure for the sparsity
    * pattern of the matrix. Therefore reallocation is necessary.
    *
    * If there are no rows to change, i.e. if entries.size()==0, nothing is·
    * done.
    *·
    * This will create a new structure for this matrix. The structure·
    * previously belonging to this matrix is not changed. So other matrices·
    * are not affected.
    *
    * @param entries for every row a map of columns-to-entries map
    */
    void changeRows(const std::map<int,std::map<int,double>>& entries);
    ///////////////// ///////////////// ///////////////// /////////////////

    /** @brief return a new TMatrix which is the transposed of this matrix
     *
     * If this is an object of a derived class (e.g. TMatrix2D, TSquareMatrix),
     * then the number of active degrees of freedom is not taken into account.
     * The returned TMatrix is really the algebraic transposed matrix.
     * */
    TMatrix* get_transposed() const;

    /** @brief compute y = A*x   (Matrix-Vector-Multiplication)
     *
     * Note that 'y' is created here and it is up to the user to delete it.
     */
    friend double* operator*(const TMatrix & A, const double* x);

    /** @brief add another matrix to this one
     *
     * This is of course only possible if the corresponding structures are the
     * same.
     */
    virtual TMatrix & operator+=(const TMatrix * A);
    /** @brief substract another matrix to this one
     *
     * This is of course only possible if the corresponding structures are the
     * same.
     */
    virtual TMatrix & operator-=(const TMatrix * A);
    /** @brief add another matrix to this one
     *
     * This is of course only possible if the corresponding structures are the
     * same. This method exists only for convenience and uses the same method
     * with a pointer to A instead of the reference.
     */
    virtual TMatrix & operator+=(const TMatrix & A)
    { *this += &A; return *this; }
    /** @brief scale matrix by a factor */
    virtual TMatrix & operator*=(const double a);

    /** @brief compute y += a * A*x
     *
     * 'A' is this TMatrix and 'x', and 'y' are given vectors. The scalar 'a'
     * is a scaling factor.
     *
     * @param x array representing the vector which is multiplied by this matrix
     * @param y array representing the vector to which a * A*x is added
     * @param a scaling factor, default is 1.0
     */
    void multiply(const double * const x, double * const y, double a = 1.0) const;

    /** @brief compute y += a * A^T * x
     *
     * 'A^T' is the transposed of this TMatrix and 'x', and 'y' are given
     * vectors. The scalar 'a' is a scaling factor.
     *
     * @param x array representing the vector which is multiplied by this matrix
     * @param y array representing the vector to which a * A*x is added
     * @param a scaling factor, default is 1.0
     */
    void transpose_multiply(const double * const x, double * const y, double a = 1.0)
      const;

    /**
     * @brief compute matrix-matrix product C = a*A*B,
     *
     * 'A' is this matrix, 'a' is a scalar factor, 'B' is given. Then matrix
     * 'C' is created during this function and the user is responsible to
     * delete C.
     *
     * Note that this is rather slow.
     *
     * @param B matrix to be multiplied (from right) to this matrix
     * @param a scaling factor, default is 1.0
     */
    TMatrix* multiply(const TMatrix * const B, double a = 1.0) const;

    /**
     * @brief compute matrix-matrix product C = a * A * op(B),
     *
     * 'A' is this matrix, 'a' is a scalar factor, 'B' is a given DenseMatrix,
     * op(B) = B if transpose is false, or else op(B) = B^T (transpose B).
     * The DenseMatrix 'C' is created during this function.
     *
     * @param[in] B         matrix to be multiplied (from right) to this matrix
     * @param[in] transpose boolean which specifies the form of op(B),
     *                      default is false
     * @param[in] a         scaling factor, default is 1.0
     * @return A pointer to the product matrix a * A * op(B).
     */
    std::shared_ptr<DenseMatrix> multiply(const DenseMatrix* const B,
                                          bool transpose = false,
                                          double a = 1.0) const;

    /**
        * @brief compute matrix-matrix product C = A * diag[d] * B,
        *
        * 'A' is this matrix, 'd' is a vector, 'B' is given. Then matrix
        * 'C' is created during this function and the user is responsible to
        * delete C.
        *
        * Note that this is rather slow.
        *
        * @param B matrix to be multiplied (from right) to this matrix
        * @param d vector
        */
       TMatrix* multiply(const TMatrix * const B,
                         const std::vector<double>& d) const;

    /**
     * @brief multiply this matrix B with its transposed B^T from the right
     *
     * @return A pointer to the product matrix B*B^T.
     */
    TMatrix* multiply_with_transpose_from_right(
#ifdef _MPI
      const std::vector<const TParFECommunicator3D*>& test_comms,
      const std::vector<const TParFECommunicator3D*>& ansatz_comms,
      bool additive_storage = true
#endif
    ) const;

    /** @brief multiply this matrix B with its transposed B^T from the right
     * and scale with a diagonal matrix in between.
     *
     * Compute the product B*D*B^T where B is this matrix and D is a diagonal
     * matrix.
     *
     * @param[in] diagonalScaling The diagonal scaling matrix D as a vector.
     * Must have as many entries as B has columns.
     * @return A pointer to the product matrix B*D*B^T.
     */
    TMatrix* multiply_with_transpose_from_right(
        const std::vector<double>& diagonalScaling
#ifdef _MPI
        , const std::vector<const TParFECommunicator3D*>& test_comms,
        const std::vector<const TParFECommunicator3D*>& ansatz_comms,
        bool additive_storage = true
#endif
    ) const;

    /** @brief multiply this matrix B with its transposed B^T from the right
     * and scale with a diagonal matrix in between.
     *
     * Computing the structure of the product is computationally most
     * intensive. If it is known from similiar former computations, pass it as
     * an argument.
     *
     * @param[in] structure The known structure.
     * @param[in] diagonalScaling The diagonal scaling matrix D as a vector.
     * Must have as many entries as B has columns.
     * @return A pointer to the product matrix B*D*B^T.
     *
     * @note in MPI mode: the resulting matrix is by default in additive
     * storage. If 'additive_storage' is set to false, the resulting matrix is
     * not correctly (ie, inconsistently) stored and this is your risk.

     */
    TMatrix* multiply_with_transpose_from_right(
        const std::vector<double>& diagonalScaling, const TStructure& structure
#ifdef _MPI
        , const std::vector<const TParFECommunicator3D*>& test_comms
        , const std::vector<const TParFECommunicator3D*>& ansatz_comms
        , bool additive_storage = true
#endif
    )
    const;

    /** @brief multiply this matrix A with its transposed A^T from the right
     * and multiply with matrix B in between.
     *
     * @param[in] B the matrix to be multiplied in between
     * @return A pointer to the product matrix A*B*A^T.
     */
    std::shared_ptr< TMatrix > multiply_with_transpose_from_right(
        const TMatrix& B
#ifdef _MPI
        , const std::vector<const TParFECommunicator3D*>& test_comms
        , const std::vector<const TParFECommunicator3D*>& ansatz_comms
#endif
    ) const;

    /**
     * @brief multiply this matrix B with B', structured like its transpose
     * B^T, from the right
     *
     * @return A pointer to the product matrix B B'.
     */
    TMatrix* multiply_with_pseudotranspose_from_right(
      const TMatrix& BT
#ifdef _MPI
      , const std::vector<const TParFECommunicator3D*>& test_comms
      , const std::vector<const TParFECommunicator3D*>& ansatz_comms
      , bool additive_storage = true
#endif
    ) const;

    /** @brief multiply this matrix B with B', structured like its transpose
     * B^T, from the right and scale with a diagonal matrix in between.
     *
     * Compute the product B D B' where B is this matrix, D is a diagonal
     * matrix, and B' is structured like B^T.
     *
     * @param[in] diagonalScaling The diagonal scaling matrix D as a vector.
     * Must have as many entries as B has columns.
     *
     * @return A pointer to the product matrix B D B'.
     */
    TMatrix* multiply_with_pseudotranspose_from_right(
      const TMatrix& BT,
      const std::vector<double>& diagonalScaling
#ifdef _MPI
      , const std::vector<const TParFECommunicator3D*>& test_comms,
      const std::vector<const TParFECommunicator3D*>& ansatz_comms,
      bool additive_storage = true
#endif
    ) const;

    /** @brief multiply this matrix B with B', structured like its transpose
     * B^T, from the right and scale with a diagonal matrix in between.
     *
     * Computing the structure of the product is computationally most
     * intensive. If it is known from similiar former computations, pass it as
     * an argument.
     *
     * @param[in] structure The known structure.
     * @param[in] diagonalScaling The diagonal scaling matrix D as a vector.
     * Must have as many entries as B has columns.
     *
     * @return A pointer to the product matrix B D B'.
     *
     * @note in MPI mode: the resulting matrix is by default in additive
     * storage. If 'additive_storage' is set to false, the resulting matrix is
     * not correctly (ie, inconsistently) stored and this is your risk.

     */
    TMatrix* multiply_with_pseudotranspose_from_right(
      const TMatrix& BT,
      const std::vector<double>& diagonalScaling, const TStructure& structure
#ifdef _MPI
      , const std::vector<const TParFECommunicator3D*>& test_comms,
      const std::vector<const TParFECommunicator3D*>& ansatz_comms,
      bool additive_storage = true
#endif
    )
    const;

    /** @brief multiply this matrix A with A', structured like its transpose
     * A^T, from the right and multiply with matrix B in between.
     *
     * @param[in] B the matrix to be multiplied in between
     * @return A pointer to the product matrix A B A'.
     */
    std::shared_ptr< TMatrix > multiply_with_pseudotranspose_from_right(
      const TMatrix& AT,
      const TMatrix& B
#ifdef _MPI
      , const std::vector<const TParFECommunicator3D*>& test_comms,
      const std::vector<const TParFECommunicator3D*>& ansatz_comms
#endif
    ) const;

    /// @brief perform one successive overrelaxation (sor) sweep.
    /// The flag can be either 0(forward sweep), 1(backward sweep), or
    /// 2(forward followed by backward sweep).
    /// @param[in] b right hand side
    /// @param[in,out] x solution (this is updated)
    /// @param[in] omega relaxation parameter
    /// @param[in] flag either 0 (forward), 1(backward), or 2(both)
    /// @param[in] par_strat The parallel strategy (MPI only). Use "all_cells" or
    ///            "own_cells".
    /// @param[in] comm MPI FE communicator which belongs to this matrix - MPI only!
    /// (might even be stored therein, if it actually is an FEMatrix)
    void sor_sweep(const double* b, double* x, double omega, size_t flag
#ifdef _MPI
                   , const std::string& par_strat
                   , const TParFECommunicator3D& comm
#endif
    ) const;

    /** @brief adding a scaled matrix to this matrix
     *
     * The summation is index-wise, i.e. A(i,j) += factor*m(i.j), where A is
     * this matrix.
     *
     * Note that this only works if the sparsity structure is the same for this
     * matrix and m.
     */
    void add_scaled(const TMatrix& m, double factor = 1.0);

#ifdef _MPI
    void add_scaled(const TMatrix& m, double factor,
        const TParFECommunicator3D &comm, bool masters_only);
#endif

    /// @brief scale all entries of this matrix
    void scale(double factor);

    /**
     * @brief scale a matrix using a vector
     *
     * think of this as multipling this matrix with a diagonal matrix from the
     * left (if the second argument is true). The parameter 'factor' are the
     * entries of that diagonal matrix.
     *
     * The i-th row (colum if from_left is false) is scaled by the i-th entry in
     * 'factor'.
     *
     * The array 'factor' must be of size this->get_n_rows() if 'from_left' is
     * true or of size this->get_n_columns() otherwise.
     *
     * If all entries in 'factor' are the same, you can use operator*= as well.
     *
     * @param factor array of scaling factors
     * @param from_left scale rows (true) or columns (false)
     */
    void scale(const double * const factor, bool from_left = true);

    /**
     * @brief remove all entries from sparsity structure where a zero is stored
     *
     * This changes the sparsity structure of this matrix. Afterwards all stored
     * entries are nonzero. This can help if a lot of zeros are explicitly
     * stored in the matrix. Afterwards matrix operations should be faster.
     *
     * if tol is greater or equal to zero, all entries with magnitude smaller
     * than tol are removed. If tol is smaller than zero, tol is set such that
     * the ratio of the largest over the smallest entry in the resulting matrix
     * is smaller than 10^{15}.
     */
    void remove_zeros(double tol = 0.0);

    /** @brief get/set a specific matrix entry
     *
     * This will give an error if that entry is not in the sparsity structure
     */
    double & operator()(const int i, const int j);
    /** @brief get a specific matrix entry
     *
     * This will give an error if that entry is not in the sparsity structure
     */
    const double & operator()(const int i, const int j) const;

    /** @brief shift all indices in the structure of this matrix by 1
     *
     * See also the documentation of TStructure::fortran_shift().
     *
     * Note that this affects all matrices sharing this structure. This is only
     * needed to call fortran code. Don't use it in other situation as it may
     * break the functionality of other functions.
     */
    void fortran_shift()
    { this->structure->fortran_shift(); }

    /// @brief print some information on this TMatrix
    void info(size_t verbose) const;

    ///! @return true if this matrix is square
    bool is_square() const
    {
      return structure->is_square();
    }

#ifdef _MPI
    /// Get the sparsity type of the matrix.
    SparsityType get_sparse_type() const{return sparse_type;};
#endif
};

#endif
