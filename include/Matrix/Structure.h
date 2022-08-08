/**
 * @class TStructure
 * @brief represent the sparsity structure of a matrix
 * 
 * Matrices in ParMooN are stored in the compressed row storage (CRS) format.
 * The structure essentially stores the number of rows (TStructure::nRows), the 
 * number of columns (TStructure::nColumns), the number of entries 
 * (TStructure::nEntries) and two vectors (TStructure::rows and 
 * TStructure::columns). 
 * 
 * The first vector TStructure::rows has length TStructure::nRows+1, indicating
 * how many entries there are in all previous rows. It always is
 *     TStructure::rows[0] = 0;
 *     TStructure::rows[TStructure::nRows] = TStructure::nEntries;
 * The number of entries in a given row `i` is given by the difference 
 *     unsigned int nEntriesInRow = TStructure::rows[i+1] - TStructure::rows[i];
 * 
 * The second vector TStructure::columns has length TStructure::nEntries and 
 * represents the column of each entry.  
 * 
 * To get the `j`-th entry in the `i`-th row (not the `j`-th entry overall) it
 * is
 *     unsigned int indexOfEntry = TStructure::rows[i] + j; 
 *     unsigned int columnOfEntry = TStructure::columns[indexOfEntry];
 * 
 * Now the entry (`i`, `columnOfEntry`) in a matrix with this structure is the
 * `indexOfEntry`-th entry in the list of entries. 
 * 
 * \todo apply naming convention to all methods
 */

#ifndef __STRUCTURE__
#define __STRUCTURE__

// forward declarations
class TFESpace;

#include <memory> // std::shared_ptr
#include <vector>

class TStructure
{
  public:
    /** @brief default constructor sets everything to 0/nullptr */
    TStructure();
    
    /// @name construct square structures using one finite element space
    /// @brief ansatz and test space is the same
    //@{
    explicit TStructure(std::shared_ptr<const TFESpace> space,
                        bool is_empty = false, bool face_integrals = false);
    //@}
    
    /// @name construct rectangular structures using two finite element spaces
    /// @brief test and ansatz space are possibly different
    ///
    /// A matrix using this structure represents a linear map from the ansatz 
    /// to the test space.
    /// @param[in] is_empty If true, the structure will not have any entries.
    /// A matrix which owns this structure thus represents the zero-map between
    /// two FE spaces.
    /// @param[in] face_integrals If true create a structure with couplings of
    /// dofs which in turn couple (in the usual sense) with a common other dof.
    //@{
    TStructure(std::shared_ptr<const TFESpace> testspace,
               std::shared_ptr<const TFESpace> ansatzspace,
               bool is_empty = false, bool face_integrals = false);
    //@}
    

    /**
     * @brief generate a square structure, all arrays are already defined
     * 
     * A deep copy of the arrays pointed to by \p col_ptr and \p row_ptr is 
     * done.
     * 
     * @param n number of rows/columns
     * @param N_entries number of entries in this structure
     * @param col_ptr the new TStructure::columns vector
     * @param row_ptr the new TStructure::rows vector
     */
    TStructure(int n, int N_entries, int *col_ptr, int *row_ptr);

    /**
     * @brief generate a structure, all arrays are already defined
     * 
     * Note that a deep copy of the arrays is performed, TStructure will not
     * take ownership of the ipnut arrays.
     * 
     * @param nRows number of rows
     * @param nCols number of columns
     * @param N_entries number of entries in this structure
     * @param col_ptr the new TStructure::columns vector
     * @param row_ptr the new TStructure::rows vector
     */
    TStructure(int nRows, int nCols, int N_entries, int *col_ptr, int *row_ptr);
    
    /** @brief Generates an empty `n`*`n` Structure with no entries */
    explicit TStructure(int n);
    
    /** @brief Generates an empty `nRows`*`nCols` Structure with no entries
     */
    TStructure(int nRows, int nCols);

    
    /** @brief Default copy constructor */
    TStructure(const TStructure&) = default;

    /// @brief Default move constructor.
    TStructure(TStructure&&) = default;

    /// @brief no copy assignment operator to avoid accidental copies
    TStructure & operator=(const TStructure& A) = delete;
    
    /// @brief no move assignment operator to avoid accidental moves
    TStructure& operator=(TStructure&&) = delete;

    /// @brief Default destructor.
    ~TStructure() = default;
    
    
    /// @brief return if this structure is square
    bool is_square() const
    { return nRows == nColumns; }

    /** @brief return number of rows */
    unsigned int get_n_rows() const
    { return nRows; }

    /** @brief return number of columns */
    unsigned int get_n_columns() const
    { return nColumns; }
    
    /** @brief return number of entries */
    unsigned int get_n_entries() const
    { return nEntries; }
    
    /** @brief return number of entries in a specified row */
    size_t get_n_entries_in_row(size_t row_index) const;
    
    /** @brief return the number of entries in the rows with indices smaller
     * than `row_index`. */
    size_t get_n_entries_up_to_row(size_t row_index) const;

    /** @brief return vector columns */
    const int *get_vector_columns() const
    { return &columns[0]; }
    /** @brief return vector columns
     * 
     * This version should never be used. It only exists because some other
     * parts of the software do not respect the const keyword (like AMG) or are
     * not well implemented (changing the structure of a matrix).
     *
     * \todo Get rid of this, as it's unsafe.
     */
    int *get_vector_columns()
    { return &columns[0]; }
    
    /** @brief return a copy of the columns vector */
    std::vector<int> get_columns() const
    { return columns; }

    /** @brief return array row pointer */
    const int *get_row_ptr() const
    { return &rows[0]; }

    /** @brief return array row pointer
     *
     * This version should never be used. It only exists because some other
     * parts of the software do not respect the const keyword (like AMG) or are
     * not well implemented (changing the structure of a matrix).
     *
     * \todo Get rid of this, as it's unsafe.
     */
    int *get_row_ptr()
    { return &rows[0]; }
    
    /** @brief return array row pointer */
    const std::vector<int>& get_row_array() const
    { return rows; }

    /**
     * Resets the number of entries to fit the row and column array. Will check
     * if the size values are in accordance. Must be called after the row array
     * was changed from the external.
     *
     */
    void reset_n_entries();

    /**
     * @brief find the index of a given entry
     * 
     * If the (i,j)-th entry is not in the sparsity pattern, -1 is returned. 
     * This is how this function can be used to check whether an entry is in the
     * sparsity pattern.
     * 
     * @param i row of entry to check
     * @param j column of entry to check
     */ 
    int get_index_of_entry(const int i, const int j) const;
    
    /** @brief shift all indices by 1
     * 
     * In fortran indices start with 1 instead of 0. Calling this function 
     * changes from one to the other, so calling it twice does not change 
     * anything.
     */
    void fortran_shift();
    
    /** @brief return if indices start with one (true) or zero (false) */
    bool is_fortran_shifted() const;
    
    /** 
     * @brief return a new structure for a transposed matrix
     * 
     * If this structure has been created using finite element spaces, then 
     * this might be different from directly constructing a structure with 
     * test and anstz space exchanged. This can happen due to the non-active 
     * degrees of freedom.
     * 
     * This function returns a true (algebraic) transposed structure.
     */
    std::shared_ptr<TStructure> get_transposed() const;
    
    /** 
     * @brief print out some information on this object
     * 
     * This depends on the current verbosity level.
     */
    void info() const;
    
    /** 
     * @brief draw a postscript picture of the sparsity pattern, similar to 
     * Matlab 'spy'
     * 
     * @param filename a file with this name will be created (overwritten)s
     */
    void draw(const std::string& filename) const;
    
    /**
     * @brief return a structure for the matrix-matrix-product A*B
     * 
     * if A and B are matrices with structures 'strucA' and 'strucB', this
     * function computes a structure for the product C = A*B
     * 
     * @param strucA structure of left factor
     * @param strucB structure of right factor
     */
    friend std::shared_ptr<TStructure> get_product_structure(
        TStructure const & strucA, TStructure const & strucB);
    
    /**
     * @brief Compute the structure of the matrix matrix product A*A^T.
     * 
     * @return A pointer to the structure of A*A^T.
     * @note Relies on sorted columns array, i.e. ColOrder should be 1.
     */
    TStructure* get_structure_of_product_with_transpose_from_right() const;
    
    /**
     * @brief Compute the structure of the matrix matrix product A*B*A^T.
     * 
     * @param B structure of the middle matrix B
     * @return A pointer to the structure of A*B*A^T.
     * @note Relies on sorted columns array, i.e. ColOrder should be 1.
     */
    TStructure* get_structure_of_product_with_transpose_from_right(
      const TStructure & B) const;
    
    /** @brief Comparision Operator 
     * 
     * All the integers and all the arrays are compared. Two structures are
     * identical iff all the values and arrays coincide.
     */
    friend bool operator==(const TStructure &lhs, const TStructure &rhs);
    friend bool operator!=(const TStructure &lhs, const TStructure &rhs);

  private:
    /** @brief number of rows */
    unsigned int nRows;

    /** @brief number columns */
    unsigned int nColumns;

    /** @brief number of matrix entries */
    unsigned int nEntries;

    /** @brief vector of column information
     * 
     * This vector has length TStructure::nEntries and stores the column index
     * for each entry.
     */
    std::vector<int> columns;
    
    /** @brief vector storing the number of entries in all previous rows
     * 
     * This vector tells you the global indices of the entries for each row.
     * The `j`-th entry in the `i`-th row is at the `TStructure[i]+j`-th 
     * position in the global entries vector. So this also gives you the right 
     * index if you are interested in the column of that entry, 
     * `TStructure::columns[TStructure::rows[i]+j]`.
     */
    std::vector<int> rows;

    /** @brief sort all rows in increasing order */
    void Sort();

    /** @brief sort one row to increasing order */
    void SortRow(int *BeginPtr, int *AfterEndPtr);
};

#endif
