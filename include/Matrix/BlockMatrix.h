#ifndef USER_PROJECTS_BLOCKMATRIX_H_
#define USER_PROJECTS_BLOCKMATRIX_H_

#include <Matrix.h>
#include <CompositeOperator.h>
#include <memory>
#include <tuple>

#include <vector>

class BlockVector;

/** ************************************************************************
 *
 * @class      BlockMatrix
 * @brief      represent a matrix consisting of blocks which are sparse matrices
 *
 *             This is a purely algebraic object. And one which does frequently
 *             appear in FEM. Identifying and exploiting a certain block structure
 *             of matrices is a fundamental concept of solvers in CFD.
 *
 *             So far the BlockMatrix does not check, whether it
 *             already holds a block which you want to assign to one of its
 *             cells, but makes its own copy of whatever you give to it.
 *
 *             The BlockMatrix offers the possibility to store matrices,
 *             which appear multiple time in the block system, only once.
 *             It is also supported to store matrices in transposed state.
 *
 *             A BlockMatrix basically consists of its 'cells', which
 *             forms the lattice where matrices are stored, and its blocks.
 *             The basic concept to realize single storage of multiple blocks
 *             is the following: Each cell has along a shared pointer to its
 *             block a 'color', which is nothing but a non-negative integer number.
 *             Cells which hold the same block have the same color, cells
 *             which hold different blocks have different colors. All colors
 *             from zero to the current number of different stored blocks must
 *             be assigned.
 *
 *             The BlockMatrix supports changes in its color structure
 *             which are necessary due to changes which affect only some of its blocks,
 *             or affect different blocks differently (e.g. storing new blocks,
 *             adding to some blocks only, scaling some blocks only) and figures
 *             out independently, to which TMatrices it has to delegate blockwise
 *             tasks in order to not perform them multiple times falsely.
 *
 *             An auxiliary vector "color_count" is held, which keeps track
 *             of the number of cells of each color. This additional structure
 *             helps in quickly deciding whether a change leads to a change in
 *             coloring.
 *
 *             Alongside the color, each cell knows whether its block is held
 *             in transposed or non-transposed state.
 *
 *             There is one principle when working with a BlockMatrix:
 *             it is only as clever as you make it! In particular: it can only
 *             maintain a 'good' block structure (that is not to say: a valid one)
 *             if you perform your actions on it in a way fo the program to maintain
 *             them. If you do differently, it will print a warning (but
 *             only from verbosity level 2 on), to give
 *             you the possibility to correct your approach.
 *
 *
 *             TODO Implement methods which take an entire block matrix as argument.
 *             (Especially summation, probably multiplication)
 *             So far the input is always only ONE TMatrix.
 *             The implementation of such methods should not be a big trouble,
 *             because one can easily break them down to modifications with only one
 *             TMatrix affecting multiple cells in the BlockMatrix.
 *
 *
 * @author     Clemens Bartsch, Ulrich Wilbrandt
 * @date       2015/12
 *
 * @ruleof0
 *
 ****************************************************************************/
class BlockMatrix
{
  public:

    /**
     * Default constructor. Creates an empty 0x0 object.
     */
    BlockMatrix();

    // implicit cast - throws an error
    BlockMatrix(CompositeOperator<BlockVector>);

    /**
     * @brief Creates a BlockMatrix which is filled with fitting zero blocks.
     *
     * The number of block rows is the length of cell_row_numbers,
     * the number of block columns the length of cell_column_numbers.
     *
     * If e.g. cell_row_numbers[3] == 100, all cells in cell row 3 will have
     * 100 matrix rows.
     *
     * @param cell_row_numbers holds the number of rows for all matrices in one cell row
     * @param cell_column_numbers holds the number of columns for all matrices in one cell column
     */
    BlockMatrix(const std::vector<size_t>& cell_row_numbers,
                const std::vector<size_t>& cell_column_numbers);

    /**
     * @brief Creates a nRows times nCols BlockMatrix filled with blocks
     *
     * The blocks are given row wise. That means block (i,j) in the resulting
     * BlockMatrix will be the (i*nCols+j)-th block in \c blocks. In other
     * words this creates a BlockMatrix where each block has its own color and
     * is not stored as transposed.
     *
     * The caller has to make sure all blocks are appropriate, otherwise this
     * constructor will throw an exception.
     *
     * @param nRows - number of blocks per column
     * @param nCols - number of blocks per row
     * @param blocks - the given blocks, must be of length nRows*nCols
     */
    BlockMatrix(int nRows, int nCols,
                const std::vector<std::shared_ptr<TMatrix>>& blocks);

    /**
     * Add a given TMatrix to the blocks in a bunch of given cells at once.
     * Of course the TMatrix must fit in dimension to the cells where it is
     * supposed to be added.
     * The input will be checked for its validity. Duplicates are removed,
     * is the input corrupted in other ways the program will quit, informing
     * you what exactly it does not like.
     *
     * @param[summand] The TMatrix to be added to certain blocks.
     * @param[in] cell_positions The places in the cell grid where to add.
     * @param[in] transposed_states If the TMatrix goes to a specific place
     * in transposed ('true') or non-transposed ('false') state.
     */
    void add_unscaled_matrix(
        const TMatrix& summand,
        const std::vector<std::vector<size_t>>& cell_positions,
        const std::vector<bool>& transposed_states);

    /**
     * Performs as add_matrix_to_blocks but with the possibility
     * to scale the summand by
     * @param[in] scaling_factor A scaling factor for the summand.
     */
    void add_matrix(
        const TMatrix& summand, double scaling_factor,
        const std::vector<std::vector<size_t>>& cell_positions,
        const std::vector<bool>& transposed_states);

    /** @brief compute y = Ax
     *
     * Write the matrix-vector product "Ax" to y if "A" is this matrix.
     * This method is purely algebraical. It does the multiplication
     * with exactly those values which it finds in the blocks of "this",
     * regardless of any non-active dofs - that's what BlockFEMatrix is for.
     *
     * @param[in] x the BlockVector which is multiplied by this matrix
     * @param[out] y result of matrix-vector-multiplication
     */
    virtual void apply(const BlockVector & x, BlockVector & y) const;

    /** @brief compute y = y + a * Ax
     *
     * Add the matrix-vector product "Ax", scaled by "a", to y if "A" is this
     * matrix.
     *
     * This function can be used to compute the residual r = b - Ax, for example
     * BlockVector r(b);  // structural copy (no values copied yet)
     * r = b;             // copy values from b to r
     * A.apply_scaled_add(x, r, -1.0);
     * cout << "Norm of residual " << r.norm() << endl;
     *
     *
     * @param x the BlockVector which is multiplied by this matrix
     * @param y result of matrix-vector-multiplication
     * @param a optional factor, defaults to 1.0
     */
    virtual void apply_scaled_add(const BlockVector & x, BlockVector & y,
                          double a = 1.0) const;

    virtual void apply_transpose(const BlockVector & x, BlockVector & y) const;
    virtual void apply_transpose_scaled_add(const BlockVector & x, BlockVector & y,
                          double a = 1.0) const;

    /// @brief perform one successive overrelaxation (sor) sweep.
    /// The flag can be either 0(forward sweep), 1(backward sweep), or
    /// 2(forward followed by backward sweep).
    ///
    /// @param[in] b right hand side
    /// @param[in,out] x solution (this is updated)
    /// @param[in] omega relaxation parameter
    /// @param[in] flag either 0 (forward), 1(backward), or 2(both)
    /// @param[in] par_strat The chosen parallelization strategy (MPI only).
    ///            Choose between "all_cells", "halo_0" and "own_cells".
    ///            Note that this does only affect the algorithm in TMatrix::sor_sweep.
    ///            There won't be any difference in here - although one could
    ///            change TMatrix' multiply and transpose_multiply to skip all
    ///            Halo rows in order to save some flops.
#ifdef _MPI
    void sor_sweep(const BlockVector& b, BlockVector& x, double omega,
                   size_t flag, const std::string& par_strat) const;
#else
    void sor_sweep(const BlockVector& b, BlockVector& x, double omega,
                   size_t flag) const;
#endif

    /**
     * @brief checks whether the coloring is correct - use in tests only
     *
     * The method checks, whether the coloring of the cells is correct, e.g.
     *
     *  - the set of assigned colors equals the set {0,...,n_colors_ -1}
     *
     *  - on the set of cells the equivalence relation "has the same color as"
     *    is equivalent to "holds a pointer to the same matrix as"
     *
     *  - the colors are ordered in such a way that from left to right, top to bottom
     *    the colors of first appearing matrices are in ascending order.
     *
     *    Does nothing but throw if one of the rules stated above is broken.
     *    Is written with comprehensibility in mind, not performance. Use
     *    for testing purpose only.
     *
     *    @note Please note that this method scales quadratically in the number of blocks
     *    as it employs checkEquivalenceOfRelations. I used it for testing
     *    of other methods in this class.
     */
    void check_coloring() const;

    ///! Check whether a BlockVector b is fit to be the rhs b of the equation Ax=b.
    virtual void check_vector_fits_image(const BlockVector& b) const;

    ///! Check whether a BlockVector x is fit to be to be factor x in the equation Ax=b.
    virtual void check_vector_fits_pre_image(const BlockVector& x) const;

    /** @brief return this BlockMatrix as one TMatrix
     *
     * This returns a merged version of this matrix. Note that the merged
     * matrix does not get stored internally, for it cannot easily be kept
     * up to date, but is recreated on every call.
     *
     * Usually this is used to pass this matrix to a solver.
     *
     * @return A shared pointer to the block matrix, merged together to
     * a TMatrix.
     *
     * @todo This method is performance critical when using a direct solver, but
     * it is rather slow - profile and speed up!
     * Make sure that zero entries are not put into the combined matrix - this
     * might already speed things up.
     */
    virtual std::shared_ptr<TMatrix> get_combined_matrix() const;


    /// Combines a rectangular submatrix specified by its upper-leftmost
    /// and lower-rightmost block into a TMatrix.
    /// @param upper_left
    /// @param lower_right
    /// @return into a TMatrix.
    virtual std::shared_ptr<TMatrix> get_combined_submatrix(
        std::pair<size_t,size_t> upper_left,
        std::pair<size_t,size_t> lower_right) const;

    // Getter.

    //! @return The number of different colors. Should be needed for testing only.
    size_t get_n_colors() const
    {
      return color_count_.size();
    }

    //! @return The number of cell rows.
    size_t get_n_cell_rows() const
    {
      return n_cell_rows_;
    }

    //! @return the number of cell columns.
    size_t get_n_cell_columns() const
    {
      return n_cell_columns_;
    }

    /// @return The number of columns a certain cell has (same for all cells in one cell columns)
    size_t get_n_columns_in_cell(size_t cell_row, size_t cell_column) const
    {
      std::vector<grid_place_and_mode>index {std::make_tuple(cell_row, cell_column, true)};
      check_indices(index);
      return cell_grid_[cell_row][cell_column].n_columns_;
    }

    /// @return The number of rows a certain cell has (same for all cells in one cell row)
    size_t get_n_rows_in_cell(size_t cell_row, size_t cell_column) const
    {
      std::vector<grid_place_and_mode>index {std::make_tuple(cell_row, cell_column, true)};
      check_indices(index);
      return cell_grid_[cell_row][cell_column].n_rows_;
    }


    /// @brief total number of columns (added over all cells in a cell row)
    size_t get_n_total_columns() const;

    /// @brief total number of entries (added over all blocks)
    size_t get_n_total_entries() const;

    /// @brief total number of rows (added over all cells in a cell column)
    size_t get_n_total_rows() const;

    /// @brief return the number of square blocks in this block matrix
    size_t get_n_square_blocks() const;

    /**
     * Spawn a new BlockMatrix which is a rectangular submatrix of this matrix.
     *
     * @param upper_left The upper leftmost block to include.
     * @param lower_right The lower rightmost block to include.
     * @return
     */
    BlockMatrix get_sub_blockmatrix(
        const std::pair<size_t,size_t>& upper_left,
        const std::pair<size_t,size_t>& lower_right) const;

    /**
     * Prints matrix coloring pattern and color_count_,
     * then runs the check. Thus one can visualize errors.
     *
     * @param[in] matrix_name A name for the matrix. Will be printed in the heading.
     */
    void print_and_check(const std::string& matrix_name = "unnamed") const;

    /**
     * Print out the current color count of the matrix.
     * Does not perform a check on consistency - use it for debugging!
     *
     * @param[in] matrix_name A name for the matrix. Will be printed in the heading.
     */
    void print_coloring_count(const std::string& matrix_name = "unnamed") const;

    /**
     * Print out a little picture of the current coloring pattern of the matrix.
     * Does not perform a check on consistency - use it for debugging!
     * A matrix stored in transposed is marked with "^T", a matrix
     * without entries (as happens in our NSE-Types) with brackets.
     *
     * @param[in] matrix_name A name for the matrix. Will be printed in the heading.
     */
    void print_coloring_pattern(const std::string& matrix_name = "unnamed",
                                bool print_adress = false) const;

    /**
     * Replace the blocks at some cells with a new one. Updates the
     * coloring accordingly.
     * If the input does not make sense it will throw an error
     * and print a sensible message.
     *
     * @param[in] new_block The new block. Will be deep copied.
     * @param[in] cell_positions In which cells to put the new block.
     * @param[in] transposed_states In which transposed state to put it there
     * respectively.
     */
    virtual void replace_blocks(
        const TMatrix& new_block,
        const std::vector<std::vector<size_t>>& cell_positions,
        const std::vector<bool>& transposed_states);

    /**
     * @brief scale the entire matrix
     *
     * That means for each block all entries are scaled.
     *
     * @param[in] factor The scaling factor.
     */
    void scale(double factor);

    /**
     * Scales some blocks, if necessary updates the coloring scheme.
     *
     * @param[in] scaling_factor The scaling factor.
     * @param[in] cell_positions The cells whose blocks are to be scaled.
     */
    void scale_blocks(
        double scaling_factor,
        const std::vector<std::vector<size_t>>& cell_positions );

    /// @brief read an individual entry
    ///
    /// @note this will not fail if the desired entry is not in the sparsity
    /// structure. In that case it will simply return 0.
    double get(unsigned int i, unsigned int j) const;

    /// @brief return the diagonal entries of this BlockMatrix
    std::vector<double> get_diagonal() const;

    std::vector<double> get_row_sums() const;
    std::vector<double> get_col_sums() const;

    /**
     * @brief Get a shared pointer to a constant version of one of the blocks.
     * Note that the non-active rows of that block might not be what you expect
     * and will have to be read with care.
     *
     * @param[in] cell_row The cell row of the desired block.
     * @param[in] cell_col The cell column of the desired block.
     * @param[out] is_transposed A flag which shows true if the block is stored
     * in transposed state, false if not so.
     *
     * @return A shared pointer to a block.
     */
    std::shared_ptr<const TMatrix> get_block(
        size_t cell_row, size_t cell_col, bool& is_transposed) const;

    // Special member functions.

    /** @brief copy constructor
     * Performs a deep copy of the stored TMatrices,
     * so that no two instances of BlockMatrix share the same blocks.
     */
    BlockMatrix(const BlockMatrix&);

    ///! Default move constructor does the job.
    BlockMatrix(BlockMatrix&&) = default;

    /** Swap function used for copy-and swap in copy assignment.
     * @param[in,out] first The object to be swapped with second.
     * @param[in,out] second The object to be swapped with first.
     */
    friend void swap(BlockMatrix& first, BlockMatrix& second);

    /** @brief Unified assignment operator
     * Performs a deep copy of the stored TMatrices, using copy-and-swap,
     * so that no two instances of BlockMatrix share the same blocks.
     */
    BlockMatrix& operator=(BlockMatrix);


    /// @brief Default destructor. Tidies up nice and clean.
    virtual ~BlockMatrix() = default;

    /// @brief Set all submatrices to zero
    void reset();

    int clean_denormals();

  protected:

    //! Store information of a certain grid cell in the block matrix.
    //! Will by default perform a shallow copy when copied, which is what we want.
    struct CellInfo
    {
      /// the number of rows a matrix stored here must possess
      size_t n_rows_;

      /// the number of columns a matrix stored here must possess
      size_t n_columns_;

      /// The TMatrix currently associated with this cell.
      std::shared_ptr<TMatrix> block_;

      /**
       * The color of the cell - std::numeric_limits<size_t>::max()
       * stands for 'uncolored'.
       */
      size_t color_;

      //! whether or not the stored block is viewed as transposed
      bool is_transposed_;

      /** A flag used in algorithms which change the coloring scheme
       *  of the owning BlockMatrix.
       */
      enum class ReColoringFlag {SPLIT, KEEP} re_color_flag_;

      /*! Default constructor, will be called implicitely. Constructs
       *  a non-transposed, uncolored 0x0 object holding a null pointer
       *  and marked with ReColoringFlag:KEEP.
       */
      CellInfo();

      /*! Constructor which initializes n_rows and n_columns,
       *  Behaves like default constructed object apart from that.
       */
      CellInfo(size_t n_rows, size_t n_columns);

    };

    //! The number of cell rows. Equals the number of blocks/cells in each cell column.
    size_t n_cell_rows_;

    //! The number of cell columns. Equals the number of blocks/cells in each cell row.
    size_t n_cell_columns_;

    /*!
     *  A tableau of cell information. Has dimension
     *   n_block_rows x n_block_columns, and stores cell information.
     *   I.e. cell_grid_[i][j] holds information on the cell in block row i,
     *   block column j.
     */
    std::vector<std::vector<CellInfo>> cell_grid_;

    /*!
     * The size of the vector is the number of colors currently
     * existing in the block matrix, i.e. the number of physically
     * stored TMatrices.
     *
     * Storing these numbers enables quick checking whether a
     * modification of the blocks will require a modification
     * of the coloring scheme.
     *
     * color_count_[i] gives you the current number of cells
     * with color i in the block matrix
     */
    std::vector<size_t> color_count_;

    /** A datatype which stores the block_row the block_column
     * and the mode where and how a cell should be modified.
     * 'Mode' can have to values: true - 'transposed' and
     * false - 'not transposed'.
     * Used internally to keep information tied together.
     */
    typedef std::tuple<size_t, size_t, bool> grid_place_and_mode;

  protected:

    /**
     * Add one matrix to several blocks at once.
     *
     * @param[in]  summand the matrix to be added.
     * @param[out] row_column_transpose_tuples A vector of tuples,
     * each of which encodes information on
     * one cell where summand is supposed to be added to,
     * and which way - "true" is for transposed,
     * "false" for not-transposed.
     *
     */
    void add_matrix_to_blocks(
        const TMatrix& summand,
        std::vector<grid_place_and_mode> row_column_transpose_tuples);

    /**
     * Add a scaled coyp of one matrix to several blocks at once.
     *
     * @param[in]  summand the matrix whose scaled copy is to be added.
     * @param[in] scaling_factor the factor by which to scale summand before adding
     * @param[out] row_column_transpose_tuples A vector of tuples,
     * each of which encodes information on one cell where summand
     * is supposed to be added to, and which way - "true" is for transposed,
     * "false" for not-transposed.
     */
    void add_scaled_matrix_to_blocks(
        const TMatrix& summand, double scaling_factor,
        std::vector<grid_place_and_mode> row_column_transpose_tuples);

    /**
     * Checks an input vector of places and edits it:
     *   - sorts in alphanumeric order of grid places
     *   - removes duplicates
     *   - throws if the input contains duplicates of the type
     *     (a,b,true) , (a, b, false) - mode at place has to
     *     be unique!
     *
     * @param[in,out] The input vector of grid places and modes.
     *
     */
    static void check_and_edit_input(
        std::vector<grid_place_and_mode>& row_column_transpose_tuples);


    /**
     * Checks and transforms input vectors to vectors of tuples
     * which the working methods expect.
     *
     * cell_positions and transposed_states must have the same lengths,
     * and all entries of cell_positions must be of size 2 themselves.
     */
    static std::vector<grid_place_and_mode> check_and_tupelize_vector_input(
        const std::vector<std::vector<size_t>>& cell_positions,
        const std::vector<bool>& transposed_states );

    /**
     * Check whether the set of assigned colors equals the set
     * {0,...,n_colors_ -1}. If not so - just throw.
     */
    void check_color_count() const;

    /**
     * This method checks whether the coloring fulfils the internal ordering.
     * If not so it just throws.
     *
     * The method relies on conditions 1 and 2 being established.
     */
    void check_coloring_order() const;

    /**
     * Check whether on the set of cells the equivalence relation
     * "has the same color as" is equivalent to "holds a pointer
     * to the same matrix as". If not so - just throw.
     *
     * The result is independent of the check which
     * checkColorAssignment performs.
     *
     * @note That this method scales quadratically in the number of blocks,
     * so by n^4 for a nxn block matrix - so please do only use it if you
     * a) want to see your computer failing or
     * b) test other methods in this class.
     */
    void check_equivalence_of_relations() const;

    /**
     * Checks whether a matrix fits to the cells and mode specified in
     * row_column_transpose_tuples.
     * Throws if not so.
     *
     * @param[in] matrix The matrix whose fit into the grid shall be checked.
     *
     * @param[in] row_column_transpose_tuples The grid places and modes to check.
     * The method check_and_edit_input is called upon it to make sure it is
     * valid input.
     */
    void check_grid_fit(
        const TMatrix& matrix,
        std::vector<grid_place_and_mode>& row_column_transpose_tuples) const;

    /**
     * Check if there is an index out of bound issue in the given grid_place_and_mode
     * vector and throw if so.
     */
    void check_indices(
             const std::vector<grid_place_and_mode>& row_column_transpose_tuples
         ) const;

    /**
     * This method checks whether all cells have their re-coloring flags set to KEEP.
     * If not so it just throws.
     */
    void check_re_coloring_flags() const;

    /**
     * Check if input transposed state fits to the transposed state of the
     * cell where the input should go. If not so: throw.
     *
     * TODO This is not perfect yet - actually the block matrix should be
     * able to differ three cases per affected color:
     *  - 1) all input transp states match the cell transp states - just add!
     *  - 2) all input transp states do not match the cell transp states - transp and add!
     *  - 3) some input tranps states match, some don't - throw! for this will only work out in
     *    symmetric structure cases on the TMatrix anyway and seems just not worth the effort...
     *    and if there's different occurences of symmetric TMatrices, why store them
     *    in different transposed states from the beginning?
     * So far this only handles only case 1 and 3.
     */
    void compare_transposed_mode(
        std::vector<grid_place_and_mode>& row_column_transpose_tuples) const;

    /** Copies a given TMatrix and wraps a smart pointer around it, which is returned.
     *  Is overridden in derrived class BlockFEMatrix, where a pointer cast has to
     *  be performed due to the storing of FEMatrices.
     */
    virtual std::shared_ptr<TMatrix> create_block_shared_pointer(const TMatrix& block) const;

    /*!
     * Check if a given block fits into a given cell in the given transposed state.
     *
     */
    static bool does_block_fit_cell(
        const TMatrix& block, const CellInfo& cell, bool transposed );

    /*!
     * Check if for two cells the relations "has the same color as"
     * and "store a pointer to the same matrix" produce the same result.
     *
     * @param[in] first  the first CellInfo object to compare
     * @param[in] second the second CellInfo object to compare
     */
    static bool does_color_match_block(const CellInfo& first, const CellInfo& second);

    /**
     * Checks if a modification which affects a set of cells requires at least
     * one color class to be split in two.
     *
     * @return True if the modification requires to split at least one color in two.
     *
     * @param[out] color_to_split The lowest color which has to be split in two.
     * @param[in] row_column_tuples The index pairs of the cells affected by
     * the modification. Must not contain the same index pair twice! Need not be sorted.
     *
     * When returning true the parameters of the method should be handed over to
     * mark_for_color_split and then to split_color
     */
    bool does_modification_require_color_split(
      size_t& color_to_split,
      const std::vector<grid_place_and_mode> row_column_transposed_tuples) const;

    /**
     * @return True if a block emplacement leads to the merging of at least two colors.
     *
     * @param[out] color_a one color which enters the merge
     * @param[out] color_b other color which enters the merge
     *
     * If no merge is needed, color_a and color_b both have the value of the color class,
     * to which all matrices in row_column_transposed_tuples belong.
     *
     * @note: To maintain the coloring conditions the method may only be called,
     * when it is ensured, that a replacement does not require color splits of any kind!
     * So call it only after does_modification_require_color_split returned false!
     */
    bool does_replace_require_color_merge (
        size_t& color_a, size_t& color_b,
        std::vector<grid_place_and_mode> row_column_transposed_tuples) const;

    /*!
     * Find the first place in the matrix (ordered from left to right, top to bottom)
     * where a certain color appears.
     *
     * @param[in] color_to_find the color to find
     *
     * @param[out] block_row the row where the first such block was found
     * @param[out] block_column the column where the first such block was found
     */
    void find_first_appearance_of_color(size_t color_to_find,
                                        size_t& block_row , size_t& block_column
                                        ) const;

    /*!
     * Find the first place in the matrix (ordered from left to right, top to bottom)
     * where a certain color appears in a certain mode (transposed or non-transposed).
     * The method will fire an exception if no block of the color in  such state was found.
     *
     * @param[in] color_to_find the color to find
     * @param[in] transposed true if you want to find the first transposed appearance,
     *                       false for the first non-transposed appearance
     *
     * @param[out] block_row the row where the first such block was found
     * @param[out] block_column the column where the first such block was found
     */
    void find_first_appearance_of_color_and_mode(size_t color_to_find, bool transposed,
                                        size_t& block_row , size_t& block_column
                                        ) const;

    /**
     * Get the next place in the grid after block_row, block_column.
     * Throws an std::logic_error exception if the last index pair
     * is given as input.
     *
     * @param[in, out] block_row the current block row
     * @param[in, out] block_column the current block column
     *
     * @throws std::logic_error exception if the last index pair
     * is given as input
     */
    void get_next_cell_grid_index(size_t& block_row, size_t& block_column) const;

    /** @brief this method is used to compare the number of actives in a block vector
     * to the number of actives in test space
     *  @param nActive number of actives
     *  @param spaceNumber number of the test space to compare the actives
     */
    virtual void handle_discovery_of_vector_non_actives(
      const int nActive, const int spaceNumber) const;
    /**
     * Check if a given index pair is the last one in the cell_grid_.
     *
     * @param[in] block_row the current block row
     * @param[in] block_column the current block column
     *
     * @return true if this is the index pair of the lower right corner
     */
    bool is_last_index_pair(size_t block_row, size_t block_column) const;

    /**
     * Mark every cell of color_to_mark whose position appears in row_column_transposed_tuples
     * with the recoloring flag "SPLIT". After that, call split_color to perform the
     * actual splitting of the color.
     *
     * @param[in] color_to_mark Cells of color color_to_mark are marked when their
     * position belongs to row_column_tuples.
     * @param[in] row_column_transposed_tuples Cells of color color_to_mark are marked when their
     * position belongs to row_column_tuples.
     */
    void mark_for_color_split(
      size_t color_to_mark,
      const std::vector<grid_place_and_mode>& row_column_transposed_tuples);

    /**
     * Merges two colors marked with ReColoringFlag::MERGE into one.
     * Assumes conditions 1, 2 and 3 to hold and maintains them.
     *
     */
    void merge_colors(size_t color_a, size_t color_b);

    /**
     * Replaces the blocks whose positions are given in grid_places by a copy of
     * new_block, as long as that fits into the given grid postions
     *
     * Expects conditions 1, 2 and 3 to hold and maintains them.
     *
     * @param[in] new_block The new block to be inserted.
     *
     * @param[in] row_column_transpose_tuples The places where the new block should to be inserted.
     */
    void replace_blocks(
        const TMatrix& new_block,
        std::vector<grid_place_and_mode> row_column_transpose_tuples);

    /**
     * Scale several blocks in the matrix at once by a factor.
     *
     * @param[in]  scaling_factor the factor by which to scale the blocks
     * @param[out] A vector of tuples, each of which encodes information on
     *             one block to be scaled. The third entry of each tuple
     *             (transposed or not) is not relevant here (except for the
     *             input consistency check!).
     */
    void scale_blocks( double scaling_factor,
                       std::vector<grid_place_and_mode> row_column_transpose_tuples);

    /**
     * Splits a color in two, assuming that conditions 1,2 and 3 hold and maintaining them.
     * Also relies on one part of the cells of the color to be split are marked with
     * ReColoringFlag::SPLIT. (Make sure to call mark_for_color_split in advance).
     *
     * @param[in] color_to_split The color which has to be split in two
     */
    void split_color(size_t color_to_split);


    // not yet adapted members of BlockMatrix - serves as a TODO list

  public:

    /**
     * @brief adding a scaled matrix to this matrix
     *
     * The summation is index-wise, i.e. M(i,j) += factor*A(i.j), where M is
     * this matrix.
     *
     * Note that this only works if the sparsity structure is the same for this
     * matrix and A.
     *
     * Possibly existing special matrices are not changed.
     */
    void add_scaled(const BlockMatrix &A, double factor = 1.0);

};



#endif /* USER_PROJECTS_BLOCKMATRIX_H_ */
