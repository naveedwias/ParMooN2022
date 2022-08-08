/** ************************************************************************
 *
 * @class      BlockFEMatrix
 * @brief      extends BlockMatrix by handling active degrees of freedom
 *
 *             A BlockMatrix of subtype BlockFEMatrix stores FEMatrices
 *             instead of simple algebraic TMatrices. Thus it can access
 *             information on finite element spaces and active degrees of
 *             freedoms and exploits this algorithmically. The majority of
 *             BlockMatrices in ParMooN will in fact be BlockFEMatrices.
 *
 *             Each cell row is associated with a certain FE test space
 *             and each cell column with a certain FE ansatz space.
 *
 *             The possibility to store matrices as transposed
 *             leads to a confusion of Test- and Ansatzspace between the
 *             BlockFEMatrix and the stored FEMatrices. If an FEMatrix is
 *             stored as transposed somewhere and one requests the testspace
 *             of the matrix stored there, it gives a different space than
 *             when requiring the testspace of the cell where it is stored.
 *             Because only the BlockFEMatrix knows, whether a block is stored
 *             as tranposed or not, its "answer" is the correct one.
 *             For the moment the only advice we can give is: do not request
 *             test- or ansatzspace of an FEMatrix stored as a block in a
 *             BlockFEMatrix, but instead request the test- or ansatzspace
 *             of the cell where it is stored.
 *
 *             The current implementation of BlockFE Matrix has restrictions
 *             due to the way ParMooN handles non-active degrees of freedom.
 *             It is restricted to:
 *               - a symmetrical test- and ansatzraum structure - i.e. the
 *                 rowwise sequence of testspaces is THE SAME as the columnwise
 *                 sequence of ansatzspaces
 *               - whose non-active dofs stem only from Dirichlet boundaries,
 *                 i.e. no hanging nodes! (which is the second reason for non-
 *                 active degrees of freedom in ParMooN)
 *             The reason for these restriction is in the specific handling of,
 *             non-active dofs, i.e. dofs from Dirichlet boundary conditions and
 *             hanging nodes. The issue is, that the corresponding rows are set in
 *             each FEMatrix seperately, depending on its testspace. The FEMatrices
 *             do not know about the space dimension they belong to, which leads
 *             to a discrepancy between the BlockFEMatrix' global non-active row
 *             and the FEMatrices' local non-active rows. We do not know how to entangle
 *             that relation at the moment, only in the case given above it is quite
 *             easily done: of all blocks in one row only the non-active entries of the
 *             block on the diagonal are of relevance. This circumstance is paid
 *             respect in the methods
 *               - apply (only indirectly)
 *               - apply_scaled_add
 *               - check_vector_fits_pre_image (actives check)
 *               - check_vector_fits_image (actives check)
 *               - get_combined_matrix.
 *             Should you intend to extend the class to the handling of hanging
 *             nodes and/or non-symmetric test- and ansatzspaces, these are the
 *             methods to changes (and of course implementing a new constructor
 *             and/or remove the "no-hanging-nodes" invariant.)
 *
 * In the case of (Navier-) Stokes where the velocity space has only Dirichlet
 * boundaries, the matrix is typically singular with the constant pressure
 * function being in the kernel of the matrix. Usually one restricts the
 * pressure space by one dimension, which is done here by setting the first
 * pressure row to be zero, except on the diagonal (where one has an entry 1.0).
 * This 'pressure correction' is mainly needed for direct solvers. So it is
 * implemented in the methods get_combined_matrix. The bool
 * `pressure_correction` determines if this should be done or not. We are aware
 * of the fact that this solution is not very nice.
 *
 * @author     Clemens Bartsch
 * @date       2015/12/08
 *
 * @ruleof0
 *
 ****************************************************************************/

#ifndef USER_PROJECTS_BLOCKFEMATRIX_H_
#define USER_PROJECTS_BLOCKFEMATRIX_H_

#include <BlockMatrix.h>
#include <FEMatrix.h>
#include <FESpace.h>
#include <templateNames.h>
#ifdef __3D__
using FESpace = TFESpace3D;
#else
using FESpace = TFESpace2D;
#endif

#ifdef _MPI
class TParFECommunicator3D;
#endif

class BlockFEMatrix : public BlockMatrix
{
  public:
    //constructors
    /**
     * @brief Creates a BlockFEMatrix which is filled with fitting zero blocks.
     * The TStructure of each of these zero blocks is the empty zero-map TStructure.
     *
     * The number of block rows and the number of block columns
     * is the length of spaces.
     *
     * @param[in] spaces FE spaces, which represent the test spaces rowwise
     * as well as the ansatz spaces columnwise. Each space
     * applies to all blocks of a particular row as testspace and
     * to all blocks of a particular column as ansatz space.
     */
#ifdef __2D__
    explicit BlockFEMatrix(const std::vector<std::shared_ptr<const TFESpace2D>>& spaces);
#elif __3D__
    explicit BlockFEMatrix(const std::vector<std::shared_ptr<const TFESpace3D>>& spaces);
#endif // 3D

    /**
     * Default constructor. Constructs emtpy matrix which maps zero space to zero space.
     */
    BlockFEMatrix();

    /**
     * @brief Creates a nRows times nCols BlockFEMatrix filled with blocks
     *
     * The blocks are given row wise. That means block (i,j) in the resulting
     * BlockFEMatrix will be the (i*nCols+j)-th block in \c blocks. In other
     * words this creates a BlockFEMatrix where each block has its own color and
     * is not stored as transposed.
     *
     * The caller has to make sure all blocks are appropriate, otherwise this
     * constructor will throw an exception.
     *
     * @param nRows - number of blocks per column
     * @param nCols - number of blocks per row
     * @param blocks - the given blocks, must be of length nRows*nCols
     */
    BlockFEMatrix(int nRows, int nCols,
                  const std::vector<std::shared_ptr<FEMatrix>>& blocks);

    // named constructors for block fe matrices often used in ParMooN
    //TODO All named constructors should further reduce the number of TMatrix-Copies made
    // in their body!

#ifdef __2D__
    /**
     * @brief Named constructor for a block Matrix used in 2D convection-
     * diffusion problems.
     *
     * Constructs a 1x1 block matrix which holds one single block with
     * fitting TStructure.
     *
     * How to use a named constructor? Have a look at the test file!
     *
     * @param space Ansatz- equals testspace.
     * @return A newly constructed BlockFEMatrix for CD2D problems.
     */
    static BlockFEMatrix CD2D(std::shared_ptr<const TFESpace2D> space);

    /**
     * Named constructor for a matrix of ParMooN-specific NSE Type 1.
     * The matrix takes the block structure
     *
     * ( A  0  B1T )
     * ( 0  A  B2T )
     * ( B1 B2 0   )
     *
     * where B1T and B2T are not explicitly stored (and marked non-transposed).
     *
     * How to use a named constructor? Have a look at the test file!
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE2D problems,
     * whose block structure is of NSE Type 1.
     */
    static BlockFEMatrix NSE2D_Type1(std::shared_ptr<const TFESpace2D> velocity,
                                     std::shared_ptr<const TFESpace2D> pressure);

    /**
     * Named constructor for a matrix of ParMooN-specific NSE Type 2.
     * The matrix takes the block structure
     *
     * ( A  0  B1T )
     * ( 0  A  B2T )
     * ( B1 B2 0   )
     *
     * where B1T and B2T are explicitly stored (and marked non-transposed).
     *
     * How to use a named constructor? Have a look at the test file!
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE2D problems,
     * whose block structure is of NSE Type 2.
     */
    static BlockFEMatrix NSE2D_Type2(std::shared_ptr<const TFESpace2D> velocity,
                                     std::shared_ptr<const TFESpace2D> pressure);

    /**
     * Named constructor for a matrix of ParMooN-specific NSE Type 3.
     * The matrix takes the block structure
     *
     * ( A11 A12 B1^T )
     * ( A21 A22 B2^T )
     * ( B1  B2  0    )
     *
     * where B1^T and B2^T are are not explicitly stored.
     *
     * How to use a named constructor? Have a look at the test file!
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE2D problems,
     * whose block structure is of NSE Type 3.
     */
    static BlockFEMatrix NSE2D_Type3(std::shared_ptr<const TFESpace2D> velocity,
                                     std::shared_ptr<const TFESpace2D> pressure);

    /**
     * Named constructor for a matrix of ParMooN-specific NSE Type 4.
     * The matrix takes the block structure
     *
     * ( A11 A12 B1T )
     * ( A21 A22 B2T )
     * ( B1  B2  0   )
     *
     * where B1^T and B2^T are explicitly stored (and marked non-transposed)..
     *
     * How to use a named constructor? Have a look at the test file!
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE2D problems,
     * whose block structure is of NSE Type 4.
     */
    static BlockFEMatrix NSE2D_Type4(std::shared_ptr<const TFESpace2D>,
                                     std::shared_ptr<const TFESpace2D> pressure);

    /**
     * Named constructor for a matrix of ParMooN-specific NSE Type 14.
     * The matrix takes the block structure
     *
     * ( A11 A12 B1^T )
     * ( A21 A22 B2^T )
     * ( B1  B2  C    ),
     *
     * where B1^T and B2^T are explicitly stored (and marked non-transposed).
     *
     * How to use a named constructor? Have a look at the test file!
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE2D problems,
     * whose block structure is of NSE Type 14.
     */
    static BlockFEMatrix NSE2D_Type14(std::shared_ptr<const TFESpace2D>,
                                      std::shared_ptr<const TFESpace2D> pressure);

    /**
     * Named constructor for a matrix for Darcy type problems in 2D.
     * Creates a 2x2 block matrix of the form
     *  ( A  B1' )
     *  ( B2 C   ),
     * where A is velocity-velocity coupling, B1' velocity-pressure,
     * B2 pressure-velocity and C pressure-pressure coupling.
     *
     * How to use a named constructor? Have a look at the test file!
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for Darcy problems in 2D.
     */
    static BlockFEMatrix Darcy2D(std::shared_ptr<const TFESpace2D> velocity,
                                 std::shared_ptr<const TFESpace2D> pressure);

    /**
     * Named constructor for a Mass matrix of ParMooN-specific NSE type 1 & 2
     *
     * @param velocity The velocity finite element space
     * @return A newly constructed BlockFEMatrix for 2D NSE problems
     */
    static BlockFEMatrix Mass_NSE2D(std::shared_ptr<const TFESpace2D> velocity);

    static BlockFEMatrix Mass_Matrix_NSE2D(
      std::shared_ptr<const TFESpace2D> velocity,
      std::shared_ptr<const TFESpace2D> pressure);

    /**
     * Named constructor for a mass matrix type 1, 2
     * The matrix takes the block structure
     *
     * ( M11  0  0 )
     * ( 0  M22  0 )
     * ( 0   0   0   )
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE2D problems,
     * whose block structure is of NSE Type 1, 2.
     */
    static BlockFEMatrix Mass_NSE2D_Type1(
      std::shared_ptr<const TFESpace2D> velocity,
      std::shared_ptr<const TFESpace2D> pressure);
    static BlockFEMatrix Mass_NSE2D_Type2(
      std::shared_ptr<const TFESpace2D> velocity,
      std::shared_ptr<const TFESpace2D> pressure);
    /**
     * Named constructor for a mass matrix type 3, 4
     * The matrix takes the block structure
     *
     * ( M11  M12  0 )
     * ( M21  M22  0 )
     * ( 0   0     0 )
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE2D problems,
     * whose block structure is of NSE Type 3, 4.
     */
    static BlockFEMatrix Mass_NSE2D_Type3(
      std::shared_ptr<const TFESpace2D> velocity,
      std::shared_ptr<const TFESpace2D> pressure);
    static BlockFEMatrix Mass_NSE2D_Type4(
      std::shared_ptr<const TFESpace2D> velocity,
      std::shared_ptr<const TFESpace2D> pressure);

    /**
     * Named constructor for a mass matrix type 14
     * The matrix takes the block structure
     *
     * ( M11  M12  MQ1 )
     * ( M21  M22  MQ2 )
     * ( QM1  QM2    0 )
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE2D problems,
     * whose block structure is of NSE Type 14.
     */
    static BlockFEMatrix Mass_NSE2D_Type14(
      std::shared_ptr<const TFESpace2D> velocity,
      std::shared_ptr<const TFESpace2D> pressure);

#elif __3D__
    /**
     * @brief Named constructor for a block Matrix used in 3D convection-
     * diffusion problems.
     *
     * Constructs a 1x1 block matrix which holds one single block with
     * fitting TStructure.
     *
     * How to use a named constructor? Have a look at the test file!
     *
     * @param space Ansatz- equals testspace.
     * @return A newly constructed BlockFEMatrix for CD3D problems.
     */
    static BlockFEMatrix CD3D(std::shared_ptr<const TFESpace3D> space );

    /**
     * Named constructor for a matrix for Darcy type problems in 3D.
     * Creates a 2x2 block matrix of the form
     *  ( A  B1' )
     *  ( B2 C   ),
     * where A is velocity-velocity coupling, B1' velocity-pressure,
     * B2 pressure-velocity and C pressure-pressure coupling.
     *
     * How to use a named constructor? Have a look at the test file!
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for Darcy problems in 2D.
     */
    static BlockFEMatrix Darcy3D(std::shared_ptr<const TFESpace3D> velocity,
                                 std::shared_ptr<const TFESpace3D> pressure);

    /**
     * Named constructor for a matrix of ParMooN-specific NSE Type 1.
     * The matrix takes the block structure
     *
     * ( A  0  0  B1T )
     * ( 0  A  0  B2T )
     * ( 0  0  A  B3T )
     * ( B1 B2 B3 0  )
     *
     * where B1T B2T and B3T are not explicitly stored (and marked non-transposed).
     *
     * How to use a named constructor? Have a look at the test file!
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE3D problems,
     * whose block structure is of NSE Type 1.
     */
    static BlockFEMatrix NSE3D_Type1(std::shared_ptr<const TFESpace3D> velocity,
                                     std::shared_ptr<const TFESpace3D> pressure);

    /**
     * Named constructor for a matrix of ParMooN-specific NSE Type 2.
     * The matrix takes the block structure
     *
     * ( A  0  0  B1T )
     * ( 0  A  0  B2T )
     * ( 0  0  A  B3T )
     * ( B1 B2 B3 0  )
     *
     * where B1T, B2T and B3T are explicitly stored (and marked non-transposed).
     *
     * How to use a named constructor? Have a look at the test file!
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE3D problems,
     * whose block structure is of NSE Type 2.
     */
    static BlockFEMatrix NSE3D_Type2(std::shared_ptr<const TFESpace3D> velocity,
                                     std::shared_ptr<const TFESpace3D> pressure);

    /**
     * Named constructor for a matrix of ParMooN-specific NSE Type 3.
     * The matrix takes the block structure
     *
     * ( A11  A12  A13  B1T )
     * ( A21  A22  A23  B2T )
     * ( A31  A32  A33  B3T )
     * ( B1   B2   B3   0   )
     *
     * where B1T, B2T and B3T are not explicitly stored (and marked non-transposed)..
     *
     * How to use a named constructor? Have a look at the test file!
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE3D problems,
     * whose block structure is of NSE Type 3.
     */
    static BlockFEMatrix NSE3D_Type3(std::shared_ptr<const TFESpace3D> velocity,
                                     std::shared_ptr<const TFESpace3D> pressure);

    /**
     * Named constructor for a matrix of ParMooN-specific NSE Type 4.
     * The matrix takes the block structure
     *
     * ( A11  A12  A13  B1T )
     * ( A21  A22  A23  B2T )
     * ( A31  A32  A33  B3T )
     * ( B1   B2   B3   0   )
     *
     * where B1T, B2T and B3T are explicitly stored (and marked non-transposed)..
     *
     * How to use a named constructor? Have a look at the test file!
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE3D problems,
     * whose block structure is of NSE Type 4.
     */
    static BlockFEMatrix NSE3D_Type4(std::shared_ptr<const TFESpace3D> velocity,
                                     std::shared_ptr<const TFESpace3D> pressure);

    /**
     * Named constructor for a matrix of ParMooN-specific NSE Type 14.
     * The matrix takes the block structure
     *
     * ( A11  A12  A13  B1T )
     * ( A21  A22  A23  B2T )
     * ( A31  A32  A33  B3T )
     * ( B1   B2   B3   C   )
     *
     * where B1T, B2T and B3T are explicitly stored (and marked non-transposed).
     *
     * How to use a named constructor? Have a look at the test file!
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE3D problems,
     * whose block structure is of NSE Type 14.
     */
    static BlockFEMatrix NSE3D_Type14(std::shared_ptr<const TFESpace3D> velocity,
                                      std::shared_ptr<const TFESpace3D> pressure);

    /**
     * Named constructor for a Mass matrix of ParMooN-specific NSE type 1 & 2
     *
     * @param velocity The velocity finite element space
     * @return A newly constructed BlockFEMatrix for
     * 3D NSE time dependent problems
     */
    static BlockFEMatrix Mass_NSE3D(std::shared_ptr<const TFESpace3D> velocity);

    /**
     * Named constructor for a mass matrix type 1, 2
     * The matrix takes the block structure
     *
     * ( M11  0  0   0)
     * ( 0  M22  0   0)
     * ( 0   0   M33 0)
     * ( 0   0   0   0)
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE2D problems,
     * whose block structure is of NSE Type 1, 2.
     */
    static BlockFEMatrix Mass_NSE3D_Type1(
      std::shared_ptr<const TFESpace3D> velocity,
      std::shared_ptr<const TFESpace3D> pressure);
    static BlockFEMatrix Mass_NSE3D_Type2(
      std::shared_ptr<const TFESpace3D> velocity,
      std::shared_ptr<const TFESpace3D> pressure);

    /**
     * Named constructor for a mass matrix type 3, 4
     * The matrix takes the block structure
     *
     * ( M11  M12  M13  0)
     * ( M21  M22  M23  0)
     * ( M31  M32  M33  0)
     * ( 0    0    0    0)
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE2D problems,
     * whose block structure is of NSE Type 3, 4.
     */
    static BlockFEMatrix Mass_NSE3D_Type3(
      std::shared_ptr<const TFESpace3D> velocity,
      std::shared_ptr<const TFESpace3D> pressure);
    static BlockFEMatrix Mass_NSE3D_Type4(
      std::shared_ptr<const TFESpace3D> velocity,
      std::shared_ptr<const TFESpace3D> pressure);

    /**
     * Named constructor for a mass matrix type 14
     * The matrix takes the block structure
     *
     * ( M11  M12  M13  MQ1 )
     * ( M21  M22  M23  MQ2 )
     * ( M31  M32  M33  MQ3 )
     * ( QM1  QM2  QM3    0 )
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE2D problems,
     * whose block structure is of NSE Type 14.
     */
    static BlockFEMatrix Mass_NSE3D_Type14(
      std::shared_ptr<const TFESpace3D> velocity,
      std::shared_ptr<const TFESpace3D> pressure);
#endif


    /**
     * Named constructor for a matrix for NSE type problems using a DG
     * discretization in 2D and 3D.
     * Creates a 2x2 block matrix of the form
     *  ( A  B1' )
     *  ( B2 C   ),
     * where A is velocity-velocity coupling, B1' velocity-pressure,
     * B2 pressure-velocity and C pressure-pressure coupling.
     *
     * How to use a named constructor? Have a look at the test file!
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE problems.
     */
    static BlockFEMatrix NSE_DG(std::shared_ptr<const FESpace> velocity,
                                 std::shared_ptr<const FESpace> pressure);


    /**
     * Add the active rows of a certain FEMatrix to several blocks at once.
     * This just figures out whether the adding will work, whether the
     * coloring scheme must be adapted and then delegates to FEMatrix
     * to perform the actual adding of active rows.
     *
     * @param[in] summand The FEMatrix to be added.
     * @param[in] factor The factor by which to scale it.
     * @param[in] cell_positions Where to add the matrix.
     * @param[in] transposed_states In which transposed state to do it.
     */
    void add_matrix_actives(
        const FEMatrix& summand, double factor,
        const std::vector<std::vector<size_t>>& cell_positions,
        const std::vector<bool>& transposed_states);

    /** @brief compute y = Ax, respecting non-active degrees of freedom.
     *
     * @param[in] x the BlockVector which is multiplied by this matrix
     * @param[out] y result of matrix-vector-multiplication
     */
    virtual void apply(const BlockVector& x, BlockVector& y) const override;

    /** @brief Compute y = y + a * Ax
     *
     * Add the matrix-vector product "Ax", scaled by "a", to y. "A" is this
     * matrix. Non-Active degrees of freedom are taken into account, i.e. the
     * matrix is treated just the way it would be treated as if the non-
     * active rows of the blocks would have been assembled globally correct.
     *
     * @param x the BlockVector which is multiplied by this matrix
     * @param y Gets added the result of the scaled matrix-vector-multiplication "aAx"
     * @param a optional factor, defaults to 1.0
     */
    virtual void apply_scaled_add(const BlockVector & x, BlockVector & y,
                                  double a = 1.0) const override;

    /**
     * Compute y = y + a * Ax, but leave y's non-active entries untouched.
     * Methods of this type are needed in the assembling process due to ParMooN's
     * handling of Dirchlet dofs.
     *
     * @param x the BlockVector which is multiplied by this matrix
     * @param y Gets added the result of the scaled matrix-vector-multiplication "aAx"
     * @param a optional factor, defaults to 1.0
     */
    void apply_scaled_add_actives(const BlockVector & x, BlockVector & y,
                                  double a = 1.0) const;

    /** @brief Compute y = y + a * Ax
     *
     * TODO This method is at the moment used in class Time_NSE2D. There it is
     * necessary at one point to multiply a vector  with a sub-matrix (velo-velo
     * blocks) only. There might a be a better approach to the issue, it would
     * be nice to get rid of this somewhat strange method.
     * THE METHOD IS NOT CONTAINED IN THE UNIT TEST.
     *
     * Add the matrix-vector product of the sub-block-matrix which consists of
     * block 0 to sub_row X 0 to sub_col, "Ax", scaled by "a", to y.
     * Only active rows are taken into account, non-actives in y are left
     * unchanged.
     *
     * @param x the BlockVector which is multiplied by this matrix
     * @param y Gets added the result of the scaled matrix-vector-multiplication "aAx"
     * @param sub_row number or rows of the submatrix
     * @param sub_col number of cols of the submatrix
     * @param a optional factor, defaults to 1.0
     */
    void apply_scaled_submatrix(const BlockVector & x, BlockVector & y,
                                size_t sub_row, size_t sub_col,
                                double a = 1.0) const;

    /**
     * Used as a developmental tool to discover slicing,
     * there should be no reason to use it anymore when the class is finished.
     * Checks if all TMatrix smart pointers stored in the
     * base class can be painlessly casted into smart pointers
     * to FEMatrix. Complains if not so, but does not throw.
     */
    void check_pointer_types();

    /*! Check whether a BlockVector b is fit to be the
     * rhs b of the equation Ax=b. (Including actives check.)
     * @param[in] b Rhs b of the equation Ax=b.
     */
    virtual void check_vector_fits_image(const BlockVector& b) const override;

    /*! Check whether a BlockVector x is fit to be the
     * factor x of the equation Ax=b. (Including actives check.)
     * @param[in] x Factor x in the equation Ax=b.
     */
    virtual void check_vector_fits_pre_image(const BlockVector& x) const override;

    /**
     * Get the ansatz space of a certain grid cell. Since the ansatzspace is
     * identical for all cells in a column, the input cell_row is not actually needed.
     * It is only requested so that the user does not have to think about
     * whether it's rowwise or columnwise constant every time she calls this method...
     *
     * @param[in] cell_row The cell row of the grid cell.
     * @param[in] cell_column The cell column of the grid cell.
     * @return The ansatzspace, which is the same for the entire column.
     */
#ifdef __2D__
    std::shared_ptr<const TFESpace2D> get_ansatz_space(size_t cell_row,
                                                       size_t cell_column) const;
#elif __3D__
    std::shared_ptr<const TFESpace3D> get_ansatz_space(size_t cell_row,
                                                       size_t cell_column) const;
#endif

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
    std::shared_ptr<const FEMatrix> get_block(
        size_t cell_row, size_t cell_col, bool& is_transposed) const;

    /**
     * @brief this is a non-const version of the method with the same name.
     */
    std::shared_ptr<FEMatrix> get_block(size_t cell_row, size_t cell_col,
                                        bool& is_transposed);

    /**
     * This method is the main interface to ParMooN solving procedures which
     * depend on the block structure.
     * It traverses the whole matrix from left to right, top to bottom, and fills
     * the return vector with pointers to all blocks, regardless of the storage
     * structure.
     *
     * The classes which make use of this method must decide how to proceed
     * with the gained blocks. They must make use of a priori knowledge
     * about the block structure of this matrix, or of methods which check
     * the block structure.
     *
     * @note The method returns shared_pointers which is not optimal,
     * for the solvers should not share ownership. But since most solvers
     * in ParMooN still expect raw pointer and one can't get those from weak
     * pointers, we chose sharedp pointers as return values.
     *
     * @return A list of pointers to the blocks. Left to right, top to bottom,
     * all blocks appear as often as they appear in cells.
     */
    std::vector<std::shared_ptr<const FEMatrix>> get_blocks() const;

    /**
     * This method is the main interface to ParMooN assembling procedures.
     * It traverses the whole matrix from left to right, top to bottom, and fills
     * the return vector with each new block it finds.
     *
     * This results in an output vector which contains each stored matrix
     * exactly once. That is the only thing the method cares about: that
     * no stored block is handed out twice.
     *
     * The classes which make use of this method must decide how to proceed
     * with the gained blocks. They must make use of a priori knowledge
     * about the block structure of this matrix, or of methods which check
     * the block structure.
     *
     * By default the method skips all those blocks which have a zero-map
     * Structure. To include them, give "true" as an input parameter.
     *
     * @param[in] include_zeroes Whether to include blocks with zero-map
     * Structure or skip them. Defaults to "false" (skipping zero-map blocks).
     * @return A list of pointers to the blocks. Each (relevant) blocks appears once.
     */
    std::vector<std::shared_ptr<FEMatrix>> get_blocks_uniquely(bool include_zeroes = false);

    /**
     * Acts just as get_blocks_uniquely(bool), but on only those
     * matrix cells whose coordinates are given in input vector "cells".
     *
     * @param[in] cells The positions of the cells whose blocks we are interested in.
     * should be a vector of double pairs, obviously, and not out of range.
     * @param[in] include_zeroes Whether to include blocks with zero-map
     * Structure or skip them. Defaults to "false" (skipping zero-map blocks).
     * @return A list of pointers to the matrix blocks. Each (relevant) blocks appears once.
     */
    std::vector<std::shared_ptr<FEMatrix>> get_blocks_uniquely(
      const std::vector<std::vector<size_t>>& cells,
      bool include_zeroes = false);

    /// @return The column (means: ansatz-)space of a certain cell column.
#ifdef __2D__
    std::shared_ptr<const TFESpace2D> get_column_space(size_t cell_column) const
#elif __3D__
    std::shared_ptr<const TFESpace3D> get_column_space(size_t cell_column) const
#endif
    {
      return ansatz_spaces_columnwise_.at(cell_column);
    }

    /** @brief return this BlockMatrix as one TMatrix
     *
     * This returns a merged version of this matrix. Note that the merged
     * matrix does not get stored internally, for it cannot easily be kept
     * up to date, but is recreated on every call.
     *
     * Treats Dirichlet rows correctly and globally, regardless of
     * what the particular blocks hold in their Dirichlet rows.
     *
     * Usually this is used to pass this matrix to a solver.
     */
    virtual std::shared_ptr<TMatrix> get_combined_matrix() const override;

    /// Returns a TMatrix representing a rectangular submatrix.
    /// Contains Dirichlet handling and pressure correction - the whole shebang.
    virtual std::shared_ptr<TMatrix> get_combined_submatrix(
        std::pair<size_t,size_t> upper_left,
        std::pair<size_t,size_t> lower_right) const override;

    double estimate_spectral_radius(
      const double* left_multipliers = nullptr,
      const double* right_multipliers = nullptr) const;

#ifdef _MPI
    /// Return a list of the FE communicators belonging to the FESpaces of
    /// the rows columns.
    std::vector<const TParFECommunicator3D*> get_communicators() const;

    bool is_additive() const { return additive_storage; };
    void set_additive(bool additive) { additive_storage = additive; };
#endif

    /**
     * Spawns a new BlockFEMatrix taken from the Block diagonal of this,
     * maintaining the coloring pattern.
     *
     * E.g. if we have a block fe matrix structure like
     *  A A A B^T
     *  A A A B^T
     *  A A A B^T
     *  B B B C,
     *  calling get_sub_blockfematrix(0,1) will give
     *
     *  A A
     *  A A
     *
     *  and calling get_sub_blockfematrix(1,3) will give
     *
     *  A A B^T
     *  A A B^T
     *  B B C.
     *
     * @param first The upper-leftmost diagonal block to include.
     * @param last The lower-rightmost diagonal block to include.
     */
    BlockFEMatrix get_sub_blockfematrix(
        size_t first, size_t last) const;

    /**
     * This method returns the number of actives of a certain cell column's
     * ansatz-space. It is needed only for the templated constructor of
     * BlockVector, and will be removed as soon as ParMooN has a new
     * handling of non-actives.
     *
     * @param cell_column The cell column.
     * @return The number of actives of the column's ansatzspace.
     */
    size_t get_n_column_actives(size_t cell_column) const;

    /**
     * This method returns the number of actives of a certain cell row's
     * test-space. It is needed only for the templated constructor of
     * BlockVector, and wil be removed as soon as ParMooN has a new
     * handling of non-actives.
     *
     * @param cell_row The cell row.
     * @return The number of actives of the row's testspace.
     */
    size_t get_n_row_actives(size_t cell_row) const;

    size_t get_n_column_inner(size_t cell_column) const;
    size_t get_n_row_inner(size_t cell_row) const;


    /// @return The row (means: test-)space of a certain cell row.
#ifdef __2D__
    std::shared_ptr<const TFESpace2D> get_row_space(size_t cell_row) const
#elif __3D__
    std::shared_ptr<const TFESpace3D> get_row_space(size_t cell_row) const
#endif
    {
      return test_spaces_rowwise_.at(cell_row);
    }

    /**
     * Get the test space of a certain grid cell. Since the testspace is
     * identical for all cells in a row, the input cell_column is not needed.
     * It is only requested so that the user does not have to think about
     * whether it's rowwise or columnwise constant every time she calls this method...
     *
     * @param[in] cell_row The cell row of the grid cell.
     * @param[in] cell_column The cell column of the grid cell.
     * @return The testspace, which is the same for the entire row.
     */
#ifdef __2D__
    std::shared_ptr<const TFESpace2D> get_test_space(size_t cell_row,
                                                     size_t cell_column) const;
#elif __3D__
    std::shared_ptr<const TFESpace3D> get_test_space(size_t cell_row,
                                                     size_t cell_column) const;
#endif

    /** @brief this method is used to compare the number of actives in a block vector
     * to the number of actives in test space
     *  @param nActive number of actives
     *  @param spaceNumber number of the test space to compare the actives
     */
    virtual void handle_discovery_of_vector_actives(const int nActive,
                                                    const int spaceNumber) const;

    /**
     * Print information on the matrix. Works in MPI case, too.
     * So far it prints only the block dimensions and the total number of d.o.f.
     * (MPI case: global number of d.o.f.)
     *
     * TODO Might be extended by further info and depend on verbosity.
     */
    void print_matrix_info(const std::string& name) const;

    /**
     * Overrides the method from the base class
     * and does nothing but print an error, when called. This ensures, that
     * no TMatrix is put into a BlockFEMatrix.
     */
    virtual void replace_blocks(
        const TMatrix& new_block,
        const std::vector<std::vector<size_t>>& cell_positions,
        const std::vector<bool>& transposed_states) override;

    /**
     * Replace the blocks in a set of cells with an FEMatrix.
     * The new block's test- and ansatzspace must fit those of
     * the grid positions where it should be placed (resp. transposed
     * mode).
     * Also, if the testspace of one of the cells has non-active
     * degrees of freedom, the matrix may not be stored as transposed
     * there.
     * If either of these conditions is violated, the method will quit
     * the program with an error.
     *
     * @param[in] new_block The new block, a copy will be made.
     * @param[in] cell_positions The positions where the new block is supposed to go.
     * @param[in] transposed_states The transposed mode of the block in the cell positions.
     */
    virtual void replace_blocks(
        const FEMatrix& new_block,
        const std::vector<std::vector<size_t>>& cell_positions,
        const std::vector<bool>& transposed_states);

    /*!
     *  Scale active degrees of freedom in the blocks at the cells
     *  whose positions are given by cell_positions. The method just figures
     *  out which blocks are to be treated and whether the coloring
     *  scheme has to be adapted, then delegates the active-scaling to the
     *  FEMatrix class.
     *
     *  @param[in] factor The scaling factor
     *  @param[in] cell_positions Those cells whose blocks are to scale.
     */
    void scale_blocks_actives(
        double factor,
        const std::vector<std::vector<size_t>>& cell_positions );

    /// @brief turn on pressure projection manually
    void enable_pressure_projection() const;
    /// @brief turn off pressure projection manually
    void disable_pressure_projection() const;
    /// @brief find out if pressure correction is currently enabled
    bool pressure_projection_enabled() const;

    // Special member functions.

    /*! Copy constructor.
     *  @note Is a bit expensive because it copies each stored
     *  matrix twice - once in the base class copy constructor
     *  and once in its body to get the cast to FEMatrix right.
     */
    BlockFEMatrix(const BlockFEMatrix&);

    /// Move constructor.
    BlockFEMatrix(BlockFEMatrix&&);

    /// Unified assignment operator. Uses copy-and-swap.
    BlockFEMatrix& operator=(BlockFEMatrix);

    /** Swap function used for copy-and swap in copy assignment.
     * @param[in,out] first The object to be swapped with second.
     * @param[in,out] second The object to be swapped with first.
     */
    friend void swap(BlockFEMatrix& first, BlockFEMatrix& second);

    /// @brief Destructor. Tidies up nice and clean.
    virtual ~BlockFEMatrix() = default;

    /// @brief add the blockfematrix to this
    void add_blockfe_matrix(const BlockFEMatrix& Matrix, double factor = 1.);

  protected:
#ifdef __2D__
    /// Store pointers to the testspaces rowwise.
    std::vector<std::shared_ptr<const TFESpace2D>> test_spaces_rowwise_;
    /// Store pointers to the ansatzspaces columnwise.
    std::vector<std::shared_ptr<const TFESpace2D>> ansatz_spaces_columnwise_;
#elif __3D__
    /// Store pointers to the testspaces rowwise.
    std::vector<std::shared_ptr<const TFESpace3D>> test_spaces_rowwise_;
    /// Store pointers to the ansatzspaces columnwise.
    std::vector<std::shared_ptr<const TFESpace3D>> ansatz_spaces_columnwise_;
#endif

    /// @brief modify the first pressure row for (Navier-) Stokes matrices
    ///
    /// This is sometimes necessary if there are only Dirichlet boundaries to
    /// avoid a singular matrix. In the methods `get_combined_matrix` this will
    /// change the firsrt pressure row to be zero except on the diagonal.
    ///
    /// It's not nice that this is mutable, but usually the call to the method
    /// 'get_combined_matrix' is made on a const BlockFEMatrix. So the setter
    /// methods `enable_pressure_correction`/`disable_pressure_projection` must
    /// be const as well.
    mutable bool use_pressure_projection_;

#ifdef _MPI
    bool additive_storage;
#endif

  private:
    /**
     * Actual implementation of add scaled actives method, whose interface is given
     * in the public part.
     */
    void add_scaled_actives(
        const FEMatrix& summand, double scaling_factor,
        std::vector<grid_place_and_mode>& row_column_transpose_tuples);

    /**
     * Try if a given TMatrix can be cast to an FEMatrix (this meaning the object actually is
     * an FEMatrix) and if so, make an FEMatrix copy of it and return a smart pointer to that,
     * so that it can be stored as a shared pointer to TMatrix in the baseclass.
     *
     * This quirky method is needed due to the base class storing TMatrices only,
     * but this class dealing with FEMatrices.
     */
    virtual std::shared_ptr<TMatrix> create_block_shared_pointer(const TMatrix& block) const override;

    /**
     * Actual implementation of the scale actives method, whose interface is given
     * in the public part.
     */
    void scale_blocks_actives(
      double scaling_factor,
      std::vector<grid_place_and_mode>& row_column_transpose_tuples);
};

#endif /* USER_PROJECTS_BLOCKFEMATRIX_H_ */
