/**
 * @class FEMatrix
 * @brief matrix which represents the mapping between two finite element spaces
 * 
 * 
 * This class essentially implements constructors which take finite element
 * spaces. Also it saves which spaces were used during construction. Otherwise
 * it behaves just like its base class Matrix.
 * 
 * Note that there is no constructor taking a TStructure together with the 
 * spaces it was created with. If you want to create multiple matrices with the
 * same structure object, use the copy constructor which will do just that.
 * 
 * @ruleof0
 * 
 * @todo write a test for this class
 */

#ifndef __FEMATRIX__
#define __FEMATRIX__

#include <FESpace1D.h>
#include <FESpace2D.h>
#include <FESpace3D.h>
#include <Matrix.h>

class FEMatrix : public TMatrix
{
  public:
    /// @brief there is no default constructor, this makes no sense
    FEMatrix() = delete;
    
    /// @name construct square matrix using one finite element space
    /// @brief ansatz and test space is the same 
    //@{
    explicit FEMatrix(std::shared_ptr<const TFESpace1D> space);
    explicit FEMatrix(std::shared_ptr<const TFESpace2D> space);
    #ifdef __3D__
    explicit FEMatrix(std::shared_ptr<const TFESpace3D> space);
    #endif // 3D
    //@}
    
    /// @name construct square matrix using one finite element space with a 
    ///       given TMatrix
    /// @brief ansatz and test space is the same 
    //@{
    explicit FEMatrix(std::shared_ptr<const TFESpace2D> space, const TMatrix&);
    #ifdef __3D__
    explicit FEMatrix(std::shared_ptr<const TFESpace3D> space, const TMatrix&);
    #endif // 3D
    //@}
    
    /// @name construct rectangular matrix using two finite element spaces
    /// @brief test and ansatz space are possibly different
    /// 
    /// A matrix using this structure represents a linear map from the ansatz 
    /// to the test space.
    /// @param[in] is_empty If true, this TMatrix is created with an empty structure.
    //@{
    FEMatrix(std::shared_ptr<const TFESpace2D> testspace,
             std::shared_ptr<const TFESpace2D> ansatzspace,
             bool is_empty = false);
    #ifdef __3D__
    FEMatrix(std::shared_ptr<const TFESpace3D> testspace,
             std::shared_ptr<const TFESpace3D> ansatzspace,
             bool is_empty = false);
    #endif // 3D
    //@}
    
    /// @name construct matrix with given structure and finite element space
    /// @brief ansatz and test space is the same. The structure must be square.
    //@{
    FEMatrix(std::shared_ptr<const TFESpace1D> space,
             std::shared_ptr<TStructure> structure);
    FEMatrix(std::shared_ptr<const TFESpace2D> space,
             std::shared_ptr<TStructure> structure);
    #ifdef __3D__
    FEMatrix(std::shared_ptr<const TFESpace3D> space,
             std::shared_ptr<TStructure> structure);
    #endif // 3D
    //@}
    
    //! Default copy constructor. Performs shallow copy.
    FEMatrix(const FEMatrix&) = default;
    
    //! Default move constructor.
    FEMatrix(FEMatrix&&) = default;
    
    //! Default copy assignment operator. Performs shallow copy.
    FEMatrix& operator=(const FEMatrix&) = default;
    
    //! Default move assignment operator
    FEMatrix& operator=(FEMatrix&&) = default;
    
    //! Default destructor.
    ~FEMatrix() = default;
    
    
    /** @brief reset all entries in active rows to zero */
    void resetActive();

    /** @brief set zeros in nonactive rows (including hanging rows). 
     * 
     * This is e.g. for the off-diagonal blocks in a Stokes matrix 
     */
    void resetNonActive();

    /** @brief scale this matrix by a factor
     * 
     * Only rows corresponding to active d.o.f are scaled. Other rows remain
     * unscaled.
     */
    void scaleActive(double factor = 1.0);
    
    /// @brief scale the diagonal entries in the nonactive rows
    void scale_non_active_diagonals(double factor);
    
    /// @brief set the diagonal entries in all Dirichlet rows to 1.0.
    void set_dirichlet_diagonals();

    /** @brief adding a scaled matrix to this matrix
     * 
     * This is only done for those rows which correspond to active degrees of 
     * freedom.
     * 
     * The summation is index-wise, i.e. A(i,j) += factor*m(i.j), where A is 
     * this matrix. 
     * 
     * Note that this only works if the sparsity structure is the same for this
     * matrix and m.
     */
    void addActive(const FEMatrix& m, double factor = 1.0);
    
    /** @brief compute y = y + a * Ax 
     *
     * add the matrix-vector product "Ax", scaled by "a", to y: only active
     * "A" is this matrix.
     * 
     * This function can be used to compute the residual r = b - Ax.
     *
     * @param x the vector which is multiplied by this matrix
     * @param y result of matrix-vector-multiplication and scaling
     * @param factor optional scaling   factor, default to 1.0
     */
    void multiplyActive(const double *x, double *y, double factor = 1.0)
      const;
      
    
    /** @brief Add the coupling from Hanging rows to the system matrix
     * 
     * Rows corresponding to the coupeld dofs with hanging dofs are modified.
     * Earlier it was part of Assemble. Helpful in implementation of internal
     * full matrix structure.
     * Input argument added for AFC schemes and grids with hanging nodes
     * 
     */  
    void ModifyMatrixAccordingToCoupling(bool assemble_dirichlet_rows);

    
    /** @brief Add the coupling from Hanging rows to the system matrix
     * 
     * This is a specific routine made for AFC schemes with hanging nodes.
     * Currently it is a work in progress and follows the idea from ongoing
     * work with Petr Knobloch.
     * 
     */  
    void ModifyMatrixAccordingToCouplingAFC();   
    
    /**
     * @brief consistent update to the hanging rows.
     * 
     * This is needed to ensure that the pre-image space consists of continuous
     * functions even in the case of hanging nodes.
     * 
     * The hanging rows should contain specific values which only depend on the
     * polynomial degree, the hanging node type (and possibly more), but not on
     * assembled values.
     */
    void correct_hanging_rows();
    
    /**
     *  compute y = y + a*A^T x
     *
     *  add the matrix-vector product "A^T x", scaled by "a", to y: only active
     * "A" is this matrix.
     * @param x
     * @param y
     * @param factor
     */
    void multiplyTransposedActive(const double *x, double *y, double factor = 1.0)
      const;

    
    /** @brief return the number of active rows */
    int get_n_active_rows() const;
    
    /** @brief return 1D test space */
    std::shared_ptr<const TFESpace1D> GetTestSpace1D() const;
    
    /** @brief return 1D ansatz space */
    std::shared_ptr<const TFESpace1D> GetAnsatzSpace1D() const;
    
    /** @brief return 2D test space */
    std::shared_ptr<const TFESpace2D> GetTestSpace2D() const;
    
    /** @brief return 2D ansatz space */
    std::shared_ptr<const TFESpace2D> GetAnsatzSpace2D() const;
    
    #ifdef __3D__
    /** @brief return 3D test space */
    std::shared_ptr<const TFESpace3D> GetTestSpace3D() const;
    
    /** @brief return 3D ansatz space */
    std::shared_ptr<const TFESpace3D> GetAnsatzSpace3D() const;
    #endif // 3D
    
    /** @brief return test space */
    std::shared_ptr<const TFESpace> GetTestSpace() const;
    
    /** @brief return ansatz space */
    std::shared_ptr<const TFESpace> GetAnsatzSpace() const;
    
    /** return FESpace */
    std::shared_ptr<const TFESpace1D> GetFESpace1D() const;
    /** return FESpace */
    std::shared_ptr<const TFESpace2D> GetFESpace2D() const;
    #ifdef __3D__
    /** return FESpace */
    std::shared_ptr<const TFESpace3D> GetFESpace3D() const;
    #endif // 3D
    
  private:
    /// @name ansatz spaces
    /// @brief the ansatz space (pre-image space)
    /// @details Exactly one of these pointers is not a nullptr.
    /// @todo make this a share_ptr
    //@{
    std::shared_ptr<const TFESpace1D> AnsatzSpace1D;
    std::shared_ptr<const TFESpace2D> AnsatzSpace2D;
    std::shared_ptr<const TFESpace3D> AnsatzSpace3D;
    //@}
    
    /// @name test spaces
    /// @brief the test space (image space)
    /// @details Exactly one of these pointers is not a nullptr.
    /// @todo make this a share_ptr
    //@{
    std::shared_ptr<const TFESpace1D> TestSpace1D;
    std::shared_ptr<const TFESpace2D> TestSpace2D;
    std::shared_ptr<const TFESpace3D> TestSpace3D;
    //@}
};

#endif //__FEMATRIX__
