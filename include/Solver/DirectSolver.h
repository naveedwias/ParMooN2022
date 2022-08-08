/** ************************************************************************ 
*
* @class DirectSolver
* @brief solving a linear system using a direct solver
*
* Given a BlockMatrix the constructor of this class computes a factorization of
* that matrix. Then the method DirectSolver::solve enables the actual solving
* step which can be called multiple time, reusing the computed factorization.
* 
* Note that there is no way to update the matrix, i.e. to recompute a 
* factorization. In case your matrix changes, you should also create a new
* object of this class DirectSolver.
* 
* @ruleof0
*
****************************************************************************/

#ifndef __DIRECTSOLVER__
#define __DIRECTSOLVER__
#include <memory>
#include <vector>

class TMatrix;

// forward declaration
class BlockMatrix;
class BlockFEMatrix;
class BlockVector;
constexpr size_t pardiso_options_array_length = 64;

template <class V>
class CompositeOperator;

class DirectSolver
{
  public:
    enum class DirectSolverTypes {umfpack, pardiso};
    /**
     * @brief compute the factorization of a matrix, ready to call solve
     * 
     * This calls the other (private) constructor with 
     * matrix.get_combined_matrix().
     *
     * @param  matrix the matrix A where Ax=b
     */
    DirectSolver(const BlockMatrix& matrix, DirectSolverTypes type);
    DirectSolver(const BlockFEMatrix& matrix, DirectSolverTypes type);
    
    /** @brief compute the factorization of a matrix, ready to call solve
     * 
     * This makes a copy of the matrix and calls the private constructor.
     */
    DirectSolver(const TMatrix& matrix, DirectSolverTypes type);

    DirectSolver(const CompositeOperator<BlockVector>& matrix,
      DirectSolverTypes type);

    /** @brief This class is not copy constructible */
    DirectSolver(const DirectSolver&) = delete;

    /** @brief move constructor */
    DirectSolver(DirectSolver&&);

    /** @brief This class is not copy assignable */
    DirectSolver& operator=(const DirectSolver&) = delete;

    /** @brief move assignment operator */
    DirectSolver& operator=(DirectSolver&&);
    
    /** @brief release all memory */
    ~DirectSolver();

    /**
     * @brief Solves the equation A*(solution)=rhs for solution.
     * 
     * This calls the other method solve.
     *
     * @param rhs the right-hand side of the problem Ax=b
     * @param solution vector to store the solution into
     */
    void solve(const BlockVector& rhs, BlockVector& solution);
    
    /**
     * @brief Solves the equation A*(solution)=rhs for solution.
     * 
     * The computed solution is stored in the provided array solution.
     * 
     * @note Please use the other method solve taking BlockVectors. Here it can
     * not be checked if the pointers really point to valid arrays.
     *
     * @param rhs the right-hand side of the problem Ax=b
     * @param solution vector to store the solution into
     */
    void solve(const double* rhs, double* solution);
    
  private:
    /** @brief type of direct solver used */
    DirectSolverTypes type;
    
    /** @brief the matrix of the linear equation A*x=b */
    std::shared_ptr<TMatrix> matrix;
    
    /** @brief temporary variables for the matrix indices 
     * 
     * These are only used if the matrix size exceeds a certain threshold value
     * and umfpack is used. In that case we call the umfpack routines which use
     * long int instead of int. This avoids out of memory errors in umfpack.
     * 
     * The DirectSolver::cols and DirectSolver::rows are replacements of the 
     * respective matrix structure members, see Structure.h.
     * 
     * This should not be the final solution. We need the structure to also
     * handle long. Implementing this is quite intrusive, and up to now it is
     * only needed for umfpack.
     */
    std::vector<long int> cols;
    std::vector<long int> rows;
    
    /** @brief storage for umfpack direct solver */
    //@{
    void* symbolic;
    void* numeric;
    //@}
    
    /** @brief storage for pardiso direct solver */
    //@{
    void *pt[pardiso_options_array_length]; /// memory pointers for pardiso_
    int maxfct;    /// maximum number of numerical factorization to perform
    int mnum;      /// number of matrices
    int mtype;     /// type of the matrix (ref. @b PARDISO_MTYPE_*)
    int perm;      /// User supplied permutation to apply in advance
    int nrhs;      /// number of right-hand sides
    int iparm[pardiso_options_array_length]; /// contains the pardiso settings
    int msglvl;    /// message level (0 - no messages, 1 - verbose)
    //@}
    
    /**
     * @brief compute the factorization of a matrix, ready to call solve
     * 
     * The indices are shifted to conform with Fortran in case of the Pardiso 
     * solver. This is why this method has to change the matrix and therefore
     * does not take a const TMatrix.
     *
     * @param  matrix  the matrix A where Ax=b
     */
    DirectSolver(std::shared_ptr<TMatrix> matrix, DirectSolverTypes type);
    
    /** @brief compute symbolic factorization */
    void symbolic_factorize();
    /** @brief compute numeric factorization (requires symbolic factorization) 
     */
    void numeric_factorize();
};

#endif // __DIRECTSOLVER__
