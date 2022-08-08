#ifndef __ITERATION_SOR__
#define __ITERATION_SOR__

#include <IterativeMethod.h>
#include <vector>

template <class LinearOperator, class Vector>
class Iteration_sor : public IterativeMethod<LinearOperator, Vector>,
                      public Preconditioner<Vector>
{
  public:
    /** @brief constructor, this uses NoPreconditioner<Vector> ! 
     * 
     * The sor method is basically a forward or backward sweep for triangular
     * parts of the linear operator. The flag determines forward (0), 
     * backward(1), or both(2). 
     * 
     * The overrelaxation parameter omega is between 0 and 2, in practice even 
     * between 1 and 2.
     */
    Iteration_sor(const LinearOperator & mat, int flag = 0, double omega = 1.0);
    
    /// @brief apply this object as a preconditioner
    void apply(const Vector & z, Vector & r) const override final;
    
    /// @brief Apply this object within a smoother for multigrid
    void apply_smoother(const Vector & z, Vector & r) const;

    /** @brief destructor */
    virtual ~Iteration_sor() = default;
    
    /** @brief iterate routine, called if this is used as a solver */
    std::pair<unsigned int, double> iterate(const LinearOperator & A, 
                                            const Vector & rhs,
                                            Vector & solution) override final;
    
    /// @brief return the used linear operator
    const LinearOperator& get_operator() const;

#ifdef _MPI
    /// @brief Choose a parallelization strategy for SOR. Still experimental.
    void set_parallel_strategy(const std::string& parallel_strategy);
#endif

  protected:
    /// @brief the type of sor: forward(0), backward(1), both(2)
    int sor_type;
    /// @brief the relaxation parameter
    double omega;
    /// @brief the linear operator (usually the matrix)
    const LinearOperator & linear_operator;

#ifdef _MPI
    /// One of the strings "all_cells", "halo_0", "own_cells", each of which
    /// corresponds with a different parallelization strategy of the SOR.
    /// This feature is still in an experimental state.
    std::string parallel_strategy_;
#endif
};

#endif // __ITERATION_JACOBI__
