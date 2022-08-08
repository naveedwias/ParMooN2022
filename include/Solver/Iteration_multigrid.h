/**
 * @file Interface for using multigrid as an iterative method or preconditioner
 * within the new Solver framework.
 *
 * @date 2016/05/10
 * @author Clemens Bartsch
 */

#ifndef INCLUDE_SOLVER_ITERATION_MULTIGRID_H_
#define INCLUDE_SOLVER_ITERATION_MULTIGRID_H_

#include <IterativeMethod.h>
#include <Multigrid.h>
#include <list>

//Forward declarations.
class Multigrid;

template <class LinearOperator, class Vector>
class Iteration_multigrid : public IterativeMethod<LinearOperator, Vector>,
                            public Preconditioner<Vector>
{
  public:
    /**
     * Constructor.
     * @param mg A multigrid object with all matrices ready.
     */
    explicit Iteration_multigrid(std::shared_ptr<Multigrid> mg);

    /**
     * For usage as preconditioner: do just one multigrid cycle.
     *
     * Should only be used in iterative methods which do not rely
     * on preconditioning being the same each iteration.
     *
     * @param[in] z Right hand side vector of the preconditioning equation.
     * @param[in, out] r Solution vector of the equation. Will be used as
     * start solution on input.
     */
    void apply(const Vector & z, Vector & r) const override final;

    /** iterate routine */
    std::pair<unsigned int, double> iterate(const LinearOperator & A,
                                            const Vector & rhs,
                                            Vector & solution) override final;
    
    void update(const LinearOperator&) override final
    { mg_->update(); };
    void update() override final
    { mg_->update(); };
    
    /** destructor */
    virtual ~Iteration_multigrid() = default;

  private:
    //The multigrid object which holds necessary grid information.
    mutable std::shared_ptr<Multigrid> mg_;

};




#endif /* INCLUDE_SOLVER_ITERATION_MULTIGRID_H_ */
