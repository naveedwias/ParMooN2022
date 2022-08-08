#ifndef __ITERATION_RICHARDSON__
#define __ITERATION_RICHARDSON__

#include <IterativeMethod.h>

template <class LinearOperator, class Vector>
class Iteration_richardson : public IterativeMethod<LinearOperator, Vector>
{
  public:
    /** constructor */
    explicit Iteration_richardson(std::shared_ptr<Preconditioner<Vector>> p);
    
    /** destructor */
    virtual ~Iteration_richardson() = default;
    
    /** iterate routine */
    std::pair<unsigned int, double> iterate(const LinearOperator & A, 
                                            const Vector & rhs,
                                            Vector & solution) override final;
};

#endif // __ITERATION_RICHARDSON__
