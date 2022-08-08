#ifndef __ITERATION_CGS__
#define __ITERATION_CGS__

#include <IterativeMethod.h>

template <class LinearOperator, class Vector>
class Iteration_cgs : public IterativeMethod<LinearOperator, Vector>
{
  public:
    /** constructor */
    explicit Iteration_cgs(std::shared_ptr<Preconditioner<Vector>> p);
    
    /** destructor */
    virtual ~Iteration_cgs() = default;
    
    /** iterate routine */
    std::pair<unsigned int, double> iterate(const LinearOperator & A, 
                                            const Vector & rhs,
                                            Vector & solution) override final;
};

#endif // __ITERATION_CGS__
