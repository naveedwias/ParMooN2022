#ifndef __ITERATION_CG__
#define __ITERATION_CG__

#include <IterativeMethod.h>

template <class LinearOperator, class Vector>
class Iteration_cg : public IterativeMethod<LinearOperator, Vector>
{
  public:
    /** constructor */
    explicit Iteration_cg(std::shared_ptr<Preconditioner<Vector>> p);
    
    /** destructor */
    virtual ~Iteration_cg() = default;
    
    /** iterate routine */
    std::pair<unsigned int, double> iterate(const LinearOperator & A, 
                                            const Vector & rhs,
                                            Vector & solution) override final;
};

#endif // __ITERATION_CG__
