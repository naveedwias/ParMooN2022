#ifndef __ITERATION_BICG__
#define __ITERATION_BICG__

#include <IterativeMethod.h>

template <class LinearOperator, class Vector>
class Iteration_bicgstab : public IterativeMethod<LinearOperator, Vector>
{
  public:
    /** constructor */
    explicit Iteration_bicgstab(std::shared_ptr<Preconditioner<Vector>> p);

    /** destructor */
    virtual ~Iteration_bicgstab() = default;

    /** iterate routine */
    std::pair<unsigned int, double> iterate(const LinearOperator & A, 
                                            const Vector & rhs,
                                            Vector & solution) override final;
};

#endif // __ITERATION_BICG__
