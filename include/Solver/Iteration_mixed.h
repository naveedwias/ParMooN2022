#ifndef __ITERATION_MIXED__
#define __ITERATION_MIXED__

#include <IterativeMethod.h>
#include <Iteration_bicgstab.h>
#include <Iteration_gmres.h>

// For now this alternates this->restart FGMRES and BiCGSTAB iterations.
// TODO: make this configurable once nested databases are available
template <class LinearOperator, class Vector>
class Iteration_mixed : public IterativeMethod<LinearOperator, Vector>
{
  public:
    /** constructor */
    explicit Iteration_mixed(std::shared_ptr<Preconditioner<Vector>> p);

    /** destructor */
    virtual ~Iteration_mixed() = default;

    /** iterate routine */
    std::pair<unsigned int, double> iterate(const LinearOperator & A,
                                            const Vector & rhs,
                                            Vector & solution) override final;
  private:
    std::shared_ptr<Iteration_bicgstab<LinearOperator, Vector>> bicgstab_iterator;
    std::shared_ptr<Iteration_gmres<LinearOperator, Vector>> gmres_iterator;
};

#endif // __ITERATION_MIXED__