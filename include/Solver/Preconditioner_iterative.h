#ifndef __PRECONDITIONER_ITERATIVE__
#define __PRECONDITIONER_ITERATIVE__

#include <IterativeMethod.h>

template <class LinearOperator, class Vector>
class Preconditioner_iterative : public Preconditioner<Vector>
{
  public:
    /** constructor */
    explicit Preconditioner_iterative(const LinearOperator &A,
      std::shared_ptr<IterativeMethod<LinearOperator, Vector>> solver);

    /** destructor */
    virtual ~Preconditioner_iterative() = default;

    void apply(const Vector & z, Vector & r) const override final;

    void update() override final;

  private:
    const LinearOperator& A;
    std::shared_ptr<IterativeMethod<LinearOperator, Vector>> solver;
};

#endif // __PRECONDITIONER_ITERATIVE__
