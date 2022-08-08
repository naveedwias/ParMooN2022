#ifndef __PRECONDITIONER_OPTIMIZED__
#define __PRECONDITIONER_OPTIMIZED__

#include <Preconditioner.h>
#include <memory>

template <class LinearOperator, class Vector>
class Preconditioner_optimized : public Preconditioner<Vector>
{
  public:
    /** constructor */
    explicit Preconditioner_optimized(const LinearOperator &A,
      std::shared_ptr<Preconditioner<Vector>> prec, double min_d, double max_d);

    /** destructor */
    virtual ~Preconditioner_optimized() = default;

    void apply(const Vector & z, Vector & r) const override final;

    void update() override final;

  private:
    const LinearOperator& A;
    std::shared_ptr<Preconditioner<Vector>> prec;
    double min_d, max_d;
};

#endif // __PRECONDITIONER_OPTIMIZED__
