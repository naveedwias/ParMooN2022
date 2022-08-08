/**
 * @file JacobiSmoother.h
 * Wrap up a Jacobi step as a smoother for multigrid method.
 *
 * @date 2016/05/13
 * @author Clemens Bartsch
 */

#ifndef INCLUDE_MULTIGRID_JACOBISMOOTHER_H_
#define INCLUDE_MULTIGRID_JACOBISMOOTHER_H_

#include <Smoother.h>

#include <memory>

class BlockVector;
class BlockFEMatrix;
template <class LinearOperator, class Vector> class Iteration_jacobi;

/**
 * To use the Jacobi iteration as a smoother in multigrid framework, this
 * class wraps up an "Iteration_jacobi" object as a smoother.
 * Note that it will throw if you try to use it with a matrix with zeroes on
 * the diagonals.
 */
class JacobiSmoother : public Smoother
{
  public:

    /// Constructor takin only a damping parameter.
    /// @param[in] damp Damping for update. Must lie in (0,1].
    explicit JacobiSmoother(double damp);

    /// Apply one step of Jacobi iteration with the stored operator
    /// to the iterate "solution"
    void smooth(const BlockVector& rhs, BlockVector& solution) override;

    /// Reset the stored Iteration_jacobi object. TODO As soon as system classes
    /// storeshared pointers to their BlockFEMatrices, so should the Jacobi_iteration
    /// class!
    void update(const BlockFEMatrix& matrix) override;

    /* ************* *
     * Special member functions. Declared, but not all are defined yet.
     * ************* */
    //! Default constructor.
    JacobiSmoother();

    //! Copy constructor.
    JacobiSmoother( const JacobiSmoother& );

    //! Move constructor.
    JacobiSmoother( JacobiSmoother&& );

    //! Copy assignment operator.
    JacobiSmoother& operator=( const JacobiSmoother& );

    //! Move assignment operator.
    JacobiSmoother& operator=( JacobiSmoother&& );

    ~JacobiSmoother() = default;


  private:
    /// The iterative method one step of which is applied in "smooth".
    std::shared_ptr<Iteration_jacobi<BlockFEMatrix, BlockVector>> jacobi_;

    /// A damping factor in (0,1], damping the update.
    double damp_;
};

#endif /* INCLUDE_MULTIGRID_JACOBISMOOTHER_H_ */
