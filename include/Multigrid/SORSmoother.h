#ifndef INCLUDE_MULTIGRID_SORSMOOTHER_H_
#define INCLUDE_MULTIGRID_SORSMOOTHER_H_

/**
 * @file SORSmoother.h
 * Wrap up a SOR sweep as a smoother for multigrid methods.
 * SOR smoother can be used in MPI. It is, without any modifications, just
 * a processorwise SOR sweep.
 *
 * @date 2016/09/09
 * @author Clemens Bartsch
 */
#include <Smoother.h>

#include <memory>
#include <string>

class BlockVector;
class BlockFEMatrix;
template <class LinearOperator, class Vector> class Iteration_sor;

/**
 * To use the SOR iteration as a smoother in multigrid framework, this
 * class wraps up an "Iteration_sor" object.
 * Note that it will throw if you try to use it with a matrix with zeroes on
 * the diagonals.
 */
class SORSmoother : public Smoother
{
  public:

    /// Constructor, stores two interesting values which will then be handed
    /// over to each new sor_ object.
#ifdef _MPI
    SORSmoother(double omega, size_t sor_strat, const std::string& sor_par_strat);
#else
    SORSmoother(double omega, size_t sor_strat);
#endif

    /// Apply one sweep of SOR iteration with the stored operator
    /// to the iterate "solution".
    void smooth(const BlockVector& rhs, BlockVector& solution) override;

    /// Reset the stored Iteration_sor object. TODO As soon as system classes
    /// storeshared pointers to their BlockFEMatrices, so should the Iteration_sor
    /// class!
    void update(const BlockFEMatrix& matrix) override;

    /* ************* *
     * Special member functions. Declared, but not defined yet.
     * ************* */
    //! Default constructor.
    SORSmoother();

    //! Copy constructor.
    SORSmoother( const SORSmoother& );

    //! Move constructor.
    SORSmoother( SORSmoother&& );

    //! Copy assignment operator.
    SORSmoother& operator=( const SORSmoother& );

    //! Move assignment operator.
    SORSmoother& operator=( SORSmoother&& );

    ~SORSmoother() = default;


  private:
    /// The iterative method one step of which is applied in "smooth".
    std::shared_ptr<Iteration_sor<BlockFEMatrix, BlockVector>> sor_;
    ///
    double omega_;
    /// (0,1,2 - check it out in Iteration_sor)
    size_t sor_strat_;

#ifdef _MPI
    /// One of the strings "all_cells", "halo_0", "own_cells", each of which
    /// corresponds with a different parallelization strategy of the SOR.
    /// This feature is still in an experimental state.
    std::string parallel_strategy_;
#endif
};

#endif /* INCLUDE_MULTIGRID_SORSMOOTHER_H_ */
