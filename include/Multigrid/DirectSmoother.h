/**
 * @file DirectSmoother.h
 *
 * @date 2016/10/11
 * @author Clemens Bartsch
 */

#ifndef INCLUDE_MULTIGRID_DIRECTSMOOTHER_H_
#define INCLUDE_MULTIGRID_DIRECTSMOOTHER_H_

#ifdef _MPI
#include <MumpsWrapper.h>
#else
#include <DirectSolver.h>
#endif
#include <Smoother.h>

#include <memory>

/**
 * A direct solver object wrapped up as a Smoother. Use only on the coarsest
 * grid - if this is enforced on a fine grid everything which happens on a
 * coarser grid is entirely useless, as the direct solver erases the
 * solution - instead of using the prolongated coarse grid solution as an initial
 * iterate, as a proper smoother does.
 */
class DirectSmoother : public Smoother
{
  public:
    // Solve the system Ax=b with the stored matrix and direct solver.
    virtual void smooth(const BlockVector& rhs, BlockVector& solution ) override;

    // Reset the stored "DirectSolver" object. Must be called whenever the
    // matrix with which to solve has changed.
    virtual void update(const BlockFEMatrix&) override;

    /* ************* *
     * Special member functions. Declared but not defined, since it
     * is not yet clear whether to shallow or deep copy here.
     * ************* */
    //! Default constructor.
    DirectSmoother();

    //! Copy constructor.
    DirectSmoother( const DirectSmoother& );

    //! Move constructor.
    DirectSmoother( DirectSmoother&& );

    //! Copy assignment operator.
    DirectSmoother& operator=( const DirectSmoother& );

    //! Move assignment operator.
    DirectSmoother& operator=( DirectSmoother&& );

    ~DirectSmoother() = default;

  private:
#ifdef _MPI
    std::shared_ptr<MumpsWrapper> solver_;
    //This has to be stored due to a deprecated feature of the MUMPS Wrapper,
    // which requires the comms to be passed to every solve.
    std::vector<const TParFECommunicator3D*> comms_;
#else
    std::shared_ptr<DirectSolver> solver_;
#endif

};



#endif /* INCLUDE_MULTIGRID_DIRECTSMOOTHER_H_ */
