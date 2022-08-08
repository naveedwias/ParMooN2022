/**
 * @file New declaration of a level in a multigrid method.
 *
 * @date 2016/05/10
 * @author Clemens Bartsch
 */

#ifndef INCLUDE_MULTIGRID_MULTIGRIDLEVEL_H_
#define INCLUDE_MULTIGRID_MULTIGRIDLEVEL_H_

#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <Smoother.h>

#include <memory>

class ParameterDatabase;
class Smoother;

enum class SmootherCode{DIRECT_SOLVE, JACOBI, SOR, SSOR,
    NODAL_VANKA, CELL_VANKA, PATCH_VANKA,
    CELL_VANKA_JACOBI,
    NODAL_VANKA_STORE, CELL_VANKA_STORE, PATCH_VANKA_STORE,
};

/**
 * The multigrid level class. Its members are private, except for copy and move,
 * since its contents are entirely managed by the Multigrid class.
 *
 * FIXME It follows its matrix with a raw pointer. This should be changed to
 * a weak_ptr as soon as the system classes hold their matrices as shared_ptr.
 */
class MultigridLevel
{
    friend class Multigrid;
  public:
    /* ************* *
     * Special member functions. Will perform shallow copies.
     * ************* */
    //! Copy constructor.
    MultigridLevel( const MultigridLevel& ) = default;

    //! Move constructor.
    MultigridLevel( MultigridLevel&& ) = default;

    //! Copy assignment operator.
    MultigridLevel& operator=( const MultigridLevel& ) = default;

    //! Move assignment operator.
    MultigridLevel& operator=( MultigridLevel&& ) = default;

    ~MultigridLevel() = default;

  private:
    /// Constructor which takes a pointer to a BlockFEMatrix and a SmootherCode,
    /// that determines which smoother will be applied on this here level.
    MultigridLevel(BlockFEMatrix* matrix, SmootherCode sm, const ParameterDatabase& db);

    /// Call the stored smoother ones - this should mean one smoothing step.
    void apply_smoother();

    /// Ask the level to compute and store its current defect and residual.
    void calculate_defect();

    /// Ask the level to update its smoother.
    /// This should be called whenever the matrix pointed to by "matrix_"
    /// has changed. The update method of the multigrid object does this for
    /// all its levels.
    void update_smoother();

    BlockFEMatrix* matrix_; // TODO This is a dangerous pointer so far -
                            // TODO change to weak_ptr as soon as the system classes store shared pointers to the matrices

    /// The defect of the equation solved (smoothed) on this level.
    /// Gets updated by calling calculate_defect.
    BlockVector defect_;

    /// The defect of the equation solved (or better: smoothed) on this level.
    /// Gets updated by calling calculate_defect.
    double residual_;

    /// The right hand side of the equation to solve (or rather: smooth) on this level.
    BlockVector rhs_;

    /// The solution of the equation to solve (or rather: smooth) on this level.
    BlockVector solution_;

    /// The smoother object. Since Smoother is an abstract base class,
    /// we can only manage it as a pointer here.
    std::shared_ptr<Smoother> smoother_;

};



#endif /* INCLUDE_MULTIGRID_MULTIGRIDLEVEL_H_ */
