/**
 * @file Pure virtual base class for multigrid smoothers.
 *
 * @date 2016/05/10
 * @author Clemens Bartsch
 */

#ifndef INCLUDE_MULTIGRID_SMOOTHER_H_
#define INCLUDE_MULTIGRID_SMOOTHER_H_

class BlockVector;
class BlockFEMatrix;

/**
 * Pure virtual base class for multigrid smoothers.
 */
class Smoother
{
  public:
    /**
     * Pure virtual smooth method - applies the smoother once.
     * @param rhs[in] The right hand side of the system to be smoothed.
     * @param solution[in, out]  The solution of the system to be smoothed.
     */
    virtual void smooth(const BlockVector& rhs, BlockVector& solution ) = 0;

    /**
     * Pure virtual update method.
     * Update the smoother. Should happen whenever the matrix of the system
     * to be smoothed has changed.
     *
     * TODO I'm not sure whether this tiny interface will suffice for all smoothers!
     */
    virtual void update(const BlockFEMatrix&) = 0;

    /// Default destructor - class does not manage resources.
    virtual ~Smoother() = default;
};


#endif /* INCLUDE_MULTIGRID_SMOOTHER_H_ */
