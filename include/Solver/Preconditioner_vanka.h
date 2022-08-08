/*
 * Preconditioner_vanka.h
 *
 *  Created on: Aug 23, 2016
 *      Author: bartsch
 */

#ifndef INCLUDE_SOLVER_PRECONDITIONER_VANKA_H_
#define INCLUDE_SOLVER_PRECONDITIONER_VANKA_H_

#include <Preconditioner.h>
#include <VankaSmoother.h>

//forward decalaration
class BlockFEMatrix;
class BlockVector;

template<class V>
class CompositeOperator;

/**
 * Wrap up a Vanka method as known from multigrid smoothing to be used as a
 * preconditioner. I do not expect this to be a very good preconditioner, but
 * it might be useful for debugging purpose, as it allows one to use Vanka
 * without the multigrid method wrapped around it.
 */
template <class Vector>
class Preconditioner_vanka : public Preconditioner<Vector>
{
  public:
    /// Construct a Vanka preconditioner.
    Preconditioner_vanka(const BlockFEMatrix& matrix,
                         VankaType type, double damp_factor);

    /** @brief don't use this constuctor. It is here only for compatability
     * in the Solver class.
     *
     * @warning Do not use this constructor. You need a BlockFEMatrix instead.
     */
    Preconditioner_vanka(const BlockMatrix & matrix,
                         VankaType type, double damp_factor);

    Preconditioner_vanka(const CompositeOperator<BlockVector> &matrix,
                         VankaType type, double damp_factor);

    // Perform one step of Vanka smoothing as preconditioning.
    virtual void apply(const Vector & z, Vector & r) const override;

    // Update the smoother. Call this whenever the matrix has changed.
    virtual void update() override;

  private:
    /// The Vanka object which does all the work.
    mutable VankaSmoother vanka_object_;
    /// The BlockFEMatrix which this is a preconditioner for.
    const BlockFEMatrix* matrix_;
};


#endif /* INCLUDE_SOLVER_PRECONDITIONER_VANKA_H_ */
