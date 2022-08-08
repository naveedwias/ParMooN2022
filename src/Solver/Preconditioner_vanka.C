/*
 * Preconditioner_vanka.C
 *
 *  Created on: Aug 23, 2016
 *      Author: bartsch
 */

#include <BlockVector.h>
#include <BlockFEMatrix.h>
#include <MooNMD_Io.h>
#include <Preconditioner_vanka.h>

template <class Vector>
Preconditioner_vanka<Vector>::Preconditioner_vanka(const BlockFEMatrix& matrix,
                         VankaType type, double damp_factor)
                         : vanka_object_(type, damp_factor, false),
                           matrix_(&matrix)
{
   // It seems like Solver.C requests a preconditioner object to
   // be fully functional directly after construction, i.e. it has to call
    // its own update method here.
   this->update();
}

// Constructor for compatibility with solver class. Not to be used.
template <class Vector>
Preconditioner_vanka<Vector>::Preconditioner_vanka(const BlockMatrix&,
                                                   VankaType type,
                                                   double damp_factor)
 : vanka_object_(type, damp_factor, false), matrix_(nullptr)
{
  ErrThrow("Creating a Preconditioner_vanka with a BlockMatrix is not "
           "possible, you need a BlockFEMatrix.");
  // And the reason is: The Vanka's need for FE spaces is inherent.
}

template <class Vector>
Preconditioner_vanka<Vector>::Preconditioner_vanka(
  const CompositeOperator<BlockVector>&, VankaType type,
  double damp_factor)
: vanka_object_(type, damp_factor, false), matrix_(nullptr)
{
  ErrThrow("Cannot use Preconditioner_vanka for a composite operator!");
}

template <class Vector>
void Preconditioner_vanka<Vector>::apply(const Vector & z, Vector & r) const
{
    // One complete smoothing sweep - that's it!
    vanka_object_.smooth(z, r);
}

template <class Vector>
void Preconditioner_vanka<Vector>::update()
{
  if(matrix_ != nullptr)
  {
    // Let the Vanka object update itself with the stored matrix.
    vanka_object_.update(*matrix_);
  }
}

// explicit instantiation
template class Preconditioner_vanka<BlockVector>;
