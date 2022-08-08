#include <Preconditioner_iterative.h>
#include <BlockMatrix.h>
#include <BlockFEMatrix.h>
#include <CompositeOperator.h>
#include <BlockVector.h>
#include <memory>

#ifdef _MPI
#include <ParFECommunicator3D.h>
#include <mpi.h>
#endif

template <class L, class V>
Preconditioner_iterative<L, V>::Preconditioner_iterative(
  const L &A,
  std::shared_ptr<IterativeMethod<L, V>> solver)
 : Preconditioner<V>(),
 A(A), solver(solver)
{
}

template <class L, class V>
void Preconditioner_iterative<L, V>::apply(const V &z, V &r) const
{
  r = z;
  solver->iterate(A, z, r);
}

template <class L, class V>
void Preconditioner_iterative<L, V>::update()
{
  solver->update(A);
}

/* ************************************************************************** */
// explicit instantiations
template class Preconditioner_iterative<BlockFEMatrix, BlockVector>;
template class Preconditioner_iterative<CompositeOperator<BlockVector>, BlockVector>;

// In MPI case we are so dependent on the connection of Matrix and FESpace, that
// it does not make sense to instantiate the function for BlockMatrix.
#ifndef _MPI
template class Preconditioner_iterative<BlockMatrix,   BlockVector>;
#endif