/**
 * @file Implementation of multigrid iteration.
 *
 * @date 2016/05/10
 * @author Clemens Bartsch
 */

#include <BlockFEMatrix.h>
#include <CompositeOperator.h>
#include <BlockVector.h>
#include <Iteration_multigrid.h>


template <class L, class V>
Iteration_multigrid<L,V>::Iteration_multigrid(std::shared_ptr<Multigrid> mg)
:  IterativeMethod<L, V>(std::make_shared<NoPreconditioner<V>>(), "Multigrid"),
   Preconditioner<V>()
{
  mg_ = mg;
}

template <class L, class V>
void Iteration_multigrid<L,V>::apply(const V & z, V & r) const
{
  //Set right hand side and solution on finest level.
  mg_->set_finest_rhs(z);
  mg_->set_finest_sol(r);

  //Apply one multigrid cycle.
  mg_->cycle();

  // Copy solution into solution output.
  r = mg_->get_finest_sol();
}

template <class L, class V>
std::pair<unsigned int, double> Iteration_multigrid<L,V>::iterate(
    const L &, const V &, V &)
{
  ErrThrow("Using Multigrid as an iterative method is a TODO");
}


/* ************************************************************************** */
// explicit instantiations
template class Iteration_multigrid<BlockFEMatrix, BlockVector>;
template class Iteration_multigrid<CompositeOperator<BlockVector>, BlockVector>;
template class Iteration_multigrid<BlockMatrix, BlockVector>; //must be instantiated, but won't work...
