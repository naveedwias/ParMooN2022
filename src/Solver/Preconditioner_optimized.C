#include <Preconditioner_optimized.h>
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
Preconditioner_optimized<L, V>::Preconditioner_optimized(
  const L &A,
  std::shared_ptr<Preconditioner<V>> prec,
  double min_d, double max_d)
 : Preconditioner<V>(),
 A(A), prec(prec), min_d(min_d), max_d(max_d)
{
}

template <class L, class V>
void Preconditioner_optimized<L, V>::apply(const V &z, V &r) const
{
  r = z;

#ifdef _MPI
  std::vector<const TParFECommunicator3D*> comms = A.get_communicators();

  for (size_t bl = 0; bl < comms.size(); ++bl)
  {
    comms[bl]->queue_consistency_update(r.block(bl), 1);
  }
  TParFECommunicator3D::flush_consistency_updates();
#endif

  V old_r(r);

  prec->apply(z, r);

  double res_0, res_1, res_2;

#ifdef _MPI
  for (size_t bl = 0; bl < comms.size(); ++bl)
  {
    comms[bl]->queue_consistency_update(r.block(bl), 1);
  }
  TParFECommunicator3D::flush_consistency_updates();
#endif

  // resid = z - A(old_r)
  V resid(z);
  A.apply_scaled_add(old_r, resid, -1.0);
#ifdef _MPI
  res_1 = resid.norm(comms);
#else
  res_1 = resid.norm();
#endif

  // resid = z - Ar
  resid = z;
  A.apply_scaled_add(r, resid, -1.0);
#ifdef _MPI
  res_2 = resid.norm(comms);
#else
  res_2 = resid.norm();
#endif

  // resid = z - A(old_r - (r - old_r))
  //       = z - 2A(old_r) + Ar
  resid = z;
  A.apply_scaled_add(old_r, resid, -2.0);
  A.apply_scaled_add(r, resid, 1.0);
#ifdef _MPI
  res_0 = resid.norm(comms);
#else
  res_0 = resid.norm();
#endif

  res_0 *= res_0;
  res_1 *= res_1;
  res_2 *= res_2;

  // res = at² + bt + c
  // res[0] = a - b + c
  // res[1] = c
  // res[2] = a + b + c
  // a + b = res[2] - res[1]
  // a - b = res[0] - res[1]
  // 2a = res[0] - 2 res[1] + res[2]

  double a = 0.5 * (res_0 + res_2) - res_1;
  double b = res_2 - res_1 - a;
  // double c = res_1;

  if (a > 0.0)
  {
    double d = -0.5 * b / a;

    d = std::max(d, min_d);
    d = std::min(d, max_d);

    if (d != 1.0)
    {
      // r <- (1 - d) old_r + d r
      r.scale(d);
      r.add_scaled(old_r, 1.0 - d);
    }

    /*Output::root_info("Preconditioner_optimized", "Optimized residual: ",
      std::sqrt(a * d * d + b * d + c),
      " / ", std::sqrt(res_2), " / ", std::sqrt(res_1));*/
  }
}

template <class L, class V>
void Preconditioner_optimized<L, V>::update()
{
  prec->update();
}

/* ************************************************************************** */
// explicit instantiations
template class Preconditioner_optimized<BlockFEMatrix, BlockVector>;
template class Preconditioner_optimized<CompositeOperator<BlockVector>, BlockVector>;

// In MPI case we are so dependent on the connection of Matrix and FESpace, that
// it does not make sense to instantiate the function for BlockMatrix.
#ifndef _MPI
template class Preconditioner_optimized<BlockMatrix,   BlockVector>;
#endif