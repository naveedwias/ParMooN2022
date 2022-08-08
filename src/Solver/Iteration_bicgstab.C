#include <Iteration_bicgstab.h>
#include <BlockFEMatrix.h>
#include <CompositeOperator.h>
#include <BlockVector.h>
#include <MooNMD_Io.h>
#include <typeinfo>

#ifdef _MPI
#include <ParFECommunicator3D.h>
#include <mpi.h>
#endif

// L - LinearOperator, V - Vector
template <class L, class V>
Iteration_bicgstab<L, V>::Iteration_bicgstab(
  std::shared_ptr<Preconditioner<V>> prec)
 : IterativeMethod<L, V>(prec, "bi-cgstab")
{
  if (!prec)
  {
    ErrThrow("No preconditioner specified. Choose NoPreconditioner<Vector> if "
             "you don't need one");
  }
}

/* ************************************************************************** */
// L - LinearOperator, V - Vector
template <class L, class Vector>
std::pair<unsigned int, double> Iteration_bicgstab<L, Vector>::iterate(
  const L & A, const Vector & rhs, Vector & solution)
{
  // MPI: rhs and solution in consistency level 0 for computation of global norm
  // MPI Environment initialization for parallel iterative solver
#ifdef _MPI
  int size, my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  std::vector<const TParFECommunicator3D*> comms = A.get_communicators();
#endif

  double rho_1, rho_2 = 1.0, alpha = 1.0, beta, omega = 1.0;
  Vector p, phat(rhs), s, shat(rhs);
  Vector t(rhs), v(rhs); // this does a copy. We don't need the entries though.

#ifdef _MPI
  // MPI: solution in consistency level 3 for matrix.vector multiplication "A.solution"
  for (size_t bl = 0; bl < comms.size(); ++bl)
  {
    comms[bl]->queue_consistency_update(solution.block(bl), 3);
  }
  TParFECommunicator3D::flush_consistency_updates();
#endif

  // Vector r = rhs - A * solution;
  Vector r(rhs); // copy values

  A.apply_scaled_add(solution, r, -1.0); // now r = rhs - A*solution

  Vector rtilde(r);

  // Computation of norms (rhs and r in consistency level 0)

double normb, rnorm;

#ifdef _MPI
  rhs.queue_norm(normb, comms);
  r.queue_norm(rnorm, comms);

  flush_norm_queue();
#else
  normb = rhs.norm();
  rnorm = r.norm();
#endif

  if (normb == 0.0)
  {
    // rhs is zero

    normb = 1.0;
    rnorm = 0.0;
    solution = 0.0;
  }

  double resid = rnorm / normb;

  // save initial residual, used to check stopping criteria later
  if (this->converged(resid, 0))
  {
    return std::pair<unsigned int, double>(0, resid);
  }

  for (unsigned int i = 1; i <= this->max_n_iterations; i++)
  {
    // compute dot-product of 2 vectors (no consistency update needed)
#ifndef _MPI
    rho_1 = dot(rtilde, r);
#elif _MPI
    rho_1 = dot(rtilde, r, comms);
#endif

    if (rho_1 == 0.0) // this should not ever happen
    {
#ifndef _MPI
      double abnormal_residual = r.norm() / normb;
#elif _MPI
      double abnormal_residual = r.norm(comms) / normb;
#endif
      return std::pair<unsigned int, double>(i, abnormal_residual);
    }

    if (i == 1)
    {
      p = r;
    }
    else
    {
      beta = (rho_1 / rho_2) * (alpha / omega);

      // p = r + beta * (p - omega * v);
      p.add_scaled(v, -omega);
      p *= beta;
      p += r;
    }

#ifdef _MPI
    // MPI: p in consistency level 1 (which the preconditioner can't do on its own)
    // TODO For a Vanka which contains toxic systems, we will need level 3 consistency here!
    for (size_t bl = 0; bl < comms.size(); ++bl)
    {
      comms[bl]->queue_consistency_update(p.block(bl), 1);
    }
    TParFECommunicator3D::flush_consistency_updates();
#endif

    // phat = M.solve(p);
    this->prec->apply(p, phat);

#ifdef _MPI
    for (size_t bl = 0; bl < comms.size(); ++bl)
    {
      comms[bl]->queue_consistency_update(phat.block(bl), 3);
    }
    TParFECommunicator3D::flush_consistency_updates();
#endif

    // v = A * phat;
    A.apply(phat, v);

    // compute dot-product of 2 vectors
#ifndef _MPI
    alpha = rho_1 / dot(rtilde, v);
#elif _MPI
    alpha = rho_1 / dot(rtilde, v, comms);
#endif

    // s = r - alpha * v;
    s = r;
    s.add_scaled(v, -alpha);

#ifndef _MPI
    resid = s.norm();
#elif _MPI
    resid = s.norm(comms);
#endif

    if (this->converged(resid, i))
    {
      // solution += alpha * phat;
      solution.add_scaled(phat, alpha);
      return std::pair<unsigned int, double>(i, resid);
    }

#ifdef _MPI
    // MPI: s in consistency level 1 (which the preconditioner can't do on its own)
    // TODO For a Vanka which contains toxic systems, we will need level 3 consistency here!
    for (size_t bl = 0; bl < comms.size(); ++bl)
    {
      comms[bl]->queue_consistency_update(s.block(bl), 1);
    }
    TParFECommunicator3D::flush_consistency_updates();
#endif

    // shat = M.solve(s);
    this->prec->apply(s, shat);

#ifdef _MPI
    for (size_t bl = 0; bl < comms.size(); ++bl)
    {
      comms[bl]->queue_consistency_update(shat.block(bl), 3);
    }
    TParFECommunicator3D::flush_consistency_updates();
#endif

    // t = A * shat;
    A.apply(shat, t);

#ifndef _MPI
    omega = dot(t,s) / dot(t,t);
#elif _MPI
    double dot_ts, dot_tt;

    queue_dot_global(t, s, comms, dot_ts);
    queue_dot_global(t, t, comms, dot_tt);
    flush_dot_queue();

    omega = dot_ts / dot_tt;

#endif

    // solution += alpha * phat + omega * shat;
    solution.add_scaled(phat, alpha);
    solution.add_scaled(shat, omega);

    // r = s - omega * t;
    r = s;
    r.add_scaled(t, -omega);

    rho_2 = rho_1;

#ifndef _MPI
    resid = r.norm();
#elif _MPI
    resid = r.norm(comms);
#endif

    if (this->converged(resid, i))
    {
      return std::pair<unsigned int, double>(i, resid);
    }

    if (omega == 0)
    {
      ErrThrow("unhappy breakdown in Bi-CGStab. iteration ", i, ", residual ",
               resid);
    }
  }

  // did not converge
  return std::pair<unsigned int, double>(this->max_n_iterations, resid);
}

/* ************************************************************************** */
// explicit instantiations
template class Iteration_bicgstab<BlockFEMatrix, BlockVector>;
template class Iteration_bicgstab<CompositeOperator<BlockVector>, BlockVector>;

// In MPI case we are so dependent on the connection of Matrix and FESpace, that
// it does not make sense to instantiate the function for BlockMatrix.
#ifndef _MPI
template class Iteration_bicgstab<BlockMatrix,   BlockVector>;
#endif