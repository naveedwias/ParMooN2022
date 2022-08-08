#include <Iteration_cgs.h>
#include <BlockFEMatrix.h>
#include <CompositeOperator.h>
#include <BlockVector.h>
#include <Iteration_jacobi.h>
#include <MooNMD_Io.h>

#ifdef _MPI
#include <ParFECommunicator3D.h>
#include <mpi.h>
#endif

// L - LinearOperator, V - Vector
template <class L, class V>
Iteration_cgs<L, V>::Iteration_cgs(std::shared_ptr<Preconditioner<V>> prec)
 : IterativeMethod<L, V>(prec, "cgs")
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
std::pair<unsigned int, double> Iteration_cgs<L, Vector>::iterate(
  const L & A, const Vector & rhs, Vector & solution)
{
#ifdef _MPI
  int size, my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  std::vector<const TParFECommunicator3D*> comms = A.get_communicators();

  for (size_t bl = 0; bl < comms.size(); ++bl)
  {
    comms[bl]->queue_consistency_update(solution.block(bl), 3);
  }
  TParFECommunicator3D::flush_consistency_updates();
#endif

  double alpha, beta, rho_1, rho_2;
  Vector p, phat(rhs), q, u, uhat(rhs);
  Vector qhat(rhs), vhat(rhs); // this does a copy. We don't need the entries.

  //Vector r = rhs - A*solution;
  Vector r(rhs); // copy values
  A.apply_scaled_add(solution, r, -1.0); // now r = rhs - A*solution

  //Vector rtilde = r;
  Vector rtilde(r); // copy values

#ifdef _MPI
  double normb = rhs.norm(comms);
#else
  double normb = rhs.norm();
#endif

  if (normb == 0.0)
  {
    // rhs is zero
    normb = 1.0;
    r = 0.0;
    solution = 0.0;
  }

#ifdef _MPI
  double resid = r.norm(comms) / normb;
#else
  double resid = r.norm() / normb;
#endif

  // save initial residual, used to check stopping criteria later
  if(this->converged(resid, 0)) 
  {
    return std::pair<unsigned int, double>(0, resid);
  }

  for (unsigned int i = 1; i <= this->max_n_iterations; i++)
  {
#ifdef _MPI
    rho_1 = dot(rtilde, r, comms);
#else
    rho_1 = dot(rtilde, r);
#endif

    if (rho_1 == 0) // this should not ever happen 
    {
      return std::pair<unsigned int, double>(i, r.norm() / normb);
    }

    if (i == 1)
    {
      u = r;
      p = u;
    }
    else
    {
      beta = rho_1 / rho_2;

      // u = r + beta * q;
      u = q;
      u *= beta;
      u += r;

      // p = u + beta * (q + beta * p);
      p *= beta;
      p += q;
      p *= beta;
      p += u;
    }

#ifdef _MPI
    for (size_t bl = 0; bl < comms.size(); ++bl)
    {
      comms[bl]->queue_consistency_update(p.block(bl), 3);
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

    //vhat = A*phat;
    A.apply(phat, vhat);

#ifdef _MPI
    alpha = rho_1 / dot(rtilde, vhat, comms);
#else
    alpha = rho_1 / dot(rtilde, vhat);
#endif

    //q = u - alpha * vhat;
    q = u;
    q.add_scaled(vhat, -alpha);

    //uhat = M.solve(u + q);
    u += q;

#ifdef _MPI
    for (size_t bl = 0; bl < comms.size(); ++bl)
    {
      comms[bl]->queue_consistency_update(u.block(bl), 3);
    }
    TParFECommunicator3D::flush_consistency_updates();
#endif

    this->prec->apply(u, uhat);

#ifdef _MPI
    for (size_t bl = 0; bl < comms.size(); ++bl)
    {
      comms[bl]->queue_consistency_update(uhat.block(bl), 3);
    }
    TParFECommunicator3D::flush_consistency_updates();
#endif

    //solution += alpha * uhat;
    solution.add_scaled(uhat, alpha);

    //qhat = A * uhat;
    A.apply(uhat, qhat);

    //r -= alpha * qhat;
    qhat *= -alpha;
    r += qhat;
    rho_2 = rho_1;

#ifdef _MPI
    resid = r.norm(comms) / normb;
#else
    resid = r.norm() / normb;
#endif

    if(this->converged(resid, i))
    {
      return std::pair<unsigned int, double>(i, resid);
    }
  }

  // did not converge
  return std::pair<unsigned int, double>(this->max_n_iterations, resid);
}

/* ************************************************************************** */
// explicit instantiations

#ifndef _MPI
template class Iteration_cgs<BlockMatrix,   BlockVector>;
#endif

template class Iteration_cgs<CompositeOperator<BlockVector>, BlockVector>;
template class Iteration_cgs<BlockFEMatrix, BlockVector>;