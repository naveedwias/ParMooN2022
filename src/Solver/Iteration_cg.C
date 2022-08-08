#include <Iteration_cg.h>
#include <BlockFEMatrix.h>
#include <CompositeOperator.h>
#include <BlockVector.h>
#include <MooNMD_Io.h>

#ifdef _MPI
#include <ParFECommunicator3D.h>
#include <mpi.h>
#endif

// L - LinearOperator, V - Vector
template <class L, class V>
Iteration_cg<L, V>::Iteration_cg(std::shared_ptr<Preconditioner<V>> prec)
 : IterativeMethod<L, V>(prec, "cg")
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
std::pair<unsigned int, double> Iteration_cg<L, Vector>::iterate(
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

  Vector p, z;
  Vector q(rhs); // this does a copy. We don't need the entries though.
  double alpha, beta, rho;
  double rho_1 = 0.0; // pour le compilateur.

  //Vector r = rhs - A*solution;
  Vector r(rhs); // copy values
  A.apply_scaled_add(solution, r, -1.0);

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

  for (unsigned int i = 1; i <= this->max_n_iterations; ++i) 
  {
#ifdef _MPI
    for (size_t bl = 0; bl < comms.size(); ++bl)
    {
      comms[bl]->queue_consistency_update(r.block(bl), 3);
    }
    TParFECommunicator3D::flush_consistency_updates();
#endif

    //z = M.solve(r);
    this->prec->apply(r, z);

#ifdef _MPI
    rho = dot(r, z, comms);
#else
    rho = dot(r, z);
#endif

    if (i == 1)
    {
      p = z;
    }
    else
    {
      beta = rho / rho_1;
      //p = z + beta * p;
      p *= beta;
      p += z;
    }

#ifdef _MPI
    for (size_t bl = 0; bl < comms.size(); ++bl)
    {
      comms[bl]->queue_consistency_update(p.block(bl), 3);
    }
    TParFECommunicator3D::flush_consistency_updates();
#endif

    //q = A*p;
    A.apply(p, q);

#ifdef _MPI
    alpha = rho / dot(p, q, comms);
#else
    alpha = rho / dot(p, q);
#endif

    //solution += alpha * p;
    solution.add_scaled(p, alpha);

    //r -= alpha * q;
    r.add_scaled(q, -alpha);

#ifdef _MPI
    resid = r.norm(comms) / normb;
#else
    resid = r.norm() / normb;
#endif

    if (this->converged(resid, i))
    {
      return std::pair<unsigned int, double>(i, resid);
    }
    rho_1 = rho;
  }

  // did not converge
  return std::pair<unsigned int, double>(this->max_n_iterations, resid);
}

/* ************************************************************************** */
// explicit instantiations

#ifndef _MPI
template class Iteration_cg<BlockMatrix,   BlockVector>;
#endif

template class Iteration_cg<BlockFEMatrix, BlockVector>;
template class Iteration_cg<CompositeOperator<BlockVector>, BlockVector>;