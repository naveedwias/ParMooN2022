#include <Iteration_richardson.h>
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
Iteration_richardson<L, V>::Iteration_richardson(
  std::shared_ptr<Preconditioner<V>> prec)
 : IterativeMethod<L, V>(prec, "Richardson")
{
  if(!prec)
  {
    ErrThrow("No preconditioner specified. Choose NoPreconditioner<Vector> if "
             "you don't need one");
  }
}

/* ************************************************************************** */
// L - LinearOperator, V - Vector
template <class L, class Vector>
std::pair<unsigned int, double> Iteration_richardson<L, Vector>::iterate(
  const L & A, const Vector & rhs, Vector & solution)
{

  //MPI: rhs and solution in consistency level 0 for computation of global norm

  Vector z(rhs);

#ifdef _MPI
  int size, my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  std::vector<const TParFECommunicator3D*> comms = A.get_communicators();
#endif

#ifndef _MPI
  double normb = rhs.norm();
#elif _MPI
  double normb = rhs.norm(comms);
#endif

  if (normb == 0.0)
  {
    // rhs is zero

    normb = 1.0;
    solution = 0.0;
  }

#ifdef _MPI
  //MPI: solution in consistency level 3 for computation of global norm
  for (size_t bl = 0; bl < comms.size(); ++bl)
  {
    comms[bl]->queue_consistency_update(solution.block(bl), 3);
  }
  TParFECommunicator3D::flush_consistency_updates();
#endif

  //Vector r = rhs - A*solution;
  Vector r(rhs); // copy values
  A.apply_scaled_add(solution, r, -1.0);

#ifdef _MPI
  double resid = r.norm(comms) / normb;
#else
  double resid = r.norm() / normb;
#endif

  // save initial residual, used to check stopping criteria later
  if (this->converged(resid, 0))
  {
    return std::pair<unsigned int, double>(0, resid);
  }

  for (unsigned int i = 1; i <= this->max_n_iterations; i++)
  {
    //z = M.solve(r);
#ifdef _MPI
    //MPI: r in consistency level 1 (which the preconditioner can't do on its own)
    // TODO For a Vanka which contains toxic systems, we will need level 3 consistency here!
    for (size_t bl = 0; bl < comms.size(); ++bl)
    {
      comms[bl]->queue_consistency_update(r.block(bl), 1);
    }
    TParFECommunicator3D::flush_consistency_updates();
#endif

    this->prec->apply(r, z);

    // update solution by w*z (w = damping factor)
    solution.add_scaled(z, this->damping);

    //r = b - A * x;
    r = rhs;

#ifdef _MPI
    //MPI: solution in consistency level 3 for computation of global norm
    for (size_t bl = 0; bl < comms.size(); ++bl)
    {
      comms[bl]->queue_consistency_update(solution.block(bl), 3);
    }
    TParFECommunicator3D::flush_consistency_updates();
#endif

    A.apply_scaled_add(solution, r, -1.0);

#ifdef _MPI
    resid = r.norm(comms)/ normb;
#else
    resid = r.norm() / normb;
#endif

    if (this->converged(resid, i))
    {
      return std::pair<unsigned int, double>(i, resid);
    }
  }

  // did not converge
  return std::pair<unsigned int, double>(this->max_n_iterations, resid);
}

/* ************************************************************************** */
// explicit instantiations
template class Iteration_richardson<BlockFEMatrix, BlockVector>;
template class Iteration_richardson<CompositeOperator<BlockVector>, BlockVector>;

// In MPI case we are so dependent on the connection of Matrix and FESpace, that
// it does not make sense to instantiate the function for BlockMatrix.
#ifndef _MPI
template class Iteration_richardson<BlockMatrix,   BlockVector>;
#endif
