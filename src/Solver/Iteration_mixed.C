#include <Iteration_mixed.h>
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
Iteration_mixed<L, V>::Iteration_mixed(
  std::shared_ptr<Preconditioner<V>> prec)
  : IterativeMethod<L, V>(prec, "mixed")
{
  if (!prec)
  {
    ErrThrow("No preconditioner specified. Choose NoPreconditioner<Vector> if "
             "you don't need one");
  }

  gmres_iterator = std::make_shared<Iteration_gmres<L, V>>(prec,
    gmres_type::flexible);
  bicgstab_iterator = std::make_shared<Iteration_bicgstab<L, V>>(prec);
}

/* ************************************************************************** */
// L - LinearOperator, V - Vector
template <class L, class Vector>
std::pair<unsigned int, double> Iteration_mixed<L, Vector>::iterate(
  const L &A, const Vector &rhs, Vector &solution)
{
  // MPI: rhs and solution in consistency level 0 for computation of global norm
  // MPI Environment initialization for parallel iterative solver
#ifdef _MPI
  int size, my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  std::vector<const TParFECommunicator3D*> comms = A.get_communicators();

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

  double rhs_norm, r_norm;

#ifdef _MPI
  rhs.queue_norm(rhs_norm, comms);
  r.queue_norm(r_norm, comms);

  flush_norm_queue();
#else
  rhs_norm = rhs.norm();
  r_norm = r.norm();
#endif

  if (rhs_norm == 0.0)
  {
    // rhs is zero

    solution = 0.0;
    rhs_norm = 1.0;
    r_norm = 0.0;
  }

  double initial_residual = r_norm / rhs_norm;

  if (this->converged(initial_residual, 0))
  {
    return std::make_pair(0u, initial_residual);
  }

  double gmres_res_tol = std::max(initial_residual * this->residual_reduction,
    this->residual_tolerance) * rhs_norm;
  double bicgstab_res_tol = std::max(initial_residual * this->residual_reduction,
    this->residual_tolerance);

  double residual = initial_residual;

  for (unsigned int i = 1; i <= this->max_n_iterations; )
  {
    int gmres_max_it = this->restart;
    if (i + this->restart > this->max_n_iterations)
    {
      gmres_max_it = this->max_n_iterations - i;
    }

    gmres_iterator->set_stopping_parameters(gmres_max_it, 0,
      gmres_res_tol, 0.0, this->divergence_factor,
      this->damping, this->restart);

    auto gmres_result = gmres_iterator->iterate(A, rhs, solution);

    i += gmres_result.first;
    residual = gmres_result.second / rhs_norm;

    if (this->converged(residual, i))
    {
      return std::make_pair(i, residual);
    }

    int bicgstab_max_it = this->restart;
    if (i + this->restart > this->max_n_iterations)
    {
      bicgstab_max_it = this->max_n_iterations - i;
    }

    bicgstab_iterator->set_stopping_parameters(bicgstab_max_it, 0,
      bicgstab_res_tol, 0.0, this->divergence_factor,
      this->damping, this->restart);

    auto bicgstab_result = bicgstab_iterator->iterate(A, rhs, solution);

    i += bicgstab_result.first;
    residual = bicgstab_result.second;

    if (this->converged(residual, i))
    {
      return std::make_pair(i, residual);
    }
  }

  // did not converge
  return std::make_pair(this->max_n_iterations, residual);
}

/* ************************************************************************** */
// explicit instantiations
template class Iteration_mixed<BlockFEMatrix, BlockVector>;
template class Iteration_mixed<CompositeOperator<BlockVector>, BlockVector>;

// In MPI case we are so dependent on the connection of Matrix and FESpace, that
// it does not make sense to instantiate the function for BlockMatrix.
#ifndef _MPI
template class Iteration_mixed<BlockMatrix, BlockVector>;
#endif