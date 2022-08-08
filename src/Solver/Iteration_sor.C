#include <Iteration_sor.h>
#include <BlockFEMatrix.h>
#include <CompositeOperator.h>
#include <BlockVector.h>
#include <MooNMD_Io.h>
#include <algorithm>

#ifdef _MPI
#include <ParFECommunicator3D.h>
#include <mpi.h>
#endif

// L - LinearOperator, V - Vector
template <class L, class V>
Iteration_sor<L, V>::Iteration_sor(const L& mat, int flag, double w)
  : IterativeMethod<L, V>(std::make_shared<NoPreconditioner<V>>(),
    flag == 2 ? "ssor" : "sor"),
    Preconditioner<V>(),
    sor_type(flag), omega(w), linear_operator(mat)
{
#ifdef _MPI
  parallel_strategy_ = std::string("own_cells"); //default parallel strategy
#endif
}

/* ************************************************************************** */
// L - LinearOperator, V - Vector
template <class L, class Vector>
std::pair<unsigned int, double> Iteration_sor<L, Vector>::iterate(
  const L & A, const Vector & rhs, Vector & solution)
{
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

  // (D + wL) x^{k+1} = w b - (wU + (w-1)D)x^k
  // D is the diagonal, L is the lower triangular part of A, U is the upper 
  // triangular part of A, w is the relaxation parameter
  // Since we can only multiply with all of A, we do the following:
  // -> x^{k+1} = w (D+wL)^{-1}(b-Ax^k) + x^k

  Vector r(rhs);
  A.apply_scaled_add(solution, r, -1.0);
  // now: r = rhs - A*solution

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

  // check for convergence
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

  for (unsigned int i = 1; i <= this->max_n_iterations; ++i)
  {
    // one (s)sor step
#ifdef _MPI
    A.sor_sweep(rhs, solution, omega, sor_type, parallel_strategy_);

    for (size_t bl = 0; bl < comms.size(); ++bl)
    {
      comms[bl]->queue_consistency_update(solution.block(bl), 3);
    }
    TParFECommunicator3D::flush_consistency_updates();
#else
    A.sor_sweep(rhs, solution, this->omega, this->sor_type);
#endif

    //solution.print(std::string("sol_") + std::to_string(i));
    // compute new residual
    r = rhs;
    A.apply_scaled_add(solution, r, -1.);

    // check for convergence

#ifdef _MPI
    resid = r.norm(comms) / normb;
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
// L - LinearOperator, V - Vector
template <class L, class V>
void Iteration_sor<L, V>::apply(const V & z, V & r) const
{
  // initialize r, in case it has not been done already
  r.copy_structure(z);
  r.reset();

#ifdef _MPI
  linear_operator.sor_sweep(z, r, omega, sor_type, parallel_strategy_);
#else
  this->linear_operator.sor_sweep(z, r, this->omega, this->sor_type);
#endif
}

#ifdef _MPI
/* ************************************************************************** */
// L - LinearOperator, V - Vector
template <class L, class V>
void Iteration_sor<L, V>::set_parallel_strategy(const std::string& parallel_strategy)
{
  parallel_strategy_= parallel_strategy;
}
#endif

/* ************************************************************************** */
// L - LinearOperator, V - Vector
template <class L, class V>
void Iteration_sor<L, V>::apply_smoother(const V & z, V & r) const
{ // (not setting the start iterate zero!)
#ifdef _MPI
  linear_operator.sor_sweep(z, r, omega, sor_type, parallel_strategy_);
#else
  this->linear_operator.sor_sweep(z, r, this->omega, this->sor_type);
#endif
}

/* ************************************************************************** */
// L - LinearOperator
template <class L, class V>
const L& Iteration_sor<L, V>::get_operator() const
{
  return linear_operator;
}

/* ************************************************************************** */
// explicit instantiations
template class Iteration_sor<BlockFEMatrix, BlockVector>;
template class Iteration_sor<CompositeOperator<BlockVector>, BlockVector>;

#ifndef _MPI
template class Iteration_sor<BlockMatrix,   BlockVector>;
#endif