#include <Iteration_jacobi.h>
#include <BlockFEMatrix.h>
#include <CompositeOperator.h>
#include <BlockVector.h>
#include <MooNMD_Io.h>
#include <algorithm>

#ifdef _MPI
#include <mpi.h>
#include <ParFECommunicator3D.h>
#endif

template <class L>
void extract_diagonal_entries(std::vector<double> &diag_entries, const L &mat)
{
  diag_entries = mat.get_diagonal();

  // make sure there are no zeros on the diagonal
  if (std::find_if(diag_entries.begin(), diag_entries.end(),
    [] (const double& d)
    {
      return d == 0.0;
    })
    != diag_entries.end())
  {
    ErrThrow("zero entry on the diagonal using Jacobi solver/preconditioner");
  }
}

// L - LinearOperator, V - Vector
template <class L, class V>
Iteration_jacobi<L, V>::Iteration_jacobi(const L& mat)
 : IterativeMethod<L, V>(std::make_shared<NoPreconditioner<V>>(), "Jacobi"),
   Preconditioner<V>(),
   diagonal_entries(mat.get_n_total_rows(), 0.0), linear_operator(mat)
{
  // extract the diagonal entries of the matrix
  extract_diagonal_entries(this->diagonal_entries, mat);
}

/* ************************************************************************** */
// L - LinearOperator, V - Vector
template <class L, class Vector>
std::pair<unsigned int, double> Iteration_jacobi<L, Vector>::iterate(
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

  // x^{k+1} = w D^{-1}(b - Rx^k) + (1-w)x^k
  // D is the diagonal, R is A-D, w is the damping parameter
  // Since we can only multiply with all of A, we do the following:
  // -> x^{k+1} = w D^{-1}(b-Ax^k) + x^k

  // extract the diagonal entries of the matrix
  extract_diagonal_entries(this->diagonal_entries, A);
  size_t length = this->diagonal_entries.size();

  Vector r(rhs);
  A.apply_scaled_add(solution, r, -1.);
  // now: r = rhs - A*solution

  double rnorm, normb;

#ifdef _MPI
  rhs.queue_norm(normb, comms);
  r.queue_norm(rnorm, comms);

  flush_norm_queue();
#else
  rnorm = r.norm();
  normb = rhs.norm();
#endif

  if (normb == 0.0)
  {
    // rhs is zero

    normb = 1.0;
    rnorm = 0.0;
    solution = 0.0;
  }

  // check for convergence
  double resid = rnorm / normb;

  // safe initial residual, used to check stopping criteria later
  if(this->converged(resid, 0)) 
  {
    return std::pair<unsigned int, double>(0, resid);
  }

  for (unsigned int i = 1; i <= this->max_n_iterations; ++i)
  {
    // now: r = rhs - A*solution
    for (size_t j = 0; j < length; ++j)
    {
      r[j] *= this->damping / this->diagonal_entries[j];
    }

    // now: r = w D^{-1}(rhs - A*solution)
    r += solution;

    // now: r = w D^{-1}(rhs - A*solution) + solution
    std::swap(solution, r); // solution is now the new solution

#ifdef _MPI
    for (size_t bl = 0; bl < comms.size(); ++bl)
    {
      comms[bl]->queue_consistency_update(solution.block(bl), 3);
    }
    TParFECommunicator3D::flush_consistency_updates();
#endif

    // compute residual
    r = rhs;
    A.apply_scaled_add(solution, r, -1.);

#ifdef _MPI
    rnorm = r.norm(comms);
#else
    rnorm = r.norm();
#endif

    // check for convergence
    resid = rnorm / normb;

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
void Iteration_jacobi<L, V>::update(const L& A)
{
  if (&this->linear_operator != &A)
  {
    ErrThrow("You are trying to update a Iteration_jacobi object with a matrix "
             "which is not the one you constructed it with. Please create a "
             "new Iteration_jacobi object");
  }

  // update the diagonal entries
  extract_diagonal_entries(this->diagonal_entries, A);
  this->IterativeMethod<L, V>::update(A);
}

/* ************************************************************************** */
// L - LinearOperator, V - Vector
template <class L, class V>
void Iteration_jacobi<L, V>::apply(const V & z, V & r) const
{
  // we don't have to do anything special in the MPI case here - unless the
  // computation is seriously overdistributed anyway, the extra overhead of
  // picking out master DOFs is not worth it

  // as if no preconditioner is used (this initialized r, in case it has not 
  // been done already)
  r = z;

  for (size_t j = 0, length = this->diagonal_entries.size(); j < length; ++j)
  {
    r[j] /= this->diagonal_entries[j];
  }
}

/* ************************************************************************** */
// L - LinearOperator
template <class L, class V>
const L& Iteration_jacobi<L, V>::get_operator() const
{
  return linear_operator;
}

/* ************************************************************************** */
// explicit instantiations
template class Iteration_jacobi<BlockFEMatrix, BlockVector>;
template class Iteration_jacobi<CompositeOperator<BlockVector>, BlockVector>;

#ifndef _MPI
template class Iteration_jacobi<BlockMatrix,   BlockVector>;
#endif