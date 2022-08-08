#include <Iteration_gmres.h>
#include <BlockFEMatrix.h>
#include <CompositeOperator.h>
#include <BlockVector.h>
#include <cmath> // std::abs
#include <MooNMD_Io.h>

#ifdef _MPI
#include <ParFECommunicator3D.h>
#include <mpi.h>
#include <vector>
#endif

#include <Utilities_gmres.h>

std::string type_name(gmres_type t)
{
  switch (t)
  {
    case gmres_type::left:
      return "left gmres";
      break;
    case gmres_type::right:
      return "right gmres";
      break;
    case gmres_type::flexible:
      return "flexible gmres";
      break;
    default:
      ErrThrow("unknown gmres type ");
      break;
  }
}

template <class LinearOperator, class Vector>
Iteration_gmres<LinearOperator, Vector>::Iteration_gmres(
  std::shared_ptr<Preconditioner<Vector>> prec, gmres_type t)
 : IterativeMethod<LinearOperator, Vector>(prec, type_name(t)), type(t), s(),
   cs(), sn(), v(), z()
{
  if (!prec)
  {
    ErrThrow("No preconditioner specified. Choose NoPreconditioner<Vector> if "
             "you don't need one");
  }
}

/* ************************************************************************** */
template <class LinearOperator, class Vector>
std::pair<unsigned int, double>
Iteration_gmres<LinearOperator, Vector>::iterate(const LinearOperator & A,
                                                 const Vector & rhs,
                                                 Vector & solution)
{
  switch (this->type)
  {
    case gmres_type::left:
      return left_gmres(A, rhs, solution);
      break;
    case gmres_type::right:
      return right_gmres(A, rhs, solution);
      break;
    case gmres_type::flexible:
      return flexible_gmres(A, rhs, solution);
      break;
    default:
      ErrThrow("unknown gmres type ");
      break;
  }
}

/* ************************************************************************** */
template <class LinearOperator, class Vector>
std::pair<unsigned int, double>
Iteration_gmres<LinearOperator, Vector>::left_gmres(const LinearOperator & A,
                                                    const Vector & rhs,
                                                    Vector & solution)
{
// MPI: rhs and solution in consistency level 0 for computation of global norm
#ifdef _MPI
  int size, my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  std::vector<const TParFECommunicator3D*> comms = A.get_communicators();
#endif

  this->s.resize(this->restart + 1);
  this->cs.resize(this->restart + 1);
  this->sn.resize(this->restart + 1);
  Vector w;

  HessenbergMatrix H(this->restart + 1);

  Vector r;
  Vector a = rhs; // intermediate vector

#ifdef _MPI
  // MPI: solution in consistency level 3 for vector.matrix multiplication
  for (size_t bl = 0; bl < comms.size(); ++bl)
  {
    comms[bl]->queue_consistency_update(solution.block(bl), 3);
  }

  TParFECommunicator3D::flush_consistency_updates();
#endif

  A.apply_scaled_add(solution, a, -1.0); // now a = rhs - A*solution
  this->prec->apply(a, r);

#ifdef _MPI
  double resid = r.norm(comms); // compute initial residual
#else
  double resid = r.norm();
#endif

  double beta = resid; // initialize beta as initial residual

  // save initial residual, used to check stopping criteria later
  if (this->converged(resid, 0))
  {
    return std::pair<unsigned int, double>(0, resid);
  }

  this->v.resize(this->restart + 1);

  unsigned int j = 1; // iteration index
  while (j <= this->max_n_iterations)
  {
    v[0] = r;
    v[0] *= 1.0 / beta;
    std::fill(s.begin(), s.end(), 0.0); // set all entries in s to zero
    s[0] = beta;

    for (unsigned int i = 0;
      i < this->restart && j <= this->max_n_iterations;
      i++, j++)
    {
#ifdef _MPI
      // MPI: solution in consistency level 3 for computation of global norm
      for (size_t bl = 0; bl < comms.size(); ++bl)
      {
        comms[bl]->queue_consistency_update(v[i].block(bl), 3);
      }

      TParFECommunicator3D::flush_consistency_updates();
#endif

      // w = M.solve(A * v[i]);
      A.apply(v[i], a); // a = A*v[i] // reuse Vector 'a'

      this->prec->apply(a, w);

#ifdef _MPI
      for (unsigned int k = 0; k <= i; k++)
      {
        queue_dot_global(w, v[k], comms, H(k, i));
      }

      flush_dot_queue();
#else
      for (unsigned int k = 0; k <= i; k++)
      {
        H(k, i) = dot(w, v[k]);
      }
#endif

      for (unsigned int k = 0; k <= i; k++)
      {
        w.add_scaled(v[k], -H(k, i));
      }

#ifdef _MPI
      H(i + 1, i) = w.norm(comms);
#else
      H(i + 1, i) = w.norm();
#endif

      if (H(i + 1, i) != 0.0)
      {
        v[i + 1] = w;
        v[i + 1] *= 1.0 / H(i + 1, i);
      }
      else
      {
        ErrThrow("Unhappy breakdown in (left) gmres at iteration ", j);
      }

      for (unsigned int k = 0; k < i; k++)
      {
        ApplyGivensRotation(H(k, i), H(k + 1, i), cs[k], sn[k]);
      }

      GenerateGivensRotation(H(i, i), H(i + 1, i), cs[i], sn[i]);
      ApplyGivensRotation(H(i, i), H(i + 1, i), cs[i], sn[i]);
      ApplyGivensRotation(s[i], s[i + 1], cs[i], sn[i]);

      resid = std::abs(s[i + 1]);
      if (this->converged(resid, j))
      {
        Update(solution, i, H, s, v);
        return std::pair<unsigned int, double>(j, resid);
      }
    }

    Update(solution, this->restart - 1, H, s, v);

#ifdef _MPI
    // MPI: solution in consistency level 3 for computation of global norm
    for (size_t bl = 0; bl < comms.size(); ++bl)
    {
      comms[bl]->queue_consistency_update(solution.block(bl), 3);
    }

    TParFECommunicator3D::flush_consistency_updates();
#endif

    // r = M.solve(b - A * x);
    a = rhs;
    A.apply_scaled_add(solution, a, -1.0);
    this->prec->apply(a, r);

#ifdef _MPI
    beta = r.norm(comms);
#else
    beta = r.norm();
#endif

    if (std::abs(beta - resid) > 0.01 * beta)
    {
      Output::root_info<1>("LGMRES", "restart residual changed ", beta, "  ", resid);
    }

    if (this->converged(resid, j))
    {
      return std::pair<unsigned int, double>(j, resid);
    }
    // else restart
  }

  // did not converge
  return std::pair<unsigned int, double>(this->max_n_iterations, resid);
}

/* ************************************************************************** */
template <class LinearOperator, class Vector>
std::pair<unsigned int, double>
Iteration_gmres<LinearOperator, Vector>::right_gmres(const LinearOperator & A,
                                                     const Vector & rhs,
                                                     Vector & solution)
{
  // MPI: rhs and solution in consistency level 0 for computation of global norm
#ifdef _MPI
  int size, my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  std::vector<const TParFECommunicator3D*> comms = A.get_communicators();
#endif

  this->s.resize(this->restart + 1);
  this->cs.resize(this->restart + 1);
  this->sn.resize(this->restart + 1);
  Vector w(rhs);

  HessenbergMatrix H(this->restart + 1);

  // Vector r = b - A * solution;
  Vector r;
  r = rhs; // reuse Vector 'a' to compute residual

#ifdef _MPI
  // MPI: solution in consistency level 3 for vector.matrix multiplication
  for (size_t bl = 0; bl < comms.size(); ++bl)
  {
    comms[bl]->queue_consistency_update(solution.block(bl), 3);
  }

  TParFECommunicator3D::flush_consistency_updates();
#endif

  A.apply_scaled_add(solution, r, -1.0); // now r = rhs - Ax

  // compute initial residual
#ifdef _MPI
  double resid = r.norm(comms);
#else
  double resid = r.norm();
#endif

  double beta = resid; // initialize beta as initial residual

  // save initial residual, used to check stopping criteria later
  if (this->converged(resid, 0))
  {
    return std::pair<unsigned int, double>(0, resid);
  }

  this->v.resize(this->restart + 1);

  unsigned int j = 1; // iteration index
  while (j <= this->max_n_iterations) 
  {
    // v[0] = r * (1.0 / beta);
    v[0] = r;
    v[0] *= 1.0 / beta;
    std::fill(s.begin(), s.end(), 0.0); // set all entries in s to zero
    s[0] = beta;

    for (unsigned int i = 0; i < this->restart && j <= this->max_n_iterations;
        i++, j++)
    {
#ifdef _MPI
      // MPI: solution in consistency level 3 for computation of global norm
      for (size_t bl = 0; bl < comms.size(); ++bl)
      {
        comms[bl]->queue_consistency_update(v[i].block(bl), 3);
      }

      TParFECommunicator3D::flush_consistency_updates();
#endif

      // r = A * M.solve(v[i]);
      w = 0.0;
      this->prec->apply(v[i], w); // w = M.solve(v[i])
      A.apply(w, r); // r = A * w

#ifdef _MPI
      for (unsigned int k = 0; k <= i; k++) 
      {
        queue_dot_global(r, v[k], comms, H(k, i));
      }

      flush_dot_queue();
#else
      for (unsigned int k = 0; k <= i; k++) 
      {
        H(k, i) = dot(r, v[k]);
      }
#endif

      for (unsigned int k = 0; k <= i; k++) 
      {
        r.add_scaled(v[k], -H(k, i));
      }

#ifdef _MPI
      H(i + 1, i) = r.norm(comms);
#else
      H(i + 1, i) = r.norm();
#endif

      if (H(i + 1, i) != 0.0)
      {
        // v[i+1] = r * (1.0 / H(i+1, i));
        v[i + 1] = r;
        v[i + 1] *= 1.0 / H(i + 1, i);
      }
      else
      {
        ErrThrow("Unhappy breakdown in (right) gmres at iteration ", j);
      }

      for (unsigned int k = 0; k < i; k++)
      {
        ApplyGivensRotation(H(k, i), H(k + 1, i), cs[k], sn[k]);
      }

      GenerateGivensRotation(H(i, i), H(i + 1, i), cs[i], sn[i]);
      ApplyGivensRotation(H(i, i), H(i + 1, i), cs[i], sn[i]);
      ApplyGivensRotation(s[i], s[i + 1], cs[i], sn[i]);

      resid = std::abs(s[i + 1]);
      if (this->converged(resid, j)) 
      {
        w = 0.0; // reuse w for update of solution
        Update(w, i, H, s, v);
        this->prec->apply(w, r); // r = M.solve(w);
        solution += r;

        return std::pair<unsigned int, double>(j, resid);
      }
    }

    w = 0.0;
    Update(w, this->restart - 1, H, s, v);

    this->prec->apply(w, r); // r = M.solve(w)

    solution += r; // update solution

    // compute new residual
    r = rhs;

#ifdef _MPI
    // MPI: solution in consistency level 3 for computation of global norm
    for (size_t bl = 0; bl < comms.size(); ++bl)
    {
      comms[bl]->queue_consistency_update(solution.block(bl), 3);
    }

    TParFECommunicator3D::flush_consistency_updates();
#endif

    A.apply_scaled_add(solution, r, -1.0); // r = rhs - A * solution;

#ifdef _MPI
    beta = r.norm(comms);
#else
    beta = r.norm();
#endif

    if (std::abs(beta - resid) > 0.01 * beta)
    {
      Output::root_info<1>("RGMRES", "restart residual changed ", beta, "  ", resid);
    }

    if(this->converged(resid, j))
    {
      return std::pair<unsigned int, double>(j, resid);
    }
    // else restart
  }

  // not converged
  return std::pair<unsigned int, double>(this->max_n_iterations, resid);
}

/* ************************************************************************** */
template <class LinearOperator, class Vector>
std::pair<unsigned int, double>
Iteration_gmres<LinearOperator, Vector>::flexible_gmres(const LinearOperator& A,
                                                        const Vector & rhs,
                                                        Vector & solution)
{
  // MPI: rhs and solution in consistency level 0 for computation of global norm
#ifdef _MPI
  int size, my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  std::vector<const TParFECommunicator3D*> comms = A.get_communicators();
#endif

  this->s.resize(this->restart + 1);
  this->cs.resize(this->restart + 1);
  this->sn.resize(this->restart + 1);

  HessenbergMatrix H(this->restart + 1);

  // Vector r = rhs - A * solution;
  Vector r(rhs); // copy values

#ifdef _MPI
  // MPI: solution in consistency level 3 for vector.matrix multiplication
  for (size_t bl = 0; bl < comms.size(); ++bl)
  {
    comms[bl]->queue_consistency_update(solution.block(bl), 3);
  }

  TParFECommunicator3D::flush_consistency_updates();
#endif

  A.apply_scaled_add(solution, r, -1.0); // now r = rhs - A*solution

  // compute initial residual
#ifdef _MPI
  double resid = r.norm(comms);
#else
  double resid = r.norm();
#endif

  double beta = resid; // initialize beta as initial residual

  // save initial residual, used to check stopping criteria later
  if (this->converged(resid, 0))
  {
    return std::pair<unsigned int, double>(0, resid);
  }

  this->v.resize(this->restart + 1);

  // array to store the outputs of the preconditioning processes
  this->z.resize(this->restart + 1);

  unsigned int j = 1;
  while (j <= this->max_n_iterations)
  {
    v[0] = r;
    v[0] *= 1.0 / beta;
    std::fill(s.begin(), s.end(), 0.0); // set all entries in s to zero
    s[0] = beta;

    for(unsigned int i = 0;
      i < this->restart && j <= this->max_n_iterations;
      i++, j++)
    {
      // z[i] = A * M.solve(v[i]); // where M.solve(v[i]) is replaced by application of a preconditioning strategy
      z[i] = r; // copy structure of r
      z[i] = 0.0;

      // apply a preconditioning strategy with rhs v[i] to obtain z[i]
      this->prec->apply(i, j, v[i], z[i]);

#ifdef _MPI
      // MPI: solution in consistency level 3 for computation of global norm

      for (size_t bl = 0; bl < comms.size(); ++bl)
      {
        comms[bl]->queue_consistency_update(z[i].block(bl), 3);
      }

      TParFECommunicator3D::flush_consistency_updates();
#endif

      A.apply(z[i], r); // r = A * z[i]

#ifdef _MPI
      // v are orthogonal, so we can separate these loops to save sync latency
      for (unsigned int k = 0; k <= i; k++)
      {
        queue_dot_global(r, v[k], comms, H(k, i));
      }

      flush_dot_queue();

      for (unsigned int k = 0; k <= i; k++)
      {
        r.add_scaled(v[k], -H(k, i));
      }
#else
      for (unsigned int k = 0; k <= i; k++)
      {
        H(k, i) = dot(r, v[k]);
        r.add_scaled(v[k], -H(k, i));
      }
#endif

#ifndef _MPI
      H(i + 1, i) = r.norm();
#elif _MPI
      H(i + 1, i) = r.norm(comms);
#endif

      // can lead to "unhappy" crash in FGMRES if H(i+1, i) == 0
      // v[i+1] = r * (1.0 / H(i+1, i));

      if (H(i + 1, i) != 0.0)
      {
        v[i + 1] = r;
        v[i + 1] *= 1.0 / H(i + 1, i);
      }
      else
      {
        ErrThrow("Unhappy breakdown in flexible gmres at iteration ", j);
      }

      for (unsigned int k = 0; k < i; k++)
      {
        ApplyGivensRotation(H(k, i), H(k + 1, i), cs[k], sn[k]);
      }

      GenerateGivensRotation(H(i, i), H(i + 1, i), cs[i], sn[i]);
      ApplyGivensRotation(H(i, i), H(i + 1, i), cs[i], sn[i]);
      ApplyGivensRotation(s[i], s[i + 1], cs[i], sn[i]);

      resid = std::abs(s[i + 1]);

      if (this->converged(resid, j))
      {
        // use r as auxiliary here
        r = 0.0; // set entries to zero

        Update(r, i, H, s, z); // solve the least squares problem by backsolve
        solution += r; // and update the solution
        return std::pair<unsigned int, double>(j, resid);
      }
    }

    // use r as auxiliary here
    r = 0.0; // set entries to zero - just in case

    // solve the least squares problem by backsolve
    Update(r, this->restart - 1, H, s, z);
    solution += r; // and update the solution

#ifdef _MPI
    // MPI: solution in consistency level 3 for computation of global norm

    for (size_t bl = 0; bl < comms.size(); ++bl)
    {
      comms[bl]->queue_consistency_update(solution.block(bl), 3);
    }

    TParFECommunicator3D::flush_consistency_updates();
#endif

    // compute new residual r
    r = rhs;
    A.apply_scaled_add(solution, r, -1.0); // r = rhs - A * solution;

    // and new defect beta
#ifndef _MPI
    beta = r.norm();
#elif _MPI
    beta = r.norm(comms);
#endif

    if (std::abs(beta - resid) > 0.01 * beta)
    {
      Output::root_info<1>("FGMRES", "restart residual changed ", beta, "  ", resid);
    }

    if (this->converged(resid, j))
    {
      return std::pair<unsigned int, double>(j, resid);
    }
    // else restart
  }

  // not converged
  return std::pair<unsigned int, double>(this->max_n_iterations, resid);
}

/* ************************************************************************** */
// explicit instantiations
template class Iteration_gmres<BlockFEMatrix,   BlockVector>;
template class Iteration_gmres<CompositeOperator<BlockVector>, BlockVector>;

// In MPI case we are so dependent on the connection of Matrix and FESpace, that
// it does not make sense to instantiate the function for BlockMatrix.
#ifndef _MPI
template class Iteration_gmres<BlockMatrix,   BlockVector>;
#endif