#include <CompositeOperator.h>
#include <BlockVector.h>
#include <MooNMD_Io.h>

#ifdef _MPI
#include <mpi.h>
#include <ParFECommunicator3D.h>
#endif

template<class V>
void ZeroOperator(const V&, V &y, double)
{
  y.reset();
}

template<class V>
CompositeOperator<V>::CompositeOperator(
  std::function<void(const V&, V&, double)> op
#ifdef _MPI
  , const std::vector<const TParFECommunicator3D*> comms
#endif
  )
: op(op)
#ifdef _MPI
  , comms(comms)
#endif
{
}

template<class V>
CompositeOperator<V>::CompositeOperator()
: op(ZeroOperator<V>)
#ifdef _MPI
  , comms()
#endif
{
}

template<class V>
void CompositeOperator<V>::set_operator(
  std::function<void(const V&, V&, double)> new_op)
{
  op = new_op;
}

#ifdef _MPI
template<class V>
void CompositeOperator<V>::set_comms(
  const std::vector<const TParFECommunicator3D*> new_comms)
{
  comms = new_comms;
}
#endif

template<class V>
void CompositeOperator<V>::set_approximate_diagonal(const std::vector<double> &diag)
{
  approximate_diagonal = diag;
}

template<class V>
void CompositeOperator<V>::set_approximate_row_sums(const std::vector<double> &row_sums)
{
  approximate_row_sums = row_sums;
}

template<class V>
void CompositeOperator<V>::set_approximate_col_sums(const std::vector<double> &col_sums)
{
  approximate_col_sums = col_sums;
}

template<class V>
void CompositeOperator<V>::apply(const V &x, V &y) const
{
  y.reset();
  op(x, y, 1.0);
}

template<class V>
void CompositeOperator<V>::apply_scaled_add(const V &x, V &y, double a) const
{
  op(x, y, a);
}

template<class V>
size_t CompositeOperator<V>::get_n_total_rows() const
{
#ifdef _MPI
  size_t n = 0;

  for (auto* c: comms)
  {
    n += c->GetNDof();
  }

  return n;
#else
  return 0;
#endif
}

template<class V>
size_t CompositeOperator<V>::get_n_total_columns() const
{
  return get_n_total_rows();
}

template<class V>
void CompositeOperator<V>::sor_sweep(const V&, V&, double, size_t
#ifdef _MPI
  , const std::string&
#endif
  ) const
{
  ErrThrow("Cannot apply SOR to a composite operator!");
}

template<class V>
std::vector<double> CompositeOperator<V>::get_diagonal() const
{
  if (approximate_diagonal.size() == 0)
  {
    ErrThrow("Cannot get a composite operator's diagonal! An approximation "
      "must be set manually.");
  }

  return approximate_diagonal;
}

template<class V>
std::vector<double> CompositeOperator<V>::get_row_sums() const
{
  if (approximate_row_sums.size() == 0)
  {
    ErrThrow("Cannot get a composite operator's row sums! An approximation "
      "must be set manually.");
  }

  return approximate_row_sums;
}

template<class V>
std::vector<double> CompositeOperator<V>::get_col_sums() const
{
  if (approximate_col_sums.size() == 0)
  {
    ErrThrow("Cannot get a composite operator's column sums! An approximation "
      "must be set manually.");
  }

  return approximate_col_sums;
}

template class CompositeOperator<BlockVector>;
template void ZeroOperator<BlockVector>(const BlockVector&, BlockVector&, double);