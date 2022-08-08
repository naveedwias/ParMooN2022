#ifndef __COMPOSITE_OPERATOR__
#define __COMPOSITE_OPERATOR__

#ifdef _MPI
class TParFECommunicator3D;
#endif

#include <vector>
#include <string>
#include <functional>

// Square operator only.
template<class V>
class CompositeOperator
{
public:
  explicit CompositeOperator(
    std::function<void(const V&, V&, double)> op
#ifdef _MPI
    , const std::vector<const TParFECommunicator3D*> comms
#endif
    );

  CompositeOperator();

  CompositeOperator(const CompositeOperator& other) = default;

  // currently only functions in MPI
  size_t get_n_total_rows() const;

  // currently only functions in MPI
  size_t get_n_total_columns() const;

  void set_operator(std::function<void(const V&, V&, double)> new_op);

#ifdef _MPI
  void set_comms(const std::vector<const TParFECommunicator3D*> new_comms);
#endif

  void set_approximate_diagonal(const std::vector<double> &diag);
  void set_approximate_row_sums(const std::vector<double> &row_sums);
  void set_approximate_col_sums(const std::vector<double> &col_sums);

  void apply(const V &x, V &y) const;
  void apply_scaled_add(const V &x, V &y,
    double a = 1.0) const;

#ifdef _MPI
  std::vector<const TParFECommunicator3D*> get_communicators() const
  {
    return comms;
  }
#endif

  // inappropriate for CompositeOperator - will throw an error
#ifdef _MPI
  void sor_sweep(const V& b, V& x, double omega,
    size_t flag, const std::string& par_strat) const;
#else
  void sor_sweep(const V& b, V& x, double omega,
    size_t flag) const;
#endif

  // will throw an error if approximate diagonal has not been set
  std::vector<double> get_diagonal() const;

  // will throw an error if approximate row sums have not been set
  std::vector<double> get_row_sums() const;

  // will throw an error if approximate column sums have not been set
  std::vector<double> get_col_sums() const;

protected:
  std::function<void(const V&, V&, double)> op;

#ifdef _MPI
  std::vector<const TParFECommunicator3D*> comms;
#endif

  std::vector<double> approximate_diagonal;
  std::vector<double> approximate_row_sums;
  std::vector<double> approximate_col_sums;
};

#endif