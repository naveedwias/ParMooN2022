#ifndef LPS_SCOTT_ZHANG_H
#define LPS_SCOTT_ZHANG_H

#include <memory>
class FEMatrix;

struct LPS_parameter_set
{
  size_t lps_coeff_type;
  double delta0;
  double delta1;
};

// this replaces the matrix C!
std::shared_ptr<FEMatrix> LPS_for_pressure_Scott_Zhang(
  const std::shared_ptr<const FEMatrix>& C, bool velocity, double nu,
  const LPS_parameter_set& lps_ps);

#endif // LPS_SCOTT_ZHANG_H
