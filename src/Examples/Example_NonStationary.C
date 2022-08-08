#include <Example_NonStationary.h>

Example_NonStationary::Example_NonStationary(
  std::vector<BoundaryCondition> && bc, std::vector<BoundaryData> && bd,
  PDECoefficients && coeffs, std::vector<AnalyticalFunction> && exact,
  std::vector<AnalyticalFunction> && initial_conditions)
 : Example(std::move(bc), std::move(bd), std::move(coeffs), std::move(exact)),
   initial(initial_conditions)
{
  if(this->initial.size() == 0)
  {
    ErrThrow("No initial condition given for Example_NonStationary "
             "constructor");
  }
  for(auto ic : this->initial)
  {
    if(!ic.exists())
    {
      ErrThrow("No explicit initial condition given. They must not be default "
               "constructed");
    }
    if(ic.depends_on_time)
    {
      Output::warn("Example_NonStationary", "The given initial condition "
                   "depends on time. This may be an error");
    }
  }
}
