#ifndef __EXAMPLE_NONSTATIONARY__
#define __EXAMPLE_NONSTATIONARY__

#include <Example.h>

class Example_NonStationary : public Example
{
 public:
  /** @brief define a time dependent convection-diffusion example
   *
   * The relevant parameter is called 'example' and switches between different
   * examples.
   */
  static Example_NonStationary Time_ConvDiff(const ParameterDatabase&);
  
  /** @brief define a time dependent Navier--Stokes example
   *
   * The relevant parameter is called 'example' and switches between different
   * examples.
   */
  static Example_NonStationary Time_NavierStokes(const ParameterDatabase&);
  
 protected:
  /** @brief create an example with all members given
   * @details This is used in the factory methods above.
   */
  Example_NonStationary(std::vector<BoundaryCondition>&& bc,
                        std::vector<BoundaryData>&& bd,
                        PDECoefficients&& coeffs,
                        std::vector<AnalyticalFunction>&& exact,
                        std::vector<AnalyticalFunction>&& initial_conditions);

 private:
  /** @brief the initial solution to this Example
   *
   * For scalar problems there is only one entry in this vector, for e.g. Stokes
   * in 2D there are three.
   *
   * The entries in this vector must exists, i.e., they must not be default 
   * constructed.
   */
  std::vector<AnalyticalFunction> initial;
};

#endif // __EXAMPLE_NONSTATIONARY__