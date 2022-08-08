#ifndef __EXAMPLE__
#define __EXAMPLE__

#include <BoundaryCondition.h>
#include <BoundaryData.h>
#include <PDECoefficients.h>
#include <AnalyticalFunction.h>
#include <ParameterDatabase.h>
#include <Domain.h>

/** ************************************************************************
*
* @class Example
* @brief Handle everything needed to set up an example both stationary and time
* dependent.
*
* Essentially this stores four objects: BoundaryCondition, BoundaryData,
* PDECoefficients, and AnalyticalFunction. The first three are mandatory.
*
* This class implements some examples for each type of Problem (e.g.
* Convection--Diffusion, Stokes, Darcy,...). If you want a different example
* you can implement it in this class (and the four classes mentioned above).
*
* We use the concept of named constructors to create objects for the different
* problem types (convection--diffusion, Stokes, Darcy, ...). These are
* implemented in different source files to keep their sizes smaller. There is
* currently no other way to construct an example.
*
* @author     Ulrich Wilbrandt
* @date       04.11.2014
*
****************************************************************************/
class Example
{
 public:
  /** @brief define a convection-diffusion example
   *
   * The relevant parameter is called 'example' and switches between different
   * examples.
   */
  static Example ConvDiff(const ParameterDatabase&);

  /** @brief define a darcy example
   *
   * The relevant parameter is called 'example' and switches between different
   * examples.
   */
  static Example Darcy(const ParameterDatabase&);

  /** @brief define a Navier--Stokes example
   *
   * The relevant parameter is called 'example' and switches between different
   * examples.
   */
  static Example NavierStokes(const ParameterDatabase&);

  /** @brief move constructor */
  Example(Example&&) = default;

  /** @brief return the boundary condition object */
  const BoundaryCondition& get_boundary_condition(unsigned int i = 0) const
  {
    return bc.at(i);
  };

  /** @brief return the boundary data object */
  const BoundaryData& get_boundary_data(unsigned int i = 0) const
  {
    return bd.at(i);
  }
  
  /** @brief return the vector of all boundary data objects */
  const std::vector<BoundaryData>& get_boundary_data_vector() const
  {
    return bd;
  }

  /** @brief return the coefficient object */
  const PDECoefficients& get_coefficients() const { return coeffs; };

  /** @brief return the exact solution object */
  const AnalyticalFunction& get_exact_solution(unsigned int i = 0) const
  {
    return exact.at(i);
  }

  /** @brief check if the domain described by a BoundaryDescription is suitable
   *         for this example.
   *
   * It is recommended to call this function even though a program will run fine
   * without it as well. If this function returns false it is no longer safe to
   * continue. For example if the domain has five boundary components while the
   * example expects only four, this is probably an error.
   */
  bool check_suitable_domain(const TDomain&) const;

  
  /// @brief extract the value of the Parameter `example` from a 
  ///        ParameterDatabase
  ///
  /// Basically this is equivalent to `param_db["example"]`.
  /// This function only exists to give a meaningful error message in case the
  /// Parameter `example` is missing. This is called in all the named
  /// constructors.
  static int get_example_from_database(const ParameterDatabase& param_db);
  
 protected:
  /** @brief create an example with all members given
   * @details This is used in the factory methods above.
   */
  Example(std::vector<BoundaryCondition>&& bc, std::vector<BoundaryData>&& bd,
          PDECoefficients&& coeffs, std::vector<AnalyticalFunction>&& exact);

  /** @brief boundary conditions
   *
   * For scalar problems there is only one entry in this vector, for e.g. Stokes
   * in 2D there are three.
   */
  std::vector<BoundaryCondition> bc;

  /** @brief boundary data
   *
   * For scalar problems there is only one entry in this vector, for e.g. Stokes
   * in 2D there are three.
   */
  std::vector<BoundaryData> bd;

  /** @brief the coefficients of the pde
   *
   * This vector should contain only one element. Maybe in the future there
   * might be examples with more than one element here for coupled problems.
   */
  PDECoefficients coeffs;

  /** @brief the exact solution to this Example
   *
   * For scalar problems there is only one entry in this vector, for e.g. Stokes
   * in 2D there are three.
   *
   * If no exact solution exists, you should still put an object of the class
   * AnalyticalFunction here. Use the standard constructor of that class which
   * indicates that no exact solution is available.
   */
  std::vector<AnalyticalFunction> exact;

  /** @brief a function which checks if the computational domain works with this
   *         example
   *
   * This function is not necessary for the program to run. However it is
   * recommended to define it to ensure that this example is not used on an
   * unsuitable domain. This function is called through the member method
   * check_suitable_domain. If it is not defined that method will always return
   * true.
   */
  std::function<bool(const TDomain&)> check_domain;
};

#endif // __EXAMPLE__
