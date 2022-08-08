#ifndef __ANALYTICALFUNCTION__
#define __ANALYTICALFUNCTION__

#include <functional> // std::function
#include <Point.h>
#include <FunctionEvaluation.h>
#include <MooNMD_Io.h>

using namespace parmoon; 

/** ************************************************************************
*
* @class      AnalyticalFunction
* @brief      handle an exact solution or initial data for an example
*
* Essentially this stores a function which can be evaluated at coordinates in
* the computational domain. There is a constructor for constant functions.
*
* For time dependent problems, this class stores a second function. There is no
* way to store both a stationary and in-stationary function in one element.
*
* This class provides a scalar exact solution. You can create a vector of
* objects of this class to handle vector problems.
*
* @author     Ulrich Wilbrandt
* @date       03.11.2014
*
****************************************************************************/
class AnalyticalFunction
{
 public:
  /** @brief create a dummy AnalyticalFunction object
   *
   * This is used to indicate that an exact solution is not known.
   */
  AnalyticalFunction();

  /** @brief construct a constant solution
   *
   * This is implemented only for simplicity.
   */
  explicit AnalyticalFunction(double a);

  /** @brief create an AnalyticalFunction object with a given (stationary)
   * function
   *
   * @note: if you use this for a time dependent solution it means this
   * AnalyticalFunction is constant in time.
   */
  explicit AnalyticalFunction(
    std::function<void(const Point&, FunctionEvaluation&)> f);

  /** @brief create an AnalyticalFunction object with a given in-stationary
   * function
   *
   * The double in std::function argument is the point in time.
   */
  explicit AnalyticalFunction(
      std::function<void(const Point&, double, FunctionEvaluation&)> f);

  /** @brief create a AnalyticalFunction which does not depend on time from 
   * another AnalyticalFunction at a particular point in time.
   * 
   * This can be used as an initial solution for a time dependent problem with
   * known solution. Then `time` should be the initial time and `other` that
   * known solution.
   * 
   * Note that if `other.exists() == false` this constructor makes no sense and
   * basically does what the default constructor does.
   */
  AnalyticalFunction(const AnalyticalFunction& other, double time);
  

  /** @brief evaluate this AnalyticalFunction at a point (stationary case)
   *
   * If this is called even if this object has been constructed using an
   * in-stationary std::function, then an exception is thrown.
   *
   * This fills the first function in the FunctionEvaluation object.
   */
  void get(const Point&, FunctionEvaluation&) const;

  /** @brief evaluate this AnalyticalFunction at a point in space and time
   *
   * If you call this method even if this object has been created using a
   * stationary std::function, then the time `t` is simply ignored
   *
   * This fills the first function in the FunctionEvaluation object.
   */
  void get(const Point&, double t, FunctionEvaluation&) const;

  /** @brief check if this exact function object really does describe a function
   *
   * If the default constructor was called to create this object, the following
   * will return false, otherwise true.
   * 
   * Note that if this returns false, this object otherwise behaves like a zero
   * function.
   */
  bool exists() const
  {
    return function.operator bool();
  }

  /// @brief indicate if this AnalyticalFunction is constant (in space and time)
  bool is_constant() const
  {
    if(!this->exists())
      Output::warn("calling AnalyticalFunction::is_constant on object which "
                   "does not have an analytical representation");
    return constant;
  }
  
  /// @brief indicate if this AnalyticalFunction depends on time
  bool depends_on_time = true;
 private:
  /** @brief the function evaluating the exact solution and its derivatives
   * at a particular point in time
   *
   * All derivatives up to order 2 are evaluated.
   */
  std::function<void(const Point&, double t, FunctionEvaluation&)> function;
  
  /** @brief indicate if the constructor for a constant AnalyticalFunction has 
   * been used.
   * 
   * This is helpful to find out proper BoundarData using this
   * AnalyticalFunction.
   */
  bool constant;
};

#endif // __ANALYTICALFUNCTION__
