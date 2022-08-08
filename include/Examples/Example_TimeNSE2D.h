/** ************************************************************************ 
*
* @class Example_TimeNSE2D
* @brief store all functions needed to describe a (Navier-)Stokes example 
* 
* The standard
* constructor of this class will fill the vectors (in Example2D_TimeProblem2D) with pointers
* to the functions needed to fully describe a particular example.
* 
* @author 
* @date 
* 
* @ruleof0
* 
 ************************************************************************  */
#ifndef _Example_TimeNSE2D_
#define _Example_TimeNSE2D_

#include <Example_NonStationary2D.h>
#include <functional>

template <int d> class TimeNavierStokes; //forward declaration

class Example_TimeNSE2D : public Example_NonStationary2D
{
public:
  /** @brief default constructor
   * This intializes a Time dependent (Navier-)Stokes example in 2D. 
   * It is chosen according to example_code.
   */
  explicit Example_TimeNSE2D(const ParameterDatabase& user_input_parameter_db);
  /** @brief initialize your own example
   * 
   * Create an example with all vectors already defined.
   */
  Example_TimeNSE2D(const std::vector<DoubleFunct2D*>& exact,
                    const std::vector<BoundCondFunct2D*>& bc,
                    const std::vector<BoundValueFunct2D*>& bd, 
                    const CoeffFct2D& coeffs,
                    bool timedependentrhs, bool timedependentcoeffs,
                    const std::vector<DoubleFunct2D*>& init_cond);

  /// Apply the function stored as post processing routine.
  void do_post_processing(TimeNavierStokes<2>& tnse2d, double& val) const;
  void do_post_processing(TimeNavierStokes<2>& tnse2d) const;

  /// Return kinematic viscosity, if set.
  double get_nu() const;

  /// Return permeability, if set.
  double get_inverse_permeability() const;

  /// Return neumann boundary ids, if set
  std::vector<size_t> get_neumann_id() const;

  /// Return nitsche boundary ids, if set
  std::vector<size_t> get_nitsche_id() const;

  /// Return windkessel boundary ids, if set
  std::vector<size_t> get_windkessel_id() const;

  
  private:
  /// Function doing the post processing for a stationary example.
  /// TODO put Time_NSE2D argument const as soon as FEFunctions can be copied properly!
  std::function<void(TimeNavierStokes<2>&, double& val)> post_processing_stat;
  
  std::function<void(TimeNavierStokes<2>&)> post_processing_stat_old;
  /// TODO Function doing the post processing for a time dependent example.

};
#endif // _Example_TimeNSE2D_
