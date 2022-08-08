/** ************************************************************************ 
*
* @class Example_TimeNSE3D
* @brief store all functions needed to describe a (Navier-)Stokes example 
* 
* The standard constructor of this class will fill the vectors (in Example2D) 
* with pointers to the functions needed to fully describe a particular example.
* 
* @author 
* @date   
* 
* @ruleof0
* 
 ************************************************************************  */

#ifndef _EXAMPLE_TimeNSE3D_
#define _EXAMPLE_TimeNSE3D_

#include <Example_NonStationary3D.h>
#include <functional>

template <int d> class TimeNavierStokes; //forward declaration

class Example_TimeNSE3D : public Example_NonStationary3D
{
public:
  /** @brief default constructor
   * 
   * This intializes a (Navier-)Stokes example in 2D. It is chosen according
   * to example_code.
   */
  explicit Example_TimeNSE3D(const ParameterDatabase& user_input_parameter_db);
  
  /** @brief initialize your own example
   * 
   * Create an example with all vectors already defined.
   */
  Example_TimeNSE3D(const std::vector<DoubleFunct3D*>& exact,
                    const std::vector<BoundCondFunct3D*>& bc,
                    const std::vector<BoundValueFunct3D*>& bd, 
                    const CoeffFct3D& coeffs,                    
                    bool timedependentrhs, bool timedependentcoeffs, 
                    const std::vector<DoubleFunct3D*>& init_cond)
  : Example_NonStationary3D(exact, bc, bd, coeffs,  timedependentrhs, 
                          timedependentcoeffs, init_cond) 
  {
  };
  
  /// Apply the function stored as post processing routine.
  void do_post_processing(TimeNavierStokes<3>& tnse3d, double& val) const;

  /// Return reynolds number, if set.
  double get_reynolds() const;

  /// Return kinematic viscosity, if set.
  double get_nu() const;
  double get_effective_viscosity() const;

  /// return Neumann boundary labels
  std::vector<size_t> get_neumann_id() const;

  /// return Nitsche boundary labels
  std::vector<size_t> get_nitsche_id() const;

  /// return Windkessel boundary labels
  std::vector<size_t> get_windkessel_id() const;

  //Declaration of special member functions - rule of zero
    
  //! Default copy constructor. Performs deep copy.
  Example_TimeNSE3D(const Example_TimeNSE3D&) = default;
  
  //! Default move constructor.
  Example_TimeNSE3D(Example_TimeNSE3D&&) = default;
  
  //! Default copy assignment operator. Performs deep copy.
  Example_TimeNSE3D& operator=(const Example_TimeNSE3D&) = default;
  
  //! Default move assignment operator
  Example_TimeNSE3D& operator=(Example_TimeNSE3D&&) = default;
  
  //! Default destructor.
  ~Example_TimeNSE3D() = default;

private:
  /// Function doing the post processing for a stationary example.
  /// TODO @ULRICH: 
  // put NSE3D argument const as soon as FEFunctions can be copied properly!
  std::function<void(TimeNavierStokes<3> &)> post_processing_stat;
  /// TODO @ULRICH Function doing the post processing for a time dependent example.
  
};
#endif // _EXAMPLE_TimeNSE3D_
