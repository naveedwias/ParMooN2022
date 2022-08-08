/** ************************************************************************ 
*
* @class Example_TimeCD3D
* @brief store all functions needed to describe a convection--diffusion example 
* 
* Depending on the value of TDatabase::ParamDB->EXAMPLE, the standard 
* constructor of this class will fill the vectors (in Example2D_TimeProblem2D) with pointers
* to the functions needed to fully describe a particular example.
* 
* So far the following examples are implemented and enabled:
* 	0 - 
* 	1 - 
* 	2 - 
*
*
* @author  
* @date    
* 
* @ruleof0
* 
 ************************************************************************  */
#ifndef _Example_TimeCD3D_
#define _Example_TimeCD3D_

#include <Example_NonStationary3D.h>
#include <functional>  // std::function

class Time_CD3D; // forward declaration

class Example_TimeCD3D : public Example_NonStationary3D
{
public:
    /**
   * @brief default constructor
   * This intializes a convection-diffusion example in 2D. It is chosen 
   * according to example_code.
   */
  explicit Example_TimeCD3D(const ParameterDatabase& user_input_parameter_db);
  
  /** @brief initialize your own example
   * 
   * Create an example with all vectors already defined.
   */
  Example_TimeCD3D(const std::vector<DoubleFunct3D*>& exact,
                   const std::vector<BoundCondFunct3D*>& bc,
                   const std::vector<BoundValueFunct3D*>& bd, 
                   const CoeffFct3D& coeffs,                    
                   bool timedependentrhs, bool timedependentcoeffs, 
                   const std::vector<DoubleFunct3D*>& init_cond)
  : Example_NonStationary3D(exact, bc, bd, coeffs,  timedependentrhs, 
                          timedependentcoeffs, init_cond) {};
  
  /// Apply the function stored as post processing routine.
  void do_post_processing(Time_CD3D& tcd3d) const;

  /// Return diffusion coefficient, if set.
  double getDiffCoeff() const;

  //Declaration of special member functions - rule of zero
  //! Default copy constructor. Performs deep copy.
  Example_TimeCD3D(const Example_TimeCD3D&) = default;
  
  //! Default move constructor.
  Example_TimeCD3D(Example_TimeCD3D&&) = default;
  
  //! Default copy assignment operator. Performs deep copy.
  Example_TimeCD3D& operator=(const Example_TimeCD3D&) = default;
  
  //! Default move assignment operator
  Example_TimeCD3D& operator=(Example_TimeCD3D&&) = default;
  
  //! Default destructor.
  ~Example_TimeCD3D() = default;  

private:
  /// Function doing the post processing for a stationary example.
  /// TODO put TCD3D argument const as soon as FEFunctions can be copied properly!
  std::function<void(Time_CD3D &)> post_processing_stat;
  /// TODO Function doing the post processing for a time dependent example.
};
#endif // _Example_TimeCD3D_
