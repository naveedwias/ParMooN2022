/** ************************************************************************ 
*
* @class Example_TimeCD2D
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
#ifndef _Example_TimeCD2D_
#define _Example_TimeCD2D_

#include <Example_NonStationary2D.h>
#include <functional>

class Time_CD2D; //forward declaration


class Example_TimeCD2D : public Example_NonStationary2D
{
public: 
  /**
   * @brief default constructor
   * This intializes a convection-diffusion example in 2D. It is chosen 
   * according to example_code.
   */
  explicit Example_TimeCD2D(const ParameterDatabase& user_input_parameter_db);
  /** @brief initialize your own example
   * 
   * Create an example with all vectors already defined.
   */
  Example_TimeCD2D(const std::vector<DoubleFunct2D*>& exact,
                   const std::vector<BoundCondFunct2D*>& bc,
                   const std::vector<BoundValueFunct2D*>& bd, 
                   const CoeffFct2D& coeffs,                    
                   bool timedependentrhs, bool timedependentcoeffs, 
                   std::vector <DoubleFunct2D*> init_cond)
  : Example_NonStationary2D(exact, bc, bd, coeffs,  timedependentrhs, 
                          timedependentcoeffs, init_cond)
  {

  };

  /// Apply the function stored as post processing routine.
  void do_post_processing(Time_CD2D& tcd2d) const;

  /// Return diffusion coefficient, if set.
  double getDiffCoeff() const;

  //Declaration of special member functions - rule of zero
  //! Default copy constructor. Performs deep copy.
  Example_TimeCD2D(const Example_TimeCD2D&) = default;
  
  //! Default move constructor.
  Example_TimeCD2D(Example_TimeCD2D&&) = default;
  
  //! Default copy assignment operator. Performs deep copy.
  Example_TimeCD2D& operator=(const Example_TimeCD2D&) = default;
  
  //! Default move assignment operator
  Example_TimeCD2D& operator=(Example_TimeCD2D&&) = default;
  
  //! Default destructor.
  ~Example_TimeCD2D() = default;  

  private:
  /// Function doing the post processing for a stationary example.
  /// TODO put Time_CD2D argument const as soon as FEFunctions can be copied properly!
  std::function<void(Time_CD2D &)> post_processing_stat;
  /// TODO Function doing the post processing for a time dependent example.

};
#endif
