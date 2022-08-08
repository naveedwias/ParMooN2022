/** ************************************************************************ 
*
* @class Example_Darcy2D
* @brief store all functions needed to describe a Darcy example using vector
*        valued basis functions.
* 
* The standard
* constructor of this class will fill the vectors (in Example2D) with pointers
* to the functions needed to fully describe a particular example.
* 
* @author    Ulrich Wilbrandt, 
* @date      14.03.15
* 
* @ruleof0
* 
 ************************************************************************  */


#ifndef __EXAMPLE_DARCY2D__
#define __EXAMPLE_DARCY2D__

#include <Example2D.h>
#include <functional>

template <int d> class Darcy; //forward declaration

class Example_Darcy2D : public Example2D 
{
  public:
    /** @brief default constructor
     * 
     * This intializes a convection-diffusion example in 2D. It is chosen 
     * according to example_code
     */
    explicit Example_Darcy2D(const ParameterDatabase& user_input_parameter_db);
    
    /** @brief initialize your own example
     * 
     * Create an example with all vectors already defined.
     */
    Example_Darcy2D(const std::vector<DoubleFunct2D*>& exact,
                    const std::vector<BoundCondFunct2D*>& bc,
                    const std::vector<BoundValueFunct2D*>& bd,
                    const CoeffFct2D& coeffs)
    : Example2D(exact, bc, bd, coeffs) {};

    /// Apply the function stored as post processing routine.
    void do_post_processing(Darcy<2>& darcy2d) const;

    /// Return kinematic viscosity, if set.
    double get_nu() const;

    //Declaration of special member functions - rule of zero

    //! Default copy constructor. Performs deep copy.
    Example_Darcy2D(const Example_Darcy2D&) = default;

    //! Default move constructor.
    Example_Darcy2D(Example_Darcy2D&&) = default;

    //! Default copy assignment operator. Performs deep copy.
    Example_Darcy2D& operator=(const Example_Darcy2D&) = default;

    //! Default move assignment operator
    Example_Darcy2D& operator=(Example_Darcy2D&&) = default;

    //! Default destructor.
    ~Example_Darcy2D() = default;

  private:
    /// Function doing the post processing for a stationary example.
    /// TODO put Darcy2D argument const as soon as FEFunctions can be copied properly!
    std::function<void(Darcy<2> &)> post_processing_stat;
    /// TODO Function doing the post processing for a time dependent example.

};


#endif // __EXAMPLE_DARCY2D__
