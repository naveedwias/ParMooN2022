/** ************************************************************************ 
*
* @class Example_Darcy3D
* @brief store all functions needed to describe a Darcy example using vector
*        valued basis functions.
* 
* The standard
* constructor of this class will fill the vectors (in Example3D) with pointers
* to the functions needed to fully describe a particular example.
* 
* @ruleof0
* 
 ************************************************************************  */


#ifndef __EXAMPLE_DARCY3D__
#define __EXAMPLE_DARCY3D__

#include<Example3D.h>
#include <functional>

template <int d> class Darcy; //forward declaration

class Example_Darcy3D : public Example3D 
{
  public:
    /** @brief default constructor
     * 
     * This intializes a convection-diffusion example in 3D. It is chosen 
     * according to example_code
     */
    explicit Example_Darcy3D(const ParameterDatabase& user_input_parameter_db);
    
    /** @brief initialize your own example
     * 
     * Create an example with all vectors already defined.
     */
    Example_Darcy3D(const std::vector<DoubleFunct3D*>& exact,
                    const std::vector<BoundCondFunct3D*>& bc,
                    const std::vector<BoundValueFunct3D*>& bd,
                    const CoeffFct3D& coeffs)
    : Example3D(exact, bc, bd, coeffs) {};

    /// Apply the function stored as post processing routine.
    void do_post_processing(Darcy<3>& darcy2d) const;

    /// Return kinematic viscosity, if set.
    double get_nu() const;

    //Declaration of special member functions - rule of zero

    //! Default copy constructor. Performs deep copy.
    Example_Darcy3D(const Example_Darcy3D&) = default;

    //! Default move constructor.
    Example_Darcy3D(Example_Darcy3D&&) = default;

    //! Default copy assignment operator. Performs deep copy.
    Example_Darcy3D& operator=(const Example_Darcy3D&) = default;

    //! Default move assignment operator
    Example_Darcy3D& operator=(Example_Darcy3D&&) = default;

    //! Default destructor.
    ~Example_Darcy3D() = default;

  private:
    /// Function doing the post processing for a stationary example.
    /// TODO put Darcy3D argument const as soon as FEFunctions can be copied properly!
    std::function<void(Darcy<3> &)> post_processing_stat;
    /// TODO Function doing the post processing for a time dependent example.
};


#endif // __EXAMPLE_DARCY3D__
