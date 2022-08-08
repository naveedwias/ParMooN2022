/** ************************************************************************ 
*
* @class Example_CD2D
* @brief store all functions needed to describe a convection--diffusion example 
* 
* The standard
* constructor of this class will fill the vectors (in Example2D) with pointers
* to the functions needed to fully describe a particular example.
* 
* So far the following examples are implemented and enabled:
* Stationary examples
* 0  -  Simple Sine Laplace (analytical example)
* 1  -  Two interior Layers  (Knopp, Lube, Rapin, CMAME 2002)
* 2  -  Hemker Rotating Body (Hemker 1996)
* 3  -  Sharp Boundary Layer (Kuzmin & Moeller 2005)
* 4  -  Smooth Solution (BJK SIAM, 2016)
* 5  -  Hughes, Mallet, Mizukami 1986
* 6  -  Inner Layer Hump, (John, Knobloch, Savescu, CMAME 2011)
* 7  -  Boundary Layer Known, (John, Knobloch, Savescu, CMAME 2011)
* 8  -  Known solution with a Boundary layer (Allendes, Barrenchea, and Rankin,
*       SISC 2017)
* 9  -  Rotating convection field (Allendes, Barrenchea, and Rankin, SISC 2017)
* 10 -  Polynomial of modifiable degree
* 11 -  Pure convection-reaction: 1+sin(...) (Houston, Schwab, Sueli, SIAM 2002)
* 12 -  Pure convection-reaction: Polynomial with jump discontinuity at 0.5 of
*       modifiable degree
* 13 -  Hughes, Mallet, Mizukami 1986, modified
*
* Time dependent examples
* 	101 - Exponential Function (analytical example)
* 	102 - Sine Sine Sine (analytical example)
* 	103 - Sine Cosine (analytical example)
*
*
* @author    Ulrich Wilbrandt, 
* @date      13.03.15
* 
* @ruleof0
* 
 ************************************************************************  */


#ifndef __EXAMPLE_CD2D__
#define __EXAMPLE_CD2D__

#include<Example2D.h>
#include <functional>

//class CD2D; //forward declaration
template<int d> class ConvectionDiffusion; // forward declaration


class Example_CD2D : public Example2D 
{
  public:
    /** @brief default constructor
     * 
     * This intializes a convection-diffusion example in 2D. It is chosen 
     * according to example_code.
     */
    explicit Example_CD2D(const ParameterDatabase& user_input_parameter_db);
    
    /** @brief initialize your own example
     * 
     * Create an example with all vectors already defined.
     */
    Example_CD2D(const std::vector<DoubleFunct2D*>& exact,
                 const std::vector<BoundCondFunct2D*>& bc,
                 const std::vector<BoundValueFunct2D*>& bd,
                 const CoeffFct2D& coeffs, double nu = 1.);

    /// Apply the function stored as post processing routine.
    void do_post_processing(ConvectionDiffusion<2>& cd2d) const;

    /// Return kinematic viscosity, if set.
    double get_nu() const;

    //Declaration of special member functions - rule of zero

    //! Default copy constructor. Performs deep copy.
    Example_CD2D(const Example_CD2D&) = default;

    //! Default move constructor.
    Example_CD2D(Example_CD2D&&) = default;

    //! Default copy assignment operator. Performs deep copy.
    Example_CD2D& operator=(const Example_CD2D&) = default;

    //! Default move assignment operator
    Example_CD2D& operator=(Example_CD2D&&) = default;

    //! Default destructor.
    ~Example_CD2D() = default;

  private:
    /// Function doing the post processing for a stationary example.
    /// TODO put CD2D argument const as soon as FEFunctions can be copied properly!
    std::function<void(ConvectionDiffusion<2> &)> post_processing_stat;
    /// TODO Function doing the post processing for a time dependent example.

};


#endif // __EXAMPLE_CD2D__
