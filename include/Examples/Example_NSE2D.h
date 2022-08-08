/** ************************************************************************ 
*
* @class Example_NSE2D
* @brief store all functions needed to describe a (Navier-)Stokes example 
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

#ifndef __EXAMPLE_NSE2D__
#define __EXAMPLE_NSE2D__

#include<Example2D.h>
#include <functional>

template <int d> class NavierStokes; //forward declaration

class Example_NSE2D : public Example2D
{
  public:
    /** @brief default constructor
     * 
     * This intializes a (Navier-)Stokes example in 2D. It is chosen according
     * to example_code.
     */
    explicit Example_NSE2D(const ParameterDatabase& user_input_parameter_db);

    /** @brief initialize your own example
     * 
     * Create an example with all vectors already defined.
     */
    Example_NSE2D(const std::vector<DoubleFunct2D*>& exact,
                  const std::vector<BoundCondFunct2D*>& bc,
                  const std::vector<BoundValueFunct2D*>& bd,
                  const CoeffFct2D& coeffs, double nu = 1.);
  
    /// Apply the function stored as post processing routine.
    void do_post_processing(NavierStokes<2>& nse2d) const;

    /// Return kinematic or effective viscosity, if set (depending on example code)
    double get_nu() const;

 
    /// Return permeability, if set.
    double get_inverse_permeability() const;

    /// Return neumann boundary ids, if set
    std::vector<size_t> get_neumann_id() const;

    /// Return nitsche boundary ids, if set
    std::vector<size_t> get_nitsche_id() const;

    /// Return windkessel boundary ids, if set
    std::vector<size_t> get_windkessel_id() const;

    //Declaration of special member functions - rule of zero

    //! Default copy constructor. Performs deep copy.
    Example_NSE2D(const Example_NSE2D&) = default;

    //! Default move constructor.
    Example_NSE2D(Example_NSE2D&&) = default;

    //! Default copy assignment operator. Performs deep copy.
    Example_NSE2D& operator=(const Example_NSE2D&) = default;

    //! Default move assignment operator
    Example_NSE2D& operator=(Example_NSE2D&&) = default;

    //! Default destructor.
    ~Example_NSE2D() = default;

  private:
  /// Function doing the post processing for a stationary example.
  /// TODO put NSE2D argument const as soon as FEFunctions can be copied properly!
  std::function<void(NavierStokes<2>&)> post_processing_stat;
  /// TODO Function doing the post processing for a time dependent example.
};


#endif // __EXAMPLE_NSE2D__
