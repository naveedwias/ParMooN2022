/** ******************************************************************
* @class Example_NonStationary2D
* @brief 
*
* @author 
* @date: 18.05.2016
*/

#ifndef __EXAMPLE_TIMEPROB2D__
#define __EXAMPLE_TIMEPROB2D__

#include <Constants.h>
#include <Example2D.h>
#include <vector>

class Example_NonStationary2D : public Example2D
{
  protected:
    /**
     * @brief default constructor
     * 
     */
    explicit Example_NonStationary2D(const ParameterDatabase &);
    
    /** @brief right hand side vector depends on time or not
     * default is true
     */
    bool timeDependentRhs;
    /** @brief problem coefficient dependency on time
     * default is true
     */
    bool timeDependentCoeffs;
    
    /** @brief function representing the initial data
     */
    std::vector<DoubleFunct2D*> initialCondition;
  public:
    
    //Declaration of special member functions - rule of zero

    //! Default copy constructor. Performs deep copy.
    Example_NonStationary2D(const Example_NonStationary2D&) = default;

    //! Default move constructor.
    Example_NonStationary2D(Example_NonStationary2D&&) = default;

    //! Default copy assignment operator. Performs deep copy.
    Example_NonStationary2D& operator=(const Example_NonStationary2D&) = default;

    //! Default move assignment operator
    Example_NonStationary2D& operator=(Example_NonStationary2D&&) = default;

    //! Default destructor.
    ~Example_NonStationary2D() = default;
    
    /**
     * @brief constructor
     */
    Example_NonStationary2D(
      const std::vector<DoubleFunct2D*>& exact,
      const std::vector<BoundCondFunct2D*>& bc,
      const std::vector<BoundValueFunct2D*>& bd, const CoeffFct2D& coeffs,
      bool timedependentrhs = true, bool timedependentcoeffs = true,
      std::vector<DoubleFunct2D*> init_cond = std::vector<DoubleFunct2D*>());
    
    
    // getters
    //
    bool get_rhs_depends_on_time(){return timeDependentRhs;};
    
    bool get_coefficients_depend_on_time(){return timeDependentCoeffs;};
    
    DoubleFunct2D* get_initial_cond(unsigned int i)const
    { return initialCondition.at(i); }
};

#endif
