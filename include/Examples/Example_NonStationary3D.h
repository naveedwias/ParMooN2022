/** ************************************************************************ 
*
* @class     Example_NonStationary3D.h
* @brief     store all functions needed to describe an example
* 
* Mainly this stores the exact solution (if available), boundary condition,
* boundary data, and the coefficients of the pde. Note that you almost always
* want to create an object of a derived class, rather than of type Example3D 
* itself.
* 
* Below we use std::vector in case there are multiple solution components (e.g.
* velocity components and pressure).
* 
* @date  
* 
* @ruleof0
* 
 ************************************************************************  */

#ifndef _Example_NonStationary3D_
#define _Example_NonStationary3D_

#include <Example3D.h>
#include <Constants.h>
#include <vector>

class Example_NonStationary3D : public Example3D
{
protected:
  /** @brief default constructor
   */
  explicit Example_NonStationary3D(const ParameterDatabase &);

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
  std::vector<DoubleFunct3D*> initialCondtion;
  
public:
  //Declaration of special member functions - rule of zero
  
  //! Default copy constructor. Performs deep copy.
  Example_NonStationary3D(const Example_NonStationary3D&) = default;

  //! Default move constructor.
  Example_NonStationary3D(Example_NonStationary3D&&) = default;
  
  //! Default copy assignment operator. Performs deep copy.
  Example_NonStationary3D& operator=(const Example_NonStationary3D&) = default;
  
  //! Default move assignment operator
  Example_NonStationary3D& operator=(Example_NonStationary3D&&) = default;
  
  //! Default destructor.
  ~Example_NonStationary3D() = default;
  /**
   * @brief constructor
   */
  Example_NonStationary3D(
    const std::vector<DoubleFunct3D*>& exact,
    const std::vector<BoundCondFunct3D*>& bc,
    const std::vector<BoundValueFunct3D*>& bd, const CoeffFct3D& coeffs,
    bool timedependentrhs = true, bool timedependentcoeffs = true,
    const std::vector<DoubleFunct3D*>& init_cond = std::vector<DoubleFunct3D*>());
  // getters
  //TODO
  bool rhs_depends_on_time();
  
  bool coefficients_depend_on_time();
  
  DoubleFunct3D* get_initial_cond(unsigned int i)const
    { return initialCondtion.at(i); }  
};

#endif // _Example_NonStationary3D_
