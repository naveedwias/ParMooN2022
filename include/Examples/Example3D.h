/** ************************************************************************ 
*
* @class     Example3D
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
* @date      01.06.2015
* 
* @ruleof0
* 
 ************************************************************************  */

#ifndef __EXAMPLE3D__
#define __EXAMPLE3D__

#include <ParameterDatabase.h>
#include <MooNMD_Io.h>
#include <Constants.h>
#include <vector>


class Example3D 
{
  protected:
    /**
    * @brief default constructor 
    * 
    * This is used only by the classes derived from this class.
    */
    explicit Example3D(const ParameterDatabase &);

    /** @brief a local parameter database which controls this class
     *
     * The database given to the constructor will be merged into this one. Only
     * parameters which are of interest to this class are stored (and the
     * default ParMooN parameters).
     */
    ParameterDatabase example_database;

  public:
    /** @brief initialize your own example
     * 
     * Create an example with all vectors already defined.
     */
    Example3D(const std::vector<DoubleFunct3D*>& exact,
              const std::vector<BoundCondFunct3D*>& bc,
              const std::vector<BoundValueFunct3D*>& bd,
              const CoeffFct3D& coeffs);

    /* functions representing the exact solution */
    std::vector <DoubleFunct3D*> exact_solution;
    /* functions representing the boundary conditions */
    std::vector <BoundCondFunct3D*> boundary_conditions;
    /* functions representing the boundary data */
    std::vector <BoundValueFunct3D*> boundary_data;
    /* functions representing the coefficients of the pde */
    CoeffFct3D problem_coefficients;
    
    //Declaration of special member functions - rule of zero

    //! Default copy constructor. Performs deep copy.
    Example3D(const Example3D&) = default;

    //! Default move constructor.
    Example3D(Example3D&&) = default;

    //! Default copy assignment operator. Performs deep copy.
    Example3D& operator=(const Example3D&) = default;

    //! Default move assignment operator
    Example3D& operator=(Example3D&&) = default;

    //! Default destructor.
    ~Example3D() = default;

    // Initialize example database, called with the constructor
    static ParameterDatabase default_example_database();

    // Getter functions

    const std::vector <DoubleFunct3D*> & get_exact() const 
    { return exact_solution; }

    DoubleFunct3D* get_exact(unsigned int i) const
    { return exact_solution.at(i); }

    BoundCondFunct3D* const * get_bc() const
    { return &boundary_conditions[0]; }

    BoundCondFunct3D* get_bc(unsigned int i) const
    { return boundary_conditions.at(i); }

    BoundValueFunct3D* const * get_bd() const
    { return &boundary_data[0]; }

    BoundValueFunct3D* get_bd(unsigned int i) const
    { return boundary_data.at(i); }

    const CoeffFct3D& get_coeffs() const
    { return problem_coefficients; }

    const ParameterDatabase & get_database() const
    { return example_database; }
    
    ParameterDatabase & get_database()
    { return example_database; }
    
    int get_n_bd_fcts() const
    { return boundary_data.size(); }
};

#endif // __EXAMPLE3D__

