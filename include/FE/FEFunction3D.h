#ifndef __FEFUNCTION3D__
#define __FEFUNCTION3D__

class TAuxParam3D;
#include "FESpace3D.h"
#include <memory>
#include <vector>

/** @brief a function from a finite element space 
 * 
 * This class essentially consists of two pointers: the finite element space and
 * the finite element vector. It represents the isomorphism between 
 * \f$\mathbb{R}^N\f$ and the finite element (function) space. Therefore, this
 * class offers methods to evaluate this function at particular points, to
 * compute integrals, and to do interpolations.
 */
class TFEFunction3D
{
  public:
    typedef void ErrorMethod(int, std::array<const double*, 3>, const double *,
                             const double *, double, const double *const*,
                             const double *const*, const double *const*,
                             double *);
  protected:
    /** @brief name of the function */
    std::string Name;

    /** @brief finite element space to which this function belongs to */
    std::shared_ptr<const TFESpace3D> FESpace3D;

    /** @brief pointer to finite element vector according to FE isomorphism
     * 
     * So this vector consists of the weights for each basis function of the
     * corresponding finite element space.
     */
    double *Values;

  public:
    /// @brief Default constructor. All empty object.
    TFEFunction3D();

    /** @brief constructor with vector initialization 
     * 
     * The pointer `values` must point to an array of length (at least) 
     * `fespace3D->get_n_dof()`.
     */
    TFEFunction3D(std::shared_ptr<const TFESpace3D> fespace3D,
                  const std::string& name, double *values);

    /** @brief (shallow) copy assignment operator. 
     * 
     * Shallow copy, as the FEFunction does not take any memory responsibility.
     * That means no data is copied (except the `name`), only the pointers are
     * overwritten.
     */
    TFEFunction3D& operator=( const TFEFunction3D & );

    /** @brief return name of this finite element function */
    std::string GetName() const
    { return Name; }

    /** @brief reset name */
    void set_name(std::string name)
    { Name = name; }

    /** @brief return 3D finite element space */
    std::shared_ptr<const TFESpace3D> GetFESpace3D() const
    { return FESpace3D; }

    /** @brief return finite element space */
    std::shared_ptr<const TFESpace> GetFESpace() const
    { return FESpace3D; }

    /** @brief return length of finite element vector (number of dofs) */
    int GetLength() const
    { return FESpace3D->get_n_dof(); }

    /** @brief return pointer to finite element vector (const version) */
    const double *GetValues() const
    { return Values; }

    /** @brief return pointer to finite element vector */
    double *GetValues()
    { return Values; }


    /** @brief calculate errors to given function 
     * 
     * @warning The array \p errors has to be of length at least `N_Errors+1` !!
     */
    void GetErrors(DoubleFunct3D *Exact, int N_Derivatives,
                   MultiIndex3D *NeededDerivatives,
                   int N_Errors, ErrorMethod *ErrorMeth, 
                   CoeffFct3D Coeff, TAuxParam3D *Aux,
                   int n_fespaces, const TFESpace3D **fespaces,
                   double *errors,
                   std::function<bool(const TBaseCell*, int)>funct =
                   [](const TBaseCell*, int){return false;}) const;

    /** @brief calculate errors to given function
     * 
     * Use this for vector valued basis functions (Raviart-Thomas (RT) or 
     * Brezzi-Douglas-Marini (BDM) elements) 
     */
    void GetErrorsForVectorValuedFunction(DoubleFunct3D * const * const Exact,
                                          ErrorMethod * const ErrMeth,
                                          double * const errors);

    /** calculate errors to given function */
    void GetMeshCellParams(DoubleFunct3D *Exact, int N_Derivatives,
                   MultiIndex3D *NeededDerivatives,
                   int N_Errors, ErrorMethod *ErrorMeth, 
                   CoeffFct3D Coeff, TAuxParam3D *Aux,
                   int n_fespaces, const TFESpace3D **fespaces,
                   double *errors, double *cell_parameters);

    /**
     * @brief Calculates value of this FEFunction at a given point and the
     * gradient recovery (arithmetic mean of the gradient value in all
     * containing mesh cells).
     *
     * @param[in] x x value of the point at which to evaluate
     * @param[in] y y value of the point at which to evaluate
     * @param[in] z z value of the point at which to evaluate
     *
     * @param[out] values A vector of length 4 (checked)! Will be filled with
     *            function value
     *            gradient recovery dx
     *            gradient recovery dy
     *            gradient recovery dz
     *
     * @return True if the point was found in the mesh cell collection
     * of this function, false if not so.
     * 
     * \note This method first searches for the cell in which the point is. 
     * Therefore it is rather expensive. If possible, use `FindGradientLocal`.
     */
    bool FindGradient(double x, double y, double z, std::vector<double>& values) const;

    /** @brief Determines the value of function and its first derivatives at
     * the given point lying WITHIN this cell (NOT on its boundary) */
    void FindGradientLocal(int cell_no, 
                           double x, double y, double z, 
                           double *values) const;

    /** @brief determine the value of function at the given point 
     * 
     * The `cell_no` must be the index of the `cell`.
     * The point (`x`, `y`) must be in the given cell.
     * 
     * @param[out] values values[0] contains the desired value. In case of 
     * vector-valued basis functions, the vector (values[0],values[1]) is the
     * desired values of this finite element function.
     */
    void FindValueLocal(const TBaseCell *cell, int cell_no, 
                        double x, double y, double z, 
                        double *values) const;

    /** @brief calculate the interpolation of an exact function */
    void Interpolate(DoubleFunct3D *Exact);

    /** @brief interpolate the old mesh fe function values to the new fe function
     * 
     * Note that this is rather slow, because no further information is 
     * required. The function 'OldFeFunction' could even live on a larger domain.
     */
    void Interpolate(TFEFunction3D *funct);
    
    /** @brief interpolate a vector valued function
     *
     * @param[in] Exact must be of length 3 (= space dimension)
     * @warning EvalAll must be correctly implemented for the used finite element
     */
    void Interpolate_vector_valued_function(
      const std::vector<DoubleFunct3D*>& Exact);
    
    typedef std::function<double(const TBaseCell* cell, int cell_index,
                                 std::array<double, 3> xyz)> AnalyticFunction;
    /**
     * @brief add a given function f to this fe function
     * 
     * The cell and cell_index are from the collection of this TFEFunction3D
     * and the point (x,y,z) is in that cell.
     * 
     * Note that this is similar to creating a second TFEFunction3D, 
     * interpolating f on it and then adding it via operator+=. Here, no second
     * TFEFunction3D is required.
     */
    void add(AnalyticFunction f);

    /**
     * @brief add a constant to the FEFunction3D (mainly for the pressure)
     *
     * @param[in] b is the constant to add to this FEFunction3D
     */
    void add_constant(double b = 0.0);

    /**
     * @brief project this functions into the space L20 (having zero mean value)
     *
     * After a call to this function the mean value (integral of this function
     * devided by the measure of its domain) has the value a. This is for
     * example needed for the pressure in a Stokes problem with Dirichlet
     * conditions on all boundaries.
     *
     * @param[in] a set mean value of this FEFunction3D to a
     */
    void project_into_L20(double a = 0.0);

    /**
     * @brief find the integral of this function and the measure of its domain
     *
     * In MPI case it returns the global integral and measure, summed up over
     * the own domains of all processes.
     *
     * @param[out] integral double value for the integral of this TFEFunction3D
     * @param[out] measure double value for the measure of its domain
     */
    void compute_integral_and_measure(double& integral, double& measure) const;

    /**
     * @brief compute the mean oscillation
     *
     * \f[ \textrm{osc}_{\textrm{mean}}(u_h) =
     *    \frac{1}{\left| \mathcal{T}_h \right|} \sum_{K \in \mathcal{T}_h}
     *    \max\left\{ 0, \max_K u_h - u_{\max} \right\}
     *   +\max\left\{ 0, u_{\min} - \min_K u_h \right\}.
     * \f]
     *
     * See also equation (15) in:
     * Derk Frerichs, Volker John: On reducing spurious oscillations 
     * in discontinuous Galerkin (DG)methods for steady-state 
     * convection-diffusion-reaction equations. 2021
     *
     * If \p use_pk_nodal_fctn is false (default), then an approximation of the
     * minimum and maximum of uh is computed using many function evaluations of
     * uh on all cells. If this parameter is set to true, it first transforms
     * the local FE-vector to using Pk elements of the same order and then takes
     * takes the minimal and maximal value of this transformed FE-vector, i.e.
     * minimum and maximum of the Pk nodal functionals evaluated at uh.
     *
     * @warning The case that \p use_pk_nodal_fctn is true works correctly only
     * if the reference transformation preserves values, e.g. for DG elements,
     * but not Hdiv elements.
     */
    double compute_mean_oscillation(double u_min, double u_max, const bool&
        use_pk_nodal_fctn = false) const;

    /// @brief compute the L^2 norm.
    double get_L2_norm() const;
    
    /// @brief compute the L^2 norm and the H1-semi norm.
    std::pair<double,double> get_L2_and_H1_norm() const;

    /** @brief Set Dirichlet values according to boundary conditions */
    void SetDirichletBC(BoundCondFunct3D *BoundaryCondition,
                        BoundValueFunct3D *BoudaryValue);

    /** @brief computes the largest and smallest value of the FE function in
     * a given cell
     *
     * If \p use_pk_nodal_fctn is set to false, it computes an approximation of
     * the largest and smallest value of the FE function by evaluating the FE
     * Function at around 900 points in triangles or around 1800 points in
     * quads.
     * If the argument is set to true, it uses the nodal functionals of
     * the Pk element with the corresponding order, i.e. point evaluations at
     * the degrees of freedom of the corresponding Pk element.
     *
     * To speed up the computations if \p use_pk_nodal_fctn is false a cache
     * function is used that updates the values in the vector \bf_values if \p
     * current_type has changed, which is why these two inputs are needed for
     * this method.
     *
     * @warning The case that \p use_pk_nodal_fctn is true works correctly only
     * if the reference transformation preserves values, e.g. for DG elements,
     * but not Hdiv elements.
     */
    std::pair<double, double> compute_cell_min_max(int cell_nr,
        std::vector<double>& bf_values, BaseFunction_type& current_type, bool
        use_pk_nodal_fctn = false) const;

    /**
     * @brief find the smallest and largest value of the FE function
     *
     * If \p use_pk_nodal_fctn is set to false, it computes an approximation of
     * the largest and smallest value of the FE function by evaluating the FE
     * Function at around 900 points in triangles or around 1800 points in
     * quads.
     * If the argument is set to true, it uses the nodal functionals of
     * the Pk element with the corresponding order, i.e. point evaluations at
     * the degrees of freedom of the corresponding Pk element.
     *
     * @warning The case that \p use_pk_nodal_fctn is true works correctly only
     * if the reference transformation preserves values, e.g. for DG elements,
     * but not Hdiv elements.
     */
    void MinMax(double & min, double & max, const bool& use_pk_nodal_fctn
        = false) const;

    /**
     * @brief print a largest and a smallest value of this FE function depending
     * on the arguments given
     *
     * This function calls TFEFunction2D::MinMax and prints the given \p name
     * together with a minimum and maximum value. If the string is empty
     * (default) the used name will be TFEFunction2D::Name.
     *
     * If \p use_pk_nodal_fctn is set to false (default), it prints an
     * approximation of the largest and smallest value of the FE function. If
     * this parameter is set to true, it uses the nodal functionals of the Pk
     * element with the same order, i.e. only point evaluations at specific
     * points.
     *
     * @warning The case that \p use_pk_nodal_fctn is true works correctly only
     * if the reference transformation preserves values, e.g. for DG elements,
     * but not Hdiv elements.
     *
     * See also TFEFunction2D::MinMax.
     */
    void PrintMinMax(const std::string& name = "", const bool& use_pk_nodal_fctn
        = false) const;
   
   protected:
     
     /// @brief Routine used for the unification of the implementation
     /// of the routines Find*Gradient*. This is its only purpose.
     void getValuesFromOriginalOnes(double* retValues, int valuesOffset,
                                    int N_BaseFunct, int cellNumber,
                                    double* uorig, double* uxorig,
                                    double* uyorig, double* uzorig) const;
     
     /// @brief Routine used for the unification of the implementation
     /// of the routines Find*Gradient*. This is its only purpose.
     void getOrigValuesForCell(int cellNumber, double x, double y, double z,
                               double* uorig, double* uxorig, double* uyorig,
                               double* uzorig) const;
};

#endif
