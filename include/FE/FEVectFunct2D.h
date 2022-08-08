#ifndef __FEVECTFUNCT2D__
#define __FEVECTFUNCT2D__

#include <FEFunction2D.h>

/** a function from a finite element space */
class TFEVectFunct2D : public TFEFunction2D
{
  protected:
    /** number of components */
    int N_Components;

  public:

    /// Default constructor. Constructs an empty object.
    TFEVectFunct2D();

    /** constructor with vector initialization */
    TFEVectFunct2D(std::shared_ptr<const TFESpace2D> fespace2D,
                   const std::string& name, double *values,
                   int n_components);

    /// Copy assignment operator. Shallow copy, as the
    /// FEFunction does not take any memory responsibility.
    TFEVectFunct2D& operator=( const TFEVectFunct2D & );

    /** return number of components */
    int GetN_Components() const
    { return N_Components; }

    /** return i-th component as FEFunction2D */
    std::unique_ptr<TFEFunction2D> GetComponent(int i) const
    {
      // the name of the component will include the index i
      std::string fct_name(Name);
      fct_name += std::to_string(i);
      int length = FESpace2D->get_n_dof();
      auto ptr = std::unique_ptr<TFEFunction2D>(
        new TFEFunction2D(FESpace2D, fct_name, Values+i*length));
      return ptr;
    }

    /** calculate errors to given vector function */
    void GetDeformationTensorErrors( 
        DoubleFunct2D *Exact, DoubleFunct2D *Exact1,
        int N_Derivatives,
        MultiIndex2D *NeededDerivatives,
        int N_Errors, TFEFunction2D::ErrorMethod *ErrorMeth, 
        const CoeffFct2D& Coeff, TAuxParam2D *Aux,
        int n_fespaces, TFESpace2D **fespaces,
        double *errors);

    /// @brief calculate L2-norm of divergence and curl
    std::pair<double,double> get_L2_norm_divergence_curl() const;
    
    /// @brief compute the integral of (almost) arbitrary functionals for this
    /// fe vector function.
    /// the length of the std::vector 'values' determines the number of 
    /// functional values you wish to compute. You can compute functionals of
    /// the form \f$ \| S(u) \|_{L^2(\Omega)} \f$ with \f$S\f$ mapping the 
    /// vector \f$u\f$ to some scalar function whose \f$L^2\f$-norm is computed.
    /// The provided 'functional' represents \f$S\f$ and computes the local 
    /// contributions for each quadrature point. Its arguments are a 
    /// std::vector with the same size as 'values' and a std::array which 
    /// consists of the coordinates of the quadrature points and the values and
    /// derivatives of this fe vector function (u1,u2), the order is: 
    /// x, y, u1, u2, u1_x, u2_x, u1_y, u2_y.
    /// Check the implementation of TFEVectFunct2D::get_L2_norm_divergence_curl
    /// to see an example.
    void get_functional_value(std::vector<double>& values,
                              const std::function<void(std::vector<double>&, std::array<double, 8>)>& functional) const;

    /** calculate L2-norm of divergence error */
    double GetL2NormDivergenceError(DoubleFunct2D *Exact_u1, DoubleFunct2D *Exact_u2);

    /** calculate L2-norm of (u \cdot n)-error at the boundary */
    double GetL2NormNormalComponentError(BoundValueFunct2D *Exact_u1, BoundValueFunct2D *Exact_u2, bool rescale_by_h_E = false);

    /** calculate L2-norm of (u \cdot n)-error at the boundary */
    double GetL2NormNormalComponentError(BoundValueFunct2D *Exact_u1, BoundValueFunct2D *Exact_u2,  int boundary_component_id, bool rescale_by_h_E = false);
    

    /** calculate the integral int_{Gamma_id} u.n (flux) */
    double compute_flux(int _id) const;

    /** @brief Determines the value of vector function and its first derivatives at
     * the given point */
    void FindVectGradient(double x, double y, double *val1, double *val2) const;

    /** @brief Determines the value of vector function and its first derivatives at
     * the given point lying WITHIN this cell (NOT on its boundary) */
    void FindVectGradientLocal(int cell_no, double x, double y,
                               double *val1, double *val2) const;

    /** @brief Interpolates the old vector valued function to the new function.
     * Note that this is rather slow because no further information is 
     * required.
     */
    void Interpolate(TFEVectFunct2D *OldVectFunct);
    
    /** determine the value of function at
    the given point componentwise */
    void FindValueLocal(const TBaseCell *cell, int cell_no, double x, double y, 
        double *values) const;

    /** @brief multiply function with a scalar alpha. Only non-Dirichlet dofs are 
     *         multiplied! */
    TFEVectFunct2D& operator*=(double alpha);

    /** @brief add one TFEVectFunct2D to another one. Both have to be defined on
     *         the same space. Only non-Dirichlet dofs are added!  */
    TFEVectFunct2D & operator+=(const TFEVectFunct2D & rhs);
};

#endif

// #ifdef __2D__
// # endif // #ifdef __2D__
