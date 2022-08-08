#ifndef __FEVECTFUNCT3D__
#define __FEVECTFUNCT3D__

#include <FEFunction3D.h>

/** a function from a finite element space */
class TFEVectFunct3D : public TFEFunction3D
{
  protected:
    /** number of components */
    int N_Components;

  public:

    /// Default constructor. Constructs an empty object.
    TFEVectFunct3D();

    /** constructor with vector initialization */
    TFEVectFunct3D(std::shared_ptr<const TFESpace3D> fespace3D,
                   const std::string& name, double *values, int n_components);

    /// Copy assignment operator. Shallow copy, as the
    /// FEFunction does not take any memory responsibility.
    TFEVectFunct3D& operator=( const TFEVectFunct3D & );

    int GetN_Components() const
    { return N_Components; }
    
    /** return number of components */
    std::unique_ptr<TFEFunction3D> GetComponent(int i) const
    {
      // the name of the component will include the index i
      std::string fct_name(Name);
      fct_name += std::to_string(i);
      int length = FESpace3D->get_n_dof();
      auto ptr = std::unique_ptr<TFEFunction3D>(
        new TFEFunction3D(FESpace3D, fct_name, Values+i*length));
      return ptr;
    }

    /**
     * @brief compute the integral \int_{Gamma_i} u.n (flux) over a given surface
     *
     * In MPI case it returns the global flux, summed up over
     * the own domains of all processes.
     * To use this function it is necessary to mark the surface cell first 
     * (e.g., reading a .mesh file)
     * @todo give a list of boundary faces as input (instead of an int)
     *
     */
    double compute_flux(int surface_id) const;

    
    /** calculate errors to given vector function */
    void GetDeformationTensorErrors( 
        DoubleFunct3D *Exact, DoubleFunct3D *Exact1,
        DoubleFunct3D *Exact2,
        int N_Derivatives,
        MultiIndex3D *NeededDerivatives,
        int N_Errors, TFEFunction3D::ErrorMethod *ErrorMeth, 
        const CoeffFct3D& Coeff, TAuxParam3D *Aux,
        int n_fespaces, TFESpace3D **fespaces,
        double *errors);
    
    /// @brief calculate L2-norm of divergence error */
    /// @note this returns the correct value on all processes when using mpi.
    double GetL2NormDivergenceError(DoubleFunct3D *Exact_u1,
                                    DoubleFunct3D *Exact_u2,
                                    DoubleFunct3D *Exact_u3);

    /// @brief calculate L2-norm of divergence and (l^2-norm of) curl
    /// @note this returns the correct value on all processes when using mpi.
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
    /// derivatives of this fe vector function (u1,u2,u3), the order is: 
    /// x, y, z, u1, u2, u3, u1_x, u2_x, u3_x, u1_y, u2_y, u3_y, u1_z, u2_z, u3_z.
    /// (note: with dimension d: d*(d+2)
    /// Check the implementation of get_L2_norm_divergence_curl to see an
    /// example.
    /// @note this returns the correct value on all processes when using mpi.
    void get_functional_value(std::vector<double>& values,
                              const std::function<void(std::vector<double>&,
                                                       std::array<double, 15>)>&
                                                       functional) const;
    
        /** determine the value of function
        the given point */
    void FindValueLocal(const TBaseCell *cell, int cell_no, 
                        double x, double y, double z, 
                        double *values) const;

    /** @brief Determines the value of a vect function and its 
first derivatives at
      a given point; returns true if the gradient at (x,y,z) was found,
      and false if not (because, e.g., (x,y,z) belongs to another MPI process). */
    bool FindVectGradient(double x, double y, double z,
                          std::vector<double>& val1,
                          std::vector<double>& val2,
                          std::vector<double>& val3) const;
    
    /** @brief Determines the value of function and its first derivatives at
     * the given point lying WITHIN this cell (NOT on its boundary) */
    void FindVectGradientLocal(int cell_no, double x, double y, double z,
                               double* val1, double* val2, double* val3) const;
                          
    /** @brief Interpolates the old vector valued function to the new function.
     * Note that this is rather slow because no further information is 
     * required.
     */
    void Interpolate(TFEVectFunct3D *OldVectFunct);
   
};

#endif
