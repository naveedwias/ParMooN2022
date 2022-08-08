#ifndef __REFTRANS2D__
#define __REFTRANS2D__

#include <vector>

// forward declaration
class TBaseCell;
class TQuadFormula;

/** @brief reference transformations for 2D geometric objects.
 * 
 * This is the abstract base class, all members methods are implemented in the
 * derived classes.
 */
class TRefTrans2D
{
  protected:
    /** @brief the cell to which this reference transformation maps. This is
     * often referred to as the 'original' cell. */
    const TBaseCell *Cell;

  public:
    /** @brief constuctor, you need to call `SetCell` to use this object */
    TRefTrans2D() = default;
    
    virtual ~TRefTrans2D() = default;

    /** @brief transfer a point \f$(\xi,\eta)\f$ form reference element to
     * original element */
    virtual void GetOrigFromRef(double xi, double eta, double &x, double &y)
      const = 0;
    
    /** @brief transfer a set of points \f$(\xi_i,\eta_i)\f$ form reference to 
     * original element */
    virtual void GetOrigFromRef(int N_Points, const double *xi,
                                const double *eta, double *x, double *y) const 
      = 0;
    
    /** @brief transfer a quadrature formula `qf_ref` on the reference cell into
     * a quadrature formula on the original cell. This includes the quadrature
     * weights as well as the quadrature points.
     */
    virtual void GetOrigFromRef(const TQuadFormula& qf_ref,
                                TQuadFormula& qf_orig) const = 0;

    /** @brief transfer a point \f$(x,y)\f$ from original element to reference 
     * element. This is the inverse of `GetOrigFromRef(xi, eta, x, y)` */
    virtual void GetRefFromOrig(double x, double y, double &xi, double &eta)
      const = 0;

    /** @brief transform functions and derivatives up to first order from
     * reference element to original element */
    virtual void GetOrigValues(double xi, double eta, int N_BaseFunct,
                               const double *uref, const double *uxiref,
                               const double *uetaref, double *uorig,
                               double *uxorig, double *uyorig,
                               int BaseVectDim = 1) const = 0;

    /** @brief transform functions and derivatives up to second order from
     * reference element to original element */
    virtual void GetOrigAllDerivatives(double xi, double eta, int N_BaseFunct,
                                       const double *refD00,
                                       const double *refD10,
                                       const double *refD01,
                                       const double *refD20,
                                       const double *refD11,
                                       const double *refD02,
                                       double *origD00, double *origD10,
                                       double *origD01, double *origD20,
                                       double *origD11, double *origD02,
                                       int BaseVectDim = 1) const = 0;

    /** @brief transform functions and derivatives up to first order from a
     * given joint (edge) in the reference element to the corresponding edge in
     * the original element */
    virtual void GetOrigValues(int joint, double zeta, int N_BaseFunct,
                               const double *uref, const double *uxiref,
                               const double *uetaref, double *uorig,
                               double *uxorig, double *uyorig,
                               int BaseVectDim = 1) const = 0;

    /** @brief set up this transformation to map from the reference cell to
     * the given cell. */
    virtual void SetCell(const TBaseCell *cell) = 0;

    /// @brief column-major Jacobian of the reference transformation
    virtual void GetTransformationDerivatives(double xi, double eta,
      double* matrix) const;

    /// @brief column-major inverse Jacobian of the reference transformation
    virtual void GetInverseTransformationDerivatives(double xi, double eta,
      double* matrix) const;
};

#endif
