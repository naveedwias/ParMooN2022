#ifndef __REFTRANS3D__
#define __REFTRANS3D__

#include <vector>

// forward declaration
class TBaseCell;
class TQuadFormula;

/** @brief reference transformations for 3D geometric objects */
class TRefTrans3D
{
  protected:
    const TBaseCell *Cell;

  public:
    /** @brief constuctor */
    TRefTrans3D() = default;
    
    virtual ~TRefTrans3D() = default;

    /** @brief transfer form reference element to original element */
    virtual void GetOrigFromRef(double xi, double eta, double zeta,
                        double &x, double &y, double &z) const = 0;
    
    virtual void GetOrigFromRef(int N_Points, const double *eta,
                                const double *xi, const double *zeta,
                                double *x, double *y, double *z) const = 0;
    
    /** @brief transfer form reference element to original element */
    virtual void GetOrigFromRef(const TQuadFormula& qf_ref,
                                TQuadFormula& qf_orig) const = 0;

    /** @brief transfer from original element to reference element */
    virtual void GetRefFromOrig(double x, double y, double z,
                                double &xi, double &eta, double &zeta) const =0;
    
    /** calculate functions and derivatives from reference element
        to original element */
    virtual void GetOrigValues(double xi, double eta, double zeta,
                               int N_BaseFunct, const double *uref,
                               const double *uxiref, const double *uetaref,
                               const double *uzetaref, double *uorig,
                               double *uxorig, double *uyorig, double *uzorig,
                               int BaseVectDim = 1) const = 0;

    /** @brief transform functions and derivatives up to second order from
     * reference element to original element */
    virtual void GetOrigAllDerivatives(double xi, double eta, double zeta,
                                       int N_BaseFunct,
                                       const double *refD000,
                                       const double *refD100,
                                       const double *refD010,
                                       const double *refD001,
                                       const double *refD200,
                                       const double *refD110,
                                       const double *refD101,
                                       const double *refD020,
                                       const double *refD011,
                                       const double *refD002,
                                       double *origD000, double *origD100,
                                       double *origD010, double *origD001,
                                       double *origD200, double *origD110,
                                       double *origD101, double *origD020,
                                       double *origD011, double *origD002,
                                       int BaseVectDim = 1) const = 0;

    /** @brief set original element to cell */
    virtual void SetCell(const TBaseCell *cell)
    {  Cell = cell; }

    /// @brief column-major Jacobian of the reference transformation
    virtual void GetTransformationDerivatives(double xi, double eta,
      double zeta, double* matrix) const;

    /// @brief column-major inverse Jacobian of the reference transformation
    virtual void GetInverseTransformationDerivatives(double xi, double eta,
      double zeta, double* matrix) const;
};

#endif
