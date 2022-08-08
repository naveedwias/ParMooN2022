#ifndef __QUADBilinear__
#define __QUADBilinear__

#include "RefTrans2D.h"

/** @brief bilinear reference transformation for quadrangles */
class TQuadBilinear : public TRefTrans2D
{
  protected:
    /** @brief x coordinates of vertices */
    double x0, x1, x2, x3;

    /** @brief y coordinate of vertices */
    double y0, y1, y2, y3;

    /** @brief x parameters for reference transformation */
    double xc0, xc1, xc2, xc3;

    /** @brief y parameters for reference transformation */
    double yc0, yc1, yc2, yc3;
    
    /** @brief transformation matrix, needed to transform second derivatives */
    mutable double Eye[5][5];
    /** @brief flag to indicate if Eye is filled with correct values */
    mutable bool ready_for_transforming_2nd_derivatives;

    /** @brief Piola transformation for vector valued basis functions */
    void PiolaMapOrigFromRef(double xi, double eta, int N_Functs,
                             const double *refD00, double *origD00) const;
    /** @brief Piola transformation for the derivatives of vector valued basis 
     *         functions */
    void PiolaMapOrigFromRef(double xi, double eta, int N_Functs,
                             const double *refD00, const double *refD10,
                             const double *refD01,
                             double *origD10, double *origD01) const;

  public:
    /** @brief constuctor, you need to call `SetCell` to use this object */
    TQuadBilinear();

    /** @brief transfer a point \f$(\xi,\eta)\f$ form reference element to
     * original element */
    void GetOrigFromRef(double xi, double eta, double &x, double &y) const
      override final;

    /** @brief transfer a set of points \f$(\xi_i,\eta_i)\f$ form reference to
     * original element */
    void GetOrigFromRef(int N_Points, const double *xi, const double *eta,
                        double *x, double *y) const override final;

     /** @brief transfer a quadrature formula `qf_ref` on the reference cell into
     * a quadrature formula on the original cell. This includes the quadrature
     * weights as well as the quadrature points.
     */
    void GetOrigFromRef(const TQuadFormula& qf_ref, TQuadFormula& qf_orig) const
      override final;

    /** @brief transfer a point \f$(x,y)\f$ from original element to reference
     * element. This is the inverse of `GetOrigFromRef(xi, eta, x, y)` */
    void GetRefFromOrig(double x, double y, double &xi, double &eta) const
      override final;

    /** @brief transform functions and derivatives up to first order from
     * reference element to original element */
    void GetOrigValues(double xi, double eta, int N_BaseFunct,
                       const double *uref, const double *uxiref,
                       const double *uetaref,
                       double *uorig, double *uxorig, double *uyorig,
                       int _BaseVectDim = 1) const override final;

    /** @brief transform functions and derivatives up to second order from
     * reference element to original element */
    void GetOrigAllDerivatives(double xi, double eta, int N_BaseFunct,
                               const double *refD00, const double *refD10,
                               const double *refD01, const double *refD20,
                               const double *refD11, const double *refD02,
                               double *origD00, double *origD10,
                               double *origD01, double *origD20,
                               double *origD11, double *origD02,
                               int BaseVectDim = 1) const override final;

    /** @brief transform functions and derivatives up to first order from a
     * given joint (edge) in the reference element to the corresponding edge in
     * the original element */
    void GetOrigValues(int joint, double zeta, int N_BaseFunct,
                       const double *uref, const double *uxiref,
                       const double *uetaref,
                       double *uorig, double *uxorig, double *uyorig,
                       int _BaseVectDim = 1) const override final;

    /** @brief set up this transformation to map from the reference cell to
     * the given cell. */
    void SetCell(const TBaseCell * cell) override final;
};

#endif
