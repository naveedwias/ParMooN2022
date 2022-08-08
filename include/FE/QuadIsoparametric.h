#ifndef __QUADISOPARAMETRIC__
#define __QUADISOPARAMETRIC__

#include "RefTrans2D.h"
#include "QuadFormula.h"
#include "BaseFunctions.h"
#include <memory>

/** @brief isoparametric reference transformation for quadrangles */
class TQuadIsoparametric : public TRefTrans2D
{
  protected:
    /** @brief x coordinates of vertices */
    double x[4];

    /** @brief y coordinate of vertices */
    double y[4];

    /** @brief x parameters for reference transformation */
    double xc0, xc1, xc2, xc3;

    /** @brief y parameters for reference transformation */
    double yc0, yc1, yc2, yc3;

    /** @brief number of additional points */
    int N_AuxPoints;

    /** @brief distance in x direction between real auxiliary point and 
        its position after a affine mapping */
    double XDistance[MaxN_BaseFunctions2D];

    /** @brief distance in y direction between real auxiliary point and 
        its position after a affine mapping */
    double YDistance[MaxN_BaseFunctions2D];

    /** @brief order of approximation */
    int ApproximationOrder;

    /** @brief values of corresponding base function at quadpoints */
    double FctValues[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];

    /** @brief xi-derivatives of corresponding base function at quadpoints */
    double XiDerValues[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];

    /** @brief eta-derivatives of corresponding base function at quadpoints */
    double EtaDerValues[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];

    /** @brief auxiliary array */
    double DoubleAux[MaxN_BaseFunctions2D];

    /** @brief auxiliary array */
    int IntAux[MaxN_BaseFunctions2D];

    /** @brief used quadrature rule */
    std::unique_ptr<TQuadFormula> quadrature_formula;

    /** @brief Piola map, needed for vector values basis functions such as 
     *         Raviart-Thomas (RT) or Brezzi-Douglas-Marini (BDM).
     */
    void PiolaMapOrigFromRef(double xi, double eta, int N_Functs,
                             const double *refD00, double *origD00) const;
    /** @brief Piola transformation for the derivatives of vector valued basis 
     *         functions */
    void PiolaMapOrigFromRef(double xi, double eta, int N_Functs,
                             const double *refD00, const double *refD10,
                             const double *refD01, double *origD10,
                             double *origD01) const;

  public:
    /** @brief constuctor, you need to call `SetCell` to use this object */
    TQuadIsoparametric();

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
                       int BaseVectDim = 1) const override final;

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
                       int BaseVectDim = 1) const override final;

    /** @brief set up this transformation to map from the reference cell to
     * the given cell. */
    void SetCell(const TBaseCell * cell) override final;

    /** @brief set order of approximation */
    void SetApproximationOrder(int order)
    {
      if(order <= 0)
        ApproximationOrder = 1;
       else
        ApproximationOrder = order;
    }

    /** @brief set used quadrature formula */
    void SetQuadFormula(QuadratureFormula_type formula);
};

#endif
