#ifndef __HexaIsoparametric__
#define __HexaIsoparametric__

#include "RefTrans3D.h"
#include "BaseFunctions.h"
#include "QuadFormula.h"
#include <memory>

/** reference transformations for Hexahedron */
class THexaIsoparametric : public TRefTrans3D
{
  protected:
    /** x coordinate */
    double x0, x1, x2, x3, x4, x5, x6, x7;

    /** y coordinate */
    double y0, y1, y2, y3, y4, y5, y6, y7;

    /** z coordinate */
    double z0, z1, z2, z3, z4, z5, z6, z7;

    /** x parameters for reference transformation */
    double xc0, xc1, xc2, xc3, xc4, xc5, xc6, xc7;

    /** y parameters for reference transformation */
    double yc0, yc1, yc2, yc3, yc4, yc5, yc6, yc7;

    /** z parameters for reference transformation */
    double zc0, zc1, zc2, zc3, zc4, zc5, zc6, zc7;

    /** number of additional points */
    int N_AuxPoints;

    /** distance in x direction between real auxiliary point and 
        its position after a trilinear mapping */
    double XDistance[MaxN_BaseFunctions3D];

    /** distance in y direction between real auxiliary point and 
        its position after a trilinear mapping */
    double YDistance[MaxN_BaseFunctions3D];

    /** distance in z direction between real auxiliary point and 
        its position after a trilinear mapping */
    double ZDistance[MaxN_BaseFunctions3D];

    /** order of approximation */
    int ApproximationOrder;

    /** values of corresponding base function at quadpoints */
    double FctValues[MaxN_QuadPoints_3D][MaxN_BaseFunctions3D];

    /** xi-derivatives of corresponding base function at quadpoints */
    double XiDerValues[MaxN_QuadPoints_3D][MaxN_BaseFunctions3D];

    /** eta-derivatives of corresponding base function at quadpoints */
    double EtaDerValues[MaxN_QuadPoints_3D][MaxN_BaseFunctions3D];

    /** zeta-derivatives of corresponding base function at quadpoints */
    double ZetaDerValues[MaxN_QuadPoints_3D][MaxN_BaseFunctions3D];

    /** auxiliary array */
    double DoubleAux[MaxN_BaseFunctions3D];

    /** auxiliary array */
    int IntAux[MaxN_BaseFunctions3D];

    /** used quadrature rule */
    QuadratureFormula_type QuadFormula_type;
    std::unique_ptr<TQuadFormula> quadrature_formula;

    /** @brief Piola map, needed for vector values basis functions such as 
     *         Raviart-Thomas (RT) or Brezzi-Douglas-Marini (BDM).
     */
    void PiolaMapOrigFromRef(double xi, double eta, double zeta, int N_Functs,
                             const double *refD00, double *origD00) const;

  public:
    /** constuctor */
    THexaIsoparametric();

    /** transfer from reference element to original element */
    void GetOrigFromRef(double eta, double xi, double zeta,
                        double &x, double &y, double &z) const override;

    /** transfer a set of points form reference to original element */
    void GetOrigFromRef(int N_Points, const double *xi, const double *eta,
                        const double *zeta, double *x, double *y, double *z)
      const override;
    void GetOrigFromRef(const TQuadFormula& qf_ref, TQuadFormula& qf_orig) const
      override;

    /** transfer from original element to reference element */
    void GetRefFromOrig(double x, double y, double z,
                        double &eta, double &xi, double &zeta) const override;
    
    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(double xi, double eta, double zeta, int N_BaseFunct,
                const double *uref, const double *uxiref, const double *uetaref,
                const double *uzetaref,
                double *uorig, double *uxorig, double *uyorig, double *uzorig,
                int BaseVectDim = 1) const override;

    void GetOrigAllDerivatives(double xi, double eta, double zeta,
                               int N_BaseFunct,
                               const double *refD000, const double *refD100,
                               const double *refD010, const double *refD001,
                               const double *refD200, const double *refD110,
                               const double *refD101, const double *refD020,
                               const double *refD011, const double *refD002,
                               double *origD000, double *origD100,
                               double *origD010, double *origD001,
                               double *origD200, double *origD110,
                               double *origD101, double *origD020,
                               double *origD011, double *origD002,
                               int BaseVectDim = 1) const override;

    /** set element to cell */
    void SetCell(const TBaseCell * cell) override;

    /** set order of approximation */
    void SetApproximationOrder(int order)
    {
      if(order <= 0)
        ApproximationOrder = 1;
       else
        ApproximationOrder = order;
    }

    /** set used quadrature formula */
    void SetQuadFormula(QuadratureFormula_type formula);
};

#endif
