#ifndef __TETRAAFFIN__
#define __TETRAAFFIN__

#include "RefTrans3D.h"

/** reference transformations for tetrahedron */
class TTetraAffin : public TRefTrans3D
{
  protected:
    /** x coordinate */
    double x0, x1, x2, x3;

    /** y coordinate */
    double y0, y1, y2, y3;

    /** z coordinate */
     double z0, z1, z2, z3;   

    /** x parameters for reference transformation */
    double xc0, xc1, xc2, xc3;

    /** y parameters for reference transformation */
    double yc0, yc1, yc2, yc3;

    /** z parameters for reference transformation */
    double zc0, zc1, zc2, zc3;

    /** detjk */
    double detjk;
    
    /** @brief transformation matrix, needed to transform second derivatives */
    mutable double Eye[6][6];
    /** @brief flag to indicate if Eye is filled with correct values */
    mutable bool ready_for_transforming_2nd_derivatives;
    
    /** @brief Piola transformation for vector basis */
    void PiolaMapOrigFromRef(double xi, double eta, double zeta, int N_Functs,
                             const double *refD000, const double *refD100,
                             const double *refD010, const double *refD001,
                             double *origD000, double *origD100,
                             double *origD010, double *origD001) const;

  public:
    /** constuctor */
    TTetraAffin();
    
    /** transfer form reference element to original element */
    void GetOrigFromRef(double xi, double eta, double zeta, double &x,
                        double &y, double &z) const override;

    /** transfer a set of points form reference to original element */
    void GetOrigFromRef(int N_Points, const double *eta, const double *xi,
                        const double *zeta, double *x, double *y, double *z)
      const override;
    void GetOrigFromRef(const TQuadFormula& qf_ref, TQuadFormula& qf_orig) const
      override;

    /** transfer from original element to reference element */
    void GetRefFromOrig(double x, double y, double z, double &eta, double &xi,
                        double &zeta) const override;

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(double xi, double eta, double zeta, int N_BaseFunct,
                       const double *uref, const double *uxiref,
                       const double *uetaref,  const double *zetaref,
                       double *uorig, double *uxorig, double *uyorig, 
                       double *uzorig, int _BaseVectDim = 1) const override;

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

    virtual void GetTransformationDerivatives(double xi, double eta,
        double zeta, double* matrix) const override;
};

#endif
