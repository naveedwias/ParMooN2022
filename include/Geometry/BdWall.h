// =======================================================================
// @(#)BdWall.h        1.1 07/16/99
//
// Class:       TBdWall
// Superclass:  TBoundComp
// Purpose:     a Wall as a component of a boundary part
//
// Author:      Volker Behns  05.07.99
//
// =======================================================================

#ifndef __BDWALL__
#define __BDWALL__

#include <BoundComp2D.h>
#include <BoundComp3D.h>

/** a Wall as a component of a boundary part */
class TBdWall : public TBoundComp3D
{
  protected:
    /** 2d boundary component */
    TBoundComp2D * BdComp2D;

    /** drift in x-direction */
    double DriftX;
    /** drift in y-direction */
    double DriftY;
    /** drift in z-direction */
    double DriftZ;
    
    /** conic scale (in x- and y-direction) */
    double ConicScale;

  public:
    // Constructor
    TBdWall(int id, TBoundComp2D *bdcomp2d);
    
    // Methods
    int SetParams(double drx, double dry, double drz, double conic_sc);

    /** return the coordinates of parameter value T, S */
    virtual int GetXYZofTS(double T, double S,
                           double &X, double &Y, double &Z) const override;

    /** return the parameter value T, S of coordinates */
    virtual int GetTSofXYZ(double X, double Y, double Z,
                           double &T, double &S) const override;

    /** return parameters and coordinates of a given linear
        combination of vertices */
    virtual int GetXYZandTS(int N_Points, double *LinComb,
                            double *xp, double *yp, double *zp,
                            double *tp, double *sp,
                            double &X, double &Y, double &Z,
                            double &T, double &S) const override;

    /** read parameter from input file */
    virtual int ReadIn(std::istream &dat) override;

    virtual void get_normal_vector(double, double, double, double&, double&,
                                   double&) const override
    {
      ErrThrow(" ** ERROR: get_normal_vector() not yet implemented for BdWall ");
    };
    
    /** return BdComp2D */
    TBoundComp2D* GetBdComp2D() const
    { return BdComp2D; }
};

#endif
