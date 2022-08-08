// =======================================================================
// @(#)BdCircle.h        1.2 07/16/99
//
// Class:       TBdCircle
// Superclass:  TBoundComp
// Purpose:     a part of a circle as a component of a boundary part
//
// Author:      Volker Behns  18.06.97
//
// =======================================================================

#ifndef __BDCIRCLE__
#define __BDCIRCLE__

#include <BoundComp2D.h>

/** a part of a circle as a component of a boundary part */
class TBdCircle : public TBoundComp2D
{
  protected:
    /** x coordinate of midpoint */
    double Xmid;
    /** y coordinate of midpoint */
    double Ymid;
    /** radii of the arc */
    double Radius_a, Radius_b;
    /** begin angle */
    double Phi1;
    /** end angle */
    double Phi2;

  public:
    // Constructor
    explicit TBdCircle(int id);

    // Methods
    /** set all parameters to the given values */
    void SetParams (double xmid, double ymid, double radius_a, 
                    double radius_b, double phi1, double phi2);

    /** return the coordinates of parameter value T */
    virtual int GetXYofT(double T, double &X, double &Y) const override;

    /** return the parameter value T of coordinates (X, Y) */
    virtual int GetTofXY(double X, double Y, double &T) const override;

    /** read parameter from input stream */
    virtual int ReadIn(std::istream &dat) override;

    /** get number of initial vertices on this component */
    virtual int GetN_InitVerts() override;
    virtual int GenInitVerts(double *&points, int I_points,
                             int *&edges, int I_edges) override;

    bool is_full_circle() const;

  protected:
    int GetN_InitVertsSub(double Phi_a, double Phi_b, int Level);
    int GenInitVertsSub(double Phi_a, double Phi_b, int Level,
                        double *&points, int &I_points,
                        int *&edges, int &I_edges);
};

#endif
