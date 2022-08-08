// =======================================================================
// @(#)BdPolygon.h        1.1 08/12/99
//
// Class:       TBdPolygon
// Superclass:  TBoundComp
// Purpose:     component is polygon
//
// Author:      Gunar Matthies 04.08.1999
//
// =======================================================================

#ifndef __BdPolygon__
#define __BdPolygon__

#include <BoundComp2D.h>

/** free boundary */
class TBdPolygon : public TBoundComp2D
{
  protected:
    /** number of points */
    int N_Points;
    /** array of coordinates of points  */
    double *Coords;

  public:
    // Constuctor
    /** constructor initializes the parameter array */
    TBdPolygon (int id, int n_points);

    // Methods
    /** set all parameters */
    void SetParams (int n_points, double *coords);

    /** get number of vertices */
    int GetN_Points();

    /** return the coordinates of parameter value T */
    virtual int GetXYofT(double T, double &X, double &Y) const override;

    /** return the parameter value T of coordinates (X, Y) */
    virtual int GetTofXY(double X, double Y, double &T) const override;

    /** read parameter from input file */
    virtual int ReadIn(std::istream &dat) override;

    /** get number of initial vertices on this component */
    virtual int GetN_InitVerts() override
    { return 2; }

    virtual int GenInitVerts(double *&, int, int *&, int) override
    { return 0; }
};

#endif
