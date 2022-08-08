// =======================================================================
// @(#)Triangle.h        1.1 10/30/98
//
// Class:       TTriangle
// Purpose:     shape descriptor of a triangle
//
// Author:      Volker Behns  16.07.97
//
// =======================================================================

#ifndef __TRIANGLE__
#define __TRIANGLE__

#include <ShapeDesc.h>

/** shape descriptor of a triangle */
class TTriangle : public TShapeDesc
{
  public:
    // Constructor
    /** build the shape descriptor for a triangle */
    TTriangle();

    // Methods
    /** return diameter of a cell */
    virtual double GetDiameter(TVertex **Verts) const override;

    /** return shortest edge of a cell */
    virtual double GetShortestEdge(TVertex **Verts) const override;

    /** return the length of the cell defined with the reference map */
    virtual double GetLengthWithReferenceMap(TVertex **Verts) const override;

    /** return measure of a cell */
    virtual double GetMeasure(TVertex **Verts) const override;
};

#endif
