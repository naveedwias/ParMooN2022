// =======================================================================
// @(#)Tetrahedron.h        1.2 10/18/99
//
// Class:       TTetrahedron
// Purpose:     shape descriptor of a tetrahedron
//
// Author:      Volker Behns  16.07.97
//
// =======================================================================

#ifndef __TETRAHEDRON__
#define __TETRAHEDRON__

#include <ShapeDesc.h>

/** shape descriptor of a tetrahedron */
class TTetrahedron : public TShapeDesc
{
  public:
    // Constructor
    /** build the shape descriptor for a tetrahedron */
    TTetrahedron();

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
