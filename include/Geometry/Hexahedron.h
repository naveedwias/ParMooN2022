// =======================================================================
// @(#)Hexahedron.h        1.2 10/18/99
//
// Class:       THexahedron
// Purpose:     shape descriptor of a hexahedron
//
// Author:      Volker Behns  16.07.97
//
// =======================================================================

#ifndef __HEXAHEDRON__
#define __HEXAHEDRON__

#include <ShapeDesc.h>

/** shape descriptor of a hexahedron */
class THexahedron : public TShapeDesc
{
  public:
    // Constructor
    THexahedron();

    // Methods
    /** return diameter of a cell */
    virtual double GetDiameter(TVertex **Verts) const override;

    /** return shortest edge of a cell */
    virtual double GetShortestEdge(TVertex **Verts) const override;

    /** return the length of the cell defined with the reference map */
    virtual double GetLengthWithReferenceMap(TVertex **Verts) const override;

    /** return measure of a cell */
    virtual double GetMeasure(TVertex **Verts) const override;

    /** check whether this hexahedron is a brick */
    Shapes CheckHexa(const TVertex * const * Vertices) const;
};

#endif
