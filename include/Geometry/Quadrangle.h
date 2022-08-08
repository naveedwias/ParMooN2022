// =======================================================================
// @(#)Quadrangle.h        1.1 10/30/98
//
// Class:       TQuadrangle
// Purpose:     shape descriptor of a quadrangle
//
// Author:      Volker Behns  16.07.97
//
// =======================================================================

#ifndef __QUADRANGLE__
#define __QUADRANGLE__

#include <ShapeDesc.h>

/** shape descriptor of a quadrangle */
class TQuadrangle : public TShapeDesc
{
  public:
    // Constructor
    /** initialize all description paramters */
    TQuadrangle();

    // Methods
    /** return diameter of a cell */
    virtual double GetDiameter(TVertex **Verts) const override;

    /** return shortest edge of a cell */
    virtual double GetShortestEdge(TVertex **Verts) const override;

    /** return the length of the cell defined with the reference map */
    virtual double GetLengthWithReferenceMap(TVertex **Verts) const override;

    /** return measure of a cell */
    virtual double GetMeasure(TVertex **Verts) const override;

    /** check a special quadrangle whether it is a parallelogram
        or even a square */
    Shapes CheckQuad(const TVertex * const * Vertices) const;
};

#endif
