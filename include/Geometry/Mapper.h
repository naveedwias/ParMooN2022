#ifndef __MAPPER__
#define __MAPPER__

#include "Enumerations_geometry.h"

/** mapper for faces of 3D geometric objects */
class TMapper
{
  protected:
    Mapper Type;

    /** mapping of original vertices */
    const int *MapOrigVerts;
    /** mapping of original edges */
    const int *MapOrigEdges;

    /** mapping of refined vertices */
    const int *MapRefVerts;
    /** mapping of refined edges */
    const int *MapRefEdges;
    /** mapping of refined faces */
    const int *MapRefFaces;

  public:
    //Constructor
    explicit TMapper(Mapper which);

    //Methods
    /** return type of mapper */
    Mapper GetType() const
    { return Type; }

    /** return mapper of vertices and faces */
    void GetMapperRef(const int *&MapVerts, const int *&MapFaces) const
    {
      MapVerts = MapRefVerts;
      MapFaces = MapRefFaces;
    }

    /** return mapper of original vertices and faces */
    void GetMapperOrig(const int *&MapVerts, const int *&MapEdges) const
    {
      MapVerts = MapOrigVerts;
      MapEdges = MapOrigEdges;
    }
};

#endif
