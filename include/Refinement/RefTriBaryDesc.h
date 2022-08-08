#ifndef __REFTRIBARYDESC__
#define __REFTRIBARYDESC__

#include <RefDesc.h>

/** refinement descriptor for barycentric refinement of a triangle */
class TRefTriBaryDesc : public TRefDesc
{
  public:
    explicit TRefTriBaryDesc(const TShapeDesc *shape);
};

#endif // __REFTRIBARYDESC__
