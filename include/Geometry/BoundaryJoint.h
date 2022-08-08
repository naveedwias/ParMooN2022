#ifndef BOUNDARYJOINT_H
#define BOUNDARYJOINT_H

#include "Joint.h"

class BoundaryJoint : public TJoint
{
  public:
    /** @brief return whether this is an interior joint */
    virtual bool InnerJoint() const
    { return false; }

    /** @brief check whether the refinement pattern on both side patch,
        dummy here: there is no neighbour */
    virtual int CheckMatchingRef(TBaseCell *, int, struct StoreGeom &Tmp)
    {
      Tmp.Filled = false;
      return 0;
    }

  protected:
    /// @brief protected constructor, only used from derived classes
    BoundaryJoint() = default;
};

#endif // BOUNDARYJOINT_H
