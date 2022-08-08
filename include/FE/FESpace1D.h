#ifndef __FESPACE1D__
#define __FESPACE1D__

#include "FESpace.h"
#include "Enumerations_fe.h"
#include "FiniteElement.h"
#include <vector>
#include <map>

/** class for all 1D finite element spaces */
class TFESpace1D : public TFESpace
{
  protected:

    /** @brief construct space */
    void ConstructSpace();
    
  public:
    /** @brief constructor for building a space with elements of order k */
    TFESpace1D(const TCollection *coll, const std::string& name, int k);

    /** @brief constructor for building a space with the given elements */
    TFESpace1D(const TCollection *coll, const std::string& name, FE_type *fes);

    ~TFESpace1D() = default;
};

#endif // __FESPACE1D__
