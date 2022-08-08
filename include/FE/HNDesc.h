#ifndef __HNDESC__
#define __HNDESC__

#include "Enumerations_fe.h"
#include <vector>

/** @brief describe fixed information for a hanging node */
class THNDesc
{
  protected:

    /** @brief type of hanging node */
    HNDesc type;

    /** @brief coefficient of other nodes in this coupling */
    std::vector<double> coeff;

  public:
    /** @brief constructor, filling all data */
    THNDesc(HNDesc type);

    // Methods
    /// @brief return the type of the hanging node descriptor
    HNDesc get_type() const
    { return type; }

    /** @brief return number of nodes in coupling */
    int GetN_Nodes() const
    { return coeff.size(); }

    /** @brief return coefficients of degrees of freedom in coupling */
    const double *GetCoeff() const
    { return coeff.data(); }

};
#endif
