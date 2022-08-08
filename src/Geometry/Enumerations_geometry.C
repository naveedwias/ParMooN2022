#include "Enumerations_geometry.h"
#include "MooNMD_Io.h"

std::ostream &operator<<(std::ostream &os, RefinementStrategyType type)
{
  switch (type)
  {
    case RefinementStrategyType::MaximalEstimatedLocalError:
      os << "Maximal estimated local error";
      break;
    case RefinementStrategyType::PortionOfEstimatedGlobalError:
      os << "Prescribed portion of estimated global error";
      break;
    case RefinementStrategyType::EquilibrationStrategy:
      os << "Equilibration strategy";
      break;
    default:
      ErrThrow("unknown RefinementStrategyType");
  }
  int typeIndex = int(type);
  os << " (database value = " << typeIndex << ")";
  return os;
}
