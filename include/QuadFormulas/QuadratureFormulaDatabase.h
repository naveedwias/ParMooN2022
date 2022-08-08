#ifndef QUADRATUREFORMULADATABASE_h
#define QUADRATUREFORMULADATABASE_h

#include "Enumerations_fe.h"
#include "Enumerations_geometry.h"
#include "Enumerations_quadrature_formula.h"

class TQuadFormula;

namespace QuadratureFormulaDatabase
{
  /** @brief initialize the cache
   *
   * @note call this before any of the other methods in this namespace
   */
  void create();
  
  /** @brief quadrature formula on given reference element */
  const TQuadFormula *qf_from_degree(int accuracy, BFRefElements RefElem);

  /** @brief quadrature formula on given shape */
  const TQuadFormula *qf_from_degree(int accuracy, Shapes shape);
  
  /** @brief free memory, delete cached quadrature formulas */
  void destroy();
}

#endif // QUADRATUREFORMULADATABASE_h
