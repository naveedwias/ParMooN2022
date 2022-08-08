#ifndef INCLUDE_GENERAL_TEMPLATENAMES_H_
#define INCLUDE_GENERAL_TEMPLATENAMES_H_

class Example_CD2D;
class Example_CD3D;
class Example_Darcy2D;
class Example_Darcy3D;
class Example_NSE2D;
class Example_NSE3D;
class Example_TimeCD2D;
class Example_TimeCD3D;
class Example_TimeNSE2D;
class Example_TimeNSE3D;
class TFEFunction2D;
class TFEFunction3D;
class TFESpace2D;
class TFESpace3D;
class TFEVectFunct2D;
class TFEVectFunct3D;
class TSquareMatrix2D;
class TSquareMatrix3D;
class TMatrix2D;
class TMatrix3D;
enum class MultiIndex2D : int;
enum class MultiIndex3D : int;


#include "Enumerations_fe.h"
#include "Constants.h"
#include <vector>

/**
 * @brief this struct gives us the opportunity to write code which does not
 * depend on the space dimension.
 *
 * Basically, classes which are implemented separately for 2D and 3D can be 
 * used in the same code via 
 *     using FEFunction = typename Template_names<d>::FEFunction;
 * where `d` is either 2 or 3.
 * 
 */
template <int d>
struct Template_names
{
};
template <>
struct Template_names<2>
{
  typedef Example_CD2D Example_CD;
  typedef Example_Darcy2D Example_Darcy;
  typedef Example_NSE2D Example_NSE;
  typedef Example_TimeCD2D Example_TimeCD;
  typedef Example_TimeNSE2D Example_TimeNSE;
  typedef TFEFunction2D FEFunction;
  typedef TFESpace2D FESpace;
  typedef TFEVectFunct2D FEVectFunct;
  typedef std::vector<MultiIndex2D> MultiIndex_vector;
  typedef CoeffFct2D CoeffFct;
  typedef BoundCondFunct2D BoundaryConditionFunction;
  typedef BoundValueFunct2D BoundaryValuesFunction;
  typedef TSquareMatrix2D SquareMatrixD;
  typedef TMatrix2D MatrixD;
};
template <>
struct Template_names<3>
{
  typedef Example_CD3D Example_CD;
  typedef Example_Darcy3D Example_Darcy;
  typedef Example_NSE3D Example_NSE;
  typedef Example_TimeCD3D Example_TimeCD;
  typedef Example_TimeNSE3D Example_TimeNSE;
  typedef TFEFunction3D FEFunction;
  typedef TFESpace3D FESpace;
  typedef TFEVectFunct3D FEVectFunct;
  typedef std::vector<MultiIndex3D> MultiIndex_vector;
  typedef CoeffFct3D CoeffFct;
  typedef BoundCondFunct3D BoundaryConditionFunction;
  typedef BoundValueFunct3D BoundaryValuesFunction;
  typedef TSquareMatrix3D SquareMatrixD;
  typedef TMatrix3D MatrixD;
};

#endif // INCLUDE_GENERAL_TEMPLATENAMES_H_
