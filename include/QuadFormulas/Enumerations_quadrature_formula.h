#ifndef ENUMERATIONS_QUADRATURE_FORMULA_H
#define ENUMERATIONS_QUADRATURE_FORMULA_H

#include <iosfwd>

constexpr int N_QuadratureFormula_type = 58;
enum class QuadratureFormula_type
{
  Gauss1Line,
  Gauss2Line,
  Gauss3Line,
  Gauss4Line,
  Gauss5Line,
  Gauss6Line,
  Gauss7Line,
  Gauss8Line,
  Gauss9Line,
  Gauss10Line,
  Gauss11Line,
  Gauss12Line,
  Gauss2W1Line,
  Gauss4W1Line,
  Gauss6W1Line,
  Gauss8W1Line,
  Gauss16W2Line,
  /////////////////////////////////////////////////////////////////////////////
  BaryCenterTria,
  MidPointTria,
  SevenPointTria,
  Gauss3Tria,
  VertexTria,
  Degree8Tria,
  Degree9Tria,
  Degree11Tria,
  Degree19Tria,
  Gauss2Quad,
  Gauss3Quad,
  Gauss4Quad,
  Gauss5Quad,
  Gauss6Quad,
  Gauss7Quad,
  Gauss8Quad,
  Gauss9Quad,
  VertexQuad,
  SimpsonQuad,
  CompGauss3Tria,
  CompGauss4Tria,
  Gauss_Degree8Tria,
  /////////////////////////////////////////////////////////////////////////////
  BaryCenterTetra,
  VertexTetra,
  P2Tetra,
  P4Tetra,
  P5Tetra,
  P8Tetra,
  VertexHexa,
  Gauss2Hexa,
  Gauss3Hexa,
  Gauss4Hexa,
  Gauss5Hexa,
  Gauss6Hexa,
  Gauss7Hexa,
  Gauss8Hexa,
  Gauss9Hexa,
  VerticesAndOrigin,
  VerticesAndOrigin15,
  VerticesAndOrigin57,
  Degree7_Points38
};

std::ostream& operator<<(std::ostream& out, const QuadratureFormula_type t);

#endif // ENUMERATIONS_QUADRATURE_FORMULA_H
