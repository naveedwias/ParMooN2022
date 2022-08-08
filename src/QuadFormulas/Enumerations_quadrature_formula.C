#include "Enumerations_quadrature_formula.h"
#include "MooNMD_Io.h"

std::ostream& operator<<(std::ostream& out, const QuadratureFormula_type t)
{
  const char* s = 0;
#define PROCESS_VAL(p)                                                         \
  case(QuadratureFormula_type::p):                                             \
    s = #p;                                                                    \
    break;
  switch(t)
  {
    PROCESS_VAL(Gauss1Line);
    PROCESS_VAL(Gauss2Line);
    PROCESS_VAL(Gauss3Line);
    PROCESS_VAL(Gauss4Line);
    PROCESS_VAL(Gauss5Line);
    PROCESS_VAL(Gauss6Line);
    PROCESS_VAL(Gauss7Line);
    PROCESS_VAL(Gauss8Line);
    PROCESS_VAL(Gauss9Line);
    PROCESS_VAL(Gauss10Line);
    PROCESS_VAL(Gauss11Line);
    PROCESS_VAL(Gauss12Line);
    PROCESS_VAL(Gauss2W1Line);
    PROCESS_VAL(Gauss4W1Line);
    PROCESS_VAL(Gauss6W1Line);
    PROCESS_VAL(Gauss8W1Line);
    PROCESS_VAL(Gauss16W2Line);
    PROCESS_VAL(BaryCenterTria);
    PROCESS_VAL(MidPointTria);
    PROCESS_VAL(SevenPointTria);
    PROCESS_VAL(Gauss3Tria);
    PROCESS_VAL(VertexTria);
    PROCESS_VAL(Degree8Tria);
    PROCESS_VAL(Degree9Tria);
    PROCESS_VAL(Degree11Tria);
    PROCESS_VAL(Degree19Tria);
    PROCESS_VAL(Gauss2Quad);
    PROCESS_VAL(Gauss3Quad);
    PROCESS_VAL(Gauss4Quad);
    PROCESS_VAL(Gauss5Quad);
    PROCESS_VAL(Gauss6Quad);
    PROCESS_VAL(Gauss7Quad);
    PROCESS_VAL(Gauss8Quad);
    PROCESS_VAL(Gauss9Quad);
    PROCESS_VAL(VertexQuad);
    PROCESS_VAL(SimpsonQuad);
    PROCESS_VAL(CompGauss3Tria);
    PROCESS_VAL(CompGauss4Tria);
    PROCESS_VAL(Gauss_Degree8Tria);
    PROCESS_VAL(BaryCenterTetra);
    PROCESS_VAL(VertexTetra);
    PROCESS_VAL(P2Tetra);
    PROCESS_VAL(P4Tetra);
    PROCESS_VAL(P5Tetra);
    PROCESS_VAL(P8Tetra);
    PROCESS_VAL(VertexHexa);
    PROCESS_VAL(Gauss2Hexa);
    PROCESS_VAL(Gauss3Hexa);
    PROCESS_VAL(Gauss4Hexa);
    PROCESS_VAL(Gauss5Hexa);
    PROCESS_VAL(Gauss6Hexa);
    PROCESS_VAL(Gauss7Hexa);
    PROCESS_VAL(Gauss8Hexa);
    PROCESS_VAL(Gauss9Hexa);
    PROCESS_VAL(VerticesAndOrigin);
    PROCESS_VAL(VerticesAndOrigin15);
    PROCESS_VAL(VerticesAndOrigin57);
    PROCESS_VAL(Degree7_Points38);
    default:
      s = "unknown QuadratureFormula_type type";
      break;
  }
#undef PROCESS_VAL
  return out << s;
}
