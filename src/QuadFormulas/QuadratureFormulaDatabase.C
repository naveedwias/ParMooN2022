#include "QuadratureFormulaDatabase.h"
#include "QuadFormula.h"
#include <array>


/** @brief quadrature formulas */
std::array<TQuadFormula *, N_QuadratureFormula_type> QuadFormulas;

/** @brief mapping the accuracy to a quadrature formula type for each 
 * reference element
 */
std::array<std::array<QuadratureFormula_type, MAXDEGREE>, N_BFRefElements>
    quad_formula_from_degree;

/** @brief mapping the reference element to the highest accuracy of a 
 * quadrature formula we have implemented
 */
std::array<int, N_BFRefElements> highest_accuracy;

///////////////////////////////////////////////////////////////////////////////

void QuadratureFormulaDatabase::create()
{
  QuadFormulas.fill(nullptr);
  using QFt = QuadratureFormula_type; // only for less typing
  // =====================================================================
  // generate accuracy->qf_type maps
  // =====================================================================
  std::array<QFt, MAXDEGREE> quad_formulas;
  quad_formulas[0] = QFt::P2Tetra;
  quad_formulas[1] = QFt::P2Tetra;
  quad_formulas[2] = QFt::P2Tetra;
  quad_formulas[3] = QFt::P4Tetra;
  quad_formulas[4] = QFt::P4Tetra;
  quad_formulas[5] = QFt::P5Tetra;
  quad_formulas[6] = QFt::P5Tetra; /// @TODO This should be P6Tetra
  quad_formulas[7] = QFt::P8Tetra;
  quad_formulas[8] = QFt::P8Tetra;
  quad_formula_from_degree[static_cast<int>(BFRefElements::BFUnitTetrahedron)]
    = quad_formulas;
  highest_accuracy[static_cast<int>(BFRefElements::BFUnitTetrahedron)] = 8;
  
  quad_formulas[0] = QFt::Gauss2Hexa;
  quad_formulas[1] = QFt::Gauss2Hexa;
  quad_formulas[2] = QFt::Gauss2Hexa;
  quad_formulas[3] = QFt::Gauss3Hexa;
  quad_formulas[4] = QFt::Gauss3Hexa;
  quad_formulas[5] = QFt::Gauss3Hexa;
  quad_formulas[6] = QFt::Gauss4Hexa;
  quad_formulas[7] = QFt::Gauss4Hexa;
  quad_formulas[8] = QFt::Gauss5Hexa;
  quad_formulas[9] = QFt::Gauss5Hexa;
  quad_formulas[10] = QFt::Gauss6Hexa;
  quad_formulas[11] = QFt::Gauss6Hexa;
  quad_formulas[12] = QFt::Gauss7Hexa;
  quad_formulas[13] = QFt::Gauss7Hexa;
  quad_formulas[14] = QFt::Gauss8Hexa;
  quad_formulas[15] = QFt::Gauss8Hexa;
  quad_formulas[16] = QFt::Gauss9Hexa;
  quad_formulas[17] = QFt::Gauss9Hexa;
  quad_formula_from_degree[static_cast<int>(BFRefElements::BFUnitHexahedron)]
    = quad_formulas;
  highest_accuracy[static_cast<int>(BFRefElements::BFUnitHexahedron)] = 17;
  
  quad_formulas[0] = QFt::Gauss2Quad;
  quad_formulas[1] = QFt::Gauss2Quad;
  quad_formulas[2] = QFt::Gauss2Quad;
  quad_formulas[3] = QFt::Gauss3Quad;
  quad_formulas[4] = QFt::Gauss3Quad;
  quad_formulas[5] = QFt::Gauss3Quad;
  quad_formulas[6] = QFt::Gauss4Quad;
  quad_formulas[7] = QFt::Gauss4Quad;
  quad_formulas[8] = QFt::Gauss5Quad;
  quad_formulas[9] = QFt::Gauss5Quad;
  quad_formulas[10] = QFt::Gauss6Quad;
  quad_formulas[11] = QFt::Gauss6Quad;
  quad_formulas[12] = QFt::Gauss7Quad;
  quad_formulas[13] = QFt::Gauss7Quad;
  quad_formulas[14] = QFt::Gauss8Quad;
  quad_formulas[15] = QFt::Gauss8Quad;
  quad_formulas[16] = QFt::Gauss9Quad;
  quad_formulas[17] = QFt::Gauss9Quad;
  quad_formula_from_degree[static_cast<int>(BFRefElements::BFUnitSquare)]
    = quad_formulas;
  highest_accuracy[static_cast<int>(BFRefElements::BFUnitSquare)] = 17;

  quad_formulas[0] = QFt::Gauss3Tria;
  quad_formulas[1] = QFt::Gauss3Tria;
  quad_formulas[2] = QFt::Gauss3Tria;
  quad_formulas[3] = QFt::Gauss3Tria;
  quad_formulas[4] = QFt::Gauss3Tria;
  quad_formulas[5] = QFt::Gauss3Tria;
  quad_formulas[6] = QFt::Degree9Tria;
  quad_formulas[7] = QFt::Degree9Tria;
  quad_formulas[8] = QFt::Degree9Tria;
  quad_formulas[9] = QFt::Degree9Tria;
  quad_formulas[10] = QFt::Degree19Tria;
  quad_formulas[11] = QFt::Degree19Tria;
  quad_formulas[12] = QFt::Degree19Tria;
  quad_formulas[13] = QFt::Degree19Tria;
  quad_formulas[14] = QFt::Degree19Tria;
  quad_formulas[15] = QFt::Degree19Tria;
  quad_formulas[16] = QFt::Degree19Tria;
  quad_formulas[17] = QFt::Degree19Tria;
  quad_formulas[18] = QFt::Degree19Tria;
  quad_formulas[19] = QFt::Degree19Tria;
  quad_formula_from_degree[static_cast<int>(BFRefElements::BFUnitTriangle)]
    = quad_formulas;
  highest_accuracy[static_cast<int>(BFRefElements::BFUnitTriangle)] = 19;

  quad_formulas[0] = QFt::Gauss1Line;
  quad_formulas[1] = QFt::Gauss1Line;
  quad_formulas[2] = QFt::Gauss2Line;
  quad_formulas[3] = QFt::Gauss2Line;
  quad_formulas[4] = QFt::Gauss3Line;
  quad_formulas[5] = QFt::Gauss3Line;
  quad_formulas[6] = QFt::Gauss4Line;
  quad_formulas[7] = QFt::Gauss4Line;
  quad_formulas[8] = QFt::Gauss5Line;
  quad_formulas[9] = QFt::Gauss5Line;
  quad_formulas[10] = QFt::Gauss6Line;
  quad_formulas[11] = QFt::Gauss6Line;
  quad_formulas[12] = QFt::Gauss7Line;
  quad_formulas[13] = QFt::Gauss7Line;
  quad_formulas[14] = QFt::Gauss8Line;
  quad_formulas[15] = QFt::Gauss8Line;
  quad_formulas[16] = QFt::Gauss9Line;
  quad_formulas[17] = QFt::Gauss9Line;
  quad_formulas[18] = QFt::Gauss10Line;
  quad_formulas[19] = QFt::Gauss10Line;
  quad_formulas[20] = QFt::Gauss11Line;
  quad_formulas[21] = QFt::Gauss11Line;
  quad_formulas[22] = QFt::Gauss12Line;
  quad_formulas[23] = QFt::Gauss12Line;
  quad_formula_from_degree[static_cast<int>(BFRefElements::BFUnitLine)]
    = quad_formulas;
  highest_accuracy[static_cast<int>(BFRefElements::BFUnitLine)] = 23;
}

const TQuadFormula *QuadratureFormulaDatabase::qf_from_degree(
    int accuracy, BFRefElements RefElem)
{
  // make sure accuracy is at most the respective highest possible accuracy
  accuracy = std::min(accuracy, highest_accuracy[static_cast<int>(RefElem)]);
  auto qf_type = quad_formula_from_degree[static_cast<int>(RefElem)][accuracy];
  if(QuadFormulas[static_cast<int>(qf_type)] == nullptr)
  {
    QuadFormulas[static_cast<int>(qf_type)] = new TQuadFormula(qf_type);
  }
  return QuadFormulas[static_cast<int>(qf_type)];
}

const TQuadFormula *QuadratureFormulaDatabase::qf_from_degree(int accuracy,
                                                              Shapes shape)
{
  // avoid compiler warning 'maybe-uninitialized'
  BFRefElements ref_element = BFRefElements::BFUnitLine;
  switch(shape)
  {
    case S_Line:
      ref_element = BFRefElements::BFUnitLine;
      break;
    case Triangle:
      ref_element = BFRefElements::BFUnitTriangle;
      break;
    case Quadrangle:
    case Parallelogram:
    case Rectangle:
      ref_element = BFRefElements::BFUnitSquare;
      break;
    case Tetrahedron:
      ref_element = BFRefElements::BFUnitTetrahedron;
      break;
    case Hexahedron:
    case Brick:
      ref_element = BFRefElements::BFUnitHexahedron;
      break;
    // no default, so that a compiler warning is generated for missing cases
  }
  return qf_from_degree(accuracy, ref_element);
}

void QuadratureFormulaDatabase::destroy()
{
  for(auto& qf : QuadFormulas)
  {
    delete qf;
    qf = nullptr;
  }
}
