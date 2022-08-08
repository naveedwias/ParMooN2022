#include "Enumerations_fe.h"
#include "MooNMD_Io.h"

std::ostream& operator<<(std::ostream& out, const BFRefElements t)
{
  switch(t)
  {
    case BFRefElements::BFUnitLine: out << "BFUnitLine"; break;
    case BFRefElements::BFUnitTriangle: out << "BFUnitTriangle"; break;
    case BFRefElements::BFUnitSquare: out << "BFUnitSquare"; break;
    case BFRefElements::BFUnitTetrahedron: out << "BFUnitTetrahedron"; break;
    case BFRefElements::BFUnitHexahedron: out << "BFUnitHexahedron"; break;
    default: out << "unknown reference element"; break;
  }
  return out;
}

std::ostream& operator<<(std::ostream& out, const HNDesc t)
{
  const char* s = 0;
#define PROCESS_VAL(p)                                                         \
  case(HNDesc::p):                                                             \
    s = #p;                                                                    \
    break;
  switch(t)
  {
    PROCESS_VAL(HN_C_P1_2D_0);
    PROCESS_VAL(HN_C_P2_2D_0);
    PROCESS_VAL(HN_C_P2_2D_1);
    PROCESS_VAL(HN_C_P3_2D_0);
    PROCESS_VAL(HN_C_P3_2D_1);
    PROCESS_VAL(HN_C_P3_2D_2);
    PROCESS_VAL(HN_C_P4_2D_0);
    PROCESS_VAL(HN_C_P4_2D_1);
    PROCESS_VAL(HN_C_P4_2D_2);
    PROCESS_VAL(HN_C_P4_2D_3);
    PROCESS_VAL(HN_C_P5_2D_0);
    PROCESS_VAL(HN_C_P5_2D_1);
    PROCESS_VAL(HN_C_P5_2D_2);
    PROCESS_VAL(HN_C_P5_2D_3);
    PROCESS_VAL(HN_C_P5_2D_4);
    PROCESS_VAL(HN_C_P6_2D_0);
    PROCESS_VAL(HN_C_P6_2D_1);
    PROCESS_VAL(HN_C_P6_2D_2);
    PROCESS_VAL(HN_C_P6_2D_3);
    PROCESS_VAL(HN_C_P6_2D_4);
    PROCESS_VAL(HN_C_P6_2D_5);
    PROCESS_VAL(HN_C_P7_2D_0);
    PROCESS_VAL(HN_C_P7_2D_1);
    PROCESS_VAL(HN_C_P7_2D_2);
    PROCESS_VAL(HN_C_P7_2D_3);
    PROCESS_VAL(HN_C_P7_2D_4);
    PROCESS_VAL(HN_C_P7_2D_5);
    PROCESS_VAL(HN_C_P7_2D_6);
    PROCESS_VAL(HN_C_P8_2D_0);
    PROCESS_VAL(HN_C_P8_2D_1);
    PROCESS_VAL(HN_C_P8_2D_2);
    PROCESS_VAL(HN_C_P8_2D_3);
    PROCESS_VAL(HN_C_P8_2D_4);
    PROCESS_VAL(HN_C_P8_2D_5);
    PROCESS_VAL(HN_C_P8_2D_6);
    PROCESS_VAL(HN_C_P8_2D_7);
    PROCESS_VAL(HN_C_P9_2D_0);
    PROCESS_VAL(HN_C_P9_2D_1);
    PROCESS_VAL(HN_C_P9_2D_2);
    PROCESS_VAL(HN_C_P9_2D_3);
    PROCESS_VAL(HN_C_P9_2D_4);
    PROCESS_VAL(HN_C_P9_2D_5);
    PROCESS_VAL(HN_C_P9_2D_6);
    PROCESS_VAL(HN_C_P9_2D_7);
    PROCESS_VAL(HN_C_P9_2D_8);
    PROCESS_VAL(HN_N_P1_2D_0);
    PROCESS_VAL(HN_N_P2_2D_0);
    PROCESS_VAL(HN_N_P3_2D_0);
    PROCESS_VAL(HN_N_P4_2D_0);
    PROCESS_VAL(HN_N_P5_2D_0);
    PROCESS_VAL(HN_C_Q1_3D_E);
    PROCESS_VAL(HN_C_Q1_3D_F);
    PROCESS_VAL(HN_C_Q2_3D_E);
    PROCESS_VAL(HN_C_Q2_3D_F);
    PROCESS_VAL(HN_C_Q3_3D_1);
    PROCESS_VAL(HN_C_Q3_3D_2);
    PROCESS_VAL(HN_C_Q3_3D_3);
    PROCESS_VAL(HN_C_Q3_3D_4);
    PROCESS_VAL(HN_C_Q3_3D_5);
    PROCESS_VAL(HN_C_P1_3D_E);
    PROCESS_VAL(HN_C_P2_3D_E);
    PROCESS_VAL(HN_C_P2_3D_F);
    PROCESS_VAL(HN_C_P3_3D_E);
    PROCESS_VAL(HN_C_P3_3D_M);
    PROCESS_VAL(HN_C_P3_3D_F);
    PROCESS_VAL(HN_C_P3_3D_G);
    PROCESS_VAL(HN_N_P1_3D_E);
    PROCESS_VAL(HN_N_P2_3D_0);
    PROCESS_VAL(HN_N_P2_3D_1);
    PROCESS_VAL(HN_N_P2_3D_2);
    default:
      s = "unknown HNDesc type";
      break;
  }
#undef PROCESS_VAL
  return out << s;
}
