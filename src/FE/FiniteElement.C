#ifdef _MPI
#  include "mpi.h"
#endif

#include <vector>
#include <cmath>
#include "FiniteElement.h"
#include "MooNMD_Io.h"


std::tuple<BaseFunction_type, NodalFunctional_type, FEDescriptor_type, ReferenceTransformation_type> get_ids(FE_type id)
{
  // for shorter code
  auto mt = std::make_tuple<BaseFunction_type, NodalFunctional_type, FEDescriptor_type,
                            ReferenceTransformation_type>;
  switch(id)
  {
    case C_P0_1D_L_A:
      return mt(BF_C_L_P0_1D, NF_C_L_P0_1D, FE_C_L_P0_1D,
                ReferenceTransformation_type::LineAffin);
      break;
    case C_P1_1D_L_A:
      return mt(BF_C_L_P1_1D, NF_C_L_P1_1D, FE_C_L_P1_1D,
                ReferenceTransformation_type::LineAffin);
      break;
    case C_P2_1D_L_A:
      return mt(BF_C_L_P2_1D, NF_C_L_P2_1D, FE_C_L_P2_1D,
                ReferenceTransformation_type::LineAffin);
      break;
    case C_P3_1D_L_A:
      return mt(BF_C_L_P3_1D, NF_C_L_P3_1D, FE_C_L_P3_1D,
                ReferenceTransformation_type::LineAffin);
      break;
    case N_P0_1D_L_A:
      return mt(BF_C_L_P0_1D, NF_C_L_P0_1D, FE_N_L_P0_1D,
                ReferenceTransformation_type::LineAffin);
      break;
    case D_P1_1D_L_A:
      return mt(BF_D_L_P1_1D, NF_D_L_P1_1D, FE_D_L_P1_1D,
                ReferenceTransformation_type::LineAffin);
      break;
    case D_P2_1D_L_A:
      return mt(BF_D_L_P2_1D, NF_D_L_P2_1D, FE_D_L_P2_1D,
                ReferenceTransformation_type::LineAffin);
      break;
// 2D
    case C_P00_2D_T_A:
      return mt(BF_C_T_P00_2D, NF_C_T_P00_2D, FE_C_T_P00_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case C_P0_2D_T_A:
      return mt(BF_C_T_P0_2D, NF_C_T_P0_2D, FE_C_T_P0_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case C_P1_2D_T_A:
      return mt(BF_C_T_P1_2D, NF_C_T_P1_2D, FE_C_T_P1_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case C_P2_2D_T_A:
      return mt(BF_C_T_P2_2D, NF_C_T_P2_2D, FE_C_T_P2_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case C_P3_2D_T_A:
      return mt(BF_C_T_P3_2D, NF_C_T_P3_2D, FE_C_T_P3_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case C_P4_2D_T_A:
      return mt(BF_C_T_P4_2D, NF_C_T_P4_2D, FE_C_T_P4_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case C_P5_2D_T_A:
      return mt(BF_C_T_P5_2D, NF_C_T_P5_2D, FE_C_T_P5_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case C_P6_2D_T_A:
      return mt(BF_C_T_P6_2D, NF_C_T_P6_2D, FE_C_T_P6_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case C_P7_2D_T_A:
      ErrThrow("finite element of type C_P7_2D_T_A not implemented")
      break;
    case C_P8_2D_T_A:
      ErrThrow("finite element of type C_P8_2D_T_A not implemented")
      break;
    case C_P9_2D_T_A:
      ErrThrow("finite element of type C_P9_2D_T_A not implemented")
      break;
    case N_P1_2D_T_A:
      return mt(BF_N_T_P1_2D, NF_N_T_P1_2D, FE_N_T_P1_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case C_Q00_2D_Q_A:
      return mt(BF_C_Q_Q00_2D, NF_C_Q_Q00_2D, FE_C_Q_Q00_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_Q0_2D_Q_A:
      return mt(BF_C_Q_Q0_2D, NF_C_Q_Q0_2D, FE_C_Q_Q0_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_Q1_2D_Q_A:
      return mt(BF_C_Q_Q1_2D, NF_C_Q_Q1_2D, FE_C_Q_Q1_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_Q2_2D_Q_A:
      return mt(BF_C_Q_Q2_2D, NF_C_Q_Q2_2D, FE_C_Q_Q2_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_Q3_2D_Q_A:
      return mt(BF_C_Q_Q3_2D, NF_C_Q_Q3_2D, FE_C_Q_Q3_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_Q4_2D_Q_A:
      return mt(BF_C_Q_Q4_2D, NF_C_Q_Q4_2D, FE_C_Q_Q4_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_Q5_2D_Q_A:
      return mt(BF_C_Q_Q5_2D, NF_C_Q_Q5_2D, FE_C_Q_Q5_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_Q6_2D_Q_A:
      return mt(BF_C_Q_Q6_2D, NF_C_Q_Q6_2D, FE_C_Q_Q6_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_Q7_2D_Q_A:
      return mt(BF_C_Q_Q7_2D, NF_C_Q_Q7_2D, FE_C_Q_Q7_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_Q8_2D_Q_A:
      return mt(BF_C_Q_Q8_2D, NF_C_Q_Q8_2D, FE_C_Q_Q8_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_Q9_2D_Q_A:
      return mt(BF_C_Q_Q9_2D, NF_C_Q_Q9_2D, FE_C_Q_Q9_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case N_Q1_2D_Q_A:
      return mt(BF_N_Q_Q1_2D, NF_N_Q_Q1_2D, FE_N_Q_Q1_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_Q00_2D_Q_M:
      return mt(BF_C_Q_Q00_2D, NF_C_Q_Q00_2D, FE_C_Q_Q00_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_Q0_2D_Q_M:
      return mt(BF_C_Q_Q0_2D, NF_C_Q_Q0_2D, FE_C_Q_Q0_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_Q1_2D_Q_M:
      return mt(BF_C_Q_Q1_2D, NF_C_Q_Q1_2D, FE_C_Q_Q1_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_Q2_2D_Q_M:
      return mt(BF_C_Q_Q2_2D, NF_C_Q_Q2_2D, FE_C_Q_Q2_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_Q3_2D_Q_M:
      return mt(BF_C_Q_Q3_2D, NF_C_Q_Q3_2D, FE_C_Q_Q3_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_Q4_2D_Q_M:
      return mt(BF_C_Q_Q4_2D, NF_C_Q_Q4_2D, FE_C_Q_Q4_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_Q5_2D_Q_M:
      return mt(BF_C_Q_Q5_2D, NF_C_Q_Q5_2D, FE_C_Q_Q5_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_Q6_2D_Q_M:
      return mt(BF_C_Q_Q6_2D, NF_C_Q_Q6_2D, FE_C_Q_Q6_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_Q7_2D_Q_M:
      return mt(BF_C_Q_Q7_2D, NF_C_Q_Q7_2D, FE_C_Q_Q7_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_Q8_2D_Q_M:
      return mt(BF_C_Q_Q8_2D, NF_C_Q_Q8_2D, FE_C_Q_Q8_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_Q9_2D_Q_M:
      return mt(BF_C_Q_Q9_2D, NF_C_Q_Q9_2D, FE_C_Q_Q9_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case N_Q1_2D_Q_M:
      return mt(BF_N_Q_Q1_2D, NF_N_Q_Q1_2D, FE_N_Q_Q1_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case D_P1_2D_Q_A:
      return mt(BF_D_Q_P1_2D, NF_D_Q_P1_2D, FE_D_Q_P1_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case D_P2_2D_Q_A:
      return mt(BF_D_Q_P2_2D, NF_D_Q_P2_2D, FE_D_Q_P2_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case D_P3_2D_Q_A:
      return mt(BF_D_Q_P3_2D, NF_D_Q_P3_2D, FE_D_Q_P3_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case D_P1_2D_Q_M:
      return mt(BF_D_Q_P1_2D, NF_D_Q_P1_2D, FE_D_Q_P1_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case D_P2_2D_Q_M:
      return mt(BF_D_Q_P2_2D, NF_D_Q_P2_2D, FE_D_Q_P2_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case D_P3_2D_Q_M:
      return mt(BF_D_Q_P3_2D, NF_D_Q_P3_2D, FE_D_Q_P3_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_B2_2D_T_A:
      return mt(BF_C_T_B2_2D, NF_C_T_B2_2D, FE_C_T_B2_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case C_B3_2D_T_A:
      return mt(BF_C_T_B3_2D, NF_C_T_B3_2D, FE_C_T_B3_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case C_SV2_2D_T_A:
      return mt(BF_C_T_SV2_2D, NF_C_T_SV2_2D, FE_C_T_SV2_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case D_P1_2D_T_A:
      return mt(BF_D_T_P1_2D, NF_D_T_P1_2D, FE_D_T_P1_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case D_P2_2D_T_A:
      return mt(BF_D_T_P2_2D, NF_D_T_P2_2D, FE_D_T_P2_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case D_SV1_2D_T_A:
      return mt(BF_D_T_SV1_2D, NF_D_T_SV1_2D, FE_D_T_SV1_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case N_Q2_2D_Q_A:
      return mt(BF_N_Q_Q2_2D, NF_N_Q_Q2_2D, FE_N_Q_Q2_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case N_Q3_2D_Q_A:
      return mt(BF_N_Q_Q3_2D, NF_N_Q_Q3_2D, FE_N_Q_Q3_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case N_Q4_2D_Q_A:
      return mt(BF_N_Q_Q4_2D, NF_N_Q_Q4_2D, FE_N_Q_Q4_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case N_Q5_2D_Q_A:
      return mt(BF_N_Q_Q5_2D, NF_N_Q_Q5_2D, FE_N_Q_Q5_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case N_Q2_2D_Q_M:
      return mt(BF_N_Q_Q2_2D, NF_N_Q_Q2_2D, FE_N_Q_Q2_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case N_Q3_2D_Q_M:
      return mt(BF_N_Q_Q3_2D, NF_N_Q_Q3_2D, FE_N_Q_Q3_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case N_Q4_2D_Q_M:
      return mt(BF_N_Q_Q4_2D, NF_N_Q_Q4_2D, FE_N_Q_Q4_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case N_Q5_2D_Q_M:
      return mt(BF_N_Q_Q5_2D, NF_N_Q_Q5_2D, FE_N_Q_Q5_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case D_P4_2D_Q_A:
      return mt(BF_D_Q_P4_2D, NF_D_Q_P4_2D, FE_D_Q_P4_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case D_P5_2D_Q_A:
      return mt(BF_D_Q_P5_2D, NF_D_Q_P5_2D, FE_D_Q_P5_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case D_P6_2D_Q_A:
      return mt(BF_D_Q_P6_2D, NF_D_Q_P6_2D, FE_D_Q_P6_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case D_P7_2D_Q_A:
      return mt(BF_D_Q_P7_2D, NF_D_Q_P7_2D, FE_D_Q_P7_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case D_P4_2D_Q_M:
      return mt(BF_D_Q_P4_2D, NF_D_Q_P4_2D, FE_D_Q_P4_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case D_P5_2D_Q_M:
      return mt(BF_D_Q_P5_2D, NF_D_Q_P5_2D, FE_D_Q_P5_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case D_P6_2D_Q_M:
      return mt(BF_D_Q_P6_2D, NF_D_Q_P6_2D, FE_D_Q_P6_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case D_P7_2D_Q_M:
      return mt(BF_D_Q_P7_2D, NF_D_Q_P7_2D, FE_D_Q_P7_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case N_P1MOD_2D_T_A:
      return mt(BF_N_T_P1MOD_2D, NF_N_T_P1MOD_2D, FE_N_T_P1MOD_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case N_P2_2D_T_A:
      return mt(BF_N_T_P2_2D, NF_N_T_P2_2D, FE_N_T_P2_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case N_P3_2D_T_A:
      return mt(BF_N_T_P3_2D, NF_N_T_P3_2D, FE_N_T_P3_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case N_P4_2D_T_A:
      return mt(BF_N_T_P4_2D, NF_N_T_P4_2D, FE_N_T_P4_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case N_P5_2D_T_A:
      return mt(BF_N_T_P5_2D, NF_N_T_P5_2D, FE_N_T_P5_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case D_P3_2D_T_A:
      return mt(BF_D_T_P3_2D, NF_D_T_P3_2D, FE_D_T_P3_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case D_P4_2D_T_A:
      return mt(BF_D_T_P4_2D, NF_D_T_P4_2D, FE_D_T_P4_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case C_P1MINI_2D_T_A:
      return mt(BF_C_T_P1MINI_2D, NF_C_T_P1MINI_2D, FE_C_T_P1MINI_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case B_IB2_2D_Q_A:
      return mt(BF_B_Q_IB2_2D, NF_B_Q_IB2_2D, FE_B_Q_IB2_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case B_IB2_2D_Q_M:
      return mt(BF_B_Q_IB2_2D, NF_B_Q_IB2_2D, FE_B_Q_IB2_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case D_Q1_2D_Q_A:
      return mt(BF_D_Q_Q1_2D, NF_D_Q_Q1_2D, FE_D_Q_Q1_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case D_Q2_2D_Q_A:
      return mt(BF_D_Q_Q2_2D, NF_D_Q_Q2_2D, FE_D_Q_Q2_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case D_Q3_2D_Q_A:
      return mt(BF_D_Q_Q3_2D, NF_D_Q_Q3_2D, FE_D_Q_Q3_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case D_Q4_2D_Q_A:
      return mt(BF_D_Q_Q4_2D, NF_D_Q_Q4_2D, FE_D_Q_Q4_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case D_Q1_2D_Q_M:
      return mt(BF_D_Q_Q1_2D, NF_D_Q_Q1_2D, FE_D_Q_Q1_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case D_Q2_2D_Q_M:
      return mt(BF_D_Q_Q2_2D, NF_D_Q_Q2_2D, FE_D_Q_Q2_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case D_Q3_2D_Q_M:
      return mt(BF_D_Q_Q3_2D, NF_D_Q_Q3_2D, FE_D_Q_Q3_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case D_Q4_2D_Q_M:
      return mt(BF_D_Q_Q4_2D, NF_D_Q_Q4_2D, FE_D_Q_Q4_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_B4_2D_T_A:
      return mt(BF_C_T_B4_2D, NF_C_T_B4_2D, FE_C_T_B4_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case D_D2_2D_Q_A:
      return mt(BF_D_Q_D2_2D, NF_D_Q_D2_2D, FE_D_Q_D2_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case D_D2_2D_Q_M:
      return mt(BF_D_Q_D2_2D, NF_D_Q_D2_2D, FE_D_Q_D2_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_UL1_2D_Q_A:
      return mt(BF_C_Q_UL1_2D, NF_C_Q_UL1_2D, FE_C_Q_UL1_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_UL2_2D_Q_A:
      return mt(BF_C_Q_UL2_2D, NF_C_Q_UL2_2D, FE_C_Q_UL2_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_UL3_2D_Q_A:
      return mt(BF_C_Q_UL3_2D, NF_C_Q_UL3_2D, FE_C_Q_UL3_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_UL4_2D_Q_A:
      return mt(BF_C_Q_UL4_2D, NF_C_Q_UL4_2D, FE_C_Q_UL4_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_UL5_2D_Q_A:
      return mt(BF_C_Q_UL5_2D, NF_C_Q_UL5_2D, FE_C_Q_UL5_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_UL1_2D_Q_M:
      return mt(BF_C_Q_UL1_2D, NF_C_Q_UL1_2D, FE_C_Q_UL1_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_UL2_2D_Q_M:
      return mt(BF_C_Q_UL2_2D, NF_C_Q_UL2_2D, FE_C_Q_UL2_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_UL3_2D_Q_M:
      return mt(BF_C_Q_UL3_2D, NF_C_Q_UL3_2D, FE_C_Q_UL3_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_UL4_2D_Q_M:
      return mt(BF_C_Q_UL4_2D, NF_C_Q_UL4_2D, FE_C_Q_UL4_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_UL5_2D_Q_M:
      return mt(BF_C_Q_UL5_2D, NF_C_Q_UL5_2D, FE_C_Q_UL5_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_UL1_2D_T_A:
      return mt(BF_C_T_UL1_2D, NF_C_T_UL1_2D, FE_C_T_UL1_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case C_UL2_2D_T_A:
      return mt(BF_C_T_UL2_2D, NF_C_T_UL2_2D, FE_C_T_UL2_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case C_UL3_2D_T_A:
      return mt(BF_C_T_UL3_2D, NF_C_T_UL3_2D, FE_C_T_UL3_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case C_UL4_2D_T_A:
      return mt(BF_C_T_UL4_2D, NF_C_T_UL4_2D, FE_C_T_UL4_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case C_UL5_2D_T_A:
      return mt(BF_C_T_UL5_2D, NF_C_T_UL5_2D, FE_C_T_UL5_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case C_UL2S_2D_Q_A:
      return mt(BF_C_Q_UL2S_2D, NF_C_Q_UL2S_2D, FE_C_Q_UL2S_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_UL3S_2D_Q_A:
      return mt(BF_C_Q_UL3S_2D, NF_C_Q_UL3S_2D, FE_C_Q_UL3S_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_UL4S_2D_Q_A:
      return mt(BF_C_Q_UL4S_2D, NF_C_Q_UL4S_2D, FE_C_Q_UL4S_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_UL5S_2D_Q_A:
      return mt(BF_C_Q_UL5S_2D, NF_C_Q_UL5S_2D, FE_C_Q_UL5S_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_UL6S_2D_Q_A:
      return mt(BF_C_Q_UL6S_2D, NF_C_Q_UL6S_2D, FE_C_Q_UL6S_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_UL7S_2D_Q_A:
      return mt(BF_C_Q_UL7S_2D, NF_C_Q_UL7S_2D, FE_C_Q_UL7S_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_UL8S_2D_Q_A:
      return mt(BF_C_Q_UL8S_2D, NF_C_Q_UL8S_2D, FE_C_Q_UL8S_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_UL9S_2D_Q_A:
      return mt(BF_C_Q_UL9S_2D, NF_C_Q_UL9S_2D, FE_C_Q_UL9S_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_UL2S_2D_Q_M:
      return mt(BF_C_Q_UL2S_2D, NF_C_Q_UL2S_2D, FE_C_Q_UL2S_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_UL3S_2D_Q_M:
      return mt(BF_C_Q_UL3S_2D, NF_C_Q_UL3S_2D, FE_C_Q_UL3S_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_UL4S_2D_Q_M:
      return mt(BF_C_Q_UL4S_2D, NF_C_Q_UL4S_2D, FE_C_Q_UL4S_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_UL5S_2D_Q_M:
      return mt(BF_C_Q_UL5S_2D, NF_C_Q_UL5S_2D, FE_C_Q_UL5S_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_UL6S_2D_Q_M:
      return mt(BF_C_Q_UL6S_2D, NF_C_Q_UL6S_2D, FE_C_Q_UL6S_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_UL7S_2D_Q_M:
      return mt(BF_C_Q_UL7S_2D, NF_C_Q_UL7S_2D, FE_C_Q_UL7S_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_UL8S_2D_Q_M:
      return mt(BF_C_Q_UL8S_2D, NF_C_Q_UL8S_2D, FE_C_Q_UL8S_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_UL9S_2D_Q_M:
      return mt(BF_C_Q_UL9S_2D, NF_C_Q_UL9S_2D, FE_C_Q_UL9S_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_UL2SE_2D_Q_A:
      return mt(BF_C_Q_UL2SE_2D, NF_C_Q_UL2SE_2D, FE_C_Q_UL2SE_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_UL3SE_2D_Q_A:
      return mt(BF_C_Q_UL3SE_2D, NF_C_Q_UL3SE_2D, FE_C_Q_UL3SE_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_UL4SE_2D_Q_A:
      return mt(BF_C_Q_UL4SE_2D, NF_C_Q_UL4SE_2D, FE_C_Q_UL4SE_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_UL5SE_2D_Q_A:
      return mt(BF_C_Q_UL5SE_2D, NF_C_Q_UL5SE_2D, FE_C_Q_UL5SE_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_UL6SE_2D_Q_A:
      return mt(BF_C_Q_UL6SE_2D, NF_C_Q_UL6SE_2D, FE_C_Q_UL6SE_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_UL7SE_2D_Q_A:
      return mt(BF_C_Q_UL7SE_2D, NF_C_Q_UL7SE_2D, FE_C_Q_UL7SE_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_UL8SE_2D_Q_A:
      return mt(BF_C_Q_UL8SE_2D, NF_C_Q_UL8SE_2D, FE_C_Q_UL8SE_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_UL9SE_2D_Q_A:
      return mt(BF_C_Q_UL9SE_2D, NF_C_Q_UL9SE_2D, FE_C_Q_UL9SE_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_UL2SE_2D_Q_M:
      return mt(BF_C_Q_UL2SE_2D, NF_C_Q_UL2SE_2D, FE_C_Q_UL2SE_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_UL3SE_2D_Q_M:
      return mt(BF_C_Q_UL3SE_2D, NF_C_Q_UL3SE_2D, FE_C_Q_UL3SE_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_UL4SE_2D_Q_M:
      return mt(BF_C_Q_UL4SE_2D, NF_C_Q_UL4SE_2D, FE_C_Q_UL4SE_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_UL5SE_2D_Q_M:
      return mt(BF_C_Q_UL5SE_2D, NF_C_Q_UL5SE_2D, FE_C_Q_UL5SE_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_UL6SE_2D_Q_M:
      return mt(BF_C_Q_UL6SE_2D, NF_C_Q_UL6SE_2D, FE_C_Q_UL6SE_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_UL7SE_2D_Q_M:
      return mt(BF_C_Q_UL7SE_2D, NF_C_Q_UL7SE_2D, FE_C_Q_UL7SE_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_UL8SE_2D_Q_M:
      return mt(BF_C_Q_UL8SE_2D, NF_C_Q_UL8SE_2D, FE_C_Q_UL8SE_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_UL9SE_2D_Q_M:
      return mt(BF_C_Q_UL9SE_2D, NF_C_Q_UL9SE_2D, FE_C_Q_UL9SE_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_M2_2D_Q_A:
      return mt(BF_C_Q_M2_2D, NF_C_Q_M2_2D, FE_C_Q_M2_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_M3_2D_Q_A:
      return mt(BF_C_Q_M3_2D, NF_C_Q_M3_2D, FE_C_Q_M3_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_M4_2D_Q_A:
      return mt(BF_C_Q_M4_2D, NF_C_Q_M4_2D, FE_C_Q_M4_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_M5_2D_Q_A:
      return mt(BF_C_Q_M5_2D, NF_C_Q_M5_2D, FE_C_Q_M5_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_M6_2D_Q_A:
      return mt(BF_C_Q_M6_2D, NF_C_Q_M6_2D, FE_C_Q_M6_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_M7_2D_Q_A:
      return mt(BF_C_Q_M7_2D, NF_C_Q_M7_2D, FE_C_Q_M7_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_M8_2D_Q_A:
      return mt(BF_C_Q_M8_2D, NF_C_Q_M8_2D, FE_C_Q_M8_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_M9_2D_Q_A:
      return mt(BF_C_Q_M9_2D, NF_C_Q_M9_2D, FE_C_Q_M9_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_M2_2D_Q_M:
      return mt(BF_C_Q_M2_2D, NF_C_Q_M2_2D, FE_C_Q_M2_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_M3_2D_Q_M:
      return mt(BF_C_Q_M3_2D, NF_C_Q_M3_2D, FE_C_Q_M3_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_M4_2D_Q_M:
      return mt(BF_C_Q_M4_2D, NF_C_Q_M4_2D, FE_C_Q_M4_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_M5_2D_Q_M:
      return mt(BF_C_Q_M5_2D, NF_C_Q_M5_2D, FE_C_Q_M5_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_M6_2D_Q_M:
      return mt(BF_C_Q_M6_2D, NF_C_Q_M6_2D, FE_C_Q_M6_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_M7_2D_Q_M:
      return mt(BF_C_Q_M7_2D, NF_C_Q_M7_2D, FE_C_Q_M7_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_M8_2D_Q_M:
      return mt(BF_C_Q_M8_2D, NF_C_Q_M8_2D, FE_C_Q_M8_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_M9_2D_Q_M:
      return mt(BF_C_Q_M9_2D, NF_C_Q_M9_2D, FE_C_Q_M9_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case C_EL1_2D_Q_A:
      return mt(BF_C_Q_EL1_2D, NF_C_Q_EL1_2D, FE_C_Q_EL1_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case C_EL1_2D_Q_M:
      return mt(BF_C_Q_EL1_2D, NF_C_Q_EL1_2D, FE_C_Q_EL1_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case N_RT0_2D_Q_A:
      return mt(BF_N_Q_RT0_2D, NF_N_Q_RT0_2D, FE_N_Q_RT0_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case N_RT0_2D_Q_M:
      return mt(BF_N_Q_RT0_2D, NF_N_Q_RT0_2D, FE_N_Q_RT0_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case N_RT1_2D_Q_A:
      return mt(BF_N_Q_RT1_2D, NF_N_Q_RT1_2D, FE_N_Q_RT1_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case N_RT1_2D_Q_M:
      return mt(BF_N_Q_RT1_2D, NF_N_Q_RT1_2D, FE_N_Q_RT1_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case N_RT2_2D_Q_A:
      return mt(BF_N_Q_RT2_2D, NF_N_Q_RT2_2D, FE_N_Q_RT2_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case N_RT2_2D_Q_M:
      return mt(BF_N_Q_RT2_2D, NF_N_Q_RT2_2D, FE_N_Q_RT2_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case N_RT3_2D_Q_A:
      return mt(BF_N_Q_RT3_2D, NF_N_Q_RT3_2D, FE_N_Q_RT3_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case N_RT3_2D_Q_M:
      return mt(BF_N_Q_RT3_2D, NF_N_Q_RT3_2D, FE_N_Q_RT3_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case N_RT0_2D_T_A:
      return mt(BF_N_T_RT0_2D, NF_N_T_RT0_2D, FE_N_T_RT0_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case N_RT1_2D_T_A:
      return mt(BF_N_T_RT1_2D, NF_N_T_RT1_2D, FE_N_T_RT1_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case N_RT2_2D_T_A:
      return mt(BF_N_T_RT2_2D, NF_N_T_RT2_2D, FE_N_T_RT2_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case N_RT3_2D_T_A:
      return mt(BF_N_T_RT3_2D, NF_N_T_RT3_2D, FE_N_T_RT3_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case N_BDM1_2D_Q_A:
      return mt(BF_N_Q_BDM1_2D, NF_N_Q_BDM1_2D, FE_N_Q_BDM1_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case N_BDM2_2D_Q_A:
      return mt(BF_N_Q_BDM2_2D, NF_N_Q_BDM2_2D, FE_N_Q_BDM2_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case N_BDM3_2D_Q_A:
      return mt(BF_N_Q_BDM3_2D, NF_N_Q_BDM3_2D, FE_N_Q_BDM3_2D, ReferenceTransformation_type::QuadAffin);
      break;
    case N_BDM1_2D_Q_M:
      return mt(BF_N_Q_BDM1_2D, NF_N_Q_BDM1_2D, FE_N_Q_BDM1_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case N_BDM2_2D_Q_M:
      return mt(BF_N_Q_BDM2_2D, NF_N_Q_BDM2_2D, FE_N_Q_BDM2_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case N_BDM3_2D_Q_M:
      return mt(BF_N_Q_BDM3_2D, NF_N_Q_BDM3_2D, FE_N_Q_BDM3_2D, ReferenceTransformation_type::QuadBilinear);
      break;
    case N_BDM1_2D_T_A:
      return mt(BF_N_T_BDM1_2D, NF_N_T_BDM1_2D, FE_N_T_BDM1_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case N_BDM2_2D_T_A:
      return mt(BF_N_T_BDM2_2D, NF_N_T_BDM2_2D, FE_N_T_BDM2_2D, ReferenceTransformation_type::TriaAffin);
      break;
    case N_BDM3_2D_T_A:
      return mt(BF_N_T_BDM3_2D, NF_N_T_BDM3_2D, FE_N_T_BDM3_2D, ReferenceTransformation_type::TriaAffin);
      break;
// 3D
    case C_P00_3D_T_A:
      return mt(BF_C_T_P00_3D, NF_C_T_P00_3D, FE_C_T_P00_3D, ReferenceTransformation_type::TetraAffin);
      break;
    case C_P0_3D_T_A:
      return mt(BF_C_T_P0_3D, NF_C_T_P0_3D, FE_C_T_P0_3D, ReferenceTransformation_type::TetraAffin);
      break;
    case C_P1_3D_T_A:
      return mt(BF_C_T_P1_3D, NF_C_T_P1_3D, FE_C_T_P1_3D, ReferenceTransformation_type::TetraAffin);
      break;
    case C_P2_3D_T_A:
      return mt(BF_C_T_P2_3D, NF_C_T_P2_3D, FE_C_T_P2_3D, ReferenceTransformation_type::TetraAffin);
      break;
    case C_P3_3D_T_A:
      return mt(BF_C_T_P3_3D, NF_C_T_P3_3D, FE_C_T_P3_3D, ReferenceTransformation_type::TetraAffin);
      break;
    case N_P1_3D_T_A:
      return mt(BF_N_T_P1_3D, NF_N_T_P1_3D, FE_N_T_P1_3D, ReferenceTransformation_type::TetraAffin);
      break;
    case C_Q00_3D_H_A:
      return mt(BF_C_H_Q00_3D, NF_C_H_Q00_3D, FE_C_H_Q00_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case C_Q0_3D_H_A:
      return mt(BF_C_H_Q0_3D, NF_C_H_Q0_3D, FE_C_H_Q0_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case C_Q1_3D_H_A:
      return mt(BF_C_H_Q1_3D, NF_C_H_Q1_3D, FE_C_H_Q1_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case C_Q2_3D_H_A:
      return mt(BF_C_H_Q2_3D, NF_C_H_Q2_3D, FE_C_H_Q2_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case C_Q3_3D_H_A:
      return mt(BF_C_H_Q3_3D, NF_C_H_Q3_3D, FE_C_H_Q3_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case C_Q4_3D_H_A:
      return mt(BF_C_H_Q4_3D, NF_C_H_Q4_3D, FE_C_H_Q4_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case C_Q00_3D_H_M:
      return mt(BF_C_H_Q00_3D, NF_C_H_Q00_3D, FE_C_H_Q00_3D, ReferenceTransformation_type::HexaTrilinear);
      break;
    case C_Q0_3D_H_M:
      return mt(BF_C_H_Q0_3D, NF_C_H_Q0_3D, FE_C_H_Q0_3D, ReferenceTransformation_type::HexaTrilinear);
      break;
    case C_Q1_3D_H_M:
      return mt(BF_C_H_Q1_3D, NF_C_H_Q1_3D, FE_C_H_Q1_3D, ReferenceTransformation_type::HexaTrilinear);
      break;
    case C_Q2_3D_H_M:
      return mt(BF_C_H_Q2_3D, NF_C_H_Q2_3D, FE_C_H_Q2_3D, ReferenceTransformation_type::HexaTrilinear);
      break;
    case C_Q3_3D_H_M:
      return mt(BF_C_H_Q3_3D, NF_C_H_Q3_3D, FE_C_H_Q3_3D, ReferenceTransformation_type::HexaTrilinear);
      break;
    case C_Q4_3D_H_M:
      return mt(BF_C_H_Q4_3D, NF_C_H_Q4_3D, FE_C_H_Q4_3D, ReferenceTransformation_type::HexaTrilinear);
      break;
    case N_Q1_3D_H_A:
      return mt(BF_N_H_Q1_3D, NF_N_H_Q1_3D, FE_N_H_Q1_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case N_Q1_3D_H_M:
      return mt(BF_N_H_Q1_3D, NF_N_H_Q1_3D, FE_N_H_Q1_3D, ReferenceTransformation_type::HexaTrilinear);
      break;
    case D_P1_3D_H_A:
      return mt(BF_D_H_P1_3D, NF_D_H_P1_3D, FE_D_H_P1_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case D_P2_3D_H_A:
      return mt(BF_D_H_P2_3D, NF_D_H_P2_3D, FE_D_H_P2_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case D_P3_3D_H_A:
      return mt(BF_D_H_P3_3D, NF_D_H_P3_3D, FE_D_H_P3_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case D_P1_3D_H_M:
      return mt(BF_D_H_P1_3D, NF_D_H_P1_3D, FE_D_H_P1_3D, ReferenceTransformation_type::HexaTrilinear);
      break;
    case D_P2_3D_H_M:
      return mt(BF_D_H_P2_3D, NF_D_H_P2_3D, FE_D_H_P2_3D, ReferenceTransformation_type::HexaTrilinear);
      break;
    case D_P3_3D_H_M:
      return mt(BF_D_H_P3_3D, NF_D_H_P3_3D, FE_D_H_P3_3D, ReferenceTransformation_type::HexaTrilinear);
      break;
    case D_Q1_3D_H_A:
      return mt(BF_D_H_Q1_3D, NF_D_H_Q1_3D, FE_D_H_Q1_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case D_Q1_3D_H_M:
      return mt(BF_D_H_Q1_3D, NF_D_H_Q1_3D, FE_D_H_Q1_3D, ReferenceTransformation_type::HexaTrilinear);
      break;
    case D_Q2_3D_H_A:
      return mt(BF_D_H_Q2_3D, NF_D_H_Q2_3D, FE_D_H_Q2_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case D_Q2_3D_H_M:
      return mt(BF_D_H_Q2_3D, NF_D_H_Q2_3D, FE_D_H_Q2_3D, ReferenceTransformation_type::HexaTrilinear);
      break;
    case B_IB2_3D_H_A:
      return mt(BF_B_H_IB2_3D, NF_B_H_IB2_3D, FE_B_H_IB2_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case B_IB2_3D_H_M:
      return mt(BF_B_H_IB2_3D, NF_B_H_IB2_3D, FE_B_H_IB2_3D, ReferenceTransformation_type::HexaTrilinear);
      break;
    case N_P2_3D_T_A:
      return mt(BF_N_T_P2_3D, NF_N_T_P2_3D, FE_N_T_P2_3D, ReferenceTransformation_type::TetraAffin);
      break;
    case N_P3_3D_T_A:
      ErrThrow("finite element of type N_P3_3D_T_A not implemented")
      break;
    case N_P4_3D_T_A:
      ErrThrow("finite element of type N_P4_3D_T_A not implemented")
      break;
    case N_P5_3D_T_A:
      ErrThrow("finite element of type N_P5_3D_T_A not implemented")
      break;
    case N_Q2_3D_H_A:
      return mt(BF_N_H_Q2_3D, NF_N_H_Q2_3D, FE_N_H_Q2_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case N_Q3_3D_H_A:
      return mt(BF_N_H_Q3_3D, NF_N_H_Q3_3D, FE_N_H_Q3_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case N_Q4_3D_H_A:
      return mt(BF_N_H_Q4_3D, NF_N_H_Q4_3D, FE_N_H_Q4_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case N_Q5_3D_H_A:
      ErrThrow("finite element of type N_Q5_3D_H_A not implemented")
      break;
    case N_Q2_3D_H_M:
      return mt(BF_N_H_Q2_3D, NF_N_H_Q2_3D, FE_N_H_Q2_3D, ReferenceTransformation_type::HexaTrilinear);
      break;
    case N_Q3_3D_H_M:
      return mt(BF_N_H_Q3_3D, NF_N_H_Q3_3D, FE_N_H_Q3_3D, ReferenceTransformation_type::HexaTrilinear);
      break;
    case N_Q4_3D_H_M:
      return mt(BF_N_H_Q4_3D, NF_N_H_Q4_3D, FE_N_H_Q4_3D, ReferenceTransformation_type::HexaTrilinear);
      break;
    case N_Q5_3D_H_M:
      ErrThrow("finite element of type N_Q5_3D_H_M not implemented")
      break;
    case C_B2_3D_T_A:
      return mt(BF_C_T_B2_3D, NF_C_T_B2_3D, FE_C_T_B2_3D, ReferenceTransformation_type::TetraAffin);
      break;
    case D_P1_3D_T_A:
      return mt(BF_D_T_P1_3D, NF_D_T_P1_3D, FE_D_T_P1_3D, ReferenceTransformation_type::TetraAffin);
      break;
    case D_P2_3D_T_A:
      return mt(BF_D_T_P2_3D, NF_D_T_P2_3D, FE_D_T_P2_3D, ReferenceTransformation_type::TetraAffin);
      break;
    case D_P3_3D_T_A:
      return mt(BF_D_T_P3_3D, NF_D_T_P3_3D, FE_D_T_P3_3D, ReferenceTransformation_type::TetraAffin);
      break;
    case C_UL1_3D_H_A:
      return mt(BF_C_H_UL1_3D, NF_C_H_UL1_3D, FE_C_H_UL1_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case C_UL2_3D_H_A:
      return mt(BF_C_H_UL2_3D, NF_C_H_UL2_3D, FE_C_H_UL2_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case C_UL3_3D_H_A:
      return mt(BF_C_H_UL3_3D, NF_C_H_UL3_3D, FE_C_H_UL3_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case C_UL4_3D_H_A:
      ErrThrow("finite element of type C_UL4_3D_H_A not implemented")
      break;
    case C_UL5_3D_H_A:
      ErrThrow("finite element of type C_UL5_3D_H_A not implemented")
      break;
    case C_UL1_3D_H_M:
      return mt(BF_C_H_UL1_3D, NF_C_H_UL1_3D, FE_C_H_UL1_3D, ReferenceTransformation_type::HexaTrilinear);
      break;
    case C_UL2_3D_H_M:
      return mt(BF_C_H_UL2_3D, NF_C_H_UL2_3D, FE_C_H_UL2_3D, ReferenceTransformation_type::HexaTrilinear);
      break;
    case C_UL3_3D_H_M:
      return mt(BF_C_H_UL3_3D, NF_C_H_UL3_3D, FE_C_H_UL3_3D, ReferenceTransformation_type::HexaTrilinear);
      break;
    case C_UL4_3D_H_M:
      ErrThrow("finite element of type C_UL4_3D_H_M not implemented")
      break;
    case C_UL5_3D_H_M:
      ErrThrow("finite element of type C_UL5_3D_H_M not implemented")
      break;
    case C_UL1_3D_T_A:
      return mt(BF_C_H_UL1_3D, NF_C_H_UL1_3D, FE_C_H_UL1_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case C_UL2_3D_T_A:
      return mt(BF_C_H_UL2_3D, NF_C_H_UL2_3D, FE_C_H_UL2_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case C_UL3_3D_T_A:
      return mt(BF_C_H_UL3_3D, NF_C_H_UL3_3D, FE_C_H_UL3_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case C_UL4_3D_T_A:
      ErrThrow("finite element of type C_UL4_3D_T_A not implemented")
      break;
    case C_UL5_3D_T_A:
      ErrThrow("finite element of type C_UL5_3D_T_A not implemented")
      break;
    case N_RT0_3D_T_A:
      return mt(BF_N_T_RT0_3D, NF_N_T_RT0_3D, FE_N_T_RT0_3D, ReferenceTransformation_type::TetraAffin);
      break;
    case N_RT1_3D_T_A:
      return mt(BF_N_T_RT1_3D, NF_N_T_RT1_3D, FE_N_T_RT1_3D, ReferenceTransformation_type::TetraAffin);
      break;
    case N_RT2_3D_T_A:
      return mt(BF_N_T_RT2_3D, NF_N_T_RT2_3D, FE_N_T_RT2_3D, ReferenceTransformation_type::TetraAffin);
      break;
    case N_RT3_3D_T_A:
      return mt(BF_N_T_RT3_3D, NF_N_T_RT3_3D, FE_N_T_RT3_3D, ReferenceTransformation_type::TetraAffin);
      break;
    case N_RT0_3D_H_A:
      return mt(BF_N_H_RT0_3D, NF_N_H_RT0_3D, FE_N_H_RT0_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case N_RT1_3D_H_A:
      return mt(BF_N_H_RT1_3D, NF_N_H_RT1_3D, FE_N_H_RT1_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case N_RT2_3D_H_A:
      return mt(BF_N_H_RT2_3D, NF_N_H_RT2_3D, FE_N_H_RT2_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case N_BDDF1_3D_T_A:
      return mt(BF_N_T_BDDF1_3D, NF_N_T_BDDF1_3D, FE_N_T_BDDF1_3D, ReferenceTransformation_type::TetraAffin);
      break;
    case N_BDDF2_3D_T_A:
      return mt(BF_N_T_BDDF2_3D, NF_N_T_BDDF2_3D, FE_N_T_BDDF2_3D, ReferenceTransformation_type::TetraAffin);
      break;
    case N_BDDF3_3D_T_A:
      return mt(BF_N_T_BDDF3_3D, NF_N_T_BDDF3_3D, FE_N_T_BDDF3_3D, ReferenceTransformation_type::TetraAffin);
      break;
    case N_BDDF1_3D_H_A:
      return mt(BF_N_H_BDDF1_3D, NF_N_H_BDDF1_3D, FE_N_H_BDDF1_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case N_BDDF2_3D_H_A:
      return mt(BF_N_H_BDDF2_3D, NF_N_H_BDDF2_3D, FE_N_H_BDDF2_3D, ReferenceTransformation_type::HexaAffin);
      break;
    case N_BDDF3_3D_H_A:
      return mt(BF_N_H_BDDF3_3D, NF_N_H_BDDF3_3D, FE_N_H_BDDF3_3D, ReferenceTransformation_type::HexaAffin);
      break;
    default:
      ErrThrow("unknown FE3D ", id);
      break;
  }
}

FiniteElement::FiniteElement ( FE_type id)
:  FiniteElement (id, get_ids(id))
{
}

FiniteElement::FiniteElement (
  FE_type id, std::tuple<BaseFunction_type, NodalFunctional_type, FEDescriptor_type, ReferenceTransformation_type> ids)
 : type(id), base_function(std::get<0>(ids)),
   nodal_functional(std::get<1>(ids)), fe_descriptor(std::get<2>(ids)),
   reference_transformation_type(std::get<3>(ids))
{
}

void FiniteElement::CheckNFandBF() const
{
#ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank != 0)
    return; // only do this on one processor
#endif
  Output::print<3>("CheckNFandBF: BaseFunct_ID: ", this->GetBaseFunct_ID(),
                   " NodalFunctional_ID: ", this->GetNodalFunctional_ID());
  int N_Points;
  const double *xi, *eta, *zeta;
  int d = this->get_space_dim();
  if(d == 3)
    nodal_functional.GetPointsForAll(N_Points, xi, eta, zeta);
  else if(d == 2)
    nodal_functional.GetPointsForAll(N_Points, xi, eta);
  else
    ErrThrow("The method 'FiniteElement::CheckNFandBF()' is not implemented "
             "for a 1D element");
  
  int n_basis_functions = base_function.GetDimension();
  int baseVectDim = base_function.GetBaseVectDim();
  std::vector<std::vector<double>> AllPointValues(
    N_Points, std::vector<double>(n_basis_functions * baseVectDim));
  for(int k=0;k<N_Points;k++)
  {
    if(d == 2)
      base_function.GetDerivatives(MultiIndex2D::D00, xi[k], eta[k],
                                   AllPointValues[k].data());
    else
      base_function.GetDerivatives(MultiIndex3D::D000, xi[k], eta[k], zeta[k],
                                   AllPointValues[k].data());
  }
  std::vector<double> PointValues(N_Points*baseVectDim);
  std::vector<double> FunctionalValues(n_basis_functions);
  for(int k = 0; k < n_basis_functions; k++)
  {
    for(int l = 0; l < N_Points; l++)
    {
      for(int i = 0; i< baseVectDim; i++)
      {
        PointValues[l + N_Points * i]
          = AllPointValues[l][k+n_basis_functions*i];
      }
    }
    nodal_functional.GetAllFunctionals(nullptr, nullptr, PointValues.data(),
                                       FunctionalValues.data());

    for(int i = 0; i < n_basis_functions; i++)
    {
      if(std::abs(FunctionalValues[i])<1e-10)
      {
        FunctionalValues[i] = 0;
      }
      if( i == k && std::abs(FunctionalValues[i]-1) > 1e-8 )
      {
        Output::print<3>("BF: ", k, " NF: ", i, " ", FunctionalValues[i]);
      }
      if( i != k && std::abs(FunctionalValues[i]-0) > 1e-8 )
      {
        Output::print<3>("BF: ", k, " NF: ", i, " ", FunctionalValues[i]);
      }
    }
  }
}

int FiniteElement::get_space_dim() const
{
  if(this->type < N_FEs1D)
    return 1;
  else if(this->type < N_FEs1D + N_FEs2D)
    return 2;
  else if(this->type < N_FEs1D + N_FEs2D + N_FEs3D)
    return 3;
  else
    ErrThrow("do not know the space dimension of this finite element, type ",
             this->type);
}
