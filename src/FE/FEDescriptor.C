#include "FEDescriptor.h"
#include "MooNMD_Io.h"
#include "FE_C_L_P0_1D.h"
#include "FE_C_L_P1_1D.h"
#include "FE_C_L_P2_1D.h"
#include "FE_C_L_P3_1D.h"
#include "FE_N_L_P0_1D.h"
#include "FE_D_L_P1_1D.h"
#include "FE_D_L_P2_1D.h"
#include "AllFEDescs2D.h"
#include "AllFEDescs3D.h"

FEDescriptor::FEDescriptor (FEDescriptor_type id)
{
  type = id;
  N_OuterDOF = 0;
  OuterDOF = nullptr;
  N_EdgeDOF = 0;
  EdgeDOF = nullptr;
  N_VertDOF = 0;
  VertDOF = nullptr;
  EdgeVertData_Filled = false;
  switch(type)
  {
    case FE_C_L_P0_1D:
      N_DOF = C_L_P0_1D_NDOF;
      N_JointDOF = C_L_P0_1D_JointDOF;
      JointDOF = C_L_P0_1D_J;
      N_InnerDOF = C_L_P0_1D_NInner;
      InnerDOF = C_L_P0_1D_Inner;
      break;
    case FE_C_L_P1_1D:
      N_DOF = C_L_P1_1D_NDOF;
      N_JointDOF = C_L_P1_1D_JointDOF;
      JointDOF = C_L_P1_1D_J;
      N_InnerDOF = C_L_P1_1D_NInner;
      InnerDOF = C_L_P1_1D_Inner;
      break;
    case FE_C_L_P2_1D:
      N_DOF = C_L_P2_1D_NDOF;
      N_JointDOF = C_L_P2_1D_JointDOF;
      JointDOF = C_L_P2_1D_J;
      N_InnerDOF = C_L_P2_1D_NInner;
      InnerDOF = C_L_P2_1D_Inner;
      break;
    case FE_C_L_P3_1D:
      N_DOF = C_L_P3_1D_NDOF;
      N_JointDOF = C_L_P3_1D_JointDOF;
      JointDOF = C_L_P3_1D_J;
      N_InnerDOF = C_L_P3_1D_NInner;
      InnerDOF = C_L_P3_1D_Inner;
      break;
    case FE_N_L_P0_1D:
      N_DOF = N_L_P0_1D_NDOF;
      N_JointDOF = N_L_P0_1D_JointDOF;
      JointDOF = N_L_P0_1D_J;
      N_InnerDOF = N_L_P0_1D_NInner;
      InnerDOF = N_L_P0_1D_Inner;
      break;
    case FE_D_L_P1_1D:
      N_DOF = D_L_P1_1D_NDOF;
      N_JointDOF = D_L_P1_1D_JointDOF;
      JointDOF = D_L_P1_1D_J;
      N_InnerDOF = D_L_P1_1D_NInner;
      InnerDOF = D_L_P1_1D_Inner;
      break;
    case FE_D_L_P2_1D:
      N_DOF = D_L_P2_1D_NDOF;
      N_JointDOF = D_L_P2_1D_JointDOF;
      JointDOF = D_L_P2_1D_J;
      N_InnerDOF = D_L_P2_1D_NInner;
      InnerDOF = D_L_P2_1D_Inner;
      break;
    case FE_C_T_P00_2D:
      N_DOF = C_T_P00_2D_NDOF;
      N_JointDOF = C_T_P00_2D_JointDOF;
      JointDOF = C_T_P00_2D_J;
      N_InnerDOF = C_T_P00_2D_NInner;
      InnerDOF = C_T_P00_2D_Inner;
      break;
    case FE_C_T_P0_2D:
      N_DOF = C_T_P0_2D_NDOF;
      N_JointDOF = C_T_P0_2D_JointDOF;
      JointDOF = C_T_P0_2D_J;
      N_InnerDOF = C_T_P0_2D_NInner;
      InnerDOF = C_T_P0_2D_Inner;
      break;
    case FE_C_T_P1_2D:
      N_DOF = C_T_P1_2D_NDOF;
      N_JointDOF = C_T_P1_2D_JointDOF;
      JointDOF = C_T_P1_2D_J;
      N_InnerDOF = C_T_P1_2D_NInner;
      InnerDOF = C_T_P1_2D_Inner;
      break;
    case FE_C_T_P2_2D:
      N_DOF = C_T_P2_2D_NDOF;
      N_JointDOF = C_T_P2_2D_JointDOF;
      JointDOF = C_T_P2_2D_J;
      N_InnerDOF = C_T_P2_2D_NInner;
      InnerDOF = C_T_P2_2D_Inner;
      break;
    case FE_C_T_P3_2D:
      N_DOF = C_T_P3_2D_NDOF;
      N_JointDOF = C_T_P3_2D_JointDOF;
      JointDOF = C_T_P3_2D_J;
      N_InnerDOF = C_T_P3_2D_NInner;
      InnerDOF = C_T_P3_2D_Inner;
      break;
    case FE_C_T_P4_2D:
      N_DOF = C_T_P4_2D_NDOF;
      N_JointDOF = C_T_P4_2D_JointDOF;
      JointDOF = C_T_P4_2D_J;
      N_InnerDOF = C_T_P4_2D_NInner;
      InnerDOF = C_T_P4_2D_Inner;
      N_OuterDOF = C_T_P4_2D_NOuter;
      OuterDOF = C_T_P4_2D_Outer;
      break;
    case FE_C_T_P5_2D:
      N_DOF = C_T_P5_2D_NDOF;
      N_JointDOF = C_T_P5_2D_JointDOF;
      JointDOF = C_T_P5_2D_J;
      N_InnerDOF = C_T_P5_2D_NInner;
      InnerDOF = C_T_P5_2D_Inner;
      N_OuterDOF = C_T_P5_2D_NOuter;
      OuterDOF = C_T_P5_2D_Outer;
      break;
    case FE_C_T_P6_2D:
      N_DOF = C_T_P6_2D_NDOF;
      N_JointDOF = C_T_P6_2D_JointDOF;
      JointDOF = C_T_P6_2D_J;
      N_InnerDOF = C_T_P6_2D_NInner;
      InnerDOF = C_T_P6_2D_Inner;
      N_OuterDOF = C_T_P6_2D_NOuter;
      OuterDOF = C_T_P6_2D_Outer;
      break;
    case FE_C_T_P7_2D:
      ErrThrow("FE descriptor of type FE_C_T_P7_2D not implemented");
      break;
    case FE_C_T_P8_2D:
      ErrThrow("FE descriptor of type FE_C_T_P8_2D not implemented");
      break;
    case FE_C_T_P9_2D:
      ErrThrow("FE descriptor of type FE_C_T_P9_2D not implemented");
      break;
    case FE_N_T_P1_2D:
      N_DOF = N_T_P1_2D_NDOF;
      N_JointDOF = N_T_P1_2D_JointDOF;
      JointDOF = N_T_P1_2D_J;
      N_InnerDOF = N_T_P1_2D_NInner;
      InnerDOF = N_T_P1_2D_Inner;
      N_OuterDOF = N_T_P1_2D_NOuter;
      OuterDOF = N_T_P1_2D_Outer;
      break;
    case FE_C_Q_Q00_2D:
      N_DOF = C_Q_Q00_2D_NDOF;
      N_JointDOF = C_Q_Q00_2D_JointDOF;
      JointDOF = C_Q_Q00_2D_J;
      N_InnerDOF = C_Q_Q00_2D_NInner;
      InnerDOF = C_Q_Q00_2D_Inner;
      break;
    case FE_C_Q_Q0_2D:
      N_DOF = C_Q_Q0_2D_NDOF;
      N_JointDOF = C_Q_Q0_2D_JointDOF;
      JointDOF = C_Q_Q0_2D_J;
      N_InnerDOF = C_Q_Q0_2D_NInner;
      InnerDOF = C_Q_Q0_2D_Inner;
      break;
    case FE_C_Q_Q1_2D:
      N_DOF = C_Q_Q1_2D_NDOF;
      N_JointDOF = C_Q_Q1_2D_JointDOF;
      JointDOF = C_Q_Q1_2D_J;
      N_InnerDOF = C_Q_Q1_2D_NInner;
      InnerDOF = C_Q_Q1_2D_Inner;
      N_OuterDOF = C_Q_Q1_2D_NOuter;
      OuterDOF = C_Q_Q1_2D_Outer;
      break;
    case FE_C_Q_Q2_2D:
      N_DOF = C_Q_Q2_2D_NDOF;
      N_JointDOF = C_Q_Q2_2D_JointDOF;
      JointDOF = C_Q_Q2_2D_J;
      N_InnerDOF = C_Q_Q2_2D_NInner;
      InnerDOF = C_Q_Q2_2D_Inner;
      N_OuterDOF = C_Q_Q2_2D_NOuter;
      OuterDOF = C_Q_Q2_2D_Outer;
      break;
    case FE_C_Q_Q3_2D:
      N_DOF = C_Q_Q3_2D_NDOF;
      N_JointDOF = C_Q_Q3_2D_JointDOF;
      JointDOF = C_Q_Q3_2D_J;
      N_InnerDOF = C_Q_Q3_2D_NInner;
      InnerDOF = C_Q_Q3_2D_Inner;
      N_OuterDOF = C_Q_Q3_2D_NOuter;
      OuterDOF = C_Q_Q3_2D_Outer;
      break;
    case FE_C_Q_Q4_2D:
      N_DOF = C_Q_Q4_2D_NDOF;
      N_JointDOF = C_Q_Q4_2D_JointDOF;
      JointDOF = C_Q_Q4_2D_J;
      N_InnerDOF = C_Q_Q4_2D_NInner;
      InnerDOF = C_Q_Q4_2D_Inner;
      N_OuterDOF = C_Q_Q4_2D_NOuter;
      OuterDOF = C_Q_Q4_2D_Outer;
      break;
    case FE_C_Q_Q5_2D:
      N_DOF = C_Q_Q5_2D_NDOF;
      N_JointDOF = C_Q_Q5_2D_JointDOF;
      JointDOF = C_Q_Q5_2D_J;
      N_InnerDOF = C_Q_Q5_2D_NInner;
      InnerDOF = C_Q_Q5_2D_Inner;
      N_OuterDOF = C_Q_Q5_2D_NOuter;
      OuterDOF = C_Q_Q5_2D_Outer;
      break;
    case FE_C_Q_Q6_2D:
      N_DOF = C_Q_Q6_2D_NDOF;
      N_JointDOF = C_Q_Q6_2D_JointDOF;
      JointDOF = C_Q_Q6_2D_J;
      N_InnerDOF = C_Q_Q6_2D_NInner;
      InnerDOF = C_Q_Q6_2D_Inner;
      N_OuterDOF = C_Q_Q6_2D_NOuter;
      OuterDOF = C_Q_Q6_2D_Outer;
      break;
    case FE_C_Q_Q7_2D:
      N_DOF = C_Q_Q7_2D_NDOF;
      N_JointDOF = C_Q_Q7_2D_JointDOF;
      JointDOF = C_Q_Q7_2D_J;
      N_InnerDOF = C_Q_Q7_2D_NInner;
      InnerDOF = C_Q_Q7_2D_Inner;
      N_OuterDOF = C_Q_Q7_2D_NOuter;
      OuterDOF = C_Q_Q7_2D_Outer;
      break;
    case FE_C_Q_Q8_2D:
      N_DOF = C_Q_Q8_2D_NDOF;
      N_JointDOF = C_Q_Q8_2D_JointDOF;
      JointDOF = C_Q_Q8_2D_J;
      N_InnerDOF = C_Q_Q8_2D_NInner;
      InnerDOF = C_Q_Q8_2D_Inner;
      N_OuterDOF = C_Q_Q8_2D_NOuter;
      OuterDOF = C_Q_Q8_2D_Outer;
      break;
    case FE_C_Q_Q9_2D:
      N_DOF = C_Q_Q9_2D_NDOF;
      N_JointDOF = C_Q_Q9_2D_JointDOF;
      JointDOF = C_Q_Q9_2D_J;
      N_InnerDOF = C_Q_Q9_2D_NInner;
      InnerDOF = C_Q_Q9_2D_Inner;
      N_OuterDOF = C_Q_Q9_2D_NOuter;
      OuterDOF = C_Q_Q9_2D_Outer;
      break;
    case FE_N_Q_Q1_2D:
      N_DOF = N_Q_Q1_2D_NDOF;
      N_JointDOF = N_Q_Q1_2D_JointDOF;
      JointDOF = N_Q_Q1_2D_J;
      N_InnerDOF = N_Q_Q1_2D_NInner;
      InnerDOF = N_Q_Q1_2D_Inner;
      N_OuterDOF = N_Q_Q1_2D_NOuter;
      OuterDOF = N_Q_Q1_2D_Outer;
      break;
    case FE_D_Q_P1_2D:
      N_DOF = D_Q_P1_2D_NDOF;
      N_JointDOF = D_Q_P1_2D_JointDOF;
      JointDOF = D_Q_P1_2D_J;
      N_InnerDOF = D_Q_P1_2D_NInner;
      InnerDOF = D_Q_P1_2D_Inner;
      break;
    case FE_D_Q_P2_2D:
      N_DOF = D_Q_P2_2D_NDOF;
      N_JointDOF = D_Q_P2_2D_JointDOF;
      JointDOF = D_Q_P2_2D_J;
      N_InnerDOF = D_Q_P2_2D_NInner;
      InnerDOF = D_Q_P2_2D_Inner;
      break;
    case FE_D_Q_P3_2D:
      N_DOF = D_Q_P3_2D_NDOF;
      N_JointDOF = D_Q_P3_2D_JointDOF;
      JointDOF = D_Q_P3_2D_J;
      N_InnerDOF = D_Q_P3_2D_NInner;
      InnerDOF = D_Q_P3_2D_Inner;
      break;
    case FE_C_T_B2_2D:
      N_DOF = C_T_B2_2D_NDOF;
      N_JointDOF = C_T_B2_2D_JointDOF;
      JointDOF = C_T_B2_2D_J;
      N_InnerDOF = C_T_B2_2D_NInner;
      InnerDOF = C_T_B2_2D_Inner;
      N_OuterDOF = C_T_B2_2D_NOuter;
      OuterDOF = C_T_B2_2D_Outer;
      break;
    case FE_C_T_B3_2D:
      N_DOF = C_T_B3_2D_NDOF;
      N_JointDOF = C_T_B3_2D_JointDOF;
      JointDOF = C_T_B3_2D_J;
      N_InnerDOF = C_T_B3_2D_NInner;
      InnerDOF = C_T_B3_2D_Inner;
      N_OuterDOF = C_T_B3_2D_NOuter;
      OuterDOF = C_T_B3_2D_Outer;
      break;
    case FE_C_T_SV2_2D:
      N_DOF = C_T_SV2_2D_NDOF;
      N_JointDOF = C_T_SV2_2D_JointDOF;
      JointDOF = C_T_SV2_2D_J;
      N_InnerDOF = C_T_SV2_2D_NInner;
      InnerDOF = C_T_SV2_2D_Inner;
      break;
    case FE_D_T_P1_2D:
      N_DOF = D_T_P1_2D_NDOF;
      N_JointDOF = D_T_P1_2D_JointDOF;
      JointDOF = D_T_P1_2D_J;
      N_InnerDOF = D_T_P1_2D_NInner;
      InnerDOF = D_T_P1_2D_Inner;
      break;
    case FE_D_T_P2_2D:
      N_DOF = D_T_P2_2D_NDOF;
      N_JointDOF = D_T_P2_2D_JointDOF;
      JointDOF = D_T_P2_2D_J;
      N_InnerDOF = D_T_P2_2D_NInner;
      InnerDOF = D_T_P2_2D_Inner;
      break;
    case FE_D_T_SV1_2D:
      N_DOF = D_T_SV1_2D_NDOF;
      N_JointDOF = D_T_SV1_2D_JointDOF;
      JointDOF = D_T_SV1_2D_J;
      N_InnerDOF = D_T_SV1_2D_NInner;
      InnerDOF = D_T_SV1_2D_Inner;
      break;
    case FE_N_Q_Q2_2D:
      N_DOF = N_Q_Q2_2D_NDOF;
      N_JointDOF = N_Q_Q2_2D_JointDOF;
      JointDOF = N_Q_Q2_2D_J;
      N_InnerDOF = N_Q_Q2_2D_NInner;
      InnerDOF = N_Q_Q2_2D_Inner;
      N_OuterDOF = N_Q_Q2_2D_NOuter;
      OuterDOF = N_Q_Q2_2D_Outer;
      break;
    case FE_N_Q_Q3_2D:
      N_DOF = N_Q_Q3_2D_NDOF;
      N_JointDOF = N_Q_Q3_2D_JointDOF;
      JointDOF = N_Q_Q3_2D_J;
      N_InnerDOF = N_Q_Q3_2D_NInner;
      InnerDOF = N_Q_Q3_2D_Inner;
      N_OuterDOF = N_Q_Q3_2D_NOuter;
      OuterDOF = N_Q_Q3_2D_Outer;
      break;
    case FE_N_Q_Q4_2D:
      N_DOF = N_Q_Q4_2D_NDOF;
      N_JointDOF = N_Q_Q4_2D_JointDOF;
      JointDOF = N_Q_Q4_2D_J;
      N_InnerDOF = N_Q_Q4_2D_NInner;
      InnerDOF = N_Q_Q4_2D_Inner;
      N_OuterDOF = N_Q_Q4_2D_NOuter;
      OuterDOF = N_Q_Q4_2D_Outer;
      break;
    case FE_N_Q_Q5_2D:
      N_DOF = N_Q_Q5_2D_NDOF;
      N_JointDOF = N_Q_Q5_2D_JointDOF;
      JointDOF = N_Q_Q5_2D_J;
      N_InnerDOF = N_Q_Q5_2D_NInner;
      InnerDOF = N_Q_Q5_2D_Inner;
      N_OuterDOF = N_Q_Q5_2D_NOuter;
      OuterDOF = N_Q_Q5_2D_Outer;
      break;
    case FE_D_Q_P4_2D:
      N_DOF = D_Q_P4_2D_NDOF;
      N_JointDOF = D_Q_P4_2D_JointDOF;
      JointDOF = D_Q_P4_2D_J;
      N_InnerDOF = D_Q_P4_2D_NInner;
      InnerDOF = D_Q_P4_2D_Inner;
      break;
    case FE_D_Q_P5_2D:
      N_DOF = D_Q_P5_2D_NDOF;
      N_JointDOF = D_Q_P5_2D_JointDOF;
      JointDOF = D_Q_P5_2D_J;
      N_InnerDOF = D_Q_P5_2D_NInner;
      InnerDOF = D_Q_P5_2D_Inner;
      break;
    case FE_D_Q_P6_2D:
      N_DOF = D_Q_P6_2D_NDOF;
      N_JointDOF = D_Q_P6_2D_JointDOF;
      JointDOF = D_Q_P6_2D_J;
      N_InnerDOF = D_Q_P6_2D_NInner;
      InnerDOF = D_Q_P6_2D_Inner;
      break;
    case FE_D_Q_P7_2D:
      N_DOF = D_Q_P7_2D_NDOF;
      N_JointDOF = D_Q_P7_2D_JointDOF;
      JointDOF = D_Q_P7_2D_J;
      N_InnerDOF = D_Q_P7_2D_NInner;
      InnerDOF = D_Q_P7_2D_Inner;
      break;
    case FE_N_T_P1MOD_2D:
      N_DOF = N_T_P1MOD_2D_NDOF;
      N_JointDOF = N_T_P1MOD_2D_JointDOF;
      JointDOF = N_T_P1MOD_2D_J;
      N_InnerDOF = N_T_P1MOD_2D_NInner;
      InnerDOF = N_T_P1MOD_2D_Inner;
      break;
    case FE_N_T_P2_2D:
      N_DOF = N_T_P2_2D_NDOF;
      N_JointDOF = N_T_P2_2D_JointDOF;
      JointDOF = N_T_P2_2D_J;
      N_InnerDOF = N_T_P2_2D_NInner;
      InnerDOF = N_T_P2_2D_Inner;
      N_OuterDOF = N_T_P2_2D_NOuter;
      OuterDOF = N_T_P2_2D_Outer;
      break;
    case FE_N_T_P3_2D:
      N_DOF = N_T_P3_2D_NDOF;
      N_JointDOF = N_T_P3_2D_JointDOF;
      JointDOF = N_T_P3_2D_J;
      N_InnerDOF = N_T_P3_2D_NInner;
      InnerDOF = N_T_P3_2D_Inner;
      N_OuterDOF = N_T_P3_2D_NOuter;
      OuterDOF = N_T_P3_2D_Outer;
      break;
    case FE_N_T_P4_2D:
      N_DOF = N_T_P4_2D_NDOF;
      N_JointDOF = N_T_P4_2D_JointDOF;
      JointDOF = N_T_P4_2D_J;
      N_InnerDOF = N_T_P4_2D_NInner;
      InnerDOF = N_T_P4_2D_Inner;
      N_OuterDOF = N_T_P4_2D_NOuter;
      OuterDOF = N_T_P4_2D_Outer;
      break;
    case FE_N_T_P5_2D:
      N_DOF = N_T_P5_2D_NDOF;
      N_JointDOF = N_T_P5_2D_JointDOF;
      JointDOF = N_T_P5_2D_J;
      N_InnerDOF = N_T_P5_2D_NInner;
      InnerDOF = N_T_P5_2D_Inner;
      N_OuterDOF = N_T_P5_2D_NOuter;
      OuterDOF = N_T_P5_2D_Outer;
      break;
    case FE_D_T_P3_2D:
      N_DOF = D_T_P3_2D_NDOF;
      N_JointDOF = D_T_P3_2D_JointDOF;
      JointDOF = D_T_P3_2D_J;
      N_InnerDOF = D_T_P3_2D_NInner;
      InnerDOF = D_T_P3_2D_Inner;
      break;
    case FE_D_T_P4_2D:
      N_DOF = D_T_P4_2D_NDOF;
      N_JointDOF = D_T_P4_2D_JointDOF;
      JointDOF = D_T_P4_2D_J;
      N_InnerDOF = D_T_P4_2D_NInner;
      InnerDOF = D_T_P4_2D_Inner;
      break;
    case FE_C_T_P1MINI_2D:
      N_DOF = C_T_P1MINI_2D_NDOF;
      N_JointDOF = C_T_P1MINI_2D_JointDOF;
      JointDOF = C_T_P1MINI_2D_J;
      N_InnerDOF = C_T_P1MINI_2D_NInner;
      InnerDOF = C_T_P1MINI_2D_Inner;
      break;
    case FE_B_Q_IB2_2D:
      N_DOF = B_Q_IB2_2D_NDOF;
      N_JointDOF = B_Q_IB2_2D_JointDOF;
      JointDOF = B_Q_IB2_2D_J;
      N_InnerDOF = B_Q_IB2_2D_NInner;
      InnerDOF = B_Q_IB2_2D_Inner;
      break;
    case FE_D_Q_Q1_2D:
      N_DOF = D_Q_Q1_2D_NDOF;
      N_JointDOF = D_Q_Q1_2D_JointDOF;
      JointDOF = D_Q_Q1_2D_J;
      N_InnerDOF = D_Q_Q1_2D_NInner;
      InnerDOF = D_Q_Q1_2D_Inner;
      break;
    case FE_D_Q_Q2_2D:
      N_DOF = D_Q_Q2_2D_NDOF;
      N_JointDOF = D_Q_Q2_2D_JointDOF;
      JointDOF = D_Q_Q2_2D_J;
      N_InnerDOF = D_Q_Q2_2D_NInner;
      InnerDOF = D_Q_Q2_2D_Inner;
      break;
    case FE_D_Q_Q3_2D:
      N_DOF = D_Q_Q3_2D_NDOF;
      N_JointDOF = D_Q_Q3_2D_JointDOF;
      JointDOF = D_Q_Q3_2D_J;
      N_InnerDOF = D_Q_Q3_2D_NInner;
      InnerDOF = D_Q_Q3_2D_Inner;
      break;
    case FE_D_Q_Q4_2D:
      N_DOF = D_Q_Q4_2D_NDOF;
      N_JointDOF = D_Q_Q4_2D_JointDOF;
      JointDOF = D_Q_Q4_2D_J;
      N_InnerDOF = D_Q_Q4_2D_NInner;
      InnerDOF = D_Q_Q4_2D_Inner;
      break;
    case FE_C_T_B4_2D:
      N_DOF = C_T_B4_2D_NDOF;
      N_JointDOF = C_T_B4_2D_JointDOF;
      JointDOF = C_T_B4_2D_J;
      N_InnerDOF = C_T_B4_2D_NInner;
      InnerDOF = C_T_B4_2D_Inner;
      N_OuterDOF = C_T_B4_2D_NOuter;
      OuterDOF = C_T_B4_2D_Outer;
      break;
    case FE_D_Q_D2_2D:
      N_DOF = D_Q_D2_2D_NDOF;
      N_JointDOF = D_Q_D2_2D_JointDOF;
      JointDOF = D_Q_D2_2D_J;
      N_InnerDOF = D_Q_D2_2D_NInner;
      InnerDOF = D_Q_D2_2D_Inner;
      break;
    case FE_C_T_UL1_2D:
      N_DOF = C_T_UL1_2D_NDOF;
      N_JointDOF = C_T_UL1_2D_JointDOF;
      JointDOF = C_T_UL1_2D_J;
      N_InnerDOF = C_T_UL1_2D_NInner;
      InnerDOF = C_T_UL1_2D_Inner;
      N_OuterDOF = C_T_UL1_2D_NOuter;
      OuterDOF = C_T_UL1_2D_Outer;
      break;
    case FE_C_T_UL2_2D:
      N_DOF = C_T_UL2_2D_NDOF;
      N_JointDOF = C_T_UL2_2D_JointDOF;
      JointDOF = C_T_UL2_2D_J;
      N_InnerDOF = C_T_UL2_2D_NInner;
      InnerDOF = C_T_UL2_2D_Inner;
      N_OuterDOF = C_T_UL2_2D_NOuter;
      OuterDOF = C_T_UL2_2D_Outer;
      break;
    case FE_C_T_UL3_2D:
      N_DOF = C_T_UL3_2D_NDOF;
      N_JointDOF = C_T_UL3_2D_JointDOF;
      JointDOF = C_T_UL3_2D_J;
      N_InnerDOF = C_T_UL3_2D_NInner;
      InnerDOF = C_T_UL3_2D_Inner;
      N_OuterDOF = C_T_UL3_2D_NOuter;
      OuterDOF = C_T_UL3_2D_Outer;
      break;
    case FE_C_T_UL4_2D:
      N_DOF = C_T_UL4_2D_NDOF;
      N_JointDOF = C_T_UL4_2D_JointDOF;
      JointDOF = C_T_UL4_2D_J;
      N_InnerDOF = C_T_UL4_2D_NInner;
      InnerDOF = C_T_UL4_2D_Inner;
      N_OuterDOF = C_T_UL4_2D_NOuter;
      OuterDOF = C_T_UL4_2D_Outer;
      break;
    case FE_C_T_UL5_2D:
      N_DOF = C_T_UL5_2D_NDOF;
      N_JointDOF = C_T_UL5_2D_JointDOF;
      JointDOF = C_T_UL5_2D_J;
      N_InnerDOF = C_T_UL5_2D_NInner;
      InnerDOF = C_T_UL5_2D_Inner;
      N_OuterDOF = C_T_UL5_2D_NOuter;
      OuterDOF = C_T_UL5_2D_Outer;
      break;
    case FE_C_Q_UL1_2D:
      N_DOF = C_Q_UL1_2D_NDOF;
      N_JointDOF = C_Q_UL1_2D_JointDOF;
      JointDOF = C_Q_UL1_2D_J;
      N_InnerDOF = C_Q_UL1_2D_NInner;
      InnerDOF = C_Q_UL1_2D_Inner;
      N_OuterDOF = C_Q_UL1_2D_NOuter;
      OuterDOF = C_Q_UL1_2D_Outer;
      break;
    case FE_C_Q_UL2_2D:
      N_DOF = C_Q_UL2_2D_NDOF;
      N_JointDOF = C_Q_UL2_2D_JointDOF;
      JointDOF = C_Q_UL2_2D_J;
      N_InnerDOF = C_Q_UL2_2D_NInner;
      InnerDOF = C_Q_UL2_2D_Inner;
      N_OuterDOF = C_Q_UL2_2D_NOuter;
      OuterDOF = C_Q_UL2_2D_Outer;
      break;
    case FE_C_Q_UL3_2D:
      N_DOF = C_Q_UL3_2D_NDOF;
      N_JointDOF = C_Q_UL3_2D_JointDOF;
      JointDOF = C_Q_UL3_2D_J;
      N_InnerDOF = C_Q_UL3_2D_NInner;
      InnerDOF = C_Q_UL3_2D_Inner;
      N_OuterDOF = C_Q_UL3_2D_NOuter;
      OuterDOF = C_Q_UL3_2D_Outer;
      break;
    case FE_C_Q_UL4_2D:
      N_DOF = C_Q_UL4_2D_NDOF;
      N_JointDOF = C_Q_UL4_2D_JointDOF;
      JointDOF = C_Q_UL4_2D_J;
      N_InnerDOF = C_Q_UL4_2D_NInner;
      InnerDOF = C_Q_UL4_2D_Inner;
      N_OuterDOF = C_Q_UL4_2D_NOuter;
      OuterDOF = C_Q_UL4_2D_Outer;
      break;
    case FE_C_Q_UL5_2D:
      N_DOF = C_Q_UL5_2D_NDOF;
      N_JointDOF = C_Q_UL5_2D_JointDOF;
      JointDOF = C_Q_UL5_2D_J;
      N_InnerDOF = C_Q_UL5_2D_NInner;
      InnerDOF = C_Q_UL5_2D_Inner;
      N_OuterDOF = C_Q_UL5_2D_NOuter;
      OuterDOF = C_Q_UL5_2D_Outer;
      break;
    case FE_C_Q_UL2S_2D:
      N_DOF = C_Q_UL2S_2D_NDOF;
      N_JointDOF = C_Q_UL2S_2D_JointDOF;
      JointDOF = C_Q_UL2S_2D_J;
      N_InnerDOF = C_Q_UL2S_2D_NInner;
      InnerDOF = C_Q_UL2S_2D_Inner;
      N_OuterDOF = C_Q_UL2S_2D_NOuter;
      OuterDOF = C_Q_UL2S_2D_Outer;
      break;
    case FE_C_Q_UL3S_2D:
      N_DOF = C_Q_UL3S_2D_NDOF;
      N_JointDOF = C_Q_UL3S_2D_JointDOF;
      JointDOF = C_Q_UL3S_2D_J;
      N_InnerDOF = C_Q_UL3S_2D_NInner;
      InnerDOF = C_Q_UL3S_2D_Inner;
      N_OuterDOF = C_Q_UL3S_2D_NOuter;
      OuterDOF = C_Q_UL3S_2D_Outer;
      break;
    case FE_C_Q_UL4S_2D:
      N_DOF = C_Q_UL4S_2D_NDOF;
      N_JointDOF = C_Q_UL4S_2D_JointDOF;
      JointDOF = C_Q_UL4S_2D_J;
      N_InnerDOF = C_Q_UL4S_2D_NInner;
      InnerDOF = C_Q_UL4S_2D_Inner;
      N_OuterDOF = C_Q_UL4S_2D_NOuter;
      OuterDOF = C_Q_UL4S_2D_Outer;
      break;
    case FE_C_Q_UL5S_2D:
      N_DOF = C_Q_UL5S_2D_NDOF;
      N_JointDOF = C_Q_UL5S_2D_JointDOF;
      JointDOF = C_Q_UL5S_2D_J;
      N_InnerDOF = C_Q_UL5S_2D_NInner;
      InnerDOF = C_Q_UL5S_2D_Inner;
      N_OuterDOF = C_Q_UL5S_2D_NOuter;
      OuterDOF = C_Q_UL5S_2D_Outer;
      break;
    case FE_C_Q_UL6S_2D:
      N_DOF = C_Q_UL6S_2D_NDOF;
      N_JointDOF = C_Q_UL6S_2D_JointDOF;
      JointDOF = C_Q_UL6S_2D_J;
      N_InnerDOF = C_Q_UL6S_2D_NInner;
      InnerDOF = C_Q_UL6S_2D_Inner;
      N_OuterDOF = C_Q_UL6S_2D_NOuter;
      OuterDOF = C_Q_UL6S_2D_Outer;
      break;
    case FE_C_Q_UL7S_2D:
      N_DOF = C_Q_UL7S_2D_NDOF;
      N_JointDOF = C_Q_UL7S_2D_JointDOF;
      JointDOF = C_Q_UL7S_2D_J;
      N_InnerDOF = C_Q_UL7S_2D_NInner;
      InnerDOF = C_Q_UL7S_2D_Inner;
      N_OuterDOF = C_Q_UL7S_2D_NOuter;
      OuterDOF = C_Q_UL7S_2D_Outer;
      break;
    case FE_C_Q_UL8S_2D:
      N_DOF = C_Q_UL8S_2D_NDOF;
      N_JointDOF = C_Q_UL8S_2D_JointDOF;
      JointDOF = C_Q_UL8S_2D_J;
      N_InnerDOF = C_Q_UL8S_2D_NInner;
      InnerDOF = C_Q_UL8S_2D_Inner;
      N_OuterDOF = C_Q_UL8S_2D_NOuter;
      OuterDOF = C_Q_UL8S_2D_Outer;
      break;
    case FE_C_Q_UL9S_2D:
      N_DOF = C_Q_UL9S_2D_NDOF;
      N_JointDOF = C_Q_UL9S_2D_JointDOF;
      JointDOF = C_Q_UL9S_2D_J;
      N_InnerDOF = C_Q_UL9S_2D_NInner;
      InnerDOF = C_Q_UL9S_2D_Inner;
      N_OuterDOF = C_Q_UL9S_2D_NOuter;
      OuterDOF = C_Q_UL9S_2D_Outer;
      break;
    case FE_C_Q_UL2SE_2D:
      N_DOF = C_Q_UL2SE_2D_NDOF;
      N_JointDOF = C_Q_UL2SE_2D_JointDOF;
      JointDOF = C_Q_UL2SE_2D_J;
      N_InnerDOF = C_Q_UL2SE_2D_NInner;
      InnerDOF = C_Q_UL2SE_2D_Inner;
      N_OuterDOF = C_Q_UL2SE_2D_NOuter;
      OuterDOF = C_Q_UL2SE_2D_Outer;
      break;
    case FE_C_Q_UL3SE_2D:
      N_DOF = C_Q_UL3SE_2D_NDOF;
      N_JointDOF = C_Q_UL3SE_2D_JointDOF;
      JointDOF = C_Q_UL3SE_2D_J;
      N_InnerDOF = C_Q_UL3SE_2D_NInner;
      InnerDOF = C_Q_UL3SE_2D_Inner;
      N_OuterDOF = C_Q_UL3SE_2D_NOuter;
      OuterDOF = C_Q_UL3SE_2D_Outer;
      break;
    case FE_C_Q_UL4SE_2D:
      N_DOF = C_Q_UL4SE_2D_NDOF;
      N_JointDOF = C_Q_UL4SE_2D_JointDOF;
      JointDOF = C_Q_UL4SE_2D_J;
      N_InnerDOF = C_Q_UL4SE_2D_NInner;
      InnerDOF = C_Q_UL4SE_2D_Inner;
      N_OuterDOF = C_Q_UL4SE_2D_NOuter;
      OuterDOF = C_Q_UL4SE_2D_Outer;
      break;
    case FE_C_Q_UL5SE_2D:
      N_DOF = C_Q_UL5SE_2D_NDOF;
      N_JointDOF = C_Q_UL5SE_2D_JointDOF;
      JointDOF = C_Q_UL5SE_2D_J;
      N_InnerDOF = C_Q_UL5SE_2D_NInner;
      InnerDOF = C_Q_UL5SE_2D_Inner;
      N_OuterDOF = C_Q_UL5SE_2D_NOuter;
      OuterDOF = C_Q_UL5SE_2D_Outer;
      break;
    case FE_C_Q_UL6SE_2D:
      N_DOF = C_Q_UL6SE_2D_NDOF;
      N_JointDOF = C_Q_UL6SE_2D_JointDOF;
      JointDOF = C_Q_UL6SE_2D_J;
      N_InnerDOF = C_Q_UL6SE_2D_NInner;
      InnerDOF = C_Q_UL6SE_2D_Inner;
      N_OuterDOF = C_Q_UL6SE_2D_NOuter;
      OuterDOF = C_Q_UL6SE_2D_Outer;
      break;
    case FE_C_Q_UL7SE_2D:
      N_DOF = C_Q_UL7SE_2D_NDOF;
      N_JointDOF = C_Q_UL7SE_2D_JointDOF;
      JointDOF = C_Q_UL7SE_2D_J;
      N_InnerDOF = C_Q_UL7SE_2D_NInner;
      InnerDOF = C_Q_UL7SE_2D_Inner;
      N_OuterDOF = C_Q_UL7SE_2D_NOuter;
      OuterDOF = C_Q_UL7SE_2D_Outer;
      break;
    case FE_C_Q_UL8SE_2D:
      N_DOF = C_Q_UL8SE_2D_NDOF;
      N_JointDOF = C_Q_UL8SE_2D_JointDOF;
      JointDOF = C_Q_UL8SE_2D_J;
      N_InnerDOF = C_Q_UL8SE_2D_NInner;
      InnerDOF = C_Q_UL8SE_2D_Inner;
      N_OuterDOF = C_Q_UL8SE_2D_NOuter;
      OuterDOF = C_Q_UL8SE_2D_Outer;
      break;
    case FE_C_Q_UL9SE_2D:
      N_DOF = C_Q_UL9SE_2D_NDOF;
      N_JointDOF = C_Q_UL9SE_2D_JointDOF;
      JointDOF = C_Q_UL9SE_2D_J;
      N_InnerDOF = C_Q_UL9SE_2D_NInner;
      InnerDOF = C_Q_UL9SE_2D_Inner;
      N_OuterDOF = C_Q_UL9SE_2D_NOuter;
      OuterDOF = C_Q_UL9SE_2D_Outer;
      break;
    case FE_C_Q_M2_2D:
      N_DOF = C_Q_M2_2D_NDOF;
      N_JointDOF = C_Q_M2_2D_JointDOF;
      JointDOF = C_Q_M2_2D_J;
      N_InnerDOF = C_Q_M2_2D_NInner;
      InnerDOF = C_Q_M2_2D_Inner;
      N_OuterDOF = C_Q_M2_2D_NOuter;
      OuterDOF = C_Q_M2_2D_Outer;
      break;
    case FE_C_Q_M3_2D:
      N_DOF = C_Q_M3_2D_NDOF;
      N_JointDOF = C_Q_M3_2D_JointDOF;
      JointDOF = C_Q_M3_2D_J;
      N_InnerDOF = C_Q_M3_2D_NInner;
      InnerDOF = C_Q_M3_2D_Inner;
      N_OuterDOF = C_Q_M3_2D_NOuter;
      OuterDOF = C_Q_M3_2D_Outer;
      break;
    case FE_C_Q_M4_2D:
      N_DOF = C_Q_M4_2D_NDOF;
      N_JointDOF = C_Q_M4_2D_JointDOF;
      JointDOF = C_Q_M4_2D_J;
      N_InnerDOF = C_Q_M4_2D_NInner;
      InnerDOF = C_Q_M4_2D_Inner;
      N_OuterDOF = C_Q_M4_2D_NOuter;
      OuterDOF = C_Q_M4_2D_Outer;
      break;
    case FE_C_Q_M5_2D:
      N_DOF = C_Q_M5_2D_NDOF;
      N_JointDOF = C_Q_M5_2D_JointDOF;
      JointDOF = C_Q_M5_2D_J;
      N_InnerDOF = C_Q_M5_2D_NInner;
      InnerDOF = C_Q_M5_2D_Inner;
      N_OuterDOF = C_Q_M5_2D_NOuter;
      OuterDOF = C_Q_M5_2D_Outer;
      break;
    case FE_C_Q_M6_2D:
      N_DOF = C_Q_M6_2D_NDOF;
      N_JointDOF = C_Q_M6_2D_JointDOF;
      JointDOF = C_Q_M6_2D_J;
      N_InnerDOF = C_Q_M6_2D_NInner;
      InnerDOF = C_Q_M6_2D_Inner;
      N_OuterDOF = C_Q_M6_2D_NOuter;
      OuterDOF = C_Q_M6_2D_Outer;
      break;
    case FE_C_Q_M7_2D:
      N_DOF = C_Q_M7_2D_NDOF;
      N_JointDOF = C_Q_M7_2D_JointDOF;
      JointDOF = C_Q_M7_2D_J;
      N_InnerDOF = C_Q_M7_2D_NInner;
      InnerDOF = C_Q_M7_2D_Inner;
      N_OuterDOF = C_Q_M7_2D_NOuter;
      OuterDOF = C_Q_M7_2D_Outer;
      break;
    case FE_C_Q_M8_2D:
      N_DOF = C_Q_M8_2D_NDOF;
      N_JointDOF = C_Q_M8_2D_JointDOF;
      JointDOF = C_Q_M8_2D_J;
      N_InnerDOF = C_Q_M8_2D_NInner;
      InnerDOF = C_Q_M8_2D_Inner;
      N_OuterDOF = C_Q_M8_2D_NOuter;
      OuterDOF = C_Q_M8_2D_Outer;
      break;
    case FE_C_Q_M9_2D:
      N_DOF = C_Q_M9_2D_NDOF;
      N_JointDOF = C_Q_M9_2D_JointDOF;
      JointDOF = C_Q_M9_2D_J;
      N_InnerDOF = C_Q_M9_2D_NInner;
      InnerDOF = C_Q_M9_2D_Inner;
      N_OuterDOF = C_Q_M9_2D_NOuter;
      OuterDOF = C_Q_M9_2D_Outer;
      break;
    case FE_C_Q_EL1_2D:
      N_DOF = C_Q_EL1_2D_NDOF;
      N_JointDOF = C_Q_EL1_2D_JointDOF;
      JointDOF = C_Q_EL1_2D_J;
      N_InnerDOF = C_Q_EL1_2D_NInner;
      InnerDOF = C_Q_EL1_2D_Inner;
      N_OuterDOF = C_Q_EL1_2D_NOuter;
      OuterDOF = C_Q_EL1_2D_Outer;
      break;
    case FE_N_Q_RT0_2D:
      N_DOF = N_Q_RT0_2D_NDOF;
      N_JointDOF = N_Q_RT0_2D_JointDOF;
      JointDOF = N_Q_RT0_2D_J;
      N_InnerDOF = N_Q_RT0_2D_NInner;
      InnerDOF = N_Q_RT0_2D_Inner;
      N_OuterDOF = N_Q_RT0_2D_NOuter;
      OuterDOF = N_Q_RT0_2D_Outer;
      break;
    case FE_N_Q_RT1_2D:
      N_DOF = N_Q_RT1_2D_NDOF;
      N_JointDOF = N_Q_RT1_2D_JointDOF;
      JointDOF = N_Q_RT1_2D_J;
      N_InnerDOF = N_Q_RT1_2D_NInner;
      InnerDOF = N_Q_RT1_2D_Inner;
      N_OuterDOF = N_Q_RT1_2D_NOuter;
      OuterDOF = N_Q_RT1_2D_Outer;
      break;
    case FE_N_Q_RT2_2D:
      N_DOF = N_Q_RT2_2D_NDOF;
      N_JointDOF = N_Q_RT2_2D_JointDOF;
      JointDOF = N_Q_RT2_2D_J;
      N_InnerDOF = N_Q_RT2_2D_NInner;
      InnerDOF = N_Q_RT2_2D_Inner;
      N_OuterDOF = N_Q_RT2_2D_NOuter;
      OuterDOF = N_Q_RT2_2D_Outer;
      break;
    case FE_N_Q_RT3_2D:
      N_DOF = N_Q_RT3_2D_NDOF;
      N_JointDOF = N_Q_RT3_2D_JointDOF;
      JointDOF = N_Q_RT3_2D_J;
      N_InnerDOF = N_Q_RT3_2D_NInner;
      InnerDOF = N_Q_RT3_2D_Inner;
      N_OuterDOF = N_Q_RT3_2D_NOuter;
      OuterDOF = N_Q_RT3_2D_Outer;
      break;
    case FE_N_T_RT0_2D:
      N_DOF = N_T_RT0_2D_NDOF;
      N_JointDOF = N_T_RT0_2D_JointDOF;
      JointDOF = N_T_RT0_2D_J;
      N_InnerDOF = N_T_RT0_2D_NInner;
      InnerDOF = N_T_RT0_2D_Inner;
      N_OuterDOF = N_T_RT0_2D_NOuter;
      OuterDOF = N_T_RT0_2D_Outer;
      break;
    case FE_N_T_RT1_2D:
      N_DOF = N_T_RT1_2D_NDOF;
      N_JointDOF = N_T_RT1_2D_JointDOF;
      JointDOF = N_T_RT1_2D_J;
      N_InnerDOF = N_T_RT1_2D_NInner;
      InnerDOF = N_T_RT1_2D_Inner;
      N_OuterDOF = N_T_RT1_2D_NOuter;
      OuterDOF = N_T_RT1_2D_Outer;
      break;
    case FE_N_T_RT2_2D:
      N_DOF = N_T_RT2_2D_NDOF;
      N_JointDOF = N_T_RT2_2D_JointDOF;
      JointDOF = N_T_RT2_2D_J;
      N_InnerDOF = N_T_RT2_2D_NInner;
      InnerDOF = N_T_RT2_2D_Inner;
      N_OuterDOF = N_T_RT2_2D_NOuter;
      OuterDOF = N_T_RT2_2D_Outer;
      break;
    case FE_N_T_RT3_2D:
      N_DOF = N_T_RT3_2D_NDOF;
      N_JointDOF = N_T_RT3_2D_JointDOF;
      JointDOF = N_T_RT3_2D_J;
      N_InnerDOF = N_T_RT3_2D_NInner;
      InnerDOF = N_T_RT3_2D_Inner;
      N_OuterDOF = N_T_RT3_2D_NOuter;
      OuterDOF = N_T_RT3_2D_Outer;
      break;
    case FE_N_Q_BDM1_2D:
      N_DOF = N_Q_BDM1_2D_NDOF;
      N_JointDOF = N_Q_BDM1_2D_JointDOF;
      JointDOF = N_Q_BDM1_2D_J;
      N_InnerDOF = N_Q_BDM1_2D_NInner;
      InnerDOF = N_Q_BDM1_2D_Inner;
      N_OuterDOF = N_Q_BDM1_2D_NOuter;
      OuterDOF = N_Q_BDM1_2D_Outer;
      break;
    case FE_N_Q_BDM2_2D:
      N_DOF = N_Q_BDM2_2D_NDOF;
      N_JointDOF = N_Q_BDM2_2D_JointDOF;
      JointDOF = N_Q_BDM2_2D_J;
      N_InnerDOF = N_Q_BDM2_2D_NInner;
      InnerDOF = N_Q_BDM2_2D_Inner;
      N_OuterDOF = N_Q_BDM2_2D_NOuter;
      OuterDOF = N_Q_BDM2_2D_Outer;
      break;
    case FE_N_Q_BDM3_2D:
      N_DOF = N_Q_BDM3_2D_NDOF;
      N_JointDOF = N_Q_BDM3_2D_JointDOF;
      JointDOF = N_Q_BDM3_2D_J;
      N_InnerDOF = N_Q_BDM3_2D_NInner;
      InnerDOF = N_Q_BDM3_2D_Inner;
      N_OuterDOF = N_Q_BDM3_2D_NOuter;
      OuterDOF = N_Q_BDM3_2D_Outer;
      break;
    case FE_N_T_BDM1_2D:
      N_DOF = N_T_BDM1_2D_NDOF;
      N_JointDOF = N_T_BDM1_2D_JointDOF;
      JointDOF = N_T_BDM1_2D_J;
      N_InnerDOF = N_T_BDM1_2D_NInner;
      InnerDOF = N_T_BDM1_2D_Inner;
      N_OuterDOF = N_T_BDM1_2D_NOuter;
      OuterDOF = N_T_BDM1_2D_Outer;
      break;
    case FE_N_T_BDM2_2D:
      N_DOF = N_T_BDM2_2D_NDOF;
      N_JointDOF = N_T_BDM2_2D_JointDOF;
      JointDOF = N_T_BDM2_2D_J;
      N_InnerDOF = N_T_BDM2_2D_NInner;
      InnerDOF = N_T_BDM2_2D_Inner;
      N_OuterDOF = N_T_BDM2_2D_NOuter;
      OuterDOF = N_T_BDM2_2D_Outer;
      break;
    case FE_N_T_BDM3_2D:
      N_DOF = N_T_BDM3_2D_NDOF;
      N_JointDOF = N_T_BDM3_2D_JointDOF;
      JointDOF = N_T_BDM3_2D_J;
      N_InnerDOF = N_T_BDM3_2D_NInner;
      InnerDOF = N_T_BDM3_2D_Inner;
      N_OuterDOF = N_T_BDM3_2D_NOuter;
      OuterDOF = N_T_BDM3_2D_Outer;
      break;
    case FE_C_T_P00_3D:
      N_DOF = C_T_P00_3D_NDOF;
      N_JointDOF = C_T_P00_3D_JointDOF;
      JointDOF = C_T_P00_3D_J;
      N_InnerDOF = C_T_P00_3D_NInner;
      InnerDOF = C_T_P00_3D_Inner;
      break;
    case FE_C_T_P0_3D:
      N_DOF = C_T_P0_3D_NDOF;
      N_JointDOF = C_T_P0_3D_JointDOF;
      JointDOF = C_T_P0_3D_J;
      N_InnerDOF = C_T_P0_3D_NInner;
      InnerDOF = C_T_P0_3D_Inner;
      N_EdgeDOF = C_T_P0_3D_EdgeDOF;
      EdgeDOF = C_T_P0_3D_E;
      N_VertDOF = C_T_P0_3D_VertDOF;
      VertDOF = C_T_P0_3D_Vert; 
      EdgeVertData_Filled = true;
      break;
    case FE_C_T_P1_3D:
      N_DOF = C_T_P1_3D_NDOF;
      N_JointDOF = C_T_P1_3D_JointDOF;
      JointDOF = C_T_P1_3D_J;
      N_InnerDOF = C_T_P1_3D_NInner;
      InnerDOF = C_T_P1_3D_Inner;
      N_EdgeDOF = C_T_P1_3D_EdgeDOF;
      EdgeDOF = C_T_P1_3D_E;
      N_VertDOF = C_T_P1_3D_VertDOF;
      VertDOF = C_T_P1_3D_Vert; 
      EdgeVertData_Filled = true;
      break;
    case FE_C_T_P2_3D:
      N_DOF = C_T_P2_3D_NDOF;
      N_JointDOF = C_T_P2_3D_JointDOF;
      JointDOF = C_T_P2_3D_J;
      N_InnerDOF = C_T_P2_3D_NInner;
      InnerDOF = C_T_P2_3D_Inner;
      N_EdgeDOF = C_T_P2_3D_EdgeDOF;
      EdgeDOF = C_T_P2_3D_E;
      N_VertDOF = C_T_P2_3D_VertDOF;
      VertDOF = C_T_P2_3D_Vert; 
      EdgeVertData_Filled = true;
      break;
    case FE_C_T_P3_3D:
      N_DOF = C_T_P3_3D_NDOF;
      N_JointDOF = C_T_P3_3D_JointDOF;
      JointDOF = C_T_P3_3D_J;
      N_InnerDOF = C_T_P3_3D_NInner;
      InnerDOF = C_T_P3_3D_Inner;
      break;
    case FE_N_T_P1_3D:
      N_DOF = N_T_P1_3D_NDOF;
      N_JointDOF = N_T_P1_3D_JointDOF;
      JointDOF = N_T_P1_3D_J;
      N_InnerDOF = N_T_P1_3D_NInner;
      InnerDOF = N_T_P1_3D_Inner;
      N_EdgeDOF = N_T_P1_3D_EdgeDOF;
      EdgeDOF = N_T_P1_3D_E;
      N_VertDOF = N_T_P1_3D_VertDOF;
      VertDOF = N_T_P1_3D_Vert;
      EdgeVertData_Filled = true;
      break;
    case FE_C_H_Q00_3D:
      N_DOF = C_H_Q00_3D_NDOF;
      N_JointDOF = C_H_Q00_3D_JointDOF;
      JointDOF = C_H_Q00_3D_J;
      N_InnerDOF = C_H_Q00_3D_NInner;
      InnerDOF = C_H_Q00_3D_Inner;
      break;
    case FE_C_H_Q0_3D:
      N_DOF = C_H_Q0_3D_NDOF;
      N_JointDOF = C_H_Q0_3D_JointDOF;
      JointDOF = C_H_Q0_3D_J;
      N_InnerDOF = C_H_Q0_3D_NInner;
      InnerDOF = C_H_Q0_3D_Inner;
      N_EdgeDOF = C_H_Q0_3D_EdgeDOF;
      EdgeDOF = C_H_Q0_3D_E;
      N_VertDOF = C_H_Q0_3D_VertDOF;
      VertDOF = C_H_Q0_3D_Vert; 
      EdgeVertData_Filled = true;
      break;
    case FE_C_H_Q1_3D:
      N_DOF = C_H_Q1_3D_NDOF;
      N_JointDOF = C_H_Q1_3D_JointDOF;
      JointDOF = C_H_Q1_3D_J;
      N_InnerDOF = C_H_Q1_3D_NInner;
      InnerDOF = C_H_Q1_3D_Inner;
      N_EdgeDOF = C_H_Q1_3D_EdgeDOF;
      EdgeDOF = C_H_Q1_3D_E;
      N_VertDOF = C_H_Q1_3D_VertDOF;
      VertDOF = C_H_Q1_3D_Vert;
      EdgeVertData_Filled = true;
      break;
    case FE_C_H_Q2_3D:
      N_DOF = C_H_Q2_3D_NDOF;
      N_JointDOF = C_H_Q2_3D_JointDOF;
      JointDOF = C_H_Q2_3D_J;
      N_InnerDOF = C_H_Q2_3D_NInner;
      InnerDOF = C_H_Q2_3D_Inner;
      N_OuterDOF = C_H_Q2_3D_NOuter;
      OuterDOF = C_H_Q2_3D_Outer;
      N_EdgeDOF = C_H_Q2_3D_EdgeDOF;
      EdgeDOF = C_H_Q2_3D_E;
      N_VertDOF = C_H_Q2_3D_VertDOF;
      VertDOF = C_H_Q2_3D_Vert;
      EdgeVertData_Filled = true;
      break;
    case FE_C_H_Q3_3D:
      N_DOF = C_H_Q3_3D_NDOF;
      N_JointDOF = C_H_Q3_3D_JointDOF;
      JointDOF = C_H_Q3_3D_J;
      N_InnerDOF = C_H_Q3_3D_NInner;
      InnerDOF = C_H_Q3_3D_Inner;
      N_OuterDOF = C_H_Q3_3D_NOuter;
      OuterDOF = C_H_Q3_3D_Outer;
      break;
    case FE_C_H_Q4_3D:
      N_DOF = C_H_Q4_3D_NDOF;
      N_JointDOF = C_H_Q4_3D_JointDOF;
      JointDOF = C_H_Q4_3D_J;
      N_InnerDOF = C_H_Q4_3D_NInner;
      InnerDOF = C_H_Q4_3D_Inner;
      break;
    case FE_N_H_Q1_3D:
      N_DOF = N_H_Q1_3D_NDOF;
      N_JointDOF = N_H_Q1_3D_JointDOF;
      JointDOF = N_H_Q1_3D_J;
      N_InnerDOF = N_H_Q1_3D_NInner;
      InnerDOF = N_H_Q1_3D_Inner;
      N_EdgeDOF = C_N_Q1_3D_EdgeDOF;
      EdgeDOF = C_N_Q1_3D_E;
      N_VertDOF = C_N_Q1_3D_VertDOF;
      VertDOF = C_N_Q1_3D_Vert;
      EdgeVertData_Filled = true;
      break;
    case FE_D_H_P1_3D:
      N_DOF = D_H_P1_3D_NDOF;
      N_JointDOF = D_H_P1_3D_JointDOF;
      JointDOF = D_H_P1_3D_J;
      N_InnerDOF = D_H_P1_3D_NInner;
      InnerDOF = D_H_P1_3D_Inner;
      N_EdgeDOF = D_H_P1_3D_EdgeDOF;
      EdgeDOF = D_H_P1_3D_E;
      N_VertDOF = D_H_P1_3D_VertDOF;
      VertDOF = D_H_P1_3D_Vert;
      EdgeVertData_Filled = true;
      break;
    case FE_D_H_P2_3D:
      N_DOF = D_H_P2_3D_NDOF;
      N_JointDOF = D_H_P2_3D_JointDOF;
      JointDOF = D_H_P2_3D_J;
      N_InnerDOF = D_H_P2_3D_NInner;
      InnerDOF = D_H_P2_3D_Inner;
      break;
    case FE_D_H_P3_3D:
      N_DOF = D_H_P3_3D_NDOF;
      N_JointDOF = D_H_P3_3D_JointDOF;
      JointDOF = D_H_P3_3D_J;
      N_InnerDOF = D_H_P3_3D_NInner;
      InnerDOF = D_H_P3_3D_Inner;
      break;
    case FE_D_H_Q1_3D:
      N_DOF = D_H_Q1_3D_NDOF;
      N_JointDOF = D_H_Q1_3D_JointDOF;
      JointDOF = D_H_Q1_3D_J;
      N_InnerDOF = D_H_Q1_3D_NInner;
      InnerDOF = D_H_Q1_3D_Inner;
      N_EdgeDOF = D_H_Q1_3D_EdgeDOF;
      EdgeDOF = D_H_Q1_3D_E;
      N_VertDOF = D_H_Q1_3D_VertDOF;
      VertDOF = D_H_Q1_3D_Vert;
      EdgeVertData_Filled = true;
      break;
    case FE_D_H_Q2_3D:
      N_DOF = D_H_Q2_3D_NDOF;
      N_JointDOF = D_H_Q2_3D_JointDOF;
      JointDOF = D_H_Q2_3D_J;
      N_InnerDOF = D_H_Q2_3D_NInner;
      InnerDOF = D_H_Q2_3D_Inner;
      N_EdgeDOF = D_H_Q2_3D_EdgeDOF;
      EdgeDOF = D_H_Q2_3D_E;
      N_VertDOF = D_H_Q2_3D_VertDOF;
      VertDOF = D_H_Q2_3D_Vert;
      EdgeVertData_Filled = true;
      break;
    case FE_B_H_IB2_3D:
      N_DOF = B_H_IB2_3D_NDOF;
      N_JointDOF = B_H_IB2_3D_JointDOF;
      JointDOF = B_H_IB2_3D_J;
      N_InnerDOF = B_H_IB2_3D_NInner;
      InnerDOF = B_H_IB2_3D_Inner;
      break;
    case FE_N_T_P2_3D:
      N_DOF = N_T_P2_3D_NDOF;
      N_JointDOF = N_T_P2_3D_JointDOF;
      JointDOF = N_T_P2_3D_J;
      N_InnerDOF = N_T_P2_3D_NInner;
      InnerDOF = N_T_P2_3D_Inner;
      break;
    case FE_N_T_P3_3D:
      ErrThrow("FE descriptor of type FE_N_T_P3_3D not implemented");
      break;
    case FE_N_T_P4_3D:
      ErrThrow("FE descriptor of type FE_N_T_P4_3D not implemented");
      break;
    case FE_N_T_P5_3D:
      ErrThrow("FE descriptor of type FE_N_T_P5_3D not implemented");
      break;
    case FE_N_H_Q2_3D:
      N_DOF = N_H_Q2_3D_NDOF;
      N_JointDOF = N_H_Q2_3D_JointDOF;
      JointDOF = N_H_Q2_3D_J;
      N_InnerDOF = N_H_Q2_3D_NInner;
      InnerDOF = N_H_Q2_3D_Inner;
      N_OuterDOF = N_H_Q2_3D_NOuter;
      OuterDOF = N_H_Q2_3D_Outer;
      break;
    case FE_N_H_Q3_3D:
      N_DOF = N_H_Q3_3D_NDOF;
      N_JointDOF = N_H_Q3_3D_JointDOF;
      JointDOF = N_H_Q3_3D_J;
      N_InnerDOF = N_H_Q3_3D_NInner;
      InnerDOF = N_H_Q3_3D_Inner;
      break;
    case FE_N_H_Q4_3D:
      N_DOF = N_H_Q4_3D_NDOF;
      N_JointDOF = N_H_Q4_3D_JointDOF;
      JointDOF = N_H_Q4_3D_J;
      N_InnerDOF = N_H_Q4_3D_NInner;
      InnerDOF = N_H_Q4_3D_Inner;
      break;
    case FE_N_H_Q5_3D:
      ErrThrow("FE descriptor of type FE_N_H_Q5_3D not implemented");
      break;
    case FE_C_T_B2_3D:
      N_DOF = C_T_B2_3D_NDOF;
      N_JointDOF = C_T_B2_3D_JointDOF;
      JointDOF = C_T_B2_3D_J;
      N_InnerDOF = C_T_B2_3D_NInner;
      InnerDOF = C_T_B2_3D_Inner;
      N_EdgeDOF = C_T_B2_3D_EdgeDOF;
      EdgeDOF = C_T_B2_3D_E;
      N_VertDOF = C_T_B2_3D_VertDOF;
      VertDOF = C_T_B2_3D_Vert; 
      EdgeVertData_Filled = true;
      break;
    case FE_D_T_P1_3D:
      N_DOF = D_T_P1_3D_NDOF;
      N_JointDOF = D_T_P1_3D_JointDOF;
      JointDOF = D_T_P1_3D_J;
      N_InnerDOF = D_T_P1_3D_NInner;
      InnerDOF = D_T_P1_3D_Inner;
      N_EdgeDOF = D_T_P1_3D_EdgeDOF;
      EdgeDOF = D_T_P1_3D_E;
      N_VertDOF = D_T_P1_3D_VertDOF;
      VertDOF = D_T_P1_3D_Vert;
      EdgeVertData_Filled = true;
      break;
    case FE_D_T_P2_3D:
      N_DOF = D_T_P2_3D_NDOF;
      N_JointDOF = D_T_P2_3D_JointDOF;
      JointDOF = D_T_P2_3D_J;
      N_InnerDOF = D_T_P2_3D_NInner;
      InnerDOF = D_T_P2_3D_Inner;
      N_EdgeDOF = D_T_P2_3D_EdgeDOF;
      EdgeDOF = D_T_P2_3D_E;
      N_VertDOF = D_T_P2_3D_VertDOF;
      VertDOF = D_T_P2_3D_Vert;
      EdgeVertData_Filled = true;
      break;
    case FE_D_T_P3_3D:
      N_DOF = D_T_P3_3D_NDOF;
      N_JointDOF = D_T_P3_3D_JointDOF;
      JointDOF = D_T_P3_3D_J;
      N_InnerDOF = D_T_P3_3D_NInner;
      InnerDOF = D_T_P3_3D_Inner;
      N_EdgeDOF = D_T_P3_3D_EdgeDOF;
      EdgeDOF = D_T_P3_3D_E;
      N_VertDOF = D_T_P3_3D_VertDOF;
      VertDOF = D_T_P3_3D_Vert;
      EdgeVertData_Filled = true;
      break;
    case FE_C_H_UL1_3D:
      N_DOF = C_H_UL1_3D_NDOF;
      N_JointDOF = C_H_UL1_3D_JointDOF;
      JointDOF = C_H_UL1_3D_J;
      N_InnerDOF = C_H_UL1_3D_NInner;
      InnerDOF = C_H_UL1_3D_Inner;
      N_OuterDOF = C_H_UL1_3D_NOuter;
      OuterDOF = C_H_UL1_3D_Outer;
      break;
    case FE_C_H_UL2_3D:
      N_DOF = C_H_UL2_3D_NDOF;
      N_JointDOF = C_H_UL2_3D_JointDOF;
      JointDOF = C_H_UL2_3D_J;
      N_InnerDOF = C_H_UL2_3D_NInner;
      InnerDOF = C_H_UL2_3D_Inner;
      N_OuterDOF = C_H_UL2_3D_NOuter;
      OuterDOF = C_H_UL2_3D_Outer;
      break;
    case FE_C_H_UL3_3D:
      N_DOF = C_H_UL3_3D_NDOF;
      N_JointDOF = C_H_UL3_3D_JointDOF;
      JointDOF = C_H_UL3_3D_J;
      N_InnerDOF = C_H_UL3_3D_NInner;
      InnerDOF = C_H_UL3_3D_Inner;
      N_OuterDOF = C_H_UL3_3D_NOuter;
      OuterDOF = C_H_UL3_3D_Outer;
      break;
    case FE_C_H_UL4_3D:
      ErrThrow("FE descriptor of type FE_C_H_UL4_3D not implemented");
      break;
    case FE_C_H_UL5_3D:
      ErrThrow("FE descriptor of type FE_C_H_UL5_3D not implemented");
      break;
    case FE_C_T_UL1_3D:
      ErrThrow("FE descriptor of type FE_C_T_UL1_3D not implemented");
      break;
    case FE_C_T_UL2_3D:
      ErrThrow("FE descriptor of type FE_C_T_UL2_3D not implemented");
      break;
    case FE_C_T_UL3_3D:
      ErrThrow("FE descriptor of type FE_C_T_UL3_3D not implemented");
      break;
    case FE_C_T_UL4_3D:
      ErrThrow("FE descriptor of type FE_C_T_UL4_3D not implemented");
      break;
    case FE_C_T_UL5_3D:
      ErrThrow("FE descriptor of type FE_C_T_UL5_3D not implemented");
      break;
    case FE_N_T_RT0_3D:
      N_DOF = N_T_RT0_3D_NDOF;
      N_JointDOF = N_T_RT0_3D_JointDOF;
      JointDOF = N_T_RT0_3D_J;
      N_InnerDOF = N_T_RT0_3D_NInner;
      InnerDOF = N_T_RT0_3D_Inner;
      N_OuterDOF = N_T_RT0_3D_NOuter;
      OuterDOF = N_T_RT0_3D_Outer;
      break;
    case FE_N_T_RT1_3D:
      N_DOF = N_T_RT1_3D_NDOF;
      N_JointDOF = N_T_RT1_3D_JointDOF;
      JointDOF = N_T_RT1_3D_J;
      N_InnerDOF = N_T_RT1_3D_NInner;
      InnerDOF = N_T_RT1_3D_Inner;
      N_OuterDOF = N_T_RT1_3D_NOuter;
      OuterDOF = N_T_RT1_3D_Outer;
      break;
    case FE_N_T_RT2_3D:
      N_DOF = N_T_RT2_3D_NDOF;
      N_JointDOF = N_T_RT2_3D_JointDOF;
      JointDOF = N_T_RT2_3D_J;
      N_InnerDOF = N_T_RT2_3D_NInner;
      InnerDOF = N_T_RT2_3D_Inner;
      N_OuterDOF = N_T_RT2_3D_NOuter;
      OuterDOF = N_T_RT2_3D_Outer;
      break;
    case FE_N_T_RT3_3D:
      N_DOF = N_T_RT3_3D_NDOF;
      N_JointDOF = N_T_RT3_3D_JointDOF;
      JointDOF = N_T_RT3_3D_J;
      N_InnerDOF = N_T_RT3_3D_NInner;
      InnerDOF = N_T_RT3_3D_Inner;
      N_OuterDOF = N_T_RT3_3D_NOuter;
      OuterDOF = N_T_RT3_3D_Outer;
      break;
    case FE_N_H_RT0_3D:
      N_DOF = N_H_RT0_3D_NDOF;
      N_JointDOF = N_H_RT0_3D_JointDOF;
      JointDOF = N_H_RT0_3D_J;
      N_InnerDOF = N_H_RT0_3D_NInner;
      InnerDOF = N_H_RT0_3D_Inner;
      N_OuterDOF = N_H_RT0_3D_NOuter;
      OuterDOF = N_H_RT0_3D_Outer;
      break;
    case FE_N_H_RT1_3D:
      N_DOF = N_H_RT1_3D_NDOF;
      N_JointDOF = N_H_RT1_3D_JointDOF;
      JointDOF = N_H_RT1_3D_J;
      N_InnerDOF = N_H_RT1_3D_NInner;
      InnerDOF = N_H_RT1_3D_Inner;
      N_OuterDOF = N_H_RT1_3D_NOuter;
      OuterDOF = N_H_RT1_3D_Outer;
      break;
    case FE_N_H_RT2_3D:
      N_DOF = N_H_RT2_3D_NDOF;
      N_JointDOF = N_H_RT2_3D_JointDOF;
      JointDOF = N_H_RT2_3D_J;
      N_InnerDOF = N_H_RT2_3D_NInner;
      InnerDOF = N_H_RT2_3D_Inner;
      N_OuterDOF = N_H_RT2_3D_NOuter;
      OuterDOF = N_H_RT2_3D_Outer;
      break;
    case FE_N_T_BDDF1_3D:
      N_DOF = N_T_BDDF1_3D_NDOF;
      N_JointDOF = N_T_BDDF1_3D_JointDOF;
      JointDOF = N_T_BDDF1_3D_J;
      N_InnerDOF = N_T_BDDF1_3D_NInner;
      InnerDOF = N_T_BDDF1_3D_Inner;
      N_OuterDOF = N_T_BDDF1_3D_NOuter;
      OuterDOF = N_T_BDDF1_3D_Outer;
      break;
    case FE_N_T_BDDF2_3D:
      N_DOF = N_T_BDDF2_3D_NDOF;
      N_JointDOF = N_T_BDDF2_3D_JointDOF;
      JointDOF = N_T_BDDF2_3D_J;
      N_InnerDOF = N_T_BDDF2_3D_NInner;
      InnerDOF = N_T_BDDF2_3D_Inner;
      N_OuterDOF = N_T_BDDF2_3D_NOuter;
      OuterDOF = N_T_BDDF2_3D_Outer;
      break;
    case FE_N_T_BDDF3_3D:
      N_DOF = N_T_BDDF3_3D_NDOF;
      N_JointDOF = N_T_BDDF3_3D_JointDOF;
      JointDOF = N_T_BDDF3_3D_J;
      N_InnerDOF = N_T_BDDF3_3D_NInner;
      InnerDOF = N_T_BDDF3_3D_Inner;
      N_OuterDOF = N_T_BDDF3_3D_NOuter;
      OuterDOF = N_T_BDDF3_3D_Outer;
      break;
    case FE_N_H_BDDF1_3D:
      N_DOF = N_H_BDDF1_3D_NDOF;
      N_JointDOF = N_H_BDDF1_3D_JointDOF;
      JointDOF = N_H_BDDF1_3D_J;
      N_InnerDOF = N_H_BDDF1_3D_NInner;
      InnerDOF = N_H_BDDF1_3D_Inner;
      N_OuterDOF = N_H_BDDF1_3D_NOuter;
      OuterDOF = N_H_BDDF1_3D_Outer;
      break;
    case FE_N_H_BDDF2_3D:
      N_DOF = N_H_BDDF2_3D_NDOF;
      N_JointDOF = N_H_BDDF2_3D_JointDOF;
      JointDOF = N_H_BDDF2_3D_J;
      N_InnerDOF = N_H_BDDF2_3D_NInner;
      InnerDOF = N_H_BDDF2_3D_Inner;
      N_OuterDOF = N_H_BDDF2_3D_NOuter;
      OuterDOF = N_H_BDDF2_3D_Outer;
      break;
    case FE_N_H_BDDF3_3D:
      N_DOF = N_H_BDDF3_3D_NDOF;
      N_JointDOF = N_H_BDDF3_3D_JointDOF;
      JointDOF = N_H_BDDF3_3D_J;
      N_InnerDOF = N_H_BDDF3_3D_NInner;
      InnerDOF = N_H_BDDF3_3D_Inner;
      N_OuterDOF = N_H_BDDF3_3D_NOuter;
      OuterDOF = N_H_BDDF3_3D_Outer;
      break;
    default:
      ErrThrow("unknown FEDesc3D ", type);
      break;
  }
}

int FEDescriptor::GetJointOfThisDOF(int localDOF) const
{
  bool is_DOF_on_edge = false;
  for (int i = 0; i < N_OuterDOF; i++)
  {
    if(OuterDOF[i]==localDOF)
    {
      is_DOF_on_edge=true;
      break;
    }
  }
  if(!is_DOF_on_edge)
    return -1;
  //else // continue to find the edge
  int i = 0;
  while (true)
  {
    // this must terminate, since we already know localDOF is a dof on an edge
    for (int j = 0; j < N_JointDOF; j++)
    {
      if(JointDOF[i][j] == localDOF)
        return i;
    }
    i++;
  }
}
