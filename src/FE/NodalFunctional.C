#include "NodalFunctional.h"
#include "MooNMD_Io.h"
#include "NF_C_L_P0_1D.h"
#include "NF_C_L_P1_1D.h"
#include "NF_C_L_P2_1D.h"
#include "NF_C_L_P3_1D.h"
#include "NF_D_L_P1_1D.h"
#include "NF_D_L_P2_1D.h"
#include "AllNodalFunctionals2D.h"
#include "AllNodalFunctionals3D.h"

NodalFunctional::NodalFunctional ( NodalFunctional_type id)
: type(id), N_AllFunctionals{}, N_EdgeFunctionals{}, N_PointsAll{}, Xi{nullptr},
  Eta{nullptr}, EvalAll{nullptr}, N_PointsEdge{}, T{nullptr}, EvalJoint {nullptr},
  N_FaceFunctionals{}, N_PointsFace{}, XiArray{nullptr}, EtaArray{nullptr},
  ZetaArray{nullptr}, S{nullptr}
{
  switch(type)
  {
    case NF_C_L_P0_1D:
      N_AllFunctionals = 1;
      N_EdgeFunctionals = 0;
      N_PointsAll = 1;
      N_PointsEdge = 0;
      Xi = NF_C_L_P0_1D_Xi;
      Eta = NF_C_L_P0_1D_Eta;
      T = NF_C_L_P0_1D_T;
      EvalAll = NF_C_L_P0_1D_EvalAll;
      EvalJoint = NF_C_L_P0_1D_EvalEdge;
      break;
    case NF_C_L_P1_1D:
      N_AllFunctionals = 2;
      N_EdgeFunctionals = 2;
      N_PointsAll = 2;
      N_PointsEdge = 2;
      Xi = NF_C_L_P1_1D_Xi;
      Eta = NF_C_L_P1_1D_Eta;
      T = NF_C_L_P1_1D_T;
      EvalAll = NF_C_L_P1_1D_EvalAll;
      EvalJoint = NF_C_L_P1_1D_EvalEdge;
      break;
    case NF_C_L_P2_1D:
      N_AllFunctionals = 3;
      N_EdgeFunctionals = 2;
      N_PointsAll = 3;
      N_PointsEdge = 2;
      Xi = NF_C_L_P2_1D_Xi;
      Eta = NF_C_L_P2_1D_Eta;
      T = NF_C_L_P2_1D_T;
      EvalAll = NF_C_L_P2_1D_EvalAll;
      EvalJoint = NF_C_L_P2_1D_EvalEdge;
      break;
    case NF_C_L_P3_1D:
      N_AllFunctionals = 4;
      N_EdgeFunctionals = 2;
      N_PointsAll = 4;
      N_PointsEdge = 2;
      Xi = NF_C_L_P3_1D_Xi;
      Eta = NF_C_L_P3_1D_Eta;
      T = NF_C_L_P3_1D_T;
      EvalAll = NF_C_L_P3_1D_EvalAll;
      EvalJoint = NF_C_L_P3_1D_EvalEdge;
      break;
    case NF_D_L_P1_1D:
      N_AllFunctionals = 2;
      N_EdgeFunctionals = 0;
      N_PointsAll = 2;
      N_PointsEdge = 0;
      Xi = NF_D_L_P1_1D_Xi;
      Eta = NF_D_L_P1_1D_Eta;
      T = NF_D_L_P1_1D_T;
      EvalAll = NF_D_L_P1_1D_EvalAll;
      EvalJoint = NF_D_L_P1_1D_EvalEdge;
      break;
    case NF_D_L_P2_1D:
      N_AllFunctionals = 3;
      N_EdgeFunctionals = 0;
      N_PointsAll = 3;
      N_PointsEdge = 0;
      Xi = NF_D_L_P2_1D_Xi;
      Eta = NF_D_L_P2_1D_Eta;
      T = NF_D_L_P2_1D_T;
      EvalAll = NF_D_L_P2_1D_EvalAll;
      EvalJoint = NF_D_L_P2_1D_EvalEdge;
      break;
    case NF_C_T_P00_2D:
      N_AllFunctionals = 1;
      N_EdgeFunctionals = 0;
      N_PointsAll = 1;
      N_PointsEdge = 0;
      Xi = NF_C_T_P00_2D_Xi;
      Eta = NF_C_T_P00_2D_Eta;
      T = NF_C_T_P00_2D_T;
      EvalAll = NF_C_T_P00_2D_EvalAll;
      EvalJoint = NF_C_T_P00_2D_EvalEdge;
      break;
    case NF_C_T_P0_2D:
      N_AllFunctionals = 1;
      N_EdgeFunctionals = 0;
      N_PointsAll = 1;
      N_PointsEdge = 0;
      Xi = NF_C_T_P0_2D_Xi;
      Eta = NF_C_T_P0_2D_Eta;
      T = NF_C_T_P0_2D_T;
      EvalAll = NF_C_T_P0_2D_EvalAll;
      EvalJoint = NF_C_T_P0_2D_EvalEdge;
      break;
    case NF_C_T_P1_2D:
      N_AllFunctionals = 3;
      N_EdgeFunctionals = 2;
      N_PointsAll = 3;
      N_PointsEdge = 2;
      Xi = NF_C_T_P1_2D_Xi;
      Eta = NF_C_T_P1_2D_Eta;
      T = NF_C_T_P1_2D_T;
      EvalAll = NF_C_T_P1_2D_EvalAll;
      EvalJoint = NF_C_T_P1_2D_EvalEdge;
      break;
    case NF_C_T_P2_2D:
      N_AllFunctionals = 6;
      N_EdgeFunctionals = 3;
      N_PointsAll = 6;
      N_PointsEdge = 3;
      Xi = NF_C_T_P2_2D_Xi;
      Eta = NF_C_T_P2_2D_Eta;
      T = NF_C_T_P2_2D_T;
      EvalAll = NF_C_T_P2_2D_EvalAll;
      EvalJoint = NF_C_T_P2_2D_EvalEdge;
      break;
    case NF_C_T_P3_2D:
      N_AllFunctionals = 10;
      N_EdgeFunctionals = 4;
      N_PointsAll = 10;
      N_PointsEdge = 4;
      Xi = NF_C_T_P3_2D_Xi;
      Eta = NF_C_T_P3_2D_Eta;
      T = NF_C_T_P3_2D_T;
      EvalAll = NF_C_T_P3_2D_EvalAll;
      EvalJoint = NF_C_T_P3_2D_EvalEdge;
      break;
    case NF_C_T_P4_2D:
      N_AllFunctionals = 15;
      N_EdgeFunctionals = 5;
      N_PointsAll = 15;
      N_PointsEdge = 5;
      Xi = NF_C_T_P4_2D_Xi;
      Eta = NF_C_T_P4_2D_Eta;
      T = NF_C_T_P4_2D_T;
      EvalAll = NF_C_T_P4_2D_EvalAll;
      EvalJoint = NF_C_T_P4_2D_EvalEdge;
      break;
    case NF_C_T_P5_2D:
      N_AllFunctionals = 21;
      N_EdgeFunctionals = 6;
      N_PointsAll = 21;
      N_PointsEdge = 6;
      Xi = NF_C_T_P5_2D_Xi;
      Eta = NF_C_T_P5_2D_Eta;
      T = NF_C_T_P5_2D_T;
      EvalAll = NF_C_T_P5_2D_EvalAll;
      EvalJoint = NF_C_T_P5_2D_EvalEdge;
      break;
    case NF_C_T_P6_2D:
      N_AllFunctionals = 28;
      N_EdgeFunctionals = 7;
      N_PointsAll = 28;
      N_PointsEdge = 7;
      Xi = NF_C_T_P6_2D_Xi;
      Eta = NF_C_T_P6_2D_Eta;
      T = NF_C_T_P6_2D_T;
      EvalAll = NF_C_T_P6_2D_EvalAll;
      EvalJoint = NF_C_T_P6_2D_EvalEdge;
      break;
    case NF_C_T_P7_2D:
      ErrThrow("Nodal functional of type NF_C_T_P7_2D not implemented");
      break;
    case NF_C_T_P8_2D:
      ErrThrow("Nodal functional of type NF_C_T_P8_2D not implemented");
      break;
    case NF_C_T_P9_2D:
      ErrThrow("Nodal functional of type NF_C_T_P9_2D not implemented");
      break;
    case NF_N_T_P1_2D:
      N_AllFunctionals = 3;
      N_EdgeFunctionals = 1;
      N_PointsAll = 9;
      N_PointsEdge = 3;
      Xi = NF_N_T_P1_2D_Xi;
      Eta = NF_N_T_P1_2D_Eta;
      T = NF_N_T_P1_2D_T;
      EvalAll = NF_N_T_P1_2D_EvalAll;
      EvalJoint = NF_N_T_P1_2D_EvalEdge;
      break;
    case NF_C_Q_Q00_2D:
      N_AllFunctionals = 1;
      N_EdgeFunctionals = 0;
      N_PointsAll = 1;
      N_PointsEdge = 0;
      Xi = NF_C_Q_Q00_2D_Xi;
      Eta = NF_C_Q_Q00_2D_Eta;
      T = NF_C_Q_Q00_2D_T;
      EvalAll = NF_C_Q_Q00_2D_EvalAll;
      EvalJoint = NF_C_Q_Q00_2D_EvalEdge;
      break;
    case NF_C_Q_Q0_2D:
      N_AllFunctionals = 1;
      N_EdgeFunctionals = 0;
      N_PointsAll = 1;
      N_PointsEdge = 0;
      Xi = NF_C_Q_Q0_2D_Xi;
      Eta = NF_C_Q_Q0_2D_Eta;
      T = NF_C_Q_Q0_2D_T;
      EvalAll = NF_C_Q_Q0_2D_EvalAll;
      EvalJoint = NF_C_Q_Q0_2D_EvalEdge;
      break;
    case NF_C_Q_Q1_2D:
      N_AllFunctionals = 4;
      N_EdgeFunctionals = 2;
      N_PointsAll = 4;
      N_PointsEdge = 2;
      Xi = NF_C_Q_Q1_2D_Xi;
      Eta = NF_C_Q_Q1_2D_Eta;
      T = NF_C_Q_Q1_2D_T;
      EvalAll = NF_C_Q_Q1_2D_EvalAll;
      EvalJoint = NF_C_Q_Q1_2D_EvalEdge;
      break;
    case NF_C_Q_Q2_2D:
      N_AllFunctionals = 9;
      N_EdgeFunctionals = 3;
      N_PointsAll = 9;
      N_PointsEdge = 3;
      Xi = NF_C_Q_Q2_2D_Xi;
      Eta = NF_C_Q_Q2_2D_Eta;
      T = NF_C_Q_Q2_2D_T;
      EvalAll = NF_C_Q_Q2_2D_EvalAll;
      EvalJoint = NF_C_Q_Q2_2D_EvalEdge;
      break;
    case NF_C_Q_Q3_2D:
      N_AllFunctionals = 16;
      N_EdgeFunctionals = 4;
      N_PointsAll = 16;
      N_PointsEdge = 4;
      Xi = NF_C_Q_Q3_2D_Xi;
      Eta = NF_C_Q_Q3_2D_Eta;
      T = NF_C_Q_Q3_2D_T;
      EvalAll = NF_C_Q_Q3_2D_EvalAll;
      EvalJoint = NF_C_Q_Q3_2D_EvalEdge;
      break;
    case NF_C_Q_Q4_2D:
      N_AllFunctionals = 25;
      N_EdgeFunctionals = 5;
      N_PointsAll = 25;
      N_PointsEdge = 5;
      Xi = NF_C_Q_Q4_2D_Xi;
      Eta = NF_C_Q_Q4_2D_Eta;
      T = NF_C_Q_Q4_2D_T;
      EvalAll = NF_C_Q_Q4_2D_EvalAll;
      EvalJoint = NF_C_Q_Q4_2D_EvalEdge;
      break;
    case NF_C_Q_Q5_2D:
      N_AllFunctionals = 36;
      N_EdgeFunctionals = 6;
      N_PointsAll = 36;
      N_PointsEdge = 6;
      Xi = NF_C_Q_Q5_2D_Xi;
      Eta = NF_C_Q_Q5_2D_Eta;
      T = NF_C_Q_Q5_2D_T;
      EvalAll = NF_C_Q_Q5_2D_EvalAll;
      EvalJoint = NF_C_Q_Q5_2D_EvalEdge;
      break;
    case NF_C_Q_Q6_2D:
      N_AllFunctionals = 49;
      N_EdgeFunctionals = 7;
      N_PointsAll = 49;
      N_PointsEdge = 49;
      Xi = NF_C_Q_Q6_2D_Xi;
      Eta = NF_C_Q_Q6_2D_Eta;
      T = NF_C_Q_Q6_2D_T;
      EvalAll = NF_C_Q_Q6_2D_EvalAll;
      EvalJoint = NF_C_Q_Q6_2D_EvalEdge;
      break;
    case NF_C_Q_Q7_2D:
      N_AllFunctionals = 64;
      N_EdgeFunctionals = 8;
      N_PointsAll = 64;
      N_PointsEdge = 8;
      Xi = NF_C_Q_Q7_2D_Xi;
      Eta = NF_C_Q_Q7_2D_Eta;
      T = NF_C_Q_Q7_2D_T;
      EvalAll = NF_C_Q_Q7_2D_EvalAll;
      EvalJoint = NF_C_Q_Q7_2D_EvalEdge;
      break;
    case NF_C_Q_Q8_2D:
      N_AllFunctionals = 81;
      N_EdgeFunctionals = 9;
      N_PointsAll = 81;
      N_PointsEdge = 9;
      Xi = NF_C_Q_Q8_2D_Xi;
      Eta = NF_C_Q_Q8_2D_Eta;
      T = NF_C_Q_Q8_2D_T;
      EvalAll = NF_C_Q_Q8_2D_EvalAll;
      EvalJoint = NF_C_Q_Q8_2D_EvalEdge;
      break;
    case NF_C_Q_Q9_2D:
      N_AllFunctionals = 10;
      N_EdgeFunctionals = 100;
      N_PointsAll = 10;
      N_PointsEdge = 100;
      Xi = NF_C_Q_Q9_2D_Xi;
      Eta = NF_C_Q_Q9_2D_Eta;
      T = NF_C_Q_Q9_2D_T;
      EvalAll = NF_C_Q_Q9_2D_EvalAll;
      EvalJoint = NF_C_Q_Q9_2D_EvalEdge;
      break;
    case NF_N_Q_Q1_2D:
      N_AllFunctionals = 4;
      N_EdgeFunctionals = 1;
      N_PointsAll = 12;
      N_PointsEdge = 3;
      Xi = NF_N_Q_Q1_2D_Xi;
      Eta = NF_N_Q_Q1_2D_Eta;
      T = NF_N_Q_Q1_2D_T;
      EvalAll = NF_N_Q_Q1_2D_EvalAll;
      EvalJoint = NF_N_Q_Q1_2D_EvalEdge;
      break;
    case NF_D_Q_P1_2D:
      N_AllFunctionals = 3;
      N_EdgeFunctionals = 0;
      N_PointsAll = 9;
      N_PointsEdge = 0;
      Xi = NF_D_Q_P1_2D_Xi;
      Eta = NF_D_Q_P1_2D_Eta;
      T = NF_D_Q_P1_2D_t;
      EvalAll = NF_D_Q_P1_2D_EvalAll;
      EvalJoint = nullptr;
      break;
    case NF_D_Q_P2_2D:
      N_AllFunctionals = 6;
      N_EdgeFunctionals = 0;
      N_PointsAll = 9;
      N_PointsEdge = 0;
      Xi = NF_D_Q_P2_2D_Xi;
      Eta = NF_D_Q_P2_2D_Eta;
      T = NF_D_Q_P2_2D_t;
      EvalAll = NF_D_Q_P2_2D_EvalAll;
      EvalJoint = nullptr;
      break;
    case NF_D_Q_P3_2D:
      N_AllFunctionals = 10;
      N_EdgeFunctionals = 0;
      N_PointsAll = 16;
      N_PointsEdge = 0;
      Xi = NF_D_Q_P3_2D_Weight1;
      Eta = NF_D_Q_P3_2D_Weight2;
      T = NF_D_Q_P3_2D_T;
      EvalAll = NF_D_Q_P3_2D_EvalAll;
      EvalJoint = nullptr;
      break;
    case NF_C_T_B2_2D:
      N_AllFunctionals = 7;
      N_EdgeFunctionals = 3;
      N_PointsAll = 7;
      N_PointsEdge = 3;
      Xi = NF_C_T_B2_2D_Xi;
      Eta = NF_C_T_B2_2D_Eta;
      T = NF_C_T_B2_2D_T;
      EvalAll = NF_C_T_B2_2D_EvalAll;
      EvalJoint = NF_C_T_B2_2D_EvalEdge;
      break;
    case NF_C_T_B3_2D:
      N_AllFunctionals = 12;
      N_EdgeFunctionals = 4;
      N_PointsAll = 16;
      N_PointsEdge = 4;
      Xi = NF_C_T_B3_2D_Xi;
      Eta = NF_C_T_B3_2D_Eta;
      T = NF_C_T_B3_2D_T;
      EvalAll = NF_C_T_B3_2D_EvalAll;
      EvalJoint = NF_C_T_B3_2D_EvalEdge;
      break;
    case NF_C_T_SV2_2D:
      N_AllFunctionals = 10;
      N_EdgeFunctionals = 3;
      N_PointsAll = 10;
      N_PointsEdge = 3;
      Xi = NF_C_T_SV2_2D_Xi;
      Eta = NF_C_T_SV2_2D_Eta;
      T = NF_C_T_SV2_2D_T;
      EvalAll = NF_C_T_SV2_2D_EvalAll;
      EvalJoint = NF_C_T_SV2_2D_EvalEdge;
      break;
    case NF_D_T_P1_2D:
      N_AllFunctionals = 3;
      N_EdgeFunctionals = 0;
      N_PointsAll = 3;
      N_PointsEdge = 0;
      Xi = NF_D_T_P1_2D_Xi;
      Eta = NF_D_T_P1_2D_Eta;
      T = NF_D_T_P1_2D_T_P;
      EvalAll = NF_D_T_P1_2D_EvalAll;
      EvalJoint = NF_D_T_P1_2D_EvalEdge;
      break;
    case NF_D_T_P2_2D:
      N_AllFunctionals = 6;
      N_EdgeFunctionals = 0;
      N_PointsAll = 7;
      N_PointsEdge = 0;
      Xi = NF_D_T_P2_2D_Xi;
      Eta = NF_D_T_P2_2D_Eta;
      T = NF_D_T_P2_2D_T_P;
      EvalAll = NF_D_T_P2_2D_EvalAll;
      EvalJoint = NF_D_T_P2_2D_EvalEdge;
      break;
    case NF_D_T_SV1_2D:
      N_AllFunctionals = 9;
      N_EdgeFunctionals = 0;
      N_PointsAll = 9;
      N_PointsEdge = 0;
      Xi = NF_D_T_SV1_2D_Xi;
      Eta = NF_D_T_SV1_2D_Eta;
      T = NF_D_T_SV1_2D_T;
      EvalAll = NF_D_T_SV1_2D_EvalAll;
      EvalJoint = NF_D_T_SV1_2D_EvalEdge;
      break;
    case NF_N_Q_Q2_2D:
      N_AllFunctionals = 9;
      N_EdgeFunctionals = 2;
      N_PointsAll = 21;
      N_PointsEdge = 3;
      Xi = NF_N_Q_Q2_2D_Xi;
      Eta = NF_N_Q_Q2_2D_Eta;
      T = NF_N_Q_Q2_2D_T;
      EvalAll = NF_N_Q_Q2_2D_EvalAll;
      EvalJoint = NF_N_Q_Q2_2D_EvalEdge;
      break;
    case NF_N_Q_Q3_2D:
      N_AllFunctionals = 15;
      N_EdgeFunctionals = 3;
      N_PointsAll = 32;
      N_PointsEdge = 4;
      Xi = NF_N_Q_Q3_2D_Xi;
      Eta = NF_N_Q_Q3_2D_Eta;
      T = NF_N_Q_Q3_2D_T;
      EvalAll = NF_N_Q_Q3_2D_EvalAll;
      EvalJoint = NF_N_Q_Q3_2D_EvalEdge;
      break;
    case NF_N_Q_Q4_2D:
      N_AllFunctionals = 22;
      N_EdgeFunctionals = 4;
      N_PointsAll = 45;
      N_PointsEdge = 5;
      Xi = NF_N_Q_Q4_2D_Xi;
      Eta = NF_N_Q_Q4_2D_Eta;
      T = NF_N_Q_Q4_2D_T;
      EvalAll = NF_N_Q_Q4_2D_EvalAll;
      EvalJoint = NF_N_Q_Q4_2D_EvalEdge;
      break;
    case NF_N_Q_Q5_2D:
      N_AllFunctionals = 30;
      N_EdgeFunctionals = 5;
      N_PointsAll = 60;
      N_PointsEdge = 6;
      Xi = NF_N_Q_Q5_2D_Xi;
      Eta = NF_N_Q_Q5_2D_Eta;
      T = NF_N_Q_Q5_2D_T;
      EvalAll = NF_N_Q_Q5_2D_EvalAll;
      EvalJoint = NF_N_Q_Q5_2D_EvalEdge;
      break;
    case NF_D_Q_P4_2D:
      N_AllFunctionals = 15;
      N_EdgeFunctionals = 0;
      N_PointsAll = 25;
      N_PointsEdge = 0;
      Xi = NF_D_Q_P4_2D_Xi;
      Eta = NF_D_Q_P4_2D_Eta;
      T = NF_D_Q_P4_2D_T;
      EvalAll = NF_D_Q_P4_2D_EvalAll;
      EvalJoint = nullptr;
      break;
    case NF_D_Q_P5_2D:
      N_AllFunctionals = 21;
      N_EdgeFunctionals = 0;
      N_PointsAll = 36;
      N_PointsEdge = 0;
      Xi = NF_D_Q_P5_2D_Xi;
      Eta = NF_D_Q_P5_2D_Eta;
      T = NF_D_Q_P5_2D_T;
      EvalAll = NF_D_Q_P5_2D_EvalAll;
      EvalJoint = nullptr;
      break;
    case NF_D_Q_P6_2D:
      N_AllFunctionals = 28;
      N_EdgeFunctionals = 0;
      N_PointsAll = 49;
      N_PointsEdge = 0;
      Xi = NF_D_Q_P6_2D_Xi;
      Eta = NF_D_Q_P6_2D_Eta;
      T = NF_D_Q_P6_2D_T;
      EvalAll = NF_D_Q_P6_2D_EvalAll;
      EvalJoint = nullptr;
      break;
    case NF_D_Q_P7_2D:
      N_AllFunctionals = 36;
      N_EdgeFunctionals = 0;
      N_PointsAll = 64;
      N_PointsEdge = 0;
      Xi = NF_D_Q_P7_2D_Xi;
      Eta = NF_D_Q_P7_2D_Eta;
      T = NF_D_Q_P7_2D_T;
      EvalAll = NF_D_Q_P7_2D_EvalAll;
      EvalJoint = nullptr;
      break;
    case NF_N_T_P1MOD_2D:
      N_AllFunctionals = 6;
      N_EdgeFunctionals = 2;
      N_PointsAll = 9;
      N_PointsEdge = 3;
      Xi = NF_N_T_P1MOD_2D_Xi;
      Eta = NF_N_T_P1MOD_2D_Eta;
      T = NF_N_T_P1MOD_2D_T;
      EvalAll = NF_N_T_P1MOD_2D_EvalAll;
      EvalJoint = NF_N_T_P1MOD_2D_EvalEdge;
      break;
    case NF_N_T_P2_2D:
      N_AllFunctionals = 7;
      N_EdgeFunctionals = 2;
      N_PointsAll = 16;
      N_PointsEdge = 3;
      Xi = NF_N_T_P2_2D_Xi;
      Eta = NF_N_T_P2_2D_Eta;
      T = NF_N_T_P2_2D_T;
      EvalAll = NF_N_T_P2_2D_EvalAll;
      EvalJoint = NF_N_T_P2_2D_EvalEdge;
      break;
    case NF_N_T_P3_2D:
      N_AllFunctionals = 12;
      N_EdgeFunctionals = 3;
      N_PointsAll = 19;
      N_PointsEdge = 4;
      Xi = NF_N_T_P3_2D_Xi;
      Eta = NF_N_T_P3_2D_Eta;
      T = NF_N_T_P3_2D_T;
      EvalAll = NF_N_T_P3_2D_EvalAll;
      EvalJoint = NF_N_T_P3_2D_EvalEdge;
      break;
    case NF_N_T_P4_2D:
      N_AllFunctionals = 18;
      N_EdgeFunctionals = 4;
      N_PointsAll = 30;
      N_PointsEdge = 5;
      Xi = NF_N_T_P4_2D_Xi;
      Eta = NF_N_T_P4_2D_Eta;
      T = NF_N_T_P4_2D_T;
      EvalAll = NF_N_T_P4_2D_EvalAll;
      EvalJoint = NF_N_T_P4_2D_EvalEdge;
      break;
    case NF_N_T_P5_2D:
      N_AllFunctionals = 25;
      N_EdgeFunctionals = 5;
      N_PointsAll = 45;
      N_PointsEdge = 6;
      Xi = NF_N_T_P5_2D_Xi;
      Eta = NF_N_T_P5_2D_Eta;
      T = NF_N_T_P5_2D_T;
      EvalAll = NF_N_T_P5_2D_EvalAll;
      EvalJoint = NF_N_T_P5_2D_EvalEdge;
      break;
    case NF_D_T_P3_2D:
      N_AllFunctionals = 10;
      N_EdgeFunctionals = 0;
      N_PointsAll = 15;
      N_PointsEdge = 0;
      Xi = NF_D_T_P3_2D_Xi;
      Eta = NF_D_T_P3_2D_Eta;
      T = nullptr;
      EvalAll = NF_D_T_P3_2D_EvalAll;
      EvalJoint = NF_D_T_P3_2D_EvalEdge;
      break;
    case NF_D_T_P4_2D:
      N_AllFunctionals = 15;
      N_EdgeFunctionals = 0;
      N_PointsAll = 15;
      N_PointsEdge = 0;
      Xi = NF_D_T_P4_2D_Xi;
      Eta = NF_D_T_P4_2D_Eta;
      T = nullptr;
      EvalAll = NF_D_T_P4_2D_EvalAll;
      EvalJoint = NF_D_T_P4_2D_EvalEdge;
      break;
    case NF_C_T_P1MINI_2D:
      N_AllFunctionals = 4;
      N_EdgeFunctionals = 2;
      N_PointsAll = 7;
      N_PointsEdge = 2;
      Xi = NF_C_T_P1MINI_2D_Xi;
      Eta = NF_C_T_P1MINI_2D_Eta;
      T = NF_C_T_P1MINI_2D_T;
      EvalAll = NF_C_T_P1MINI_2D_EvalAll;
      EvalJoint = NF_C_T_P1MINI_2D_EvalEdge;
      break;
    case NF_B_Q_IB2_2D:
      N_AllFunctionals = 1;
      N_EdgeFunctionals = 0;
      N_PointsAll = 1;
      N_PointsEdge = 0;
      Xi = NF_B_Q_IB2_2D_Xi;
      Eta = NF_B_Q_IB2_2D_Eta;
      T = NF_B_Q_IB2_2D_T;
      EvalAll = NF_B_Q_IB2_2D_EvalAll;
      EvalJoint = NF_B_Q_IB2_2D_EvalEdge;
      break;
    case NF_S_Q_Q2_2D:
      N_AllFunctionals = 9;
      N_EdgeFunctionals = 3;
      N_PointsAll = 25;
      N_PointsEdge = 5;
      Xi = NF_S_Q_Q2_2D_Xi;
      Eta = NF_S_Q_Q2_2D_Eta;
      T = NF_S_Q_Q2_2D_T;
      EvalAll = NF_S_Q_Q2_2D_EvalAll;
      EvalJoint = NF_S_Q_Q2_2D_EvalEdge;
      break;
    case NF_D_Q_Q1_2D:
      N_AllFunctionals = 4;
      N_EdgeFunctionals = 0;
      N_PointsAll = 9;
      N_PointsEdge = 0;
      Xi = NF_D_Q_Q1_2D_Xi;
      Eta = NF_D_Q_Q1_2D_Eta;
      T = NF_D_Q_Q1_2D_t;
      EvalAll = NF_D_Q_Q1_2D_EvalAll;
      EvalJoint = nullptr;
      break;
    case NF_D_Q_Q2_2D:
      N_AllFunctionals = 9;
      N_EdgeFunctionals = 0;
      N_PointsAll = 9;
      N_PointsEdge = 0;
      Xi = NF_D_Q_Q2_2D_Xi;
      Eta = NF_D_Q_Q2_2D_Eta;
      T = NF_D_Q_Q2_2D_t;
      EvalAll = NF_D_Q_Q2_2D_EvalAll;
      EvalJoint = nullptr;
      break;
    case NF_D_Q_Q3_2D:
      N_AllFunctionals = 16;
      N_EdgeFunctionals = 0;
      N_PointsAll = 16;
      N_PointsEdge = 0;
      Xi = NF_D_Q_Q3_2D_Xi;
      Eta = NF_D_Q_Q3_2D_Eta;
      T = NF_D_Q_Q3_2D_T;
      EvalAll = NF_D_Q_Q3_2D_EvalAll;
      EvalJoint = nullptr;
      break;
    case NF_D_Q_Q4_2D:
      N_AllFunctionals = 25;
      N_EdgeFunctionals = 0;
      N_PointsAll = 25;
      N_PointsEdge = 0;
      Xi = NF_D_Q_Q4_2D_Xi;
      Eta = NF_D_Q_Q4_2D_Eta;
      T = NF_D_Q_Q4_2D_T;
      EvalAll = NF_D_Q_Q4_2D_EvalAll;
      EvalJoint = nullptr;
      break;
    case NF_C_T_B4_2D:
      N_AllFunctionals = 18;
      N_EdgeFunctionals = 5;
      N_PointsAll = 18;
      N_PointsEdge = 5;
      Xi = NF_C_T_B4_2D_Xi;
      Eta = NF_C_T_B4_2D_Eta;
      T = NF_C_T_B4_2D_T;
      EvalAll = NF_C_T_B4_2D_EvalAll;
      EvalJoint = NF_C_T_B4_2D_EvalEdge;
      break;
    case NF_D_Q_D2_2D:
      N_AllFunctionals = 9;
      N_EdgeFunctionals = 0;
      N_PointsAll = 16;
      N_PointsEdge = 0;
      Xi = NF_D_Q_D2_2D_Xi;
      Eta = NF_D_Q_D2_2D_Eta;
      T = NF_D_Q_D2_2D_t;
      EvalAll = NF_D_Q_D2_2D_EvalAll;
      EvalJoint = nullptr;
      break;
    case NF_C_Q_UL1_2D:
      N_AllFunctionals = 5;
      N_EdgeFunctionals = 2;
      N_PointsAll = 8;
      N_PointsEdge = 2;
      Xi = NF_C_Q_UL1_2D_Xi;
      Eta = NF_C_Q_UL1_2D_Eta;
      T = NF_C_Q_UL1_2D_T;
      EvalAll = NF_C_Q_UL1_2D_EvalAll;
      EvalJoint = NF_C_Q_UL1_2D_EvalEdge;
      break;
    case NF_C_Q_UL2_2D:
      N_AllFunctionals = 11;
      N_EdgeFunctionals = 3;
      N_PointsAll = 17;
      N_PointsEdge = 3;
      Xi = NF_C_Q_UL2_2D_Xi;
      Eta = NF_C_Q_UL2_2D_Eta;
      T = NF_C_Q_UL2_2D_T;
      EvalAll = NF_C_Q_UL2_2D_EvalAll;
      EvalJoint = NF_C_Q_UL2_2D_EvalEdge;
      break;
    case NF_C_Q_UL3_2D:
      N_AllFunctionals = 18;
      N_EdgeFunctionals = 4;
      N_PointsAll = 28;
      N_PointsEdge = 4;
      Xi = NF_C_Q_UL3_2D_Xi;
      Eta = NF_C_Q_UL3_2D_Eta;
      T = NF_C_Q_UL3_2D_T;
      EvalAll = NF_C_Q_UL3_2D_EvalAll;
      EvalJoint = NF_C_Q_UL3_2D_EvalEdge;
      break;
    case NF_C_Q_UL4_2D:
      N_AllFunctionals = 27;
      N_EdgeFunctionals = 5;
      N_PointsAll = 41;
      N_PointsEdge = 5;
      Xi = NF_C_Q_UL4_2D_Xi;
      Eta = NF_C_Q_UL4_2D_Eta;
      T = NF_C_Q_UL4_2D_T;
      EvalAll = NF_C_Q_UL4_2D_EvalAll;
      EvalJoint = NF_C_Q_UL4_2D_EvalEdge;
      break;
    case NF_C_Q_UL5_2D:
      N_AllFunctionals = 38;
      N_EdgeFunctionals = 6;
      N_PointsAll = 56;
      N_PointsEdge = 6;
      Xi = NF_C_Q_UL5_2D_Xi;
      Eta = NF_C_Q_UL5_2D_Eta;
      T = NF_C_Q_UL5_2D_T;
      EvalAll = NF_C_Q_UL5_2D_EvalAll;
      EvalJoint = NF_C_Q_UL5_2D_EvalEdge;
      break;
    case NF_C_T_UL1_2D:
      N_AllFunctionals = 4;
      N_EdgeFunctionals = 2;
      N_PointsAll = 10;
      N_PointsEdge = 2;
      Xi = NF_C_T_UL1_2D_Xi;
      Eta = NF_C_T_UL1_2D_Eta;
      T = NF_C_T_UL1_2D_T;
      EvalAll = NF_C_T_UL1_2D_EvalAll;
      EvalJoint = NF_C_T_UL1_2D_EvalEdge;
      break;
    case NF_C_T_UL2_2D:
      N_AllFunctionals = 9;
      N_EdgeFunctionals = 3;
      N_PointsAll = 13;
      N_PointsEdge = 3;
      Xi = NF_C_T_UL2_2D_Xi;
      Eta = NF_C_T_UL2_2D_Eta;
      T = NF_C_T_UL2_2D_T;
      EvalAll = NF_C_T_UL2_2D_EvalAll;
      EvalJoint = NF_C_T_UL2_2D_EvalEdge;
      break;
    case NF_C_T_UL3_2D:
      N_AllFunctionals = 15;
      N_EdgeFunctionals = 4;
      N_PointsAll = 24;
      N_PointsEdge = 4;
      Xi = NF_C_T_UL3_2D_Xi;
      Eta = NF_C_T_UL3_2D_Eta;
      T = NF_C_T_UL3_2D_T;
      EvalAll = NF_C_T_UL3_2D_EvalAll;
      EvalJoint = NF_C_T_UL3_2D_EvalEdge;
      break;
    case NF_C_T_UL4_2D:
      N_AllFunctionals = 22;
      N_EdgeFunctionals = 5;
      N_PointsAll = 39;
      N_PointsEdge = 5;
      Xi = NF_C_T_UL4_2D_Xi;
      Eta = NF_C_T_UL4_2D_Eta;
      T = NF_C_T_UL4_2D_T;
      EvalAll = NF_C_T_UL4_2D_EvalAll;
      EvalJoint = NF_C_T_UL4_2D_EvalEdge;
      break;
    case NF_C_T_UL5_2D:
      N_AllFunctionals = 30;
      N_EdgeFunctionals = 6;
      N_PointsAll = 42;
      N_PointsEdge = 6;
      Xi = NF_C_T_UL5_2D_Xi;
      Eta = NF_C_T_UL5_2D_Eta;
      T = NF_C_T_UL5_2D_T;
      EvalAll = NF_C_T_UL5_2D_EvalAll;
      EvalJoint = NF_C_T_UL5_2D_EvalEdge;
      break;
    case NF_C_Q_UL2S_2D:
      N_AllFunctionals = 9;
      N_EdgeFunctionals = 3;
      N_PointsAll = 9;
      N_PointsEdge = 3;
      Xi = NF_C_Q_UL2S_2D_Xi;
      Eta = NF_C_Q_UL2S_2D_Eta;
      T = NF_C_Q_UL2S_2D_T;
      EvalAll = NF_C_Q_UL2S_2D_EvalAll;
      EvalJoint = NF_C_Q_UL2S_2D_EvalEdge;
      break;
    case NF_C_Q_UL3S_2D:
      N_AllFunctionals = 15;
      N_EdgeFunctionals = 4;
      N_PointsAll = 15;
      N_PointsEdge = 4;
      Xi = NF_C_Q_UL3S_2D_Xi;
      Eta = NF_C_Q_UL3S_2D_Eta;
      T = NF_C_Q_UL3S_2D_T;
      EvalAll = NF_C_Q_UL3S_2D_EvalAll;
      EvalJoint = NF_C_Q_UL3S_2D_EvalEdge;
      break;
    case NF_C_Q_UL4S_2D:
      N_AllFunctionals = 22;
      N_EdgeFunctionals = 5;
      N_PointsAll = 22;
      N_PointsEdge = 5;
      Xi = NF_C_Q_UL4S_2D_Xi;
      Eta = NF_C_Q_UL4S_2D_Eta;
      T = NF_C_Q_UL4S_2D_T;
      EvalAll = NF_C_Q_UL4S_2D_EvalAll;
      EvalJoint = NF_C_Q_UL4S_2D_EvalEdge;
      break;
    case NF_C_Q_UL5S_2D:
      N_AllFunctionals = 31;
      N_EdgeFunctionals = 6;
      N_PointsAll = 31;
      N_PointsEdge = 6;
      Xi = NF_C_Q_UL5S_2D_Xi;
      Eta = NF_C_Q_UL5S_2D_Eta;
      T = NF_C_Q_UL5S_2D_T;
      EvalAll = NF_C_Q_UL5S_2D_EvalAll;
      EvalJoint = NF_C_Q_UL5S_2D_EvalEdge;
      break;
    case NF_C_Q_UL6S_2D:
      N_AllFunctionals = 42;
      N_EdgeFunctionals = 7;
      N_PointsAll = 42;
      N_PointsEdge = 7;
      Xi = NF_C_Q_UL6S_2D_Xi;
      Eta = NF_C_Q_UL6S_2D_Eta;
      T = NF_C_Q_UL6S_2D_T;
      EvalAll = NF_C_Q_UL6S_2D_EvalAll;
      EvalJoint = NF_C_Q_UL6S_2D_EvalEdge;
      break;
    case NF_C_Q_UL7S_2D:
      N_AllFunctionals = 55;
      N_EdgeFunctionals = 8;
      N_PointsAll = 55;
      N_PointsEdge = 8;
      Xi = NF_C_Q_UL7S_2D_Xi;
      Eta = NF_C_Q_UL7S_2D_Eta;
      T = NF_C_Q_UL7S_2D_T;
      EvalAll = NF_C_Q_UL7S_2D_EvalAll;
      EvalJoint = NF_C_Q_UL7S_2D_EvalEdge;
      break;
    case NF_C_Q_UL8S_2D:
      N_AllFunctionals = 70;
      N_EdgeFunctionals = 9;
      N_PointsAll = 70;
      N_PointsEdge = 9;
      Xi = NF_C_Q_UL8S_2D_Xi;
      Eta = NF_C_Q_UL8S_2D_Eta;
      T = NF_C_Q_UL8S_2D_T;
      EvalAll = NF_C_Q_UL8S_2D_EvalAll;
      EvalJoint = NF_C_Q_UL8S_2D_EvalEdge;
      break;
    case NF_C_Q_UL9S_2D:
      N_AllFunctionals = 87;
      N_EdgeFunctionals = 10;
      N_PointsAll = 87;
      N_PointsEdge = 10;
      Xi = NF_C_Q_UL9S_2D_Xi;
      Eta = NF_C_Q_UL9S_2D_Eta;
      T = NF_C_Q_UL9S_2D_T;
      EvalAll = NF_C_Q_UL9S_2D_EvalAll;
      EvalJoint = NF_C_Q_UL9S_2D_EvalEdge;
      break;
    case NF_C_Q_UL2SE_2D:
      N_AllFunctionals = 8;
      N_EdgeFunctionals = 3;
      N_PointsAll = 8;
      N_PointsEdge = 3;
      Xi = NF_C_Q_UL2SE_2D_Xi;
      Eta = NF_C_Q_UL2SE_2D_Eta;
      T = NF_C_Q_UL2SE_2D_T;
      EvalAll = NF_C_Q_UL2SE_2D_EvalAll;
      EvalJoint = NF_C_Q_UL2SE_2D_EvalEdge;
      break;
    case NF_C_Q_UL3SE_2D:
      N_AllFunctionals = 13;
      N_EdgeFunctionals = 4;
      N_PointsAll = 13;
      N_PointsEdge = 4;
      Xi = NF_C_Q_UL3SE_2D_Xi;
      Eta = NF_C_Q_UL3SE_2D_Eta;
      T = NF_C_Q_UL3SE_2D_T;
      EvalAll = NF_C_Q_UL3SE_2D_EvalAll;
      EvalJoint = NF_C_Q_UL3SE_2D_EvalEdge;
      break;
    case NF_C_Q_UL4SE_2D:
      N_AllFunctionals = 20;
      N_EdgeFunctionals = 5;
      N_PointsAll = 20;
      N_PointsEdge = 5;
      Xi = NF_C_Q_UL4SE_2D_Xi;
      Eta = NF_C_Q_UL4SE_2D_Eta;
      T = NF_C_Q_UL4SE_2D_T;
      EvalAll = NF_C_Q_UL4SE_2D_EvalAll;
      EvalJoint = NF_C_Q_UL4SE_2D_EvalEdge;
      break;
    case NF_C_Q_UL5SE_2D:
      N_AllFunctionals = 29;
      N_EdgeFunctionals = 6;
      N_PointsAll = 29;
      N_PointsEdge = 6;
      Xi = NF_C_Q_UL5SE_2D_Xi;
      Eta = NF_C_Q_UL5SE_2D_Eta;
      T = NF_C_Q_UL5SE_2D_T;
      EvalAll = NF_C_Q_UL5SE_2D_EvalAll;
      EvalJoint = NF_C_Q_UL5SE_2D_EvalEdge;
      break;
    case NF_C_Q_UL6SE_2D:
      N_AllFunctionals = 40;
      N_EdgeFunctionals = 7;
      N_PointsAll = 40;
      N_PointsEdge = 7;
      Xi = NF_C_Q_UL6SE_2D_Xi;
      Eta = NF_C_Q_UL6SE_2D_Eta;
      T = NF_C_Q_UL6SE_2D_T;
      EvalAll = NF_C_Q_UL6SE_2D_EvalAll;
      EvalJoint = NF_C_Q_UL6SE_2D_EvalEdge;
      break;
    case NF_C_Q_UL7SE_2D:
      N_AllFunctionals = 53;
      N_EdgeFunctionals = 8;
      N_PointsAll = 53;
      N_PointsEdge = 8;
      Xi = NF_C_Q_UL7SE_2D_Xi;
      Eta = NF_C_Q_UL7SE_2D_Eta;
      T = NF_C_Q_UL7SE_2D_T;
      EvalAll = NF_C_Q_UL7SE_2D_EvalAll;
      EvalJoint = NF_C_Q_UL7SE_2D_EvalEdge;
      break;
    case NF_C_Q_UL8SE_2D:
      N_AllFunctionals = 68;
      N_EdgeFunctionals = 9;
      N_PointsAll = 68;
      N_PointsEdge = 9;
      Xi = NF_C_Q_UL8SE_2D_Xi;
      Eta = NF_C_Q_UL8SE_2D_Eta;
      T = NF_C_Q_UL8SE_2D_T;
      EvalAll = NF_C_Q_UL8SE_2D_EvalAll;
      EvalJoint = NF_C_Q_UL8SE_2D_EvalEdge;
      break;
    case NF_C_Q_UL9SE_2D:
      N_AllFunctionals = 85;
      N_EdgeFunctionals = 10;
      N_PointsAll = 85;
      N_PointsEdge = 10;
      Xi = NF_C_Q_UL9SE_2D_Xi;
      Eta = NF_C_Q_UL9SE_2D_Eta;
      T = NF_C_Q_UL9SE_2D_T;
      EvalAll = NF_C_Q_UL9SE_2D_EvalAll;
      EvalJoint = NF_C_Q_UL9SE_2D_EvalEdge;
      break;
    case NF_C_Q_M2_2D:
      N_AllFunctionals = 8;
      N_EdgeFunctionals = 3;
      N_PointsAll = 8;
      N_PointsEdge = 3;
      Xi = NF_C_Q_M2_2D_Xi;
      Eta = NF_C_Q_M2_2D_Eta;
      T = NF_C_Q_M2_2D_T;
      EvalAll = NF_C_Q_M2_2D_EvalAll;
      EvalJoint = NF_C_Q_M2_2D_EvalEdge;
      break;
    case NF_C_Q_M3_2D:
      N_AllFunctionals = 12;
      N_EdgeFunctionals = 4;
      N_PointsAll = 12;
      N_PointsEdge = 4;
      Xi = NF_C_Q_M3_2D_Xi;
      Eta = NF_C_Q_M3_2D_Eta;
      T = NF_C_Q_M3_2D_T;
      EvalAll = NF_C_Q_M3_2D_EvalAll;
      EvalJoint = NF_C_Q_M3_2D_EvalEdge;
      break;
    case NF_C_Q_M4_2D:
      N_AllFunctionals = 17;
      N_EdgeFunctionals = 5;
      N_PointsAll = 17;
      N_PointsEdge = 5;
      Xi = NF_C_Q_M4_2D_Xi;
      Eta = NF_C_Q_M4_2D_Eta;
      T = NF_C_Q_M4_2D_T;
      EvalAll = NF_C_Q_M4_2D_EvalAll;
      EvalJoint = NF_C_Q_M4_2D_EvalEdge;
      break;
    case NF_C_Q_M5_2D:
      N_AllFunctionals = 23;
      N_EdgeFunctionals = 6;
      N_PointsAll = 23;
      N_PointsEdge = 6;
      Xi = NF_C_Q_M5_2D_Xi;
      Eta = NF_C_Q_M5_2D_Eta;
      T = NF_C_Q_M5_2D_T;
      EvalAll = NF_C_Q_M5_2D_EvalAll;
      EvalJoint = NF_C_Q_M5_2D_EvalEdge;
      break;
    case NF_C_Q_M6_2D:
      N_AllFunctionals = 30;
      N_EdgeFunctionals = 7;
      N_PointsAll = 30;
      N_PointsEdge = 7;
      Xi = NF_C_Q_M6_2D_Xi;
      Eta = NF_C_Q_M6_2D_Eta;
      T = NF_C_Q_M6_2D_T;
      EvalAll = NF_C_Q_M6_2D_EvalAll;
      EvalJoint = NF_C_Q_M6_2D_EvalEdge;
      break;
    case NF_C_Q_M7_2D:
      N_AllFunctionals = 38;
      N_EdgeFunctionals = 8;
      N_PointsAll = 38;
      N_PointsEdge = 8;
      Xi = NF_C_Q_M7_2D_Xi;
      Eta = NF_C_Q_M7_2D_Eta;
      T = NF_C_Q_M7_2D_T;
      EvalAll = NF_C_Q_M7_2D_EvalAll;
      EvalJoint = NF_C_Q_M7_2D_EvalEdge;
      break;
    case NF_C_Q_M8_2D:
      N_AllFunctionals = 47;
      N_EdgeFunctionals = 9;
      N_PointsAll = 47;
      N_PointsEdge = 9;
      Xi = NF_C_Q_M8_2D_Xi;
      Eta = NF_C_Q_M8_2D_Eta;
      T = NF_C_Q_M8_2D_T;
      EvalAll = NF_C_Q_M8_2D_EvalAll;
      EvalJoint = NF_C_Q_M8_2D_EvalEdge;
      break;
    case NF_C_Q_M9_2D:
      N_AllFunctionals = 57;
      N_EdgeFunctionals = 10;
      N_PointsAll = 57;
      N_PointsEdge = 10;
      Xi = NF_C_Q_M9_2D_Xi;
      Eta = NF_C_Q_M9_2D_Eta;
      T = NF_C_Q_M9_2D_T;
      EvalAll = NF_C_Q_M9_2D_EvalAll;
      EvalJoint = NF_C_Q_M9_2D_EvalEdge;
      break;
    case NF_C_Q_EL1_2D:
      N_AllFunctionals = 5;
      N_EdgeFunctionals = 2;
      N_PointsAll = 8;
      N_PointsEdge = 2;
      Xi = NF_C_Q_EL1_2D_Xi;
      Eta = NF_C_Q_EL1_2D_Eta;
      T = NF_C_Q_EL1_2D_T;
      EvalAll = NF_C_Q_EL1_2D_EvalAll;
      EvalJoint = NF_C_Q_EL1_2D_EvalEdge;
      break;
    case NF_N_Q_RT0_2D:
      N_AllFunctionals = 4;
      N_EdgeFunctionals = 1;
      N_PointsAll = 8;
      N_PointsEdge = 2;
      Xi = NF_N_Q_RT0_2D_Xi;
      Eta = NF_N_Q_RT0_2D_Eta;
      T = NF_N_Q_RT0_2D_T;
      EvalAll = NF_N_Q_RT0_2D_EvalAll;
      EvalJoint = NF_N_Q_RT0_2D_EvalEdge;
      break;
    case NF_N_Q_RT1_2D:
      N_AllFunctionals = 12;
      N_EdgeFunctionals = 2;
      N_PointsAll = 37;
      N_PointsEdge = 3;
      Xi = NF_N_Q_RT1_2D_Xi;
      Eta = NF_N_Q_RT1_2D_Eta;
      T = NF_N_Q_RT1_2D_T;
      EvalAll = NF_N_Q_RT1_2D_EvalAll;
      EvalJoint = NF_N_Q_RT1_2D_EvalEdge;
      break;
    case NF_N_Q_RT2_2D:
      N_AllFunctionals = 24;
      N_EdgeFunctionals = 3;
      N_PointsAll = 32;
      N_PointsEdge = 4;
      Xi = NF_N_Q_RT2_2D_Xi;
      Eta = NF_N_Q_RT2_2D_Eta;
      T = NF_N_Q_RT2_2D_T;
      EvalAll = NF_N_Q_RT2_2D_EvalAll;
      EvalJoint = NF_N_Q_RT2_2D_EvalEdge;
      break;
    case NF_N_Q_RT3_2D:
      N_AllFunctionals = 40;
      N_EdgeFunctionals = 4;
      N_PointsAll = 45;
      N_PointsEdge = 5;
      Xi = NF_N_Q_RT3_2D_Xi;
      Eta = NF_N_Q_RT3_2D_Eta;
      T = NF_N_Q_RT3_2D_q;
      EvalAll = NF_N_Q_RT3_2D_EvalAll;
      EvalJoint = NF_N_Q_RT3_2D_EvalEdge;
      break;
    case NF_N_T_RT0_2D:
      N_AllFunctionals = 3;
      N_EdgeFunctionals = 1;
      N_PointsAll = 6;
      N_PointsEdge = 2;
      Xi = NF_N_T_RT0_2D_Xi;
      Eta = NF_N_T_RT0_2D_Eta;
      T = NF_N_T_RT0_2D_T;
      EvalAll = NF_N_T_RT0_2D_EvalAll;
      EvalJoint = NF_N_T_RT0_2D_EvalEdge;
      break;
    case NF_N_T_RT1_2D:
      N_AllFunctionals = 8;
      N_EdgeFunctionals = 2;
      N_PointsAll = 16;
      N_PointsEdge = 3;
      Xi = NF_N_T_RT1_2D_Xi;
      Eta = NF_N_T_RT1_2D_Eta;
      T = NF_N_T_RT1_2D_T;
      EvalAll = NF_N_T_RT1_2D_EvalAll;
      EvalJoint = NF_N_T_RT1_2D_EvalEdge;
      break;
    case NF_N_T_RT2_2D:
      N_AllFunctionals = 15;
      N_EdgeFunctionals = 3;
      N_PointsAll = 31;
      N_PointsEdge = 4;
      Xi = NF_N_T_RT2_2D_Xi;
      Eta = NF_N_T_RT2_2D_Eta;
      T = NF_N_T_RT2_2D_T;
      EvalAll = NF_N_T_RT2_2D_EvalAll;
      EvalJoint = NF_N_T_RT2_2D_EvalEdge;
      break;
    case NF_N_T_RT3_2D:
      N_AllFunctionals = 24;
      N_EdgeFunctionals = 4;
      N_PointsAll = 34;
      N_PointsEdge = 5;
      Xi = NF_N_T_RT3_2D_Xi;
      Eta = NF_N_T_RT3_2D_Eta;
      T = NF_N_T_RT3_2D_T;
      EvalAll = NF_N_T_RT3_2D_EvalAll;
      EvalJoint = NF_N_T_RT3_2D_EvalEdge;
      break;
    case NF_N_Q_BDM1_2D:
      N_AllFunctionals = 8;
      N_EdgeFunctionals = 2;
      N_PointsAll = 12;
      N_PointsEdge = 3;
      Xi = NF_N_Q_BDM1_2D_Xi;
      Eta = NF_N_Q_BDM1_2D_Eta;
      T = NF_N_Q_BDM1_2D_T;
      EvalAll = NF_N_Q_BDM1_2D_EvalAll;
      EvalJoint = NF_N_Q_BDM1_2D_EvalEdge;
      break;
    case NF_N_Q_BDM2_2D:
      N_AllFunctionals = 14;
      N_EdgeFunctionals = 3;
      N_PointsAll = 32;
      N_PointsEdge = 4;
      Xi = NF_N_Q_BDM2_2D_Xi;
      Eta = NF_N_Q_BDM2_2D_Eta;
      T = NF_N_Q_BDM2_2D_q;
      EvalAll = NF_N_Q_BDM2_2D_EvalAll;
      EvalJoint = NF_N_Q_BDM2_2D_EvalEdge;
      break;
    case NF_N_Q_BDM3_2D:
      N_AllFunctionals = 22;
      N_EdgeFunctionals = 4;
      N_PointsAll = 45;
      N_PointsEdge = 5;
      Xi = NF_N_Q_BDM3_2D_Xi;
      Eta = NF_N_Q_BDM3_2D_Eta;
      T = NF_N_Q_BDM3_2D_q;
      EvalAll = NF_N_Q_BDM3_2D_EvalAll;
      EvalJoint = NF_N_Q_BDM3_2D_EvalEdge;
      break;
    case NF_N_T_BDM1_2D:
      N_AllFunctionals = 6;
      N_EdgeFunctionals = 2;
      N_PointsAll = 9;
      N_PointsEdge = 3;
      Xi = NF_N_T_BDM1_2D_Xi;
      Eta = NF_N_T_BDM1_2D_Eta;
      T = NF_N_T_BDM1_2D_T;
      EvalAll = NF_N_T_BDM1_2D_EvalAll;
      EvalJoint = NF_N_T_BDM1_2D_EvalEdge;
      break;
    case NF_N_T_BDM2_2D:
      N_AllFunctionals = 12;
      N_EdgeFunctionals = 3;
      N_PointsAll = 31;
      N_PointsEdge = 4;
      Xi = NF_N_T_BDM2_2D_Xi;
      Eta = NF_N_T_BDM2_2D_Eta;
      T = NF_N_T_BDM2_2D_T;
      EvalAll = NF_N_T_BDM2_2D_EvalAll;
      EvalJoint = NF_N_T_BDM2_2D_EvalEdge;
      break;
    case NF_N_T_BDM3_2D:
      N_AllFunctionals = 20;
      N_EdgeFunctionals = 4;
      N_PointsAll = 34;
      N_PointsEdge = 5;
      Xi = NF_N_T_BDM3_2D_Xi;
      Eta = NF_N_T_BDM3_2D_Eta;
      T = NF_N_T_BDM3_2D_T;
      EvalAll = NF_N_T_BDM3_2D_EvalAll;
      EvalJoint = NF_N_T_BDM3_2D_EvalEdge;
      break;
    case NF_C_T_P00_3D:
      N_AllFunctionals = NF_C_T_P00_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_C_T_P00_3D_N_FaceFunctionals;
      N_PointsAll = NF_C_T_P00_3D_N_PointsAll;
      N_PointsFace = NF_C_T_P00_3D_N_PointsFace;
      Xi = NF_C_T_P00_3D_Xi;
      Eta = NF_C_T_P00_3D_Eta;
      Zeta = NF_C_T_P00_3D_Zeta;
      XiArray = NF_C_T_P00_3D_XiArray;
      EtaArray = NF_C_T_P00_3D_EtaArray;
      ZetaArray = NF_C_T_P00_3D_ZetaArray;
      T = NF_C_T_P00_3D_T;
      S = NF_C_T_P00_3D_S;
      EvalAll = NF_C_T_P00_3D_EvalAll;
      EvalJoint = NF_C_T_P00_3D_EvalFace;
      break;
    case NF_C_T_P0_3D:
      N_AllFunctionals = NF_C_T_P0_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_C_T_P0_3D_N_FaceFunctionals;
      N_PointsAll = NF_C_T_P0_3D_N_PointsAll;
      N_PointsFace = NF_C_T_P0_3D_N_PointsFace;
      Xi = NF_C_T_P0_3D_Xi;
      Eta = NF_C_T_P0_3D_Eta;
      Zeta = NF_C_T_P0_3D_Zeta;
      XiArray = NF_C_T_P0_3D_XiArray;
      EtaArray = NF_C_T_P0_3D_EtaArray;
      ZetaArray = NF_C_T_P0_3D_ZetaArray;
      T = NF_C_T_P0_3D_T;
      S = NF_C_T_P0_3D_S;
      EvalAll = NF_C_T_P0_3D_EvalAll;
      EvalJoint = NF_C_T_P0_3D_EvalFace;
      break;
    case NF_C_T_P1_3D:
      N_AllFunctionals = NF_C_T_P1_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_C_T_P1_3D_N_FaceFunctionals;
      N_PointsAll = NF_C_T_P1_3D_N_PointsAll;
      N_PointsFace = NF_C_T_P1_3D_N_PointsFace;
      Xi = NF_C_T_P1_3D_Xi;
      Eta = NF_C_T_P1_3D_Eta;
      Zeta = NF_C_T_P1_3D_Zeta;
      XiArray = NF_C_T_P1_3D_XiArray;
      EtaArray = NF_C_T_P1_3D_EtaArray;
      ZetaArray = NF_C_T_P1_3D_ZetaArray;
      T = NF_C_T_P1_3D_T;
      S = NF_C_T_P1_3D_S;
      EvalAll = NF_C_T_P1_3D_EvalAll;
      EvalJoint = NF_C_T_P1_3D_EvalFace;
      break;
    case NF_C_T_P2_3D:
      N_AllFunctionals = NF_C_T_P2_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_C_T_P2_3D_N_FaceFunctionals;
      N_PointsAll = NF_C_T_P2_3D_N_PointsAll;
      N_PointsFace = NF_C_T_P2_3D_N_PointsFace;
      Xi = NF_C_T_P2_3D_Xi;
      Eta = NF_C_T_P2_3D_Eta;
      Zeta = NF_C_T_P2_3D_Zeta;
      XiArray = NF_C_T_P2_3D_XiArray;
      EtaArray = NF_C_T_P2_3D_EtaArray;
      ZetaArray = NF_C_T_P2_3D_ZetaArray;
      T = NF_C_T_P2_3D_T;
      S = NF_C_T_P2_3D_S;
      EvalAll = NF_C_T_P2_3D_EvalAll;
      EvalJoint = NF_C_T_P2_3D_EvalFace;
      break;
    case NF_C_T_P3_3D:
      N_AllFunctionals = NF_C_T_P3_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_C_T_P3_3D_N_FaceFunctionals;
      N_PointsAll = NF_C_T_P3_3D_N_PointsAll;
      N_PointsFace = NF_C_T_P3_3D_N_PointsFace;
      Xi = NF_C_T_P3_3D_Xi;
      Eta = NF_C_T_P3_3D_Eta;
      Zeta = NF_C_T_P3_3D_Zeta;
      XiArray = NF_C_T_P3_3D_XiArray;
      EtaArray = NF_C_T_P3_3D_EtaArray;
      ZetaArray = NF_C_T_P3_3D_ZetaArray;
      T = NF_C_T_P3_3D_T;
      S = NF_C_T_P3_3D_S;
      EvalAll = NF_C_T_P3_3D_EvalAll;
      EvalJoint = NF_C_T_P3_3D_EvalFace;
      break;
    case NF_N_T_P1_3D:
      N_AllFunctionals = NF_N_T_P1_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_N_T_P1_3D_N_FaceFunctionals;
      N_PointsAll = NF_N_T_P1_3D_N_T_PointsAll;
      N_PointsFace = NF_N_T_P1_3D_N_T_PointsFace;
      Xi = NF_N_T_P1_3D_Xi;
      Eta = NF_N_T_P1_3D_Eta;
      Zeta = NF_N_T_P1_3D_Zeta;
      XiArray = NF_N_T_P1_3D_XiArray;
      EtaArray = NF_N_T_P1_3D_EtaArray;
      ZetaArray = NF_N_T_P1_3D_ZetaArray;
      T = NF_N_T_P1_3D_T;
      S = NF_N_T_P1_3D_S;
      EvalAll = NF_N_T_P1_3D_EvalAll;
      EvalJoint = NF_N_T_P1_3D_EvalFace;
      break;
    case NF_C_H_Q00_3D:
      N_AllFunctionals = NF_C_H_Q00_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_C_H_Q00_3D_N_FaceFunctionals;
      N_PointsAll = NF_C_H_Q00_3D_N_PointsAll;
      N_PointsFace = NF_C_H_Q00_3D_N_PointsFace;
      Xi = NF_C_H_Q00_3D_Xi;
      Eta = NF_C_H_Q00_3D_Eta;
      Zeta = NF_C_H_Q00_3D_Zeta;
      XiArray = NF_C_H_Q00_3D_XiArray;
      EtaArray = NF_C_H_Q00_3D_EtaArray;
      ZetaArray = NF_C_H_Q00_3D_ZetaArray;
      T = NF_C_H_Q00_3D_T;
      S = NF_C_H_Q00_3D_S;
      EvalAll = NF_C_H_Q00_3D_EvalAll;
      EvalJoint = NF_C_H_Q00_3D_EvalFace;
      break;
    case NF_C_H_Q0_3D:
      N_AllFunctionals = NF_C_H_Q0_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_C_H_Q0_3D_N_FaceFunctionals;
      N_PointsAll = NF_C_H_Q0_3D_N_PointsAll;
      N_PointsFace = NF_C_H_Q0_3D_N_PointsFace;
      Xi = NF_C_H_Q0_3D_Xi;
      Eta = NF_C_H_Q0_3D_Eta;
      Zeta = NF_C_H_Q0_3D_Zeta;
      XiArray = NF_C_H_Q0_3D_XiArray;
      EtaArray = NF_C_H_Q0_3D_EtaArray;
      ZetaArray = NF_C_H_Q0_3D_ZetaArray;
      T = NF_C_H_Q0_3D_T;
      S = NF_C_H_Q0_3D_S;
      EvalAll = NF_C_H_Q0_3D_EvalAll;
      EvalJoint = NF_C_H_Q0_3D_EvalFace;
      break;
    case NF_C_H_Q1_3D:
      N_AllFunctionals = NF_C_H_Q1_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_C_H_Q1_3D_N_FaceFunctionals;
      N_PointsAll = NF_C_H_Q1_3D_N_PointsAll;
      N_PointsFace = NF_C_H_Q1_3D_N_PointsFace;
      Xi = NF_C_H_Q1_3D_Xi;
      Eta = NF_C_H_Q1_3D_Eta;
      Zeta = NF_C_H_Q1_3D_Zeta;
      XiArray = NF_C_H_Q1_3D_XiArray;
      EtaArray = NF_C_H_Q1_3D_EtaArray;
      ZetaArray = NF_C_H_Q1_3D_ZetaArray;
      T = NF_C_H_Q1_3D_T;
      S = NF_C_H_Q1_3D_S;
      EvalAll = NF_C_H_Q1_3D_EvalAll;
      EvalJoint = NF_C_H_Q1_3D_EvalFace;
      break;
    case NF_C_H_Q2_3D:
      N_AllFunctionals = NF_C_H_Q2_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_C_H_Q2_3D_N_FaceFunctionals;
      N_PointsAll = NF_C_H_Q2_3D_N_PointsAll;
      N_PointsFace = NF_C_H_Q2_3D_N_PointsFace;
      Xi = NF_C_H_Q2_3D_Xi;
      Eta = NF_C_H_Q2_3D_Eta;
      Zeta = NF_C_H_Q2_3D_Zeta;
      XiArray = NF_C_H_Q2_3D_XiArray;
      EtaArray = NF_C_H_Q2_3D_EtaArray;
      ZetaArray = NF_C_H_Q2_3D_ZetaArray;
      T = NF_C_H_Q2_3D_T;
      S = NF_C_H_Q2_3D_S;
      EvalAll = NF_C_H_Q2_3D_EvalAll;
      EvalJoint = NF_C_H_Q2_3D_EvalFace;
      break;
    case NF_C_H_Q3_3D:
      N_AllFunctionals = NF_C_H_Q3_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_C_H_Q3_3D_N_FaceFunctionals;
      N_PointsAll = NF_C_H_Q3_3D_N_PointsAll;
      N_PointsFace = NF_C_H_Q3_3D_N_PointsFace;
      Xi = NF_C_H_Q3_3D_Xi;
      Eta = NF_C_H_Q3_3D_Eta;
      Zeta = NF_C_H_Q3_3D_Zeta;
      XiArray = NF_C_H_Q3_3D_XiArray;
      EtaArray = NF_C_H_Q3_3D_EtaArray;
      ZetaArray = NF_C_H_Q3_3D_ZetaArray;
      T = NF_C_H_Q3_3D_T;
      S = NF_C_H_Q3_3D_S;
      EvalAll = NF_C_H_Q3_3D_EvalAll;
      EvalJoint = NF_C_H_Q3_3D_EvalFace;
      break;
    case NF_C_H_Q4_3D:
      N_AllFunctionals = NF_C_H_Q4_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_C_H_Q4_3D_N_FaceFunctionals;
      N_PointsAll = NF_C_H_Q4_3D_N_PointsAll;
      N_PointsFace = NF_C_H_Q4_3D_N_PointsFace;
      Xi = NF_C_H_Q4_3D_Xi;
      Eta = NF_C_H_Q4_3D_Eta;
      Zeta = NF_C_H_Q4_3D_Zeta;
      XiArray = NF_C_H_Q4_3D_XiArray;
      EtaArray = NF_C_H_Q4_3D_EtaArray;
      ZetaArray = NF_C_H_Q4_3D_ZetaArray;
      T = NF_C_H_Q4_3D_T;
      S = NF_C_H_Q4_3D_S;
      EvalAll = NF_C_H_Q4_3D_EvalAll;
      EvalJoint = NF_C_H_Q4_3D_EvalFace;
      break;
    case NF_N_H_Q1_3D:
      N_AllFunctionals = NF_N_H_Q1_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_N_H_Q1_3D_N_FaceFunctionals;
      N_PointsAll = NF_N_H_Q1_3D_N_PointsAll;
      N_PointsFace = NF_N_H_Q1_3D_N_PointsFace;
      Xi = NF_N_H_Q1_3D_Xi;
      Eta = NF_N_H_Q1_3D_Eta;
      Zeta = NF_N_H_Q1_3D_Zeta;
      XiArray = NF_N_H_Q1_3D_XiArray;
      EtaArray = NF_N_H_Q1_3D_EtaArray;
      ZetaArray = NF_N_H_Q1_3D_ZetaArray;
      T = NF_N_H_Q1_3D_T;
      S = NF_N_H_Q1_3D_S;
      EvalAll = NF_N_H_Q1_3D_EvalAll;
      EvalJoint = NF_N_H_Q1_3D_EvalFace;
      break;
    case NF_D_H_P1_3D:
      N_AllFunctionals = NF_D_H_P1_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_D_H_P1_3D_N_FaceFunctionals;
      N_PointsAll = NF_D_H_P1_3D_N_PointsAll;
      N_PointsFace = NF_D_H_P1_3D_N_PointsFace;
      Xi = NF_D_H_P1_3D_Xi;
      Eta = NF_D_H_P1_3D_Eta;
      Zeta = NF_D_H_P1_3D_Zeta;
      XiArray = NF_D_H_P1_3D_XiArray;
      EtaArray = NF_D_H_P1_3D_EtaArray;
      ZetaArray = NF_D_H_P1_3D_ZetaArray;
      T = NF_D_H_P1_3D_T;
      S = NF_D_H_P1_3D_S;
      EvalAll = NF_D_H_P1_3D_EvalAll;
      EvalJoint = NF_D_H_P1_3D_EvalFace;
      break;
    case NF_D_H_P2_3D:
      N_AllFunctionals = NF_D_H_P2_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_D_H_P2_3D_N_FaceFunctionals;
      N_PointsAll = NF_D_H_P2_3D_N_PointsAll;
      N_PointsFace = NF_D_H_P2_3D_N_PointsFace;
      Xi = NF_D_H_P2_3D_Xi;
      Eta = NF_D_H_P2_3D_Eta;
      Zeta = NF_D_H_P2_3D_Zeta;
      XiArray = NF_D_H_P2_3D_XiArray;
      EtaArray = NF_D_H_P2_3D_EtaArray;
      ZetaArray = NF_D_H_P2_3D_ZetaArray;
      T = NF_D_H_P2_3D_T;
      S = NF_D_H_P2_3D_S;
      EvalAll = NF_D_H_P2_3D_EvalAll;
      EvalJoint = NF_D_H_P2_3D_EvalFace;
      break;
    case NF_D_H_P3_3D:
      N_AllFunctionals = NF_D_H_P3_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_D_H_P3_3D_N_FaceFunctionals;
      N_PointsAll = NF_D_H_P3_3D_N_PointsAll;
      N_PointsFace = NF_D_H_P3_3D_N_PointsFace;
      Xi = NF_D_H_P3_3D_Array2;
      Eta = NF_D_H_P3_3D_Array3;
      Zeta = NF_D_H_P3_3D_Array4;
      XiArray = NF_D_H_P3_3D_XiArray;
      EtaArray = NF_D_H_P3_3D_EtaArray;
      ZetaArray = NF_D_H_P3_3D_ZetaArray;
      T = NF_D_H_P3_3D_T;
      S = NF_D_H_P3_3D_S;
      EvalAll = NF_D_H_P3_3D_EvalAll;
      EvalJoint = NF_D_H_P3_3D_EvalFace;
      break;
    case NF_D_H_Q1_3D:
      N_AllFunctionals = NF_D_H_Q1_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_D_H_Q1_3D_N_FaceFunctionals;
      N_PointsAll = NF_D_H_Q1_3D_N_PointsAll;
      N_PointsFace = NF_D_H_Q1_3D_N_PointsFace;
      Xi = NF_D_H_Q1_3D_Xi;
      Eta = NF_D_H_Q1_3D_Eta;
      Zeta = NF_D_H_Q1_3D_Zeta;
      XiArray = NF_D_H_Q1_3D_XiArray;
      EtaArray = NF_D_H_Q1_3D_EtaArray;
      ZetaArray = NF_D_H_Q1_3D_ZetaArray;
      T = NF_D_H_Q1_3D_T;
      S = NF_D_H_Q1_3D_S;
      EvalAll = NF_D_H_Q1_3D_EvalAll;
      EvalJoint = NF_D_H_Q1_3D_EvalFace;
      break;
    case NF_D_H_Q2_3D:
      N_AllFunctionals = NF_D_H_Q2_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_D_H_Q2_3D_N_FaceFunctionals;
      N_PointsAll = NF_D_H_Q2_3D_N_PointsAll;
      N_PointsFace = NF_D_H_Q2_3D_N_PointsFace;
      Xi = NF_D_H_Q2_3D_Xi;
      Eta = NF_D_H_Q2_3D_Eta;
      Zeta = NF_D_H_Q2_3D_Zeta;
      XiArray = NF_D_H_Q2_3D_XiArray;
      EtaArray = NF_D_H_Q2_3D_EtaArray;
      ZetaArray = NF_D_H_Q2_3D_ZetaArray;
      T = NF_D_H_Q2_3D_T;
      S = NF_D_H_Q2_3D_S;
      EvalAll = NF_D_H_Q2_3D_EvalAll;
      EvalJoint = NF_D_H_Q2_3D_EvalFace;
      break;
    case NF_B_H_IB2_3D:
      N_AllFunctionals = NF_B_H_IB2_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_B_H_IB2_3D_N_FaceFunctionals;
      N_PointsAll = NF_B_H_IB2_3D_N_PointsAll;
      N_PointsFace = NF_B_H_IB2_3D_N_PointsFace;
      Xi = NF_B_H_IB2_3D_Xi;
      Eta = NF_B_H_IB2_3D_Eta;
      Zeta = NF_B_H_IB2_3D_Zeta;
      XiArray = NF_B_H_IB2_3D_XiArray;
      EtaArray = NF_B_H_IB2_3D_EtaArray;
      ZetaArray = NF_B_H_IB2_3D_ZetaArray;
      T = NF_B_H_IB2_3D_T;
      S = NF_B_H_IB2_3D_S;
      EvalAll = NF_B_H_IB2_3D_EvalAll;
      EvalJoint = NF_B_H_IB2_3D_EvalFace;
      break;
    case NF_S_H_Q2_3D:
      N_AllFunctionals = NF_S_H_Q2_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_S_H_Q2_3D_N_FaceFunctionals;
      N_PointsAll = NF_S_H_Q2_3D_N_PointsAll;
      N_PointsFace = NF_S_H_Q2_3D_N_PointsFace;
      Xi = NF_S_H_Q2_3D_Xi;
      Eta = NF_S_H_Q2_3D_Eta;
      Zeta = NF_S_H_Q2_3D_Zeta;
      XiArray = NF_S_H_Q2_3D_XiArray;
      EtaArray = NF_S_H_Q2_3D_EtaArray;
      ZetaArray = NF_S_H_Q2_3D_ZetaArray;
      T = NF_S_H_Q2_3D_T;
      S = NF_S_H_Q2_3D_S;
      EvalAll = NF_S_H_Q2_3D_EvalAll;
      EvalJoint = NF_S_H_Q2_3D_EvalFace;
      break;
    case NF_N_T_P2_3D:
      N_AllFunctionals = NF_N_T_P2_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_N_T_P2_3D_N_FaceFunctionals;
      N_PointsAll = NF_N_T_P2_3D_N_T_PointsAll;
      N_PointsFace = NF_N_T_P2_3D_N_T_PointsFace;
      Xi = NF_N_T_P2_3D_Xi;
      Eta = NF_N_T_P2_3D_Eta;
      Zeta = NF_N_T_P2_3D_Zeta;
      XiArray = NF_N_T_P2_3D_XiArray;
      EtaArray = NF_N_T_P2_3D_EtaArray;
      ZetaArray = NF_N_T_P2_3D_ZetaArray;
      T = NF_N_T_P2_3D_T;
      S = NF_N_T_P2_3D_S;
      EvalAll = NF_N_T_P2_3D_EvalAll;
      EvalJoint = NF_N_T_P2_3D_EvalFace;
      break;
    case NF_N_T_P3_3D:
      ErrThrow("Nodal functional of type NF_N_T_P3_3D not implemented");
      break;
    case NF_N_T_P4_3D:
      ErrThrow("Nodal functional of type NF_N_T_P4_3D not implemented");
      break;
    case NF_N_T_P5_3D:
      ErrThrow("Nodal functional of type NF_N_T_P5_3D not implemented");
      break;
    case NF_N_H_Q2_3D:
      N_AllFunctionals = NF_N_H_Q2_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_N_H_Q2_3D_N_FaceFunctionals;
      N_PointsAll = NF_N_H_Q2_3D_N_PointsAll;
      N_PointsFace = NF_N_H_Q2_3D_N_PointsFace;
      Xi = NF_N_H_Q2_3D_Xi;
      Eta = NF_N_H_Q2_3D_Eta;
      Zeta = NF_N_H_Q2_3D_Zeta;
      XiArray = NF_N_H_Q2_3D_XiArray;
      EtaArray = NF_N_H_Q2_3D_EtaArray;
      ZetaArray = NF_N_H_Q2_3D_ZetaArray;
      T = NF_N_H_Q2_3D_T;
      S = NF_N_H_Q2_3D_S;
      EvalAll = NF_N_H_Q2_3D_EvalAll;
      EvalJoint = NF_N_H_Q2_3D_EvalFace;
      break;
    case NF_N_H_Q3_3D:
      N_AllFunctionals = NF_N_H_Q3_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_N_H_Q3_3D_N_FaceFunctionals;
      N_PointsAll = NF_N_H_Q3_3D_N_PointsAll;
      N_PointsFace = NF_N_H_Q3_3D_N_PointsFace;
      Xi = NF_N_H_Q3_3D_Xi;
      Eta = NF_N_H_Q3_3D_Eta;
      Zeta = NF_N_H_Q3_3D_Zeta;
      XiArray = NF_N_H_Q3_3D_XiArray;
      EtaArray = NF_N_H_Q3_3D_EtaArray;
      ZetaArray = NF_N_H_Q3_3D_ZetaArray;
      T = NF_N_H_Q3_3D_T;
      S = NF_N_H_Q3_3D_S;
      EvalAll = NF_N_H_Q3_3D_EvalAll;
      EvalJoint = NF_N_H_Q3_3D_EvalFace;
      break;
    case NF_N_H_Q4_3D:
      N_AllFunctionals = NF_N_H_Q4_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_N_H_Q4_3D_N_FaceFunctionals;
      N_PointsAll = NF_N_H_Q4_3D_N_PointsAll;
      N_PointsFace = NF_N_H_Q4_3D_N_PointsFace;
      Xi = NF_N_H_Q4_3D_Xi;
      Eta = NF_N_H_Q4_3D_Eta;
      Zeta = NF_N_H_Q4_3D_Zeta;
      XiArray = NF_N_H_Q4_3D_XiArray;
      EtaArray = NF_N_H_Q4_3D_EtaArray;
      ZetaArray = NF_N_H_Q4_3D_ZetaArray;
      T = NF_N_H_Q4_3D_T;
      S = NF_N_H_Q4_3D_S;
      EvalAll = NF_N_H_Q4_3D_EvalAll;
      EvalJoint = NF_N_H_Q4_3D_EvalFace;
      break;
    case NF_N_H_Q5_3D:
      ErrThrow("Nodal functional of type NF_N_H_Q5_3D not implemented");
      break;
    case NF_C_T_B2_3D:
      N_AllFunctionals = NF_C_T_B2_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_C_T_B2_3D_N_FaceFunctionals;
      N_PointsAll = NF_C_T_B2_3D_N_PointsAll;
      N_PointsFace = NF_C_T_B2_3D_N_PointsFace;
      Xi = NF_C_T_B2_3D_Xi;
      Eta = NF_C_T_B2_3D_Eta;
      Zeta = NF_C_T_B2_3D_Zeta;
      XiArray = NF_C_T_B2_3D_XiArray;
      EtaArray = NF_C_T_B2_3D_EtaArray;
      ZetaArray = NF_C_T_B2_3D_ZetaArray;
      T = NF_C_T_B2_3D_T;
      S = NF_C_T_B2_3D_S;
      EvalAll = NF_C_T_B2_3D_EvalAll;
      EvalJoint = NF_C_T_B2_3D_EvalFace;
      break;
    case NF_D_T_P1_3D:
      N_AllFunctionals = NF_D_T_P1_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_D_T_P1_3D_N_FaceFunctionals;
      N_PointsAll = NF_D_T_P1_3D_N_PointsAll;
      N_PointsFace = NF_D_T_P1_3D_N_PointsFace;
      Xi = NF_D_T_P1_3D_Xi;
      Eta = NF_D_T_P1_3D_Eta;
      Zeta = NF_D_T_P1_3D_Zeta;
      XiArray = NF_D_T_P1_3D_XiArray;
      EtaArray = NF_D_T_P1_3D_EtaArray;
      ZetaArray = NF_D_T_P1_3D_ZetaArray;
      T = NF_D_T_P1_3D_T;
      S = NF_D_T_P1_3D_S;
      EvalAll = NF_D_T_P1_3D_EvalAll;
      EvalJoint = NF_D_T_P1_3D_EvalFace;
      break;
    case NF_D_T_P2_3D:
      N_AllFunctionals = NF_D_T_P2_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_D_T_P2_3D_N_FaceFunctionals;
      N_PointsAll = NF_D_T_P2_3D_N_PointsAll;
      N_PointsFace = NF_D_T_P2_3D_N_PointsFace;
      Xi = NF_D_T_P2_3D_Xi;
      Eta = NF_D_T_P2_3D_Eta;
      Zeta = NF_D_T_P2_3D_Zeta;
      XiArray = NF_D_T_P2_3D_XiArray;
      EtaArray = NF_D_T_P2_3D_EtaArray;
      ZetaArray = NF_D_T_P2_3D_ZetaArray;
      T = NF_D_T_P2_3D_T;
      S = NF_D_T_P2_3D_S;
      EvalAll = NF_D_T_P2_3D_EvalAll;
      EvalJoint = NF_D_T_P2_3D_EvalFace;
      break;
    case NF_D_T_P3_3D:
      N_AllFunctionals = NF_D_T_P3_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_D_T_P3_3D_N_FaceFunctionals;
      N_PointsAll = NF_D_T_P3_3D_N_PointsAll;
      N_PointsFace = NF_D_T_P3_3D_N_PointsFace;
      Xi = NF_D_T_P3_3D_Xi;
      Eta = NF_D_T_P3_3D_Eta;
      Zeta = NF_D_T_P3_3D_Zeta;
      XiArray = NF_D_T_P3_3D_XiArray;
      EtaArray = NF_D_T_P3_3D_EtaArray;
      ZetaArray = NF_D_T_P3_3D_ZetaArray;
      T = NF_D_T_P3_3D_T;
      S = NF_D_T_P3_3D_S;
      EvalAll = NF_D_T_P3_3D_EvalAll;
      EvalJoint = NF_D_T_P3_3D_EvalFace;
      break;
    case NF_C_H_UL1_3D:
      N_AllFunctionals = NF_C_H_UL1_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_C_H_UL1_3D_N_FaceFunctionals;
      N_PointsAll = NF_C_H_UL1_3D_N_PointsAll;
      N_PointsFace = NF_C_H_UL1_3D_N_PointsFace;
      Xi = NF_C_H_UL1_3D_Xi;
      Eta = NF_C_H_UL1_3D_Eta;
      Zeta = NF_C_H_UL1_3D_Zeta;
      XiArray = NF_C_H_UL1_3D_XiArray;
      EtaArray = NF_C_H_UL1_3D_EtaArray;
      ZetaArray = NF_C_H_UL1_3D_ZetaArray;
      T = NF_C_H_UL1_3D_T;
      S = NF_C_H_UL1_3D_S;
      EvalAll = NF_C_H_UL1_3D_EvalAll;
      EvalJoint = NF_C_H_UL1_3D_EvalFace;
      break;
    case NF_C_H_UL2_3D:
      N_AllFunctionals = NF_C_H_UL2_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_C_H_UL2_3D_N_FaceFunctionals;
      N_PointsAll = NF_C_H_UL2_3D_N_PointsAll;
      N_PointsFace = NF_C_H_UL2_3D_N_PointsFace;
      Xi = NF_C_H_UL2_3D_Xi;
      Eta = NF_C_H_UL2_3D_Eta;
      Zeta = NF_C_H_UL2_3D_Zeta;
      XiArray = NF_C_H_UL2_3D_XiArray;
      EtaArray = NF_C_H_UL2_3D_EtaArray;
      ZetaArray = NF_C_H_UL2_3D_ZetaArray;
      T = NF_C_H_UL2_3D_T;
      S = NF_C_H_UL2_3D_S;
      EvalAll = NF_C_H_UL2_3D_EvalAll;
      EvalJoint = NF_C_H_UL2_3D_EvalFace;
      break;
    case NF_C_H_UL3_3D:
      N_AllFunctionals = NF_C_H_UL3_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_C_H_UL3_3D_N_FaceFunctionals;
      N_PointsAll = NF_C_H_UL3_3D_N_PointsAll;
      N_PointsFace = NF_C_H_UL3_3D_N_PointsFace;
      Xi = NF_C_H_UL3_3D_Xi;
      Eta = NF_C_H_UL3_3D_Eta;
      Zeta = NF_C_H_UL3_3D_Zeta;
      XiArray = NF_C_H_UL3_3D_XiArray;
      EtaArray = NF_C_H_UL3_3D_EtaArray;
      ZetaArray = NF_C_H_UL3_3D_ZetaArray;
      T = NF_C_H_UL3_3D_T;
      S = NF_C_H_UL3_3D_S;
      EvalAll = NF_C_H_UL3_3D_EvalAll;
      EvalJoint = NF_C_H_UL3_3D_EvalFace;
      break;
    case NF_C_H_UL4_3D:
      ErrThrow("Nodal functional of type NF_C_H_UL4_3D not implemented");
      break;
    case NF_C_H_UL5_3D:
      ErrThrow("Nodal functional of type NF_C_H_UL5_3D not implemented");
      break;
    case NF_C_T_UL1_3D:
      ErrThrow("Nodal functional of type NF_C_T_UL1_3D not implemented");
      break;
    case NF_C_T_UL2_3D:
      ErrThrow("Nodal functional of type NF_C_T_UL2_3D not implemented");
      break;
    case NF_C_T_UL3_3D:
      ErrThrow("Nodal functional of type NF_C_T_UL3_3D not implemented");
      break;
    case NF_C_T_UL4_3D:
      ErrThrow("Nodal functional of type NF_C_T_UL4_3D not implemented");
      break;
    case NF_C_T_UL5_3D:
      ErrThrow("Nodal functional of type NF_C_T_UL5_3D not implemented");
      break;
    case NF_N_T_RT0_3D:
      N_AllFunctionals = NF_N_T_RT0_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_N_T_RT0_3D_N_FaceFunctionals;
      N_PointsAll = NF_N_T_RT0_3D_N_PointsAll;
      N_PointsFace = NF_N_T_RT0_3D_N_PointsFace;
      Xi = NF_N_T_RT0_3D_Xi;
      Eta = NF_N_T_RT0_3D_Eta;
      Zeta = NF_N_T_RT0_3D_Zeta;
      XiArray = NF_N_T_RT0_3D_XiArray;
      EtaArray = NF_N_T_RT0_3D_EtaArray;
      ZetaArray = NF_N_T_RT0_3D_ZetaArray;
      T = NF_N_T_RT0_3D_T;
      S = NF_N_T_RT0_3D_S;
      EvalAll = NF_N_T_RT0_3D_EvalAll;
      EvalJoint = NF_N_T_RT0_3D_EvalFace;
      break;
    case NF_N_T_RT1_3D:
      N_AllFunctionals = NF_N_T_RT1_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_N_T_RT1_3D_N_FaceFunctionals;
      N_PointsAll = NF_N_T_RT1_3D_N_PointsAll;
      N_PointsFace = NF_N_T_RT1_3D_N_PointsFace;
      Xi = NF_N_T_RT1_3D_Xi;
      Eta = NF_N_T_RT1_3D_Eta;
      Zeta = NF_N_T_RT1_3D_Zeta;
      XiArray = NF_N_T_RT1_3D_XiArray;
      EtaArray = NF_N_T_RT1_3D_EtaArray;
      ZetaArray = NF_N_T_RT1_3D_ZetaArray;
      T = NF_N_T_RT1_3D_T;
      S = NF_N_T_RT1_3D_S;
      EvalAll = NF_N_T_RT1_3D_EvalAll;
      EvalJoint = NF_N_T_RT1_3D_EvalFace;
      break;
    case NF_N_T_RT2_3D:
      N_AllFunctionals = NF_N_T_RT2_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_N_T_RT2_3D_N_FaceFunctionals;
      N_PointsAll = NF_N_T_RT2_3D_N_PointsAll;
      N_PointsFace = NF_N_T_RT2_3D_N_PointsFace;
      Xi = NF_N_T_RT2_3D_Xi;
      Eta = NF_N_T_RT2_3D_Eta;
      Zeta = NF_N_T_RT2_3D_Zeta;
      XiArray = NF_N_T_RT2_3D_XiArray;
      EtaArray = NF_N_T_RT2_3D_EtaArray;
      ZetaArray = NF_N_T_RT2_3D_ZetaArray;
      T = NF_N_T_RT2_3D_T;
      S = NF_N_T_RT2_3D_S;
      EvalAll = NF_N_T_RT2_3D_EvalAll;
      EvalJoint = NF_N_T_RT2_3D_EvalFace;
      break;
    case NF_N_T_RT3_3D:
      N_AllFunctionals = NF_N_T_RT3_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_N_T_RT3_3D_N_FaceFunctionals;
      N_PointsAll = NF_N_T_RT3_3D_N_PointsAll;
      N_PointsFace = NF_N_T_RT3_3D_N_PointsFace;
      Xi = NF_N_T_RT3_3D_Xi;
      Eta = NF_N_T_RT3_3D_Eta;
      Zeta = NF_N_T_RT3_3D_Zeta;
      XiArray = NF_N_T_RT3_3D_XiArray;
      EtaArray = NF_N_T_RT3_3D_EtaArray;
      ZetaArray = NF_N_T_RT3_3D_ZetaArray;
      T = NF_N_T_RT3_3D_T;
      S = NF_N_T_RT3_3D_S;
      EvalAll = NF_N_T_RT3_3D_EvalAll;
      EvalJoint = NF_N_T_RT3_3D_EvalFace;
      break;
    case NF_N_H_RT0_3D:
      N_AllFunctionals = NF_N_H_RT0_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_N_H_RT0_3D_N_FaceFunctionals;
      N_PointsAll = NF_N_H_RT0_3D_N_PointsAll;
      N_PointsFace = NF_N_H_RT0_3D_N_PointsFace;
      Xi = NF_N_H_RT0_3D_Xi;
      Eta = NF_N_H_RT0_3D_Eta;
      Zeta = NF_N_H_RT0_3D_Zeta;
      XiArray = NF_N_H_RT0_3D_XiArray;
      EtaArray = NF_N_H_RT0_3D_EtaArray;
      ZetaArray = NF_N_H_RT0_3D_ZetaArray;
      T = NF_N_H_RT0_3D_T;
      S = NF_N_H_RT0_3D_S;
      EvalAll = NF_N_H_RT0_3D_EvalAll;
      EvalJoint = NF_N_H_RT0_3D_EvalFace;
      break;
    case NF_N_H_RT1_3D:
      N_AllFunctionals = NF_N_H_RT1_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_N_H_RT1_3D_N_FaceFunctionals;
      N_PointsAll = NF_N_H_RT1_3D_N_PointsAll;
      N_PointsFace = NF_N_H_RT1_3D_N_PointsFace;
      Xi = NF_N_H_RT1_3D_Xi;
      Eta = NF_N_H_RT1_3D_Eta;
      Zeta = NF_N_H_RT1_3D_Zeta;
      XiArray = NF_N_H_RT1_3D_XiArray;
      EtaArray = NF_N_H_RT1_3D_EtaArray;
      ZetaArray = NF_N_H_RT1_3D_ZetaArray;
      T = NF_N_H_RT1_3D_T;
      S = NF_N_H_RT1_3D_S;
      EvalAll = NF_N_H_RT1_3D_EvalAll;
      EvalJoint = NF_N_H_RT1_3D_EvalFace;
      break;
    case NF_N_H_RT2_3D:
      N_AllFunctionals = NF_N_H_RT2_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_N_H_RT2_3D_N_FaceFunctionals;
      N_PointsAll = NF_N_H_RT2_3D_N_PointsAll;
      N_PointsFace = NF_N_H_RT2_3D_N_PointsFace;
      Xi = NF_N_H_RT2_3D_Xi;
      Eta = NF_N_H_RT2_3D_Eta;
      Zeta = NF_N_H_RT2_3D_Zeta;
      XiArray = NF_N_H_RT2_3D_XiArray;
      EtaArray = NF_N_H_RT2_3D_EtaArray;
      ZetaArray = NF_N_H_RT2_3D_ZetaArray;
      T = NF_N_H_RT2_3D_T;
      S = NF_N_H_RT2_3D_S;
      EvalAll = NF_N_H_RT2_3D_EvalAll;
      EvalJoint = NF_N_H_RT2_3D_EvalFace;
      break;
    case NF_N_T_BDDF1_3D:
      N_AllFunctionals = NF_N_T_BDDF1_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_N_T_BDDF1_3D_N_FaceFunctionals;
      N_PointsAll = NF_N_T_BDDF1_3D_N_PointsAll;
      N_PointsFace = NF_N_T_BDDF1_3D_N_PointsFace;
      Xi = NF_N_T_BDDF1_3D_Xi;
      Eta = NF_N_T_BDDF1_3D_Eta;
      Zeta = NF_N_T_BDDF1_3D_Zeta;
      XiArray = NF_N_T_BDDF1_3D_XiArray;
      EtaArray = NF_N_T_BDDF1_3D_EtaArray;
      ZetaArray = NF_N_T_BDDF1_3D_ZetaArray;
      T = NF_N_T_BDDF1_3D_T;
      S = NF_N_T_BDDF1_3D_S;
      EvalAll = NF_N_T_BDDF1_3D_EvalAll;
      EvalJoint = NF_N_T_BDDF1_3D_EvalFace;
      break;
    case NF_N_T_BDDF2_3D:
      N_AllFunctionals = NF_N_T_BDDF2_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_N_T_BDDF2_3D_N_FaceFunctionals;
      N_PointsAll = NF_N_T_BDDF2_3D_N_PointsAll;
      N_PointsFace = NF_N_T_BDDF2_3D_N_PointsFace;
      Xi = NF_N_T_BDDF2_3D_Xi;
      Eta = NF_N_T_BDDF2_3D_Eta;
      Zeta = NF_N_T_BDDF2_3D_Zeta;
      XiArray = NF_N_T_BDDF2_3D_XiArray;
      EtaArray = NF_N_T_BDDF2_3D_EtaArray;
      ZetaArray = NF_N_T_BDDF2_3D_ZetaArray;
      T = NF_N_T_BDDF2_3D_T;
      S = NF_N_T_BDDF2_3D_S;
      EvalAll = NF_N_T_BDDF2_3D_EvalAll;
      EvalJoint = NF_N_T_BDDF2_3D_EvalFace;
      break;
    case NF_N_T_BDDF3_3D:
      N_AllFunctionals = NF_N_T_BDDF3_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_N_T_BDDF3_3D_N_FaceFunctionals;
      N_PointsAll = NF_N_T_BDDF3_3D_N_PointsAll;
      N_PointsFace = NF_N_T_BDDF3_3D_N_PointsFace;
      Xi = NF_N_T_BDDF3_3D_Xi;
      Eta = NF_N_T_BDDF3_3D_Eta;
      Zeta = NF_N_T_BDDF3_3D_Zeta;
      XiArray = NF_N_T_BDDF3_3D_XiArray;
      EtaArray = NF_N_T_BDDF3_3D_EtaArray;
      ZetaArray = NF_N_T_BDDF3_3D_ZetaArray;
      T = NF_N_T_BDDF3_3D_T;
      S = NF_N_T_BDDF3_3D_S;
      EvalAll = NF_N_T_BDDF3_3D_EvalAll;
      EvalJoint = NF_N_T_BDDF3_3D_EvalFace;
      break;
    case NF_N_H_BDDF1_3D:
      N_AllFunctionals = NF_N_H_BDDF1_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_N_H_BDDF1_3D_N_FaceFunctionals;
      N_PointsAll = NF_N_H_BDDF1_3D_N_PointsAll;
      N_PointsFace = NF_N_H_BDDF1_3D_N_PointsFace;
      Xi = NF_N_H_BDDF1_3D_Xi;
      Eta = NF_N_H_BDDF1_3D_Eta;
      Zeta = NF_N_H_BDDF1_3D_Zeta;
      XiArray = NF_N_H_BDDF1_3D_XiArray;
      EtaArray = NF_N_H_BDDF1_3D_EtaArray;
      ZetaArray = NF_N_H_BDDF1_3D_ZetaArray;
      T = NF_N_H_BDDF1_3D_T;
      S = NF_N_H_BDDF1_3D_S;
      EvalAll = NF_N_H_BDDF1_3D_EvalAll;
      EvalJoint = NF_N_H_BDDF1_3D_EvalFace;
      break;
    case NF_N_H_BDDF2_3D:
      N_AllFunctionals = NF_N_H_BDDF2_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_N_H_BDDF2_3D_N_FaceFunctionals;
      N_PointsAll = NF_N_H_BDDF2_3D_N_PointsAll;
      N_PointsFace = NF_N_H_BDDF2_3D_N_PointsFace;
      Xi = NF_N_H_BDDF2_3D_Xi;
      Eta = NF_N_H_BDDF2_3D_Eta;
      Zeta = NF_N_H_BDDF2_3D_Zeta;
      XiArray = NF_N_H_BDDF2_3D_XiArray;
      EtaArray = NF_N_H_BDDF2_3D_EtaArray;
      ZetaArray = NF_N_H_BDDF2_3D_ZetaArray;
      T = NF_N_H_BDDF2_3D_T;
      S = NF_N_H_BDDF2_3D_S;
      EvalAll = NF_N_H_BDDF2_3D_EvalAll;
      EvalJoint = NF_N_H_BDDF2_3D_EvalFace;
      break;
    case NF_N_H_BDDF3_3D:
      N_AllFunctionals = NF_N_H_BDDF3_3D_N_AllFunctionals;
      N_FaceFunctionals = NF_N_H_BDDF3_3D_N_FaceFunctionals;
      N_PointsAll = NF_N_H_BDDF3_3D_N_PointsAll;
      N_PointsFace = NF_N_H_BDDF3_3D_N_PointsFace;
      Xi = NF_N_H_BDDF3_3D_Xi;
      Eta = NF_N_H_BDDF3_3D_Eta;
      Zeta = NF_N_H_BDDF3_3D_Zeta;
      XiArray = NF_N_H_BDDF3_3D_XiArray;
      EtaArray = NF_N_H_BDDF3_3D_EtaArray;
      ZetaArray = NF_N_H_BDDF3_3D_ZetaArray;
      T = NF_N_H_BDDF3_3D_T;
      S = NF_N_H_BDDF3_3D_S;
      EvalAll = NF_N_H_BDDF3_3D_EvalAll;
      EvalJoint = NF_N_H_BDDF3_3D_EvalFace;
      break;
    default:
      ErrThrow("unknown NodalFunctional2D ", type);
      break;
  }
}

void NodalFunctional::GetPointsForAll(int &n_points, const double* &xi,
                        const double* &eta, const double* &zeta) const
{ 
  n_points = N_PointsAll;
  xi = Xi; 
  eta = Eta;
  zeta = Zeta;
}

void NodalFunctional::GetPointsForFace(int j, int &n_points,
                        const double* &xi, const double* &eta,
                        const double* &zeta) const
{ 
  n_points = N_PointsFace[j];
  xi   = XiArray[j];
  eta  = EtaArray[j];
  zeta = ZetaArray[j];
}

void NodalFunctional::GetPointsForFace(int &n_points, const double* &t,
                                          const double* &s) const
{
  n_points = N_PointsFace[0],
  t = T;
  s = S;
}
