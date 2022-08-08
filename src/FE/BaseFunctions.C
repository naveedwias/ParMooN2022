#include "BaseFunctions.h"
#include "BaseCell.h"
#include "Collection.h"
#include "Point.h"
#include "QuadFormula.h"
#include "MooNMD_Io.h"
#include "BF_C_L_P0_1D.h"
#include "BF_C_L_P1_1D.h"
#include "BF_C_L_P2_1D.h"
#include "BF_C_L_P3_1D.h"
#include "BF_D_L_P1_1D.h"
#include "BF_D_L_P2_1D.h"
#include "AllBaseFunctions2D.h"
#include "AllBaseFunctions3D.h"
#include <stdlib.h>

BaseFunctions::BaseFunctions(BaseFunction_type id)
{
  BaseFunct = id;
  N_BF2Change = 0;
  BF2Change = nullptr;
  BaseVectDim = 1;
  int D00 = static_cast<int>(MultiIndex2D::D00);
  int D10 = static_cast<int>(MultiIndex2D::D10);
  int D01 = static_cast<int>(MultiIndex2D::D01);
  int D20 = static_cast<int>(MultiIndex2D::D20);
  int D11 = static_cast<int>(MultiIndex2D::D11);
  int D02 = static_cast<int>(MultiIndex2D::D02);
  int D000 = static_cast<int>(MultiIndex3D::D000);
  int D100 = static_cast<int>(MultiIndex3D::D100);
  int D010 = static_cast<int>(MultiIndex3D::D010);
  int D001 = static_cast<int>(MultiIndex3D::D001);
  int D200 = static_cast<int>(MultiIndex3D::D200);
  int D110 = static_cast<int>(MultiIndex3D::D110);
  int D101 = static_cast<int>(MultiIndex3D::D101);
  int D020 = static_cast<int>(MultiIndex3D::D020);
  int D011 = static_cast<int>(MultiIndex3D::D011);
  int D002 = static_cast<int>(MultiIndex3D::D002);
  switch(BaseFunct)
  {
    case BF_C_L_P0_1D:
      Dimension = 1;
      RefElement = BFRefElements::BFUnitLine;
      Functions1D[static_cast<int>(MultiIndex1D::D0)] = C_L_P0_1D_Funct;
      Functions1D[static_cast<int>(MultiIndex1D::D1)] = C_L_P0_1D_DeriveXi;
      Functions1D[static_cast<int>(MultiIndex1D::D2)] = C_L_P0_1D_DeriveXiXi;
      PolynomialDegree = 0;
      Accuracy = 0;
      break;
    case BF_C_L_P1_1D:
      Dimension = 2;
      RefElement = BFRefElements::BFUnitLine;
      Functions1D[static_cast<int>(MultiIndex1D::D0)] = C_L_P1_1D_Funct;
      Functions1D[static_cast<int>(MultiIndex1D::D1)] = C_L_P1_1D_DeriveXi;
      Functions1D[static_cast<int>(MultiIndex1D::D2)] = C_L_P1_1D_DeriveXiXi;
      PolynomialDegree = 1;
      Accuracy = 1;
      break;
    case BF_C_L_P2_1D:
      Dimension = 3;
      RefElement = BFRefElements::BFUnitLine;
      Functions1D[static_cast<int>(MultiIndex1D::D0)] = C_L_P2_1D_Funct;
      Functions1D[static_cast<int>(MultiIndex1D::D1)] = C_L_P2_1D_DeriveXi;
      Functions1D[static_cast<int>(MultiIndex1D::D2)] = C_L_P2_1D_DeriveXiXi;
      PolynomialDegree = 2;
      Accuracy = 2;
      break;
    case BF_C_L_P3_1D:
      Dimension = 4;
      RefElement = BFRefElements::BFUnitLine;
      Functions1D[static_cast<int>(MultiIndex1D::D0)] = C_L_P3_1D_Funct;
      Functions1D[static_cast<int>(MultiIndex1D::D1)] = C_L_P3_1D_DeriveXi;
      Functions1D[static_cast<int>(MultiIndex1D::D2)] = C_L_P3_1D_DeriveXiXi;
      PolynomialDegree = 3;
      Accuracy = 3;
      break;
    case BF_D_L_P1_1D:
      Dimension = 2;
      RefElement = BFRefElements::BFUnitLine;
      Functions1D[static_cast<int>(MultiIndex1D::D0)] = D_L_P1_1D_Funct;
      Functions1D[static_cast<int>(MultiIndex1D::D1)] = D_L_P1_1D_DeriveXi;
      Functions1D[static_cast<int>(MultiIndex1D::D2)] = D_L_P1_1D_DeriveXiXi;
      PolynomialDegree = 1;
      Accuracy = 1;
      break;
    case BF_D_L_P2_1D:
      Dimension = 3;
      RefElement = BFRefElements::BFUnitLine;
      Functions1D[static_cast<int>(MultiIndex1D::D0)] = D_L_P2_1D_Funct;
      Functions1D[static_cast<int>(MultiIndex1D::D1)] = D_L_P2_1D_DeriveXi;
      Functions1D[static_cast<int>(MultiIndex1D::D2)] = D_L_P2_1D_DeriveXiXi;
      PolynomialDegree = 2;
      Accuracy = 2;
      break;
    case BF_C_T_P00_2D:
      Dimension = 1;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = C_T_P00_2D_Funct;
      Functions[D10] = C_T_P00_2D_DeriveXi;
      Functions[D01] = C_T_P00_2D_DeriveEta;
      Functions[D20] = C_T_P00_2D_DeriveXiXi;
      Functions[D11] = C_T_P00_2D_DeriveXiEta;
      Functions[D02] = C_T_P00_2D_DeriveEtaEta;
      PolynomialDegree = 0;
      Accuracy = 0;
      break;
    case BF_C_T_P0_2D:
      Dimension = 1;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = C_T_P0_2D_Funct;
      Functions[D10] = C_T_P0_2D_DeriveXi;
      Functions[D01] = C_T_P0_2D_DeriveEta;
      Functions[D20] = C_T_P0_2D_DeriveXiXi;
      Functions[D11] = C_T_P0_2D_DeriveXiEta;
      Functions[D02] = C_T_P0_2D_DeriveEtaEta;
      PolynomialDegree = 0;
      Accuracy = 0;
      break;
    case BF_C_T_P1_2D:
      Dimension = 3;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = C_T_P1_2D_Funct;
      Functions[D10] = C_T_P1_2D_DeriveXi;
      Functions[D01] = C_T_P1_2D_DeriveEta;
      Functions[D20] = C_T_P1_2D_DeriveXiXi;
      Functions[D11] = C_T_P1_2D_DeriveXiEta;
      Functions[D02] = C_T_P1_2D_DeriveEtaEta;
      PolynomialDegree = 1;
      Accuracy = 1;
      break;
    case BF_C_T_P2_2D:
      Dimension = 6;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = C_T_P2_2D_Funct;
      Functions[D10] = C_T_P2_2D_DeriveXi;
      Functions[D01] = C_T_P2_2D_DeriveEta;
      Functions[D20] = C_T_P2_2D_DeriveXiXi;
      Functions[D11] = C_T_P2_2D_DeriveXiEta;
      Functions[D02] = C_T_P2_2D_DeriveEtaEta;
      PolynomialDegree = 2;
      Accuracy = 2;
      break;
    case BF_C_T_P3_2D:
      Dimension = 10;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = C_T_P3_2D_Funct;
      Functions[D10] = C_T_P3_2D_DeriveXi;
      Functions[D01] = C_T_P3_2D_DeriveEta;
      Functions[D20] = C_T_P3_2D_DeriveXiXi;
      Functions[D11] = C_T_P3_2D_DeriveXiEta;
      Functions[D02] = C_T_P3_2D_DeriveEtaEta;
      PolynomialDegree = 3;
      Accuracy = 3;
      break;
    case BF_N_T_P1_2D:
      Dimension = 3;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = N_T_P1_2D_Funct;
      Functions[D10] = N_T_P1_2D_DeriveXi;
      Functions[D01] = N_T_P1_2D_DeriveEta;
      Functions[D20] = N_T_P1_2D_DeriveXiXi;
      Functions[D11] = N_T_P1_2D_DeriveXiEta;
      Functions[D02] = N_T_P1_2D_DeriveEtaEta;
      PolynomialDegree = 1;
      Accuracy = 1;
      break;
    case BF_C_Q_Q00_2D:
      Dimension = 1;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_Q00_2D_Funct;
      Functions[D10] = C_Q_Q00_2D_DeriveXi;
      Functions[D01] = C_Q_Q00_2D_DeriveEta;
      Functions[D20] = C_Q_Q00_2D_DeriveXiXi;
      Functions[D11] = C_Q_Q00_2D_DeriveXiEta;
      Functions[D02] = C_Q_Q00_2D_DeriveEtaEta;
      PolynomialDegree = 0;
      Accuracy = 0;
      break;
    case BF_C_Q_Q0_2D:
      Dimension = 1;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_Q0_2D_Funct;
      Functions[D10] = C_Q_Q0_2D_DeriveXi;
      Functions[D01] = C_Q_Q0_2D_DeriveEta;
      Functions[D20] = C_Q_Q0_2D_DeriveXiXi;
      Functions[D11] = C_Q_Q0_2D_DeriveXiEta;
      Functions[D02] = C_Q_Q0_2D_DeriveEtaEta;
      PolynomialDegree = 0;
      Accuracy = 0;
      break;
    case BF_C_Q_Q1_2D:
      Dimension = 4;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_Q1_2D_Funct;
      Functions[D10] = C_Q_Q1_2D_DeriveXi;
      Functions[D01] = C_Q_Q1_2D_DeriveEta;
      Functions[D20] = C_Q_Q1_2D_DeriveXiXi;
      Functions[D11] = C_Q_Q1_2D_DeriveXiEta;
      Functions[D02] = C_Q_Q1_2D_DeriveEtaEta;
      PolynomialDegree = 1;
      Accuracy = 1;
      break;
    case BF_C_Q_Q2_2D:
      Dimension = 9;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_Q2_2D_Funct;
      Functions[D10] = C_Q_Q2_2D_DeriveXi;
      Functions[D01] = C_Q_Q2_2D_DeriveEta;
      Functions[D20] = C_Q_Q2_2D_DeriveXiXi;
      Functions[D11] = C_Q_Q2_2D_DeriveXiEta;
      Functions[D02] = C_Q_Q2_2D_DeriveEtaEta;
      PolynomialDegree = 2;
      Accuracy = 2;
      break;
    case BF_C_Q_Q3_2D:
      Dimension = 16;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_Q3_2D_Funct;
      Functions[D10] = C_Q_Q3_2D_DeriveXi;
      Functions[D01] = C_Q_Q3_2D_DeriveEta;
      Functions[D20] = C_Q_Q3_2D_DeriveXiXi;
      Functions[D11] = C_Q_Q3_2D_DeriveXiEta;
      Functions[D02] = C_Q_Q3_2D_DeriveEtaEta;
      PolynomialDegree = 3;
      Accuracy = 3;
      break;
    case BF_N_Q_Q1_2D:
      Dimension = 4;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = N_Q_Q1_2D_Funct;
      Functions[D10] = N_Q_Q1_2D_DeriveXi;
      Functions[D01] = N_Q_Q1_2D_DeriveEta;
      Functions[D20] = N_Q_Q1_2D_DeriveXiXi;
      Functions[D11] = N_Q_Q1_2D_DeriveXiEta;
      Functions[D02] = N_Q_Q1_2D_DeriveEtaEta;
      PolynomialDegree = 2;
      Accuracy = 1;
      break;
    case BF_D_Q_P1_2D:
      Dimension = 3;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = D_Q_P1_2D_Funct;
      Functions[D10] = D_Q_P1_2D_DeriveXi;
      Functions[D01] = D_Q_P1_2D_DeriveEta;
      Functions[D20] = D_Q_P1_2D_DeriveXiXi;
      Functions[D11] = D_Q_P1_2D_DeriveXiEta;
      Functions[D02] = D_Q_P1_2D_DeriveEtaEta;
      PolynomialDegree = 1;
      Accuracy = 1;
      break;
    case BF_D_Q_P2_2D:
      Dimension = 6;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = D_Q_P2_2D_Funct;
      Functions[D10] = D_Q_P2_2D_DeriveXi;
      Functions[D01] = D_Q_P2_2D_DeriveEta;
      Functions[D20] = D_Q_P2_2D_DeriveXiXi;
      Functions[D11] = D_Q_P2_2D_DeriveXiEta;
      Functions[D02] = D_Q_P2_2D_DeriveEtaEta;
      PolynomialDegree = 2;
      Accuracy = 2;
      break;
    case BF_D_Q_P3_2D:
      Dimension = 10;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = D_Q_P3_2D_Funct;
      Functions[D10] = D_Q_P3_2D_DeriveXi;
      Functions[D01] = D_Q_P3_2D_DeriveEta;
      Functions[D20] = D_Q_P3_2D_DeriveXiXi;
      Functions[D11] = D_Q_P3_2D_DeriveXiEta;
      Functions[D02] = D_Q_P3_2D_DeriveEtaEta;
      PolynomialDegree = 3;
      Accuracy = 3;
      break;
    case BF_C_T_P4_2D:
      Dimension = 15;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = C_T_P4_2D_Funct;
      Functions[D10] = C_T_P4_2D_DeriveXi;
      Functions[D01] = C_T_P4_2D_DeriveEta;
      Functions[D20] = C_T_P4_2D_DeriveXiXi;
      Functions[D11] = C_T_P4_2D_DeriveXiEta;
      Functions[D02] = C_T_P4_2D_DeriveEtaEta;
      PolynomialDegree = 4;
      Accuracy = 4;
      break;
    case BF_C_T_P5_2D:
      Dimension = 21;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = C_T_P5_2D_Funct;
      Functions[D10] = C_T_P5_2D_DeriveXi;
      Functions[D01] = C_T_P5_2D_DeriveEta;
      Functions[D20] = C_T_P5_2D_DeriveXiXi;
      Functions[D11] = C_T_P5_2D_DeriveXiEta;
      Functions[D02] = C_T_P5_2D_DeriveEtaEta;
      PolynomialDegree = 5;
      Accuracy = 5;
      break;
    case BF_C_T_P6_2D:
      Dimension = 28;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = C_T_P6_2D_Funct;
      Functions[D10] = C_T_P6_2D_DeriveXi;
      Functions[D01] = C_T_P6_2D_DeriveEta;
      Functions[D20] = C_T_P6_2D_DeriveXiXi;
      Functions[D11] = C_T_P6_2D_DeriveXiEta;
      Functions[D02] = C_T_P6_2D_DeriveEtaEta;
      PolynomialDegree = 6;
      Accuracy = 6;
      break;
    case BF_C_T_P7_2D:
      ErrThrow("Basis Functions of type BF_C_T_P7_2D not implemented");
      break;
    case BF_C_T_P8_2D:
      ErrThrow("Basis Functions of type BF_C_T_P8_2D not implemented");
      break;
    case BF_C_T_P9_2D:
      ErrThrow("Basis Functions of type BF_C_T_P9_2D not implemented");
      break;
    case BF_C_Q_Q4_2D:
      Dimension = 25;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_Q4_2D_Funct;
      Functions[D10] = C_Q_Q4_2D_DeriveXi;
      Functions[D01] = C_Q_Q4_2D_DeriveEta;
      Functions[D20] = C_Q_Q4_2D_DeriveXiXi;
      Functions[D11] = C_Q_Q4_2D_DeriveXiEta;
      Functions[D02] = C_Q_Q4_2D_DeriveEtaEta;
      PolynomialDegree = 4;
      Accuracy = 4;
      break;
    case BF_C_Q_Q5_2D:
      Dimension = 36;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_Q5_2D_Funct;
      Functions[D10] = C_Q_Q5_2D_DeriveXi;
      Functions[D01] = C_Q_Q5_2D_DeriveEta;
      Functions[D20] = C_Q_Q5_2D_DeriveXiXi;
      Functions[D11] = C_Q_Q5_2D_DeriveXiEta;
      Functions[D02] = C_Q_Q5_2D_DeriveEtaEta;
      PolynomialDegree = 5;
      Accuracy = 5;
      break;
    case BF_C_Q_Q6_2D:
      Dimension = 49;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_Q6_2D_Funct;
      Functions[D10] = C_Q_Q6_2D_DeriveXi;
      Functions[D01] = C_Q_Q6_2D_DeriveEta;
      Functions[D20] = C_Q_Q6_2D_DeriveXiXi;
      Functions[D11] = C_Q_Q6_2D_DeriveXiEta;
      Functions[D02] = C_Q_Q6_2D_DeriveEtaEta;
      PolynomialDegree = 6;
      Accuracy = 6;
      break;
    case BF_C_Q_Q7_2D:
      Dimension = 64;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_Q7_2D_Funct;
      Functions[D10] = C_Q_Q7_2D_DeriveXi;
      Functions[D01] = C_Q_Q7_2D_DeriveEta;
      Functions[D20] = C_Q_Q7_2D_DeriveXiXi;
      Functions[D11] = C_Q_Q7_2D_DeriveXiEta;
      Functions[D02] = C_Q_Q7_2D_DeriveEtaEta;
      PolynomialDegree = 7;
      Accuracy = 7;
      break;
    case BF_C_Q_Q8_2D:
      Dimension = 81;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_Q8_2D_Funct;
      Functions[D10] = C_Q_Q8_2D_DeriveXi;
      Functions[D01] = C_Q_Q8_2D_DeriveEta;
      Functions[D20] = C_Q_Q8_2D_DeriveXiXi;
      Functions[D11] = C_Q_Q8_2D_DeriveXiEta;
      Functions[D02] = C_Q_Q8_2D_DeriveEtaEta;
      PolynomialDegree = 8;
      Accuracy = 8;
      break;
    case BF_C_Q_Q9_2D:
      Dimension = 100;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_Q9_2D_Funct;
      Functions[D10] = C_Q_Q9_2D_DeriveXi;
      Functions[D01] = C_Q_Q9_2D_DeriveEta;
      Functions[D20] = C_Q_Q9_2D_DeriveXiXi;
      Functions[D11] = C_Q_Q9_2D_DeriveXiEta;
      Functions[D02] = C_Q_Q9_2D_DeriveEtaEta;
      PolynomialDegree = 9;
      Accuracy = 9;
      break;
    case BF_C_T_B2_2D:
      Dimension = 7;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = C_T_B2_2D_Funct;
      Functions[D10] = C_T_B2_2D_DeriveXi;
      Functions[D01] = C_T_B2_2D_DeriveEta;
      Functions[D20] = C_T_B2_2D_DeriveXiXi;
      Functions[D11] = C_T_B2_2D_DeriveXiEta;
      Functions[D02] = C_T_B2_2D_DeriveEtaEta;
      PolynomialDegree = 3;
      Accuracy = 2;
      break;
    case BF_C_T_B3_2D:
      Dimension = 12;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = C_T_B3_2D_Funct;
      Functions[D10] = C_T_B3_2D_DeriveXi;
      Functions[D01] = C_T_B3_2D_DeriveEta;
      Functions[D20] = C_T_B3_2D_DeriveXiXi;
      Functions[D11] = C_T_B3_2D_DeriveXiEta;
      Functions[D02] = C_T_B3_2D_DeriveEtaEta;
      PolynomialDegree = 4;
      Accuracy = 3;
      break;
    case BF_C_T_SV2_2D:
      Dimension = 10;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = C_T_SV2_2D_Funct;
      Functions[D10] = C_T_SV2_2D_DeriveXi;
      Functions[D01] = C_T_SV2_2D_DeriveEta;
      Functions[D20] = C_T_SV2_2D_DeriveXiXi;
      Functions[D11] = C_T_SV2_2D_DeriveXiEta;
      Functions[D02] = C_T_SV2_2D_DeriveEtaEta;
      PolynomialDegree = 2;
      Accuracy = 2;
      break;
    case BF_D_T_P1_2D:
      Dimension = 3;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = D_T_P1_2D_Funct;
      Functions[D10] = D_T_P1_2D_DeriveXi;
      Functions[D01] = D_T_P1_2D_DeriveEta;
      Functions[D20] = D_T_P1_2D_DeriveXiXi;
      Functions[D11] = D_T_P1_2D_DeriveXiEta;
      Functions[D02] = D_T_P1_2D_DeriveEtaEta;
      PolynomialDegree = 1;
      Accuracy = 1;
      break;
    case BF_D_T_P2_2D:
      Dimension = 6;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = D_T_P2_2D_Funct;
      Functions[D10] = D_T_P2_2D_DeriveXi;
      Functions[D01] = D_T_P2_2D_DeriveEta;
      Functions[D20] = D_T_P2_2D_DeriveXiXi;
      Functions[D11] = D_T_P2_2D_DeriveXiEta;
      Functions[D02] = D_T_P2_2D_DeriveEtaEta;
      PolynomialDegree = 2;
      Accuracy = 2;
      break;
    case BF_D_T_SV1_2D:
      Dimension = 9;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = D_T_SV1_2D_Funct;
      Functions[D10] = D_T_SV1_2D_DeriveXi;
      Functions[D01] = D_T_SV1_2D_DeriveEta;
      Functions[D20] = D_T_SV1_2D_DeriveXiXi;
      Functions[D11] = D_T_SV1_2D_DeriveXiEta;
      Functions[D02] = D_T_SV1_2D_DeriveEtaEta;
      PolynomialDegree = 1;
      Accuracy = 1;
      break;
    case BF_N_Q_Q2_2D:
      Dimension = 9;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = N_Q_Q2_2D_Funct;
      Functions[D10] = N_Q_Q2_2D_DeriveXi;
      Functions[D01] = N_Q_Q2_2D_DeriveEta;
      Functions[D20] = N_Q_Q2_2D_DeriveXiXi;
      Functions[D11] = N_Q_Q2_2D_DeriveXiEta;
      Functions[D02] = N_Q_Q2_2D_DeriveEtaEta;
      PolynomialDegree = 3;
      Accuracy = 2;
      N_BF2Change = 1;
      BF2Change = N_Q_Q2_2D_Change;
      break;
    case BF_N_Q_Q3_2D:
      Dimension = 15;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = N_Q_Q3_2D_Funct;
      Functions[D10] = N_Q_Q3_2D_DeriveXi;
      Functions[D01] = N_Q_Q3_2D_DeriveEta;
      Functions[D20] = N_Q_Q3_2D_DeriveXiXi;
      Functions[D11] = N_Q_Q3_2D_DeriveXiEta;
      Functions[D02] = N_Q_Q3_2D_DeriveEtaEta;
      PolynomialDegree = 4;
      Accuracy = 3;
      N_BF2Change = 1;
      BF2Change = N_Q_Q3_2D_Change;
      break;
    case BF_N_Q_Q4_2D:
      Dimension = 22;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = N_Q_Q4_2D_Funct;
      Functions[D10] = N_Q_Q4_2D_DeriveXi;
      Functions[D01] = N_Q_Q4_2D_DeriveEta;
      Functions[D20] = N_Q_Q4_2D_DeriveXiXi;
      Functions[D11] = N_Q_Q4_2D_DeriveXiEta;
      Functions[D02] = N_Q_Q4_2D_DeriveEtaEta;
      PolynomialDegree = 5;
      Accuracy = 4;
      N_BF2Change = 2;
      BF2Change = N_Q_Q4_2D_Change;
      break;
    case BF_N_Q_Q5_2D:
      Dimension = 30;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = N_Q_Q5_2D_Funct;
      Functions[D10] = N_Q_Q5_2D_DeriveXi;
      Functions[D01] = N_Q_Q5_2D_DeriveEta;
      Functions[D20] = N_Q_Q5_2D_DeriveXiXi;
      Functions[D11] = N_Q_Q5_2D_DeriveXiEta;
      Functions[D02] = N_Q_Q5_2D_DeriveEtaEta;
      PolynomialDegree = 6;
      Accuracy = 5;
      N_BF2Change = 2;
      BF2Change = N_Q_Q5_2D_Change;
      break;
    case BF_D_Q_P4_2D:
      Dimension = 15;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = D_Q_P4_2D_Funct;
      Functions[D10] = D_Q_P4_2D_DeriveXi;
      Functions[D01] = D_Q_P4_2D_DeriveEta;
      Functions[D20] = D_Q_P4_2D_DeriveXiXi;
      Functions[D11] = D_Q_P4_2D_DeriveXiEta;
      Functions[D02] = D_Q_P4_2D_DeriveEtaEta;
      PolynomialDegree = 4;
      Accuracy = 4;
      break;
    case BF_D_Q_P5_2D:
      Dimension = 21;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = D_Q_P5_2D_Funct;
      Functions[D10] = D_Q_P5_2D_DeriveXi;
      Functions[D01] = D_Q_P5_2D_DeriveEta;
      Functions[D20] = D_Q_P5_2D_DeriveXiXi;
      Functions[D11] = D_Q_P5_2D_DeriveXiEta;
      Functions[D02] = D_Q_P5_2D_DeriveEtaEta;
      PolynomialDegree = 5;
      Accuracy = 5;
      break;
    case BF_D_Q_P6_2D:
      Dimension = 28;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = D_Q_P6_2D_Funct;
      Functions[D10] = D_Q_P6_2D_DeriveXi;
      Functions[D01] = D_Q_P6_2D_DeriveEta;
      Functions[D20] = D_Q_P6_2D_DeriveXiXi;
      Functions[D11] = D_Q_P6_2D_DeriveXiEta;
      Functions[D02] = D_Q_P6_2D_DeriveEtaEta;
      PolynomialDegree = 6;
      Accuracy = 6;
      break;
    case BF_D_Q_P7_2D:
      Dimension = 36;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = D_Q_P7_2D_Funct;
      Functions[D10] = D_Q_P7_2D_DeriveXi;
      Functions[D01] = D_Q_P7_2D_DeriveEta;
      Functions[D20] = D_Q_P7_2D_DeriveXiXi;
      Functions[D11] = D_Q_P7_2D_DeriveXiEta;
      Functions[D02] = D_Q_P7_2D_DeriveEtaEta;
      PolynomialDegree = 7;
      Accuracy = 7;
      break;
    case BF_N_T_P1MOD_2D:
      Dimension = 6;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = N_T_P1MOD_2D_Funct;
      Functions[D10] = N_T_P1MOD_2D_DeriveXi;
      Functions[D01] = N_T_P1MOD_2D_DeriveEta;
      Functions[D20] = N_T_P1MOD_2D_DeriveXiXi;
      Functions[D11] = N_T_P1MOD_2D_DeriveXiEta;
      Functions[D02] = N_T_P1MOD_2D_DeriveEtaEta;
      PolynomialDegree = 3;
      Accuracy = 1;
      N_BF2Change = 1;
      BF2Change = N_T_P1MOD_2D_Change;
      break;
    case BF_N_T_P2_2D:
      Dimension = 7;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = N_T_P2_2D_Funct;
      Functions[D10] = N_T_P2_2D_DeriveXi;
      Functions[D01] = N_T_P2_2D_DeriveEta;
      Functions[D20] = N_T_P2_2D_DeriveXiXi;
      Functions[D11] = N_T_P2_2D_DeriveXiEta;
      Functions[D02] = N_T_P2_2D_DeriveEtaEta;
      PolynomialDegree = 3;
      Accuracy = 2;
      N_BF2Change = 1;
      BF2Change = N_T_P2_2D_Change;
      break;
    case BF_N_T_P3_2D:
      Dimension = 12;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = N_T_P3_2D_Funct;
      Functions[D10] = N_T_P3_2D_DeriveXi;
      Functions[D01] = N_T_P3_2D_DeriveEta;
      Functions[D20] = N_T_P3_2D_DeriveXiXi;
      Functions[D11] = N_T_P3_2D_DeriveXiEta;
      Functions[D02] = N_T_P3_2D_DeriveEtaEta;
      PolynomialDegree = 4;
      Accuracy = 3;
      N_BF2Change = 1;
      BF2Change = N_T_P3_2D_Change;
      break;
    case BF_N_T_P4_2D:
      Dimension = 18;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = N_T_P4_2D_Funct;
      Functions[D10] = N_T_P4_2D_DeriveXi;
      Functions[D01] = N_T_P4_2D_DeriveEta;
      Functions[D20] = N_T_P4_2D_DeriveXiXi;
      Functions[D11] = N_T_P4_2D_DeriveXiEta;
      Functions[D02] = N_T_P4_2D_DeriveEtaEta;
      PolynomialDegree = 5;
      Accuracy = 4;
      N_BF2Change = 2;
      BF2Change = N_T_P4_2D_Change;
      break;
    case BF_N_T_P5_2D:
      Dimension = 25;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = N_T_P5_2D_Funct;
      Functions[D10] = N_T_P5_2D_DeriveXi;
      Functions[D01] = N_T_P5_2D_DeriveEta;
      Functions[D20] = N_T_P5_2D_DeriveXiXi;
      Functions[D11] = N_T_P5_2D_DeriveXiEta;
      Functions[D02] = N_T_P5_2D_DeriveEtaEta;
      PolynomialDegree = 6;
      Accuracy = 5;
      N_BF2Change = 2;
      BF2Change = N_T_P5_2D_Change;
      break;
    case BF_D_T_P3_2D:
      Dimension = 10;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = D_T_P3_2D_Funct;
      Functions[D10] = D_T_P3_2D_DeriveXi;
      Functions[D01] = D_T_P3_2D_DeriveEta;
      Functions[D20] = D_T_P3_2D_DeriveXiXi;
      Functions[D11] = D_T_P3_2D_DeriveXiEta;
      Functions[D02] = D_T_P3_2D_DeriveEtaEta;
      PolynomialDegree = 3;
      Accuracy = 3;
      break;
    case BF_D_T_P4_2D:
      Dimension = 15;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = D_T_P4_2D_Funct;
      Functions[D10] = D_T_P4_2D_DeriveXi;
      Functions[D01] = D_T_P4_2D_DeriveEta;
      Functions[D20] = D_T_P4_2D_DeriveXiXi;
      Functions[D11] = D_T_P4_2D_DeriveXiEta;
      Functions[D02] = D_T_P4_2D_DeriveEtaEta;
      PolynomialDegree = 4;
      Accuracy = 4;
      break;
    case BF_C_T_P1MINI_2D:
      Dimension = 4;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = C_T_P1MINI_2D_Funct;
      Functions[D10] = C_T_P1MINI_2D_DeriveXi;
      Functions[D01] = C_T_P1MINI_2D_DeriveEta;
      Functions[D20] = C_T_P1MINI_2D_DeriveXiXi;
      Functions[D11] = C_T_P1MINI_2D_DeriveXiEta;
      Functions[D02] = C_T_P1MINI_2D_DeriveEtaEta;
      PolynomialDegree = 3;
      Accuracy = 1;
      break;
    case BF_B_Q_IB2_2D:
      Dimension = 1;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = B_Q_IB2_2D_Funct;
      Functions[D10] = B_Q_IB2_2D_DeriveXi;
      Functions[D01] = B_Q_IB2_2D_DeriveEta;
      Functions[D20] = B_Q_IB2_2D_DeriveXiXi;
      Functions[D11] = B_Q_IB2_2D_DeriveXiEta;
      Functions[D02] = B_Q_IB2_2D_DeriveEtaEta;
      PolynomialDegree = 2;
      Accuracy = 2;
      break;
    case BF_D_Q_Q1_2D:
      Dimension = 4;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = D_Q_Q1_2D_Funct;
      Functions[D10] = D_Q_Q1_2D_DeriveXi;
      Functions[D01] = D_Q_Q1_2D_DeriveEta;
      Functions[D20] = D_Q_Q1_2D_DeriveXiXi;
      Functions[D11] = D_Q_Q1_2D_DeriveXiEta;
      Functions[D02] = D_Q_Q1_2D_DeriveEtaEta;
      PolynomialDegree = 1;
      Accuracy = 1;
      break;
    case BF_D_Q_Q2_2D:
      Dimension = 9;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = D_Q_Q2_2D_Funct;
      Functions[D10] = D_Q_Q2_2D_DeriveXi;
      Functions[D01] = D_Q_Q2_2D_DeriveEta;
      Functions[D20] = D_Q_Q2_2D_DeriveXiXi;
      Functions[D11] = D_Q_Q2_2D_DeriveXiEta;
      Functions[D02] = D_Q_Q2_2D_DeriveEtaEta;
      PolynomialDegree = 2;
      Accuracy = 2;
      break;
    case BF_D_Q_Q3_2D:
      Dimension = 16;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = D_Q_Q3_2D_Funct;
      Functions[D10] = D_Q_Q3_2D_DeriveXi;
      Functions[D01] = D_Q_Q3_2D_DeriveEta;
      Functions[D20] = D_Q_Q3_2D_DeriveXiXi;
      Functions[D11] = D_Q_Q3_2D_DeriveXiEta;
      Functions[D02] = D_Q_Q3_2D_DeriveEtaEta;
      PolynomialDegree = 3;
      Accuracy = 3;
      break;
    case BF_D_Q_Q4_2D:
      Dimension = 25;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = D_Q_Q4_2D_Funct;
      Functions[D10] = D_Q_Q4_2D_DeriveXi;
      Functions[D01] = D_Q_Q4_2D_DeriveEta;
      Functions[D20] = D_Q_Q4_2D_DeriveXiXi;
      Functions[D11] = D_Q_Q4_2D_DeriveXiEta;
      Functions[D02] = D_Q_Q4_2D_DeriveEtaEta;
      PolynomialDegree = 4;
      Accuracy = 4;
      break;
    case BF_C_T_B4_2D:
      Dimension = 18;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = C_T_B4_2D_Funct;
      Functions[D10] = C_T_B4_2D_DeriveXi;
      Functions[D01] = C_T_B4_2D_DeriveEta;
      Functions[D20] = C_T_B4_2D_DeriveXiXi;
      Functions[D11] = C_T_B4_2D_DeriveXiEta;
      Functions[D02] = C_T_B4_2D_DeriveEtaEta;
      PolynomialDegree = 5;
      Accuracy = 4;
      break;
    case BF_D_Q_D2_2D:
      Dimension = 9;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = D_Q_D2_2D_Funct;
      Functions[D10] = D_Q_D2_2D_DeriveXi;
      Functions[D01] = D_Q_D2_2D_DeriveEta;
      Functions[D20] = D_Q_D2_2D_DeriveXiXi;
      Functions[D11] = D_Q_D2_2D_DeriveXiEta;
      Functions[D02] = D_Q_D2_2D_DeriveEtaEta;
      PolynomialDegree = 3;
      Accuracy = 2;
      break;
    case BF_C_Q_UL1_2D:
      Dimension = 5;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_UL1_2D_Funct;
      Functions[D10] = C_Q_UL1_2D_DeriveXi;
      Functions[D01] = C_Q_UL1_2D_DeriveEta;
      Functions[D20] = C_Q_UL1_2D_DeriveXiXi;
      Functions[D11] = C_Q_UL1_2D_DeriveXiEta;
      Functions[D02] = C_Q_UL1_2D_DeriveEtaEta;
      PolynomialDegree = 2;
      Accuracy = 1;
      break;
    case BF_C_Q_UL2_2D:
      Dimension = 11;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_UL2_2D_Funct;
      Functions[D10] = C_Q_UL2_2D_DeriveXi;
      Functions[D01] = C_Q_UL2_2D_DeriveEta;
      Functions[D20] = C_Q_UL2_2D_DeriveXiXi;
      Functions[D11] = C_Q_UL2_2D_DeriveXiEta;
      Functions[D02] = C_Q_UL2_2D_DeriveEtaEta;
      PolynomialDegree = 3;
      Accuracy = 2;
      break;
    case BF_C_Q_UL3_2D:
      Dimension = 18;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_UL3_2D_Funct;
      Functions[D10] = C_Q_UL3_2D_DeriveXi;
      Functions[D01] = C_Q_UL3_2D_DeriveEta;
      Functions[D20] = C_Q_UL3_2D_DeriveXiXi;
      Functions[D11] = C_Q_UL3_2D_DeriveXiEta;
      Functions[D02] = C_Q_UL3_2D_DeriveEtaEta;
      PolynomialDegree = 4;
      Accuracy = 3;
      break;
    case BF_C_Q_UL4_2D:
      Dimension = 27;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_UL4_2D_Funct;
      Functions[D10] = C_Q_UL4_2D_DeriveXi;
      Functions[D01] = C_Q_UL4_2D_DeriveEta;
      Functions[D20] = C_Q_UL4_2D_DeriveXiXi;
      Functions[D11] = C_Q_UL4_2D_DeriveXiEta;
      Functions[D02] = C_Q_UL4_2D_DeriveEtaEta;
      PolynomialDegree = 5;
      Accuracy = 4;
      break;
    case BF_C_Q_UL5_2D:
      Dimension = 38;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_UL5_2D_Funct;
      Functions[D10] = C_Q_UL5_2D_DeriveXi;
      Functions[D01] = C_Q_UL5_2D_DeriveEta;
      Functions[D20] = C_Q_UL5_2D_DeriveXiXi;
      Functions[D11] = C_Q_UL5_2D_DeriveXiEta;
      Functions[D02] = C_Q_UL5_2D_DeriveEtaEta;
      PolynomialDegree = 6;
      Accuracy = 5;
      break;
    case BF_C_T_UL1_2D:
      Dimension = 4;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = C_T_UL1_2D_Funct;
      Functions[D10] = C_T_UL1_2D_DeriveXi;
      Functions[D01] = C_T_UL1_2D_DeriveEta;
      Functions[D20] = C_T_UL1_2D_DeriveXiXi;
      Functions[D11] = C_T_UL1_2D_DeriveXiEta;
      Functions[D02] = C_T_UL1_2D_DeriveEtaEta;
      PolynomialDegree = 3;
      Accuracy = 1;
      break;
    case BF_C_T_UL2_2D:
      Dimension = 9;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = C_T_UL2_2D_Funct;
      Functions[D10] = C_T_UL2_2D_DeriveXi;
      Functions[D01] = C_T_UL2_2D_DeriveEta;
      Functions[D20] = C_T_UL2_2D_DeriveXiXi;
      Functions[D11] = C_T_UL2_2D_DeriveXiEta;
      Functions[D02] = C_T_UL2_2D_DeriveEtaEta;
      PolynomialDegree = 4;
      Accuracy = 2;
      break;
    case BF_C_T_UL3_2D:
      Dimension = 15;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = C_T_UL3_2D_Funct;
      Functions[D10] = C_T_UL3_2D_DeriveXi;
      Functions[D01] = C_T_UL3_2D_DeriveEta;
      Functions[D20] = C_T_UL3_2D_DeriveXiXi;
      Functions[D11] = C_T_UL3_2D_DeriveXiEta;
      Functions[D02] = C_T_UL3_2D_DeriveEtaEta;
      PolynomialDegree = 5;
      Accuracy = 3;
      break;
    case BF_C_T_UL4_2D:
      Dimension = 22;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = C_T_UL4_2D_Funct;
      Functions[D10] = C_T_UL4_2D_DeriveXi;
      Functions[D01] = C_T_UL4_2D_DeriveEta;
      Functions[D20] = C_T_UL4_2D_DeriveXiXi;
      Functions[D11] = C_T_UL4_2D_DeriveXiEta;
      Functions[D02] = C_T_UL4_2D_DeriveEtaEta;
      PolynomialDegree = 6;
      Accuracy = 4;
      break;
    case BF_C_T_UL5_2D:
      Dimension = 30;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = C_T_UL5_2D_Funct;
      Functions[D10] = C_T_UL5_2D_DeriveXi;
      Functions[D01] = C_T_UL5_2D_DeriveEta;
      Functions[D20] = C_T_UL5_2D_DeriveXiXi;
      Functions[D11] = C_T_UL5_2D_DeriveXiEta;
      Functions[D02] = C_T_UL5_2D_DeriveEtaEta;
      PolynomialDegree = 7;
      Accuracy = 5;
      break;
    case BF_C_Q_UL2S_2D:
      Dimension = 9;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_UL2S_2D_Funct;
      Functions[D10] = C_Q_UL2S_2D_DeriveXi;
      Functions[D01] = C_Q_UL2S_2D_DeriveEta;
      Functions[D20] = C_Q_UL2S_2D_DeriveXiXi;
      Functions[D11] = C_Q_UL2S_2D_DeriveXiEta;
      Functions[D02] = C_Q_UL2S_2D_DeriveEtaEta;
      PolynomialDegree = 2;
      Accuracy = 2;
      break;
    case BF_C_Q_UL3S_2D:
      Dimension = 15;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_UL3S_2D_Funct;
      Functions[D10] = C_Q_UL3S_2D_DeriveXi;
      Functions[D01] = C_Q_UL3S_2D_DeriveEta;
      Functions[D20] = C_Q_UL3S_2D_DeriveXiXi;
      Functions[D11] = C_Q_UL3S_2D_DeriveXiEta;
      Functions[D02] = C_Q_UL3S_2D_DeriveEtaEta;
      PolynomialDegree = 3;
      Accuracy = 3;
      break;
    case BF_C_Q_UL4S_2D:
      Dimension = 22;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_UL4S_2D_Funct;
      Functions[D10] = C_Q_UL4S_2D_DeriveXi;
      Functions[D01] = C_Q_UL4S_2D_DeriveEta;
      Functions[D20] = C_Q_UL4S_2D_DeriveXiXi;
      Functions[D11] = C_Q_UL4S_2D_DeriveXiEta;
      Functions[D02] = C_Q_UL4S_2D_DeriveEtaEta;
      PolynomialDegree = 4;
      Accuracy = 4;
      break;
    case BF_C_Q_UL5S_2D:
      Dimension = 31;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_UL5S_2D_Funct;
      Functions[D10] = C_Q_UL5S_2D_DeriveXi;
      Functions[D01] = C_Q_UL5S_2D_DeriveEta;
      Functions[D20] = C_Q_UL5S_2D_DeriveXiXi;
      Functions[D11] = C_Q_UL5S_2D_DeriveXiEta;
      Functions[D02] = C_Q_UL5S_2D_DeriveEtaEta;
      PolynomialDegree = 5;
      Accuracy = 5;
      break;
    case BF_C_Q_UL6S_2D:
      Dimension = 42;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_UL6S_2D_Funct;
      Functions[D10] = C_Q_UL6S_2D_DeriveXi;
      Functions[D01] = C_Q_UL6S_2D_DeriveEta;
      Functions[D20] = C_Q_UL6S_2D_DeriveXiXi;
      Functions[D11] = C_Q_UL6S_2D_DeriveXiEta;
      Functions[D02] = C_Q_UL6S_2D_DeriveEtaEta;
      PolynomialDegree = 6;
      Accuracy = 6;
      break;
    case BF_C_Q_UL7S_2D:
      Dimension = 55;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_UL7S_2D_Funct;
      Functions[D10] = C_Q_UL7S_2D_DeriveXi;
      Functions[D01] = C_Q_UL7S_2D_DeriveEta;
      Functions[D20] = C_Q_UL7S_2D_DeriveXiXi;
      Functions[D11] = C_Q_UL7S_2D_DeriveXiEta;
      Functions[D02] = C_Q_UL7S_2D_DeriveEtaEta;
      PolynomialDegree = 7;
      Accuracy = 7;
      break;
    case BF_C_Q_UL8S_2D:
      Dimension = 70;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_UL8S_2D_Funct;
      Functions[D10] = C_Q_UL8S_2D_DeriveXi;
      Functions[D01] = C_Q_UL8S_2D_DeriveEta;
      Functions[D20] = C_Q_UL8S_2D_DeriveXiXi;
      Functions[D11] = C_Q_UL8S_2D_DeriveXiEta;
      Functions[D02] = C_Q_UL8S_2D_DeriveEtaEta;
      PolynomialDegree = 8;
      Accuracy = 8;
      break;
    case BF_C_Q_UL9S_2D:
      Dimension = 87;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_UL9S_2D_Funct;
      Functions[D10] = C_Q_UL9S_2D_DeriveXi;
      Functions[D01] = C_Q_UL9S_2D_DeriveEta;
      Functions[D20] = C_Q_UL9S_2D_DeriveXiXi;
      Functions[D11] = C_Q_UL9S_2D_DeriveXiEta;
      Functions[D02] = C_Q_UL9S_2D_DeriveEtaEta;
      PolynomialDegree = 9;
      Accuracy = 9;
      break;
    case BF_C_Q_UL2SE_2D:
      Dimension = 8;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_UL2SE_2D_Funct;
      Functions[D10] = C_Q_UL2SE_2D_DeriveXi;
      Functions[D01] = C_Q_UL2SE_2D_DeriveEta;
      Functions[D20] = C_Q_UL2SE_2D_DeriveXiXi;
      Functions[D11] = C_Q_UL2SE_2D_DeriveXiEta;
      Functions[D02] = C_Q_UL2SE_2D_DeriveEtaEta;
      PolynomialDegree = 2;
      Accuracy = 2;
      break;
    case BF_C_Q_UL3SE_2D:
      Dimension = 13;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_UL3SE_2D_Funct;
      Functions[D10] = C_Q_UL3SE_2D_DeriveXi;
      Functions[D01] = C_Q_UL3SE_2D_DeriveEta;
      Functions[D20] = C_Q_UL3SE_2D_DeriveXiXi;
      Functions[D11] = C_Q_UL3SE_2D_DeriveXiEta;
      Functions[D02] = C_Q_UL3SE_2D_DeriveEtaEta;
      PolynomialDegree = 3;
      Accuracy = 3;
      break;
    case BF_C_Q_UL4SE_2D:
      Dimension = 20;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_UL4SE_2D_Funct;
      Functions[D10] = C_Q_UL4SE_2D_DeriveXi;
      Functions[D01] = C_Q_UL4SE_2D_DeriveEta;
      Functions[D20] = C_Q_UL4SE_2D_DeriveXiXi;
      Functions[D11] = C_Q_UL4SE_2D_DeriveXiEta;
      Functions[D02] = C_Q_UL4SE_2D_DeriveEtaEta;
      PolynomialDegree = 4;
      Accuracy = 4;
      break;
    case BF_C_Q_UL5SE_2D:
      Dimension = 29;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_UL5SE_2D_Funct;
      Functions[D10] = C_Q_UL5SE_2D_DeriveXi;
      Functions[D01] = C_Q_UL5SE_2D_DeriveEta;
      Functions[D20] = C_Q_UL5SE_2D_DeriveXiXi;
      Functions[D11] = C_Q_UL5SE_2D_DeriveXiEta;
      Functions[D02] = C_Q_UL5SE_2D_DeriveEtaEta;
      PolynomialDegree = 5;
      Accuracy = 5;
      break;
    case BF_C_Q_UL6SE_2D:
      Dimension = 40;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_UL6SE_2D_Funct;
      Functions[D10] = C_Q_UL6SE_2D_DeriveXi;
      Functions[D01] = C_Q_UL6SE_2D_DeriveEta;
      Functions[D20] = C_Q_UL6SE_2D_DeriveXiXi;
      Functions[D11] = C_Q_UL6SE_2D_DeriveXiEta;
      Functions[D02] = C_Q_UL6SE_2D_DeriveEtaEta;
      PolynomialDegree = 6;
      Accuracy = 6;
      break;
    case BF_C_Q_UL7SE_2D:
      Dimension = 53;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_UL7SE_2D_Funct;
      Functions[D10] = C_Q_UL7SE_2D_DeriveXi;
      Functions[D01] = C_Q_UL7SE_2D_DeriveEta;
      Functions[D20] = C_Q_UL7SE_2D_DeriveXiXi;
      Functions[D11] = C_Q_UL7SE_2D_DeriveXiEta;
      Functions[D02] = C_Q_UL7SE_2D_DeriveEtaEta;
      PolynomialDegree = 7;
      Accuracy = 7;
      break;
    case BF_C_Q_UL8SE_2D:
      Dimension = 68;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_UL8SE_2D_Funct;
      Functions[D10] = C_Q_UL8SE_2D_DeriveXi;
      Functions[D01] = C_Q_UL8SE_2D_DeriveEta;
      Functions[D20] = C_Q_UL8SE_2D_DeriveXiXi;
      Functions[D11] = C_Q_UL8SE_2D_DeriveXiEta;
      Functions[D02] = C_Q_UL8SE_2D_DeriveEtaEta;
      PolynomialDegree = 8;
      Accuracy = 8;
      break;
    case BF_C_Q_UL9SE_2D:
      Dimension = 85;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_UL9SE_2D_Funct;
      Functions[D10] = C_Q_UL9SE_2D_DeriveXi;
      Functions[D01] = C_Q_UL9SE_2D_DeriveEta;
      Functions[D20] = C_Q_UL9SE_2D_DeriveXiXi;
      Functions[D11] = C_Q_UL9SE_2D_DeriveXiEta;
      Functions[D02] = C_Q_UL9SE_2D_DeriveEtaEta;
      PolynomialDegree = 9;
      Accuracy = 9;
      break;
    case BF_C_Q_M2_2D:
      Dimension = 8;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_M2_2D_Funct;
      Functions[D10] = C_Q_M2_2D_DeriveXi;
      Functions[D01] = C_Q_M2_2D_DeriveEta;
      Functions[D20] = C_Q_M2_2D_DeriveXiXi;
      Functions[D11] = C_Q_M2_2D_DeriveXiEta;
      Functions[D02] = C_Q_M2_2D_DeriveEtaEta;
      PolynomialDegree = 2;
      Accuracy = 2;
      break;
    case BF_C_Q_M3_2D:
      Dimension = 12;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_M3_2D_Funct;
      Functions[D10] = C_Q_M3_2D_DeriveXi;
      Functions[D01] = C_Q_M3_2D_DeriveEta;
      Functions[D20] = C_Q_M3_2D_DeriveXiXi;
      Functions[D11] = C_Q_M3_2D_DeriveXiEta;
      Functions[D02] = C_Q_M3_2D_DeriveEtaEta;
      PolynomialDegree = 3;
      Accuracy = 3;
      break;
    case BF_C_Q_M4_2D:
      Dimension = 17;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_M4_2D_Funct;
      Functions[D10] = C_Q_M4_2D_DeriveXi;
      Functions[D01] = C_Q_M4_2D_DeriveEta;
      Functions[D20] = C_Q_M4_2D_DeriveXiXi;
      Functions[D11] = C_Q_M4_2D_DeriveXiEta;
      Functions[D02] = C_Q_M4_2D_DeriveEtaEta;
      PolynomialDegree = 4;
      Accuracy = 4;
      break;
    case BF_C_Q_M5_2D:
      Dimension = 23;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_M5_2D_Funct;
      Functions[D10] = C_Q_M5_2D_DeriveXi;
      Functions[D01] = C_Q_M5_2D_DeriveEta;
      Functions[D20] = C_Q_M5_2D_DeriveXiXi;
      Functions[D11] = C_Q_M5_2D_DeriveXiEta;
      Functions[D02] = C_Q_M5_2D_DeriveEtaEta;
      PolynomialDegree = 5;
      Accuracy = 5;
      break;
    case BF_C_Q_M6_2D:
      Dimension = 30;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_M6_2D_Funct;
      Functions[D10] = C_Q_M6_2D_DeriveXi;
      Functions[D01] = C_Q_M6_2D_DeriveEta;
      Functions[D20] = C_Q_M6_2D_DeriveXiXi;
      Functions[D11] = C_Q_M6_2D_DeriveXiEta;
      Functions[D02] = C_Q_M6_2D_DeriveEtaEta;
      PolynomialDegree = 6;
      Accuracy = 6;
      break;
    case BF_C_Q_M7_2D:
      Dimension = 38;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_M7_2D_Funct;
      Functions[D10] = C_Q_M7_2D_DeriveXi;
      Functions[D01] = C_Q_M7_2D_DeriveEta;
      Functions[D20] = C_Q_M7_2D_DeriveXiXi;
      Functions[D11] = C_Q_M7_2D_DeriveXiEta;
      Functions[D02] = C_Q_M7_2D_DeriveEtaEta;
      PolynomialDegree = 7;
      Accuracy = 7;
      break;
    case BF_C_Q_M8_2D:
      Dimension = 47;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_M8_2D_Funct;
      Functions[D10] = C_Q_M8_2D_DeriveXi;
      Functions[D01] = C_Q_M8_2D_DeriveEta;
      Functions[D20] = C_Q_M8_2D_DeriveXiXi;
      Functions[D11] = C_Q_M8_2D_DeriveXiEta;
      Functions[D02] = C_Q_M8_2D_DeriveEtaEta;
      PolynomialDegree = 8;
      Accuracy = 8;
      break;
    case BF_C_Q_M9_2D:
      Dimension = 57;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_M9_2D_Funct;
      Functions[D10] = C_Q_M9_2D_DeriveXi;
      Functions[D01] = C_Q_M9_2D_DeriveEta;
      Functions[D20] = C_Q_M9_2D_DeriveXiXi;
      Functions[D11] = C_Q_M9_2D_DeriveXiEta;
      Functions[D02] = C_Q_M9_2D_DeriveEtaEta;
      PolynomialDegree = 9;
      Accuracy = 9;
      break;
    case BF_C_Q_EL1_2D:
      Dimension = 5;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = C_Q_EL1_2D_Funct;
      Functions[D10] = C_Q_EL1_2D_DeriveXi;
      Functions[D01] = C_Q_EL1_2D_DeriveEta;
      Functions[D20] = C_Q_EL1_2D_DeriveXiXi;
      Functions[D11] = C_Q_EL1_2D_DeriveXiEta;
      Functions[D02] = C_Q_EL1_2D_DeriveEtaEta;
      PolynomialDegree = 2;
      Accuracy = 1;
      // note: this element is the only one which used to set a member called
      // 'spaceDeptBasis' to true. This member was then never used.
      break;
    case BF_N_Q_RT0_2D:
      Dimension = 4;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = N_Q_RT0_2D_Funct;
      Functions[D10] = N_Q_RT0_2D_DeriveXi;
      Functions[D01] = N_Q_RT0_2D_DeriveEta;
      Functions[D20] = N_Q_RT0_2D_DeriveXiXi;
      Functions[D11] = N_Q_RT0_2D_DeriveXiEta;
      Functions[D02] = N_Q_RT0_2D_DeriveEtaEta;
      PolynomialDegree = 2;
      Accuracy = 1;
      N_BF2Change = 1;
      BF2Change = N_Q_RT0_2D_Change;
      BaseVectDim = 2;
      break;
    case BF_N_Q_RT1_2D:
      Dimension = 12;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = N_Q_RT1_2D_Funct;
      Functions[D10] = N_Q_RT1_2D_DeriveXi;
      Functions[D01] = N_Q_RT1_2D_DeriveEta;
      Functions[D20] = N_Q_RT1_2D_DeriveXiXi;
      Functions[D11] = N_Q_RT1_2D_DeriveXiEta;
      Functions[D02] = N_Q_RT1_2D_DeriveEtaEta;
      PolynomialDegree = 3;
      Accuracy = 2;
      N_BF2Change = 1;
      BF2Change = N_Q_RT1_2D_Change;
      BaseVectDim = 2;
      break;
    case BF_N_Q_RT2_2D:
      Dimension = 24;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = N_Q_RT2_2D_Funct;
      Functions[D10] = N_Q_RT2_2D_DeriveXi;
      Functions[D01] = N_Q_RT2_2D_DeriveEta;
      Functions[D20] = N_Q_RT2_2D_DeriveXiXi;
      Functions[D11] = N_Q_RT2_2D_DeriveXiEta;
      Functions[D02] = N_Q_RT2_2D_DeriveEtaEta;
      PolynomialDegree = 4;
      Accuracy = 3;
      N_BF2Change = 2;
      BF2Change = N_Q_RT2_2D_Change;
      BaseVectDim = 2;
      break;
    case BF_N_Q_RT3_2D:
      Dimension = 40;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = N_Q_RT3_2D_Funct;
      Functions[D10] = N_Q_RT3_2D_DeriveXi;
      Functions[D01] = N_Q_RT3_2D_DeriveEta;
      Functions[D20] = N_Q_RT3_2D_DeriveXiXi;
      Functions[D11] = N_Q_RT3_2D_DeriveXiEta;
      Functions[D02] = N_Q_RT3_2D_DeriveEtaEta;
      PolynomialDegree = 7;
      Accuracy = 4;
      N_BF2Change = 2;
      BF2Change = N_Q_RT3_2D_Change;
      BaseVectDim = 2;
      break;
    case BF_N_T_RT0_2D:
      Dimension = 3;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = N_T_RT0_2D_Funct;
      Functions[D10] = N_T_RT0_2D_DeriveXi;
      Functions[D01] = N_T_RT0_2D_DeriveEta;
      Functions[D20] = N_T_RT0_2D_DeriveXiXi;
      Functions[D11] = N_T_RT0_2D_DeriveXiEta;
      Functions[D02] = N_T_RT0_2D_DeriveEtaEta;
      PolynomialDegree = 1;
      Accuracy = 1;
      N_BF2Change = 1;
      BF2Change = N_T_RT0_2D_Change;
      BaseVectDim = 2;
      break;
    case BF_N_T_RT1_2D:
      Dimension = 8;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = N_T_RT1_2D_Funct;
      Functions[D10] = N_T_RT1_2D_DeriveXi;
      Functions[D01] = N_T_RT1_2D_DeriveEta;
      Functions[D20] = N_T_RT1_2D_DeriveXiXi;
      Functions[D11] = N_T_RT1_2D_DeriveXiEta;
      Functions[D02] = N_T_RT1_2D_DeriveEtaEta;
      PolynomialDegree = 2;
      Accuracy = 1;
      N_BF2Change = 1;
      BF2Change = N_T_RT1_2D_Change;
      BaseVectDim = 2;
      break;
    case BF_N_T_RT2_2D:
      Dimension = 15;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = N_T_RT2_2D_Funct;
      Functions[D10] = N_T_RT2_2D_DeriveXi;
      Functions[D01] = N_T_RT2_2D_DeriveEta;
      Functions[D20] = N_T_RT2_2D_DeriveXiXi;
      Functions[D11] = N_T_RT2_2D_DeriveXiEta;
      Functions[D02] = N_T_RT2_2D_DeriveEtaEta;
      PolynomialDegree = 3;
      Accuracy = 2;
      N_BF2Change = 2;
      BF2Change = N_T_RT2_2D_Change;
      BaseVectDim = 2;
      break;
    case BF_N_T_RT3_2D:
      Dimension = 24;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = N_T_RT3_2D_Funct;
      Functions[D10] = N_T_RT3_2D_DeriveXi;
      Functions[D01] = N_T_RT3_2D_DeriveEta;
      Functions[D20] = N_T_RT3_2D_DeriveXiXi;
      Functions[D11] = N_T_RT3_2D_DeriveXiEta;
      Functions[D02] = N_T_RT3_2D_DeriveEtaEta;
      PolynomialDegree = 4;
      Accuracy = 3;
      N_BF2Change = 2;
      BF2Change = N_T_RT3_2D_Change;
      BaseVectDim = 2;
      break;
    case BF_N_Q_BDM1_2D:
      Dimension = 8;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = N_Q_BDM1_2D_Funct;
      Functions[D10] = N_Q_BDM1_2D_DeriveXi;
      Functions[D01] = N_Q_BDM1_2D_DeriveEta;
      Functions[D20] = N_Q_BDM1_2D_DeriveXiXi;
      Functions[D11] = N_Q_BDM1_2D_DeriveXiEta;
      Functions[D02] = N_Q_BDM1_2D_DeriveEtaEta;
      PolynomialDegree = 6;
      Accuracy = 2;
      N_BF2Change = 1;
      BF2Change = N_Q_BDM1_2D_Change;
      BaseVectDim = 2;
      break;
    case BF_N_Q_BDM2_2D:
      Dimension = 14;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = N_Q_BDM2_2D_Funct;
      Functions[D10] = N_Q_BDM2_2D_DeriveXi;
      Functions[D01] = N_Q_BDM2_2D_DeriveEta;
      Functions[D20] = N_Q_BDM2_2D_DeriveXiXi;
      Functions[D11] = N_Q_BDM2_2D_DeriveXiEta;
      Functions[D02] = N_Q_BDM2_2D_DeriveEtaEta;
      PolynomialDegree = 3;
      Accuracy = 3;
      N_BF2Change = 2;
      BF2Change = N_Q_BDM2_2D_Change;
      BaseVectDim = 2;
      break;
    case BF_N_Q_BDM3_2D:
      Dimension = 22;
      RefElement = BFRefElements::BFUnitSquare;
      Functions[D00] = N_Q_BDM3_2D_Funct;
      Functions[D10] = N_Q_BDM3_2D_DeriveXi;
      Functions[D01] = N_Q_BDM3_2D_DeriveEta;
      Functions[D20] = N_Q_BDM3_2D_DeriveXiXi;
      Functions[D11] = N_Q_BDM3_2D_DeriveXiEta;
      Functions[D02] = N_Q_BDM3_2D_DeriveEtaEta;
      PolynomialDegree = 4;
      Accuracy = 4;
      N_BF2Change = 2;
      BF2Change = N_Q_BDM3_2D_Change;
      BaseVectDim = 2;
      break;
    case BF_N_T_BDM1_2D:
      Dimension = 6;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = N_T_BDM1_2D_Funct;
      Functions[D10] = N_T_BDM1_2D_DeriveXi;
      Functions[D01] = N_T_BDM1_2D_DeriveEta;
      Functions[D20] = N_T_BDM1_2D_DeriveXiXi;
      Functions[D11] = N_T_BDM1_2D_DeriveXiEta;
      Functions[D02] = N_T_BDM1_2D_DeriveEtaEta;
      PolynomialDegree = 2;
      Accuracy = 1;
      N_BF2Change = 1;
      BF2Change = N_T_BDM1_2D_Change;
      BaseVectDim = 2;
      break;
    case BF_N_T_BDM2_2D:
      Dimension = 12;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = N_T_BDM2_2D_Funct;
      Functions[D10] = N_T_BDM2_2D_DeriveXi;
      Functions[D01] = N_T_BDM2_2D_DeriveEta;
      Functions[D20] = N_T_BDM2_2D_DeriveXiXi;
      Functions[D11] = N_T_BDM2_2D_DeriveXiEta;
      Functions[D02] = N_T_BDM2_2D_DeriveEtaEta;
      PolynomialDegree = 2;
      Accuracy = 2;
      N_BF2Change = 2;
      BF2Change = N_T_BDM2_2D_Change;
      BaseVectDim = 2;
      break;
    case BF_N_T_BDM3_2D:
      Dimension = 20;
      RefElement = BFRefElements::BFUnitTriangle;
      Functions[D00] = N_T_BDM3_2D_Funct;
      Functions[D10] = N_T_BDM3_2D_DeriveXi;
      Functions[D01] = N_T_BDM3_2D_DeriveEta;
      Functions[D20] = N_T_BDM3_2D_DeriveXiXi;
      Functions[D11] = N_T_BDM3_2D_DeriveXiEta;
      Functions[D02] = N_T_BDM3_2D_DeriveEtaEta;
      PolynomialDegree = 3;
      Accuracy = 3;
      N_BF2Change = 2;
      BF2Change = N_T_BDM3_2D_Change;
      BaseVectDim = 2;
      break;
    case BF_C_T_P00_3D:
      Dimension = 1;
      RefElement = BFRefElements::BFUnitTetrahedron;
      Functions3D[D000] = C_T_P00_3D_Funct;
      Functions3D[D100] = C_T_P00_3D_DeriveXi;
      Functions3D[D010] = C_T_P00_3D_DeriveEta;
      Functions3D[D001] = C_T_P00_3D_DeriveZeta;
      Functions3D[D200] = C_T_P00_3D_DeriveXiXi;
      Functions3D[D110] = C_T_P00_3D_DeriveXiEta;
      Functions3D[D101] = C_T_P00_3D_DeriveXiZeta;
      Functions3D[D020] = C_T_P00_3D_DeriveEtaEta;
      Functions3D[D011] = C_T_P00_3D_DeriveEtaZeta;
      Functions3D[D002] = C_T_P00_3D_DeriveZetaZeta;
      PolynomialDegree = 1;
      Accuracy = 1;
      break;
    case BF_C_T_P0_3D:
      Dimension = 1;
      RefElement = BFRefElements::BFUnitTetrahedron;
      Functions3D[D000] = C_T_P0_3D_Funct;
      Functions3D[D100] = C_T_P0_3D_DeriveXi;
      Functions3D[D010] = C_T_P0_3D_DeriveEta;
      Functions3D[D001] = C_T_P0_3D_DeriveZeta;
      Functions3D[D200] = C_T_P0_3D_DeriveXiXi;
      Functions3D[D110] = C_T_P0_3D_DeriveXiEta;
      Functions3D[D101] = C_T_P0_3D_DeriveXiZeta;
      Functions3D[D020] = C_T_P0_3D_DeriveEtaEta;
      Functions3D[D011] = C_T_P0_3D_DeriveEtaZeta;
      Functions3D[D002] = C_T_P0_3D_DeriveZetaZeta;
      PolynomialDegree = 0;
      Accuracy = 0;
      break;
    case BF_C_T_P1_3D:
      Dimension = 4;
      RefElement = BFRefElements::BFUnitTetrahedron;
      Functions3D[D000] = C_T_P1_3D_Funct;
      Functions3D[D100] = C_T_P1_3D_DeriveXi;
      Functions3D[D010] = C_T_P1_3D_DeriveEta;
      Functions3D[D001] = C_T_P1_3D_DeriveZeta;
      Functions3D[D200] = C_T_P1_3D_DeriveXiXi;
      Functions3D[D110] = C_T_P1_3D_DeriveXiEta;
      Functions3D[D101] = C_T_P1_3D_DeriveXiZeta;
      Functions3D[D020] = C_T_P1_3D_DeriveEtaEta;
      Functions3D[D011] = C_T_P1_3D_DeriveEtaZeta;
      Functions3D[D002] = C_T_P1_3D_DeriveZetaZeta;
      PolynomialDegree = 1;
      Accuracy = 1;
      break;
    case BF_C_T_P2_3D:
      Dimension = 10;
      RefElement = BFRefElements::BFUnitTetrahedron;
      Functions3D[D000] = C_T_P2_3D_Funct;
      Functions3D[D100] = C_T_P2_3D_DeriveXi;
      Functions3D[D010] = C_T_P2_3D_DeriveEta;
      Functions3D[D001] = C_T_P2_3D_DeriveZeta;
      Functions3D[D200] = C_T_P2_3D_DeriveXiXi;
      Functions3D[D110] = C_T_P2_3D_DeriveXiEta;
      Functions3D[D101] = C_T_P2_3D_DeriveXiZeta;
      Functions3D[D020] = C_T_P2_3D_DeriveEtaEta;
      Functions3D[D011] = C_T_P2_3D_DeriveEtaZeta;
      Functions3D[D002] = C_T_P2_3D_DeriveZetaZeta;
      PolynomialDegree = 2;
      Accuracy = 2;
      break;
    case BF_C_T_P3_3D:
      Dimension = 20;
      RefElement = BFRefElements::BFUnitTetrahedron;
      Functions3D[D000] = C_T_P3_3D_Funct;
      Functions3D[D100] = C_T_P3_3D_DeriveXi;
      Functions3D[D010] = C_T_P3_3D_DeriveEta;
      Functions3D[D001] = C_T_P3_3D_DeriveZeta;
      Functions3D[D200] = C_T_P3_3D_DeriveXiXi;
      Functions3D[D110] = C_T_P3_3D_DeriveXiEta;
      Functions3D[D101] = C_T_P3_3D_DeriveXiZeta;
      Functions3D[D020] = C_T_P3_3D_DeriveEtaEta;
      Functions3D[D011] = C_T_P3_3D_DeriveEtaZeta;
      Functions3D[D002] = C_T_P3_3D_DeriveZetaZeta;
      PolynomialDegree = 3;
      Accuracy = 3;
      break;
    case BF_N_T_P1_3D:
      Dimension = 4;
      RefElement = BFRefElements::BFUnitTetrahedron;
      Functions3D[D000] = N_T_P1_3D_Funct;
      Functions3D[D100] = N_T_P1_3D_DeriveXi;
      Functions3D[D010] = N_T_P1_3D_DeriveEta;
      Functions3D[D001] = N_T_P1_3D_DeriveZeta;
      Functions3D[D200] = N_T_P1_3D_DeriveXiXi;
      Functions3D[D110] = N_T_P1_3D_DeriveXiEta;
      Functions3D[D101] = N_T_P1_3D_DeriveXiZeta;
      Functions3D[D020] = N_T_P1_3D_DeriveEtaEta;
      Functions3D[D011] = N_T_P1_3D_DeriveEtaZeta;
      Functions3D[D002] = N_T_P1_3D_DeriveZetaZeta;
      PolynomialDegree = 1;
      Accuracy = 1;
      break;
    case BF_C_H_Q00_3D:
      Dimension = 1;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = C_H_Q00_3D_Funct;
      Functions3D[D100] = C_H_Q00_3D_DeriveXi;
      Functions3D[D010] = C_H_Q00_3D_DeriveEta;
      Functions3D[D001] = C_H_Q00_3D_DeriveZeta;
      Functions3D[D200] = C_H_Q00_3D_DeriveXiXi;
      Functions3D[D110] = C_H_Q00_3D_DeriveXiEta;
      Functions3D[D101] = C_H_Q00_3D_DeriveXiZeta;
      Functions3D[D020] = C_H_Q00_3D_DeriveEtaEta;
      Functions3D[D011] = C_H_Q00_3D_DeriveEtaZeta;
      Functions3D[D002] = C_H_Q00_3D_DeriveZetaZeta;
      PolynomialDegree = 1;
      Accuracy = 1;
      break;
    case BF_C_H_Q0_3D:
      Dimension = 1;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = C_H_Q0_3D_Funct;
      Functions3D[D100] = C_H_Q0_3D_DeriveXi;
      Functions3D[D010] = C_H_Q0_3D_DeriveEta;
      Functions3D[D001] = C_H_Q0_3D_DeriveZeta;
      Functions3D[D200] = C_H_Q0_3D_DeriveXiXi;
      Functions3D[D110] = C_H_Q0_3D_DeriveXiEta;
      Functions3D[D101] = C_H_Q0_3D_DeriveXiZeta;
      Functions3D[D020] = C_H_Q0_3D_DeriveEtaEta;
      Functions3D[D011] = C_H_Q0_3D_DeriveEtaZeta;
      Functions3D[D002] = C_H_Q0_3D_DeriveZetaZeta;
      PolynomialDegree = 1;
      Accuracy = 1;
      break;
    case BF_C_H_Q1_3D:
      Dimension = 8;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = C_H_Q1_3D_Funct;
      Functions3D[D100] = C_H_Q1_3D_DeriveXi;
      Functions3D[D010] = C_H_Q1_3D_DeriveEta;
      Functions3D[D001] = C_H_Q1_3D_DeriveZeta;
      Functions3D[D200] = C_H_Q1_3D_DeriveXiXi;
      Functions3D[D110] = C_H_Q1_3D_DeriveXiEta;
      Functions3D[D101] = C_H_Q1_3D_DeriveXiZeta;
      Functions3D[D020] = C_H_Q1_3D_DeriveEtaEta;
      Functions3D[D011] = C_H_Q1_3D_DeriveEtaZeta;
      Functions3D[D002] = C_H_Q1_3D_DeriveZetaZeta;
      PolynomialDegree = 1;
      Accuracy = 1;
      break;
    case BF_C_H_Q2_3D:
      Dimension = 27;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = C_H_Q2_3D_Funct;
      Functions3D[D100] = C_H_Q2_3D_DeriveXi;
      Functions3D[D010] = C_H_Q2_3D_DeriveEta;
      Functions3D[D001] = C_H_Q2_3D_DeriveZeta;
      Functions3D[D200] = C_H_Q2_3D_DeriveXiXi;
      Functions3D[D110] = C_H_Q2_3D_DeriveXiEta;
      Functions3D[D101] = C_H_Q2_3D_DeriveXiZeta;
      Functions3D[D020] = C_H_Q2_3D_DeriveEtaEta;
      Functions3D[D011] = C_H_Q2_3D_DeriveEtaZeta;
      Functions3D[D002] = C_H_Q2_3D_DeriveZetaZeta;
      PolynomialDegree = 2;
      Accuracy = 2;
      break;
    case BF_C_H_Q3_3D:
      Dimension = 64;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = C_H_Q3_3D_Funct;
      Functions3D[D100] = C_H_Q3_3D_DeriveXi;
      Functions3D[D010] = C_H_Q3_3D_DeriveEta;
      Functions3D[D001] = C_H_Q3_3D_DeriveZeta;
      Functions3D[D200] = C_H_Q3_3D_DeriveXiXi;
      Functions3D[D110] = C_H_Q3_3D_DeriveXiEta;
      Functions3D[D101] = C_H_Q3_3D_DeriveXiZeta;
      Functions3D[D020] = C_H_Q3_3D_DeriveEtaEta;
      Functions3D[D011] = C_H_Q3_3D_DeriveEtaZeta;
      Functions3D[D002] = C_H_Q3_3D_DeriveZetaZeta;
      PolynomialDegree = 3;
      Accuracy = 3;
      break;
    case BF_C_H_Q4_3D:
      Dimension = 125;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = C_H_Q4_3D_Funct;
      Functions3D[D100] = C_H_Q4_3D_DeriveXi;
      Functions3D[D010] = C_H_Q4_3D_DeriveEta;
      Functions3D[D001] = C_H_Q4_3D_DeriveZeta;
      Functions3D[D200] = C_H_Q4_3D_DeriveXiXi;
      Functions3D[D110] = C_H_Q4_3D_DeriveXiEta;
      Functions3D[D101] = C_H_Q4_3D_DeriveXiZeta;
      Functions3D[D020] = C_H_Q4_3D_DeriveEtaEta;
      Functions3D[D011] = C_H_Q4_3D_DeriveEtaZeta;
      Functions3D[D002] = C_H_Q4_3D_DeriveZetaZeta;
      PolynomialDegree = 4;
      Accuracy = 4;
      break;
    case BF_N_H_Q1_3D:
      Dimension = 6;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = N_H_Q1_3D_Funct;
      Functions3D[D100] = N_H_Q1_3D_DeriveXi;
      Functions3D[D010] = N_H_Q1_3D_DeriveEta;
      Functions3D[D001] = N_H_Q1_3D_DeriveZeta;
      Functions3D[D200] = N_H_Q1_3D_DeriveXiXi;
      Functions3D[D110] = N_H_Q1_3D_DeriveXiEta;
      Functions3D[D101] = N_H_Q1_3D_DeriveXiZeta;
      Functions3D[D020] = N_H_Q1_3D_DeriveEtaEta;
      Functions3D[D011] = N_H_Q1_3D_DeriveEtaZeta;
      Functions3D[D002] = N_H_Q1_3D_DeriveZetaZeta;
      PolynomialDegree = 2;
      Accuracy = 1;
      break;
    case BF_D_H_P1_3D:
      Dimension = 4;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = D_H_P1_3D_Funct;
      Functions3D[D100] = D_H_P1_3D_DeriveXi;
      Functions3D[D010] = D_H_P1_3D_DeriveEta;
      Functions3D[D001] = D_H_P1_3D_DeriveZeta;
      Functions3D[D200] = D_H_P1_3D_DeriveXiXi;
      Functions3D[D110] = D_H_P1_3D_DeriveXiEta;
      Functions3D[D101] = D_H_P1_3D_DeriveXiZeta;
      Functions3D[D020] = D_H_P1_3D_DeriveEtaEta;
      Functions3D[D011] = D_H_P1_3D_DeriveEtaZeta;
      Functions3D[D002] = D_H_P1_3D_DeriveZetaZeta;
      PolynomialDegree = 1;
      Accuracy = 1;
      break;
    case BF_D_H_P2_3D:
      Dimension = 10;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = D_H_P2_3D_Funct;
      Functions3D[D100] = D_H_P2_3D_DeriveXi;
      Functions3D[D010] = D_H_P2_3D_DeriveEta;
      Functions3D[D001] = D_H_P2_3D_DeriveZeta;
      Functions3D[D200] = D_H_P2_3D_DeriveXiXi;
      Functions3D[D110] = D_H_P2_3D_DeriveXiEta;
      Functions3D[D101] = D_H_P2_3D_DeriveXiZeta;
      Functions3D[D020] = D_H_P2_3D_DeriveEtaEta;
      Functions3D[D011] = D_H_P2_3D_DeriveEtaZeta;
      Functions3D[D002] = D_H_P2_3D_DeriveZetaZeta;
      PolynomialDegree = 2;
      Accuracy = 2;
      break;
    case BF_D_H_P3_3D:
      Dimension = 20;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = D_H_P3_3D_Funct;
      Functions3D[D100] = D_H_P3_3D_DeriveXi;
      Functions3D[D010] = D_H_P3_3D_DeriveEta;
      Functions3D[D001] = D_H_P3_3D_DeriveZeta;
      Functions3D[D200] = D_H_P3_3D_DeriveXiXi;
      Functions3D[D110] = D_H_P3_3D_DeriveXiEta;
      Functions3D[D101] = D_H_P3_3D_DeriveXiZeta;
      Functions3D[D020] = D_H_P3_3D_DeriveEtaEta;
      Functions3D[D011] = D_H_P3_3D_DeriveEtaZeta;
      Functions3D[D002] = D_H_P3_3D_DeriveZetaZeta;
      PolynomialDegree = 3;
      Accuracy = 3;
      break;
    case BF_D_H_Q1_3D:
      Dimension = 8;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = D_H_Q1_3D_Funct;
      Functions3D[D100] = D_H_Q1_3D_DeriveXi;
      Functions3D[D010] = D_H_Q1_3D_DeriveEta;
      Functions3D[D001] = D_H_Q1_3D_DeriveZeta;
      Functions3D[D200] = D_H_Q1_3D_DeriveXiXi;
      Functions3D[D110] = D_H_Q1_3D_DeriveXiEta;
      Functions3D[D101] = D_H_Q1_3D_DeriveXiZeta;
      Functions3D[D020] = D_H_Q1_3D_DeriveEtaEta;
      Functions3D[D011] = D_H_Q1_3D_DeriveEtaZeta;
      Functions3D[D002] = D_H_Q1_3D_DeriveZetaZeta;
      PolynomialDegree = 3;
      Accuracy = 1;
      break;
    case BF_D_H_Q2_3D:
      Dimension = 27;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = D_H_Q2_3D_Funct;
      Functions3D[D100] = D_H_Q2_3D_DeriveXi;
      Functions3D[D010] = D_H_Q2_3D_DeriveEta;
      Functions3D[D001] = D_H_Q2_3D_DeriveZeta;
      Functions3D[D200] = D_H_Q2_3D_DeriveXiXi;
      Functions3D[D110] = D_H_Q2_3D_DeriveXiEta;
      Functions3D[D101] = D_H_Q2_3D_DeriveXiZeta;
      Functions3D[D020] = D_H_Q2_3D_DeriveEtaEta;
      Functions3D[D011] = D_H_Q2_3D_DeriveEtaZeta;
      Functions3D[D002] = D_H_Q2_3D_DeriveZetaZeta;
      PolynomialDegree = 3;
      Accuracy = 1;
      break;
    case BF_B_H_IB2_3D:
      Dimension = 1;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = B_H_IB2_3D_Funct;
      Functions3D[D100] = B_H_IB2_3D_DeriveXi;
      Functions3D[D010] = B_H_IB2_3D_DeriveEta;
      Functions3D[D001] = B_H_IB2_3D_DeriveZeta;
      Functions3D[D200] = B_H_IB2_3D_DeriveXiXi;
      Functions3D[D110] = B_H_IB2_3D_DeriveXiEta;
      Functions3D[D101] = B_H_IB2_3D_DeriveXiZeta;
      Functions3D[D020] = B_H_IB2_3D_DeriveEtaEta;
      Functions3D[D011] = B_H_IB2_3D_DeriveEtaZeta;
      Functions3D[D002] = B_H_IB2_3D_DeriveZetaZeta;
      PolynomialDegree = 2;
      Accuracy = 2;
      break;
    case BF_N_T_P2_3D:
      Dimension = 13;
      RefElement = BFRefElements::BFUnitTetrahedron;
      Functions3D[D000] = N_T_P2_3D_Funct;
      Functions3D[D100] = N_T_P2_3D_DeriveXi;
      Functions3D[D010] = N_T_P2_3D_DeriveEta;
      Functions3D[D001] = N_T_P2_3D_DeriveZeta;
      Functions3D[D200] = N_T_P2_3D_DeriveXiXi;
      Functions3D[D110] = N_T_P2_3D_DeriveXiEta;
      Functions3D[D101] = N_T_P2_3D_DeriveXiZeta;
      Functions3D[D020] = N_T_P2_3D_DeriveEtaEta;
      Functions3D[D011] = N_T_P2_3D_DeriveEtaZeta;
      Functions3D[D002] = N_T_P2_3D_DeriveZetaZeta;
      PolynomialDegree = 3;
      Accuracy = 2;
      break;
    case BF_N_T_P3_3D:
      ErrThrow("Basis Functions of type BF_N_T_P3_3D not implemented");
      break;
    case BF_N_T_P4_3D:
      ErrThrow("Basis Functions of type BF_N_T_P4_3D not implemented");
      break;
    case BF_N_T_P5_3D:
      ErrThrow("Basis Functions of type BF_N_T_P5_3D not implemented");
      break;
    case BF_N_H_Q2_3D:
      Dimension = 19;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = N_H_Q2_3D_Funct;
      Functions3D[D100] = N_H_Q2_3D_DeriveXi;
      Functions3D[D010] = N_H_Q2_3D_DeriveEta;
      Functions3D[D001] = N_H_Q2_3D_DeriveZeta;
      Functions3D[D200] = N_H_Q2_3D_DeriveXiXi;
      Functions3D[D110] = N_H_Q2_3D_DeriveXiEta;
      Functions3D[D101] = N_H_Q2_3D_DeriveXiZeta;
      Functions3D[D020] = N_H_Q2_3D_DeriveEtaEta;
      Functions3D[D011] = N_H_Q2_3D_DeriveEtaZeta;
      Functions3D[D002] = N_H_Q2_3D_DeriveZetaZeta;
      PolynomialDegree = 3;
      Accuracy = 2;
      N_BF2Change = 1;
      BF2Change = N_H_Q2_3D_Change;
      break;
    case BF_N_H_Q3_3D:
      Dimension = 40;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = N_H_Q3_3D_Funct;
      Functions3D[D100] = N_H_Q3_3D_DeriveXi;
      Functions3D[D010] = N_H_Q3_3D_DeriveEta;
      Functions3D[D001] = N_H_Q3_3D_DeriveZeta;
      Functions3D[D200] = N_H_Q3_3D_DeriveXiXi;
      Functions3D[D110] = N_H_Q3_3D_DeriveXiEta;
      Functions3D[D101] = N_H_Q3_3D_DeriveXiZeta;
      Functions3D[D020] = N_H_Q3_3D_DeriveEtaEta;
      Functions3D[D011] = N_H_Q3_3D_DeriveEtaZeta;
      Functions3D[D002] = N_H_Q3_3D_DeriveZetaZeta;
      PolynomialDegree = 4;
      Accuracy = 3;
      N_BF2Change = 2;
      BF2Change = N_H_Q3_3D_Change;
      break;
    case BF_N_H_Q4_3D:
      Dimension = 70;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = N_H_Q4_3D_Funct;
      Functions3D[D100] = N_H_Q4_3D_DeriveXi;
      Functions3D[D010] = N_H_Q4_3D_DeriveEta;
      Functions3D[D001] = N_H_Q4_3D_DeriveZeta;
      Functions3D[D200] = N_H_Q4_3D_DeriveXiXi;
      Functions3D[D110] = N_H_Q4_3D_DeriveXiEta;
      Functions3D[D101] = N_H_Q4_3D_DeriveXiZeta;
      Functions3D[D020] = N_H_Q4_3D_DeriveEtaEta;
      Functions3D[D011] = N_H_Q4_3D_DeriveEtaZeta;
      Functions3D[D002] = N_H_Q4_3D_DeriveZetaZeta;
      PolynomialDegree = 5;
      Accuracy = 4;
      N_BF2Change = 4;
      BF2Change = N_H_Q4_3D_Change;
      break;
    case BF_N_H_Q5_3D:
      ErrThrow("Basis Functions of type BF_N_H_Q5_3D not implemented");
      break;
    case BF_C_T_B2_3D:
      Dimension = 15;
      RefElement = BFRefElements::BFUnitTetrahedron;
      Functions3D[D000] = C_T_B2_3D_Funct;
      Functions3D[D100] = C_T_B2_3D_DeriveXi;
      Functions3D[D010] = C_T_B2_3D_DeriveEta;
      Functions3D[D001] = C_T_B2_3D_DeriveZeta;
      Functions3D[D200] = C_T_B2_3D_DeriveXiXi;
      Functions3D[D110] = C_T_B2_3D_DeriveXiEta;
      Functions3D[D101] = C_T_B2_3D_DeriveXiZeta;
      Functions3D[D020] = C_T_B2_3D_DeriveEtaEta;
      Functions3D[D011] = C_T_B2_3D_DeriveEtaZeta;
      Functions3D[D002] = C_T_B2_3D_DeriveZetaZeta;
      PolynomialDegree = 3;
      Accuracy = 2;
      break;
    case BF_D_T_P1_3D:
      Dimension = 4;
      RefElement = BFRefElements::BFUnitTetrahedron;
      Functions3D[D000] = D_T_P1_3D_Funct;
      Functions3D[D100] = D_T_P1_3D_DeriveXi;
      Functions3D[D010] = D_T_P1_3D_DeriveEta;
      Functions3D[D001] = D_T_P1_3D_DeriveZeta;
      Functions3D[D200] = D_T_P1_3D_DeriveXiXi;
      Functions3D[D110] = D_T_P1_3D_DeriveXiEta;
      Functions3D[D101] = D_T_P1_3D_DeriveXiZeta;
      Functions3D[D020] = D_T_P1_3D_DeriveEtaEta;
      Functions3D[D011] = D_T_P1_3D_DeriveEtaZeta;
      Functions3D[D002] = D_T_P1_3D_DeriveZetaZeta;
      PolynomialDegree = 1;
      Accuracy = 1;
      break;
    case BF_D_T_P2_3D:
      Dimension = 10;
      RefElement = BFRefElements::BFUnitTetrahedron;
      Functions3D[D000] = D_T_P2_3D_Funct;
      Functions3D[D100] = D_T_P2_3D_DeriveXi;
      Functions3D[D010] = D_T_P2_3D_DeriveEta;
      Functions3D[D001] = D_T_P2_3D_DeriveZeta;
      Functions3D[D200] = D_T_P2_3D_DeriveXiXi;
      Functions3D[D110] = D_T_P2_3D_DeriveXiEta;
      Functions3D[D101] = D_T_P2_3D_DeriveXiZeta;
      Functions3D[D020] = D_T_P2_3D_DeriveEtaEta;
      Functions3D[D011] = D_T_P2_3D_DeriveEtaZeta;
      Functions3D[D002] = D_T_P2_3D_DeriveZetaZeta;
      PolynomialDegree = 2;
      Accuracy = 1;
      break;
    case BF_D_T_P3_3D:
      Dimension = 20;
      RefElement = BFRefElements::BFUnitTetrahedron;
      Functions3D[D000] = D_T_P3_3D_Funct;
      Functions3D[D100] = D_T_P3_3D_DeriveXi;
      Functions3D[D010] = D_T_P3_3D_DeriveEta;
      Functions3D[D001] = D_T_P3_3D_DeriveZeta;
      Functions3D[D200] = D_T_P3_3D_DeriveXiXi;
      Functions3D[D110] = D_T_P3_3D_DeriveXiEta;
      Functions3D[D101] = D_T_P3_3D_DeriveXiZeta;
      Functions3D[D020] = D_T_P3_3D_DeriveEtaEta;
      Functions3D[D011] = D_T_P3_3D_DeriveEtaZeta;
      Functions3D[D002] = D_T_P3_3D_DeriveZetaZeta;
      PolynomialDegree = 3;
      Accuracy = 1;
      break;
    case BF_C_H_UL1_3D:
      Dimension = 9;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = C_H_UL1_3D_Funct;
      Functions3D[D100] = C_H_UL1_3D_DeriveXi;
      Functions3D[D010] = C_H_UL1_3D_DeriveEta;
      Functions3D[D001] = C_H_UL1_3D_DeriveZeta;
      Functions3D[D200] = C_H_UL1_3D_DeriveXiXi;
      Functions3D[D110] = C_H_UL1_3D_DeriveXiEta;
      Functions3D[D101] = C_H_UL1_3D_DeriveXiZeta;
      Functions3D[D020] = C_H_UL1_3D_DeriveEtaEta;
      Functions3D[D011] = C_H_UL1_3D_DeriveEtaZeta;
      Functions3D[D002] = C_H_UL1_3D_DeriveZetaZeta;
      PolynomialDegree = 2;
      Accuracy = 1;
      break;
    case BF_C_H_UL2_3D:
      Dimension = 30;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = C_H_UL2_3D_Funct;
      Functions3D[D100] = C_H_UL2_3D_DeriveXi;
      Functions3D[D010] = C_H_UL2_3D_DeriveEta;
      Functions3D[D001] = C_H_UL2_3D_DeriveZeta;
      Functions3D[D200] = C_H_UL2_3D_DeriveXiXi;
      Functions3D[D110] = C_H_UL2_3D_DeriveXiEta;
      Functions3D[D101] = C_H_UL2_3D_DeriveXiZeta;
      Functions3D[D020] = C_H_UL2_3D_DeriveEtaEta;
      Functions3D[D011] = C_H_UL2_3D_DeriveEtaZeta;
      Functions3D[D002] = C_H_UL2_3D_DeriveZetaZeta;
      PolynomialDegree = 3;
      Accuracy = 2;
      break;
    case BF_C_H_UL3_3D:
      Dimension = 67;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = C_H_UL3_3D_Funct;
      Functions3D[D100] = C_H_UL3_3D_DeriveXi;
      Functions3D[D010] = C_H_UL3_3D_DeriveEta;
      Functions3D[D001] = C_H_UL3_3D_DeriveZeta;
      Functions3D[D200] = C_H_UL3_3D_DeriveXiXi;
      Functions3D[D110] = C_H_UL3_3D_DeriveXiEta;
      Functions3D[D101] = C_H_UL3_3D_DeriveXiZeta;
      Functions3D[D020] = C_H_UL3_3D_DeriveEtaEta;
      Functions3D[D011] = C_H_UL3_3D_DeriveEtaZeta;
      Functions3D[D002] = C_H_UL3_3D_DeriveZetaZeta;
      PolynomialDegree = 4;
      Accuracy = 3;
      break;
    case BF_C_H_UL4_3D:
      ErrThrow("Basis Functions of type BF_C_H_UL4_3D not implemented");
      break;
    case BF_C_H_UL5_3D:
      ErrThrow("Basis Functions of type BF_C_H_UL5_3D not implemented");
      break;
    case BF_C_T_UL1_3D:
      ErrThrow("Basis Functions of type BF_C_T_UL1_3D not implemented");
      break;
    case BF_C_T_UL2_3D:
      ErrThrow("Basis Functions of type BF_C_T_UL2_3D not implemented");
      break;
    case BF_C_T_UL3_3D:
      ErrThrow("Basis Functions of type BF_C_T_UL3_3D not implemented");
      break;
    case BF_C_T_UL4_3D:
      ErrThrow("Basis Functions of type BF_C_T_UL4_3D not implemented");
      break;
    case BF_C_T_UL5_3D:
      ErrThrow("Basis Functions of type BF_C_T_UL5_3D not implemented");
      break;
    case BF_N_T_RT0_3D:
      Dimension = 4;
      RefElement = BFRefElements::BFUnitTetrahedron;
      Functions3D[D000] = N_T_RT0_3D_Funct;
      Functions3D[D100] = N_T_RT0_3D_DeriveXi;
      Functions3D[D010] = N_T_RT0_3D_DeriveEta;
      Functions3D[D001] = N_T_RT0_3D_DeriveZeta;
      Functions3D[D200] = N_T_RT0_3D_DeriveXiXi;
      Functions3D[D110] = N_T_RT0_3D_DeriveXiEta;
      Functions3D[D101] = N_T_RT0_3D_DeriveXiZeta;
      Functions3D[D020] = N_T_RT0_3D_DeriveEtaEta;
      Functions3D[D011] = N_T_RT0_3D_DeriveEtaZeta;
      Functions3D[D002] = N_T_RT0_3D_DeriveZetaZeta;
      PolynomialDegree = 1;
      Accuracy = 1;
      N_BF2Change = 1;
      BF2Change = N_T_RT0_3D_Change;
      BaseVectDim = 3;
      break;
    case BF_N_T_RT1_3D:
      Dimension = 15;
      RefElement = BFRefElements::BFUnitTetrahedron;
      Functions3D[D000] = N_T_RT1_3D_Funct;
      Functions3D[D100] = N_T_RT1_3D_DeriveXi;
      Functions3D[D010] = N_T_RT1_3D_DeriveEta;
      Functions3D[D001] = N_T_RT1_3D_DeriveZeta;
      Functions3D[D200] = N_T_RT1_3D_DeriveXiXi;
      Functions3D[D110] = N_T_RT1_3D_DeriveXiEta;
      Functions3D[D101] = N_T_RT1_3D_DeriveXiZeta;
      Functions3D[D020] = N_T_RT1_3D_DeriveEtaEta;
      Functions3D[D011] = N_T_RT1_3D_DeriveEtaZeta;
      Functions3D[D002] = N_T_RT1_3D_DeriveZetaZeta;
      PolynomialDegree = 2;
      Accuracy = 1;
      N_BF2Change = 3;
      BF2Change = N_T_RT1_3D_Change;
      BaseVectDim = 3;
      break;
    case BF_N_T_RT2_3D:
      Dimension = 36;
      RefElement = BFRefElements::BFUnitTetrahedron;
      Functions3D[D000] = N_T_RT2_3D_Funct;
      Functions3D[D100] = N_T_RT2_3D_DeriveXi;
      Functions3D[D010] = N_T_RT2_3D_DeriveEta;
      Functions3D[D001] = N_T_RT2_3D_DeriveZeta;
      Functions3D[D200] = N_T_RT2_3D_DeriveXiXi;
      Functions3D[D110] = N_T_RT2_3D_DeriveXiEta;
      Functions3D[D101] = N_T_RT2_3D_DeriveXiZeta;
      Functions3D[D020] = N_T_RT2_3D_DeriveEtaEta;
      Functions3D[D011] = N_T_RT2_3D_DeriveEtaZeta;
      Functions3D[D002] = N_T_RT2_3D_DeriveZetaZeta;
      PolynomialDegree = 3;
      Accuracy = 1;
      N_BF2Change = 6;
      BF2Change = N_T_RT2_3D_Change;
      BaseVectDim = 3;
      break;
    case BF_N_T_RT3_3D:
      Dimension = 70;
      RefElement = BFRefElements::BFUnitTetrahedron;
      Functions3D[D000] = N_T_RT3_3D_Funct;
      Functions3D[D100] = N_T_RT3_3D_DeriveXi;
      Functions3D[D010] = N_T_RT3_3D_DeriveEta;
      Functions3D[D001] = N_T_RT3_3D_DeriveZeta;
      Functions3D[D200] = N_T_RT3_3D_DeriveXiXi;
      Functions3D[D110] = N_T_RT3_3D_DeriveXiEta;
      Functions3D[D101] = N_T_RT3_3D_DeriveXiZeta;
      Functions3D[D020] = N_T_RT3_3D_DeriveEtaEta;
      Functions3D[D011] = N_T_RT3_3D_DeriveEtaZeta;
      Functions3D[D002] = N_T_RT3_3D_DeriveZetaZeta;
      PolynomialDegree = 4;
      Accuracy = 1;
      N_BF2Change = 10;
      BF2Change = N_T_RT3_3D_Change;
      BaseVectDim = 3;
      break;
    case BF_N_H_RT0_3D:
      Dimension = 6;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = N_H_RT0_3D_Funct;
      Functions3D[D100] = N_H_RT0_3D_DeriveXi;
      Functions3D[D010] = N_H_RT0_3D_DeriveEta;
      Functions3D[D001] = N_H_RT0_3D_DeriveZeta;
      Functions3D[D200] = N_H_RT0_3D_DeriveXiXi;
      Functions3D[D110] = N_H_RT0_3D_DeriveXiEta;
      Functions3D[D101] = N_H_RT0_3D_DeriveXiZeta;
      Functions3D[D020] = N_H_RT0_3D_DeriveEtaEta;
      Functions3D[D011] = N_H_RT0_3D_DeriveEtaZeta;
      Functions3D[D002] = N_H_RT0_3D_DeriveZetaZeta;
      PolynomialDegree = 1;
      Accuracy = 1;
      N_BF2Change = 1;
      BF2Change = N_H_RT0_3D_Change;
      BaseVectDim = 3;
      break;
    case BF_N_H_RT1_3D:
      Dimension = 36;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = N_H_RT1_3D_Funct;
      Functions3D[D100] = N_H_RT1_3D_DeriveXi;
      Functions3D[D010] = N_H_RT1_3D_DeriveEta;
      Functions3D[D001] = N_H_RT1_3D_DeriveZeta;
      Functions3D[D200] = N_H_RT1_3D_DeriveXiXi;
      Functions3D[D110] = N_H_RT1_3D_DeriveXiEta;
      Functions3D[D101] = N_H_RT1_3D_DeriveXiZeta;
      Functions3D[D020] = N_H_RT1_3D_DeriveEtaEta;
      Functions3D[D011] = N_H_RT1_3D_DeriveEtaZeta;
      Functions3D[D002] = N_H_RT1_3D_DeriveZetaZeta;
      PolynomialDegree = 4;
      Accuracy = 1;
      N_BF2Change = 4;
      BF2Change = N_H_RT1_3D_Change;
      BaseVectDim = 3;
      break;
    case BF_N_H_RT2_3D:
      Dimension = 108;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = N_H_RT2_3D_Funct;
      Functions3D[D100] = N_H_RT2_3D_DeriveXi;
      Functions3D[D010] = N_H_RT2_3D_DeriveEta;
      Functions3D[D001] = N_H_RT2_3D_DeriveZeta;
      Functions3D[D200] = N_H_RT2_3D_DeriveXiXi;
      Functions3D[D110] = N_H_RT2_3D_DeriveXiEta;
      Functions3D[D101] = N_H_RT2_3D_DeriveXiZeta;
      Functions3D[D020] = N_H_RT2_3D_DeriveEtaEta;
      Functions3D[D011] = N_H_RT2_3D_DeriveEtaZeta;
      Functions3D[D002] = N_H_RT2_3D_DeriveZetaZeta;
      PolynomialDegree = 7;
      Accuracy = 1;
      N_BF2Change = 9;
      BF2Change = N_H_RT2_3D_Change;
      BaseVectDim = 3;
      break;
    case BF_N_T_BDDF1_3D:
      Dimension = 12;
      RefElement = BFRefElements::BFUnitTetrahedron;
      Functions3D[D000] = N_T_BDDF1_3D_Funct;
      Functions3D[D100] = N_T_BDDF1_3D_DeriveXi;
      Functions3D[D010] = N_T_BDDF1_3D_DeriveEta;
      Functions3D[D001] = N_T_BDDF1_3D_DeriveZeta;
      Functions3D[D200] = N_T_BDDF1_3D_DeriveXiXi;
      Functions3D[D110] = N_T_BDDF1_3D_DeriveXiEta;
      Functions3D[D101] = N_T_BDDF1_3D_DeriveXiZeta;
      Functions3D[D020] = N_T_BDDF1_3D_DeriveEtaEta;
      Functions3D[D011] = N_T_BDDF1_3D_DeriveEtaZeta;
      Functions3D[D002] = N_T_BDDF1_3D_DeriveZetaZeta;
      PolynomialDegree = 1;
      Accuracy = 1;
      N_BF2Change = 3;
      BF2Change = N_T_BDDF1_3D_Change;
      BaseVectDim = 3;
      break;
    case BF_N_T_BDDF2_3D:
      Dimension = 30;
      RefElement = BFRefElements::BFUnitTetrahedron;
      Functions3D[D000] = N_T_BDDF2_3D_Funct;
      Functions3D[D100] = N_T_BDDF2_3D_DeriveXi;
      Functions3D[D010] = N_T_BDDF2_3D_DeriveEta;
      Functions3D[D001] = N_T_BDDF2_3D_DeriveZeta;
      Functions3D[D200] = N_T_BDDF2_3D_DeriveXiXi;
      Functions3D[D110] = N_T_BDDF2_3D_DeriveXiEta;
      Functions3D[D101] = N_T_BDDF2_3D_DeriveXiZeta;
      Functions3D[D020] = N_T_BDDF2_3D_DeriveEtaEta;
      Functions3D[D011] = N_T_BDDF2_3D_DeriveEtaZeta;
      Functions3D[D002] = N_T_BDDF2_3D_DeriveZetaZeta;
      PolynomialDegree = 2;
      Accuracy = 1;
      N_BF2Change = 6;
      BF2Change = N_T_BDDF2_3D_Change;
      BaseVectDim = 3;
      break;
    case BF_N_T_BDDF3_3D:
      Dimension = 60;
      RefElement = BFRefElements::BFUnitTetrahedron;
      Functions3D[D000] = N_T_BDDF3_3D_Funct;
      Functions3D[D100] = N_T_BDDF3_3D_DeriveXi;
      Functions3D[D010] = N_T_BDDF3_3D_DeriveEta;
      Functions3D[D001] = N_T_BDDF3_3D_DeriveZeta;
      Functions3D[D200] = N_T_BDDF3_3D_DeriveXiXi;
      Functions3D[D110] = N_T_BDDF3_3D_DeriveXiEta;
      Functions3D[D101] = N_T_BDDF3_3D_DeriveXiZeta;
      Functions3D[D020] = N_T_BDDF3_3D_DeriveEtaEta;
      Functions3D[D011] = N_T_BDDF3_3D_DeriveEtaZeta;
      Functions3D[D002] = N_T_BDDF3_3D_DeriveZetaZeta;
      //polynomial degree is 3 but somehow 4 is needed in order to work properly
      PolynomialDegree = 4;
      Accuracy = 1;
      N_BF2Change = 10;
      BF2Change = N_T_BDDF3_3D_Change;
      BaseVectDim = 3;
      BaseVectDim = 3;
      break;
    case BF_N_H_BDDF1_3D:
      Dimension = 18;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = N_H_BDDF1_3D_Funct;
      Functions3D[D100] = N_H_BDDF1_3D_DeriveXi;
      Functions3D[D010] = N_H_BDDF1_3D_DeriveEta;
      Functions3D[D001] = N_H_BDDF1_3D_DeriveZeta;
      Functions3D[D200] = N_H_BDDF1_3D_DeriveXiXi;
      Functions3D[D110] = N_H_BDDF1_3D_DeriveXiEta;
      Functions3D[D101] = N_H_BDDF1_3D_DeriveXiZeta;
      Functions3D[D020] = N_H_BDDF1_3D_DeriveEtaEta;
      Functions3D[D011] = N_H_BDDF1_3D_DeriveEtaZeta;
      Functions3D[D002] = N_H_BDDF1_3D_DeriveZetaZeta;
      PolynomialDegree = 2;
      Accuracy = 1;
      N_BF2Change = 3;
      BF2Change = N_H_BDDF1_3D_Change;
      BaseVectDim = 3;
      break;
    case BF_N_H_BDDF2_3D:
      Dimension = 39;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = N_H_BDDF2_3D_Funct;
      Functions3D[D100] = N_H_BDDF2_3D_DeriveXi;
      Functions3D[D010] = N_H_BDDF2_3D_DeriveEta;
      Functions3D[D001] = N_H_BDDF2_3D_DeriveZeta;
      Functions3D[D200] = N_H_BDDF2_3D_DeriveXiXi;
      Functions3D[D110] = N_H_BDDF2_3D_DeriveXiEta;
      Functions3D[D101] = N_H_BDDF2_3D_DeriveXiZeta;
      Functions3D[D020] = N_H_BDDF2_3D_DeriveEtaEta;
      Functions3D[D011] = N_H_BDDF2_3D_DeriveEtaZeta;
      Functions3D[D002] = N_H_BDDF2_3D_DeriveZetaZeta;
      PolynomialDegree = 3;
      Accuracy = 1;
      N_BF2Change = 6;
      BF2Change = N_H_BDDF2_3D_Change;
      BaseVectDim = 3;
      break;
    case BF_N_H_BDDF3_3D:
      Dimension = 72;
      RefElement = BFRefElements::BFUnitHexahedron;
      Functions3D[D000] = N_H_BDDF3_3D_Funct;
      Functions3D[D100] = N_H_BDDF3_3D_DeriveXi;
      Functions3D[D010] = N_H_BDDF3_3D_DeriveEta;
      Functions3D[D001] = N_H_BDDF3_3D_DeriveZeta;
      Functions3D[D200] = N_H_BDDF3_3D_DeriveXiXi;
      Functions3D[D110] = N_H_BDDF3_3D_DeriveXiEta;
      Functions3D[D101] = N_H_BDDF3_3D_DeriveXiZeta;
      Functions3D[D020] = N_H_BDDF3_3D_DeriveEtaEta;
      Functions3D[D011] = N_H_BDDF3_3D_DeriveEtaZeta;
      Functions3D[D002] = N_H_BDDF3_3D_DeriveZetaZeta;
      PolynomialDegree = 4;
      Accuracy = 1;
      N_BF2Change = 10;
      BF2Change = N_H_BDDF3_3D_Change;
      BaseVectDim = 3;
      break;
    default:
      ErrThrow("unknown BaseFunct2D ", BaseFunct);
      break;
  }
}

void BaseFunctions::GetDerivatives(MultiIndex2D MultiIndex,
                                  const TQuadFormula *formula,
                                  double **values) const
{
  unsigned int N_ = formula->GetN_QuadPoints();
  for(unsigned int i = 0; i < N_; i++)
  {
    auto p = formula->get_point(i);
    GetDerivatives(MultiIndex, p.x, p.y, values[i]);
  }
}

parmoon::Point transform(BFRefElements RefElement, int joint_nr,
                                    double zeta)
{
  switch(RefElement)
  {
    case BFRefElements::BFUnitTriangle:
      switch(joint_nr)
      {
        case 0:
          return {parmoon::Point(0.5*(zeta+1), 0.)};
          break;
        case 1:
          return {parmoon::Point(0.5*(1-zeta), 0.5*(zeta+1))};
          break;
        case 2:
          return {parmoon::Point(0., 0.5*(-zeta+1))};
          break;
      } // endswitch
      break;
    case BFRefElements::BFUnitSquare:
      switch(joint_nr)
      {
        case 0:
          return {parmoon::Point(zeta, -1.)};
          break;
        case 1:
          return {parmoon::Point(1., zeta)};
          break;
        case 2:
          return {parmoon::Point(-zeta, 1.)};
          break;
        case 3:
          return {parmoon::Point(-1., -zeta)};
          break;
      } // endswitch
      break;
    default:
      ErrThrow("unknown reference element ", RefElement);
      break;
  } // endswitch
  ErrThrow("unknown joint_nr ", joint_nr, " on reference element ", RefElement);
}

void BaseFunctions::GetDerivatives(MultiIndex2D MultiIndex,
                                  const TQuadFormula* formula, int joint,
                                  double ** Values) const
{
  unsigned int N_ = formula->GetN_QuadPoints();
  for(unsigned int i = 0; i < N_; i++)
  {
    double zeta = formula->get_point(i).x;
    auto xy = transform(RefElement, joint, zeta);
    GetDerivatives(MultiIndex, xy.x, xy.y, Values[i]);
  }
}


void BaseFunctions::GetDerivatives(MultiIndex2D MultiIndex, int N_Points,
                                  const double *zeta, int joint,
                                  double **Values)
  const
{
  for(int i=0;i<N_Points;i++)
  {
    auto xy = transform(RefElement, joint, zeta[i]);
    GetDerivatives(MultiIndex, xy.x, xy.y, Values[i]);
  }
}

void BaseFunctions::GetDerivatives(MultiIndex3D MultiIndex, 
                        const TQuadFormula *formula, double **values) const
{
  unsigned int n_points = formula->GetN_QuadPoints();
  for(unsigned int i = 0; i < n_points; i++)
  {
    auto p = formula->get_point(i);
    GetDerivatives(MultiIndex, p.x, p.y, p.z, values[i]);
  }
}

parmoon::Point transform(BFRefElements RefElement,
                                             int joint_nr, double t, double s)
{
  switch(RefElement)
  {
    case BFRefElements::BFUnitTetrahedron:
      switch(joint_nr)
      {
        case 0:
          return {parmoon::Point(t, s, 0.)};
          break;
        case 1:
          return {parmoon::Point(s, 0., t)};
          break;
        case 2:
          return {parmoon::Point(t, 1.-t-s, s)};
          break;
        case 3:
          return {parmoon::Point(0., t, s)};
          break;
      } // endswitch
      break;
    case BFRefElements::BFUnitHexahedron:
      switch(joint_nr)
      {
        case 0:
          return {parmoon::Point(t, s, -1.)};
          break;
        case 1:
          return {parmoon::Point(s, -1., t)};
          break;
        case 2:
          return {parmoon::Point(1., s, t)};
          break;
        case 3:
          return {parmoon::Point(-s, 1., t)};
          break;
        case 4:
          return {parmoon::Point(-1, t, s)};
          break;
        case 5:
          return {parmoon::Point(s, t, 1.)};
          break;
      } // endswitch
      break;
    default:
      ErrThrow("unknown reference element ", RefElement);
      break;
  } // endswitch
  ErrThrow("unknown joint_nr ", joint_nr, " on reference element ", RefElement);
}

void BaseFunctions::GetDerivatives(MultiIndex3D MultiIndex,
                                  const TQuadFormula* formula, int joint,
                                  double ** values) const
{
  unsigned int N_ = formula->GetN_QuadPoints();
  for(unsigned int i = 0; i < N_; i++)
  {
    auto p = formula->get_point(i);
    auto xyz = transform(RefElement, joint, p.x, p.y);
    GetDerivatives(MultiIndex, xyz.x, xyz.y,
                   xyz.z, values[i]);
  }
}



void BaseFunctions::ChangeBF(const TCollection *Coll, const TBaseCell *Cell,
                            double *Values) const
{
  ChangeBF(Coll, Cell, 1, &Values);
}

void BaseFunctions::ChangeBF(const TCollection *Coll, const TBaseCell *Cell,
                            int N_Points, double **Values) const
{
  const int *JointArray;
  int i, j, k;
  double *Array;
  const TJoint *joint;
  const TBaseCell *neigh;
  int OwnNum, NeighNum;
  int N_Joints;
  // Hdiv elements are vector valued, in this case the array Values is longer
  int BaseVectDim = this->GetBaseVectDim();

  if(BF2Change == nullptr || Values == nullptr) return;

  // /*
  OwnNum = Coll->get_cell_index(Cell);
  N_Joints = Cell->GetN_Joints();

  for(j=0;j<N_Joints;j++)
  {
    joint = Cell->GetJoint(j);
    neigh = joint->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->get_cell_index(neigh);
      if(NeighNum < OwnNum)
      {
        int maptype = joint->GetMapType(); // always 0 in 2D
        int index = (maptype <= 2) ? 0 : 1;
        JointArray = BF2Change[index][j];
        for(k=0;k<N_Points;k++)
        {
          Array = Values[k];
          if(Array != nullptr)
          {
            for(i=0;i<N_BF2Change;i++)
            {
              for(int m = 0; m < BaseVectDim; ++m)
              {
                Array[JointArray[i] + m*Dimension] *= -1.;
              }
            }
          }
        } // endfor k
      } // endif NeighNum < OwnNum
    } // endif neigh
  } // endfor j
}
