// =======================================================================
// @(#)FESpace.C        1.5 09/15/99
// 
// Class:       TFESpace
// Purpose:     general super class for all finite element spaces
//              special spaces are implemented in subclasses
//
// Author:      Gunar Matthies (04.11.97)
//
// History:     start of implementation 04.11.97 (Gunar Matthies)
//
//              split TFESpace into TFESpacexD (15.04.1998) Volker Behns
//
// =======================================================================

#include "FESpace.h"
#include "BaseCell.h"
#include "MooNMD_Io.h"
#include "HangingNode.h"
#include <algorithm>
#include <numeric>      // std::accumulate

TFESpace::TFESpace(const TCollection *coll, const std::string& name)
:  Name(name), Collection(coll), UsedElements{}, ElementForShape{},
   AllElements{}
{
  // The following values are dummies, they are re-set by
  // the ConstructSpace method of the daughter classes.
  N_DegreesOfFreedom = 0;
  N_Dirichlet = 0;
  N_BoundaryNodes = {};
  N_Inner = 0;
  is_discontinuous_space = false;
}

TFESpace::TFESpace(const TCollection* coll, const std::string& name, int k,
                   int dim)
 : TFESpace(coll, name)
{
  ElementForShape = get_element_for_shape(k, dim);

  if(k < 0 && k > -6 && Collection->includes_hanging_vertices())
  {
    // see e.g.
    // "ON HANGING NODE CONSTRAINTS FOR NONCONFORMING FINITE ELEMENTS USING THE
    // DOUGLAS–SANTOS–SHEEN–YEELEMENT AS AN EXAMPLE", 2017
    // by WOLFGANG BANGERTH, IMBUNM KIM, DONGWOO SHEEN, AND JAERYUN YIM
    ErrThrow("non-conforming elements currently do not work with hanging nodes."
             " It could be implemented but is not trivial.");
  }
  FindUsedElements();
}


TFESpace::TFESpace(const TCollection* coll, const std::string& name,
                   SpaceType type, int k, int dim)
 : TFESpace(coll, name)
{
  // build ElementForShape array
  switch(type)
  {
    // find velo space for discontinuous pressure
    case DiscP_USpace:
    {
      switch(k)
      {
        case 1:
          Output::warn("TFESpace::TFESpace", "This makes no sense.");
          ElementForShape[Triangle] =      C_P1_2D_T_A;
          ElementForShape[Quadrangle] =    C_Q1_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q1_2D_Q_A;
          ElementForShape[Rectangle] =     C_Q1_2D_Q_A;
          ElementForShape[Tetrahedron] =   C_P1_3D_T_A;
          ElementForShape[Hexahedron] =    C_Q1_3D_H_M;
          ElementForShape[Brick] =         C_Q1_3D_H_A;
          break;
        case 2:
          ElementForShape[Triangle] =      C_B2_2D_T_A;
          ElementForShape[Quadrangle] =    C_Q2_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q2_2D_Q_A;
          ElementForShape[Rectangle] =     C_Q2_2D_Q_A;
          ElementForShape[Tetrahedron] =   C_B2_3D_T_A;
          ElementForShape[Hexahedron] =    C_Q2_3D_H_M;
          ElementForShape[Brick] =         C_Q2_3D_H_A;
          break;
        case 3:
          ElementForShape[Triangle] =      C_B3_2D_T_A;
          ElementForShape[Quadrangle] =    C_Q3_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q3_2D_Q_A;
          ElementForShape[Rectangle] =     C_Q3_2D_Q_A;
          // ElementForShape[Tetrahedron] =   C_B3_3D_T_A;
          ElementForShape[Hexahedron] =    C_Q3_3D_H_M;
          ElementForShape[Brick] =         C_Q3_3D_H_A;
        break;
        case 4:
          ElementForShape[Triangle] =      C_B4_2D_T_A;
          ElementForShape[Quadrangle] =    C_Q4_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q4_2D_Q_A;
          ElementForShape[Rectangle] =     C_Q4_2D_Q_A;
          // ElementForShape[Tetrahedron] =   C_B4_3D_T_A;
          ElementForShape[Hexahedron] =    C_Q4_3D_H_M;
          ElementForShape[Brick] =         C_Q4_3D_H_A;
        break;
        case 222: // Scott-Vogelius
          ElementForShape[Triangle] = C_SV2_2D_T_A;
          ElementForShape[Quadrangle] = C_Q2_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q2_2D_Q_A;
          ElementForShape[Rectangle] = C_Q2_2D_Q_A;
          if(dim == 3)
            ErrThrow("Scott-Vogelius not available in 3D");
          break;
        default:
          ErrThrow("Space is not available");
          break;
      }
      break;
    }
    // find pressure space for discontinuous pressure
    case DiscP_PSpace:
    {
      is_discontinuous_space = true;
      switch(k)
      {
        case 0:
          ElementForShape[Triangle] =      C_P0_2D_T_A;
          ElementForShape[Quadrangle] =    C_Q0_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q0_2D_Q_A;
          ElementForShape[Rectangle] =     C_Q0_2D_Q_A;
          ElementForShape[Tetrahedron] =   C_P0_3D_T_A;
          ElementForShape[Hexahedron] =    C_Q0_3D_H_M;
          ElementForShape[Brick] =         C_Q0_3D_H_A;
          break;
        case 1:
          ElementForShape[Triangle] =      D_P1_2D_T_A;
          ElementForShape[Quadrangle] =    D_P1_2D_Q_M;
          ElementForShape[Parallelogram] = D_P1_2D_Q_A;
          ElementForShape[Rectangle] =     D_P1_2D_Q_A;
          ElementForShape[Tetrahedron] =   D_P1_3D_T_A;
          ElementForShape[Hexahedron] =    D_P1_3D_H_M;
          ElementForShape[Brick] =         D_P1_3D_H_A;
          break;
        case 2:
          ElementForShape[Triangle] =      D_P2_2D_T_A;
          ElementForShape[Quadrangle] =    D_P2_2D_Q_M;
          ElementForShape[Parallelogram] = D_P2_2D_Q_A;
          ElementForShape[Rectangle] =     D_P2_2D_Q_A;
          if(dim == 3)
            ErrThrow("This will not work on Tetrahedrons!!");
          // ElementForShape[Tetrahedron] =   D_P2_3D_T_A;
          ElementForShape[Hexahedron] =    D_P2_3D_H_M;
          ElementForShape[Brick] =         D_P2_3D_H_A;
          break;
        case 3:
          ElementForShape[Triangle] =      D_P3_2D_T_A;
          ElementForShape[Quadrangle] =    D_P3_2D_Q_M;
          ElementForShape[Parallelogram] = D_P3_2D_Q_A;
          ElementForShape[Rectangle] =     D_P3_2D_Q_A;
          if(dim == 3)
            ErrThrow("This will not work on Tetrahedrons!!");
          // ElementForShape[Tetrahedron] =   D_P3_3D_T_A;
          ElementForShape[Hexahedron] =    D_P3_3D_H_M;
          ElementForShape[Brick] =         D_P3_3D_H_A;
          break;
        case 4:
          ElementForShape[Triangle] =      D_P4_2D_T_A;
          ElementForShape[Quadrangle] =    D_P4_2D_Q_M;
          ElementForShape[Parallelogram] = D_P4_2D_Q_A;
          ElementForShape[Rectangle] =     D_P4_2D_Q_A;
          if(dim == 3)
            ErrThrow("No 4-th order discontinuous element implemented in 3D");
          break;
        case 5:
          ElementForShape[Triangle] =      D_P4_2D_T_A;
          if(dim == 2)
            Output::print("Using P4 on triangles");
          ElementForShape[Quadrangle] =    D_P5_2D_Q_M;
          ElementForShape[Parallelogram] = D_P5_2D_Q_A;
          ElementForShape[Rectangle] =     D_P5_2D_Q_A;
          if(dim == 3)
            ErrThrow("No 5-th order discontinuous element implemented in 3D");
          break;
        case 6:
          ElementForShape[Triangle] =      D_P4_2D_T_A;
          if(dim == 2)
            Output::print("Using P4 on triangles");
          ElementForShape[Quadrangle] =    D_P6_2D_Q_M;
          ElementForShape[Parallelogram] = D_P6_2D_Q_A;
          ElementForShape[Rectangle] =     D_P6_2D_Q_A;
          if(dim == 3)
            ErrThrow("No 6-th order discontinuous element implemented in 3D");
          break;
        case 7:
          ElementForShape[Triangle] =      D_P4_2D_T_A;
          if(dim == 2)
            Output::print("Using P4 on triangles");
          ElementForShape[Quadrangle] =    D_P7_2D_Q_M;
          ElementForShape[Parallelogram] = D_P7_2D_Q_A;
          ElementForShape[Rectangle] =     D_P7_2D_Q_A;
          if(dim == 3)
            ErrThrow("No 7-th order discontinuous element implemented in 3D");
          break;
        case 92:
          ElementForShape[Triangle] =      D_P2_2D_T_A;
          ElementForShape[Quadrangle] =    D_D2_2D_Q_M;
          ElementForShape[Parallelogram] = D_D2_2D_Q_A;
          ElementForShape[Rectangle] =     D_D2_2D_Q_A;
          if(dim == 3)
            ErrThrow("element not implemented in 3D");
          break;
        case 102: // only function zero
          ElementForShape[Triangle] =      C_P00_2D_T_A;
          ElementForShape[Quadrangle] =    C_Q00_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q00_2D_Q_A;
          ElementForShape[Rectangle] =     C_Q00_2D_Q_A;
          ElementForShape[Tetrahedron] =   C_P00_3D_T_A;
          ElementForShape[Hexahedron] =    C_Q00_3D_H_M;
          ElementForShape[Brick] =         C_Q00_3D_H_A;
          break;
        case 222: // Scott-Vogelius
          ElementForShape[Triangle] =      D_SV1_2D_T_A;
          ElementForShape[Quadrangle] =    D_P1_2D_Q_M;
          ElementForShape[Parallelogram] = D_P1_2D_Q_A;
          ElementForShape[Rectangle] =     D_P1_2D_Q_A;
          if(dim == 3)
            ErrThrow("Scott-Vogelius not available in 3D");
          break;
        default:
          ErrThrow("Space is not available");
          break;
      }
      break;
    }
    // find pressure and velo space for continuous pressure
    case ContP_USpace:
    case ContP_PSpace:
    {
      switch(k)
      {
        case 1:
          ElementForShape[Triangle] =      C_P1_2D_T_A;
          ElementForShape[Quadrangle] =    C_Q1_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q1_2D_Q_A;
          ElementForShape[Rectangle] =     C_Q1_2D_Q_A;
          ElementForShape[Tetrahedron] =   C_P1_3D_T_A;
          ElementForShape[Hexahedron] =    C_Q1_3D_H_M;
          ElementForShape[Brick] =         C_Q1_3D_H_A;
          break;
        case 2:
          ElementForShape[Triangle] =      C_P2_2D_T_A;
          ElementForShape[Quadrangle] =    C_Q2_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q2_2D_Q_A;
          ElementForShape[Rectangle] =     C_Q2_2D_Q_A;
          ElementForShape[Tetrahedron] =   C_P2_3D_T_A;
          ElementForShape[Hexahedron] =    C_Q2_3D_H_M;
          ElementForShape[Brick] =         C_Q2_3D_H_A;
          break;
        case 3:
          ElementForShape[Triangle] =      C_P3_2D_T_A;
          ElementForShape[Quadrangle] =    C_Q3_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q3_2D_Q_A;
          ElementForShape[Rectangle] =     C_Q3_2D_Q_A;
          ElementForShape[Tetrahedron] =   C_P3_3D_T_A;
          ElementForShape[Hexahedron] =    C_Q3_3D_H_M;
          ElementForShape[Brick] =         C_Q3_3D_H_A;
          break;
        case 4:
          ElementForShape[Triangle] =      C_P4_2D_T_A;
          ElementForShape[Quadrangle] =    C_Q4_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q4_2D_Q_A;
          ElementForShape[Rectangle] =     C_Q4_2D_Q_A;
          ElementForShape[Tetrahedron] =   C_P3_3D_T_A;
          if(dim == 3)
            Output::warn("TFESpace3D::TFESpace3D",
                         "No elements of order 4 for Tetrahedron!");
          ElementForShape[Hexahedron] =    C_Q4_3D_H_M;
          ElementForShape[Brick] =         C_Q4_3D_H_A;
          break;
        case 5:
          ElementForShape[Triangle] =      C_P5_2D_T_A;
          ElementForShape[Quadrangle] =    C_Q5_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q5_2D_Q_A;
          ElementForShape[Rectangle] =     C_Q5_2D_Q_A;
          if(dim == 3)
            ErrThrow(k, "-th order Lagrange elements not implemented in 3D");
          break;
        case 6:
          ElementForShape[Triangle] =      C_P6_2D_T_A;
          ElementForShape[Quadrangle] =    C_Q6_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q6_2D_Q_A;
          ElementForShape[Rectangle] =     C_Q6_2D_Q_A;
          if(dim == 3)
            ErrThrow(k, "-th order Lagrange elements not implemented in 3D");
          break;
        case 7:
          ElementForShape[Triangle] =      C_P7_2D_T_A;
          ElementForShape[Quadrangle] =    C_Q7_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q7_2D_Q_A;
          ElementForShape[Rectangle] =     C_Q7_2D_Q_A;
          if(dim == 3)
            ErrThrow(k, "-th order Lagrange elements not implemented in 3D");
          break;
        case 8:
          ElementForShape[Triangle] =      C_P8_2D_T_A;
          ElementForShape[Quadrangle] =    C_Q8_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q8_2D_Q_A;
          ElementForShape[Rectangle] =     C_Q8_2D_Q_A;
          if(dim == 3)
            ErrThrow(k, "-th order Lagrange elements not implemented in 3D");
          break;
        case 9:
          ElementForShape[Triangle] =      C_P9_2D_T_A;
          ElementForShape[Quadrangle] =    C_Q9_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q9_2D_Q_A;
          ElementForShape[Rectangle] =     C_Q9_2D_Q_A;
          if(dim == 3)
            ErrThrow(k, "-th order Lagrange elements not implemented in 3D");
          break;
        case 101: // P1mini
          ElementForShape[Triangle] =      C_P1MINI_2D_T_A;
          ElementForShape[Quadrangle] =    C_Q1_2D_Q_M;
          ElementForShape[Parallelogram] = C_Q1_2D_Q_A;
          ElementForShape[Rectangle] =     C_Q1_2D_Q_A;
          if(dim == 2)
            Output::print("P1MINI works only on triangles");
          if(dim == 3)
            ErrThrow("No P1MINI finite element implemented in 3D");
          break;
        case -2:
          ElementForShape[Triangle] =      C_P9_2D_T_A;
          ElementForShape[Quadrangle] =    B_IB2_2D_Q_M;
          ElementForShape[Parallelogram] = B_IB2_2D_Q_A;
          ElementForShape[Rectangle] =     B_IB2_2D_Q_A;
          if(dim == 2)
            Output::print("Space for triangles not yet implemented !!!!!!");
          ElementForShape[Tetrahedron] =   C_P3_3D_T_A;
          if(dim == 3)
            Output::warn("TFESpace3D::TFESpace3D",
                         "No elements of order 4 for Tetrahedron!");
          ElementForShape[Hexahedron] =    B_IB2_3D_H_M;
          ElementForShape[Brick] =         B_IB2_3D_H_A;
          break;
        //========LOCALPROJECTION=============
        case 100: // Q1+bubble*P0
          ElementForShape[Triangle] =      C_UL1_2D_T_A;
          ElementForShape[Quadrangle] =    C_UL1_2D_Q_M;
          ElementForShape[Parallelogram] = C_UL1_2D_Q_A;
          ElementForShape[Rectangle] =     C_UL1_2D_Q_A;
          if(dim == 3)
            ErrThrow("Q1+bubble*P0 element not implemented in 3D");
          break;
        case 201: // Q2+bubble*P1
          ElementForShape[Triangle] =      C_UL2_2D_T_A;
          ElementForShape[Quadrangle] =    C_UL2_2D_Q_M;
          ElementForShape[Parallelogram] = C_UL2_2D_Q_A;
          ElementForShape[Rectangle] =     C_UL2_2D_Q_A;
          if(dim == 3)
            ErrThrow("Q2+bubble*P1 element not implemented in 3D");
          break;
        case 302: // Q3+bubble*P2
          ElementForShape[Triangle] =      C_UL3_2D_T_A;
          ElementForShape[Quadrangle] =    C_UL3_2D_Q_M;
          ElementForShape[Parallelogram] = C_UL3_2D_Q_A;
          ElementForShape[Rectangle] =     C_UL3_2D_Q_A;
          if(dim == 3)
            ErrThrow("Q3+bubble*P2 element not implemented in 3D");
          break;
        case 403: // Q4+bubble*P3
          ElementForShape[Triangle] =      C_UL4_2D_T_A;
          ElementForShape[Quadrangle] =    C_UL4_2D_Q_M;
          ElementForShape[Parallelogram] = C_UL4_2D_Q_A;
          ElementForShape[Rectangle] =     C_UL4_2D_Q_A;
          if(dim == 3)
            ErrThrow("Q4+bubble*P3 element not implemented in 3D");
          break;
        case 504: // Q5+bubble*P4
          ElementForShape[Triangle] =      C_UL5_2D_T_A;
          ElementForShape[Quadrangle] =    C_UL5_2D_Q_M;
          ElementForShape[Parallelogram] = C_UL5_2D_Q_A;
          ElementForShape[Rectangle] =     C_UL5_2D_Q_A;
          if(dim == 3)
            ErrThrow("Q4+bubble*P3 element not implemented in 3D");
          break;
        case 200: // enriched Q_1-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_UL2_2D_T_A;
          if(dim == 2)
            Output::print("Using usual local projection element on triangles");
          ElementForShape[Quadrangle] =    C_UL2S_2D_Q_M;
          ElementForShape[Parallelogram] = C_UL2S_2D_Q_A;
          ElementForShape[Rectangle] =     C_UL2S_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched Q_1-element not implemented in 3D");
          break;
        case 211: // enriched Q_1-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_UL2_2D_T_A;
          if(dim == 2)
            Output::print("Using usual local projection element on triangles");
          ElementForShape[Quadrangle] =    C_UL2SE_2D_Q_M;
          ElementForShape[Parallelogram] = C_UL2SE_2D_Q_A;
          ElementForShape[Rectangle] =     C_UL2SE_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched Q_1-element not implemented in 3D");
          break;
        case 221: // enriched P_2-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_UL2_2D_T_A;
          if(dim == 2)
            Output::print("Using usual local projection element on triangles");
          ElementForShape[Quadrangle] =    C_M2_2D_Q_M;
          ElementForShape[Parallelogram] = C_M2_2D_Q_A;
          ElementForShape[Rectangle] =     C_M2_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched P_2-element not implemented in 3D");
          break;
        case 301: // enriched Q_2-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_UL3_2D_T_A;
          if(dim == 2)
            Output::print("Using usual local projection element on triangles");
          ElementForShape[Quadrangle] =    C_UL3S_2D_Q_M;
          ElementForShape[Parallelogram] = C_UL3S_2D_Q_A;
          ElementForShape[Rectangle] =     C_UL3S_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched Q_2-element not implemented in 3D");
          break;
        case 312: // enriched Q_2-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_UL3_2D_T_A;
          if(dim == 2)
            Output::print("Using usual local projection element on triangles");
          ElementForShape[Quadrangle] =    C_UL3SE_2D_Q_M;
          ElementForShape[Parallelogram] = C_UL3SE_2D_Q_A;
          ElementForShape[Rectangle] =     C_UL3SE_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched Q_2-element not implemented in 3D");
          break;
        case 322: // enriched P_3-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_UL3_2D_T_A;
          if(dim == 2)
            Output::print("Using usual local projection element on triangles");
          ElementForShape[Quadrangle] =    C_M3_2D_Q_M;
          ElementForShape[Parallelogram] = C_M3_2D_Q_A;
          ElementForShape[Rectangle] =     C_M3_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched P_3-element not implemented in 3D");
          break;
        case 402: // enriched Q_3-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_UL4_2D_T_A;
          if(dim == 2)
            Output::print("Using usual local projection element on triangles");
          ElementForShape[Quadrangle] =    C_UL4S_2D_Q_M;
          ElementForShape[Parallelogram] = C_UL4S_2D_Q_A;
          ElementForShape[Rectangle] =     C_UL4S_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched Q_3-element not implemented in 3D");
          break;
        case 413: // enriched Q_3-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_UL4_2D_T_A;
          if(dim == 2)
            Output::print("Using usual local projection element on triangles");
          ElementForShape[Quadrangle] =    C_UL4SE_2D_Q_M;
          ElementForShape[Parallelogram] = C_UL4SE_2D_Q_A;
          ElementForShape[Rectangle] =     C_UL4SE_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched Q_3-element not implemented in 3D");
          break;
        case 423: // enriched P_4-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_UL4_2D_T_A;
          if(dim == 2)
            Output::print("Using usual local projection element on triangles");
          ElementForShape[Quadrangle] =    C_M4_2D_Q_M;
          ElementForShape[Parallelogram] = C_M4_2D_Q_A;
          ElementForShape[Rectangle] =     C_M4_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched P_4-element not implemented in 3D");
          break;
        case 503: // enriched Q_4-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_UL5_2D_T_A;
          if(dim == 2)
            Output::print("Using usual local projection element on triangles");
          ElementForShape[Quadrangle] =    C_UL5S_2D_Q_M;
          ElementForShape[Parallelogram] = C_UL5S_2D_Q_A;
          ElementForShape[Rectangle] =     C_UL5S_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched Q_4-element not implemented in 3D");
          break;
        case 514: // enriched Q_4-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_UL5_2D_T_A;
          if(dim == 2)
            Output::print("Using usual local projection element on triangles");
          ElementForShape[Quadrangle] =    C_UL5SE_2D_Q_M;
          ElementForShape[Parallelogram] = C_UL5SE_2D_Q_A;
          ElementForShape[Rectangle] =     C_UL5SE_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched Q_4-element not implemented in 3D");
          break;
        case 524: // enriched P_5-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_UL5_2D_T_A;
          if(dim == 2)
            Output::print("Using usual local projection element on triangles");
          ElementForShape[Quadrangle] =    C_M5_2D_Q_M;
          ElementForShape[Parallelogram] = C_M5_2D_Q_A;
          ElementForShape[Rectangle] =     C_M5_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched P_5-element not implemented in 3D");
          break;
        case 604: // enriched Q_5-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_P6_2D_T_A;
          if(dim == 2)
            Output::print("Using usual element on triangles");
          ElementForShape[Quadrangle] =    C_UL6S_2D_Q_M;
          ElementForShape[Parallelogram] = C_UL6S_2D_Q_A;
          ElementForShape[Rectangle] =     C_UL6S_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched Q_5-element not implemented in 3D");
          break;
        case 615: // enriched Q_5-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_P6_2D_T_A;
          if(dim == 2)
            Output::print("Using usual element on triangles");
          ElementForShape[Quadrangle] =    C_UL6SE_2D_Q_M;
          ElementForShape[Parallelogram] = C_UL6SE_2D_Q_A;
          ElementForShape[Rectangle] =     C_UL6SE_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched Q_5-element not implemented in 3D");
          break;
        case 625: // enriched P_6-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_P6_2D_T_A;
          if(dim == 2)
            Output::print("Using usual element on triangles");
          ElementForShape[Quadrangle] =    C_M6_2D_Q_M;
          ElementForShape[Parallelogram] = C_M6_2D_Q_A;
          ElementForShape[Rectangle] =     C_M6_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched P_6-element not implemented in 3D");
          break;
        case 705: // enriched Q_6-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_P7_2D_T_A;
          if(dim == 2)
            Output::print("Using usual element on triangles");
          ElementForShape[Quadrangle] =    C_UL7S_2D_Q_M;
          ElementForShape[Parallelogram] = C_UL7S_2D_Q_A;
          ElementForShape[Rectangle] =     C_UL7S_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched Q_6-element not implemented in 3D");
          break;
        case 716: // enriched Q_6-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_P7_2D_T_A;
          if(dim == 2)
            Output::print("Using usual element on triangles");
          ElementForShape[Quadrangle] =    C_UL7SE_2D_Q_M;
          ElementForShape[Parallelogram] = C_UL7SE_2D_Q_A;
          ElementForShape[Rectangle] =     C_UL7SE_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched Q_6-element not implemented in 3D");
          break;
        case 726: // enriched P_8-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_P8_2D_T_A;
          if(dim == 2)
            Output::print("Using usual element on triangles");
          ElementForShape[Quadrangle] =    C_M8_2D_Q_M;
          ElementForShape[Parallelogram] = C_M8_2D_Q_A;
          ElementForShape[Rectangle] =     C_M8_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched P_8-element not implemented in 3D");
          break;
        case 806: // enriched Q_7-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_P8_2D_T_A;
          if(dim == 2)
            Output::print("Using usual element on triangles");
          ElementForShape[Quadrangle] =    C_UL8S_2D_Q_M;
          ElementForShape[Parallelogram] = C_UL8S_2D_Q_A;
          ElementForShape[Rectangle] =     C_UL8S_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched Q_7-element not implemented in 3D");
          break;
        case 817: // enriched Q_7-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_P8_2D_T_A;
          if(dim == 2)
            Output::print("Using usual element on triangles");
          ElementForShape[Quadrangle] =    C_UL8SE_2D_Q_M;
          ElementForShape[Parallelogram] = C_UL8SE_2D_Q_A;
          ElementForShape[Rectangle] =     C_UL8SE_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched Q_7-element not implemented in 3D");
          break;
        case 827: // enriched P_8-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_P8_2D_T_A;
          if(dim == 2)
            Output::print("Using usual element on triangles");
          ElementForShape[Quadrangle] =    C_M8_2D_Q_M;
          ElementForShape[Parallelogram] = C_M8_2D_Q_A;
          ElementForShape[Rectangle] =     C_M8_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched P_8-element not implemented in 3D");
          break;
        case 907: // enriched Q_8-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_P9_2D_T_A;
          if(dim == 2)
            Output::print("Using usual element on triangles");
          ElementForShape[Quadrangle] =    C_UL9S_2D_Q_M;
          ElementForShape[Parallelogram] = C_UL9S_2D_Q_A;
          ElementForShape[Rectangle] =     C_UL9S_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched Q_8-element not implemented in 3D");
          break;
        case 918: // enriched Q_8-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_P9_2D_T_A;
          if(dim == 2)
            Output::print("Using usual element on triangles");
          ElementForShape[Quadrangle] =    C_UL9SE_2D_Q_M;
          ElementForShape[Parallelogram] = C_UL9SE_2D_Q_A;
          ElementForShape[Rectangle] =     C_UL9SE_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched Q_8-element not implemented in 3D");
          break;
        case 928: // enriched P_9-element, used on Shishkin meshes
          ElementForShape[Triangle] =      C_P9_2D_T_A;
          if(dim == 2)
            Output::print("Using usual element on triangles");
          ElementForShape[Quadrangle] =    C_M9_2D_Q_M;
          ElementForShape[Parallelogram] = C_M9_2D_Q_A;
          ElementForShape[Rectangle] =     C_M9_2D_Q_A;
          if(dim == 3)
            ErrThrow("enriched P_9-element not implemented in 3D");
          break;
        default:
          ErrThrow("Space is not available.");
          break;
      }
      break;
    }
    // find velo space for nonconforming fe
    case Non_USpace:
    {
      is_discontinuous_space = true;
      switch(k)
      {
        case 1:
          ElementForShape[Triangle] =      N_P1_2D_T_A;
          ElementForShape[Quadrangle] =    N_Q1_2D_Q_M;
          ElementForShape[Parallelogram] = N_Q1_2D_Q_A;
          ElementForShape[Rectangle] =     N_Q1_2D_Q_A;
          ElementForShape[Tetrahedron] =   N_P1_3D_T_A;
          ElementForShape[Hexahedron] =    N_Q1_3D_H_M;
          ElementForShape[Brick] =         N_Q1_3D_H_A;
          break;
        case 2:
          ElementForShape[Triangle] =      N_P2_2D_T_A;
          ElementForShape[Quadrangle] =    N_Q2_2D_Q_M;
          ElementForShape[Parallelogram] = N_Q2_2D_Q_A;
          ElementForShape[Rectangle] =     N_Q2_2D_Q_A;
          ElementForShape[Hexahedron] =    N_Q2_3D_H_M;
          ElementForShape[Brick] =         N_Q2_3D_H_A;
          ElementForShape[Tetrahedron] =   N_P2_3D_T_A;
          break;
        case 3:
          ElementForShape[Triangle] =      N_P3_2D_T_A;
          ElementForShape[Quadrangle] =    N_Q3_2D_Q_M;
          ElementForShape[Parallelogram] = N_Q3_2D_Q_A;
          ElementForShape[Rectangle] =     N_Q3_2D_Q_A;
          ElementForShape[Hexahedron] =    N_Q3_3D_H_M;
          ElementForShape[Brick] =         N_Q3_3D_H_A;
          ElementForShape[Tetrahedron] =   N_P3_3D_T_A;
          if(dim == 3)
            Output::warn("TFESpace::TFESpace", "No nonconforming elements of "
                        "order 3 for Tetrahedron!");
          break;
        case 4:
          ElementForShape[Triangle] =      N_P4_2D_T_A;
          ElementForShape[Quadrangle] =    N_Q4_2D_Q_M;
          ElementForShape[Parallelogram] = N_Q4_2D_Q_A;
          ElementForShape[Rectangle] =     N_Q4_2D_Q_A;
          ElementForShape[Hexahedron] =    N_Q4_3D_H_M;
          ElementForShape[Brick] =         N_Q4_3D_H_A;
          ElementForShape[Tetrahedron] =   N_P4_3D_T_A;
          if(dim == 3)
            Output::warn("TFESpace::TFESpace", "No nonconforming elements of "
                         "order 4 for Tetrahedron!");
          break;
        case 5:
          ElementForShape[Triangle] =      N_P5_2D_T_A;
          ElementForShape[Quadrangle] =    N_Q5_2D_Q_M;
          ElementForShape[Parallelogram] = N_Q5_2D_Q_A;
          ElementForShape[Rectangle] =     N_Q5_2D_Q_A;
          if(dim == 3)
            Output::warn("TFESpace::TFESpace", "No nonconforming elements of "
                         "order 5 for Tetrahedron!");
          break;
        case 101:
          ElementForShape[Triangle] =      N_P1MOD_2D_T_A;
          ElementForShape[Quadrangle] =    N_Q1_2D_Q_M;
          ElementForShape[Parallelogram] = N_Q1_2D_Q_A;
          ElementForShape[Rectangle] =     N_Q1_2D_Q_A;
          if(dim == 2)
            Output::print("P1MOD works only on triangles");
          if(dim == 3)
            ErrThrow("No P1MOD finite element implemented in 3D");
          break;
        default:
          ErrThrow("This nonconforming space (order ", k, ") is not available.");
          break;
      }
      break;
    }
    default:
      ErrThrow("Wrong space type ", type);
      break;
  }
  FindUsedElements();
}

TFESpace::TFESpace(const TCollection* coll, const std::string& name,
                   FE_type* fes)
 : TFESpace(coll, name)
{
  AllElements.resize(coll->GetN_Cells());
  std::copy(fes, fes + Collection->GetN_Cells(), AllElements.begin());
  FindUsedElements();
}

TFESpace::~TFESpace()
{
  for(auto hn : HangingNodeArray)
    delete hn;
}

std::array<FE_type, N_SHAPES> TFESpace::get_element_for_shape(int k, int dim)
{
  std::array<FE_type, N_SHAPES> ElementForShape;
  switch(k)
  {
    case 0:
      ElementForShape[Triangle] =      C_P0_2D_T_A;
      ElementForShape[Quadrangle] =    C_Q0_2D_Q_M;
      ElementForShape[Parallelogram] = C_Q0_2D_Q_A;
      ElementForShape[Rectangle] =     C_Q0_2D_Q_A;
      ElementForShape[Tetrahedron] =   C_P0_3D_T_A;
      ElementForShape[Hexahedron] =    C_Q0_3D_H_M;
      ElementForShape[Brick] =         C_Q0_3D_H_A;
      break;
    case 1:
      ElementForShape[S_Line] =        C_P1_1D_L_A;
      ElementForShape[Triangle] =      C_P1_2D_T_A;
      ElementForShape[Quadrangle] =    C_Q1_2D_Q_M;
      ElementForShape[Parallelogram] = C_Q1_2D_Q_A;
      ElementForShape[Rectangle] =     C_Q1_2D_Q_A;
      ElementForShape[Tetrahedron] =   C_P1_3D_T_A;
      ElementForShape[Brick] =         C_Q1_3D_H_A;
      ElementForShape[Hexahedron] =    C_Q1_3D_H_M;
      break;
    case 2:
      ElementForShape[S_Line] =        C_P2_1D_L_A;
      ElementForShape[Triangle] =      C_P2_2D_T_A;
      ElementForShape[Quadrangle] =    C_Q2_2D_Q_M;
      ElementForShape[Parallelogram] = C_Q2_2D_Q_A;
      ElementForShape[Rectangle] =     C_Q2_2D_Q_A;
      ElementForShape[Tetrahedron] =   C_P2_3D_T_A;
      ElementForShape[Brick] =         C_Q2_3D_H_A;
      ElementForShape[Hexahedron] =    C_Q2_3D_H_M;
      break;
    case 3:
      ElementForShape[S_Line] =        C_P3_1D_L_A;
      ElementForShape[Triangle] =      C_P3_2D_T_A;
      ElementForShape[Quadrangle] =    C_Q3_2D_Q_M;
      ElementForShape[Parallelogram] = C_Q3_2D_Q_A;
      ElementForShape[Rectangle] =     C_Q3_2D_Q_A;
      ElementForShape[Tetrahedron] =   C_P3_3D_T_A;
      ElementForShape[Brick] =         C_Q3_3D_H_A;
      ElementForShape[Hexahedron] =    C_Q3_3D_H_M;
      break;
    case 4:
      ElementForShape[Triangle] =      C_P4_2D_T_A;
      ElementForShape[Quadrangle] =    C_Q4_2D_Q_M;
      ElementForShape[Parallelogram] = C_Q4_2D_Q_A;
      ElementForShape[Rectangle] =     C_Q4_2D_Q_A;
      ElementForShape[Tetrahedron] =   C_P3_3D_T_A;
      ElementForShape[Brick] =         C_Q4_3D_H_A;
      ElementForShape[Hexahedron] =    C_Q4_3D_H_M;
      if(dim == 3)
        Output::warn("TFESpace::TFESpace",
                     "No elements of order 4 for Tetrahedron!");
      break;
    case 5:
      ElementForShape[Triangle] =      C_P5_2D_T_A;
      ElementForShape[Quadrangle] =    C_Q5_2D_Q_M;
      ElementForShape[Parallelogram] = C_Q5_2D_Q_A;
      ElementForShape[Rectangle] =     C_Q5_2D_Q_A;
      if(dim == 3)
        ErrThrow(k, "-th order Lagrange elements not implemented in 3D");
      break;
    case 6:
      ElementForShape[Triangle] =      C_P6_2D_T_A;
      ElementForShape[Quadrangle] =    C_Q6_2D_Q_M;
      ElementForShape[Parallelogram] = C_Q6_2D_Q_A;
      ElementForShape[Rectangle] =     C_Q6_2D_Q_A;
      if(dim == 3)
        ErrThrow(k, "-th order Lagrange elements not implemented in 3D");
      break;
    case 7:
      ElementForShape[Triangle] =      C_P7_2D_T_A;
      ElementForShape[Quadrangle] =    C_Q7_2D_Q_M;
      ElementForShape[Parallelogram] = C_Q7_2D_Q_A;
      ElementForShape[Rectangle] =     C_Q7_2D_Q_A;
      if(dim == 3)
        ErrThrow(k, "-th order Lagrange elements not implemented in 3D");
      break;
    case 8:
      ElementForShape[Triangle] =      C_P8_2D_T_A;
      ElementForShape[Quadrangle] =    C_Q8_2D_Q_M;
      ElementForShape[Parallelogram] = C_Q8_2D_Q_A;
      ElementForShape[Rectangle] =     C_Q8_2D_Q_A;
      if(dim == 3)
        ErrThrow(k, "-th order Lagrange elements not implemented in 3D");
      break;
    case 9:
      ElementForShape[Triangle] =      C_P9_2D_T_A;
      ElementForShape[Quadrangle] =    C_Q9_2D_Q_M;
      ElementForShape[Parallelogram] = C_Q9_2D_Q_A;
      ElementForShape[Rectangle] =     C_Q9_2D_Q_A;
      if(dim == 3)
        ErrThrow(k, "-th order Lagrange elements not implemented in 3D");
      break;
    //==========================================================================
    // conforming fe spaces with bubbles on triangles 
    case 21:
      ErrThrow("Conforming space ", k, " with bubble is not implemented");
      break;
    case 22:
      ElementForShape[Triangle] =      C_B2_2D_T_A;
      ElementForShape[Quadrangle] =    C_Q2_2D_Q_M;
      ElementForShape[Parallelogram] = C_Q2_2D_Q_A;
      ElementForShape[Rectangle] =     C_Q2_2D_Q_A;
      ElementForShape[Tetrahedron] =   C_B2_3D_T_A;
      ElementForShape[Brick] =         C_Q2_3D_H_A;
      ElementForShape[Hexahedron] =    C_Q2_3D_H_M;
      break;
    case 23:
      ElementForShape[Triangle] =      C_B3_2D_T_A;
      ElementForShape[Quadrangle] =    C_Q3_2D_Q_M;
      ElementForShape[Parallelogram] = C_Q3_2D_Q_A;
      ElementForShape[Rectangle] =     C_Q3_2D_Q_A;
      if(dim == 3)
        ErrThrow(k-20, "-th order bubble elements not implemented in 3D");
      break;
    case 24:
      ElementForShape[Triangle] =      C_B4_2D_T_A;
      ElementForShape[Quadrangle] =    C_Q4_2D_Q_M;
      ElementForShape[Parallelogram] = C_Q4_2D_Q_A;
      ElementForShape[Rectangle] =     C_Q4_2D_Q_A;
      if(dim == 3)
        ErrThrow(k-20, "-th order bubble elements not implemented in 3D");
      break;
    case 25:
      ElementForShape[Triangle] =      C_P5_2D_T_A;
      ElementForShape[Quadrangle] =    C_Q5_2D_Q_M;
      ElementForShape[Parallelogram] = C_Q5_2D_Q_A;
      ElementForShape[Rectangle] =     C_Q5_2D_Q_A;
      if(dim == 3)
        ErrThrow(k-20, "-th order bubble elements not implemented in 3D");
      Output::print("there are no elements P5 with bubbles for triangles");
      break;
//===============================================================================
    case -1: // P1//PQ1 nonconforming
      ElementForShape[S_Line] =        N_P0_1D_L_A;
      ElementForShape[Triangle] =      N_P1_2D_T_A;
      ElementForShape[Quadrangle] =    N_Q1_2D_Q_M;
      ElementForShape[Parallelogram] = N_Q1_2D_Q_A;
      ElementForShape[Rectangle] =     N_Q1_2D_Q_A;
      ElementForShape[Hexahedron] =    N_Q1_3D_H_M;
      ElementForShape[Brick] =         N_Q1_3D_H_A;
      ElementForShape[Tetrahedron] =   N_P1_3D_T_A;
      break;
    case -2: // P2/Q2 nonconforming
      ElementForShape[Triangle] =      N_P2_2D_T_A;
      ElementForShape[Quadrangle] =    N_Q2_2D_Q_M;
      ElementForShape[Parallelogram] = N_Q2_2D_Q_A;
      ElementForShape[Rectangle] =     N_Q2_2D_Q_A;
      ElementForShape[Hexahedron] =    N_Q2_3D_H_M;
      ElementForShape[Brick] =         N_Q2_3D_H_A;
      ElementForShape[Tetrahedron] =   N_P2_3D_T_A;
      break;
    case -3: // P3/Q3 nonconforming
      ElementForShape[Triangle] =      N_P3_2D_T_A;
      ElementForShape[Quadrangle] =    N_Q3_2D_Q_M;
      ElementForShape[Parallelogram] = N_Q3_2D_Q_A;
      ElementForShape[Rectangle] =     N_Q3_2D_Q_A;
      ElementForShape[Hexahedron] =    N_Q3_3D_H_M;
      ElementForShape[Brick] =         N_Q3_3D_H_A;
      ElementForShape[Tetrahedron] =   N_P3_3D_T_A;
      if(dim == 3)
         Output::warn("TFESpace::TFESpace", "No nonconforming elements of "
                      "order 3 for Tetrahedron!");
      break;
    case -4: // P4/Q4 nonconforming
      ElementForShape[Triangle] =      N_P4_2D_T_A;
      ElementForShape[Quadrangle] =    N_Q4_2D_Q_M;
      ElementForShape[Parallelogram] = N_Q4_2D_Q_A;
      ElementForShape[Rectangle] =     N_Q4_2D_Q_A;
      ElementForShape[Hexahedron] =    N_Q4_3D_H_M;
      ElementForShape[Brick] =         N_Q4_3D_H_A;
      ElementForShape[Tetrahedron] =   N_P4_3D_T_A;
      break;
    case -5: // P5/Q5 nonconforming
      ElementForShape[Triangle] =      N_P5_2D_T_A;
      ElementForShape[Quadrangle] =    N_Q5_2D_Q_M;
      ElementForShape[Parallelogram] = N_Q5_2D_Q_A;
      ElementForShape[Rectangle] =     N_Q5_2D_Q_A;
      if(dim == 3)
        ErrThrow("No 5-th order nonconforming finite element implemented in "
                 "3D");
      break;
    //==========================================================================
    case -101: // P1mod
      ElementForShape[Triangle] =      N_P1MOD_2D_T_A;
      ElementForShape[Quadrangle] =    N_Q1_2D_Q_M;
      ElementForShape[Parallelogram] = N_Q1_2D_Q_A;
      ElementForShape[Rectangle] =     N_Q1_2D_Q_A;
      Output::print("P1MOD works only on triangles");
      if(dim == 3)
        ErrThrow("No P1MOD finite element implemented in 3D");
      break;
    case 101: // P1mini
      ElementForShape[Triangle] =      C_P1MINI_2D_T_A;
      ElementForShape[Quadrangle] =    C_Q1_2D_Q_M;
      ElementForShape[Parallelogram] = C_Q1_2D_Q_A;
      ElementForShape[Rectangle] =     C_Q1_2D_Q_A;
      Output::print("P1MINI works only on triangles");
      if(dim == 3)
        ErrThrow("No P1MINI finite element implemented in 3D");
      break;
    //==========================================================================
    //LOCALPROJECTION
    case 100: // Q1+bubble*P0
      ElementForShape[Triangle] =      C_UL1_2D_T_A;
      ElementForShape[Quadrangle] =    C_UL1_2D_Q_M;
      ElementForShape[Parallelogram] = C_UL1_2D_Q_A;
      ElementForShape[Rectangle] =     C_UL1_2D_Q_A;
      ElementForShape[Hexahedron] =    C_UL1_3D_H_M;
      ElementForShape[Brick] =         C_UL1_3D_H_A;
      // ElementForShape[Tetrahedron] =   C_UL1_3D_T_A;
      break;
    case 201: // Q2+bubble*P1
      ElementForShape[Triangle] =      C_UL2_2D_T_A;
      ElementForShape[Quadrangle] =    C_UL2_2D_Q_M;
      ElementForShape[Parallelogram] = C_UL2_2D_Q_A;
      ElementForShape[Rectangle] =     C_UL2_2D_Q_A;
      ElementForShape[Hexahedron] =    C_UL2_3D_H_M;
      ElementForShape[Brick] =         C_UL2_3D_H_A;
      // ElementForShape[Tetrahedron] =   C_UL2_3D_T_A;
      break;
    case 302: // Q3+bubble*P2
      ElementForShape[Triangle] =      C_UL3_2D_T_A;
      ElementForShape[Quadrangle] =    C_UL3_2D_Q_M;
      ElementForShape[Parallelogram] = C_UL3_2D_Q_A;
      ElementForShape[Rectangle] =     C_UL3_2D_Q_A;
      ElementForShape[Hexahedron] =    C_UL3_3D_H_M;
      ElementForShape[Brick] =         C_UL3_3D_H_A;
      // ElementForShape[Tetrahedron] =   C_UL3_3D_T_A;
      break;
    case 403: // Q4+bubble*P3
      ElementForShape[Triangle] =      C_UL4_2D_T_A;
      ElementForShape[Quadrangle] =    C_UL4_2D_Q_M;
      ElementForShape[Parallelogram] = C_UL4_2D_Q_A;
      ElementForShape[Rectangle] =     C_UL4_2D_Q_A;
      if(dim == 3)
        ErrThrow("No Q4+bubble*P3 finite element implemented in 3D");
      break;
    case 504: // Q5+bubble*P4
      ElementForShape[Triangle] =      C_UL5_2D_T_A;
      ElementForShape[Quadrangle] =    C_UL5_2D_Q_M;
      ElementForShape[Parallelogram] = C_UL5_2D_Q_A;
      ElementForShape[Rectangle] =     C_UL5_2D_Q_A;
      if(dim == 3)
        ErrThrow("No Q5+bubble*P4 finite element implemented in 3D");
      break; 
    //==========================================================================
    // elements used for Shishkin meshes
    case 200: // enriched Q_1-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_UL2_2D_T_A;
      Output::print("Using usual local projection element on triangles");
      ElementForShape[Quadrangle] =    C_UL2S_2D_Q_M;
      ElementForShape[Parallelogram] = C_UL2S_2D_Q_A;
      ElementForShape[Rectangle] =     C_UL2S_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched Q_1-element not implemented in 3D");
      break;
    case 211: // enriched Q_1-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_UL2_2D_T_A;
      Output::print("Using usual local projection element on triangles");
      ElementForShape[Quadrangle] =    C_UL2SE_2D_Q_M;
      ElementForShape[Parallelogram] = C_UL2SE_2D_Q_A;
      ElementForShape[Rectangle] =     C_UL2SE_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched Q_1-element not implemented in 3D");
      break;
    case 221: // enriched P_2-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_UL2_2D_T_A;
      Output::print("Using usual local projection element on triangles");
      ElementForShape[Quadrangle] =    C_M2_2D_Q_M;
      ElementForShape[Parallelogram] = C_M2_2D_Q_A;
      ElementForShape[Rectangle] =     C_M2_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched P_2-element not implemented in 3D");
      break;
    case 301: // enriched Q_2-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_UL3_2D_T_A;
      Output::print("Using usual local projection element on triangles");
      ElementForShape[Quadrangle] =    C_UL3S_2D_Q_M;
      ElementForShape[Parallelogram] = C_UL3S_2D_Q_A;
      ElementForShape[Rectangle] =     C_UL3S_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched Q_2-element not implemented in 3D");
      break;
    case 312: // enriched Q_2-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_UL3_2D_T_A;
      Output::print("Using usual local projection element on triangles");
      ElementForShape[Quadrangle] =    C_UL3SE_2D_Q_M;
      ElementForShape[Parallelogram] = C_UL3SE_2D_Q_A;
      ElementForShape[Rectangle] =     C_UL3SE_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched Q_2-element not implemented in 3D");
      break;
    case 322: // enriched P_3-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_UL3_2D_T_A;
      Output::print("Using usual local projection element on triangles");
      ElementForShape[Quadrangle] =    C_M3_2D_Q_M;
      ElementForShape[Parallelogram] = C_M3_2D_Q_A;
      ElementForShape[Rectangle] =     C_M3_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched P_3-element not implemented in 3D");
      break;
    case 402: // enriched Q_3-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_UL4_2D_T_A;
      Output::print("Using usual local projection element on triangles");
      ElementForShape[Quadrangle] =    C_UL4S_2D_Q_M;
      ElementForShape[Parallelogram] = C_UL4S_2D_Q_A;
      ElementForShape[Rectangle] =     C_UL4S_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched Q_3-element not implemented in 3D");
      break;
    case 413: // enriched Q_3-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_UL4_2D_T_A;
      Output::print("Using usual local projection element on triangles");
      ElementForShape[Quadrangle] =    C_UL4SE_2D_Q_M;
      ElementForShape[Parallelogram] = C_UL4SE_2D_Q_A;
      ElementForShape[Rectangle] =     C_UL4SE_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched Q_3-element not implemented in 3D");
      break;
    case 423: // enriched P_4-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_UL4_2D_T_A;
      Output::print("Using usual local projection element on triangles");
      ElementForShape[Quadrangle] =    C_M4_2D_Q_M;
      ElementForShape[Parallelogram] = C_M4_2D_Q_A;
      ElementForShape[Rectangle] =     C_M4_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched P_4-element not implemented in 3D");
      break;
    case 503: // enriched Q_4-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_UL5_2D_T_A;
      Output::print("Using usual local projection element on triangles");
      ElementForShape[Quadrangle] =    C_UL5S_2D_Q_M;
      ElementForShape[Parallelogram] = C_UL5S_2D_Q_A;
      ElementForShape[Rectangle] =     C_UL5S_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched Q_4-element not implemented in 3D");
      break;
    case 514: // enriched Q_4-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_UL5_2D_T_A;
      Output::print("Using usual local projection element on triangles");
      ElementForShape[Quadrangle] =    C_UL5SE_2D_Q_M;
      ElementForShape[Parallelogram] = C_UL5SE_2D_Q_A;
      ElementForShape[Rectangle] =     C_UL5SE_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched Q_4-element not implemented in 3D");
      break;
    case 524: // enriched P_5-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_UL5_2D_T_A;
      Output::print("Using usual local projection element on triangles");
      ElementForShape[Quadrangle] =    C_M5_2D_Q_M;
      ElementForShape[Parallelogram] = C_M5_2D_Q_A;
      ElementForShape[Rectangle] =     C_M5_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched P_5-element not implemented in 3D");
      break;
    case 604: // enriched Q_5-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_P6_2D_T_A;
      Output::print("Using usual element on triangles");
      ElementForShape[Quadrangle] =    C_UL6S_2D_Q_M;
      ElementForShape[Parallelogram] = C_UL6S_2D_Q_A;
      ElementForShape[Rectangle] =     C_UL6S_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched Q_5-element not implemented in 3D");
      break;
    case 615: // enriched Q_5-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_P6_2D_T_A;
      Output::print("Using usual element on triangles");
      ElementForShape[Quadrangle] =    C_UL6SE_2D_Q_M;
      ElementForShape[Parallelogram] = C_UL6SE_2D_Q_A;
      ElementForShape[Rectangle] =     C_UL6SE_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched Q_5-element not implemented in 3D");
      break;
    case 625: // enriched P_6-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_P6_2D_T_A;
      Output::print("Using usual element on triangles");
      ElementForShape[Quadrangle] =    C_M6_2D_Q_M;
      ElementForShape[Parallelogram] = C_M6_2D_Q_A;
      ElementForShape[Rectangle] =     C_M6_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched P_6-element not implemented in 3D");
      break;
    case 705: // enriched Q_6-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_P7_2D_T_A;
      Output::print("Using usual element on triangles");
      ElementForShape[Quadrangle] =    C_UL7S_2D_Q_M;
      ElementForShape[Parallelogram] = C_UL7S_2D_Q_A;
      ElementForShape[Rectangle] =     C_UL7S_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched Q_6-element not implemented in 3D");
      break;
    case 716: // enriched Q_6-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_P7_2D_T_A;
      Output::print("Using usual element on triangles");
      ElementForShape[Quadrangle] =    C_UL7SE_2D_Q_M;
      ElementForShape[Parallelogram] = C_UL7SE_2D_Q_A;
      ElementForShape[Rectangle] =     C_UL7SE_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched Q_6-element not implemented in 3D");
      break;
    case 726: // enriched P_8-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_P8_2D_T_A;
      Output::print("Using usual element on triangles");
      ElementForShape[Quadrangle] =    C_M8_2D_Q_M;
      ElementForShape[Parallelogram] = C_M8_2D_Q_A;
      ElementForShape[Rectangle] =     C_M8_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched P_8-element not implemented in 3D");
      break;
    case 806: // enriched Q_7-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_P8_2D_T_A;
      Output::print("Using usual element on triangles");
      ElementForShape[Quadrangle] =    C_UL8S_2D_Q_M;
      ElementForShape[Parallelogram] = C_UL8S_2D_Q_A;
      ElementForShape[Rectangle] =     C_UL8S_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched Q_7-element not implemented in 3D");
      break;
    case 817: // enriched Q_7-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_P8_2D_T_A;
      Output::print("Using usual element on triangles");
      ElementForShape[Quadrangle] =    C_UL8SE_2D_Q_M;
      ElementForShape[Parallelogram] = C_UL8SE_2D_Q_A;
      ElementForShape[Rectangle] =     C_UL8SE_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched Q_7-element not implemented in 3D");
      break;
    case 827: // enriched P_8-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_P8_2D_T_A;
      Output::print("Using usual element on triangles");
      ElementForShape[Quadrangle] =    C_M8_2D_Q_M;
      ElementForShape[Parallelogram] = C_M8_2D_Q_A;
      ElementForShape[Rectangle] =     C_M8_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched P_8-element not implemented in 3D");
      break;
    case 907: // enriched Q_8-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_P9_2D_T_A;
      Output::print("Using usual element on triangles");
      ElementForShape[Quadrangle] =    C_UL9S_2D_Q_M;
      ElementForShape[Parallelogram] = C_UL9S_2D_Q_A;
      ElementForShape[Rectangle] =     C_UL9S_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched Q_8-element not implemented in 3D");
      break;
    case 918: // enriched Q_8-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_P9_2D_T_A;
      Output::print("Using usual element on triangles");
      ElementForShape[Quadrangle] =    C_UL9SE_2D_Q_M;
      ElementForShape[Parallelogram] = C_UL9SE_2D_Q_A;
      ElementForShape[Rectangle] =     C_UL9SE_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched Q_8-element not implemented in 3D");
      break;
    case 928: // enriched P_9-element, used on Shishkin meshes
      ElementForShape[Triangle] =      C_P9_2D_T_A;
      Output::print("Using usual element on triangles");
      ElementForShape[Quadrangle] =    C_M9_2D_Q_M;
      ElementForShape[Parallelogram] = C_M9_2D_Q_A;
      ElementForShape[Rectangle] =     C_M9_2D_Q_A;
      if(dim == 3)
        ErrThrow("enriched P_9-element not implemented in 3D");
      break;
    case 102: // only function zero
      ElementForShape[Triangle] =      C_P00_2D_T_A;
      ElementForShape[Quadrangle] =    C_Q00_2D_Q_M;
      ElementForShape[Parallelogram] = C_Q00_2D_Q_A;
      ElementForShape[Rectangle] =     C_Q00_2D_Q_A;
      ElementForShape[Hexahedron] =    C_Q00_3D_H_M;
      ElementForShape[Brick] =         C_Q00_3D_H_A;
      ElementForShape[Tetrahedron] =   C_P00_3D_T_A;
      break;
    //==========================================================================
    // discontinous elements
    case -11:
      ElementForShape[S_Line] =        D_P1_1D_L_A;
      ElementForShape[Triangle] =      D_P1_2D_T_A;
      ElementForShape[Quadrangle] =    D_Q1_2D_Q_M;
      ElementForShape[Parallelogram] = D_Q1_2D_Q_A;
      ElementForShape[Rectangle] =     D_Q1_2D_Q_A;
      ElementForShape[Hexahedron] =    D_Q1_3D_H_M;
      ElementForShape[Brick] =         D_Q1_3D_H_A;
      ElementForShape[Tetrahedron] =   D_P1_3D_T_A;
      break;
    case -12:
      ElementForShape[S_Line] =        D_P2_1D_L_A;
      ElementForShape[Triangle] =      D_P2_2D_T_A;
      ElementForShape[Quadrangle] =    D_Q2_2D_Q_M;
      ElementForShape[Parallelogram] = D_Q2_2D_Q_A;
      ElementForShape[Rectangle] =     D_Q2_2D_Q_A;
      ElementForShape[Hexahedron] =    D_Q2_3D_H_M;
      ElementForShape[Brick] =         D_Q2_3D_H_A;
      ElementForShape[Tetrahedron] =   D_P2_3D_T_A;
      break;
    case -13:
      ElementForShape[Triangle] =      D_P3_2D_T_A;
      ElementForShape[Quadrangle] =    D_Q3_2D_Q_M;
      ElementForShape[Parallelogram] = D_Q3_2D_Q_A;
      ElementForShape[Rectangle] =     D_Q3_2D_Q_A;
      //ElementForShape[Hexahedron] =    D_Q2_3D_H_M;
      //ElementForShape[Brick] =         D_Q2_3D_H_A;
      ElementForShape[Tetrahedron] =   D_P3_3D_T_A;
      if(dim == 3)
         Output::warn("TFESpace::TFESpace", "No discontinuous elements of "
                      "order 3 for Hexahedra!");
      break;
    case -14:
      ElementForShape[Triangle] =      D_P4_2D_T_A;
      ElementForShape[Quadrangle] =    D_Q4_2D_Q_M;
      ElementForShape[Parallelogram] = D_Q4_2D_Q_A;
      ElementForShape[Rectangle] =     D_Q4_2D_Q_A;
      if(dim == 3)
        ErrThrow("No discontinuous elements of order 4 implemented in 3D");
      break;
    //==========================================================================
    case -110: // discontionuous P-elements on quadrangles
      ElementForShape[Triangle] =      D_P1_2D_T_A;
      ElementForShape[Quadrangle] =    D_P1_2D_Q_M;
      ElementForShape[Parallelogram] = D_P1_2D_Q_A;
      ElementForShape[Rectangle] =     D_P1_2D_Q_A;
      ElementForShape[Hexahedron] =    D_P1_3D_H_M;
      ElementForShape[Brick] =         D_P1_3D_H_A;
      ElementForShape[Tetrahedron] =   D_P1_3D_T_A;
      break;
    case -120:
      ElementForShape[Triangle] =      D_P2_2D_T_A;
      ElementForShape[Quadrangle] =    D_P2_2D_Q_M;
      ElementForShape[Parallelogram] = D_P2_2D_Q_A;
      ElementForShape[Rectangle] =     D_P2_2D_Q_A;
      ElementForShape[Hexahedron] =    D_P2_3D_H_M;
      ElementForShape[Brick] =         D_P2_3D_H_A;
      ElementForShape[Tetrahedron] =   D_P2_3D_T_A;
      break;
    case -130:
      ElementForShape[Triangle] =      D_P3_2D_T_A;
      ElementForShape[Quadrangle] =    D_P3_2D_Q_M;
      ElementForShape[Parallelogram] = D_P3_2D_Q_A;
      ElementForShape[Rectangle] =     D_P3_2D_Q_A;
      if(dim == 3)
            ErrThrow("This will not work on Tetrahedrons!!");
      // ElementForShape[Tetrahedron] =   D_P3_3D_T_A;
      ElementForShape[Hexahedron] =    D_P3_3D_H_M;
      ElementForShape[Brick] =         D_P3_3D_H_A;
      break;
    case -140:
      ElementForShape[Triangle] =      D_P4_2D_T_A;
      ElementForShape[Quadrangle] =    D_P4_2D_Q_M;
      ElementForShape[Parallelogram] = D_P4_2D_Q_A;
      ElementForShape[Rectangle] =     D_P4_2D_Q_A;
      if(dim == 3)
        ErrThrow("No discontionuous P-elements of order 4 implemented in 3D");
      break;
    //==========================================================================
    //========Vector basis Raviart-Thomas  element=============
    case 1000:
      ElementForShape[Triangle] =      N_RT0_2D_T_A;
      ElementForShape[Quadrangle] =    N_RT0_2D_Q_M;
      ElementForShape[Parallelogram] = N_RT0_2D_Q_A;
      ElementForShape[Rectangle] =     N_RT0_2D_Q_A;
      // Piola map not yet implemented in THexaTrilinear
      //ElementForShape[Hexahedron] =    N_RT0_3D_H_M;
      ElementForShape[Brick] =         N_RT0_3D_H_A;
      ElementForShape[Tetrahedron] =   N_RT0_3D_T_A;
      break;
    case 1001:
      ElementForShape[Triangle] =      N_RT1_2D_T_A;
      ElementForShape[Quadrangle] =    N_RT1_2D_Q_M;
      ElementForShape[Parallelogram] = N_RT1_2D_Q_A;
      ElementForShape[Rectangle] =     N_RT1_2D_Q_A;
      // Piola map not yet implemented in THexaTrilinear
      //ElementForShape[Hexahedron] =    N_RT1_3D_H_M;
      ElementForShape[Brick] =         N_RT1_3D_H_A;
      ElementForShape[Tetrahedron] =   N_RT1_3D_T_A;
      break;
    case 1002:
      ElementForShape[Triangle] =      N_RT2_2D_T_A;
      ElementForShape[Quadrangle] =    N_RT2_2D_Q_M;
      ElementForShape[Parallelogram] = N_RT2_2D_Q_A;
      ElementForShape[Rectangle] =     N_RT2_2D_Q_A;
      // Piola map not yet implemented in THexaTrilinear
      //ElementForShape[Hexahedron] =    N_RT2_3D_H_M;
      ElementForShape[Brick] =         N_RT2_3D_H_A;
      ElementForShape[Tetrahedron] =   N_RT2_3D_T_A;
      break;
    case 1003:
      ElementForShape[Triangle] =      N_RT3_2D_T_A;
      ElementForShape[Quadrangle] =    N_RT3_2D_Q_M;
      ElementForShape[Parallelogram] = N_RT3_2D_Q_A;
      ElementForShape[Rectangle] =     N_RT3_2D_Q_A;
      ElementForShape[Tetrahedron] =   N_RT3_3D_T_A;
      if(dim == 3)
        Output::warn("TFESpace::TFESpace", "No 3rd order Raviart-Thomas "
                     "elements for Hexahedra implemented!");
      break;
    //==========================================================================
      //========Vector basis BDM  element=============
    case 1011:
      ElementForShape[Triangle] =      N_BDM1_2D_T_A;
      ElementForShape[Quadrangle] =    N_BDM1_2D_Q_M;
      ElementForShape[Parallelogram] = N_BDM1_2D_Q_A;
      ElementForShape[Rectangle] =     N_BDM1_2D_Q_A;
      // Piola map not yet implemented in THexaTrilinear
      //ElementForShape[Hexahedron] =    N_BDM1_3D_H_M;
      ElementForShape[Brick] =         N_BDDF1_3D_H_A;
      ElementForShape[Tetrahedron] =   N_BDDF1_3D_T_A;
      break;
    case 1012:
      ElementForShape[Triangle] =      N_BDM2_2D_T_A;
      ElementForShape[Quadrangle] =    N_BDM2_2D_Q_M;
      ElementForShape[Parallelogram] = N_BDM2_2D_Q_A;
      ElementForShape[Rectangle] =     N_BDM2_2D_Q_A;
      // Piola map not yet implemented in THexaTrilinear
      //ElementForShape[Hexahedron] =    N_BDM2_3D_H_M;
      ElementForShape[Brick] =         N_BDDF2_3D_H_A;
      ElementForShape[Tetrahedron] =   N_BDDF2_3D_T_A;
      break;
    case 1013:
      ElementForShape[Triangle] =      N_BDM3_2D_T_A;
      ElementForShape[Quadrangle] =    N_BDM3_2D_Q_M;
      ElementForShape[Parallelogram] = N_BDM3_2D_Q_A;
      ElementForShape[Rectangle] =     N_BDM3_2D_Q_A;
      // Piola map not yet implemented in THexaTrilinear
      //ElementForShape[Hexahedron] =    N_BDM3_3D_H_M;
      ElementForShape[Brick] =         N_BDDF3_3D_H_A;
      ElementForShape[Tetrahedron] =   N_BDDF3_3D_T_A;
      break;
    //==========================================================================
    //========LOCALPROJECTION WITH EXP BUBBLE=============
    case 1100: // Q1+bubble*P0
      ElementForShape[Triangle] =      C_UL1_2D_T_A;
      Output::print("Using usual LPS bubble element on triangles");
      ElementForShape[Quadrangle] =    C_EL1_2D_Q_M;
      ElementForShape[Parallelogram] = C_EL1_2D_Q_A;
      ElementForShape[Rectangle] =     C_EL1_2D_Q_A;
      if(dim == 3)
        ErrThrow("No Q1+bubble*P0 elements implemented in 3D");
      break;

    default:
      ErrThrow("unknown order ", k);
      break;
  } // endswitch
  return ElementForShape;
}

void TFESpace::set_is_discontinuous_space()
{
  for (auto fe_type : UsedElements)
  {
    switch(fe_type)
    {
      case C_P1_1D_L_A:
      case C_P1_2D_T_A:
      case C_Q1_2D_Q_M:
      case C_Q1_2D_Q_A:
      case C_P1_3D_T_A:
      case C_Q1_3D_H_A:
      case C_Q1_3D_H_M:
      case C_P2_1D_L_A:
      case C_P2_2D_T_A:
      case C_Q2_2D_Q_M:
      case C_Q2_2D_Q_A:
      case C_P2_3D_T_A:
      case C_Q2_3D_H_A:
      case C_Q2_3D_H_M:
      case C_P3_1D_L_A:
      case C_P3_2D_T_A:
      case C_Q3_2D_Q_M:
      case C_Q3_2D_Q_A:
      case C_P3_3D_T_A:
      case C_Q3_3D_H_A:
      case C_Q3_3D_H_M:
      case C_P4_2D_T_A:
      case C_Q4_2D_Q_M:
      case C_Q4_2D_Q_A:
      case C_Q4_3D_H_A:
      case C_Q4_3D_H_M:
      case C_P5_2D_T_A:
      case C_Q5_2D_Q_M:
      case C_Q5_2D_Q_A:
      case C_P6_2D_T_A:
      case C_Q6_2D_Q_M:
      case C_Q6_2D_Q_A:
      case C_P7_2D_T_A:
      case C_Q7_2D_Q_M:
      case C_Q7_2D_Q_A:
      case C_P8_2D_T_A:
      case C_Q8_2D_Q_M:
      case C_Q8_2D_Q_A:
      case C_P9_2D_T_A:
      case C_Q9_2D_Q_M:
      case C_Q9_2D_Q_A:
      case C_B2_2D_T_A:
      case C_B2_3D_T_A:
      case C_B3_2D_T_A:
      case C_B4_2D_T_A:
      case N_P1MOD_2D_T_A:
      case N_Q1_2D_Q_M:
      case N_Q1_2D_Q_A:
      case C_P1MINI_2D_T_A:
      case C_UL1_2D_T_A:
      case C_UL1_2D_Q_M:
      case C_UL1_2D_Q_A:
      case C_UL1_3D_H_M:
      case C_UL1_3D_H_A:
      case C_UL2_2D_T_A:
      case C_UL2_2D_Q_M:
      case C_UL2_2D_Q_A:
      case C_UL2_3D_H_M:
      case C_UL2_3D_H_A:
      case C_UL3_2D_T_A:
      case C_UL3_2D_Q_M:
      case C_UL3_2D_Q_A:
      case C_UL3_3D_H_M:
      case C_UL3_3D_H_A:
      case C_UL4_2D_T_A:
      case C_UL4_2D_Q_M:
      case C_UL4_2D_Q_A:
      case C_UL5_2D_T_A:
      case C_UL5_2D_Q_M:
      case C_UL5_2D_Q_A:
      case C_UL2S_2D_Q_M:
      case C_UL2S_2D_Q_A:
      case C_UL2SE_2D_Q_M:
      case C_UL2SE_2D_Q_A:
      case C_M2_2D_Q_M:
      case C_M2_2D_Q_A:
      case C_UL3S_2D_Q_M:
      case C_UL3S_2D_Q_A:
      case C_UL3SE_2D_Q_M:
      case C_UL3SE_2D_Q_A:
      case C_M3_2D_Q_M:
      case C_M3_2D_Q_A:
      case C_UL4S_2D_Q_M:
      case C_UL4S_2D_Q_A:
      case C_UL4SE_2D_Q_M:
      case C_UL4SE_2D_Q_A:
      case C_M4_2D_Q_M:
      case C_M4_2D_Q_A:
      case C_UL5S_2D_Q_M:
      case C_UL5S_2D_Q_A:
      case C_UL5SE_2D_Q_M:
      case C_UL5SE_2D_Q_A:
      case C_M5_2D_Q_M:
      case C_M5_2D_Q_A:
      case C_UL6S_2D_Q_M:
      case C_UL6S_2D_Q_A:
      case C_UL6SE_2D_Q_M:
      case C_UL6SE_2D_Q_A:
      case C_M6_2D_Q_M:
      case C_M6_2D_Q_A:
      case C_UL7S_2D_Q_M:
      case C_UL7S_2D_Q_A:
      case C_UL7SE_2D_Q_M:
      case C_UL7SE_2D_Q_A:
      case C_M8_2D_Q_M:
      case C_M8_2D_Q_A:
      case C_UL8S_2D_Q_M:
      case C_UL8S_2D_Q_A:
      case C_UL8SE_2D_Q_M:
      case C_UL8SE_2D_Q_A:
      case C_UL9S_2D_Q_M:
      case C_UL9S_2D_Q_A:
      case C_UL9SE_2D_Q_M:
      case C_UL9SE_2D_Q_A:
      case C_M9_2D_Q_M:
      case C_M9_2D_Q_A:
      case C_P00_2D_T_A:
      case C_Q00_2D_Q_A:
      case C_Q00_3D_H_M:
      case C_Q00_3D_H_A:
      case C_P00_3D_T_A:
        is_discontinuous_space = false;
        break;
      case C_P0_2D_T_A:
      case C_Q0_2D_Q_M:
      case C_Q0_2D_Q_A:
      case C_P0_3D_T_A:
      case C_Q0_3D_H_M:
      case C_Q0_3D_H_A:
      case N_P0_1D_L_A:
      case N_P1_2D_T_A:
      case N_Q1_3D_H_M:
      case N_Q1_3D_H_A:
      case N_P1_3D_T_A:
      case N_P2_2D_T_A:
      case N_Q2_2D_Q_M:
      case N_Q2_2D_Q_A:
      case N_Q2_3D_H_M:
      case N_Q2_3D_H_A:
      case N_P2_3D_T_A:
      case N_P3_2D_T_A:
      case N_Q3_2D_Q_M:
      case N_Q3_2D_Q_A:
      case N_Q3_3D_H_M:
      case N_Q3_3D_H_A:
      case N_P3_3D_T_A:
      case N_P4_2D_T_A:
      case N_Q4_2D_Q_M:
      case N_Q4_2D_Q_A:
      case N_Q4_3D_H_M:
      case N_Q4_3D_H_A:
      case N_P4_3D_T_A:
      case N_P5_2D_T_A:
      case N_Q5_2D_Q_M:
      case N_Q5_2D_Q_A:
      case D_P1_1D_L_A:
      case D_P1_2D_T_A:
      case D_Q1_2D_Q_M:
      case D_Q1_2D_Q_A:
      case D_Q1_3D_H_M:
      case D_Q1_3D_H_A:
      case D_P1_3D_T_A:
      case D_P2_1D_L_A:
      case D_P2_2D_T_A:
      case D_Q2_2D_Q_M:
      case D_Q2_2D_Q_A:
      case D_Q2_3D_H_M:
      case D_Q2_3D_H_A:
      case D_P2_3D_T_A:
      case D_P3_2D_T_A:
      case D_Q3_2D_Q_M:
      case D_Q3_2D_Q_A:
      case D_P3_3D_T_A:
      case D_P4_2D_T_A:
      case D_Q4_2D_Q_M:
      case D_Q4_2D_Q_A:
      case D_P1_2D_Q_M:
      case D_P1_2D_Q_A:
      case D_P1_3D_H_M:
      case D_P1_3D_H_A:
      case D_P2_2D_Q_M:
      case D_P2_2D_Q_A:
      case D_P2_3D_H_M:
      case D_P2_3D_H_A:
      case D_P3_2D_Q_M:
      case D_P3_2D_Q_A:
      case D_P3_3D_H_M:
      case D_P3_3D_H_A:
      case D_P4_2D_Q_M:
      case D_P4_2D_Q_A:
      case N_RT0_2D_T_A:
      case N_RT0_2D_Q_M:
      case N_RT0_2D_Q_A:
      case N_RT0_3D_H_A:
      case N_RT0_3D_T_A:
      case N_RT1_2D_T_A:
      case N_RT1_2D_Q_M:
      case N_RT1_2D_Q_A:
      case N_RT1_3D_H_A:
      case N_RT1_3D_T_A:
      case N_RT2_2D_T_A:
      case N_RT2_2D_Q_M:
      case N_RT2_2D_Q_A:
      case N_RT2_3D_H_A:
      case N_RT2_3D_T_A:
      case N_RT3_2D_T_A:
      case N_RT3_2D_Q_M:
      case N_RT3_2D_Q_A:
      case N_RT3_3D_T_A:
      case N_BDM1_2D_T_A:
      case N_BDM1_2D_Q_M:
      case N_BDM1_2D_Q_A:
      case N_BDDF1_3D_H_A:
      case N_BDDF1_3D_T_A:
      case N_BDM2_2D_T_A:
      case N_BDM2_2D_Q_M:
      case N_BDM2_2D_Q_A:
      case N_BDDF2_3D_H_A:
      case N_BDDF2_3D_T_A:
      case N_BDM3_2D_T_A:
      case N_BDM3_2D_Q_M:
      case N_BDM3_2D_Q_A:
      case N_BDDF3_3D_H_A:
      case N_BDDF3_3D_T_A:
      case C_EL1_2D_Q_M:
      case C_EL1_2D_Q_A:
        is_discontinuous_space = true;
        break;
      default:
        ErrThrow("unknown FE_type ", fe_type);
        break;
    } // endswitch
    if (is_discontinuous_space == true)
    {
      break;
    }
  } // endfor fe_types
}



unsigned int TFESpace::get_n_local_dof(unsigned int i) const
{
  if((long int)i >= Collection->GetN_Cells())
    ErrThrow("TFESpace::get_n_local_dof, given cell number ", i,
             " is too large, maximum is ", Collection->GetN_Cells());
  return BeginIndex[i+1] - BeginIndex[i];
}

unsigned int TFESpace::get_max_n_local_dof() const
{
  unsigned int n_cells = Collection->GetN_Cells();
  unsigned int max = 0;
  for(unsigned int i = 0u; i < n_cells; ++i)
  {
    max = std::max(max, (unsigned)(BeginIndex[i+1] - BeginIndex[i]));
  }
  return max;
}

bool TFESpace::is_dof_in_cell(int dof, int cell_index) const
{
  auto dofs = GetGlobalDOF(cell_index);
  unsigned int n = get_n_local_dof(cell_index);
  return std::any_of(dofs, dofs + n, [dof](int d){ return dof == d; });
}

const THangingNode * TFESpace::get_hanging_node(unsigned int hanging_dof_index)
  const
{
  if(hanging_dof_index >= (unsigned int) get_n_hanging())
  {
    ErrThrow("there are only ", get_n_hanging(), " hanging dofs in this fe "
             "space. You wanted to access one with index ", hanging_dof_index);
  }
  return HangingNodeArray[hanging_dof_index];
}

int TFESpace::get_n_non_dirichlet_boundary() const
{
  return std::accumulate(std::begin(N_BoundaryNodes), std::end(N_BoundaryNodes),
                         0);
}

FE_type TFESpace::get_fe_type(int i) const
{
  if(!AllElements.empty())
    return AllElements[i];
  else
    return ElementForShape[Collection->GetCell(i)->GetType()];
}

const FiniteElement& TFESpace::get_fe(unsigned int cell_number) const
{
  // find corresponding cell
  unsigned int n_cells = this->Collection->GetN_Cells();
  if(cell_number >= n_cells)
    ErrThrow("unable to find the finite element for cell ", cell_number, 
             ". There are only ", n_cells, " cells");
  // find finite element id
  FE_type fe_type = this->get_fe_type(cell_number);
  // get the finite element
  auto it = elements.find(fe_type);
  if (it != elements.end())
  {
    return it->second;
  }
  else
  {
    ErrThrow("no finite element of type ", fe_type, " stored in this TFESpace");
  }
}

int TFESpace::GetBaseVectDim() const 
{
  // the desired information is stored in the BasisFunction2D object. We take 
  // the one on the first cell, on all other cells it should be the same
  return this->get_fe(0).GetBaseFunct()->GetBaseVectDim();
}

int TFESpace::getFEDegree(const TBaseCell *cell) const
{
  auto FE = this->get_fe(cell->GetCellIndex());
  return FE.GetBaseFunct()->GetPolynomialDegree();
}

/** write info on fespace into file */
int TFESpace::Write(const std::string& filename)
{
  int header[4];
  int N_LocalDOF;

  std::ofstream dat(filename);
  if(!dat)
  {
    cerr << "cannot open file '" << filename << "' for output" << endl;
    return -1;
  }

  int n_cells= Collection->GetN_Cells();
  N_LocalDOF = BeginIndex[n_cells];
  header[0] = n_cells;
  header[1] = N_DegreesOfFreedom;
  header[3] = N_LocalDOF;

  dat.write((char *)header, sizeof(int)*4);
  
  dat.write((char *)&BeginIndex[0], sizeof(int)*(n_cells +1));
  dat.write((char *)&GlobalNumbers[0], sizeof(int)*N_LocalDOF);
  
  dat.close();

  return 0;
}

void TFESpace::fill_begin_index()
{
  int N_Cells = Collection->GetN_Cells();
  BeginIndex.clear();
  BeginIndex.resize(N_Cells + 1);
  BeginIndex[0] = 0;
  for(int i = 0; i < N_Cells; i++)
  {
    BeginIndex[i+1] = BeginIndex[i] + get_fe(i).GetN_DOF();
  }
}

void TFESpace::FindUsedElements()
{
  std::array<int, n_FiniteElements> used_fe = {}; // initialize with zeros
  for(int i = 0, n_cells = Collection->GetN_Cells(); i < n_cells; i++)
  {
    used_fe[get_fe_type(i)] = 1;
  }

  // count the number of different finite elements used
  int N_UsedElements = std::accumulate(used_fe.begin(), used_fe.end(), 0);
  
  UsedElements.resize(N_UsedElements);
  for(int i = 0, j = 0; i < n_FiniteElements; i++)
  {
    if(used_fe[i])
    {
      UsedElements[j] = (FE_type)i;
      j++;
    }
  }
  
  // Construct and store all used finite elements.
  for(auto u : UsedElements)
  {
    elements.insert(std::make_pair(u, FiniteElement (u)));
  }

  // Output::print("N_UsedElements: ", N_UsedElements);
  // for(int i = 0; i < N_UsedElements; i++)
  //   Output::print("UsedElement[", i, "]: ", UsedElements[i]);
#ifdef _MPI
  // making sure all elements used here are usable in mpi:
  for(auto& p : elements)
  {
    auto FEDesc0_Obj = p.second.GetFEDesc();
    if(!(FEDesc0_Obj->IsEdgeVertData_Filled()) )
    {
      ErrThrow("FESpace Error! Edge and vertex data are not set in "
               "FEDescriptor for this FE");
    }
  }
#endif
  int N_Cells = Collection->GetN_Cells();
  this->Collection->mark_all_cells();
  this->fill_begin_index();

  GlobalNumbers.clear();
  GlobalNumbers.resize(BeginIndex[N_Cells], -1); // initialize with -1

  // Now the underlying finite elements are known and the variable
  // TFESpace::is_discontinuous_space can be set
  set_is_discontinuous_space();
}

void TFESpace::info() const
{
  using namespace Output;
  print<2>("finite element space '", Name, "' with ", this->get_n_dof(),
           " degerees of freedom on ", this->GetN_Cells(), " cells");
  print<2>(" n_active ", get_n_active(), "  n_hanging ", get_n_hanging());
  for(int c = 0; c < this->GetN_Cells(); ++c)
  {
    const int * localDofs = this->GetGlobalDOF(c);
    int nLocalDof = get_n_local_dof(c);
    std::stringstream dofsInCell;
    for(int ld = 0; ld < nLocalDof; ++ld)
    {
      dofsInCell << "  " << localDofs[ld];
    }
    auto cell = GetCollection()->GetCell(c);
    double x = 0., y = 0., z = 0.;
    int n_v = cell->GetN_Vertices();
    for(int iv = 0; iv < n_v; ++iv)
    {
      auto v = cell->GetVertex(iv);
      x += v->GetX();
      y += v->GetY();
      z += v->GetZ();
    }
    print<3>("cell ", c, " has dofs ", dofsInCell.str(), "   cell center (",
             x/n_v, ",", y/n_v, ",", z/n_v, ")");
  }
  int n_hanging = this->get_n_hanging();
  for(int i = 0; i < n_hanging; ++i)
  {
    auto hn = this->get_hanging_node(i);
    int n_coupled = hn->GetN_Nodes();
    auto coupled_dofs = hn->GetDOF();
    auto coupling_factors = hn->GetCoeff();
    std::stringstream s;
    for(int j = 0; j < n_coupled; ++j)
    {
      s << "  " << coupled_dofs[j];
    }
    s << "  and corresponding coefficients ";
    for(int j = 0; j < n_coupled; ++j)
    {
      s << "  " << coupling_factors[j];
    }
    auto partners = hn->get_partners();
    if(!partners.empty())
    {
      s << "  partner dofs:";
      for(int p : partners)
      {
        s << " " << p;
      }
    }
    print<2>("hanging dof ", i, " has coupled dofs", s.str());
  }
}


void sort_hanging_nodes(
  std::vector<std::pair<const THangingNode*, int>>& to_be_sorted,
  int ActiveBound, int HangingBound);

void TFESpace::compute_various_numbers(
  const std::array<int, N_DiffBoundNodeTypes>& BoundaryUpperBound,
  const std::vector<THangingNode *>& VHN, const std::vector<int>& HNNumbers,
  bool flag_2d)
{
  // find global numbers
  int total_boundary_upper_bound = std::accumulate(BoundaryUpperBound.begin(),
                                                   BoundaryUpperBound.end(), 0);
  int FIRSTMARK = -10; // todo: this magic number needs to go away.

  int N_Cells = this->Collection->GetN_Cells();
  int SumLocDOF = BeginIndex[N_Cells];
  int n = VHN.size();
  for(int i = 0, m = -SumLocDOF - total_boundary_upper_bound + FIRSTMARK; i < n; i++)
  {
    int j=HNNumbers.at(i);
    int k;
    while( (k=GlobalNumbers[j]) > -1 )
      j=k;

    GlobalNumbers[j] = m;
    m--;
  }

  std::vector<int> BoundMark(N_DiffBoundNodeTypes);
  int l = FIRSTMARK;
  for(int i=0;i<N_DiffBoundNodeTypes;i++)
  {
    l -= BoundaryUpperBound[i];
    BoundMark[i] = l;
  }
  int InnerMark = l;
  int SlaveMark = l - SumLocDOF;

  //for(i=0;i<N_DiffBoundNodeTypes;i++)
  //    Output::print(i, "   ", BoundMark[i]);
  // Output::print("InnerMark: ", InnerMark);
  // Output::print("SlaveMark: ", SlaveMark);

  int DirichletCounter=0;
  l = 0;
  std::vector<int> BoundCounter(N_DiffBoundNodeTypes);
  for(int i=0;i<N_DiffBoundNodeTypes;i++)
  {
    l += BoundaryUpperBound[i];
    BoundCounter[i] = l;
  }
  int count = l;

  N_Dirichlet = 0;
  N_BoundaryNodes.fill(0);
  N_Inner = 0;
  int N_Slave = 0;

  for(int i=0;i<N_Cells;i++)
  {
    // Output::print("cell: ", i);
    int * J_K0 = &GlobalNumbers[BeginIndex[i]];
    int k = this->get_n_local_dof(i);
    for(int j=0;j<k;j++)
    {
      int l = J_K0[j];
      if (l < -1)
      {
        // Output::print("new node");
        if(l<=SlaveMark)
        {
          // hanging node
          J_K0[j] = -l + (flag_2d ?  0 : SumLocDOF);
          N_Slave++;
          // Output::print("slave");
          continue;
        }

        if(l<=InnerMark)
        {
          // inner node
          J_K0[j] = count;
          count++;
          N_Inner++;
          // Output::print("inner");
          continue;
        }

        int m=N_DiffBoundNodeTypes-1;
        for(; m >= 0; m--)
        {
          if(l<=BoundMark[m])
          {
            J_K0[j] = BoundCounter[m];
            BoundCounter[m]++;
            N_BoundaryNodes[m]++;
            // Output::print("type " << m);
            m = -2;
            break;
          }
        }

        if(m!=-2 && l<=FIRSTMARK) // no match in loop above
        {
          // Dirichlet nodes
          J_K0[j] = DirichletCounter;
          DirichletCounter++;
          N_Dirichlet++;
          // Output::print("Dirichlet");
          continue;
        }

      } // l < -1
      else
      {
        if (l >= 0)
        {
          J_K0[j] = GlobalNumbers[l];
        }
        else
        {
          Output::print("J_K0[j]==-1 at locdof: ", j, " in ele: ", i);
        }
      } // l >= -1
    } // endfor j
  } // endfor i
  
  /*
  Output::print("N_Inner: ", N_Inner);
  Output::print("N_Slave: ", N_Slave);
  Output::print("N_Dirichlet: ", N_Dirichlet);
  for(int i=0;i<N_DiffBoundNodeTypes;i++)
    Output::print(i, " N_BoundaryNodes: ", N_BoundaryNodes[i]);
  */

  // create real numbers
  l = std::accumulate(N_BoundaryNodes.begin(), N_BoundaryNodes.end(), 0);
  
  int DirichletOffset = N_Inner + N_Slave + l - 0;
  // Output::print("DirichletOffset: ", DirichletOffset);
  int SlaveOffset = N_Inner + l - (flag_2d ? 1 : 2)*SumLocDOF - total_boundary_upper_bound + FIRSTMARK;
  int InnerOffset = 0 - total_boundary_upper_bound;

  l = N_Inner;
  std::vector<int> BoundOffset(N_DiffBoundNodeTypes);
  for(int i = 0, m = 0;i<N_DiffBoundNodeTypes;i++)
  {
    m += BoundaryUpperBound[i];
    BoundOffset[i] = l-m;
    l += N_BoundaryNodes[i];
  }

  int DirichletMark = BoundaryUpperBound[0];
  InnerMark = (flag_2d ? 1 : 2)*SumLocDOF-FIRSTMARK;
  l = BoundaryUpperBound[0];
  for(int i=0;i<N_DiffBoundNodeTypes-1;i++)
  {
    l += BoundaryUpperBound[i+1];
    BoundMark[i] = l;
  }
  BoundMark[N_DiffBoundNodeTypes-1] = l;

  // Output::print("DirichletMark: ", DirichletMark);
  // Output::print("InnerMark: ", InnerMark);

  for(int i=0;i<SumLocDOF;i++)
  {
    int n=GlobalNumbers[i];
    // Output::print(i, "  ", n);
    if(n<DirichletMark)
    {
      // Dirichlet node
      GlobalNumbers[i] += DirichletOffset;
      // Output::print("Diri", GlobalNumbers[i]);
    }
    else
    {
      if(n<BoundMark[N_DiffBoundNodeTypes-1])
      {
        // non Dirichlet boundary type
        for(int m=0;m<N_DiffBoundNodeTypes;m++)
        {
          if(n<BoundMark[m])
          {
            // node of type m
            GlobalNumbers[i] += BoundOffset[m];
            // Output::print("type: " << m);
            break;
          }
        } // endfor m
      }
      else
      {
        // no boundary node
        if(n<InnerMark)
        {
          // inner node
          GlobalNumbers[i] += InnerOffset;
          // Output::print("inner");
        }
        else
        {
          // slave node
          GlobalNumbers[i] += SlaveOffset;
          // Output::print("slave: ");
        }
      }
    } // non Dirichlet node
  } // endfor i

  // fill in information for hanging nodes
  HangingNodeArray.resize(N_Slave);

  for(int i=0;i<N_Slave;i++)
  {
    auto hn=VHN.at(i);
    int k=hn->GetN_Nodes();
    int * v=hn->GetDOF();
    for(int j=0;j<k;j++)
      v[j] = GlobalNumbers[v[j]];
    for(int& p : hn->get_partners())
      p = GlobalNumbers[p];

    HangingNodeArray[i] = hn;

  }

  l = N_Inner;
  for(int i=0;i<N_DiffBoundNodeTypes;i++)
  {
    l += N_BoundaryNodes[i];
  }
  int ActiveBound = l;
  N_DegreesOfFreedom = ActiveBound + N_Slave + N_Dirichlet;

  sorted_hanging_nodes.resize(N_Slave);
  for(int i = 0; i < N_Slave; ++i)
    sorted_hanging_nodes[i] = {HangingNodeArray[i], ActiveBound + i};
  
  sort_hanging_nodes(sorted_hanging_nodes, ActiveBound, ActiveBound + N_Slave);
}

void sort_hanging_nodes(
  std::vector<std::pair<const THangingNode*, int>>& to_be_sorted,
  int begin_hanging, int end_hanging)
{
  // the set of hanging nodes form an acyclic directed graph where an edge from
  // hn1 to hn2 means "hn2 is a coupling dof of hn1". This can not be sorted
  // using the std::sort algorithm, instead we need to do "topological sorting".
  // Wikipedia (https://en.wikipedia.org/wiki/Topological_sorting) describes
  // Kahn's algorithm which needs to know the incoming edges of a node, which we
  // do not store here, we only store the outgoing edges, i.e. a hanging node
  // stores the dofs with which it couples but not if another hanging node
  // couples with it.
  // Instead we use another algorithm which it is based on depth-first search.
  // This algorithm starts with an empty vector which is then gradually filled.
  // Some notes regarding the implementation:
  // - usually the entries are written to the beginning of the vector, but since
  //   that is not recommended for std::vectors, we append them at the end and
  //   later reverse the vector
  // - the algorithm uses a sort function which calls itself recursively
  // - this does not change the input `to_be_sorted` in place. Maybe this is
  //   possible, but most of the times there are not too many hanging nodes
  bool check_result = false; // check if sorting was successful
  using HNP = std::pair<const THangingNode*, int>;
  
  size_t n_hn = to_be_sorted.size();
  std::vector<HNP> sorted_hn;
  std::vector<bool> visited(n_hn, false);
  // recursive lambda, see https://stackoverflow.com/a/4081391
  std::function<void(HNP)> sort;
  sort = [&sorted_hn, &visited, &to_be_sorted, begin_hanging, end_hanging, &sort]
         (HNP hn)
  {
    visited.at(hn.second-begin_hanging) = true;
    int n_coupling_nodes = hn.first->GetN_Nodes();
    auto coupling_dofs = hn.first->GetDOF();
    for(int i = 0; i < n_coupling_nodes; ++i)
    {
      int c_dof = coupling_dofs[i];
      bool is_hanging = c_dof >= begin_hanging && c_dof < end_hanging;
      if(is_hanging && !visited[c_dof-begin_hanging])
        sort(to_be_sorted[c_dof-begin_hanging]);
    }
    sorted_hn.push_back(hn);
  };
  
  // the actual sorting is done here
  for(size_t i = 0; i < n_hn; ++i)
  {
    if(visited[i])
      continue;
    sort(to_be_sorted[i]);
  }
  std::reverse(sorted_hn.begin(), sorted_hn.end());
  
  if(n_hn != sorted_hn.size())
  {
    ErrThrow("not all hanging dof in sorted vector ", n_hn, " != ",
             sorted_hn.size());
  }
  to_be_sorted = sorted_hn;
  
  if(check_result)
  {
    // the following lambda detects if there is a direct edge from hn1 to hn2,
    // i.e., if hn1<hn2. This is the case whenever hn2 is a coupling dof of hn1.
    auto compare = [](HNP hn1, HNP hn2)
    {
      int n_coupling_nodes1 = hn1.first->GetN_Nodes();
      auto coupling_dofs1 = hn1.first->GetDOF();
      for(int i = 0; i < n_coupling_nodes1; ++i)
      {
        int coupling_dof1 = coupling_dofs1[i];
        if(coupling_dof1 == hn2.second)
          return true;
      }
      return false;
    };

    // check if everything is correctly ordered:
    for(size_t i = 0; i < n_hn; ++i)
    {
      for(size_t j = i+1; j < n_hn; ++j)
      {
        if(compare(to_be_sorted[j], to_be_sorted[i]))
          ErrThrow("wrong sorting ", to_be_sorted[j].second, " ",
                   to_be_sorted[i].second);
      }
    }
  }
}
