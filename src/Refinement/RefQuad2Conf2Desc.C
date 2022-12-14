// =======================================================================
// @(#)RefQuad2Conf2Desc.C        1.2 09/13/99
//
// Class:       TRefQuad2Conf2Desc
// Purpose:     refinement descriptor for conforming closure of a
//              quadrangle with 2 hanging node
//
// Author:      Matthias Ebeling  20.08.99
//
// =======================================================================

#include <RefQuad2Conf2Desc.h>

// Constructor
TRefQuad2Conf2Desc::TRefQuad2Conf2Desc(const TShapeDesc *shape) : TRefDesc(shape)
{
static const Shapes DatChildType[] = { Quadrangle, Quadrangle,
                                       Quadrangle };

static const Refinements DatEdgeType[] = { NoRef, LineReg, LineReg, NoRef};

static const int DatChildVertex[][QUAD2ConfMAXN_VpC] =
                 { {0, 4, 6, 5},  {4, 1, 2, 6},  {5, 6, 2, 3}};
static const int DatVertexChild[][QUAD2ConfMAXN_CpV] =
                 { {0},  {1},  {1,2},  {2},  {0, 1},  {0, 2},  {0, 1, 2}};
                 
static const int DatVertexChildIndex[][QUAD2ConfMAXN_CpV] =
               { {0},  {1},  {2, 2},  {3},  {1, 0},  {3, 0},  {2, 3, 1}};

static const int DatVertexChildLen[] =
               { 1,  1,  2,  1,  2,  2,  3};

static const int DatChildEdge[][QUAD2ConfMAXN_EpC] =
               { {0, 6, 8, 5},  {1, 2, 7, 6},  {8, 7, 3, 4}};
static const int DatEdgeChild[][QUAD2ConfMAXN_CpE] =
               { {0},  {1},  {1},  {2},  {2},  {0},  {0, 1},  {1, 2},
                 {0, 2}};
static const int DatEdgeChildIndex[][QUAD2ConfMAXN_CpE] = 
               { {0},  {0},  {1},  {2},  {3},  {3},  {1, 3},  {2, 1},
                 {2, 0}};
static const int DatEdgeChildLen[] =
               { 1,  1,  1,  1,  1,  1,  2,  2,  2};

static const int DatEdgeVertex[][2] =
               { {0, 4},  {4, 1},  {1, 2},  {2, 3},  {3, 5},  {5, 0},
                 {4, 6},  {2, 6},  {5, 6}};
static const int DatVertexEdge[][QUAD2ConfMAXN_EpV] = 
               { {0, 5},  {1, 2},  {2, 3},  {3, 4},  {0, 1, 6},
                 {5, 8, 4},  {6, 7, 8}};
static const int DatVertexEdgeIndex[][QUAD2ConfMAXN_EpV] =
               { {0, 1},  {1, 0},  {1, 0},  {1,0},  {1, 0, 0},
                 {0, 0, 1},  {1, 1, 1}};
static const int DatVertexEdgeLen[] =
               { 2,  2,  2,  2,  3,  3,  3};

static const int DatNewVertexEqOldVertex[] = { 0, 1, 2, 3};
static const int DatNewVertexEqOldVertexIndex[] = { 2, 3, 0, 1};

static const int DatInteriorVertexOfCell[] = { 6};
static const double DatPositionOfIntVert[][QUADN_V] =
                  { {0.25, 0.25, 0.25, 0.25}};
                                 
static const int DatInteriorEdgeOfCell[] = { 6, 7, 8};
static const int DatInteriorVertexOfEdge[][QUAD2ConfMAXN_iVpE] =
               { {-1}, {5}, {4}, {-1}};
static const int DatInteriorVertexOfEdgeLen[] = { 0, 1, 1, 0};

static const int DatOldEdgeNewVertex[][QUAD2ConfMAXN_nVpoE] =
               { {2, 3},  {3, 5, 0},  {0, 4, 1},  {1, 2}};

static const int DatOldEdgeNewEdge[][QUAD2ConfMAXN_nEpoE] =
               { {3},  {4, 5},  {0, 1},  {2}};

static const int DatOldEdgeNewLocEdge[][QUADConfN_E] =
               { {-1, 3, 0, -1}, {-1, -1, 0, 1},
                 {2, 3, -1, -1} };

static const int DatNewEdgeOldEdge[] =
               { 2,  2,  3,  0,  1,  1,  -1,  -1,  -1};

static const int DatNewEdgeEqOldEdge[] = { 2, 3};
static const int DatNewEdgeEqOldEdgeIndex[] = { 3, 0};

  Type = Quad2Conf2;

  // set all numbers
  N_Edges = 9;
  N_Vertices = 7;
  N_Children = 3;
  N_NewVertEqOldVert = 4;
  N_NewEdgeEqOldEdge = 2;
  N_InnerVertices = 1;    
  N_InnerEdges = 3;

  // initialize all dimension values
  MaxN_VpC = QUAD2ConfMAXN_VpC;
  MaxN_CpV = QUAD2ConfMAXN_CpV;
  MaxN_EpC = QUAD2ConfMAXN_EpC;
  MaxN_CpE = QUAD2ConfMAXN_CpE;
  MaxN_EpV = QUAD2ConfMAXN_EpV;
  MaxN_iVpE = QUAD2ConfMAXN_iVpE;
  MaxN_nVpoE = QUAD2ConfMAXN_nVpoE;
  MaxN_nEpoE = QUAD2ConfMAXN_nEpoE;

  // initialize all pointers
  ChildType = (const Shapes *) DatChildType;
  EdgeType = (const Refinements *) DatEdgeType;

  ChildVertex = (const int *) DatChildVertex;
  VertexChild = (const int *) DatVertexChild;
  VertexChildIndex = (const int *) DatVertexChildIndex;
  VertexChildLen = (const int *) DatVertexChildLen;

  ChildEdge = (const int *) DatChildEdge;
  EdgeChild = (const int *) DatEdgeChild;
  EdgeChildIndex = (const int *) DatEdgeChildIndex;
  EdgeChildLen = (const int *) DatEdgeChildLen;

  EdgeVertex = (const int *) DatEdgeVertex;
  VertexEdge = (const int *) DatVertexEdge;
  VertexEdgeIndex = (const int *) DatVertexEdgeIndex;
  VertexEdgeLen = (const int *) DatVertexEdgeLen;

  NewVertexEqOldVertex = (const int *) DatNewVertexEqOldVertex;
  NewVertexEqOldVertexIndex = (const int *) DatNewVertexEqOldVertexIndex;

  InteriorVertexOfCell = (const int *) DatInteriorVertexOfCell;
  PositionOfIntVert = (const double *) DatPositionOfIntVert;

  InteriorEdgeOfCell = (const int *) DatInteriorEdgeOfCell;
  InteriorVertexOfEdge = (const int *) DatInteriorVertexOfEdge;
  InteriorVertexOfEdgeLen = (const int *) DatInteriorVertexOfEdgeLen;

  OldEdgeNewVertex = (const int *) DatOldEdgeNewVertex;

  OldEdgeNewEdge = (const int *) DatOldEdgeNewEdge;
  OldEdgeNewLocEdge = (const int *) DatOldEdgeNewLocEdge;
  NewEdgeOldEdge = (const int *) DatNewEdgeOldEdge;

  NewEdgeEqOldEdge = (const int *) DatNewEdgeEqOldEdge;
  NewEdgeEqOldEdgeIndex = (const int *) DatNewEdgeEqOldEdgeIndex;
}

// Methods
