// =======================================================================
// @(#)RefQuadToTri1Desc.C        1.4 10/19/99
//
// Class:       TRefQuadToTri1Desc
// Purpose:     refinement descriptor for refinement of a quadrangle
//              into two triangles
//
// Author:      Volker Behns  16.03.98
//
// =======================================================================

#include <RefQuadToTri1Desc.h>

// Constructor
TRefQuadToTri1Desc::TRefQuadToTri1Desc(const TShapeDesc *shape) : TRefDesc(shape)
{
static const Shapes DatChildType[] = { Triangle, Triangle};

static const Refinements DatEdgeType[] = { NoRef, NoRef, NoRef, NoRef};

static const int DatChildVertex[][QUADTTMAXN_VpC] =
                 { {0, 1, 3},  {1, 2, 3}};
static const int DatVertexChild[][QUADTTMAXN_CpV] =
                 { {0},  {0, 1},  {1},  {0, 1}};
static const int DatVertexChildIndex[][QUADTTMAXN_CpV] =
               { {0},  {1, 0},  {1},  {2, 2}};
static const int DatVertexChildLen[] =
               { 1,  2,  1,  2};

static const int DatChildEdge[][QUADTTMAXN_EpC] =
               { {0, 4, 3},  {1, 2, 4}};
static const int DatEdgeChild[][QUADTTMAXN_CpE] =
               { {0},  {1},  {1},  {0},  {0, 1}};
static const int DatEdgeChildIndex[][QUADTTMAXN_CpE] = 
               { {0},  {0},  {1},  {2},  {1, 2}};
static const int DatEdgeChildLen[] =
               { 1,  1,  1,  1,  2};

static const int DatEdgeVertex[][2] =
               { {0, 1},  {1, 2},  {2, 3},  {3, 0},  {0, 2}};
static const int DatVertexEdge[][QUADTTMAXN_EpV] = 
               { {0, 3},  {0, 1, 4},  {1, 2},  {2, 3, 4}};
static const int DatVertexEdgeIndex[][QUADTTMAXN_EpV] =
               { {0, 1},  {1, 0, 0},  {1, 0},  {1, 0, 1}};
static const int DatVertexEdgeLen[] =
               { 2,  3,  2,  3};

static const int DatNewVertexEqOldVertex[] = { 0, 1, 2, 3};
static const int DatNewVertexEqOldVertexIndex[] = { 0, 1, 2, 3};

static const int DatNewEdgeEqOldEdge[] = { 0, 1, 2, 3};
static const int DatNewEdgeEqOldEdgeIndex[] = { 0, 1, 2, 3};

static const int DatInteriorEdgeOfCell[] = {4};

static const int DatInteriorVertexOfEdgeLen[] = { 0,  0,  0,  0};

static const int DatOldEdgeNewVertex[][QUADTTMAXN_nVpoE] =
               { {0, 1},  {1, 2},  {2, 3},  {3, 0}};

static const int DatOldEdgeNewEdge[][QUADTTMAXN_nEpoE] =
               { {0},  {1},  {2},  {3}};

// static const int DatOldEdgeNewLocEdge[][QUADTTN_E] =
//                { {0, -1, 3}, {1, 2, -1}};

static const int DatOldEdgeNewLocEdge[][4] =
               { {  0, -1, -1,  2 },
                 { -1,  0,  1, -1 } };


static const int DatNewEdgeOldEdge[] =
               { 0,  1,  2,  3,  -1};

  Type = QuadToTri1;

  // set all numbers
  N_Edges = 5;
  N_Vertices = 4;
  N_Children = 2;
  N_NewVertEqOldVert = 4;
  N_NewEdgeEqOldEdge = 4;
  N_InnerEdges = 1;

  // initialize all dimension values
  MaxN_VpC = QUADTTMAXN_VpC;
  MaxN_CpV = QUADTTMAXN_CpV;
  MaxN_EpC = QUADTTMAXN_EpC;
  MaxN_CpE = QUADTTMAXN_CpE;
  MaxN_EpV = QUADTTMAXN_EpV;
  MaxN_nVpoE = QUADTTMAXN_nVpoE;
  MaxN_nEpoE = QUADTTMAXN_nEpoE;

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

  NewEdgeEqOldEdge = (const int *) DatNewEdgeEqOldEdge;
  NewEdgeEqOldEdgeIndex = (const int *) DatNewEdgeEqOldEdgeIndex;

  InteriorEdgeOfCell = (const int *) DatInteriorEdgeOfCell;

  OldEdgeNewVertex = (const int *) DatOldEdgeNewVertex;

  OldEdgeNewEdge = (const int *) DatOldEdgeNewEdge;
  OldEdgeNewLocEdge = (const int *) DatOldEdgeNewLocEdge;
  NewEdgeOldEdge = (const int *) DatNewEdgeOldEdge;

  InteriorVertexOfEdgeLen = (const int *) DatInteriorVertexOfEdgeLen;
}

// Methods
