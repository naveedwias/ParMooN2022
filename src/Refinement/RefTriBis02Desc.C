// =======================================================================
// @(#)RefTriBis0Desc.C        1.3 09/13/99
//
// Class:       TRefTriBis0Desc
// Purpose:     bisection of edge 0
//              into two triangles
//
// Author:      Volker Behns  16.03.98
//              Gunar Matthies 03.09.1998
//
// =======================================================================

#include <RefTriBis02Desc.h>

// Constructor
TRefTriBis02Desc::TRefTriBis02Desc(const TShapeDesc *shape) : TRefDesc(shape)
{
static const Shapes DatChildType[] = { Triangle, Triangle, Triangle};

static const Refinements DatEdgeType[] = { LineReg, NoRef, LineReg};

static const int DatChildVertex[][TRIBI02MAXN_VpC] = {{0,3,4},{2,4,3},{1,2,3}};
static const int DatVertexChild[][TRIBI02MAXN_CpV] = {{0},{2},{1,2},{0,1,2},{0,1}};
static const int DatVertexChildIndex[][TRIBI02MAXN_CpV] = {{0},{0},{0,1},{1,2,2},{2,1}};
static const int DatVertexChildLen[] = {1,1,2,3,2};

static const int DatChildEdge[][TRIBI02MAXN_EpC] = {{0,6,4},{3,6,5},{2,5,1}};
static const int DatEdgeChild[][TRIBI02MAXN_CpE] = {{0},{2},{2},{1},{0},{1,2},{0,1}};
static const int DatEdgeChildIndex[][TRIBI02MAXN_CpE] = {{0},{2},{0},{0},{2},{2,1},{1,1}};
static const int DatEdgeChildLen[] = {1,1,1,1,1,2,2};

static const int DatEdgeVertex[][2] = {{0,3},{3,1},{1,2},{2,4},{4,0},{3,2},{3,4}};
static const int DatVertexEdge[][TRIBI02MAXN_EpV] = {{0,4},{1,2},{2,3,5},{0,1,5,6},{3,4,6}};
static const int DatVertexEdgeIndex[][TRIBI02MAXN_EpV] = {{0,1},{1,0},{1,0,1},{1,0,0,0},{1,0,1}};
static const int DatVertexEdgeLen[] = {2,2,3,4,3};

static const int DatNewVertexEqOldVertex[] = { 0, 1, 2};
static const int DatNewVertexEqOldVertexIndex[] = { 0, 1, 2};

static const int DatNewEdgeEqOldEdge[] = {2};
static const int DatNewEdgeEqOldEdgeIndex[] = {1};

static const int DatInteriorVertexOfEdge[][TRIBI02MAXN_iVpE] = { {3}, {}, {4} };
static const int DatInteriorVertexOfEdgeLen[] = { 1, 0, 1};

static const int DatInteriorEdgeOfCell[] = {5, 6};

static const int DatOldEdgeNewVertex[][TRIBI02MAXN_nVpoE] = { {0, 3, 1}, { 1, 2}, {2, 4, 0}};

static const int DatOldEdgeNewEdge[][TRIBI02MAXN_nEpoE] =
               { {0, 1}, {2}, {3, 4}};

static const int DatOldEdgeNewLocEdge[][TRIBI02N_E] =
               { {0, -1, 2}, {-1, -1, 0}, {2, 0, -1} };

static const int DatNewEdgeOldEdge[] =
               { 0, 0, 1, 2, 2, -1, -1};

  Type = TriBis02;

  // set all numbers
  N_Edges = 7;
  N_Vertices = 5;
  N_Children = 3;
  N_NewVertEqOldVert = 3;
  N_NewEdgeEqOldEdge = 1;
  N_InnerEdges = 2;

  // initialize all dimension values
  MaxN_VpC = TRIBI02MAXN_VpC;
  MaxN_CpV = TRIBI02MAXN_CpV;
  MaxN_EpC = TRIBI02MAXN_EpC;
  MaxN_CpE = TRIBI02MAXN_CpE;
  MaxN_EpV = TRIBI02MAXN_EpV;
  MaxN_nVpoE = TRIBI02MAXN_nVpoE;
  MaxN_nEpoE = TRIBI02MAXN_nEpoE;

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

  InteriorVertexOfEdge = (const int *) DatInteriorVertexOfEdge;
  InteriorVertexOfEdgeLen = (const int *) DatInteriorVertexOfEdgeLen;
}

// Methods
