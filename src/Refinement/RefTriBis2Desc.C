// =======================================================================
// @(#)RefTriBis2Desc.C        1.3 09/13/99
//
// Class:       TRefTriBis2Desc
// Purpose:     bisection of edge 0
//              into two triangles
//
// Author:      Volker Behns  16.03.98
//              Gunar Matthies 03.09.1998
//
// =======================================================================

#include <RefTriBis2Desc.h>

// Constructor
TRefTriBis2Desc::TRefTriBis2Desc(const TShapeDesc *shape) : TRefDesc(shape)
{
static const Shapes DatChildType[] = { Triangle, Triangle};

static const Refinements DatEdgeType[] = { NoRef, NoRef, LineReg };

static const int DatChildVertex[][TRIBI2MAXN_VpC] =
                 { {2, 3, 1}, {0, 1, 3}};
static const int DatVertexChild[][TRIBI2MAXN_CpV] =
                 { {1}, {0, 1}, {0}, {0, 1} };
static const int DatVertexChildIndex[][TRIBI2MAXN_CpV] =
               { {0}, {2, 1}, {0}, {1, 2} };
static const int DatVertexChildLen[] =
               { 1, 2, 1, 2 };

static const int DatChildEdge[][TRIBI2MAXN_EpC] =
               { {1, 4, 0}, {3, 4, 2}};
static const int DatEdgeChild[][TRIBI2MAXN_CpE] =
               { {0},  {0},  {1},  {1},  {0, 1}};
static const int DatEdgeChildIndex[][TRIBI2MAXN_CpE] = 
               { {2},  {0},  {2},  {0},  {1, 1}};
static const int DatEdgeChildLen[] =
               { 1,  1,  1,  1,  2};

static const int DatEdgeVertex[][2] =
               { {1, 2}, {2, 3}, {3, 0}, {0, 1}, {1, 3}};
static const int DatVertexEdge[][TRIBI2MAXN_EpV] = 
               { {2, 3}, {0, 3, 4}, {0, 1}, {1, 2, 4}};
               
static const int DatVertexEdgeIndex[][TRIBI2MAXN_EpV] =
               { {0, 1}, {0, 1, 0}, {1, 0}, {1, 0, 1}};
static const int DatVertexEdgeLen[] =
               { 2, 3, 3, 2};

static const int DatNewVertexEqOldVertex[] = { 0, 1, 2};
static const int DatNewVertexEqOldVertexIndex[] = { 0, 1, 2};

static const int DatNewEdgeEqOldEdge[] = { 0, 3};
static const int DatNewEdgeEqOldEdgeIndex[] = { 1, 0};

static const int DatInteriorVertexOfEdge[][TRIBI2MAXN_iVpE] =
                { {-1}, {-1}, {3} };
static const int DatInteriorVertexOfEdgeLen[] = { 0, 0, 1};

static const int DatInteriorEdgeOfCell[] = { 4};

static const int DatOldEdgeNewVertex[][TRIBI2MAXN_nVpoE] =
               { {0, 1}, {1, 2}, {2, 3, 0}};

static const int DatOldEdgeNewEdge[][TRIBI2MAXN_nEpoE] =
               { {3}, {0}, {1, 2}};

static const int DatOldEdgeNewLocEdge[][TRIBI2N_E] =
               { {-1, 2, 0}, { 0, -1, 2}};

static const int DatNewEdgeOldEdge[] =
               { 1, 2, 2, 0, -1};

  Type = TriBis2;

  // set all numbers
  N_Edges = 5;
  N_Vertices = 4;
  N_Children = 2;
  N_NewVertEqOldVert = 3;
  N_NewEdgeEqOldEdge = 2;
  N_InnerEdges = 1;

  // initialize all dimension values
  MaxN_VpC = TRIBI2MAXN_VpC;
  MaxN_CpV = TRIBI2MAXN_CpV;
  MaxN_EpC = TRIBI2MAXN_EpC;
  MaxN_CpE = TRIBI2MAXN_CpE;
  MaxN_EpV = TRIBI2MAXN_EpV;
  MaxN_nVpoE = TRIBI2MAXN_nVpoE;
  MaxN_nEpoE = TRIBI2MAXN_nEpoE;

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
