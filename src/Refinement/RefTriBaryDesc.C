#include <RefTriBaryDesc.h>

// Constructor
TRefTriBaryDesc::TRefTriBaryDesc(const TShapeDesc *shape) : TRefDesc(shape)
{
static const Shapes DatChildType[] = { Triangle, Triangle, Triangle};

static const Refinements DatEdgeType[] = { NoRef, NoRef, NoRef};

static const int DatChildVertex[][3] = {{3, 0, 1}, {3, 1, 2}, {3, 2, 0}};
static const int DatVertexChild[][3] = {{2, 0}, {0, 1}, {1, 2}, {2, 1, 0}};
static const int DatVertexChildIndex[][3] = {{2, 1}, {2, 1}, {2, 1}, {0, 0, 0}};
static const int DatVertexChildLen[] = {2, 2, 2, 3};

static const int DatChildEdge[][3] = {{3,0,4}, {4,1,5}, {5,2,3}};
static const int DatEdgeChild[][2] =
                {{0}, {1}, {2}, {0,2}, {0,1}, {1,2}};
static const int DatEdgeChildIndex[][2] = 
                 {{1}, {1}, {1}, {0,2}, {2,0}, {2,0}};
static const int DatEdgeChildLen[] = {1, 1, 1, 2, 2, 2};

static const int DatEdgeVertex[][2] =
                 {{0,1}, {1,2}, {2,0}, {0,3}, {1,3}, {2,3}};
static const int DatVertexEdge[][3] = 
                 {{0,2,3}, {0,1,4}, {1,2,5}, {3,4,5}};
static const int DatVertexEdgeIndex[][3] =
                 {{0,1,0}, {1,0,0}, {1,0,0}, {1,1,1}};
static const int DatVertexEdgeLen[] = {3, 3, 3, 3};

static const int DatNewVertexEqOldVertex[] = { 0, 1, 2};
static const int DatNewVertexEqOldVertexIndex[] = { 0, 1, 2};

static const int DatNewEdgeEqOldEdge[] = {0, 1, 2};
static const int DatNewEdgeEqOldEdgeIndex[] = {0, 1, 2};

static const int DatInteriorVertexOfCell[] = {3};
static const double DatPositionOfIntVert[][3] =
                  { {1./3., 1./3., 1./3.}};

static const int DatInteriorEdgeOfCell[] = { 3, 4, 5};

static const int DatInteriorVertexOfEdgeLen[] = { 0,  0,  0};

static const int DatOldEdgeNewVertex[][2] =
                 {{0, 1}, {1, 2}, {2, 0}};

static const int DatOldEdgeNewEdge[][1] = {{0}, {1}, {2}};

static const int DatOldEdgeNewLocEdge[][3] =
                 {{-1, 0, -1}, {-1, 1, -1}, {-1, 2, -1}};

static const int DatNewEdgeOldEdge[] =
                 {0, 1, 2, -1, -1, -1};

  Type = TriBary;

  // set all numbers
  N_Edges = 6;
  N_Vertices = 4;
  N_Children = 3;
  N_NewVertEqOldVert = 3;
  N_NewEdgeEqOldEdge = 3;
  N_InnerVertices = 1;    
  N_InnerEdges = 3;

  // initialize all dimension values
  MaxN_VpC = 3;
  MaxN_CpV = 3;
  MaxN_EpC = 3;
  MaxN_CpE = 2;
  MaxN_EpV = 3;
  MaxN_iVpE = 0;
  MaxN_nVpoE = 2;
  MaxN_nEpoE = 1;

  // initialize all pointers
  ChildType = (const Shapes *) DatChildType;
  EdgeType = (const Refinements *) DatEdgeType;

  ChildVertex = (const int *) &(DatChildVertex[0][0]);
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
  
  NewEdgeEqOldEdge = (const int *) DatNewEdgeEqOldEdge;
  NewEdgeEqOldEdgeIndex = (const int *) DatNewEdgeEqOldEdgeIndex;

  InteriorEdgeOfCell = (const int *) DatInteriorEdgeOfCell;
  InteriorVertexOfEdge = (const int *) nullptr;
  InteriorVertexOfEdgeLen = (const int *) DatInteriorVertexOfEdgeLen;

  OldEdgeNewVertex = (const int *) DatOldEdgeNewVertex;

  OldEdgeNewEdge = (const int *) DatOldEdgeNewEdge;
  OldEdgeNewLocEdge= (const int *) DatOldEdgeNewLocEdge;
  NewEdgeOldEdge = (const int *) DatNewEdgeOldEdge;
}


// Methods
