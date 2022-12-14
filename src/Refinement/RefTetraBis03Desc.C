#include <RefTetraBis03Desc.h>

// Constructor
TRefTetraBis03Desc::TRefTetraBis03Desc(const TShapeDesc *shape) : TRefDesc(shape)
{
static const Shapes DatChildType[] = {Tetrahedron, Tetrahedron, Tetrahedron};

static const Refinements DatEdgeType[] = {LineReg, NoRef, NoRef, LineReg, NoRef, NoRef};

static const Refinements DatFaceType [] = {TriBis0, TriBis20, NoRef, TriBis2};

/*
 * Vertex
 */
static const int DatChildVertex[][4] = {{0,4,2,5},{4,1,2,3},{5,4,2,3}};
static const int DatVertexChild[][3] = {{0},{1},{0,1,2},{1,2},{0,1,2},{0,2}};
static const int DatVertexChildIndex[][3] = {{0},{1},{2,2,2},{3,3},{1,0,1},{3,0}};
static const int DatVertexChildLen[] = {1,1,3,2,3,2};

/*
 * Faces
 */
static const int DatChildFace[][4] = {{0,2,9,6},{1,3,5,8},{9,4,8,7}};
static const int DatFaceChild[][2] = {{0},{1},{0},{1},{2},{1},{0},{2},{1,2},{0,2}};
static const int DatFaceChildIndex[][2] = {{0},{0},{1},{1},{1},{2},{3},{3},{3,2},{2,0}};
static const int DatFaceChildLen[] = {1,1,1,1,1,1,1,1,2,2};

/*
 * Edges
 */
static const int DatChildEdge[][6] = {{0,8,3,4,10,11},{1,2,8,9,6,7},{10,8,11,5,9,7}};
static const int DatEdgeChild[][3] = {{0},{1},{1},{0},{0},{2},{1},{1,2},{0,1,2},{1,2},{0,2},{0,2}};
static const int DatEdgeChildIndex[][3] = {{0},{0},{1},{2},{3},{3},{4},{5,5},{1,2,1},{3,4},{4,0},{5,2}};
static const int DatEdgeChildLen[] = {1,1,1,1,1,1,1,2,3,2,2,2};

/*
 * Edge-Vertex
 */
static const int DatEdgeVertex[][2] = {{0,4},{4,1},{1,2},{2,0},{0,5},{5,3},{1,3},{2,3},{4,2},{4,3},{4,5},{2,5}};
static const int DatVertexEdge[][5] = {{0,3,4},{1,2,6},{2,3,7,8,11},{5,6,7,9},{0,1,8,9,10},{4,5,10,11}};
static const int DatVertexEdgeIndex[][5] = {{0,1,0},{1,0,0},{1,0,0,1,0},{1,1,1,1},{1,0,0,0,0},{1,0,1,1}};
static const int DatVertexEdgeLen[] = {3,3,5,4,5,4};

/*
 * Face-Vertex
 */
static const int DatFaceVertex[][3] = {{0,4,2},{4,1,2},{0,5,4},{4,3,1},{5,3,4},{2,1,3},{0,2,5},{5,2,3},{4,2,3},{2,4,5}};
static const int DatVertexFace[][7] = {{0,2,6},{1,3,5},{0,1,5,6,7,8,9},{3,4,5,7,8},{0,1,2,3,4,8,9},{2,4,6,7,9}};
static const int DatVertexFaceIndex[][7] = {{0,0,0},{1,2,1},{2,2,0,1,1,1,0},{1,1,2,2,2},{1,0,2,0,2,0,1},{1,0,2,0,2}};
static const int DatVertexFaceLen[] = {3,3,7,5,7,5};

/*
 * Face-Edge
 */
static const int DatFaceEdge[][3] = {{0,8,3},{1,2,8},{4,10,0},{9,6,1},{5,9,10},{2,6,7},{3,11,4},{11,7,5},{8,7,9},{8,10,11}};
static const int DatEdgeFace[][4] = {{0,2},{1,3},{1,5},{0,6},{2,6},{4,7},{3,5},{5,7,8},{0,1,8,9},{3,4,8},{2,4,9},{6,7,9}};
static const int DatEdgeFaceIndex[][4] = {{0,2},{0,2},{1,0},{2,0},{0,2},{0,2},{1,1},{2,1,1},{1,2,0,0},{0,1,2},{1,2,1},{1,0,2}};
static const int DatEdgeFaceLen[] = {2,2,2,2,2,2,2,3,4,3,3,3};


static const int DatNewVertexEqOldVertex[] = {0, 1, 2, 3};
static const int DatNewVertexEqOldVertexIndex[] = {0, 1, 2, 3};

static const int DatNewFaceEqOldFace[] = {5};
static const int DatNewFaceEqOldFaceIndex[] = {2};

static const int DatOldFaceNewVertex[][REFTETRABIS03MAXN_nVpoF] =
  { {0, 1, 2, 4}, {0, 3, 1, 4, 5}, {2, 1, 3}, {0, 2, 3, 5} };
static const int DatOldFaceNewVertexLen[] =
  { 4, 5, 3, 4};
static const double DatOldFaceNewVertexPos[][REFTETRABIS03MAXN_nVpoF][REFTETRABIS03MAXN_oVpoF] =
  { { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0.5, 0} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0, 0.5}, {0.5, 0.5, 0} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0, 0.5} } };

static const int DatInteriorFaceOfCell[] = {8, 9};
static const int DatInteriorEdgeOfCell[] = {-1};

static const int DatInteriorVertexOfEdge[][REFTETRABIS03MAXN_iVpE] =
  { {4}, {}, {}, {5}, {}, {} };
static const int DatInteriorVertexOfEdgeLen[] = {1, 0, 0, 1, 0, 0};

static const int DatInteriorEdgeOfFace[][REFTETRABIS03MAXN_iEpF] =
  { {8}, {9, 10}, {}, {11} };
static const int DatInteriorEdgeOfFaceLen[] =
  { 1, 2, 0, 1};

static const int DatOldEdgeNewVertex[][REFTETRABIS03MAXN_nVpoE] =
  { {0, 4, 1}, {1, 2}, {2, 0}, {0, 5, 3}, {1, 3}, {2, 3} };

static const int DatOldEdgeNewVertexLen[] =
  {3, 2, 2, 3, 2, 2};

static const int DatOldEdgeNewEdge[][REFTETRABIS03MAXN_nEpoE] =
  { {0, 1}, {2}, {3}, {4, 5}, {6}, {7} };

static const int DatOldEdgeNewEdgeLen[] =
  {2, 1, 1, 2, 1, 1};

static const int DatOldFaceNewEdge[][REFTETRABIS03MAXN_nEpoF] =
  { {0, 1, 2, 3}, {4, 5, 6, 1, 0}, {2, 6, 7}, {3, 7, 5, 4} };

static const int DatOldFaceNewEdgeLen[] =
  {4, 5, 3, 4};

static const int DatOldFaceNewFace[][REFTETRABIS03MAXN_nFpoF] =
  { {0, 1}, {2, 4, 3}, {5}, {7, 6} };

static const int DatOldFaceNewFaceLen[] =
  {2, 3, 1, 2};

static const int DatNewEdgeOldEdge[] =
  {0, 0, 1, 2, 3, 3, 4, 5, -1, -1, -1, -1};

static const int DatNewFaceOldFace[] =
  {0, 0, 1, 1, 1, 2, 3, 3, -1, -1};

static const int DatOldFaceNewLocFace[][4] =
  { {0, 1, -1, 3}, {0, 1, 2, -1}, {-1, 1, -1, 3} };

static const int DatChildTwistIndex[] =
  {0, 1, 0, 2, 1, 0, 2, -1, -1};

  Type = TetraBis03;

  //set all numbers
  N_Vertices = 6;
  N_Edges = 12;
  N_Faces = 10;
  N_Children = 3;
  N_NewVertEqOldVert = 4;
  N_NewFaceEqOldFace = 1;
  N_InnerEdges = 0;
  N_InnerFaces = 2;

  // initialize all dimension values
  MaxN_VpC = REFTETRABIS03MAXN_VpC;
  MaxN_CpV = REFTETRABIS03MAXN_CpV;
  MaxN_EpC = REFTETRABIS03MAXN_EpC;
  MaxN_CpE = REFTETRABIS03MAXN_CpE;
  MaxN_EpV = REFTETRABIS03MAXN_EpV;
  MaxN_EpF = REFTETRABIS03MAXN_EpF;
  MaxN_FpE = REFTETRABIS03MAXN_FpE;
  MaxN_VpF = REFTETRABIS03MAXN_VpF;
  MaxN_FpV = REFTETRABIS03MAXN_FpV;
  MaxN_FpC = REFTETRABIS03MAXN_FpC;
  MaxN_CpF = REFTETRABIS03MAXN_CpF;
  MaxN_iVpE = REFTETRABIS03MAXN_iVpE;
  MaxN_iEpF = REFTETRABIS03MAXN_iEpF;
  MaxN_nVpoE = REFTETRABIS03MAXN_nVpoE;
  MaxN_nEpoE = REFTETRABIS03MAXN_nEpoE;
  MaxN_nVpoF = REFTETRABIS03MAXN_nVpoF;
  MaxN_oVpoF = REFTETRABIS03MAXN_oVpoF;
  MaxN_nEpoF = REFTETRABIS03MAXN_nEpoF;
  MaxN_nFpoF = REFTETRABIS03MAXN_nFpoF;

  // initialize all pointers
  ChildType = (const Shapes *) DatChildType;
  EdgeType = (const Refinements *) DatEdgeType;
  FaceType = (const Refinements *) DatFaceType;

  ChildVertex = (const int *) DatChildVertex;
  VertexChild = (const int *) DatVertexChild;
  VertexChildIndex = (const int *) DatVertexChildIndex;
  VertexChildLen = (const int *) DatVertexChildLen;

  ChildFace = (const int *) DatChildFace;
  FaceChild = (const int *) DatFaceChild;
  FaceChildIndex = (const int *) DatFaceChildIndex;
  FaceChildLen = (const int *) DatFaceChildLen;

  ChildEdge = (const int *) DatChildEdge;
  EdgeChild = (const int *) DatEdgeChild;
  EdgeChildIndex = (const int *) DatEdgeChildIndex;
  EdgeChildLen = (const int *) DatEdgeChildLen;

  EdgeVertex = (const int *) DatEdgeVertex;
  VertexEdge = (const int *) DatVertexEdge;
  VertexEdgeIndex = (const int *) DatVertexEdgeIndex;
  VertexEdgeLen = (const int *) DatVertexEdgeLen;

  FaceVertex = (const int *) DatFaceVertex;
  VertexFace = (const int *) DatVertexFace;
  VertexFaceIndex = (const int *) DatVertexFaceIndex;
  VertexFaceLen = (const int *) DatVertexFaceLen;

  FaceEdge = (const int *) DatFaceEdge;
  EdgeFace = (const int *) DatEdgeFace;
  EdgeFaceIndex = (const int *) DatEdgeFaceIndex;
  EdgeFaceLen = (const int *) DatEdgeFaceLen;

  NewVertexEqOldVertex = (const int *) DatNewVertexEqOldVertex;
  NewVertexEqOldVertexIndex = (const int *) DatNewVertexEqOldVertexIndex;

  NewFaceEqOldFace = (const int *) DatNewFaceEqOldFace;
  NewFaceEqOldFaceIndex = (const int *) DatNewFaceEqOldFaceIndex;

  InteriorEdgeOfCell = (const int *) DatInteriorEdgeOfCell;
  InteriorFaceOfCell = (const int *) DatInteriorFaceOfCell;
  InteriorVertexOfEdge = (const int *) DatInteriorVertexOfEdge;
  InteriorVertexOfEdgeLen = (const int *) DatInteriorVertexOfEdgeLen;
  InteriorEdgeOfFace = (const int *) DatInteriorEdgeOfFace;
  InteriorEdgeOfFaceLen = (const int *) DatInteriorEdgeOfFaceLen;

  OldEdgeNewVertex = (const int *) DatOldEdgeNewVertex;
  OldEdgeNewVertexLen = (const int *) DatOldEdgeNewVertexLen;

  OldEdgeNewEdge = (const int *) DatOldEdgeNewEdge;
  OldEdgeNewEdgeLen = (const int *) DatOldEdgeNewEdgeLen;
  NewEdgeOldEdge = (const int *) DatNewEdgeOldEdge;
  NewFaceOldFace = (const int *) DatNewFaceOldFace;

  OldFaceNewVertex = (const int *) DatOldFaceNewVertex;
  OldFaceNewVertexPos = (const double *) DatOldFaceNewVertexPos;
  OldFaceNewVertexLen = (const int *) DatOldFaceNewVertexLen;
  OldFaceNewEdge = (const int *) DatOldFaceNewEdge;
  OldFaceNewEdgeLen = (const int *) DatOldFaceNewEdgeLen;
  OldFaceNewFace = (const int *) DatOldFaceNewFace;
  OldFaceNewFaceLen = (const int *) DatOldFaceNewFaceLen;
  OldFaceNewLocFace = (const int *) DatOldFaceNewLocFace;
  ChildTwistIndex = (const int *) DatChildTwistIndex;
}

// Methods

