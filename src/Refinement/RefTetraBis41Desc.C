#include <RefTetraBis41Desc.h>

// Constructor
TRefTetraBis41Desc::TRefTetraBis41Desc(const TShapeDesc *shape) : TRefDesc(shape)
{
static const Shapes DatChildType[] = {Tetrahedron, Tetrahedron, Tetrahedron};

static const Refinements DatEdgeType[] = {NoRef, LineReg, NoRef, NoRef, LineReg, NoRef};

static const Refinements DatFaceType [] = {TriBis1, TriBis1, TriBis10, NoRef};

/*
 * Vertex
 */
static const int DatChildVertex[][4] = {{0,1,5,4},{0,5,2,4},{0,4,2,3}};
static const int DatVertexChild[][3] = {{0,1,2},{0},{1,2},{2},{0,1,2},{0,1}};
static const int DatVertexChildIndex[][3] = {{0,0,0},{1},{2,2},{3},{3,3,1},{2,1}};
static const int DatVertexChildLen[] = {3,1,2,1,3,2};

/*
 * Faces
 */
static const int DatChildFace[][4] = {{0,2,4,9},{1,9,5,8},{8,3,6,7}};
static const int DatFaceChild[][2] = {{0},{1},{0},{2},{0},{1},{2},{2},{1,2},{0,1}};
static const int DatFaceChildIndex[][2] = {{0},{0},{1},{1},{2},{2},{2},{3},{3,0},{3,1}};
static const int DatFaceChildLen[] = {1,1,1,1,1,1,1,1,2,2};

/*
 * Edges
 */
static const int DatChildEdge[][6] = {{0,1,8,9,5,11},{8,2,3,9,11,10},{9,10,3,4,6,7}};
static const int DatEdgeChild[][3] = {{0},{0},{1},{1,2},{2},{0},{2},{2},{0,1},{0,1,2},{1,2},{0,1}};
static const int DatEdgeChildIndex[][3] = {{0},{1},{1},{2,2},{3},{4},{4},{5},{2,0},{3,3,0},{5,1},{5,4}};
static const int DatEdgeChildLen[] = {1,1,1,2,1,1,1,1,2,3,2,2};

/*
 * Edge-Vertex
 */
static const int DatEdgeVertex[][2] = {{0,1},{1,5},{5,2},{2,0},{0,3},{1,4},{4,3},{2,3},{5,0},{0,4},{2,4},{5,4}};
static const int DatVertexEdge[][5] = {{0,3,4,8,9},{0,1,5},{2,3,7,10},{4,6,7},{5,6,9,10,11},{1,2,8,11}};
static const int DatVertexEdgeIndex[][5] = {{0,1,0,1,0},{1,0,0},{1,0,0,0},{1,1,1},{1,0,1,1,1},{1,0,0,0}};
static const int DatVertexEdgeLen[] = {5,3,4,3,5,4};

/*
 * Face-Vertex
 */
static const int DatFaceVertex[][3] = {{0,1,5},{0,5,2},{0,4,1},{0,3,4},{5,1,4},{2,5,4},{2,4,3},{0,2,3},{0,2,4},{0,5,4}};
static const int DatVertexFace[][7] = {{0,1,2,3,7,8,9},{0,2,4},{1,5,6,7,8},{3,6,7},{2,3,4,5,6,8,9},{0,1,4,5,9}};
static const int DatVertexFaceIndex[][7] = {{0,0,0,0,0,0,0},{1,2,1},{2,0,0,1,1},{1,2,2},{1,2,2,2,1,2,2},{2,1,0,1,1}};
static const int DatVertexFaceLen[] = {7,3,5,3,7,5};

/*
 * Face-Edge
 */
static const int DatFaceEdge[][3] = {{0,1,8},{8,2,3},{9,5,0},{4,6,9},{1,5,11},{2,11,10},{10,6,7},{3,7,4},{3,10,9},{8,11,9}};
static const int DatEdgeFace[][4] = {{0,2},{0,4},{1,5},{1,7,8},{3,7},{2,4},{3,6},{6,7},{0,1,9},{2,3,8,9},{5,6,8},{4,5,9}};
static const int DatEdgeFaceIndex[][4] = {{0,2},{1,0},{1,0},{2,0,0},{0,2},{1,1},{1,1},{2,1},{2,0,0},{0,2,2,2},{2,0,1},{2,1,1}};
static const int DatEdgeFaceLen[] = {2,2,2,3,2,2,2,2,3,4,3,3};

static const int DatNewVertexEqOldVertex[] = {0, 1, 2, 3};
static const int DatNewVertexEqOldVertexIndex[] = {0, 1, 2, 3};

static const int DatNewFaceEqOldFace[] = {7};
static const int DatNewFaceEqOldFaceIndex[] = {3};

static const int DatOldFaceNewVertex[][REFTETRABIS41MAXN_nVpoF] =
  { {0, 1, 2, 5}, {0, 3, 1, 4}, {2, 1, 3, 4, 5}, {0, 2, 3} };
static const int DatOldFaceNewVertexLen[] =
  { 4, 4, 5, 3};
static const double DatOldFaceNewVertexPos[][REFTETRABIS41MAXN_nVpoF][REFTETRABIS41MAXN_oVpoF] =
  { { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0.5, 0.5} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0.5, 0.5} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0.5, 0.5}, {0.5, 0.5, 0} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} } };

static const int DatInteriorFaceOfCell[] = {8, 9};
static const int DatInteriorEdgeOfCell[] = {-1};

static const int DatInteriorVertexOfEdge[][REFTETRABIS41MAXN_iVpE] =
  { {}, {5}, {}, {}, {4}, {} };
static const int DatInteriorVertexOfEdgeLen[] = {0, 1, 0, 0, 1, 0};

static const int DatInteriorEdgeOfFace[][REFTETRABIS41MAXN_iEpF] =
  { {8}, {9}, {10, 11}, {} };
static const int DatInteriorEdgeOfFaceLen[] =
  { 1, 1, 2, 0};

static const int DatOldEdgeNewVertex[][REFTETRABIS41MAXN_nVpoE] =
  { {0, 1}, {1, 5, 2}, {2, 0}, {0, 3}, {1, 4, 3}, {2, 3} };

static const int DatOldEdgeNewVertexLen[] =
  {2, 3, 2, 2, 3, 2};

static const int DatOldEdgeNewEdge[][REFTETRABIS41MAXN_nEpoE] =
  { {0}, {1, 2}, {3}, {4}, {5, 6}, {7} };

static const int DatOldEdgeNewEdgeLen[] =
  {1, 2, 1, 1, 2, 1};

static const int DatOldFaceNewEdge[][REFTETRABIS41MAXN_nEpoF] =
  { {0, 1, 2, 3}, {4, 6, 5, 0}, {2, 1, 5, 6, 7}, {3, 7, 4} };

static const int DatOldFaceNewEdgeLen[] =
  {4, 4, 5, 3};

static const int DatOldFaceNewFace[][REFTETRABIS41MAXN_nFpoF] =
  { {0, 1}, {3, 2}, {5, 4, 6}, {7} };

static const int DatOldFaceNewFaceLen[] =
  {2, 2, 3, 1};

static const int DatNewEdgeOldEdge[] =
  {0, 1, 1, 2, 3, 4, 4, 5, -1, -1, -1, -1};

static const int DatNewFaceOldFace[] =
  {0, 0, 1, 1, 2, 2, 2, 3, -1, -1};

static const int DatOldFaceNewLocFace[][4] =
  { {0, 1, 2, -1}, {0, -1, 2, -1}, {-1, 1, 2, 3} };

static const int DatChildTwistIndex[] =
  {1, 2, 2, 1, 1, 0, 2, 0, -1, -1};

  Type = TetraBis41;

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
  MaxN_VpC = REFTETRABIS41MAXN_VpC;
  MaxN_CpV = REFTETRABIS41MAXN_CpV;
  MaxN_EpC = REFTETRABIS41MAXN_EpC;
  MaxN_CpE = REFTETRABIS41MAXN_CpE;
  MaxN_EpV = REFTETRABIS41MAXN_EpV;
  MaxN_EpF = REFTETRABIS41MAXN_EpF;
  MaxN_FpE = REFTETRABIS41MAXN_FpE;
  MaxN_VpF = REFTETRABIS41MAXN_VpF;
  MaxN_FpV = REFTETRABIS41MAXN_FpV;
  MaxN_FpC = REFTETRABIS41MAXN_FpC;
  MaxN_CpF = REFTETRABIS41MAXN_CpF;
  MaxN_iVpE = REFTETRABIS41MAXN_iVpE;
  MaxN_iEpF = REFTETRABIS41MAXN_iEpF;
  MaxN_nVpoE = REFTETRABIS41MAXN_nVpoE;
  MaxN_nEpoE = REFTETRABIS41MAXN_nEpoE;
  MaxN_nVpoF = REFTETRABIS41MAXN_nVpoF;
  MaxN_oVpoF = REFTETRABIS41MAXN_oVpoF;
  MaxN_nEpoF = REFTETRABIS41MAXN_nEpoF;
  MaxN_nFpoF = REFTETRABIS41MAXN_nFpoF;

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





