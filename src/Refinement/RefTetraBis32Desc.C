#include <RefTetraBis32Desc.h>

// Constructor
TRefTetraBis32Desc::TRefTetraBis32Desc(const TShapeDesc *shape) : TRefDesc(shape)
{
static const Shapes DatChildType[] = {Tetrahedron, Tetrahedron, Tetrahedron};

static const Refinements DatEdgeType[] = {NoRef, NoRef, LineReg, LineReg, NoRef, NoRef};

static const Refinements DatFaceType [] = {TriBis2, TriBis0, NoRef, TriBis20};

/*
 * Vertex
 */
static const int DatChildVertex[][4] = {{0,1,5,4},{5,1,2,4},{4,1,2,3}};
static const int DatVertexChild[][3] = {{0},{0,1,2},{1,2},{2},{0,1,2},{0,1}};
static const int DatVertexChildIndex[][3] = {{0},{1,1,1},{2,2},{3},{3,3,0},{2,0}};
static const int DatVertexChildLen[] = {1,3,2,1,3,2};

/*
 * Faces
 */
static const int DatChildFace[][4] = {{0,2,9,7},{1,9,8,6},{8,3,4,5}};
static const int DatFaceChild[][2] = {{0},{1},{0},{2},{2},{2},{1},{0},{1,2},{0,1}};
static const int DatFaceChildIndex[][2] = {{0},{0},{1},{1},{2},{3},{3},{3},{2,0},{2,1}};
static const int DatFaceChildLen[] = {1,1,1,1,1,1,1,1,2,2};

/*
 * Edges
 */
static const int DatChildEdge[][6] = {{0,8,3,4,9,11},{8,1,2,11,9,10},{9,1,10,5,6,7}};
static const int DatEdgeChild[][3] = {{0},{1,2},{1},{0},{0},{2},{2},{2},{0,1},{0,1,2},{1,2},{0,1}};
static const int DatEdgeChildIndex[][3] = {{0},{1,1},{2},{2},{3},{3},{4},{5},{1,0},{4,4,0},{5,2},{5,3}};
static const int DatEdgeChildLen[] = {1,2,1,1,1,1,1,1,2,3,2,2};

/*
 * Edge-Vertex
 */
static const int DatEdgeVertex[][2] = {{0,1},{1,2},{2,5},{5,0},{0,4},{4,3},{1,3},{2,3},{1,5},{1,4},{2,4},{5,4}};
static const int DatVertexEdge[][5] = {{0,3,4},{0,1,6,8,9},{1,2,7,10},{5,6,7},{4,5,9,10,11},{2,3,8,11}};
static const int DatVertexEdgeIndex[][5] = {{0,1,0},{1,0,0,0,0},{1,0,0,0},{1,1,1},{1,0,1,1,1},{1,0,1,0}};
static const int DatVertexEdgeLen[] = {3,5,4,3,5,4};

/*
 * Face-Vertex
 */
static const int DatFaceVertex[][3] = {{0,1,5},{5,1,2},{0,4,1},{4,3,1},{2,1,3},{4,2,3},{5,2,4},{0,5,4},{2,1,4},{5,1,4}};
static const int DatVertexFace[][7] = {{0,2,7},{0,1,2,3,4,8,9},{1,4,5,6,8},{3,4,5},{2,3,5,6,7,8,9},{0,1,6,7,9}};
static const int DatVertexFaceIndex[][7] = {{0,0,0},{1,1,2,2,1,1,1},{2,0,1,1,0},{1,2,2},{1,0,0,2,2,2,2},{2,0,0,1,0}};
static const int DatVertexFaceLen[] = {3,7,5,3,7,5};

/*
 * Face-Edge
 */
static const int DatFaceEdge[][3] = {{0,8,3},{8,1,2},{4,9,0},{5,6,9},{1,6,7},{10,7,5},{2,10,11},{3,11,4},{1,9,10},{8,9,11}};
static const int DatEdgeFace[][4] = {{0,2},{1,4,8},{1,6},{0,7},{2,7},{3,5},{3,4},{4,5},{0,1,9},{2,3,8,9},{5,6,8},{6,7,9}};
static const int DatEdgeFaceIndex[][4] = {{0,2},{1,0,0},{2,0},{2,0},{0,2},{0,2},{1,1},{2,1},{1,0,0},{1,2,1,1},{0,1,2},{2,1,2}};
static const int DatEdgeFaceLen[] = {2,3,2,2,2,2,2,2,3,4,3,3};

static const int DatNewVertexEqOldVertex[] = {0, 1, 2, 3};
static const int DatNewVertexEqOldVertexIndex[] = {0, 1, 2, 3};

static const int DatNewFaceEqOldFace[] = {4};
static const int DatNewFaceEqOldFaceIndex[] = {2};

static const int DatOldFaceNewVertex[][REFTETRABIS32MAXN_nVpoF] =
  { {0, 1, 2, 5}, {0, 3, 1, 4}, {2, 1, 3}, {0, 2, 3, 4, 5} };
static const int DatOldFaceNewVertexLen[] =
  { 4, 4, 3, 5};
static const double DatOldFaceNewVertexPos[][REFTETRABIS32MAXN_nVpoF][REFTETRABIS32MAXN_oVpoF] =
  { { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0, 0.5} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0.5, 0} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0, 0.5}, {0.5, 0.5, 0} } };

static const int DatInteriorFaceOfCell[] = {8, 9};
static const int DatInteriorEdgeOfCell[] = {-1};

static const int DatInteriorVertexOfEdge[][REFTETRABIS32MAXN_iVpE] =
  { {}, {}, {5}, {4}, {}, {} };
static const int DatInteriorVertexOfEdgeLen[] = {0, 0, 1, 1, 0, 0};

static const int DatInteriorEdgeOfFace[][REFTETRABIS32MAXN_iEpF] =
  { {8}, {9}, {}, {10, 11} };
static const int DatInteriorEdgeOfFaceLen[] =
  { 1, 1, 0, 2};

static const int DatOldEdgeNewVertex[][REFTETRABIS32MAXN_nVpoE] =
  { {0, 1}, {1, 2}, {2, 5, 0}, {0, 4, 3}, {1, 3}, {2, 3} };

static const int DatOldEdgeNewVertexLen[] =
  {2, 2, 3, 3, 2, 2};

static const int DatOldEdgeNewEdge[][REFTETRABIS32MAXN_nEpoE] =
  { {0}, {1}, {2, 3}, {4, 5}, {6}, {7} };

static const int DatOldEdgeNewEdgeLen[] =
  {1, 1, 2, 2, 1, 1};

static const int DatOldFaceNewEdge[][REFTETRABIS32MAXN_nEpoF] =
  { {0, 1, 2, 3}, {4, 5, 6, 0}, {1, 6, 7}, {3, 2, 7, 5, 4} };

static const int DatOldFaceNewEdgeLen[] =
  {4, 4, 3, 5};

static const int DatOldFaceNewFace[][REFTETRABIS32MAXN_nFpoF] =
  { {1, 0}, {2, 3}, {4}, {7, 6, 5} };

static const int DatOldFaceNewFaceLen[] =
  {2, 2, 1, 3};

static const int DatNewEdgeOldEdge[] =
  {0, 1, 2, 2, 3, 3, 4, 5, -1, -1, -1, -1};

static const int DatNewFaceOldFace[] =
  {0, 0, 1, 1, 2, 3, 3, 3, -1, -1};

static const int DatOldFaceNewLocFace[][4] =
  { {0, 1, -1, 3}, {0, -1, -1, 3}, {-1, 1, 2, 3} };

static const int DatChildTwistIndex[] =
  {0, 2, 0, 1, 0, 2, 1, 0, -1, -1};

  Type = TetraBis32;

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
  MaxN_VpC = REFTETRABIS32MAXN_VpC;
  MaxN_CpV = REFTETRABIS32MAXN_CpV;
  MaxN_EpC = REFTETRABIS32MAXN_EpC;
  MaxN_CpE = REFTETRABIS32MAXN_CpE;
  MaxN_EpV = REFTETRABIS32MAXN_EpV;
  MaxN_EpF = REFTETRABIS32MAXN_EpF;
  MaxN_FpE = REFTETRABIS32MAXN_FpE;
  MaxN_VpF = REFTETRABIS32MAXN_VpF;
  MaxN_FpV = REFTETRABIS32MAXN_FpV;
  MaxN_FpC = REFTETRABIS32MAXN_FpC;
  MaxN_CpF = REFTETRABIS32MAXN_CpF;
  MaxN_iVpE = REFTETRABIS32MAXN_iVpE;
  MaxN_iEpF = REFTETRABIS32MAXN_iEpF;
  MaxN_nVpoE = REFTETRABIS32MAXN_nVpoE;
  MaxN_nEpoE = REFTETRABIS32MAXN_nEpoE;
  MaxN_nVpoF = REFTETRABIS32MAXN_nVpoF;
  MaxN_oVpoF = REFTETRABIS32MAXN_oVpoF;
  MaxN_nEpoF = REFTETRABIS32MAXN_nEpoF;
  MaxN_nFpoF = REFTETRABIS32MAXN_nFpoF;

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




