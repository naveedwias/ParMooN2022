#include <RefTetraBis4Desc.h>

// Constructor
TRefTetraBis4Desc::TRefTetraBis4Desc(const TShapeDesc *shape) : TRefDesc(shape)
{
static const Shapes DatChildType[] = {Tetrahedron, Tetrahedron};

static const Refinements DatEdgeType[] = {NoRef, NoRef, NoRef, NoRef, LineReg, NoRef};

static const Refinements DatFaceType [] = {NoRef, TriBis1, TriBis1, NoRef};

/*
 * Vertex
 */
static const int DatChildVertex[][4] = {{0,1,2,4},{0,4,2,3}};
static const int DatVertexChild[][2] = {{0,1},{0},{0,1},{1},{0,1}};
static const int DatVertexChildIndex[][2] = {{0,0},{1},{2,2},{3},{3,1}};
static const int DatVertexChildLen[] = {2,1,2,1,2};

/*
 * Faces
 */
static const int DatChildFace[][4] = {{0,1,3,6},{6,2,4,5}};
static const int DatFaceChild[][2] = {{0},{0},{1},{0},{1},{1},{0,1}};
static const int DatFaceChildIndex[][2] = {{0},{1},{1},{2},{2},{3},{3,0}};
static const int DatFaceChildLen[] = {1,1,1,1,1,1,2};

/*
 * Edges
 */
static const int DatChildEdge[][6] = {{0,1,2,7,4,8},{7,8,2,3,5,6}};
static const int DatEdgeChild[][2] = {{0},{0},{0,1},{1},{0},{1},{1},{0,1},{0,1}};
static const int DatEdgeChildIndex[][2] = {{0},{1},{2,2},{3},{4},{4},{5},{3,0},{5,1}};
static const int DatEdgeChildLen[] = {1,1,2,1,1,1,1,2,2};

/*
 * Edge-Vertex
 */
static const int DatEdgeVertex[][2] = {{0,1},{1,2},{2,0},{0,3},{1,4},{4,3},{2,3},{0,4},{2,4}};
static const int DatVertexEdge[][4] = {{0,2,3,7},{0,1,4},{1,2,6,8},{3,5,6},{4,5,7,8}};
static const int DatVertexEdgeIndex[][4] = {{0,1,0,0},{1,0,0},{1,0,0,0},{1,1,1},{1,0,1,1}};
static const int DatVertexEdgeLen[] = {4,3,4,3,4};

/*
 * Face-Vertex
 */
static const int DatFaceVertex[][3] = {{0,1,2},{0,4,1},{0,3,4},{2,1,4},{2,4,3},{0,2,3},{0,2,4}};
static const int DatVertexFace[][5] = {{0,1,2,5,6},{0,1,3},{0,3,4,5,6},{2,4,5},{1,2,3,4,6}};
static const int DatVertexFaceIndex[][5] = {{0,0,0,0,0},{1,2,1},{2,0,0,1,1},{1,2,2},{1,2,2,1,2}};
static const int DatVertexFaceLen[] = {5,3,5,3,5};

/*
 * Face-Edge
 */
static const int DatFaceEdge[][3] = {{0,1,2},{7,4,0},{3,5,7},{1,4,8},{8,5,6},{2,6,3},{2,8,7}};
static const int DatEdgeFace[][3] = {{0,1},{0,3},{0,5,6},{2,5},{1,3},{2,4},{4,5},{1,2,6},{3,4,6}};
static const int DatEdgeFaceIndex[][3] = {{0,2},{1,0},{2,0,0},{0,2},{1,1},{1,1},{2,1},{0,2,2},{2,0,1}};
static const int DatEdgeFaceLen[] = {2,2,3,2,2,2,2,3,3};

/*
 * New index equal to old index
 */
static const int DatNewVertexEqOldVertex[] = {0, 1, 2, 3};
static const int DatNewVertexEqOldVertexIndex[] = {0, 1, 2, 3};

static const int DatNewFaceEqOldFace[] = {0, 5};
static const int DatNewFaceEqOldFaceIndex[] = {0, 3};

static const int DatOldFaceNewVertex[][REFTETRABIS4MAXN_nVpoF] =
  { {0, 1, 2}, {0, 3, 1, 4}, {2, 1, 3, 4}, {0, 2, 3} };
static const int DatOldFaceNewVertexLen[] =
  { 3, 4, 4, 3};
static const double DatOldFaceNewVertexPos[][REFTETRABIS4MAXN_nVpoF][REFTETRABIS4MAXN_oVpoF] =
  { { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0.5, 0.5} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0.5, 0.5} },
    { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} } };

static const int DatInteriorFaceOfCell[] = {6};
static const int DatInteriorEdgeOfCell[] = {-1};

static const int DatInteriorVertexOfEdge[][REFTETRABIS4MAXN_iVpE] =
  { {}, {}, {}, {}, {4}, {} };

static const int DatInteriorVertexOfEdgeLen[] = {0, 0, 0, 0, 1, 0};

static const int DatInteriorEdgeOfFace[][REFTETRABIS4MAXN_iEpF] =
  { {}, {7}, {8}, {} };

static const int DatInteriorEdgeOfFaceLen[] =
  { 0, 1, 1, 0};

/*
 * Old-New Relations
 */
static const int DatOldEdgeNewVertex[][REFTETRABIS4MAXN_nVpoE] =
  { {0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 4, 3}, {2, 3} };

static const int DatOldEdgeNewVertexLen[] =
  {2, 2, 2, 2, 3, 2};

static const int DatOldEdgeNewEdge[][REFTETRABIS4MAXN_nEpoE] =
  { {0}, {1}, {2}, {3}, {4, 5}, {6} };

static const int DatOldEdgeNewEdgeLen[] =
  {1, 1, 1, 1, 2, 1};

static const int DatOldFaceNewEdge[][REFTETRABIS4MAXN_nEpoF] =
  { {0, 1, 2}, {3, 5, 4, 0}, {1, 4, 5, 6}, {2, 6, 3} };

static const int DatOldFaceNewEdgeLen[] =
  {3, 4, 4, 3};

static const int DatOldFaceNewFace[][REFTETRABIS4MAXN_nFpoF] =
  { {0}, {2, 1}, {3, 4}, {5} };

static const int DatOldFaceNewFaceLen[] =
  {1, 2, 2, 1};

/*
 * New-Old Relations
 */
static const int DatNewEdgeOldEdge[] =
  {0, 1, 2, 3, 4, 4, 5, -1, -1};

static const int DatNewFaceOldFace[] =
  {0, 1, 1, 2, 2, 3, -1};

static const int DatOldFaceNewLocFace[][4] =
  { {0, 1, 2, -1}, {-1, 1, 2, 3} };

static const int DatChildTwistIndex[] =
  {0, 2, 1, 1, 2, 0, -1};

  Type = TetraBis4;

  //set all numbers
  N_Vertices = 5;
  N_Edges = 9;
  N_Faces = 7;
  N_Children = 2;
  N_NewVertEqOldVert = 4;
  N_NewFaceEqOldFace = 2;
  N_InnerEdges = 0;
  N_InnerFaces = 1;

  // initialize all dimension values
  MaxN_VpC = REFTETRABIS4MAXN_VpC;
  MaxN_CpV = REFTETRABIS4MAXN_CpV;
  MaxN_EpC = REFTETRABIS4MAXN_EpC;
  MaxN_CpE = REFTETRABIS4MAXN_CpE;
  MaxN_EpV = REFTETRABIS4MAXN_EpV;
  MaxN_EpF = REFTETRABIS4MAXN_EpF;
  MaxN_FpE = REFTETRABIS4MAXN_FpE;
  MaxN_VpF = REFTETRABIS4MAXN_VpF;
  MaxN_FpV = REFTETRABIS4MAXN_FpV;
  MaxN_FpC = REFTETRABIS4MAXN_FpC;
  MaxN_CpF = REFTETRABIS4MAXN_CpF;
  MaxN_iVpE = REFTETRABIS4MAXN_iVpE;
  MaxN_iEpF = REFTETRABIS4MAXN_iEpF;
  MaxN_nVpoE = REFTETRABIS4MAXN_nVpoE;
  MaxN_nEpoE = REFTETRABIS4MAXN_nEpoE;
  MaxN_nVpoF = REFTETRABIS4MAXN_nVpoF;
  MaxN_oVpoF = REFTETRABIS4MAXN_oVpoF;
  MaxN_nEpoF = REFTETRABIS4MAXN_nEpoF;
  MaxN_nFpoF = REFTETRABIS4MAXN_nFpoF;

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

