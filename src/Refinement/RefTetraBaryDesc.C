#include <RefTetraBaryDesc.h>

// Constructor
TRefTetraBaryDesc::TRefTetraBaryDesc(const TShapeDesc *shape) : TRefDesc(shape)
{
static const Shapes DatChildType[] = {Tetrahedron, Tetrahedron, Tetrahedron,
                                      Tetrahedron};

static const Refinements DatEdgeType[] = 
  {NoRef, NoRef, NoRef, NoRef, NoRef, NoRef};

static const Refinements DatFaceType[] = { NoRef, NoRef, NoRef, NoRef};

static const int DatChildVertex[][4] =
           {{0,1,2,4}, {0,3,1,4}, {1,3,2,4}, {0,2,3,4}};

static const int DatVertexChild[][4] =
           {{0,1,3}, {0,1,2}, {0,2,3}, {1,2,3}, {0,1,2,3}};

static const int DatVertexChildIndex[][4] =
           {{0,0,0}, {1,2,0}, {2,2,1}, {1,1,2}, {3,3,3,3}};

static const int DatVertexChildLen[] = {3, 3, 3, 4};

static const int DatChildEdge[][6] =
           {{0,4,1,3,6,8}, {2,5,0,3,9,6}, {5,7,4,6,9,8}, {1,7,2,3,8,9}};

static const int DatEdgeChild[][3] =
           {{0,1}, {0,3}, {1,3}, {0,1,3}, {0,2}, {1,2}, {0,1,2}, {2,3},
            {0,2,3}, {1,2,3}};

static const int DatEdgeChildIndex[][3] =
           {{0,2}, {2,0}, {0,2}, {3,3,3}, {1,2}, {1,0}, {4,5,3},
            {1,1}, {5,5,4}, {4,4,5}};

static const int DatEdgeChildLen[] = {2, 2, 2, 3, 2, 2, 3, 2, 3, 3};

static const int DatChildFace[][4] =
           {{0,4,7,5}, {1,6,8,4}, {2,8,9,7}, {3,5,9,6}};

static const int DatFaceChild[][2] =
           {{0}, {1}, {2}, {3}, {0,1}, {0,3}, {1,3}, {0,2}, {1,2}, {2,3}};

static const int DatFaceChildIndex[][2] =
           {{0}, {0}, {0}, {0}, {1,3}, {3,1}, {1,3}, {2,3}, {2,1}, {2,2}};

static const int DatFaceChildLen[] =
           {1, 1, 1, 1, 2, 2, 2, 2, 2, 2};

static const int DatEdgeVertex[][2] =
           {{0,1}, {0,2}, {0,3}, {0,4}, {1,2}, {1,3}, {1,4}, {2,3},
            {2,4}, {3,4}};

static const int DatVertexEdge[][4] =
           {{0,1,2,3}, {0,4,5,6}, {1,4,7,8}, {2,5,7,9}, {3,6,8,9}};

static const int DatVertexEdgeIndex[][4] =
           {{0,0,0,0}, {1,0,0,0}, {1,1,0,0}, {1,1,1,0}, {1,1,1,1}};

static const int DatVertexEdgeLen[] = {4, 4, 4, 4, 4};

static const int DatFaceVertex[][3] =
           {{0,1,2}, {0,1,3}, {1,2,3}, {2,0,3},
            {0,1,4}, {0,2,4}, {0,3,4}, {1,2,4}, {1,3,4}, {2,3,4}};

static const int DatVertexFace[][6] =
          {{0,1,3,4,5,6}, {0,1,2,4,7,8}, {0,2,3,5,7,9},
           {1,2,3,6,8,9}, {4,5,6,7,8,9}};

static const int DatVertexFaceIndex[][6] =
           {{0,0,1,0,0,0}, {1,1,0,1,0,0}, {2,1,0,1,1,0}, 
            {2,2,2,1,1,1}, {2,2,2,2,2,2}};

static const int DatVertexFaceLen[] = {6, 6, 6, 6, 6};

static const int DatFaceEdge[][3] =
           {{0,1,4}, {0,2,5}, {4,5,7}, {1,2,7}, 
            {0,3,6}, {1,3,8}, {2,3,9}, {4,6,8}, {5,6,9}, {7,8,9}};

static const int DatEdgeFace[][3] =
           {{0,1,4}, {0,3,5}, {1,3,6}, {4,5,6}, {0,2,7}, {1,2,8}, 
            {4,7,8}, {2,3,9}, {5,7,9}, {6,8,9}};

static const int DatEdgeFaceIndex[][3] =
           {{0,0,0}, {1,0,0}, {1,1,0}, {1,1,1}, {2,0,0}, {2,1,0},
            {2,1,1}, {2,2,0}, {2,2,1}, {2,2,2}};

static const int DatEdgeFaceLen[] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3};

static const int DatNewVertexEqOldVertex[] = {0, 1, 2, 3};
static const int DatNewVertexEqOldVertexIndex[] = {0, 1, 2, 3};

static const int DatNewFaceEqOldFace[] = {0, 1, 2, 3};
static const int DatNewFaceEqOldFaceIndex[] = {0, 1, 2, 3};

static const int DatInteriorVertexOfCell[] = {4};
static const double DatPositionOfIntVert[][4] =
        { { 0.25, 0.25, 0.25, 0.25 } };

static const int DatInteriorEdgeOfCell[] = {6, 7, 8, 9};

static const int DatInteriorFaceOfCell[] = {4, 5, 6, 7, 8, 9};

static const int* DatInteriorVertexOfEdge = nullptr;
static const int DatInteriorVertexOfEdgeLen[] = {0, 0, 0, 0, 0, 0};

static const int* DatInteriorVertexOfFace = nullptr;
static const int DatInteriorVertexOfFaceLen[] = { 0, 0, 0, 0};

static const int* DatInteriorEdgeOfFace = nullptr;
static const int DatInteriorEdgeOfFaceLen[] = {0, 0, 0, 0};

static const int DatOldEdgeNewVertex[][2] =
           {{0,1}, {1,2}, {2,0}, {0,3}, {1,3}, {2,3}};

static const int DatOldEdgeNewVertexLen[] = {2, 2, 2, 2, 2, 2};

static const int DatOldEdgeNewEdge[][1] =
           {{0}, {4}, {1}, {2}, {5}, {7}};

static const int DatOldEdgeNewEdgeLen[] =
           { 1, 1, 1, 1, 1, 1};

static const int DatNewFaceOldFace[] = {0, 1, 2, 3, -1, -1, -1, -1};

static const int DatOldFaceNewVertex[][3] =
           {{0,1,2}, {0,3,1}, {2,1,3}, {0,2,3}};

static const double DatOldFaceNewVertexPos[][3][3] =
           {{{1,0,0}, {0,1,0}, {0,0,1}}, 
            {{1,0,0}, {0,1,0}, {0,0,1}}, 
            {{1,0,0}, {0,1,0}, {0,0,1}}, 
            {{1,0,0}, {0,1,0}, {0,0,1}}};
             
static const int DatOldFaceNewVertexLen[] = {3, 3, 3, 3};

static const int DatOldFaceNewEdge[][3] =
           {{0,4,1}, {2,5,0}, {4,5,7}, {1,7,2}};

static const int DatOldFaceNewEdgeLen[] = {3, 3, 3, 3};

static const int DatOldFaceNewFace[][1] = {{0}, {1}, {2}, {3}};

static const int DatOldFaceNewFaceLen[] = {1, 1, 1, 1};

static const int DatOldFaceNewLocFace[][4] =
           {{0,-1,-1,-1}, {1,-1,-1,-1}, {2,-1,-1,-1}, {3,-1,-1,-1}};

static const int DatChildTwistIndex[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

static const int DatNewEdgeOldEdge[] = {0, 2, 3, -1, 1, 4, -1, 5, -1, -1};


  Type = TetraReg0;

  //set all numbers
  N_Vertices = 5;
  N_Edges = 10;
  N_Faces = 10;
  N_Children = 4;
  N_NewVertEqOldVert = 4;
  N_NewFaceEqOldFace = 4;
  
  N_InnerVertices = 1;
  N_InnerEdges = 4;
  N_InnerFaces = 6;

  // initialize all dimension values
  MaxN_VpC = 4;
  MaxN_CpV = 4;
  MaxN_EpC = 6;
  MaxN_CpE = 3;
  MaxN_EpV = 4;
  MaxN_EpF = 3;
  MaxN_FpE = 3;
  MaxN_VpF = 3;
  MaxN_FpV = 6;
  MaxN_FpC = 4;
  MaxN_CpF = 2;
  MaxN_iVpE = 0;
  MaxN_iEpF = 0;
  MaxN_nVpoE = 2;
  MaxN_nEpoE = 1;
  MaxN_nVpoF = 3;
  MaxN_oVpoF = 3;
  MaxN_nEpoF = 3;
  MaxN_nFpoF = 1;

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
  
  InteriorVertexOfCell = (const int *) DatInteriorVertexOfCell;
  PositionOfIntVert = (const double *) DatPositionOfIntVert;
  
  NewFaceEqOldFace = (const int *) DatNewFaceEqOldFace;
  NewFaceEqOldFaceIndex = (const int *) DatNewFaceEqOldFaceIndex;

  InteriorEdgeOfCell = (const int *) DatInteriorEdgeOfCell;
  InteriorFaceOfCell = (const int *) DatInteriorFaceOfCell;
  InteriorVertexOfEdge = (const int *) DatInteriorVertexOfEdge;
  InteriorVertexOfEdgeLen = (const int *) DatInteriorVertexOfEdgeLen;
  InteriorEdgeOfFace = (const int *) DatInteriorEdgeOfFace;
  InteriorEdgeOfFaceLen = (const int *) DatInteriorEdgeOfFaceLen;
  InteriorVertexOfFace = (const int *) DatInteriorVertexOfFace;
  InteriorVertexOfFaceLen = (const int *) DatInteriorVertexOfFaceLen;

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
