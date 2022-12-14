#include <RefTetraRegDesc.h>

// Constructor
TRefTetraRegDesc::TRefTetraRegDesc(const TShapeDesc *shape) : TRefDesc(shape)
{
static const Shapes DatChildType[] = {Tetrahedron, Tetrahedron, Tetrahedron,
                                      Tetrahedron, Tetrahedron, Tetrahedron,
                                      Tetrahedron, Tetrahedron };

static const Refinements DatEdgeType[] = { LineReg, LineReg, LineReg, LineReg, 
                                           LineReg, LineReg};

static const Refinements DatFaceType[] = { TriReg, TriReg, TriReg, TriReg};

static const int DatChildVertex[][REFTETRAREGMAXN_VpC] =
           { { 0, 4, 6, 7}, { 4, 1, 5, 8}, { 6, 5, 2, 9}, { 7, 8, 9, 3}, { 4, 6, 7, 8}, { 5, 6, 4, 8}, { 6, 7, 8, 9}, { 8, 5, 6, 9}, };

static const int DatVertexChild[][REFTETRAREGMAXN_CpV] =
           { { 0}, { 1}, { 2}, { 3}, { 0, 1, 4, 5}, { 1, 2, 5, 7}, { 0, 2, 4, 5, 6, 7}, { 0, 3, 4, 6}, { 1, 3, 4, 5, 6, 7}, { 2, 3, 6, 7}, };

static const int DatVertexChildIndex[][REFTETRAREGMAXN_CpV] =
           { { 0}, { 1}, { 2}, { 3}, { 1, 0, 0, 2}, { 2, 1, 0, 1}, { 2, 0, 1, 1, 0, 2}, { 3, 0, 2, 1}, { 3, 1, 3, 3, 2, 0}, { 3, 2, 3, 3}, };

static const int DatVertexChildLen[] =
           { 1, 1, 1, 1, 4, 4, 6, 4, 6, 4};

static const int DatChildEdge[][REFTETRAREGMAXN_EpC] =
           { { 0, 12, 5, 6, 15, 22}, { 1, 2, 13, 16, 8, 18}, { 14, 3, 4, 21, 19, 10}, { 17, 20, 23, 7, 9, 11}, { 12, 22, 15, 16, 24, 17}, { 14, 12, 13, 18, 24, 16}, { 22, 17, 24, 21, 23, 20}, { 18, 14, 24, 20, 19, 21}, };

static const int DatEdgeChild[][REFTETRAREGMAXN_CpE] =
           { { 0}, { 1}, { 1}, { 2}, { 2}, { 0}, { 0}, { 3}, { 1}, { 3}, { 2}, { 3}, { 0, 4, 5}, { 1, 5}, { 2, 5, 7}, { 0, 4}, { 1, 4, 5}, { 3, 4, 6}, { 1, 5, 7}, { 2, 7}, { 3, 6, 7}, { 2, 6, 7}, { 0, 4, 6}, { 3, 6}, { 4, 5, 6, 7}, };

static const int DatEdgeChildIndex[][REFTETRAREGMAXN_CpE] =
           { { 0}, { 0}, { 1}, { 1}, { 2}, { 2}, { 3}, { 3}, { 4}, { 4}, { 5}, { 5}, { 1, 0, 1}, { 2, 2}, { 0, 0, 1}, { 4, 2}, { 3, 3, 5}, { 0, 5, 1}, { 5, 3, 0}, { 4, 4}, { 1, 5, 3}, { 3, 3, 5}, { 5, 1, 0}, { 2, 4}, { 4, 4, 2, 2}, };

static const int DatEdgeChildLen[] =
           { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2, 3, 3, 3, 2, 3, 3, 3, 2, 4};

static const int DatChildFace[][REFTETRAREGMAXN_FpC] =
           { { 0, 1, 2, 3}, { 4, 5, 6, 7}, { 8, 9, 10, 11}, { 12, 13, 14, 15}, { 2, 16, 17, 18}, { 19, 20, 16, 7}, { 17, 21, 12, 22}, { 20, 23, 9, 22}, };

static const int DatFaceChild[][REFTETRAREGMAXN_CpF] =
           { { 0}, { 0}, { 0, 4}, { 0}, { 1}, { 1}, { 1}, { 1, 5}, { 2}, { 2, 7}, { 2}, { 2}, { 3, 6}, { 3}, { 3}, { 3}, { 4, 5}, { 4, 6}, { 4}, { 5}, { 5, 7}, { 6}, { 6, 7}, { 7}, };

static const int DatFaceChildIndex[][REFTETRAREGMAXN_CpF] =
           { { 0}, { 1}, { 2, 0}, { 3}, { 0}, { 1}, { 2}, { 3, 3}, { 0}, { 1, 2}, { 2}, { 3}, { 0, 2}, { 1}, { 2}, { 3}, { 1, 2}, { 2, 0}, { 3}, { 0}, { 1, 0}, { 1}, { 3, 3}, { 1}, };

static const int DatFaceChildLen[] =
           { 1, 1, 2, 1, 1, 1, 1, 2, 1, 2, 1, 1, 2, 1, 1, 1, 2, 2, 1, 1, 2, 1, 2, 1};

static const int DatEdgeVertex[][2] =
           { { 0, 4}, { 4, 1}, { 1, 5}, { 5, 2}, { 2, 6}, { 6, 0}, { 0, 7}, { 7, 3}, { 1, 8}, { 8, 3}, { 2, 9}, { 9, 3}, { 6, 4}, { 4, 5}, { 5, 6}, { 7, 4}, { 4, 8}, { 8, 7}, { 8, 5}, { 5, 9}, { 9, 8}, { 9, 6}, { 6, 7}, { 7, 9}, { 6, 8}, };

static const int DatVertexEdge[][REFTETRAREGMAXN_EpV] =
           { { 0, 5, 6}, { 1, 2, 8}, { 3, 4, 10}, { 7, 9, 11}, { 0, 1, 12, 13, 15, 16}, { 2, 3, 13, 14, 18, 19}, { 4, 5, 12, 14, 21, 22, 24}, { 6, 7, 15, 17, 22, 23}, { 8, 9, 16, 17, 18, 20, 24}, { 10, 11, 19, 20, 21, 23}, };

static const int DatVertexEdgeIndex[][REFTETRAREGMAXN_EpV] =
           { { 0, 1, 0}, { 1, 0, 0}, { 1, 0, 0}, { 1, 1, 1}, { 1, 0, 1, 0, 1, 0}, { 1, 0, 1, 0, 1, 0}, { 1, 0, 0, 1, 1, 0, 0}, { 1, 0, 0, 1, 1, 0}, { 1, 0, 1, 0, 0, 1, 1}, { 1, 0, 1, 0, 0, 1}, };

static const int DatVertexEdgeLen[] =
           { 3, 3, 3, 3, 6, 6, 7, 6, 7, 6};

static const int DatFaceVertex[][REFTETRAREGMAXN_VpF] =
           { { 0, 4, 6}, { 0, 7, 4}, { 6, 4, 7}, { 0, 6, 7}, { 4, 1, 5}, { 4, 8, 1}, { 5, 1, 8}, { 4, 5, 8}, { 6, 5, 2}, { 6, 9, 5}, { 2, 5, 9}, { 6, 2, 9}, { 7, 8, 9}, { 7, 3, 8}, { 9, 8, 3}, { 7, 9, 3}, { 4, 8, 6}, { 7, 6, 8}, { 4, 7, 8}, { 5, 6, 4}, { 5, 8, 6}, { 6, 9, 7}, { 6, 8, 9}, { 8, 9, 5}, };

static const int DatVertexFace[][REFTETRAREGMAXN_FpV] =
           { { 0, 1, 3}, { 4, 5, 6}, { 8, 10, 11}, { 13, 14, 15}, { 0, 1, 2, 4, 5, 7, 16, 18, 19}, { 4, 6, 7, 8, 9, 10, 19, 20, 23}, { 0, 2, 3, 8, 9, 11, 16, 17, 19, 20, 21, 22}, { 1, 2, 3, 12, 13, 15, 17, 18, 21}, { 5, 6, 7, 12, 13, 14, 16, 17, 18, 20, 22, 23}, { 9, 10, 11, 12, 14, 15, 21, 22, 23}, };

static const int DatVertexFaceIndex[][REFTETRAREGMAXN_FpV] =
           { { 0, 0, 2}, { 1, 2, 0}, { 2, 2, 0}, { 1, 1, 1}, { 1, 2, 0, 0, 0, 2, 0, 2, 2}, { 2, 2, 0, 1, 2, 0, 0, 0, 2}, { 2, 2, 0, 0, 0, 2, 2, 0, 1, 2, 0, 2}, { 1, 1, 1, 0, 0, 2, 2, 0, 2}, { 1, 1, 1, 1, 2, 0, 1, 1, 1, 1, 0, 0}, { 1, 1, 1, 2, 2, 0, 1, 1, 1}, };

static const int DatVertexFaceLen[] =
           { 3, 3, 3, 3, 9, 9, 12, 9, 12, 9};

static const int DatFaceEdge[][REFTETRAREGMAXN_EpF] =
           { { 0, 12, 5}, { 6, 15, 0}, { 15, 22, 12}, { 22, 6, 5}, { 1, 2, 13}, { 16, 8, 1}, { 8, 18, 2}, { 18, 16, 13}, { 14, 3, 4}, { 21, 19, 14}, { 19, 10, 3}, { 10, 21, 4}, { 17, 20, 23}, { 7, 9, 17}, { 9, 11, 20}, { 11, 7, 23}, { 16, 24, 12}, { 24, 17, 22}, { 17, 16, 15}, { 14, 12, 13}, { 18, 24, 14}, { 21, 23, 22}, { 20, 21, 24}, { 20, 19, 18}, };

static const int DatEdgeFace[][REFTETRAREGMAXN_FpE] =
           { { 0, 1}, { 4, 5}, { 4, 6}, { 8, 10}, { 8, 11}, { 0, 3}, { 1, 3}, { 13, 15}, { 5, 6}, { 13, 14}, { 10, 11}, { 14, 15}, { 0, 2, 16, 19}, { 4, 7, 19}, { 8, 9, 19, 20}, { 1, 2, 18}, { 5, 7, 16, 18}, { 12, 13, 17, 18}, { 6, 7, 20, 23}, { 9, 10, 23}, { 12, 14, 22, 23}, { 9, 11, 21, 22}, { 2, 3, 17, 21}, { 12, 15, 21}, { 16, 17, 20, 22}, };

static const int DatEdgeFaceIndex[][REFTETRAREGMAXN_FpE] =
           { { 0, 2}, { 0, 2}, { 1, 2}, { 1, 2}, { 2, 2}, { 2, 2}, { 0, 1}, { 0, 1}, { 1, 0}, { 1, 0}, { 1, 0}, { 1, 0}, { 1, 2, 2, 1}, { 2, 2, 2}, { 0, 2, 0, 2}, { 1, 0, 2}, { 0, 1, 0, 1}, { 0, 2, 1, 0}, { 1, 0, 0, 2}, { 1, 0, 1}, { 1, 2, 0, 0}, { 0, 1, 0, 1}, { 1, 0, 2, 2}, { 2, 2, 1}, { 1, 0, 1, 2}, };

static const int DatEdgeFaceLen[] =
           { 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 3, 4, 3, 4, 4, 4, 3, 4, 4, 4, 3, 4};

static const int DatNewVertexEqOldVertex[] = { 0, 1, 2, 3};
static const int DatNewVertexEqOldVertexIndex[] = { 0, 1, 2, 3};

static const double DatOldFaceNewVertexPos[][REFTETRAREGMAXN_nVpoF]
                                            [REFTETRAREGMAXN_oVpoF] =
           { { { 1, 0, 0}, {0, 1, 0}, { 0, 0, 1}, 
               { 0.5, 0.5, 0}, { 0, 0.5, 0.5}, { 0.5, 0, 0.5}},
             { { 1, 0, 0}, {0, 1, 0}, { 0, 0, 1}, 
               { 0.5, 0.5, 0}, { 0, 0.5, 0.5}, { 0.5, 0, 0.5}},
             { { 1, 0, 0}, {0, 1, 0}, { 0, 0, 1}, 
               { 0.5, 0.5, 0}, { 0, 0.5, 0.5}, { 0.5, 0, 0.5}},
             { { 1, 0, 0}, {0, 1, 0}, { 0, 0, 1}, 
               { 0.5, 0.5, 0}, { 0, 0.5, 0.5}, { 0.5, 0, 0.5}}};

static const int DatInteriorEdgeOfCell[] =
           { 24};

static const int DatInteriorFaceOfCell[] =
           { 2, 7, 9, 12, 16, 17, 20, 22};

static const int DatInteriorVertexOfEdge[][REFTETRAREGMAXN_iVpE] =
           {{4}, {5}, {6}, {7}, {8}, {9}};

static const int DatInteriorVertexOfEdgeLen[]=
           { 1, 1, 1, 1, 1, 1};

static const int DatInteriorEdgeOfFace[][REFTETRAREGMAXN_iEpF] =
           { { 12, 13, 14}, { 15, 16, 17}, { 18, 19, 20}, { 21, 22, 23}, };

static const int DatInteriorEdgeOfFaceLen[] =
           { 3, 3, 3, 3};

static const int DatOldEdgeNewVertex[][REFTETRAREGMAXN_nVpoE] =
           { { 0, 4, 1}, { 1, 5, 2}, { 2, 6, 0}, { 0, 7, 3}, { 1, 8, 3}, { 2, 9, 3}, };

static const int DatOldEdgeNewVertexLen[] = 
           { 3, 3, 3, 3, 3, 3 };

static const int DatOldEdgeNewEdge[][REFTETRAREGMAXN_nEpoE] =
           { { 0, 1}, { 2, 3}, { 4, 5}, { 6, 7}, { 8, 9}, { 10, 11}, };

static const int DatOldEdgeNewEdgeLen[] =
           { 2, 2, 2, 2, 2, 2};

static const int DatOldFaceNewVertex[][REFTETRAREGMAXN_nVpoF] =
           { { 0, 1, 2, 4, 5, 6}, { 0, 3, 1, 7, 8, 4}, { 2, 1, 3, 5, 8, 9}, { 0, 2, 3, 6, 9, 7}, };

static const int DatOldFaceNewVertexLen[] =
           { 6, 6, 6, 6};

static const int DatOldFaceNewEdge[][REFTETRAREGMAXN_nEpoF] =
           { { 0, 1, 2, 3, 4, 5, 12, 13, 14}, { 0, 1, 6, 7, 8, 9, 15, 16, 17}, { 2, 3, 8, 9, 10, 11, 18, 19, 20}, { 4, 5, 6, 7, 10, 11, 21, 22, 23}, };

static const int DatOldFaceNewEdgeLen[] =
           { 9, 9, 9, 9};

static const int DatOldFaceNewFace[][REFTETRAREGMAXN_nFpoF] =
           { { 0, 4, 8, 19}, { 1, 13, 5, 18}, { 10, 6, 14, 23}, { 3, 11, 15, 21}, };

static const int DatOldFaceNewFaceLen[] =
           { 4, 4, 4, 4};

static const int DatNewEdgeOldEdge[] = 
           { 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

static const int DatNewFaceOldFace[] =
           { 0, 1, -1, 3, 0, 1, 2, -1, 0, -1, 2, 3, -1, 1, 2, 3, -1, -1, 1, 0, -1, 3, -1, 2};

static const int DatOldFaceNewLocFace[][TETRAN_F] =
           { { 0, 1, -1, 3}, { 0, 1, 2, -1}, { 0, -1, 2, 3}, { -1, 1, 2, 3}, { -1, 3, -1, -1}, { 0, -1, -1, -1}, { -1, -1, -1, 1}, { -1, -1, 1, -1}, };

static const int DatChildTwistIndex[] =
           { 0, 0, -1, 0, 1, 2, 1, -1, 2, -1, 0, 1, -1, 1, 2, 2, -1, -1, 1, 2, -1, 0, -1, 2};


  Type = TetraReg;

  //set all numbers
  N_Vertices = 10;
  N_Edges = 25;
  N_Faces = 24;
  N_Children = 8;
  N_NewVertEqOldVert = 4;
  N_InnerEdges = 1;
  N_InnerFaces = 8;

  // initialize all dimension values
  MaxN_VpC = REFTETRAREGMAXN_VpC;
  MaxN_CpV = REFTETRAREGMAXN_CpV;
  MaxN_EpC = REFTETRAREGMAXN_EpC;
  MaxN_CpE = REFTETRAREGMAXN_CpE;
  MaxN_EpV = REFTETRAREGMAXN_EpV;
  MaxN_EpF = REFTETRAREGMAXN_EpF;
  MaxN_FpE = REFTETRAREGMAXN_FpE;
  MaxN_VpF = REFTETRAREGMAXN_VpF;
  MaxN_FpV = REFTETRAREGMAXN_FpV;
  MaxN_FpC = REFTETRAREGMAXN_FpC;
  MaxN_CpF = REFTETRAREGMAXN_CpF;
  MaxN_iVpE = REFTETRAREGMAXN_iVpE;
  MaxN_iEpF = REFTETRAREGMAXN_iEpF;
  MaxN_nVpoE = REFTETRAREGMAXN_nVpoE;
  MaxN_nEpoE = REFTETRAREGMAXN_nEpoE;
  MaxN_nVpoF = REFTETRAREGMAXN_nVpoF;
  MaxN_oVpoF = REFTETRAREGMAXN_oVpoF;
  MaxN_nEpoF = REFTETRAREGMAXN_nEpoF;
  MaxN_nFpoF = REFTETRAREGMAXN_nFpoF;

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
