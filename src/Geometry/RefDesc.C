// =======================================================================
// @(#)RefDesc.C        1.3 11/15/99
//
// Class:       TRefDesc
// Purpose:     super class of all refinement descriptors
//
// Author:      Volker Behns  18.06.97
//
// =======================================================================

#include <RefDesc.h>

// Constructor
TRefDesc::TRefDesc(const TShapeDesc *shape)
{
  Shape = shape;

  Type = NoRef;

  N_OrigEdges = Shape->GetN_Edges();
  N_OrigVertices = Shape->GetN_Vertices();

  N_Edges = 0;
  N_Vertices = 0;
  
  N_Children = 0;

  #ifdef __3D__
    N_OrigFaces = Shape->GetN_Faces();
    N_Faces = 0;
    N_InnerFaces = 0;
  #endif

  N_NewVertEqOldVert = 0;
  N_InnerVertices = 0;
  N_NewEdgeEqOldEdge = 0;
  N_InnerEdges = 0;

  MaxN_VpC = 0;
  MaxN_CpV = 0;
  MaxN_EpC = 0;
  MaxN_CpE = 0;
  MaxN_EpV = 0;
  MaxN_iVpE = 0;
  MaxN_nVpoE = 0;
  MaxN_nEpoE = 0;

  ChildType = nullptr;
  EdgeType = nullptr;

  ChildVertex = nullptr;
  ChildVertexLen = nullptr;
  VertexChild = nullptr;
  VertexChildIndex = nullptr;
  VertexChildLen = nullptr;

  ChildEdge = nullptr;
  EdgeChild = nullptr;
  EdgeChildIndex = nullptr;
  EdgeChildLen = nullptr;

  EdgeVertex = nullptr;
  VertexEdge = nullptr;
  VertexEdgeIndex = nullptr;
  VertexEdgeLen = nullptr;

  NewVertexEqOldVertex = nullptr;
  NewVertexEqOldVertexIndex = nullptr;

  InteriorVertexOfCell = nullptr;
  PositionOfIntVert = nullptr;

  NewEdgeEqOldEdge = nullptr;
  NewEdgeEqOldEdgeIndex = nullptr;

  InteriorEdgeOfCell = nullptr;
  InteriorVertexOfEdge = nullptr;
  InteriorVertexOfEdgeLen = nullptr;

  OldEdgeNewVertex = nullptr;
  OldEdgeNewVertexLen = nullptr;
  OldEdgeNewLocEdge = nullptr;
  
  OldEdgeNewEdge = nullptr;
  OldEdgeNewEdgeLen = nullptr;
  NewEdgeOldEdge = nullptr;

  #ifdef __3D__
  MaxN_VpF = 0;
  MaxN_oVpoF = 0;
  MaxN_FpV = 0;
  MaxN_EpF = 0;
  MaxN_FpE = 0;
  MaxN_CpF = 0;
  MaxN_FpC = 0;
  MaxN_iVpF = 0;
  MaxN_iEpF = 0;
  MaxN_nEpoF = 0;
  MaxN_nVpoF = 0;
  MaxN_niVpoF = 0;
  MaxN_nFpoF = 0;

  FaceType = nullptr;

  ChildFace = nullptr;
  FaceChild = nullptr;
  FaceChildIndex = nullptr;
  FaceChildLen = nullptr;

  FaceVertex = nullptr;
  VertexFace = nullptr;
  VertexFaceIndex = nullptr;
  VertexFaceLen = nullptr;

  FaceEdge = nullptr;
  EdgeFace = nullptr;
  EdgeFaceIndex = nullptr;
  EdgeFaceLen = nullptr;

  InteriorFaceOfCell = nullptr;
  InteriorVertexOfFace = nullptr;
  InteriorVertexOfFaceLen = nullptr;
  InteriorEdgeOfFace = nullptr;
  InteriorEdgeOfFaceLen = nullptr;

  OldFaceNewInnerVertices = nullptr;
  OldFaceNewInnerVerticesLen = nullptr;
  
  N_NewFaceEqOldFace = 0;
  NewFaceEqOldFace = nullptr;
  NewFaceEqOldFaceIndex = nullptr;
  
  N_NewVertsOnOldFace = nullptr;
  NewVertsOnOldFace = nullptr;
  NewVertsOnOldFacePos = nullptr;
  
  NewFaceOldFace = nullptr;

  OldFaceNewVertex = nullptr;
  OldFaceNewVertexPos = nullptr;
  OldFaceNewVertexLen = nullptr;
  OldFaceNewEdge = nullptr;
  OldFaceNewEdgeLen = nullptr;
  OldFaceNewFace = nullptr;
  OldFaceNewFaceLen = nullptr;

  OldFaceNewLocFace = nullptr;
  ChildTwistIndex = nullptr;

  #endif
}
 

// Methods
