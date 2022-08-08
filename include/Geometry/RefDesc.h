#ifndef __REFDESC__
#define __REFDESC__

#include "Enumerations_geometry.h"
#include <ShapeDesc.h>

#ifdef __2D__
  #define MAXN_ORIGEDGES    4
  #define MAXN_NEWVERTICES  9
  #define MAXN_NEWJOINTS   12
  #define MAXN_CHILDREN     4
#else
  #define MAXN_ORIGEDGES   12
  #define MAXN_NEWVERTICES 27
  #define MAXN_NEWJOINTS   36
  #define MAXN_NEWEDGES    54
  #define MAXN_CHILDREN     8
#endif

/** super class of all refinement descriptors */
class TRefDesc
{
  protected:
    /** pointer to corresponding shape descriptor */
    const TShapeDesc *Shape;

    /** type of refinement */
    Refinements Type;

    /** number of vertices */
    int N_Vertices;
    /** number of edges */
    int N_Edges;
    /** number of children */
    int N_Children;

    /** number of vertices of base cell */
    int N_OrigVertices;
    /** number of edges of base cell */
    int N_OrigEdges;

    #ifdef __3D__
      /** number of faces */
      int N_Faces;
      /** number of faces of base cell */
      int N_OrigFaces;
      /** number of inner faces */
      int N_InnerFaces;
    #endif

    /** number of new vertices, which equal to old vertices */
    int N_NewVertEqOldVert;
    /** number of new inner vertices */
    int N_InnerVertices;
    /** number of new edges, which equal to old edges */
    int N_NewEdgeEqOldEdge;
    /** number of new inner edges */
    int N_InnerEdges;

    /** maximum number of vertices per cell */
    int MaxN_VpC;
    /** maximum number of cells per vertex */
    int MaxN_CpV;
    /** maximum number of edges per cell */
    int MaxN_EpC;
    /** maximum number of cells per edge */
    int MaxN_CpE;
    /** maximum number of edges per vertex */
    int MaxN_EpV;
    /** maximum number of inner vertices per edge */
    int MaxN_iVpE;
    /** maximum number of new vertices per old edge */
    int MaxN_nVpoE;
    /** maximum number of new edges per old edge */
    int MaxN_nEpoE;

    /** type of children after refinement */
    const Shapes *ChildType;
    /** refinement's types of the edges */
    const Refinements *EdgeType;

    /** which vertices build a child */
    const int *ChildVertex;
    /** number of vertices per child */
    const int *ChildVertexLen;
    /** which children meet on a vertex */
    const int *VertexChild;
    /** which local index has the vertex in each child */
    const int *VertexChildIndex;
    /** length of array VertexChild */
    const int *VertexChildLen;

    /** which edges build a child */
    const int *ChildEdge;
    /** which children are connectred with an edge */
    const int *EdgeChild;
    /** which local index has the edge in each child */
    const int *EdgeChildIndex;
    /** length of array EdgeChild */
    const int *EdgeChildLen;

    /** which vertices belong to an edge */
    const int *EdgeVertex;
    /** which edges meet on a vertex */
    const int *VertexEdge;
    /** which index has the vertex in each edge */
    const int *VertexEdgeIndex;
    /** length of array VertexEdge */
    const int *VertexEdgeLen;

    /** which new vertices are equal to old vertices */
    const int *NewVertexEqOldVertex;
    /** which old indices have the new vertices */
    const int *NewVertexEqOldVertexIndex;

    /** which new edges are equal to old edges */
    const int *NewEdgeEqOldEdge;
    /** which old indices have the new edges */
    const int *NewEdgeEqOldEdgeIndex;

    /** indices of new inner vertices */
    const int *InteriorVertexOfCell;
    /** position of new inner vertices as the coefficients for a linear
        combination of the old vertices */
    const double *PositionOfIntVert;

    /** indices  of new inner edges */
    const int *InteriorEdgeOfCell;
    /** indices of inner vertices on an edge */
    const int *InteriorVertexOfEdge;
    /** length of array InteriorVertexOfEdge */
    const int *InteriorVertexOfEdgeLen;

    /** which new vertices are on an old edge */
    const int *OldEdgeNewVertex;
    /** length of array NewVertex per old edge */
    const int *OldEdgeNewVertexLen;

    /** which local edge of a child lies on an old edge (for each child) */
    const int *OldEdgeNewLocEdge;

    /** which new edges lie on which old edge */
    const int *OldEdgeNewEdge;
    /** length of array NewEdge per old edge */
    const int *OldEdgeNewEdgeLen;
    /** to which old edge belongs a new edge */
    const int *NewEdgeOldEdge;

    #ifdef __3D__
      /** maximum number of vertices per face */
      int MaxN_VpF;
      /** maximum number of old vertices per old face */
      int MaxN_oVpoF;
      /** maximum number of faces per vertex */
      int MaxN_FpV;
      /** maximum number of edges per face */
      int MaxN_EpF;
      /** maximum number of faces per edge */
      int MaxN_FpE;
      /** maximum number of cells per face */
      int MaxN_CpF;
      /** maximum number of faces per cell */
      int MaxN_FpC;
      /** maximum number of inner vertices per face */
      int MaxN_iVpF;
      /** maximum number of inner edges per face */
      int MaxN_iEpF;
      /** maximum number of new edges per old face */
      int MaxN_nEpoF;
      /** maximum number of new vertices per old face */
      int MaxN_nVpoF;
      /** maximum number of new inner vertices per old face */
      int MaxN_niVpoF;
      /** maximum number of new faces per old face */
      int MaxN_nFpoF;
      /** refinement's types of the faces */
      const Refinements *FaceType;
      /** which children meet on a face */
      const int *FaceChild;
      /** which local index has the face in each child */
      const int *FaceChildIndex;
      /** length of array FaceChild */
      const int *FaceChildLen;

      /** which edges built a face */
      const int *FaceEdge;
      /** which faces are connected with edge . */
      const int *EdgeFace;
      const int *EdgeFaceLen;
      /** which local index has the edge in each face */
      const int *EdgeFaceIndex;

      /** which vertices built a face */
      const int *FaceVertex;
      /** which faces are connected with vertex . */
      const int *VertexFace;
      const int *VertexFaceLen;
      /** which local index has the vertex in each face */
      const int *VertexFaceIndex;

      /** field of new inner faces */
      const int *InteriorFaceOfCell;

      /** which vertices lie on interior of old face */
      const int *InteriorVertexOfFace;
      const int *InteriorVertexOfFaceLen;

      /** which edges lie on interior of old face */
      const int *InteriorEdgeOfFace;
      const int *InteriorEdgeOfFaceLen;

      /** which new Inner vertices lie on which old face */
      const int *OldFaceNewInnerVertices;
      /** lenght of OldFaceNewVertices entries */
      const int *OldFaceNewInnerVerticesLen;

      /** number of new faces equal old faces */
      int N_NewFaceEqOldFace;
      /** which new faces are equal to old faces */
      const int *NewFaceEqOldFace;
      /** which old indices have the new faces */
      const int *NewFaceEqOldFaceIndex;

      /** number of new vertices on an old face */
      const int *N_NewVertsOnOldFace;
      /** new vertices on an old face */
      const int *NewVertsOnOldFace;
      /** position of new vertices on an old face */
      const double *NewVertsOnOldFacePos;

      /** which faces built a child */
      const int *ChildFace;

      /** which new vertices belong to an old face */
      const int *OldFaceNewVertex;
      /** position of new vertices on old face (convex linear combination) */
      const double *OldFaceNewVertexPos;
      /** number of new vertices which belong to an old face */
      const int *OldFaceNewVertexLen;

      /** which new edges belong to an old face */
      const int *OldFaceNewEdge;
      /** number of new vertices which belong to an old face */
      const int *OldFaceNewEdgeLen;

      /** which new faces belong to an old face */
      const int *OldFaceNewFace;
      /** number of new faces which belong to an old face */
      const int *OldFaceNewFaceLen;

      /** on which old face does a new face lie */
      const int *NewFaceOldFace;

      /** ??? */
      const int *OldFaceNewLocFace;

      /** ??? */
      const int *ChildTwistIndex;
    #endif

  public:
    // Constructor
    /** initialize the whole data structure */
    explicit TRefDesc(const TShapeDesc *shape);

    virtual ~TRefDesc(){};
    
    // Methods
    /** return type of refinement */
    Refinements GetType() const
    { return Type; }
    /** return number of children */
    int GetN_Children() const
    { return N_Children; }
    /** return number of edges */
    int GetN_Edges() const
    { return N_Edges; }
    /** return number of vertices */
    int GetN_Vertices() const
    { return N_Vertices; }

    /** return number of edges on base cell */
    int GetN_OrigEdges() const
    { return N_OrigEdges; }
    /** return number of vertices on base cell */
    int GetN_OrigVertices() const
    { return N_OrigVertices; }

    #ifdef __3D__
      /** return number of faces */
      int GetN_OrigFaces() const
      { return N_OrigFaces; }
      /** return number of faces on base cell*/
      int GetN_Faces() const
      { return N_Faces; }
    #endif

    /** return a bool, whether to refine or not */
    virtual int IsToRefine() const
    { return true; }
    
    /** return shape descriptor */
    const TShapeDesc *GetShapeDesc() const
    { return Shape; }
    /** return refinement type of edge J\_i */
    Refinements GetEdgeRef(int J_i) const
    { return EdgeType[J_i]; }
    /** return type of child number C\_i */
    Shapes GetChildType(int C_i) const
    { return ChildType[C_i]; }

    /** return number of new vertices, which equal old vertices */
    int GetN_NewVertEqOldVert() const
    { return N_NewVertEqOldVert; }
    /** return number of new inner vertices */
    int GetN_InnerVertices() const
    { return N_InnerVertices; }
    /** return number of new edges, which equal old edges */
    int GetN_NewEdgeEqOldEdge() const
    { return N_NewEdgeEqOldEdge; }
    /** return number of new inner edges */
    int GetN_InnerEdges() const
    { return N_InnerEdges; }

    /** return auxilary fields in order to copy existing vertices */
    int GetNewVertEqOldVert(const int *&TmpValues, const int *&TmpIndex) const
    {
        TmpValues = NewVertexEqOldVertex;
        TmpIndex = NewVertexEqOldVertexIndex;
        return 0;
     }

    /** return auxilary fields in order to create new inner vertices */
    int GetInnerVerts(const int *&TmpValues, const double *&TmpPos,
          int &MaxLen) const
    {
      TmpValues = InteriorVertexOfCell;
      TmpPos = PositionOfIntVert;
      MaxLen = N_OrigVertices;
      return 0;
    }

    /** return auxilary fields in order to copy existing edges */
    int GetNewEdgeEqOldEdge(const int *&TmpValues, const int *&TmpIndex) const
    {
      TmpValues = NewEdgeEqOldEdge;
      TmpIndex = NewEdgeEqOldEdgeIndex;
      return 0;
    }

    /** return auxilary fields in order to create new inner edges */
    int GetInnerEdges(const int *&TmpinE, const int *&TmpEC, int &MaxLen) const
    {
       TmpinE = InteriorEdgeOfCell;
       TmpEC = EdgeChild;
       MaxLen = MaxN_CpE;
       return 0;
    }

    /** return the array OldEdgeNewEdge */
    int GetOldEdgeNewEdge(const int *&TmpoEnE, const int *&TmpLen, int &MaxLen)
      const
    {
      TmpoEnE = OldEdgeNewEdge;
      TmpLen = InteriorVertexOfEdgeLen;
      MaxLen = MaxN_nEpoE;
      return 0;
    }

    /** return the array OldEdgeNewLocEdge */
    int GetOldEdgeNewLocEdge(const int *&TmpoEnlE) const
    {
      TmpoEnlE=OldEdgeNewLocEdge;
      return 0;
    }

    /** return the array NewEdgeOldEdge */
    int GetNewEdgeOldEdge(const int *&TmpnEoE) const
    {
      TmpnEoE = NewEdgeOldEdge;
      return 0;
    }

    /** return the array EdgeChild */
    int GetEdgeChild(const int *&TmpEC, const int *&TmpLen, int &MaxLen) const
    {
        TmpEC = EdgeChild;
        TmpLen = EdgeChildLen;
        MaxLen = MaxN_CpE;
        return 0;
    }

    /** return the array EdgeChildIndex */
    int GetEdgeChildIndex(const int *&TmpECI, const int *&TmpLen, int &MaxLen) 
      const
    {
      TmpECI = EdgeChildIndex;
      TmpLen = EdgeChildLen;
      MaxLen = MaxN_CpE;
      return 0;
    }

    /** return the array OldEdgeNewVertex */
    int GetOldEdgeNewVertex(const int *&TmpoEnV, const int *&TmpLen, 
                            int &MaxLen) const
    {
      TmpoEnV = OldEdgeNewVertex;
      TmpLen = InteriorVertexOfEdgeLen;
      MaxLen = MaxN_nVpoE;
      return 0;
    }
    
    /** return the array EdgeVertex */
    int GetEdgeVertex(const int *&TmpEV) const
    {
      TmpEV = EdgeVertex;
      return 0;
    }

    /** return the array VertexEdge */
    int GetVertexEdge(const int *&TmpVE, const int *&TmpLen, int &MaxLen) const
    {
      TmpVE = VertexEdge;
      TmpLen = VertexEdgeLen;
      MaxLen = MaxN_EpV;
      return 0;
    }

    /** return the array VertexEdgeIndex */
    int GetVertexEdgeIndex(const int *&TmpVEI, const int *&TmpLen, int &MaxLen)
      const
    {
      TmpVEI = VertexEdgeIndex;
      TmpLen = VertexEdgeLen;
      MaxLen = MaxN_EpV;
      return 0;
    }

    /** return the array VertexChild */
    int GetVertexChild(const int *&TmpVC, const int *&TmpLen, int &MaxLen) const
    {
      TmpVC = VertexChild;
      TmpLen = VertexChildLen;
      MaxLen = MaxN_CpV;
      return 0;
    }

    /** return the array VertexChildIndex */
    int GetVertexChildIndex(const int *&TmpVCI, const int *&TmpLen, int &MaxLen)
      const
    {
      TmpVCI = VertexChildIndex;
      TmpLen = VertexChildLen;
      MaxLen = MaxN_CpV;
      return 0;
    }

    /** return the array ChildVertex */
    int GetChildVertex(const int *&TmpCV, int &MaxLen) const
    {
      TmpCV =  ChildVertex;
      MaxLen = MaxN_VpC;
      return 0;
    }

    /** return the array ChildEdge */
    int GetChildEdge(const int *&TmpCE, int &MaxLen) const
    {
      TmpCE =  ChildEdge;
      MaxLen = MaxN_EpC;
      return 0;
    }

    #ifdef __3D__
      /** return number of inner faces */
      int GetN_InnerFaces() const
      { return N_InnerFaces; }
    
      /** return auxilary fields in order to create new inner faces */
      int GetInnerFaces(const int *&TmpinF, const int *&TmpFC, int &MaxLen)
        const
      {
         TmpinF = InteriorFaceOfCell;
         TmpFC = FaceChild;
         MaxLen = MaxN_CpF;
         return 0;
      }

      /** return field of new vertices on an old face */
      int GetOldFaceNewInnerVertex(const int *&TmpoFniV,
                                   const int *&TmpLen, int &MaxLen) const
      {
        TmpoFniV = OldFaceNewInnerVertices;
        TmpLen = OldFaceNewInnerVerticesLen;
        MaxLen = MaxN_niVpoF;
        return 0;
      }

      /** return field of new faces on old faces */
      int GetOldFaceNewFace(const int *&TmpoFnF, const int *&TmpLen,
                            int &MaxLen) const
      {
        TmpoFnF = OldFaceNewFace;
        TmpLen = OldFaceNewFaceLen;
        MaxLen = MaxN_nFpoF;
        return 0;
      }

      /** return the refinement type of face i */
      Refinements GetFaceRef(int i) const
      { return FaceType[i]; }

      /** return number of new faces equal old faces */
      int GetN_NewFaceEqOldFace() const
      { return N_NewFaceEqOldFace; }

      /** return the array NewFaceEqOldFace */
      int GetNewFaceEqOldFace(const int *&TmpValues, const int *&TmpIndex) const
      {
        TmpValues = NewFaceEqOldFace;
        TmpIndex = NewFaceEqOldFaceIndex;
        return 0;
      }

      /** return the array FaceChild */
      int GetFaceChild(const int *&TmpFC, const int *&TmpLen, int &MaxLen) const
      {
          TmpFC = FaceChild;
          TmpLen = FaceChildLen;
          MaxLen = MaxN_CpF;
          return 0;
      }

      /** return the array FaceChildIndex */
      int GetFaceChildIndex(const int *&TmpFCI, const int *&TmpLen, int &MaxLen)
        const
      {
          TmpFCI = FaceChildIndex;
          TmpLen = FaceChildLen;
          MaxLen = MaxN_CpF;
          return 0;
      }

      /** return the array FaceEdge */
      int GetFaceEdge(const int *&TmpFE, int &MaxLen) const
      {
        TmpFE =  FaceEdge;
        MaxLen = MaxN_EpF;
        return 0;
      }

      /** return field NewVertsOnOldFace for face i */
      int GetNewVertsOnOldFace(const int *&TmpNV, const double *&TmpPos,
                               int &MaxLen) const
      {
        TmpNV = NewVertsOnOldFace;
        TmpPos = NewVertsOnOldFacePos;
        MaxLen = MaxN_nVpoF;

        return 0;
      }

      /** return the field ChildFace */
      int GetChildFace(const int *&TmpCF, int &MaxLen) const
      {
        TmpCF = ChildFace;
        MaxLen = MaxN_FpC;

        return 0;
      }

      /** return the field OldFaceNewVertex */
      int GetOldFaceNewVertex(const int *&TmpoFnV, const int *&TmpLen,
                              int &MaxLen) const
      {
        TmpoFnV = OldFaceNewVertex;
        TmpLen = OldFaceNewVertexLen;
        MaxLen = MaxN_nVpoF;

        return 0;
      }

      /** return the field OldFaceNewVertex */
      int GetOldFaceNewVertex(const int *&TmpoFnV, const double *&TmpPos,
                              const int *&TmpLen, int &MaxLen1, 
                              int &MaxLen2) const
      {
        TmpoFnV = OldFaceNewVertex;
        TmpPos = OldFaceNewVertexPos;
        TmpLen = OldFaceNewVertexLen;
        MaxLen1 = MaxN_nVpoF;
        MaxLen2 = MaxN_oVpoF;

        return 0;
      }

      /** return the field NewFaceOldFace */
      int GetNewFaceOldFace(const int *&TmpnFoF) const
      {
        TmpnFoF = NewFaceOldFace;

        return 0;
      }

      /** return the array OldFaceNewLocFace */
      int GetOldFaceNewLocFace(const int *&TmpoFnlF) const
      {
        TmpoFnlF = OldFaceNewLocFace;

        return 0;
      }

      /** return the array ChildTwistIndex */
      int GetChildTwistIndex(const int *&TmpCTI) const
      {
        TmpCTI = ChildTwistIndex;

        return 0;
      }
      
      // added 25.04.2010 for fixing refinement problem
      /** return the array FaceVertex **/
      int GetFaceVertex (const int *&TmpFV, int &TmpLen) const
      {
        TmpFV = FaceVertex;
        TmpLen = MaxN_VpF;
        
        return 0;
      }

    #endif
};

#endif
