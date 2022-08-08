#ifndef __SHAPEDESC__
#define __SHAPEDESC__

class TVertex;

#include "Enumerations_geometry.h"

constexpr int MAXN_JOINTS2D = 4;
constexpr int MAXN_JOINTS3D = 6;
constexpr int MAXN_EDGES3D = 12;

/** This is a super class of all shape descriptors.
 *
 * In ParMooN different cell geometries, e.g., triangles or hexahedra, are
 * encoded using shape descriptors. For each type of cell a different shape
 * descriptor exists describing the actual look and shape of the cell. See \ref
 * Shapes for all possible shapes a cell can have in ParMooN.
 *
 * In general, every cell consists of vertices and edges, and in 3D also of
 * faces, and different shapes differ in the number of vertices they consist of,
 * the number of edges and so on. Therefore, the description that is used in
 * ParMooN is based on exactly these information. How this is done and which
 * information are used is described here in general and in particular for
 * the hexahedron depicted in the following picture.
 *
 * \image html ShapeDesc_hexahedron-1.svg "Example hexahedron with vertex numbers"
 *
 * First of all, the number of vertices \ref N_Vertices, e.g. 3 for triangles,
 * 4 for quadrilaterals or 4 for tetrahedra, is an basic information about the
 * shape of a cell. For the hexahedron N_Vertices equals 8. Another basic
 * information is the number of edges \ref N_Edges that is used to form the
 * cell, e.g. 4 for quadrilaterals or 6 for tetrahedra. The hexahedron has 12
 * edges, hence N_Edges=12. In 3D also the number of
 * faces \ref N_Faces is of interest, e.g. 4 for tetrahedra and 6 for hexahedra.
 * ParMooN stores and needs in addition another redundant information about the
 * cells namely the number of joints \ref N_Joints, i.e. the number of (d-1)
 * dimensional facets, which equals the number of edges in 2D and the number of
 * faces in 3D, respectively. If all these information are known, one is in the
 * position to describe the particular shape using the vertices. The shape of
 * the cell that is described is stored in \ref Type.
 *
 * Assume that there is a fixed local numbering of the vertices \f$V_1, V_2,
 * ..., V_{N_{Vertices}}\f$. Each edge consists of two vertices, i.e. a starting
 * point and end point of the edge. Therefore, edges can be describes as
 * 2-tuples of vertices, e.g. \f$\{0, 1\}\f$ representing the vertices 0 and 1,
 * or \f$\{5, 3\}\f$ representing the vertices 5 and 3, and so on. These
 * information about the edges are stored in an array \ref EdgeVertex points to.
 * For the hexahedron the array is given by
 * \f[
 * EdgeVertex = \{ \{0, 1\},  \{1, 2\},  \{2, 3\},  \{3, 0\}, \{0, 4\},  \{1,
 * 5\},  \{2, 6\},  \{3, 7\}, \{4, 5\},  \{5, 6\},  \{6, 7\},  \{7, 4\}\}.
 * \f]
 * As you can see, there are N_Edges tuple of Vertices, each defining
 * a particular edge of the hexahedron. This list of edges implies also
 * a local numbering of the edges. Locally the first edge is given by the first
 * tuple in EdgeVertex, the second edge is given by the second tuple and the
 * N_Edges-th edge is given by the last tuple. For the hexahedron given above
 * the edge numbering is indicated with blue numbers in the following picture.
 *
 * \image html ShapeDesc_hexahedron-2.svg "Example hexahedron with vertex numbers (black) and local edge numbers (blue)"
 *
 * After the edges are defined also the faces can be defined using the vertices.
 * The information about which vertices build a face are stored in an array \ref
 * FaceVertex points to. It has N_Faces entries of tuples of possibly various
 * sizes. For the hexahedron the array is given by
 * \f[
 * \{ \{0, 1, 2, 3\},  \{0, 4, 5, 1\},  \{1, 5, 6, 2\}, \{2, 6, 7, 3\},  \{0, 3,
 * 7, 4\}, \{4, 7, 6, 5\}\}
 * \f]
 * consisting of 6 entries of 4-tuples. As for the edges the ordering of the
 * tuples implies a local face numbering where the first tuple corresponds to
 * the first face and so on. The local numbering for the hexahedron implied by
 * the array above can be seen below.
 * \image html ShapeDesc_hexahedron-3.svg "Exploded view of hexahedron with vertex and face numbers"
 * Often always the same number of vertices build and edge, e.g. 3 for
 * tetrahedra and 4 for hexahedra. In general this has not to be the case, e.g.
 * a pyramid with a square base. The information about the number of vertices
 * that build an face is given by an array of size N_Faces to which \ref
 * FaceVertexLen points to. The ordering is in correspondence with the numbering
 * of the faces, i.e. FaceVertexLen[0] returns the number of vertices that are
 * involved in the first face. The largest value in this array is stored in \ref
 * MaxN_VpF. For the hexahedron these arrays are \f$
 * FaceVertexLen = \{4,4,4,4,4,4\} \f$ and \f$MaxN\_VpF = 4\f$. Related to these
 * description \ref FaceType points to an array of length N_Faces that consists
 * of \href Shapes , where each entry indicates the shape of the corresponding
 * face. For the hexahedron this array is given by
 * \f[
 * FaceTypes = \{Quadrangle, Quadrangle, Quadrangle, Quadrangle, Quadrangle,
* Quadrangle\}.
 * \f]
 *
 * With all this information the shape of the cell is uniquely defined. The
 * numbering of the vertices, edges and faces of the hexahedron example can be
 * seen together in this net of the hexahedron.
 * \image html ShapeDesc_hexahedron-4.svg "Net of hexahedron with vertex (black), edge (blue) and face numbers"
 *
 * Even though with these information the shape of a cell is uniquely defined,
 * ParMooN stores more information that come in handy from time to time.
 * First of all, an array of length N_Vertices to which \ref VertexEdge points
 * to is stored. This array can be seen as something like a transposed of the
 * array EdgeVertex points to. Each entry corresponds to one vertex using the
 * numbering of the vertices and returns the local edge numbers that meet in
 * that particular vertex. Or in other words, it lists the local edge numbers
 * for which a particular vertex is the start or the end point. For the
 * hexahedron this array is given by
 * \f[
 * VertexEdge = \{\{3, 0, 4\},  \{0, 1, 5\},  \{1, 2, 6\},  \{2, 3, 7\}, \{4,
 * 8,11\},  \{5, 9, 8\},  \{6,10, 9\},  \{7,11,10\}\},
 * \f]
 * e.g. the first vertex is involved in edges 3, 0 and 4.
 * The largest number of edges that contain the same vertex, i.e. the maximum of
 * the length of these tuples, is stored in \ref MaxN_VpF.
 *
 * Analogously, it is also stored which faces contain a particular vertex. The
 * pointer \ref VertexFace points to this array where each tuple corresponds to
 * the faces that a particular vertex is a part of. This is again the transposed
 * of FaceVertex. For the hexahedron the array is
 * \f[
 * VertexFace = \{\{0, 1, 4\},   \{0, 2, 1\},   \{0, 3, 2\},   \{0, 4, 3\}, \{1,
 * 5, 4\},   \{2, 5, 1\},   \{3, 5, 2\},   \{4, 5, 3\}\}.
 * \f]
 * The length of the largest tuple in this array is stored in \ref MaxN_FpV.
 *
 * As written above, the faces are described by using the vertices. In contrast
 * to this the faces could also be described using the edges, i.e. for each face
 * which edges build the particular face. This information is stored in an array
 * \ref FaceEdge points to.
 * Again the numbering corresponds to the numbering of the faces
 * and each tuple gives information about which edges are involved in building
 * the particular face. For the hexahedron this array is given by
 * \f[
 * FaceEdge =\{ \{0, 1, 2, 3\},  \{4, 8, 5, 0\},  \{5, 9, 6, 1\}, \{6,10, 7, 2\},  \{3,
 * 7,11, 4\},  \{11,10, 9,8\}\}.
 * \f]
 * Since again the faces can have different numbers of
 * edges, an array which has information about the length of each entry of
 * FaceEdge ist stored. \ref FaceEdgeLen points exactly to this array. The
 * largest of these length of tuples is stored in \ref MaxN_EpF. For the
 * hexahedron FaceEdgeLen looks exactly as FaceVertexLen given above and
 * \f$MaxN\_EpF=4\f$.
 *
 * Last but not least again an array transposed to the previous array is stored
 * to which \ref EdgeFace points to. It has length N_Edges and returns for each
 * edge the faces that meet in the edge. The largest of the length of these
 * tuples is stored in \ref MaxN_FpE. For the hexahedron the array is given by
 * \f[
 * EdgeFace = \{ \{1, 0\},  \{2, 0\},  \{3, 0\},  \{4, 0\}, \{4, 1\},  \{1, 2\},
 * \{2, 3\},  \{3, 4\}, \{5, 1\},  \{5, 2\},  \{5, 3\},  \{5, 4\} \}
 * \f]
 * and it is \f$MaxN\_FpE = 2\f$.
 *
 */
class TShapeDesc { protected:
    /** @brief type of shape */
    Shapes Type;

    /** @brief number of vertices */
    int N_Vertices;
    /** @brief number of edges */
    int N_Edges;
    /** @brief maximum number of edges per vertex */
    int MaxN_EpV;

    /** @brief number of faces (3D) */
    int N_Faces = 0;
    /** @brief maximum number of vertices per face */
    int MaxN_VpF = 0;
    /** @brief maximum number of faces per vertex */
    int MaxN_FpV = 0;
    /** @brief maximum number of edges per face */
    int MaxN_EpF = 0;
    /** @brief maximum number of faces per edge */
    int MaxN_FpE = 0;

    /** @brief number of joints */
    int N_Joints;

    /** @brief which vertices belong to one edge */
    const int *EdgeVertex;
    /** @brief which edges meet at a vertex */
    const int *VertexEdge;

    /** @brief which vertices are on one face */
    const int *FaceVertex = nullptr;
    /** @brief number of  vertices on one face */
    const int *FaceVertexLen = nullptr;
    /** @brief which edges are on one face */
    const int *FaceEdge = nullptr;
    /** @brief number of edges on one face */
    const int *FaceEdgeLen = nullptr;
    /** @brief which shapes have the faces got */
    const Shapes *FaceType = nullptr;
    /** @brief which faces meet at a vertex */
    const int *VertexFace = nullptr;
    /** @brief which faces meet at one edge */
    const int *EdgeFace = nullptr;

  public:
    // Constructor
    virtual ~TShapeDesc(){};

    // Methods
    /** @brief return the number of vertices */
    int GetN_Vertices() const
    { return N_Vertices; }
    /** @brief return the number of edges */
    int GetN_Edges() const
    { return N_Edges; }
    /** @brief return the number of joints */
    int GetN_Joints() const
    { return N_Joints; }

    /** @brief return the number of faces */
    int GetN_Faces() const
    { return N_Faces; }

    /** @brief return the shape */
    Shapes GetType() const
    { return Type; }

    /** @brief return the EdgeVertex array */
    int GetEdgeVertex(const int *&TmpEV) const
    {
      TmpEV = EdgeVertex;
      return 0;
    }

    /** @brief return the MaxN_EpV */
    int GetMaxN_EpV() const
    { return MaxN_EpV; }

    /** @brief return the EdgeVertex array */
    int GetVertexEdge(const int *&TmpVE) const
    {
      TmpVE = VertexEdge;
      return 0;
    }

    /** @brief return the FaceVertex array */
    int GetFaceVertex(const int *&TmpFV, const int *&TmpLen, int &MaxLen)
    const
    {
      TmpFV = FaceVertex;
      TmpLen = FaceVertexLen;
      MaxLen = MaxN_VpF;
      return 0;
    }
    /** @brief return the FaceEdge array */
    int GetFaceEdge(const int *&TmpFV, const int *&TmpLen, int &MaxLen) const
    {
      TmpFV = FaceEdge;
      TmpLen = FaceEdgeLen;
      MaxLen = MaxN_EpF;
      return 0;
    }

    /** @brief return the EdgeFace array */
    int GetEdgeFace(const int *&TmpEF, int &MaxLen) const
    {
      TmpEF = EdgeFace;
      MaxLen = MaxN_FpE;
      return 0;
    }

    /** @brief return the FaceType array */
    int GetFaceType(const Shapes *&TmpFT) const
    {
      TmpFT = FaceType;
      return 0;
    }

    /** @brief return diameter of a cell */
    virtual double GetDiameter(TVertex **Verts) const = 0;

    /** @brief return shortest of a cell */
    virtual double GetShortestEdge(TVertex **Verts) const = 0;

    /** @brief return the length of the cell defined with the reference map */
    virtual double GetLengthWithReferenceMap(TVertex **Verts) const = 0;

    /** @brief return measure of a cell */
    virtual double GetMeasure(TVertex **Verts) const = 0;
};

#endif
