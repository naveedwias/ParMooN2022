/** ************************************************************************ 
*
* @class     Mesh
* @brief     stores mesh arrays and allows for conversion
*
*            Store the mesh using several lists
*            (vertices, tria, quad, etc.), which are
*            implemented as separate struct (in this file)
*            Moreover, it contains a Boundary object that
*            can store information about the boundary (i.e. a PRM file)
*            Note: the Boundary must be initialized in order to
*            write the geometry in .GEO format.
* 
* @author    Alfonso Caiazzo 
* @date      21.03.16
 ************************************************************************  */

// A lot of the internals of this class are public, and there's a great deal of
// improper entanglement with TDomain in particular.
// TODO: figure out how much of this can be reasonably cleaned up and do so!

#include <vector>
#include <array>
#include <set>
#include <unordered_map>
#include <string>
#include "Boundary.h"

#ifndef __MESH__
#define __MESH__

struct meshNode
{
  double x,y,z;
  int reference;
};

///@brief a structure for edges 
struct meshEdge
{
  std::array<int, 2> nodes;
  int reference;
};

///@brief a structure for triangles
struct meshTriangle
{
  std::array<int, 3> nodes;
  int reference;
};

///@brief a structure for quadrilaterals
struct meshQuad
{
  std::array<int, 4> nodes;
  int reference;
};

///@brief a structure for tetrahedra
struct meshTetrahedron
{
  std::array<int, 4> nodes;
  int reference;
};

///@brief a structure for hexahedra
struct meshHexahedron
{
  std::array<int, 8> nodes;
  int reference;
};



/**  @brief mesh details */
class Mesh 
{
  public:

    /// @brief dimension
    /// @attention this is not necessarily the geometrical dimension
    /// (a mesh file could be written in dimension 3 even if
    /// the domain lies on a plane)
    unsigned int dimension;

    /// @brief nodes
    std::vector<meshNode> vertex;

    /// @brief edges
    std::vector<meshEdge> edge;

    /// @brief triangles
    std::vector<meshTriangle> triangle;

    /// @brief quadrilateral
    std::vector<meshQuad> quad;

    /// @brief flag for the case of two dimensional mesh with both
    /// triangles and quadrilaterals
    bool hasBothTriaAndQuads;

    /// @brief flag for mesh with tetrahedra
    bool is_tetramesh = false;

    /// @brief flag for mesh with hexahedra
    bool is_hexamesh = false;

    /// @brief tetrahedron
    std::vector<meshTetrahedron> tetra;

    /// @brief hexahedron
    std::vector<meshHexahedron> hexa;

    /// @brief boundary handler class
    Boundary boundary;

    /// @brief number of boundary faces
    int n_boundary_faces;

    /// @brief map face to neighboring tetra
    ///
    /// faceToTetra[i][0] and faceToTetra[i][1] contains the the indices
    /// of two elements sharing face i
    /// faceToTetra[i][0]=faceToTetra[i][1]=-1 if the face is on the boundary
    std::vector< std::vector<int> > faceToTetra;

    /// @brief map face to neighboring hexa
    ///
    /// faceToHexa[i][0] and faceToHexa[i][1] contains the the indices
    /// of two elements sharing face i
    /// faceToHexa[i][0]=faceToHexa[i][1]=-1 if the face is on the boundary
    std::vector< std::vector<int> > faceToHexa;

    /// @brief boundaryFacesMarker: 0 = inner face, > 0 = boundary face
    std::vector<int> boundaryFacesMarker;

    /// @brief Emtpy constructor initialize an empty mesh,
    /// a single file initialize the mesh structures (.mesh)
    /// while passing also a boundary file initializes also
    /// the boundary description.
    Mesh();

    explicit Mesh(const std::string& f);

    Mesh(const std::string& filename, const std::string& filenameBoundary);

    ~Mesh();

    /// @brief read mesh from a file
    /// @note supported formats: .mesh
    void readFromFile(const std::string& filename);

    ///@brief write mesh to a file .mesh
    void writeToMesh(const std::string& filename, bool includeInteriorFaces = true);

    /// @brief write mesh to a file .xGEO (extended ParMooN format)
    /// @param prmfile is an input prm file describing the boundary
    /// @attention the Boundary object must have been initialized
    /// (from a  prm file consistent with the geometry)
    /// @warning it works only in 2D at the moment
    void writeToGEO(const std::string& filename);

    /// @brief initialize the Boundary class reading a PRM file 
    /// @warning it works only in 2D at the moment
    void setBoundary(const std::string& PRM);

    /// @brief compute the number of boundary faces of the mesh
    ///
    /// Note: this functions uses the vector faceToTetra.
    /// If this has not been filled yet, it will be created by the function
    void computeNumberOfBoundaryFaces();

    /// @brief create inner faces is these are not written in the mesh file
    void createInnerFaces();

    /// @brief Try to fix boundary elements with too many facets on the
    /// boundary, e.g. by splitting edges.
    /// @warning Currently only implemented for tetrahedral meshes.
    void fixBoundaryElements();

    /// @brief find the index of the triangular face with indices a,b,c
    ///
    /// Given the vertex indices (a,b,c), this function uses the hash a+b+c
    /// to find the corresponing face.
    /// Note: the meshTrifaceHash vector must have been created before. If not,
    /// the funcion creates it.
    /// Return -1 (and a warning) if no face is found.
    int findTriFace(int a, int b, int c);

    /// @brief find the index of the quadrilateral face with indices a,b,c,d
    ///
    /// Given the vertex indices (a,b,c,d), this function uses the hash a+b+c+d
    /// to find the corresponing face.
    /// Note: the meshQuadfaceHash vector must have been created before. If not,
    /// the funcion creates it.
    /// Return -1 (and a warning) if no face is found.
    int findQuadFace(int a, int b, int c, int d);

    /// @brief display some info on screen
    void info();

    /// @brief number of nodes in the mesh
    unsigned int nPoints() { return vertex.size(); };

    /// @brief number of (inner) elements in the mesh
    unsigned int nElements()
    {
      if (hexa.size() + tetra.size())
      {
        return hexa.size() + tetra.size();
      }
      else
      {
        return triangle.size() + quad.size();
      }
    };

  private:

    class IndexLookup
    {
      public:
        bool contains(int hash) const;

        const std::vector<int>& at(int hash) const;

        void reserve(int hash, int capacity);

        void add(int value, int hash);

        void remove(int value, int hash);

        void clear();

        int size() const { return map.size(); };
      private:
        std::unordered_map<int, std::vector<int>> map;
    };

    /// @brief count the faces with the same hash
    ///
    /// This function is used to speed up the search of a nodes' face
    /// Once the "hash" array has been created, given three vertices (a,b,c),
    /// we can search the ID of the corresponding face by looking up only the
    /// faces with the same hash (a+b+c).
    void hashTriFaces();

    /// @brief count the faces with the same hash
    ///
    /// This function is used to speed up the search of a nodes' face
    /// Once the "hash" array has been created, given four vertices (a,b,c,d),
    /// we can search the ID of the corresponding face by looking up only the
    /// faces with the same hash (a+b+c+d).
    void hashQuadFaces();

    /// @brief create a list mapping each face to the neighboring tetra
    void createFaceToTetrahedraMap();

    /// @brief create a list mapping each face to the neighboring hexa
    void createFaceToHexahedraMap();

    /// @brief Try to fix boundary tetrahedra with two or more
    /// faces on the boundary
    void fixBoundaryTetrahedra();

    /// @brief Split a tetrahedron into four at its centroid
    /// @param c Global cell index
    void splitTetrahedronCentroid(int c);

    /// @brief Split an edge of a tetrahedral mesh,
    /// bisecting all adjacent cells
    /// @param a Global vertex index (1-based)
    /// @param b Global vertex index (1-based)
    int splitTetrahedraAtEdge(int a, int b);

    /// @brief Split a triangle of a tetrahedral mesh and mark all adjacent
    /// tetrahedra for splitting
    /// @param f Global face index
    /// @param i Local edge index
    /// @param new_node_index Index of the edge midpoint node (1-based)
    /// @param cellsToSplit Indices of tetrahedra to split
    void splitTriangleAtEdge(int f, int i, int new_node_index, std::set<int> &cellsToSplit);

    /// @brief Assuming any relevant triangles have already been split,
    /// split a tetrahedron along one of its edges.
    /// @param c Global cell index
    /// @param a Global vertex index (1-based)
    /// @param b Global vertex index (1-based)
    /// @param new_node_index Index of the edge midpoint node (1-based)
    void splitTetrahedronAtEdge(int c, int a, int b, int new_node_index);

    /// @brief Set up the edge-to-triangle mapping.
    void buildEdgeSums();

    /// @brief Remove all edge-to-triangle entries for a given face.
    /// @param f Global face index
    void removeEdgeSums(int f);

    /// @brief Add all edge-to-triangle entries for a given face.
    /// @param f Global face index
    void updateEdgeSums(int f);

    /// @brief Remove a triangle from the vertex sum lookup.
    /// @param f Global face index
    /// @param old_hash Previous vertex sum
    void unhashTriangle(int f, int old_hash);

    /// @brief Add a triangle to the vertex sum lookup.
    /// @param f Global face index
    void hashNewTriangle(int f);

    ///@brief routine for Element2D = meshQuad, meshTriangle
    template <class Element2D>
    void correct_numbering_of_vertices_if_needed(std::vector<Element2D>& elements_2d);

    template <class Element2D>
    bool first_three_vertices_are_numbered_clockwisely(const Element2D& element_2d);

    template <class Element2D>
    void reverse_clockwise_numbering(Element2D& element_2d);

    double get_z_component_of_cross_product(double* vec1, double* vec2)
    { return vec1[0] * vec2[1] - vec1[1] * vec2[0]; }

    /// @brief For each integer, stores a list of all triangles whose vertex
    /// indices sum to that integer.
    IndexLookup meshTrifaceHash;

    /// @brief For each integer, stores a list of all quads whose vertex
    /// indices sum to that integer.
    IndexLookup meshQuadfaceHash;

    /// @brief For each integer, stores a list of all triangles with edges
    /// whose vertex indices sum to that integer.
    IndexLookup edgeSumToTriangle;
};

#endif