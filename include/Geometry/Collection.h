/** ************************************************************************ 
*
* @class TCollection 
* @date  14.10.97
* @brief store cells in an array
* @author Gunar Matthies & Sashikumaar Ganesan
* @History: MPI methods (Sashikumaar Ganesan, 08.08.14)
   
****************************************************************************/

#ifndef __COLLECTION__
#define __COLLECTION__

#include <vector>
#include <map>
#include "Point.h"
class TBaseCell;
class TBoundEdge;
class TBoundFace;

/** @brief store cells in an array, used by cell iterators */
class TCollection
{
  protected:
    /** @brief number of cells stored */
    int N_Cells;

    /** @brief array containing the pointers to the cells */
    TBaseCell **Cells;

    /** @brief map each cell to its index within this collection.
     * This enables the method get_cell_index. It serves as a cache is therefore
     * mutable.
     */
    mutable std::map<const TBaseCell*, int> cell_to_index_map;
    
#ifdef  _MPI
    /** @brief Number of own cells (excluding Halo cells) */
    int N_OwnCells;
#endif
    
    struct ElementLists
    {
      //----- used in DataWriter ----
      std::vector<double> NodesCoords;
      std::vector<unsigned int> NodesReferences;
      std::vector< std::vector<unsigned int> > ElementNodes;
      std::vector<unsigned int> ElementReferences;
      std::vector<unsigned int> BdFacesNodes;
      std::vector<int> BdFacesReferences;
      std::vector<unsigned int> DomainVertexNumbers;
      unsigned int NLocVertices;
      bool empty() const { return NodesCoords.empty(); }
    };
    /// @brief store redundant data for faster access, mostly needed for 
    /// DataWriter. It serves as a cache is therefore mutable.
    mutable ElementLists element_lists;
    /// @brief create a list of nodes, vertices, elements.
    /// This method is 'const' because it only alters the mutable member 
    /// 'element_lists'.
    void createElementLists() const;
   
    
  public:
    /** @brief constructor */
    TCollection(int n_cells, TBaseCell **cells);

    /** @brief return number of cells */
    int GetN_Cells() const
    { return N_Cells; }

    /** @brief return Cell with index i in Cells-array */
    TBaseCell *GetCell(int i) const
    { return Cells[i]; }

    /** @brief destructor: delete arrays */
    ~TCollection();

    /** @brief get maximal and minimal diameter */
    int GetHminHmax(double *hmin, double *hmax) const;
    
    /** @brief tells if the mesh is Delaunay or not
     * @note Currently implemented only in 2D, 3D returns false 
     * @note Hexahedra and Quad grids return false */
    bool IsDelaunay() const;
    
    /** @brief find out if one of the cells in this collection is a quad/hex. */
    bool has_quad_cell() const;
    bool has_hex_cell() const
    { return has_quad_cell(); }
    
    /** @brief find out if one of the cells in this collection is a tria/tetra. */
    bool has_tria_cell() const;
    bool has_tetra_cell() const
    { return has_tria_cell(); }

    /** @brief return index of cell in Cells-array */
    int get_cell_index(const TBaseCell *cell) const;
    
    /** @brief set the clipboards of all cells to the cell's index in this
     * collection. Furthermore, all neighbors which are not in this collection 
     * get a clipboard value of -1.
     */
    void mark_all_cells() const;
    
    /** @brief find out if there are hanging vertices in this grid. */
    bool includes_hanging_vertices() const;
    
#ifdef  _MPI
    void SetN_OwnCells(int n_OwnCells)
     { N_OwnCells = n_OwnCells; }

    int GetN_OwnCells() const
     { return N_OwnCells; }

    int GetN_HaloCells() const
     { return (N_Cells - N_OwnCells); }

    /**
     * Find the lowest-number process which contains the given point (x,y,z)
     * in an OwnCell.
     *
     * Throws an error if the point was not found anywhere.
     *
     * @param x x value
     * @param y y value
     * @param z z value
     * @return The lowest number process of those processes which contain
     * the given point in an "OwnCell".
     */
    int find_process_of_point(double x, double y, double z) const;
#endif
    
    /**
     * @brief find the index of the cell in which the point p is in
     * 
     * This method first checks if p lies in the given cell. If yes, nothing is
     * done. If no, all neighbors are checked. If p is not in any of the
     * neighbors, it is sought in all cells.
     * 
     * This method will alter `current_cell_nr` to make sure that the point is 
     * in the cell with that index in this collection.
     * 
     * This method is much faster if you repeatedly have to find the cell of a
     * point where you know the next point is close to the previous one. Then
     * you can call this for point p1 and the the second call with p2 (close to
     * p1) is very cheap, because one does not search through all cells.
     * 
     * For example if you want to evaluate a finite element function along a 
     * line with many points, then you know that the next point is in the same 
     * cell as the previous one or in one of its neighbors.
     */
    void find_cell(parmoon::Point p, int& current_cell_nr) const;

   ///@brief write the geometry in .mesh format
   int writeMesh(const char *meshFileName, unsigned int VERTEX_OFFSET=0) const;
    
   /**@brief Write a list of boundary edges
    */
   void get_edge_list_on_component(int i,std::vector<TBoundEdge*> &edges) const;
   ///@todo it is better to return the vector?
   // std::vector<TBoundEdge*> get_edge_list_on_component(int i);
   // ------------------------------------------------
   
   void get_boundary_edge_list(std::vector<TBoundEdge*> &edges) const;
   
   //New LB 11.10.18
#ifdef __3D__
   void get_face_list_on_component(int boundary_component_id,
                                   std::vector<TBoundFace*> &faces) const;
#endif

   //####################################################
   //--------------- used in DataWriter -----------------
   //####################################################
   
   unsigned int GetN_Vertices() const;
   
   unsigned int GetN_BdFaces() const;
   
   double GetCoord(unsigned int vert) const;
   
   /// @brief return the vector of the coordinates of all vertices.
   /// The ordering in 2D is vertex_0_x, vertex_0_y, vertex_1_x, vertex_1_y, ...
   /// In 3D similarly with a third component for each verex.
   const std::vector<double>& get_nodes_coords() const;
   
   unsigned int GetBdFacesNode(unsigned int face) const;
   
   unsigned int GetGlobalVerNo(unsigned int cell, unsigned int locvert) const;
   
   unsigned int GetNLocVertices() const;
};

#endif
