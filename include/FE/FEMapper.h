#ifndef INCLUDE_FE_FEMAPPER_H
#define INCLUDE_FE_FEMAPPER_H

constexpr int HANGINGNODE = -5;
constexpr int FIRSTMARK = -10;

#include <HangingNode.h>
#include <string>
#include <vector>

/// @brief a unique identifier for each FEMapper
enum class FEMapper_type;


/** @brief find out which of the given local degress of freedom are equivalent
 *  to the same global degree of freedom 
 * 
 * This class in its derived classes are only needed during the construction of
 * finite element spaces. Its methods are called for joints and try to identify
 * degrees of freedom which are in different cells but should be the same global
 * dof. The setting is therefore always a joint which has two neighbor cells.
 * 
 * @ruleof0
 */
class FEMapper
{
  protected:

    /** @brief the type identifier for this mapper */
    FEMapper_type type;

    /** @brief name for the mapper */
    std::string name;

    /** @brief some words describing the mapper */
    std::string description;

    /** @brief number of local degrees on the first side */
    int n_dof0;

    /** @brief number of local degrees on the second side */
    int n_dof1;

    /** @brief number of formally different degrees which are to be identified
     */
    int n_pairs;

    /** @brief which pairs of local degrees matching in an array[][][2]
     * 
     * In 2D there is only one set of pairs. In 3D there are 3 (tetrahedra) or
     * 4 (hexahedra). That's why we need the outer std::vector.
     */
    std::vector<std::vector<std::pair<int,int>>> pairs;

    /** @brief memory for internal storage.
     * 
     * Only to avoid reallocation. Its size is `n_dof0+n_dof1`.
     */
    mutable std::vector<int> aux;

    /** @brief map the two given degrees of freedom
     * 
     * This changes entries in `GlobalNumbers` such that the two dofs `dof0` and
     * `dof1` are identified.
     * 
     * @param GlobalNumbers should be TFESpace::GlobalNumbers
     * @param dof0 dof which is to be identified with `dof1`
     * @param dof1 dof which is to be identified with `dof0`
     * @param counter
     */
    void map_single_dof(int *GlobalNumbers, int dof0, int dof1, int &counter)
      const;

    /** @brief helper method which is used multiple times in routines of this 
     *         class.
     */
    int find_in_global_numbers(const int *GlobalNumbers, int w) const;
    
    /** @brief set the marks in `GlobalNumbers` for a boundary dof
     * 
     * @param GlobalNumbers should be `TFESpace::GlobalNumbers`
     * @param dof position of this dof in the `TFESpace::GlobalNumbers` vector
     */
    void map_single_boundary_dof(int *GlobalNumbers, int dof, int &BoundCounter,
                                 std::vector<int>& hanging_numbers) const;

    /** @brief protected default constructor used in derived classes. */
    FEMapper() = default;

  public:
    /** @brief construct an FEMapper of the given type. 
     * 
     * Note that certain types are not supported here and need to be constructed
     * in a derived class, in particular these are the ones with "1Reg" in their
     * type name.
     */
    FEMapper(FEMapper_type t);

    /** @brief delete all members, no memory on heap is used.*/
    virtual ~FEMapper() = default;

    /** @brief default copy constructor (currently never used) */
    FEMapper(const FEMapper&) = default;
    /** @brief default move constructor */
    FEMapper(FEMapper&&) = default;
    /** @brief default copy assignment */
    FEMapper & operator=(const FEMapper&) = default;
    /** @brief default move assignment */
    FEMapper& operator=(FEMapper&&) = default;

    /** @brief return name of mapper */
    std::string get_name() const
    { return name; }

    /** @brief return description of mapper */
    std::string get_description() const
    { return description; }

    /** @brief map the given local degrees of freedom
     * 
     * This method is used for inner joints. The neighbors are K0 and K1.
     * 
     * @param GlobalNumbers should be TFESpace::GlobalNumbers
     * @param I_K0 should be TFESpace::BeginIndex[K0]
     * @param I_K1 should be TFESpace::BeginIndex[K1]
     * @param Indices0 the cell-local indices of the dofs on this joint, seen
     *        from cell K0
     * @param Indices1 the cell-local indices of the dofs on this joint, seen
     *        from cell K1
     * @param counter
     * @param type should be left 0 in 2D, in 3D this the map type index of this
     *        joint (see `TJoint::GetMapType()`)
     */
    void map(int *GlobalNumbers, int I_K0, int I_K1, const int *Indices0,
             const int *Indices1, int &counter, int type = 0) const;

    /** @brief "map" the given dof on a boundary joint
     * 
     * similar to `FEMapper::map()` but on the boundary there is no neighbor and
     * therefore no need to identify dofs.
     */
    void map_boundary(int *GlobalNumbers, int I_K, const int *Indices,
                      int &BoundCounter,
                      std::vector<int>& hanging_numbers) const;


    /** @brief "map" the given dof on a boundary edge
     * 
     * similar to `FEMapper::map_boundary()` but only for edges. This is only
     * necessary in mpi mode, where it can happen that an edge is on the 
     * boundary, while no joint on this boundary is part of a (process') mesh.
     */
    void map_boundary_edge(int N_EdgeDOF, int *GlobalNumbers, int I_K,
                           const int *Indices, int &BoundCounter,
                           std::vector<int>& hanging_numbers) const;

    /** @brief "map" the given dof on a boundary vert 
     * 
     * similar to `FEMapper::map_boundary_edge()` but only for vertices. This is
     * only necessary in mpi mode, where it can happen that a vertex is on the
     * boundary, while no edge on this boundary is part of a (process') mesh.
     * Here we assume that there is only one dof at a given vertex.
     */
    void map_boundary_vertex(int *GlobalNumbers, int I_K, int Index,
                             int &BoundCounter,
                             std::vector<int>& hanging_numbers) const;
};

enum class FEMapper_type {
  // 2D mapper types
  C0_2_C0_2_2D,
  C1_2_C1_2_2D,
  C2_2_C2_2_2D,
  C3_2_C3_2_2D,
  C4_2_C4_2_2D,
  C5_2_C5_2_2D,
  C6_2_C6_2_2D,
  C7_2_C7_2_2D,
  C8_2_C8_2_2D,
  C9_2_C9_2_2D,
  N1_2_N1_2_2D,
  N2_2_N2_2_2D,
  N3_2_N3_2_2D,
  N4_2_N4_2_2D,
  N5_2_N5_2_2D,
  // ONE regular grid, same pattern on both sides
  C0_2_C0_2_1Reg_2D,
  C1_2_C1_2_1Reg_2D,
  C2_2_C2_2_1Reg_2D,
  C3_2_C3_2_1Reg_2D,
  C4_2_C4_2_1Reg_2D,
  C5_2_C5_2_1Reg_2D,
  N1_2_N1_2_1Reg_2D,
  N2_2_N2_2_1Reg_2D,
  N3_2_N3_2_1Reg_2D,
  N4_2_N4_2_1Reg_2D,
  N5_2_N5_2_1Reg_2D,
  // ONE regular grid, different pattern
  N2_2_N1_2_1Reg_2D, // coarse: second order, fine: first order
  N3_2_N2_2_1Reg_2D, // coarse: third order,  fine: second order
  N4_2_N3_2_1Reg_2D, // coarse: fourth order, fine: third order
  N5_2_N4_2_1Reg_2D, // coarse: fifth order,  fine: fourth order
  // 3D mapper types
  D_D_3D,
  P1_P1_3D,
  P2_P2_3D,
  P3_P3_3D,
  P2B_P2B_3D,
  Q1_Q1_3D,
  Q2_Q2_3D,
  Q3_Q3_3D,
  Q4_Q4_3D,
  N1_N1_3D,
  N2_N2_3D,
  N3_N3_3D,
  N4_N4_3D,
  NP1_NP1_3D,
  NP2_NP2_3D,
  // 1-regular grid, same pattern on both sides
  P1_P1_1Reg_3D,
  P2_P2_1Reg_3D,
  P3_P3_1Reg_3D,
  Q1_Q1_1Reg_3D,
  Q2_Q2_1Reg_3D,
  Q3_Q3_1Reg_3D,
  NP1_NP1_1Reg_3D,
  NP2_NP2_1Reg_3D
};


std::ostream& operator<<(std::ostream& out, const FEMapper_type t);

#endif // INCLUDE_FE_FEMAPPER_H
