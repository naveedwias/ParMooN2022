#ifndef __FESPACE__
#define __FESPACE__

#include "BoundaryJoint.h"
#include "Collection.h"
#include "Constants.h"
#include "Enumerations_fe.h"
#include "FiniteElement.h"
#include "HNDesc.h"
#include <string>
#include <array>
#include "ShapeDesc.h"

class THangingNode;

enum SpaceType { DiscP_USpace, DiscP_PSpace, ContP_USpace, ContP_PSpace,
                 Non_USpace };

/** @brief abstract base class for all finite element spaces, special spaces are
 * implemented in derived classes.
 * 
 * Here the mapping from local degrees of freedom (dof) to global dof is stored:
 * Let \f$n_C\f$ be the the number of cells and let each cell have an index
 * \f$i \in I_C:=\{0,...,n_C-1\}\f$. Furthermore, on a given cell \f$K_i\f$,
 * denote the set of local degrees of freedom on \f$K_i\f$ by
 * \f$\mathcal N_i\f$ with cardinality \f$ |\mathcal N_i |\f$. Each local degree
 * of freedom \f$\Phi_{i,j} \in\mathcal N_i\f$, is indexed via its cell index
 * \f$i\f$ and its local dof number 
 * \f$j \in I_i:=\{0,...,|\mathcal N_i|-1 \}\f$. Often we say local degree of
 * freedom \f$(i,j)\f$ when we mean \f$\Phi_{i,j}\f$. Finally, let
 * \f$N_{\text{dof}} \in \mathbb N\f$ be the number of global degrees of freedom.
 * 
 * The local-to-global map 
 * \f[
 *   \text{l2g} : \bigcup_{i \in I_C} \{ (i,j) | j \in I_i \} \to
 *                \{0,...,N_{\text{dof}}-1\}
 * \f]
 * assigns each local degree of freedom its global number. It is surjective.
 * 
 * Each local dof \f$(i,j)\f$ is assigned a global number \f$\text{l2g}(i,j)\f$.
 * In each cell each dof maps to exactly one global dof, i.e., the map 
 * \f$\text{l2g}(i,\cdot)\f$ is injective for each \f$i\in I_C\f$.
 * Dofs \f$(i_1,j_1)\f$ and \f$(i_2,j_2)\f$ from different cells \f$i_1\f$
 * and \f$i_2\f$ may map to the same global dof
 * \f$k = \text{l2g}(i_i,j_1) = \text{l2g}(i_2,j_2)\f$, which then means these 
 * dofs are identified in this fe space. For example think of a P1/Q1 space,
 * then the dof on a vertex is assigned to the same global dof number in all its
 * neighboring cells.
 * 
 * The map \f$\text{l2g}\f$ is represented as the member `get_global_dof()`, see
 * also `GetGlobalDOF()`.
 * 
 * This class therefore defines a (global) numbering of all degrees of freedom,
 * in particular it defines the dimension of the space, which is stored in the
 * member variable `N_DegreesOfFreedom`. Every finite element function in this
 * space can then be represented as a vector of that size.
 * 
 * Additionally the degrees of freedom are ordered:
 * 1. inner dofs
 * 2. non-Dirichlet boundary dofs (in particular Neumann, Robin)
 * 3. hanging dofs, i.e., dofs which are found near a hanging vertex
 * 4. Dirichlet dofs, i.e., dofs which are on a boundary where Dirichlet 
 *    conditions are prescribed.
 * Note that dofs which are both on a Dirichlet and on a Neumann boundary count
 * as Dirichlet dofs.
 * 
 * The number of dofs in each of the four sets of dofs can be obtained through
 * the respective member functions:
 * 1. `get_n_inner()`
 * 2. `get_n_non_dirichlet_boundary()`
 * 3. `get_n_hanging()`
 * 4. `get_n_dirichlet()`
 * 
 * Note that the sum of the four is equal to `get_n_dof()`.
 * Furthermore, there are members to get the number of dofs on Neumann and Robin
 * boundary parts, `get_n_neumann()` and `get_n_robin()`.
 * 
 * One often needs the combination of some of these sets of degrees of freedom.
 * This is what the methods `get_n_active()` and `get_n_active_non_hanging()`
 * are for.
 * 
 * If one wants to iterate over any of the four sets of dofs, use the methods
 * 1. `get_n_inner()`: indices \f$0 \le i \le\f$ `get_n_inner()` 
 *    correspond to the inner dofs
 * 2. `get_n_active_non_hanging()`, `get_n_inner()`: indices 
 *    `get_n_inner()`\f$\le i \le\f$`get_n_active_non_hanging()` are the
 *    ones on non-Dirichlet boundaries
 * 3. `get_n_active_non_hanging()`, `get_n_active()`: indices 
 *    `get_n_active_non_hanging()`\f$\le i \le\f$`get_n_active()` correspond to
 *    hanging nodes.
 * 4. `get_n_dof()`, `get_n_active()`: indices 
 *    `get_n_active()`\f$\le i \le\f$`get_n_dof()` are the Dirichlet dofs
 *
 * Here is an illustration for the degrees of freedom:
 * \image html fe_space_degrees_of_freedom.svg
 */
class TFESpace
{
  protected:
    /** @brief name of the space
     * 
     * This serves to computational purpose, instead it is mostly used to easily
     * distinguish between multiple spaces.
     */
    std::string Name;

    /** @brief Collection containing the cells used for building this space
     * 
     * This is the representation of the mesh on which this finite element space
     * is defined on.
     */
    const TCollection *Collection;

// =======================================================================
// general information on degrees of freedom and their numbers
// =======================================================================
    /** @brief number of all degrees of freedom
     * 
     * This is equal to the largest entry in `GlobalNumbers`+1 (plus 1 because
     * the starting index is 0 in c++).
     */
    int N_DegreesOfFreedom = 0;

    /// @name vectors defining the local-to-global map \f$\text{l2g}\f$
    /// @brief determine the global index of each degree of freedom.
    /// 
    /// These two vectors implement the local-to-global map \f$\text{l2g}\f$ as
    /// defined in this class' description. Here is an illustration for these
    /// two vectors \image html fe_space_local_to_global_map.svg
    ///@{
    /** @brief vector containing all global numbers of local degrees of freedom
     * for all elements.
     *
     * The length of this vector is the number of local degrees of freedom,
     * i.e., \f$ \sum_{i\in I_C} |\mathcal N_i| \ge N_{\text{dof}}\f$. Each
     * entry is the global dof of some specific local dof \f$(i,j)\f$, where
     * both \f$i\f$ and \f$j\f$ are determined through the vector `BeginIndex`.
     */
    std::vector<int> GlobalNumbers;

    /** @brief vector containing the indices in `GlobalNumbers` for each
     * element, size: n_cells+1.
     * 
     * This vector determines where in `GlobalNumbers` one can find the global
     * dof number of a given local dof \f$(i,j)\f$, namely at position 
     * `BeginIndex[i]+j`. Therefore the local-to-global mapping \f$\text{l2g}\f$
     * is given by \f$\text{l2g}(i,j) = \f$ `GlobalNumbers[BeginIndex[i]+j]`.
     * @note the number of local dofs on cell `i` is
     * `BeginIndex[i+1]-BeginIndex[i]`, see `get_n_local_dof()`
     * @note The first entry is always zero: `BeginIndex[0] = 0`;
     */
    std::vector<int> BeginIndex;
    ///@}

// =======================================================================
// storing the finite elements
// =======================================================================
    /** array containing the used elements */
    std::vector<FE_type> UsedElements;

    /** array with an element for each shape */
    std::array<FE_type, N_SHAPES> ElementForShape;
    
    /** @brief store all the elements which are used in this space */
    std::map<FE_type, FiniteElement> elements;

    /** array storing the fe for each element, if necessary */
    std::vector<FE_type> AllElements;
    
// =======================================================================
// counts and bounds for different types of degrees of freedom
// =======================================================================
    /** number of different boundary node type */
    static const int N_DiffBoundNodeTypes = N_BOUNDCOND;

    /** number of Dirichlet nodes */
    int N_Dirichlet = 0;

    /** number of nodes for each boundary node type */
    std::array<int, N_DiffBoundNodeTypes> N_BoundaryNodes = {};

    /** number of inner nodes 
     * 0 <= i < N_Inner for all inner degrees of freedom i */
    int N_Inner = 0;

    /** array of hanging nodes, the i-th entry corresponds to the global dof
     * i+get_n_active_non_hanging(). */
    std::vector<const THangingNode*> HangingNodeArray;
    
    /** array of hanging nodes in a particular order. One needs this to properly
     * adjust the matrix and right-hand side.
     * 
     * Note that a hanging dof i can as well be a coupled dof of another hanging
     * dof j. In that case the necessary matrix modifications need to be done in
     * a certain order, namely a multiple (e.g. 0.5 for P1) of the j-th row has
     * to be added to the i-th row before a multiple of the i-th row is added to
     * the coupling dofs of i. This array is constructed such that the j-th
     * hanging dof appears before the i-th.
     * The second parameter in the pair is the global index of this dof.
     */
    std::vector<std::pair<const THangingNode*, int>> sorted_hanging_nodes;
    
    /** hanging node descriptors */
    std::map<HNDesc, THNDesc *> HN_descriptors;

    /** True if this space contains discontinuous elements. Discontinuous means
     * not globally continuous, i.e also non-conforming finite elements that are
     * continuous at some points are defined as discontinuous.
     */
    bool is_discontinuous_space = false;
    
    /** @brief fill the array TFESpace::BeginIndex during construction */
    void fill_begin_index();
    
    /** find used elements */
    void FindUsedElements();
    
    /// helper method called during construction in 2D and 3D.
    /// This simply avoids double code, there still needs to be some refactoring
    /// here.
    void compute_various_numbers(
      const std::array<int, N_DiffBoundNodeTypes>& BoundaryUpperBound,
      const std::vector<THangingNode *>& VHN, const std::vector<int>& HNNumbers,
      bool flag_2d);

    /**
     * Constructor, setting name, description, cellgrid, and a lot of
     * dummy variables. Note that FESpace is more like an interface class,
     * only its daughter classes (with fixed dimension) are usable.
     *
     * @param[in] coll The cell grid (i.e. the finite element mesh).
     * @param[in] name The name of the space, used in printout etc.
     */
    TFESpace(const TCollection *coll, const std::string& name);
    
    /** @brief constructor for building a space with elements of order k
     * 
     * @param[in] coll collection of cells, i.e., the mesh, a pointer a stored
     * @param[in] name name for this finite element space
     * @param[in] k polynomial degree of finite elements
     * @param[in] dim dimension of the space (1, 2, or 3)
     */
    TFESpace(const TCollection *coll, const std::string& name, int k, int dim);
    
    /** @brief constructor for building a space with elements of order k
     * 
     * This constructor is mostly useful for Navier-Stokes problems.
     * 
     * @param[in] coll collection of cells, i.e., the mesh, a pointer a stored
     * @param[in] name name for this finite element space
     * @param[in] type specification of the type of space
     * @param[in] k polynomial degree of finite elements
     * @param[in] dim dimension of the space (1, 2, or 3)
     * 
     * @todo this constructor is only used when VMS projection is done. Can we
     * maybe get rid of it completely?
     */
    TFESpace(const TCollection *coll, const std::string& name, SpaceType type,
             int k, int dim);
    
    /** @brief constructor for building a space with the given elements 
     * 
     * Use this constructor if all elements are known explicitly.
     * 
     * @param[in] coll collection of cells, i.e., the mesh, a pointer a stored
     * @param[in] name name for this finite element space
     * @param[in] fes array of finite element types, size: `coll->GetN_Cells()`
     */
    TFESpace(const TCollection *coll, const std::string& name, FE_type *fes);

    /** @brief sets the variable TFESpace::is_discontinuous_space for given
     * finite elements
     *
     * @details
     * This function sets the variable TFESpace::is_discontinuous_space with
     * respect to the underlying finite elements. If there are more than one
     * finite elements then this variable is false if at least one of the
     * finite element types is not-continuous.
     *
     * @return TFESpace::is_discontinuous_space will be set to false if one of
     * the elements in TFESpace::UsedElements is not-continuous
     */
    void set_is_discontinuous_space();

  public:

    /** destrcutor */
    virtual ~TFESpace();

    /** return name */
    std::string GetName() const
    { return Name; }

    /** return number of cells in the triangulation used for building 
        this space */
    int GetN_Cells() const
    { return Collection->GetN_Cells(); }

    /** return the collection of this space */
    const TCollection *GetCollection() const
    { return Collection; }
    
    /** @brief get dimension of the vector basis function */
    int GetBaseVectDim() const;
    
    /** @brief return local-2-global map \f$\text{l2g}(i,\cdot)\f$ for a given 
     * cell with index `i`.
     * 
     * set `int * DOF = feSpace->GetGlobalDOF(i);` then the `j`-th local dof 
     * corresponds to the `DOF[j]`-th global degree of freedom, i.e.,
     * `DOF[j]`=\f$\text{l2g}(i,j)\f$.
     * 
     * @param i cell index
     * 
     * @warning it is not checked if `i` is a valid cell index.
    */
    const int* GetGlobalDOF(int i) const
    { return &GlobalNumbers[BeginIndex[i]];}
    
    /**
     * @brief the local-to-global map
     * 
     * return the global dof number of the `local_dof_index`-th local dof in 
     * cell `cell_index`.
     * @param cell_index index of a cell in the collection (mesh)
     * @param local_dof_index index of the dof within the cell
     * 
     * @warning for faster access the input parameters are not checked.
     */
    int get_global_dof(int cell_index, int local_dof_index) const
    { return GlobalNumbers[BeginIndex[cell_index]+local_dof_index];}
    
    /// @brief get the number of local degrees of freedom (dof) on cell `i`.
    unsigned int get_n_local_dof(unsigned int i) const;
    
    /// @brief get the maximum number of local degrees of freedom over all cells
    unsigned int get_max_n_local_dof() const;
    
    /// @brief find out if the given `dof` is in the cell with the given index
    bool is_dof_in_cell(int dof, int cell_index) const;

    /** return number of used elements */
    int GetN_UsedElements() const
    { return UsedElements.size(); }

    /** return number of all degrees of freedom */
    int get_n_dof() const
    { return N_DegreesOfFreedom; }

// =======================================================================
// counts and bounds for different types of degrees of freedom
// =======================================================================
    /** @brief number of inner dofs */
    int get_n_inner() const
    { return N_Inner; }

    /** @brief number of dofs on the boundary which are not on the Dirichlet 
     * boundary */
    int get_n_non_dirichlet_boundary() const;

    /** @brief number of hanging dofs. */
    size_t get_n_hanging() const
    { return HangingNodeArray.size(); }

    /** @brief number of dofs on the Dirichlet boundary */
    int get_n_dirichlet() const
    { return N_Dirichlet; }

    /** @return number of Neumann boundary dofs */
    int get_n_neumann() const
    { return N_BoundaryNodes[0]; }

    /** @brief number of Robin boundary dofs */
    int get_n_robin() const
    { return N_BoundaryNodes[1]; }

    /** @brief number of active degrees of freedom, including hanging dofs */
    int get_n_active() const
    { return get_n_dof() - get_n_dirichlet(); }

    /** @brief number of active degrees of freedom, excluding hanging dofs */
    int get_n_active_non_hanging() const
    { return get_n_inner() + get_n_non_dirichlet_boundary(); }


    const std::vector<std::pair<const THangingNode*, int>>
    get_sorted_hanging_nodes() const
    { return sorted_hanging_nodes; }

    /** @brief hanging node of a given index. Valid indices start with 0. */
    const THangingNode* get_hanging_node(unsigned int hanging_dof_index) const;


    /** return identifiers of used elements */
    const std::vector<FE_type> GetUsedElements() const
    { return UsedElements; }

    /** @brief return the type of finite element on a cell */
    FE_type get_fe_type(int i) const;

    /** @brief return the Finite Element on a given cell */
    const FiniteElement& get_fe(unsigned int cell_number) const;

    /** @brief get degree of basis functions on a given cell
     * 
     * Note that this only works correctly if the cell->GetCellIndex() returns
     * the correct cell index.
     */
    int getFEDegree(const TBaseCell *cell) const;
    
    /** @brief return the boundary condition on a given boundary joint. */
    virtual BoundCond get_boundary_condition(const BoundaryJoint& bd_joint) 
      const = 0;

    /** @brief write info on fespace into file */
    int Write(const std::string& filename);

    /** Check whether the space is discontinuous or not
     */
    bool is_discontinuous() const { return is_discontinuous_space; }

    /**
     * @brief print some information on this fe space
     * 
     * Depending on the verbosity level more or less is printed.
    */
    void info() const;

    /**
     * @brief returns element for shape depending on the order of the method and
     * the dimension
     *
     * @details
     * @param[in] k identifyer / number of finite element
     * @param[in] dim dimension of space
     * @return array of finite element types defining the shape of the
     * finite element space
     */
    static std::array<FE_type, N_SHAPES> get_element_for_shape(int k, int dim);
};

#endif
