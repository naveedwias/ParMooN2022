#ifndef ASSEMBLE_DG
#define ASSEMBLE_DG

class TAuxParam2D;
#include "FiniteElement.h"
#include "LocalAssembling.h"
#ifdef __2D__
#include "BoundEdge.h"
  using BoundaryFacet = TBoundEdge;
#else
#include "BoundFace.h"
  using BoundaryFacet = TBoundFace;
#endif

/** Number of local coefficients of a PDE. Increase if needed. 40 is just a
 * value that is sure to be large enough right now.
 */
constexpr int n_local_coeff = 40;

/**
 * @brief constructs a database with needed parameters for DG method.
 * @details
 * This function construct a default DG database. This includes symmetry_DG
 * specifying whether SIPG, IIPG, or NIPG shall be used, face_sigma_DG, i.e. the
 * factor in front of the penalty terms arising from the discretization of the
 * diffusive term, and eta_upwind_DG with which the amount of upwinding can be
 * controlled.
 **/
ParameterDatabase default_dg_database();

/**
 * @brief assembles the DG terms of a discontinuous Galerkin discretization for
 * CD problems.
 * @details This function assembles the edge contributions of discontinuous
 * Galerkin discretization of convection-diffusion-reaction problems, see e.g.
 * <a href="https://doi.org/10.1016/j.cam.2021.113487">this paper</a>, equations (4)-(6).
 * \note Note that only scalar cases are implemented yet, i.e. only a single
 * finite element space is allowed.
 * \note Be aware that for this method it might be good to choose the grid such
 * that on Dirichlet edges \f$\mathbf{b}\cdot \nu\f$ does not change its sign
 * along the edge, where \f$\mathbf{b}\f$ is the convection term and \f$\nu\f$
 * is the outer unit normal vector. If the sign changes, the assembled matrix is
 * inaccurate for this entry. This inaccuracy might be not significant, but
 * please have this in mind.
 * \note This code may allows for different polynomial degrees on adjacent
 * cells. However, this feature was not tested. So before using it, test it!
 *
 * \todo specify input and output of this function and all other assemble
 * routines. For now please refer to the code of
 * ConvectionDiffusion< d >::call_assembling_routine.
 */
template< int d >
void Assemble_DG(const typename Template_names<d>::CoeffFct& Coeff, int
    n_fespaces, const typename Template_names<d>::FESpace** fespaces, int
    n_sqmatrices, typename Template_names<d>::SquareMatrixD** sqmatrices, int
    n_matrices, typename Template_names<d>::MatrixD **/*matrices*/, int n_rhs,
    double** rhs, typename Template_names<d>::BoundaryConditionFunction**
    BoundaryConditions, typename Template_names<d>::BoundaryValuesFunction**
    BoundaryValues, const ParameterDatabase& param_db, 
    const std::string system_type, int matrix_flag = -1, 
    typename Template_names<d>::SquareMatrixD **sqMatA22 = nullptr);

/** @brief checks if quadrature points on a given edge are in the same order
 * seen from the two adjacent cells.
 * @details
 * This function checks whether the quadrature points on a given edge are in the
 * same order from the perspective of the two adjacent cells.
 * @param[in] cell original cell
 * @param[in] cell_neigh neighboring cell
 * @param[in] joint_nr number of facet in original cell
 * @param[in] joint_nr_neigh number of the same facet in the neighboring cell
 * @return a boolean indicating if they are in correct order or not
 */
template< int d >
bool is_correct_order_of_quad_points( const TBaseCell* cell, const TBaseCell*
    cell_neigh, const int& joint_nr, const int& joint_nr_neigh);

/** @brief computes an list of indices that maps (d-1) dimensional
 * quadrature points on a facet of a cell to the same quadrature points on the
 * same facet from the neighboring cell
 * @details
 * Assume two neighboring cells share a facet and a (d-1) dimensional quadrature
 * formula on this facet is given. The d dimensional quadrature points seen from
 * both cells with respect to this quad formula are the same but may not be in
 * the same order. This function computes an index map that map the index j of
 * the j-th quadrature point in the cell to the index that the quadrature point
 * has in the neighboring cell.
 *
 * With this map it holds
 * quad_pt[j] = quad_pt_neigh[index_map[j]]
 * if quad_pt and quad_pt_neigh are the d dimensional quad points of the facet
 * of the cell and the neighboring cell, j an index, and index_map is the map
 * that is computed by this function.
 * @param[in] cell original cell
 * @param[in] cell_neigh neighboring cell
 * @param[in] joint_nr number of facet in original cell
 * @param[in] joint_nr_neigh number of facet in neighboring cell
 * @param[in] fin_ele finite element on the cell
 * @param[in] fin_ele_neigh finite element on the neighboring cell
 * @param[in] quad_formula_facet quadrature formula on the facet
 * \return a vector of integers that are the map of the index of the
 * quadrature point from the cell to the index the quadrature point has in the
 * neighboring cell
 */
template <int d >
std::vector<int> compute_index_map_for_facet_quad_points_of_neigh( const
    TBaseCell* cell, const TBaseCell* cell_neigh, const int& joint_nr, const
    int& joint_nr_neigh, const FiniteElement& fin_ele, const FiniteElement&
    fin_ele_neigh, const TQuadFormula* quad_formula_facet );

/** @brief computes the coefficients of PDE in quadrature points on the facet
 * @details This function computes the coefficients of the PDE evaluated at
 * the quadrature points of a d-1 dimensional quadrature rule for a given facet.
 * @param[in] collection collection of cells
 * @param[in] cell cell in which the joint is a part of
 * @param[in] joint_nr local index of joint with respect to cell
 * @param[in] finite_element finite element on this cell
 * @param[in] quad_form (d-1) dimensional quadrature formula for this edge
 * @param[in] parameters parameter function
 * @param[in] coeff_fct coefficient function
 * @param[out] coefficients array sufficiently large where the evaluation of the
 * coefficient function for all quadrature points are stored in
 * @result evaluate the coefficient function in quadrature points on the given
 * edge and store in coefficients
 */
template< int d >
void compute_coefficients( const TBaseCell* cell, const int& joint_nr, const
    FiniteElement& finite_element, const TQuadFormula* quad_form, const
    typename Template_names<d>::CoeffFct& coeff_fct, double** coefficients);

/**
 * @brief get the type of the boundary condition on a given joint.
 * @details Return the type of a given boundary condition on a specific joint.
 * @param[in] joint joint on which the type of boundary condition shall be
 * determined
 * @param[in] boundary_condition boundary condition on this joint from which the
 * type shall be determined
 * @return the type of the boundary condition on this edge
 */
template< int d >
BoundCond get_type_of_boundary_condition(const BoundaryFacet* joint,
    typename Template_names<d>::BoundaryConditionFunction* boundary_condition);

/**
 * @brief computes the local contributions of the DG bilinear form for CD problems.
 * @details This function computes the edge-wise local contribution of the DG
 * bilinear form for convection-diffusion-reaction problems. It calculates both
 * contributions of boundary edges and interior edges. It returns a local matrix
 * of local evaluations of the bilinear form.
 * @param[in] param_db parameter database including values for symmetry_DG,
 * face_sigma_DG and eta_upwind_DG
 * @param[in] val_bas_fct pointer to vector of values of the basis functions of the cell
 * evaluated at the quadrature points
 * @param[in] val_bas_fct_dx pointer to vector of values of the derivative wrt. x of the
 * basis functions of the cell evaluated at the quadrature points
 * @param[in] val_bas_fct_dy pointer to vector of values of the derivative wrt. y of the
 * basis functions of the cell evaluated at the quadrature points
 * @param[in] val_bas_fct_dz pointer to vector of values of the derivative wrt. z of the
 * basis functions of the cell evaluated at the quadrature points. In 2D this
 * can be the nullptr.
 * @param[in] val_bas_fct_neigh see \p val_bas_fct but for basis functions of
 * neighbor
 * @param[in] val_bas_fct_neigh_dx see \p val_bas_fct_dx but for basis
 * functions of neighbor
 * @param[in] val_bas_fct_neigh_dy see \p val_bas_fct_dy but for basis
 * functions of neighbor
 * @param[in] val_bas_fct_neigh_dz see \p val_bas_fct_dz but for basis
 * functions of neighbor
 * @param[in] joint_normal array of components of the joint_normal vector of the underlying
 * edge
 * @param[in] joint_diam length of edge
 * @param[in] quad_form_1D 1D quadrature formula of sufficiently high degree
 * for integration along edges
 * @param[in] coefficients array of coefficients of underlying PDE evaluated at
 * the quadrature points
 * @param[in] is_boundary_edge boolean to distinguish between inner and boundary
 * edges
 * @param[out] loc_mat array of size 2*ndofs x 2*ndofs for interior edges and of
 * size ndofs x ndofs for boundary edges where the local contributions are
 * written into. The rows correspond to test functions and the columns to ansatz
 * functions. The first \p n_dofs rows / columns correlate to dofs from the
 * cell, afterwards to dofs from the neighbor.
 * @return matrix of local edge-wise contributions of the underlying edge for
 * all combination of test- and ansatz functions of dofs of cell and its
 * potential neighbor along the edge
*/
template< int d>
void compute_DG_bilinear_form_local(const ParameterDatabase& param_db, const int& base_vec_dim, const
    std::vector<double>* val_bas_fct, const std::vector<double>* val_bas_fct_dx,
    const std::vector<double>* val_bas_fct_dy, const std::vector<double>*
    val_bas_fct_dz, const std::vector<double>* val_bas_fct_neigh, const
    std::vector<double>* val_bas_fct_neigh_dx, const std::vector<double>*
    val_bas_fct_neigh_dy, const std::vector<double>* val_bas_fct_neigh_dz, const
    std::vector<int>& index_map, const double* joint_normal, const double&
    joint_diam, const TQuadFormula* quad_form_joint, const std::vector<double>&
    orig_weights, double** coefficients, const bool& is_boundary_edge, double**
    loc_mat);

/**
 * @brief is a wrapper for boundary edges to call compute_DG_bilinear_form_local()
 * @details Boundary edges don't have neighbor and therefore do not need this as
 * a input parameter. To simplify the code one can also call this wrapper
 * function. See also compute_DG_bilinear_form_local().
 */
template< int d>
void compute_DG_bilinear_form_local(const ParameterDatabase& param_db, const int& base_vec_dim, const
    std::vector<double>* val_bas_fct, const std::vector<double>* val_bas_fct_dx,
    const std::vector<double>* val_bas_fct_dy, const std::vector<double>*
    val_bas_fct_dz, const double* joint_normal, const double& joint_diam, const
    TQuadFormula* quad_form_joint, const std::vector<double>& orig_weights,
    double** coefficients, double** loc_mat);
    
/**
* @brief computes the local contributions of the DG bilinear form for Navier-Stokes problems.
*/    
 
template< int d>
void compute_DG_bilinear_form_local_NSE(const ParameterDatabase& param_db, const int& base_vec_dim, const
    std::vector<double>* val_bas_fct, const std::vector<double>* val_bas_fct_dx,
    const std::vector<double>* val_bas_fct_dy, const std::vector<double>*
    val_bas_fct_dz, const std::vector<double>* val_bas_fct_neigh, const
    std::vector<double>* val_bas_fct_neigh_dx, const std::vector<double>*
    val_bas_fct_neigh_dy, const std::vector<double>* val_bas_fct_neigh_dz, const
    std::vector<int>& index_map, const double* joint_normal, const double&
    joint_diam, const TQuadFormula* quad_form_joint, const std::vector<double>&
    orig_weights, double** coefficients, const bool& is_boundary_edge, double**
    loc_mat, const int* global_dof, const int* global_dof_neigh, const TBaseCell*
    cell, int cell_no, int egde_no);
/**
 * @brief is a wrapper for boundary edges to call compute_DG_bilinear_form_local_NSE()
 * @details Boundary edges don't have neighbor and therefore do not need this as
 * a input parameter. To simplify the code one can also call this wrapper
 * function. See also compute_DG_bilinear_form_local_NSE().
*/
    
template< int d>
void compute_DG_bilinear_form_local_NSE(const ParameterDatabase& param_db, const int& base_vec_dim, const
    std::vector<double>* val_bas_fct, const std::vector<double>* val_bas_fct_dx,
    const std::vector<double>* val_bas_fct_dy, const std::vector<double>*
    val_bas_fct_dz, const double* joint_normal, const double& joint_diam, const
    TQuadFormula* quad_form_joint, const std::vector<double>& orig_weights,
    double** coefficients, double** loc_mat, const int* global_dof);    
    
/**
 * @brief computes the local contributions of the DG right-hand side for CD
 * problems.
 * @details This function computes the edge-wise local contribution of the DG
 * linear form for convection-diffusion-reaction problems. It returns a local
 * vector of local evaluations of the discrete right-hand side for all
 * test functions of a cell.
 * @param[in] param_db parameter database including values for symmetry_DG,
 * face_sigma_DG and eta_upwind_DG
 * @param[in] val_bas_fct vector of values of the basis functions of the cell
 * evaluated at the quadrature points
 * @param[in] val_bas_fct_dx vector of values of the derivative wrt. x of the
 * basis functions of the cell evaluated at the quadrature points
 * @param[in] val_bas_fct_dy vector of values of the derivative wrt. y of the
 * basis functions of the cell evaluated at the quadrature points
 * @param[in] joint_normal array of components of the joint_normal vector of the underlying
 * edge
 * @param[in] joint_diam length of edge
 * @param[in] quad_form_1D 1D quadrature formula of sufficiently high degree
 * for integration along edges
 * @param[in] coefficients array of coefficients of underlying PDE evaluated at
 * the quadrature points
 * @param[in] boundary_value function returning the boundary value at a given
 * component at a given point
 * @param[out] loc_rhs array of size ndofs where the local contributions are
 * written into.
 * @return vector of local edge-wise contributions of the underlying edge for
 * all test functions of dofs of cell using the parameters given in param_db
 */
template< int d >
void compute_DG_rhs_local(const ParameterDatabase& param_db, const TBaseCell*
    cell, const int& joint_index, const FiniteElement& finite_element, const
    std::vector<double>* val_bas_fct, const std::vector<double>* val_bas_fct_dx,
    const std::vector<double>* val_bas_fct_dy, const std::vector<double>*
    val_bas_fct_dz, const double* joint_normal, const double& joint_diam, const TQuadFormula*
    quad_form_facet, const std::vector<double>& orig_weights, double**
    coefficients, const typename Template_names<d>::BoundaryValuesFunction*
    boundary_value, double* loc_rhs);

/**
* @brief computes the local contributions of the DG right-hand side for Navier-Stokes
 * problems.
*/
template< int d >
void compute_DG_rhs_local_NSE(const ParameterDatabase& param_db, TBaseCell*
    cell, const int& joint_index, const FiniteElement& finite_element, const
    std::vector<double>* val_bas_fct, const std::vector<double>* val_bas_fct_dx,
    const std::vector<double>* val_bas_fct_dy, const std::vector<double>*
    val_bas_fct_dz, const double* joint_normal, const double& joint_diam, const
    TQuadFormula* quad_form_facet, const std::vector<double>& orig_weights,
    double** coefficients, const typename
    Template_names<d>::BoundaryValuesFunction* boundary_value, double* loc_rhs);

/**
 * @brief adds local matrix and right-hand side to global ones.
 * @details This function adds local contributions to global contributions both
 * for the system matrix and the right-hand side. Local contributions can be
 * cell-wise or edge-wise, i.e. for dofs of a cell and a neighbor.
 * @param[in] n_dofs number of local dofs on cell
 * @param[in] n_dofs_neigh number of local dofs on neighbor cell
 * @param[in] global_dof global dof numbers of local dofs of cell, see
 * TFESpace::GetGlobalDOF()
 * @param[in] global_dof_neigh see \p global_dof but for neighbor cells. If this
 * equals the nullptr only the first \p n_dofs entries are considered
 * @param[in] loc_mat local matrix with structure from
 * compute_DG_bilinear_form_local()
 * @param[in] loc_rhs local right-hand side, see also compute_DG_rhs_local(). It
 * his is the nullptr nothing happens to the right-hand side
 * @param[out] global_matrix global matrix of the system
 * @param[out] global_rhs global right-hand side
 * @return adds the values of \p loc_mat and conditionally of \p loc_rhs to the
 * correct positions in the global matrix and global right-hand side
 */
template< int d >
void copy_local_to_global(const int& n_dofs, const int& n_dofs_neigh, const int*
    global_dof, const int* global_dof_neigh, double** loc_mat, double* loc_rhs,
    typename Template_names<d>::SquareMatrixD* global_matrix, double*
    global_rhs, double** loc_mata22=nullptr, 
    typename Template_names<d>::SquareMatrixD*
    global_matrixA22 = nullptr);

/**
 * @brief is a wrapper function to call copy_local_to_global() which can be used
 * if only the dofs of a cell and not of a neighbor cell play a role. See also
 * copy_local_to_global().
 */
template< int d >
void copy_local_to_global(const int& n_dofs, const int* global_dof, double**
    loc_mat, double* loc_rhs, typename Template_names<d>::SquareMatrixD*
    global_matrix, double* global_rhs);

/**
 * @brief computes the values of basis functions and their derivatives in
 * quadrature points on a joint of a cell
 *
 * @details In DG bilinear forms integrals of the basis functions and their
 * derivatives along facets are computed. To prepare these integrals, this
 * function computes the evaluation of the basis functions and their derivatives
 * in quadrature points on the facet of a given cell. The quadrature points come
 * from a (d-1) dimensional quadrature formula on a reference domain.
 *
 * @param[in] finite_element finite element on the given cell that determines
 * the basis functions on the cell
 * @param[in] cell pointer to a cell which contains the facet of interest
 * @param[in] joint_index number of joint in the cell for which the basis functions
 * shall be evaluated
 * @param[in] base_vec_dim the dimension of the basis functions
 * @param[in] quad_formula (d-1) dimensional quadrature formula for the facet
 * @param[out] original_values vector of values of basis functions in quadrature
 * points. The values for all basis functions in the j-th quadrature point are
 * stored in original_values[j]. To access the value of the i-th basis function
 * in the j-th quadrature point original_values[j][i] has to be read.
 * @param[out] original_values_dx vector of values of derivatives of basis
 * functions with respect to x in quadrature points. The values are ordered as
 * in original_values.
 * @param[out] original_values_dy vector of values of derivatives of basis
 * functions with respect to y in quadrature points. The values are ordered as
 * in original_values.
 * @param[out] original_values_dz vector of values of derivatives of basis
 * functions with respect to z in quadrature points. The values are ordered as
 * in original_values. In 2D still a vector is needed but it is not touched.
 * @returns for every quadrature point of a given (d-1) dimensional quadrature
 * formula on a facet of a cell the evaluation of the basis functions and their
 * derivatives in each point.
 */
template < int d >
void compute_basis_fct_values_at_joint(const FiniteElement& finite_element,
    const TBaseCell* cell, const int& joint_index, const int& base_vec_dim, const
    TQuadFormula& quad_formula,
    std::vector< std::vector<double> >& original_values,
    std::vector< std::vector<double> >& original_values_dx,
    std::vector< std::vector<double> >& original_values_dy,
    std::vector< std::vector<double> >& original_values_dz,
    const TCollection& coll
    );

template < int d >
void compute_basis_fct_values_at_joint(const FiniteElement& finite_element,
    const TBaseCell* cell, const int& joint_index, const int& base_vec_dim,
    const TQuadFormula& quad_formula,
    std::vector< std::vector<double> >& u,
    std::vector< std::vector<double> >& ux,
    std::vector< std::vector<double> >& uy,
    std::vector< std::vector<double> >& uxx,
    std::vector< std::vector<double> >& uxy,
    std::vector< std::vector<double> >& uyy,
    std::vector< std::vector<double> >& uz,
    const TCollection& coll
    );
/**
 * @brief computes the outer unit normal of a facet of a given cell in given
 * quadrature points.
 * @details This function computes a vector of parmoon::Point where each entry
 * is the outer unit normal of a given facet of a given cell in given quadrature
 * points. The quadrature points are given by an (d-1) dimensional quadrature
 * formula on the reference cell. These points are then transformed onto points
 * in the facet of the original cell and the unit normal is computed in those
 * points.
 * @param[in] cell the cell for which the unit normal of the given facet shall
 * be computed
 * @param[in] joint_index the number of the facet of the cell
 * @param[in] quad_formula A (d-1) dimensional quad formula on the reference
 * cell
 * @return A vector of parmoon::Point that are the d dimensional unit normal
 * vectors of the (d-1) dimensional facet in the quadrature points.
 */
std::vector<parmoon::Point> compute_face_normal( const TBaseCell* cell,
    const int joint_index, const TQuadFormula quad_formula );

/**
 * @brief transform quadrature points of a (d-1) dimensional quadrature rule on
 * a reference geometry to d dimensional points on the facet of a original cell
 * @details Assume one wants to integrate along a facet. To this purpose,
 * a (d-1) dimensional quadrature rule is needed, which is defined on
 * a reference geometry, i.e. a line (2D), a triangle or a quadrilateral (3D).
 * This function computes the d dimensional points that lie in a facet of the
 * original d dimensional cell with respect to this quadrature formula.
 * @param[in] cell original cell
 * @param[in] ref_element reference element of the cell
 * @param[in] ref_trans_id type of reference transformation
 * @param[in] joint_index local index of the facet
 * @param[in] quad_form (d-1) dimensional quadrature formula on the reference
 * cell
 * @returns a vector of d dimensional parmoon::Point that are the coordinates of
 * the quadrature points in the facet of the original cell
 */
std::vector<parmoon::Point> transform_quad_points( const TBaseCell* cell, const
    BFRefElements& ref_element, const ReferenceTransformation_type&
    ref_trans_id, const int& joint_index, const TQuadFormula* quad_form );

/**
 * @brief computes the weights on the original cell for a given quadrature
 * formula that is defined on a (d-1) dimensional reference facet
 * @details Assume a quadrature rule on a (d-1) dimensional reference facet is
 * given. To use this quadrature rule a transformation to the original
 * d dimensional cell is needed. This function computes the weights for the
 * quadrature points on the original cell.
 * @param[in] cell original mesh cell
 * @param[in] joint_index the facet to which the quadrature points shall be
 * transformed
 * @param[in] quad_form (d-1) dimensional quad formula on a reference facet
 * @returns a vector of weights for the quadrature points
 */
std::vector<double> compute_orig_weights( TBaseCell* cell, const int& joint_index,
    const TQuadFormula* quad_form);
#endif // ASSEMBLE_DG


template< int d>
void compute_DG_bilinear_form_local_Oseen(const ParameterDatabase& param_db, const int& base_vec_dim, const
    std::vector<double>* val_bas_fct, const std::vector<double>* val_bas_fct_dx,
    const std::vector<double>* val_bas_fct_dy, const std::vector<double>*
    val_bas_fct_dz, const std::vector<double>* val_bas_fct_neigh, const
    std::vector<double>* val_bas_fct_neigh_dx, const std::vector<double>*
    val_bas_fct_neigh_dy, const std::vector<double>* val_bas_fct_neigh_dz, const
    std::vector<int>& index_map, const double* joint_normal, const double&
    joint_diam, const TQuadFormula* quad_form_joint, const std::vector<double>&
    orig_weights, double** coefficients, const bool& is_boundary_edge, double**
    loc_mat, int matrix_flag);

#ifdef __2D__
void compute_DG_bilinear_form_local_Oseen(const ParameterDatabase& param_db, 
    const int& base_vec_dim, 
    const std::vector<double>* val_bas_fct, const std::vector<double>* val_bas_fct_dx,
    const std::vector<double>* val_bas_fct_dy,
    const std::vector<double>* val_bas_fct_dxx,
    const std::vector<double>* val_bas_fct_dxy,
    const std::vector<double>* val_bas_fct_dyy,
    const std::vector<double>* val_bas_fct_dz,
    const std::vector<double>* val_bas_fct_neigh,
    const std::vector<double>* val_bas_fct_neigh_dx, 
    const std::vector<double>* val_bas_fct_neigh_dy, 
    const std::vector<double>* val_bas_fct_neigh_dxx,
    const std::vector<double>* val_bas_fct_neigh_dxy, 
    const std::vector<double>* val_bas_fct_neigh_dyy,
    const std::vector<double>* val_bas_fct_neigh_dz,
    const std::vector<int>& index_map, const double* joint_normal, const double&
    joint_diam, const TQuadFormula* quad_form_joint, const std::vector<double>&
    orig_weights, double** coefficients, const bool& is_boundary_edge, double**
    loc_mat, int matrix_flag);
#endif


#ifdef __2D__
void compute_CIP_Stab_Oseen(const ParameterDatabase& param_db, const int& base_vec_dim, const
    std::vector<double>* val_bas_fct, 
    const std::vector<double>* val_bas_fct_dx,
    const std::vector<double>* val_bas_fct_dy,
    const std::vector<double>* val_bas_fct_dxx,
    const std::vector<double>* val_bas_fct_dxy,
    const std::vector<double>* val_bas_fct_dyy,
    const std::vector<double>* val_bas_fct_dz,
    const std::vector<double>* val_bas_fct_neigh, 
    const std::vector<double>* val_bas_fct_neigh_dx, 
    const std::vector<double>* val_bas_fct_neigh_dy, 
    const std::vector<double>* val_bas_fct_neigh_dxx,
    const std::vector<double>* val_bas_fct_neigh_dxy, 
    const std::vector<double>* val_bas_fct_neigh_dyy,
    const std::vector<double>* val_bas_fct_neigh_dz, 
    const std::vector<int>& index_map, 
    const double* joint_normal, const double&
    joint_diam, const TQuadFormula* quad_form_joint, const std::vector<double>&
    orig_weights, double** coefficients, 
    double ** parameters,
    const bool& is_boundary_edge, 
    const int* global_dof, const int* global_dof_neigh, const TBaseCell*
    cell, int cell_no, int matrix_flag, double **loc_matA11, 
    double **loc_matA12, double **loc_matA21, double **loc_matA22);

#endif
#ifdef __2D__
void AddJumpStabilizationCIP(int n_fespaces, const TFESpace2D** fespaces, int
    n_sqmatrices, TSquareMatrix2D** sqmatrices, 
    BoundCondFunct2D** BoundaryConditions, 
    BoundValueFunct2D** BoundaryValues, 
    const ParameterDatabase& param_db, 
    LocalAssembling2D& la);
#endif
