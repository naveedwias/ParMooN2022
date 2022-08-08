#include <Assemble_DG.h>
#include <SquareMatrix2D.h>
#include "BaseCell.h"
#include "BoundaryAssembling2D.h"
#include "Database.h"
#include <algorithm>
#include "HexaIsoparametric.h"
#include "InterfaceJoint.h"
#include "QuadIsoparametric.h"
#include "RefTrans3D.h"
#include "SquareMatrix3D.h"
#include "TetraIsoparametric.h"
#include "TriaIsoparametric.h"
#include "templateNames.h"
#include "Point.h"
#include "Line.h"
#include "Triangle.h"
#include "Quadrangle.h"

// added for the sign change
#include <FEFunction2D.h>
#include "QuadratureFormulaDatabase.h"

#include "FEDatabase.h"
#ifdef __2D__
using BoundaryFacet = TBoundEdge;
#else
using BoundaryFacet = TBoundFace;
#endif

ParameterDatabase default_dg_database()
{
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("Discontinuous Galerkin parameter database");

  std::string description = "This specifies whether a symmetric (1 = SIPG), an "
    "incomplete (0 = IIPG) or a non-symmetric (-1 = NIPG) interior penalty "
    "method shall be used in the DG discretization of the Laplace operator.";
  std::set<int> range = {-1, 0, 1};
  db.add<int>("symmetry_DG", 1, description, range);

  description = "The penalty parameter used in the DG discretization of the "
    "Laplacian to penalize jumps can be set using this parameter. Note that on "
    "boundary edges twice this parameter is chosen, while on interior edges the"
    " exact value is used.";
  db.add("face_sigma_DG", 1., description, 0., 1.e100);

  description = "This parameter controls the amount of upwinding "
    "used in the DG discretization of the convective term. The value 0 "
    "corresponds to a central flux, i.e. no upwinding, and 1 leads to the usual"
    " upwind flux used in finite volume methods. Any other value between 0 and "
    "1e100 is possible as well, but not a-priori meaningful.";
  db.add("eta_upwind_DG", 1., description, 0., 1.e100);
  return db;
}

template< int d >
void Assemble_DG(const typename Template_names<d>::CoeffFct& coeff_fct, int
    n_fespaces, const typename Template_names<d>::FESpace** fespaces, int
    n_sqmatrices, typename Template_names<d>::SquareMatrixD** sqmatrices, int
    n_matrices, typename Template_names<d>::MatrixD **/*matrices*/, int n_rhs,
    double** rhs, typename Template_names<d>::BoundaryConditionFunction**
    BoundaryConditions, typename Template_names<d>::BoundaryValuesFunction**
    BoundaryValues, const ParameterDatabase& param_db, 
    const std::string system_type, int matrix_flag, 
    typename Template_names<d>::SquareMatrixD **sqMatA22)
{
  // IDEA: Here the edge or face contributions of the DG method are computed. To
  // perform a loop  over the facets, a loop over cells is carried out in which
  // all the local facets are considered if they haven't been considered before.
  // On each facet the values of the basis functions of the two neighbouring
  // cells have to be computed. Afterwards, the entries for the matrix and the
  // right-hand side can be calculated according to the DG bilinear form and
  // linear form and can be added to the correct position in the matrix / vector
  // of the system.
  //
  // In principle, the code is maybe written such that different polynomial
  // degrees on adjacent cells should not make any problem. However, this is
  // untested since the constructor in the ConvectionDiffusion class only allows
  // for a single degree shared by all cells.

  // Isoparametric elements do not work. If you have time, ParMooN would be
  // happy if you fix this.
  // It is important to check that on boundary facets the values of the basis
  // functions are correct and in addition the length of the facets and the
  // normal is correct. And the integration has to be set up correctly.  Right
  // now, this is not the case.
  if (TDatabase::ParamDB->USE_ISOPARAMETRIC == 1)
  {
    ErrThrow("DG method with isoparametric elements are not implemented yet.");
  }

  if (n_sqmatrices > 1 || n_fespaces > 1 || n_matrices > 0 || n_rhs > 1)
  {
    ErrThrow("Only the scalar case is implemented yet.");
  }

  if ( d < 2 || d > 3 )
  {
    ErrThrow("Only 2D or 3D are available");
  }

  // Create needed parameter for dg bilinear form and update database
  auto db = default_dg_database();
  db.merge(param_db, false); // update this database with given values
  if ( db["eta_upwind_DG"].get<double>() > 1 )
  {
    Output::warn("Values of eta_upwind_DG larger than 1 might not be "
        "reasonable. Be sure what you do.");
  }

  // Prepare loop over cells -> mark cells to access cell numbers
  auto collection = fespaces[0]->GetCollection();
  collection->mark_all_cells(); // set cell_index
  auto n_cells = collection->GetN_Cells();  // number of cells

  // Set up quadrature rule on facets
  // Find out maximal degree needed for integration, i.e. the largest degree
  // over all cells. Note that on single cells this degree can be too large if
  // the degree differs between cell, and hence a computational overhead exists.
  // Fix it, if you think it is a substantial problem.
  int max_degree = 0;
  for (int cell_i = 0; cell_i < n_cells; ++cell_i)
  {
    for(int fe_i = 0; fe_i < n_fespaces; fe_i++)
    {
      auto& element = fespaces[fe_i]->get_fe(cell_i);
      auto degree = element.GetBaseFunct()->GetPolynomialDegree();
      max_degree = (degree > max_degree) ? degree : max_degree;
    }
  }
  // There will be integrands of the type base_fct * base_fct. Hence, the
  // quadrature rule has to be exact for at least 2 * degree.
  // This in principle allows for different polynomial degrees in two
  // adjacent cells, since max_degree is large enough, but at the price that the
  // quadrature rule may be too complicated in possibly many cells.
  BFRefElements joint_type;
  if (d == 2)
  {
    // In 2D facets are lines
    joint_type = BFRefElements::BFUnitLine;
  }
  else
  {
    // Currently in 3D mixed grids are not possible, i.e. the mesh
    // consists of either hexahedra or tetrahedra. Therefore, only one
    // quadrature rule is needed depending on the type of cells. A criterion for
    // distinguishing between tetrahedra and hexahedra is e.g. the number of
    // vertices in the first cell ( 4 for tetrahedra and 6 for hexahedra ).
    if (collection->GetCell(0)->GetN_Vertices() == 4)
    {
      joint_type = BFRefElements::BFUnitTriangle;
    }
    else if (collection->GetCell(0)->GetN_Vertices() == 8)
    {
      joint_type = BFRefElements::BFUnitSquare;
    }
    else
    {
      ErrThrow("Only tetrahedral or hexahedral grids are allowed.");
    }
  }
  auto quad_formula_joint = QuadratureFormulaDatabase::qf_from_degree(
      2 * max_degree, joint_type);
  auto n_quad_pts = quad_formula_joint->GetN_QuadPoints();

  //sqmatrices[0]->PrintFull("A_volume");
  //added to delete the terms from the volume integral in the global matrix
  //sqmatrices[0]->reset();

  // Perform loop over cells to loop over local edges
  for (int cell_i = 0; cell_i < n_cells; ++cell_i)
  { // loop over cells
    //Output::print("cell ", cell_i);
    auto finite_element = fespaces[0]->get_fe(cell_i);
    auto cell = collection->GetCell(cell_i);
    auto n_joints = cell->GetN_Joints();

    cell-> SetNormalOrientation();

    // for all edges
    for (int joint_i = 0; joint_i < n_joints; ++joint_i)
    {
      //Output::print(" joint ", joint_i);
      // Calculate the values of the basis function and their derivatives in
      // quad pts for this edge
      std::vector<std::vector<double>> val_bas_fct;
      std::vector<std::vector<double>> val_bas_fct_dx;
      std::vector<std::vector<double>> val_bas_fct_dy;
      std::vector<std::vector<double>> val_bas_fct_dz;
      std::vector<std::vector<double>> val_bas_fct_uxx;
      std::vector<std::vector<double>> val_bas_fct_uxy;
      std::vector<std::vector<double>> val_bas_fct_uyy;
      // val_bas_fct[j] contains the values of the basis functions on the
      // original cell evaluated at quadrature point j.
      // val_bas_fct_dx (_dy, _dz) are the derivatives wrt. x (y, z)
      // This is filled in the next lines
      auto base_vec_dim = finite_element.GetBaseFunct()->GetBaseVectDim();
      if(system_type == "CD" || system_type == "NSE")
      {
        compute_basis_fct_values_at_joint<d>(finite_element, cell, joint_i,
                    base_vec_dim, *quad_formula_joint,
                    val_bas_fct, val_bas_fct_dx, val_bas_fct_dy, val_bas_fct_dz, *collection);
      }
      else
      {
        compute_basis_fct_values_at_joint<d>(finite_element, cell, joint_i, 
                base_vec_dim, *quad_formula_joint, val_bas_fct, val_bas_fct_dx,
               val_bas_fct_dy, 
               val_bas_fct_uxx, val_bas_fct_uxy, val_bas_fct_uyy,
               val_bas_fct_dz, *collection);
      }

      auto ref_element = finite_element.GetBaseFunct()->GetRefElement();
      auto ref_trans_id = finite_element.GetRefTransID();
      auto transformed_quad_points = transform_quad_points(cell, ref_element,
      ref_trans_id, joint_i, quad_formula_joint);



      // Check for neighbours to distinguish between boundary and inner edges
      auto joint = cell->GetJoint(joint_i);
      auto neigh = joint->GetNeighbour(cell);

      // Not every facet has to be considered. Only Dirichlet boundary facets
      // and interior facet have to be considered. Neumann boundary facets are
      // already treated by previously called assemble routines. Each interior
      // facets only needs to be considered once. This is checked by the
      // following boolean.
      bool consider_joint = false;

      // Storage for values of neighbor basis function evaluated at the quad
      // pts. This is filled only in the case of interior edges.
      std::vector<std::vector<double>> val_bas_fct_neigh;
      std::vector<std::vector<double>> val_bas_fct_neigh_dx;
      std::vector<std::vector<double>> val_bas_fct_neigh_dy;
      std::vector<std::vector<double>> val_bas_fct_neigh_uxx;
      std::vector<std::vector<double>> val_bas_fct_neigh_uxy;
      std::vector<std::vector<double>> val_bas_fct_neigh_uyy;
      std::vector<std::vector<double>> val_bas_fct_neigh_dz;
      std::vector<int> index_map;
      int cell_nr_neigh = -1; // dummy value to remove warning
      if (neigh)
      { // interior edge
        // In MPI case unfortunately the unused cells are not deleted correctly.
        // This results in the fact that cell think they have neighbours which
        // are not in the collection. The only way I found to exclude those
        // cells is to check the clipboard even if I try to avoid this. The
        // clipboard of these "ghost cells" is set to -1 as for cells with
        // hanging nodes. But those cells have to be excluded
        if (neigh->GetClipBoard() >= 0 )
        {
          cell_nr_neigh = neigh->GetCellIndex(); // cell number of neighbor

          // Every interior edge has to be considered only once. We consider it
          // only at their first appearance as local edge, i.e. if cell_i is
          // smaller than cell_nr_neigh. In this case also the normal is defined
          // as the outer unit normal of cell i.
          if (cell_i < cell_nr_neigh)
          { // interior edges was not considered yet
            consider_joint = true;

            // find out the local joint number corresponding to joint_i
            auto facet_nr_neigh = joint->get_joint_nr_in_cell(neigh);

            // val_bas_fct_neigh[j] contains the values of the basis functions
            // on the original cell evaluated at quadrature point j.
            // val_bas_fct_neigh_dx (_y, _z) are the derivatives wrt. x (y, z)
            // This is filled in the next line
            auto fe_neigh = fespaces[0]->get_fe(cell_nr_neigh);
            int base_vec_dim_ne = fe_neigh.GetBaseFunct()->GetBaseVectDim();
            if(system_type == "CD" || system_type == "NSE")
            {
              compute_basis_fct_values_at_joint<d>( fe_neigh, neigh,
                  facet_nr_neigh, base_vec_dim_ne,
                  *quad_formula_joint, val_bas_fct_neigh, val_bas_fct_neigh_dx,
                  val_bas_fct_neigh_dy, val_bas_fct_neigh_dz, *collection);
            }
            else{
              compute_basis_fct_values_at_joint<d>(fe_neigh, neigh, facet_nr_neigh, 
                base_vec_dim_ne, *quad_formula_joint, val_bas_fct_neigh, val_bas_fct_neigh_dx,
               val_bas_fct_neigh_dy, 
               val_bas_fct_neigh_uxx, val_bas_fct_neigh_uxy, val_bas_fct_neigh_uyy,
               val_bas_fct_neigh_dz, *collection);
            }

//             auto ref_element = fe_neigh.GetBaseFunct()->GetRefElement();
//             auto ref_trans_id = fe_neigh.GetRefTransID();
//             auto transformed_quad_points = transform_quad_points(neigh, ref_element,
//             ref_trans_id, facet_nr_neigh, quad_formula_joint);

            // Compute an index map that maps the indices of the quadrature
            // points seen from the cell to the quadrature points seen form the
            // neighboring cell
            index_map = compute_index_map_for_facet_quad_points_of_neigh<d>(
                cell, neigh, joint_i, facet_nr_neigh, finite_element, fe_neigh,
                quad_formula_joint );
          } // endif cell_i < cell_nr_neigh
        }
      } // endif neigh
      else
      { // boundary edge
        // determine type of boundary condition since only Dirichlet edges
        // appear in the bilinear form
        auto condition = get_type_of_boundary_condition<d>(
            (const BoundaryFacet*) joint, BoundaryConditions[0]);
        if (condition != DIRICHLET && condition != NEUMANN)
        {
          ErrThrow("Other conditions than Dirichlet and Neumann boundary ",
              "conditions are not implemented yet.");
        }

        // Only Dirichlet boundary facets have to be considered.
        consider_joint = (condition == DIRICHLET);
      } // endif boundary edge

      if (consider_joint)
      {
        // Compute bilinear form:

        // Compute coefficients of PDE in quadrature points on edge_i
        std::vector<double*> coefficients(n_quad_pts);
        std::vector<double> aux(n_quad_pts * n_local_coeff);
        for(int quad_pt_i = 0; quad_pt_i < n_quad_pts; quad_pt_i++)
        {
          coefficients[quad_pt_i] = aux.data() + quad_pt_i * n_local_coeff;
        }
        compute_coefficients<d>(cell, joint_i, finite_element,
            quad_formula_joint, coeff_fct, coefficients.data());

        // compute geometrical information needed for integration
        auto facet_diam = cell->ComputeDiameterOfJoint(joint_i);
        auto face_normal = cell->ComputeOuterJointNormal(joint_i);

        // allocate local matrix. Upper-left square of size n_dofs^2
        // correspond to dofs from cell.
        // In case of interior edges, lower-right square of size
        // n_dofs_neigh^2 corresponds to dofs of neighbor cell, upper-right
        // and lower-left part of size n_dofs * n_dofs_neigh correspond to
        // coupling between test and ansatz dofs of cell and neighbour and
        // other way around.
        auto n_dofs = finite_element.GetN_DOF();
        int n_dofs_neigh = 0;
        if (neigh)
        {
          auto element_neigh = fespaces[0]->get_fe(cell_nr_neigh);
          n_dofs_neigh = element_neigh.GetN_DOF();
        }
        auto n_dof_total = n_dofs + n_dofs_neigh;
        std::vector<double*> loc_mat(n_dof_total);
        std::vector<double*> loc_matA22(n_dof_total);
        std::vector<double> pointer_vec(n_dof_total * n_dof_total, 0);
        std::vector<double> pointer_vec1(n_dof_total * n_dof_total, 0);
        for (int i = 0; i < n_dof_total; ++i)
        {
          loc_mat[i] = &pointer_vec[i * n_dof_total]; // local matrix
          loc_matA22[i] = &pointer_vec1[i * n_dof_total]; // local matrix
        }

        auto orig_weights = compute_orig_weights(cell, joint_i,
            quad_formula_joint);
        auto base_vec_dim = finite_element.GetBaseFunct()->GetBaseVectDim();

        auto global_dof = fespaces[0]->GetGlobalDOF(cell_i);
        if (neigh)
        {
            if (system_type == "CD")
            {   compute_DG_bilinear_form_local<d>(db, base_vec_dim, val_bas_fct.data(),
                        val_bas_fct_dx.data(), val_bas_fct_dy.data(),
                        val_bas_fct_dz.data(), val_bas_fct_neigh.data(),
                        val_bas_fct_neigh_dx.data(), val_bas_fct_neigh_dy.data(),
                        val_bas_fct_neigh_dz.data(), index_map, face_normal.data(),
                        facet_diam, quad_formula_joint, orig_weights,
                        coefficients.data(), false, loc_mat.data());
            }
            else if (system_type == "NSE")
             {
                auto global_dof_neigh = fespaces[0]->GetGlobalDOF(cell_nr_neigh);
                compute_DG_bilinear_form_local_NSE<d>(db, base_vec_dim, val_bas_fct.data(),
                        val_bas_fct_dx.data(), val_bas_fct_dy.data(),
                        val_bas_fct_dz.data(), val_bas_fct_neigh.data(),
                        val_bas_fct_neigh_dx.data(), val_bas_fct_neigh_dy.data(),
                        val_bas_fct_neigh_dz.data(), index_map, face_normal.data(),
                        facet_diam, quad_formula_joint, orig_weights,
                        coefficients.data(), false, loc_mat.data(), global_dof, global_dof_neigh, cell, cell_i, joint_i);
            }
            else if(system_type == "Oseen")
            {
#ifdef __2D__
//               compute_DG_bilinear_form_local_Oseen(db, base_vec_dim, val_bas_fct.data(),
//                         val_bas_fct_dx.data(), val_bas_fct_dy.data(),
//                         val_bas_fct_uxx.data(), val_bas_fct_uxy.data(),
//                         val_bas_fct_uyy.data(),
//                         val_bas_fct_dz.data(), val_bas_fct_neigh.data(),
//                         val_bas_fct_neigh_dx.data(), val_bas_fct_neigh_dy.data(),
//                         val_bas_fct_neigh_uxx.data(), val_bas_fct_neigh_uxy.data(),
//                         val_bas_fct_neigh_uyy.data(),
//                         val_bas_fct_neigh_dz.data(), index_map, face_normal.data(),
//                         facet_diam, quad_formula_joint, orig_weights,
//                         coefficients.data(), false, loc_mat.data(), matrix_flag);
              
              auto global_dof_neigh = fespaces[0]->GetGlobalDOF(cell_nr_neigh);
                compute_CIP_Stab_Oseen(db, base_vec_dim, 
                        val_bas_fct.data(),
                        val_bas_fct_dx.data(), 
                        val_bas_fct_dy.data(),
                        val_bas_fct_uxx.data(), 
                        val_bas_fct_uxy.data(),
                        val_bas_fct_uyy.data(),
                        val_bas_fct_dz.data(), 
                        val_bas_fct_neigh.data(),
                        val_bas_fct_neigh_dx.data(), 
                        val_bas_fct_neigh_dy.data(),
                        val_bas_fct_neigh_uxx.data(),
                        val_bas_fct_neigh_uxy.data(),
                        val_bas_fct_neigh_uyy.data(),
                        val_bas_fct_neigh_dz.data(), 
                        index_map, face_normal.data(),
                        facet_diam, quad_formula_joint, orig_weights,
                        coefficients.data(), 
                        nullptr,
                        false, 
                        global_dof, 
                        global_dof_neigh, cell, cell_i, matrix_flag, 
                        loc_mat.data(), nullptr, nullptr, loc_matA22.data());
#endif
            }
        }
        else
        {

            if (system_type == "CD")
            {
                compute_DG_bilinear_form_local<d>(db, base_vec_dim, val_bas_fct.data(),
                              val_bas_fct_dx.data(), val_bas_fct_dy.data(),
                              val_bas_fct_dz.data(), face_normal.data(), facet_diam,
                              quad_formula_joint, orig_weights, coefficients.data(),
                              loc_mat.data());
            }
            else if (system_type == "NSE")
             {
                compute_DG_bilinear_form_local_NSE<d>(db, base_vec_dim, val_bas_fct.data(),
                              val_bas_fct_dx.data(), val_bas_fct_dy.data(),
                              val_bas_fct_dz.data(), face_normal.data(), facet_diam,
                              quad_formula_joint, orig_weights, coefficients.data(),
                              loc_mat.data(), global_dof);

             }
        }

        // compute local right-hand side
        std::vector<double> loc_rhs(n_dofs, 0);
        if (!neigh)
        {
          // compute local right-hand side

         if (system_type == "CD")
            {
                compute_DG_rhs_local<d>(db, cell, joint_i, finite_element,
                              val_bas_fct.data(), val_bas_fct_dx.data(), val_bas_fct_dy.data(),
                              val_bas_fct_dz.data(), face_normal.data(), facet_diam,
                              quad_formula_joint, orig_weights, coefficients.data(),
                              BoundaryValues[0], loc_rhs.data());
            }
            else if (system_type == "NSE")
             {
                compute_DG_rhs_local_NSE<d>(db, cell, joint_i, finite_element,
                              val_bas_fct.data(), val_bas_fct_dx.data(), val_bas_fct_dy.data(),
                              val_bas_fct_dz.data(), face_normal.data(), facet_diam,
                              quad_formula_joint, orig_weights, coefficients.data(),
                              BoundaryValues[0], loc_rhs.data());
             }
        }

        // add local to global
        if (neigh && neigh->GetClipBoard() >= 0)
        {
          auto global_dof_neigh = fespaces[0]->GetGlobalDOF(cell_nr_neigh);
          copy_local_to_global<d>(n_dofs, n_dofs_neigh, global_dof,
              global_dof_neigh, loc_mat.data(), nullptr, sqmatrices[0],
              nullptr, loc_mat.data(), sqMatA22[0]);
        }
        else
        {
          if(system_type != "Oseen")
            copy_local_to_global<d>(n_dofs, global_dof, loc_mat.data(),
                                    loc_rhs.data(), sqmatrices[0], rhs[0]);
        }
      }
    } // endfor loop over edges
  } // endfor loop over cells
  //sqmatrices[0]->PrintFull("A_DG");
} // end of Assemble2D_DG
#ifdef __3D__
template void Assemble_DG<3>(const typename Template_names<3>::CoeffFct&
    Coeff, int n_fespaces, const typename Template_names<3>::FESpace**
    fespaces, int n_sqmatrices, typename Template_names<3>::SquareMatrixD**
    sqmatrices, int n_matrices, typename Template_names<3>::MatrixD
    **/*matrices*/, int n_rhs, double** rhs, typename
    Template_names<3>::BoundaryConditionFunction** BoundaryConditions, typename
    Template_names<3>::BoundaryValuesFunction** BoundaryValues, const
    ParameterDatabase& param_db, const std::string system_type, int matrix_flag, 
    typename Template_names<3>::SquareMatrixD **sqMatA22);
#else
template void Assemble_DG<2>(const typename Template_names<2>::CoeffFct&
    Coeff, int n_fespaces, const typename Template_names<2>::FESpace**
    fespaces, int n_sqmatrices, typename Template_names<2>::SquareMatrixD**
    sqmatrices, int n_matrices, typename Template_names<2>::MatrixD
    **/*matrices*/, int n_rhs, double** rhs, typename
    Template_names<2>::BoundaryConditionFunction** BoundaryConditions, typename
    Template_names<2>::BoundaryValuesFunction** BoundaryValues, const
    ParameterDatabase& param_db, const std::string system_type, int matrix_flag, 
    typename Template_names<2>::SquareMatrixD **sqMatA22);
#endif

// Check whether quadrature points on an edge are in the same order seen from
// the two adjacent cells
template <int d >
bool is_correct_order_of_quad_points( const TBaseCell* cell, const TBaseCell*
    cell_neigh, const int& joint_nr, const int& joint_nr_neigh)
{
  // Check input of function
  if (cell->GetJoint(joint_nr) != cell_neigh->GetJoint(joint_nr_neigh))
  {
    ErrThrow("The joints given by ", joint_nr, " (cell) and ", joint_nr_neigh,
        " (neighbor cell) are not the same. Either one of the cells or one of ",
        "the edge number is wrong.");
  }
  // Here the actual computation takes place
  bool is_correct_order;
  if (d == 2)
  {
    // In 2D the quadrature points on an edge can either be in the correct order
    // or exactly the other way around. This comes from the fact that the edge
    // on each cell is defined by a start and an end vertex and either the start
    // vertices from both cells (and the end vertices) are the same or they are
    // exactly the other way around
    // Therefore, we compute the coordinates of the start and end point of the
    // edge and compare if the start points are the same.
    double x0, y0, z0;
    cell->GetVertex(joint_nr)->GetCoords(x0, y0, z0);
    double x0_neigh, y0_neigh;
    cell_neigh->GetVertex(joint_nr_neigh)->GetCoords(x0_neigh, y0_neigh, z0);

    is_correct_order = (x0_neigh == x0 && y0_neigh == y0);
  }
  else if (d == 3)
  {
    // In 3D the vertices defining a face can be ordered differently. They have
    // exactly the same order if the first three vertices of a face are the
    // same. For triangles this is trivially true and for polygons with more
    // than three vertices this is true since the neighboring relationship
    // between the vertices is fixed, i.e. the numbering of the vertices can
    // either be mirrored or turned or a combination of both. But two
    // neighboring vertices in a facet are neighbors no matter if you look from
    // the cell or the neighboring cell. Using turning and mirroring it is
    // possible to let two vertices coincide and the rest not, but if three
    // consecutive vertices coincide then also the rest has to be equal.
    // Therefore, the first three vertices of the facet have to be determined
    // via the face_vertex_map, see TShapeDesc.
    const int* face_vertex_map;
    const int* face_vertex_map_length;
    int max_n_vertices_per_face;
    cell->GetShapeDesc()->GetFaceVertex(face_vertex_map, face_vertex_map_length,
        max_n_vertices_per_face);
    const int* face_vertex_map_neigh;
    const int* face_vertex_map_length_neigh;
    int max_n_vertices_per_face_neigh;
    cell_neigh->GetShapeDesc()->GetFaceVertex(face_vertex_map_neigh,
        face_vertex_map_length_neigh, max_n_vertices_per_face_neigh);

    // Check consistency of calculated maps
    if ( face_vertex_map_length[joint_nr] !=
        face_vertex_map_length_neigh[joint_nr_neigh] )
    {
      auto cell_nr = cell->GetCellIndex();
      auto cell_nr_neigh = cell_neigh->GetCellIndex();
      ErrThrow("Number of vertices of face ", joint_nr, " of cell ", cell_nr,
          " is not equal to number of vertices of face ", joint_nr_neigh, " of",
          " neighboring cell ", cell_nr_neigh, ".");
    }
    if (face_vertex_map_length[joint_nr] < 3)
    {
      ErrThrow("A face with positive area in 3D has to have at least ",
          "3 coordinates.");
    }

    // Find the index in face_vertex_map where the joint_nr-th face starts. Same
    // for neighbor. To do so we add the length all faces prior to joint_nr have
    int start_index = 0;
    for (int face_i = 0; face_i < joint_nr; ++face_i)
    {
      start_index += face_vertex_map_length[face_i];
    }
    int start_index_neigh = 0;
    for (int face_i = 0; face_i < joint_nr_neigh; ++face_i)
    {
      start_index_neigh += face_vertex_map_length_neigh[face_i];
    }

    // Determine coordinates of first three vertices
    double x0, y0, z0;
    double x1, y1, z1;
    double x2, y2, z2;
    double x0_neigh, y0_neigh, z0_neigh;
    double x1_neigh, y1_neigh, z1_neigh;
    double x2_neigh, y2_neigh, z2_neigh;

    cell->GetVertex(face_vertex_map[start_index + 0])->GetCoords(x0, y0, z0);
    cell->GetVertex(face_vertex_map[start_index + 1])->GetCoords(x1, y1, z1);
    cell->GetVertex(face_vertex_map[start_index + 2])->GetCoords(x2, y2, z2);
    cell_neigh->GetVertex(face_vertex_map_neigh[start_index_neigh + 0])
      ->GetCoords(x0_neigh, y0_neigh, z0_neigh);
    cell_neigh->GetVertex(face_vertex_map_neigh[start_index_neigh + 1])
      ->GetCoords(x1_neigh, y1_neigh, z1_neigh);
    cell_neigh->GetVertex(face_vertex_map_neigh[start_index_neigh + 2])
      ->GetCoords(x2_neigh, y2_neigh, z2_neigh);

    // If three consecutive vertices coincide all vertices coincide
    is_correct_order = (
        x0_neigh == x0 && y0_neigh == y0 && z0_neigh == z0 &&
        x1_neigh == x1 && y1_neigh == y1 && z1_neigh == z1 &&
        x2_neigh == x2 && y2_neigh == y2 && z2_neigh == z2
        );
  }
  else
  {
    ErrThrow("Only 2D and 3D are possible.");
  }
  return is_correct_order;
}
#ifdef __3D__
template bool is_correct_order_of_quad_points<3>( const TBaseCell* cell, const
    TBaseCell* cell_neigh, const int& joint_nr, const int& joint_nr_neigh);
#else
template bool is_correct_order_of_quad_points<2>( const TBaseCell* cell, const
    TBaseCell* cell_neigh, const int& joint_nr, const int& joint_nr_neigh);
#endif


template <int d >
std::vector<int> compute_index_map_for_facet_quad_points_of_neigh( const
    TBaseCell* cell, const TBaseCell* cell_neigh, const int& joint_nr, const
    int& joint_nr_neigh, const FiniteElement& fin_ele, const FiniteElement&
    fin_ele_neigh, const TQuadFormula* quad_formula_facet )
{
  if ( !(d == 2 || d==3) )
  {
    ErrThrow("Wrong dimension. Only 2D and 3D are implemented yet.");
  }
  else
  {
    // Set up an index list that maps
    // {0,1,..., n_quad_pts-1} to itself where n_quad_pts is
    // the number of quad points the quadrature formula has.
    // First of all the map is created such that the quadrature points are in
    // the same order, i.e. j |-> j for j \in {0,1,..., n_quad_pts-1}.
    auto n_quad_pts = quad_formula_facet->GetN_QuadPoints();
    std::vector<int> index_map(n_quad_pts);
    for (int index_i = 0; index_i < n_quad_pts; ++index_i)
    {
      index_map[index_i] = index_i;
    }

    // Check if the quadrature points are in the same order and if not compute
    // the map of corresponding indices
    if (!is_correct_order_of_quad_points<d>(cell, cell_neigh, joint_nr,
          joint_nr_neigh))
    {
      if (d == 2)
      {
        // in 2D a (d-1) dimensional quadrature rule is a rule for a line. This
        // line is defined by two vertices, i.e. the starting and the end point
        // of the line. The starting point form the line has in the cell is
        // either the same starting point in the neighboring cell or the end
        // point. Therefore, only two configurations of orderings exist. Here we
        // are in the case that the starting point in the neighboring cell is
        // the end point in the original cell and therefore the quadrature
        // points have to be sorted in reverse order.
        std::reverse(index_map.begin(), index_map.end());
      }
      else
      {
        // in 3D the situation is a little bit more involved. I don't know if
        // also a rule exists how the quadrature points are sorted with respect
        // to the vertices of the face. Therefore, the quadrature points are
        // sorted via searching for the same quadrature point in the list of
        // neighboring quadrature points and identifying the index.

        // We need the three dimensional quadrature points for the facet to
        // identify corresponding points. TODO DO WE REALLY NEED THIS? AREN'T 2D
        // POINTS ENOUGH?
        // Hence the 3D quadrature points in the facet are computed in the next
        // line using the respective reference element.
        auto ref_element = fin_ele.GetBaseFunct()->GetRefElement();
        auto ref_trans_id = fin_ele.GetRefTransID();
        auto ref_element_neigh = fin_ele_neigh.GetBaseFunct()->GetRefElement();
        auto ref_trans_id_neigh = fin_ele_neigh.GetRefTransID();

        auto transformed_quad_points = transform_quad_points(cell, ref_element,
            ref_trans_id, joint_nr, quad_formula_facet );
        auto transformed_quad_points_neigh = transform_quad_points(cell_neigh,
            ref_element_neigh, ref_trans_id_neigh, joint_nr_neigh,
            quad_formula_facet );

        // Now we have two lists of quadrature points and we can identify
        // corresponding points. Since the map is an injective map any
        // quadrature point that was already found can be ignored for the
        // upcoming searches. To archive this a list of possible indices is
        // created and every index that was used is deleted in this list. This
        // might be favourable for many quadrature points.
        auto list_of_possible_indices = index_map;

        // For every quadrature point the corresponding quadrature point is
        // searched for.
        for (int pt_i = 0; pt_i < n_quad_pts; ++pt_i)
        {
          int current_index = -1; // This is the index of the quadrature point
          // search all possible indices
          for ( unsigned int index_i = 0;
              index_i < list_of_possible_indices.size(); ++index_i )
          {
            auto index_candidate = list_of_possible_indices[index_i];
            // Check quad_pt[pt_i] == quad_pt_neigh[index_candidate] and if
            // so set the index and erase the entry in the list of possible
            // indices
            double tolerance = 1e-10;
            if (
                (std::abs(transformed_quad_points[pt_i].x
                          - transformed_quad_points_neigh[index_candidate].x)
                 < tolerance)
                && (std::abs(transformed_quad_points[pt_i].y
                    - transformed_quad_points_neigh[index_candidate].y)
                  < tolerance)
                && (std::abs(transformed_quad_points[pt_i].z
                    - transformed_quad_points_neigh[index_candidate].z)
                  < tolerance)
               )
            {
              current_index = index_candidate;
              list_of_possible_indices.erase(list_of_possible_indices.begin()
                  + index_i);
              break;
            }
          }

          // Check if an index was found
          if (current_index < 0)
          {
            ErrThrow("The quadrature points do not coincide.");
          }
          // Store the index
          index_map[pt_i] = current_index;
        }
      }
    }
    return index_map;
  }
}
#ifdef __3D__
template std::vector<int> compute_index_map_for_facet_quad_points_of_neigh<3>(
    const TBaseCell* cell, const TBaseCell* cell_neigh, const int& facet_nr,
    const int& facet_nr_neigh, const FiniteElement& fin_ele, const
    FiniteElement& fin_ele_neigh, const TQuadFormula* quad_formula_facet );
#else
template std::vector<int> compute_index_map_for_facet_quad_points_of_neigh<2>(
    const TBaseCell* cell, const TBaseCell* cell_neigh, const int& facet_nr,
    const int& facet_nr_neigh, const FiniteElement& fin_ele, const
    FiniteElement& fin_ele_neigh, const TQuadFormula* quad_formula_facet );
#endif


// Compute the coefficients of PDE: Find the coordinates of the  quadrature
// points on original cell and store in x, y, and evaluate the coefficient
// function at those points
template< int d>
void compute_coefficients( const TBaseCell* cell, const int& joint_nr, const
    FiniteElement& finite_element, const TQuadFormula* quad_form, const typename
    Template_names<d>::CoeffFct& coeff_fct, double** coefficients)
{
  // To compute the coefficients at the quadrature points the coeff_fct has to
  // be evaluated at those points.
  if(coeff_fct)
  {
    // Prepare transformation to original cell
    auto ref_element = finite_element.GetBaseFunct()->GetRefElement();
    auto ref_trans_id = finite_element.GetRefTransID();

    // Compute quadrature points on original cell
    auto transformed_quad_points = transform_quad_points(cell, ref_element,
        ref_trans_id, joint_nr, quad_form );

    auto n_quad_pts = quad_form->GetN_QuadPoints();
    std::vector<double> x(n_quad_pts);
    std::vector<double> y(n_quad_pts);

#ifdef __2D__
    for (int pt_i = 0; pt_i < n_quad_pts; ++pt_i)
    {
      x[pt_i] = transformed_quad_points[pt_i].x;
      y[pt_i] = transformed_quad_points[pt_i].y;
    }

    // Evaluate the coefficient function at those points
    coeff_fct(n_quad_pts, x.data(), y.data(), nullptr, coefficients);
#else
    std::vector<double> z(n_quad_pts);
    for (int pt_i = 0; pt_i < n_quad_pts; ++pt_i)
    {
      x[pt_i] = transformed_quad_points[pt_i].x;
      y[pt_i] = transformed_quad_points[pt_i].y;
      z[pt_i] = transformed_quad_points[pt_i].z;
    }
    // Evaluate the coefficient function at those points
    coeff_fct(n_quad_pts, x.data(), y.data(), z.data(), nullptr, coefficients);
#endif
  }
  else
  {
    ErrThrow("Please plug in a coefficient function.");
  }
}
#ifdef __3D__
template void compute_coefficients<3>( const TBaseCell* cell, const int&
    joint_nr, const FiniteElement& finite_element, const TQuadFormula*
    quad_form, const typename Template_names<3>::CoeffFct& coeff_fct, double**
    coefficients);
#else
template void compute_coefficients<2>( const TBaseCell* cell, const int&
    joint_nr, const FiniteElement& finite_element, const TQuadFormula*
    quad_form, const typename Template_names<2>::CoeffFct& coeff_fct, double**
    coefficients);
#endif


// Function to compute the type of a boundary condition
  template< int d>
BoundCond get_type_of_boundary_condition(const BoundaryFacet* joint,
    typename Template_names<d>::BoundaryConditionFunction* boundary_condition)
{
  auto boundary_comp = joint->GetBoundComp();
#ifdef __2D__
  // In 2D the joint aka face is parametrised using two parameters t0 and t1. We
  // therefore look for the parameters and then plug it in into
  // boundary_condition() that returns for a given point the boundary condition

  // Compute parameters of the joint
  double t0, t1;
  joint->GetParameters(t0, t1);
  // get id of the boundary component
  auto comp = boundary_comp->GetID();
  // get type of the boundary condition at the beginning
  // and at the end of the current edge
  BoundCond cond_0, cond_1;
  if (t0 < t1)
  {
    double eps = 1e-12;
    boundary_condition(comp, t0+eps, cond_0);
    boundary_condition(comp, t1-eps, cond_1);
  }
  else
  {
    double eps = 1e-12;
    boundary_condition(comp, t0-eps, cond_0);
    boundary_condition(comp, t1+eps, cond_1);
  }
  if(cond_0 != cond_1)
  {
    ErrThrow("Different boundary conditions on a single edge. This is not ",
        "implemented yet.");
  }
#else
  // The information we need comes from boundary_condition(...) that gives the
  // type of boundary condition at a given point. In general, this point can be
  // any point inside the facet, and in particular the code takes the barycentre
  // of the facet. The barycentre can be computed from the vertices of the
  // facet.  Therefore, we first need to find out the vertices that are a member
  // of this face.

  // Get a neighboring cell and find the number the facet has in the cell
  auto cell = joint->GetNeighb(0);
  auto joint_i = joint->get_joint_nr_in_cell(cell);

  // Get the vertices of the facet using the FaceVertex array from the shape
  // descriptor
  auto shape_desc = cell->GetShapeDesc();
  const int* face_vertex;
  const int* face_vertex_len;
  int max_face_vertex_len;
  shape_desc->GetFaceVertex(face_vertex, face_vertex_len,
      max_face_vertex_len);
  // The barycentre can be computed using the information about the vertices.
  // Therefore the coordinates of the vertices are stored now.

  std::vector<std::array<double, 3>> vertices_vec( face_vertex_len[joint_i] );

  // Find the position where the vertices of the joint_i-th face start. The
  // ordering in FaceVertex is such that first the vertices of the 0th face
  // are specified, then of the 1st face and so on. Therefore, the starting
  // point for our current face is the sum of all vertices that belong to
  // faces with a smaller local face number.
  int starting_index = 0;
  for (int facet_nr = 0; facet_nr < joint_i; ++facet_nr)
  {
    starting_index += face_vertex_len[facet_nr];
  }
  for (int vertex_i = 0; vertex_i < face_vertex_len[joint_i]; ++vertex_i)
  {
    cell->GetVertex( face_vertex[ starting_index + vertex_i]
        )->GetCoords(vertices_vec[vertex_i][0], vertices_vec[vertex_i][1],
          vertices_vec[vertex_i][2]);
  }


  // The barycentre is given by the sum of all vertices divided by the
  // number of vertices.
  std::array<double, 3> barycenter {0, 0, 0};
  for (int dim_i = 0; dim_i < 3; ++dim_i)
  {
    for (int vertex_i = 0; vertex_i < face_vertex_len[joint_i]; ++vertex_i)
    {
      barycenter[dim_i] += vertices_vec[vertex_i][dim_i];
    }
    barycenter[dim_i] /= face_vertex_len[joint_i];
  }

  // Determine the boundary condition in the barycentre of the face
  auto comp = boundary_comp->GetID();
  BoundCond cond_0;
  boundary_condition(comp, barycenter[0], barycenter[1], barycenter[2], cond_0);
#endif
  return cond_0;
} // end get_type_of_boundary_condition
#ifdef __3D__
template BoundCond get_type_of_boundary_condition<3>(const BoundaryFacet* joint,
    typename Template_names<3>::BoundaryConditionFunction* boundary_condition);
#else
template BoundCond get_type_of_boundary_condition<2>(const BoundaryFacet* joint,
    typename Template_names<2>::BoundaryConditionFunction* boundary_condition);
#endif


// Computes the local DG bilinear form for convection-diffusion problems
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
    loc_mat)
{

  // get some information needed for integration
  auto n_quad_pts = quad_form_joint->GetN_QuadPoints();
  auto n_dofs = val_bas_fct[0].size()/base_vec_dim;
  auto n_dofs_neigh = (is_boundary_edge) ? 0 : val_bas_fct_neigh[0].size()/base_vec_dim;
  auto n_dofs_total = n_dofs + n_dofs_neigh;
  auto nx = joint_normal[0];
  auto ny = joint_normal[1];
  auto nz = (d == 3) ? joint_normal[2] : 0.0;

  // get parameters of DG method
  auto face_sigma = param_db["face_sigma_DG"].get<double>(); // edge penalty
  auto kappa = param_db["symmetry_DG"].get<int>();  // symmetry
  auto eta = param_db["eta_upwind_DG"].get<double>(); // upwind

  if (is_boundary_edge)
  { // adapt boundary penalty, see e.g. Kanschat 2008
    face_sigma *= 2;
  }

  for (int quad_pt_i = 0; quad_pt_i < n_quad_pts; ++quad_pt_i)
  {
    // parameter of PDE
    auto diffusion = coefficients[quad_pt_i][0];
    auto beta_1 = coefficients[quad_pt_i][1];
    auto beta_2 = coefficients[quad_pt_i][2];
    auto convection_3 = (d == 3) ? coefficients[quad_pt_i][3] : 0.0;

    for (unsigned int t_dof = 0; t_dof < n_dofs_total; ++t_dof)
    { // for all test dofs

      // Define jump and averages for test functions index [0] corresponds to
      // the basis functions itself index [1/2/3] corresponds to the derivative
      // wrt. x/y/z for the basis functions.
      // No matter if we are in 2D or 3D we always assume 4 values, i.e. the
      // 4 indices mentioned above, and in 2D we set every entry corresponding
      // to z to 0.
      std::vector<double> t_jump(4 * base_vec_dim);
      std::vector<double> t_average(4 * base_vec_dim);

      // distinguish whether test function are related to the cell or to the
      // neighbor. Compute jump and average  accordingly using the fact that the
      // support of a basis function is exactly the cell it lives on.
      if (t_dof < n_dofs)
      { // test dof from cell
        t_jump[0] = val_bas_fct [quad_pt_i][t_dof];
        t_jump[1] = val_bas_fct_dx [quad_pt_i][t_dof];
        t_jump[2] = val_bas_fct_dy [quad_pt_i][t_dof];
        t_jump[3] = (d==3) ? val_bas_fct_dz [quad_pt_i][t_dof] : 0;
        t_average[0] = val_bas_fct [quad_pt_i][t_dof] / 2;
        t_average[1] = val_bas_fct_dx [quad_pt_i][t_dof] / 2;
        t_average[2] = val_bas_fct_dy [quad_pt_i][t_dof] / 2;
        t_average[3] = (d == 3) ? val_bas_fct_dz [quad_pt_i][t_dof] / 2 : 0;

      }
      else
      { // test dof from neighbor
        t_jump[0] = -val_bas_fct_neigh[index_map[quad_pt_i]][t_dof - n_dofs];
        t_jump[1] = -val_bas_fct_neigh_dx[index_map[quad_pt_i]][t_dof - n_dofs];
        t_jump[2] = -val_bas_fct_neigh_dy[index_map[quad_pt_i]][t_dof - n_dofs];
        t_jump[3] = (d == 3) ? -val_bas_fct_neigh_dz[index_map[quad_pt_i]][t_dof - n_dofs] : 0;
        t_average[0] = val_bas_fct_neigh[index_map[quad_pt_i]][t_dof - n_dofs] / 2;
        t_average[1] = val_bas_fct_neigh_dx[index_map[quad_pt_i]][t_dof - n_dofs] / 2;
        t_average[2] = val_bas_fct_neigh_dy[index_map[quad_pt_i]][t_dof - n_dofs] / 2;
        t_average[3] = (d == 3) ? val_bas_fct_neigh_dz[index_map[quad_pt_i]][t_dof - n_dofs] / 2 : 0;
      }
      if (is_boundary_edge)
      { // Correct values for boundary edges: jump and average on the boundary
        // are defined as values seen from cell
        for (int i = 0; i < 4; ++i)
        {
          t_average[i] *= 2;
        }
      }
      for (unsigned int a_dof = 0; a_dof < n_dofs_total; ++a_dof)
      { // for all ansatz dofs
        // Define jump and averages for ansatz functions index [0] corresponds
        // to the basis functions itself index [1/2/3] corresponds to the
        // derivative wrt. x/y/z for the basis functions.  No matter if we are
        // in 2D or 3D we always assume 4 values, i.e. the
        // 4 indices mentioned above, and in 2D we set every entry corresponding
        // to z to 0.
        std::vector<double> a_jump(4 * base_vec_dim);
        std::vector<double> a_average(4 * base_vec_dim);

        // distinguish whether ansatz function are related to the cell or to the
        // neighbor. Compute jump and average  accordingly using the fact that
        // the support of a basis function is exactly the cell it lives on.
        if (a_dof < n_dofs)
        { // ansatz dof from cell
          a_jump[0] = val_bas_fct [quad_pt_i][a_dof];
          a_jump[1] = val_bas_fct_dx [quad_pt_i][a_dof];
          a_jump[2] = val_bas_fct_dy [quad_pt_i][a_dof];
          a_jump[3] = (d == 3) ? val_bas_fct_dz [quad_pt_i][a_dof] : 0;
          a_average[0] = val_bas_fct [quad_pt_i][a_dof] / 2;
          a_average[1] = val_bas_fct_dx [quad_pt_i][a_dof] / 2;
          a_average[2] = val_bas_fct_dy [quad_pt_i][a_dof] / 2;
          a_average[3] = (d == 3) ? val_bas_fct_dz [quad_pt_i][a_dof] / 2 : 0;
        }
        else
        { // ansatz dof from neighbor
          a_jump[0] = -val_bas_fct_neigh [index_map[quad_pt_i]][a_dof - n_dofs];
          a_jump[1] = -val_bas_fct_neigh_dx [index_map[quad_pt_i]][a_dof - n_dofs];
          a_jump[2] = -val_bas_fct_neigh_dy [index_map[quad_pt_i]][a_dof - n_dofs];
          a_jump[3] = (d == 3) ? -val_bas_fct_neigh_dz[index_map[quad_pt_i]][a_dof - n_dofs] : 0;
          a_average[0] = val_bas_fct_neigh [index_map[quad_pt_i]][a_dof - n_dofs] / 2;
          a_average[1] = val_bas_fct_neigh_dx [index_map[quad_pt_i]][a_dof - n_dofs] / 2;
          a_average[2] = val_bas_fct_neigh_dy [index_map[quad_pt_i]][a_dof - n_dofs] / 2;
          a_average[3] = (d == 3) ? val_bas_fct_neigh_dz[index_map[quad_pt_i]][a_dof - n_dofs] / 2 : 0;
        }
        if (is_boundary_edge)
        { // Correct values for boundary edges: jump and average on the boundary
          // are defined as values seen from cell
          for (int i = 0; i < 4; ++i)
          {
            a_average[i] *= 2;
          }
        }

        double integrand = 0;
        // contribution of diffusive term
        integrand -= diffusion * ((a_average[1] * nx + a_average[2] * ny
              + a_average[3] * nz) * t_jump[0]
            // symmetry term
            + kappa * (t_average[1] * nx + t_average[2] * ny + t_average[3]
              * nz) * a_jump[0]);

        // penalty term to achieve coercivity
        integrand += face_sigma / joint_diam * a_jump[0] * t_jump[0];

        // contribution of convective term
        // This is exactly the upwind part of the bilinear form (see e.g.
        // - Kanschat: "Discontinuous Galerkin Methods for Viscous
        //   Incompressible Flow" (2008), chapter 1.4 or
        // - Houston, Schwab, Süli: SIAM J. Numer. Anal., 39(6), 2133–2163,
        //   doi = 10.1137/S0036142900374111)
        // even though it uses the point of view of stabilisation (see e.g.
        // - Brezzi, Marini, Süli: Mathematical Models and Methods in Applied
        //   Sciences, 14(12), 1893-1903, doi = 10.1142/S0218202504003866, or
        // - Di Pietro, Ern: "Mathematical Aspects of Discontinuous Galerkin
        //   Methods", doi = 10.1007/978-3-642-22980-0, chapter 2.3)
        auto inflow = beta_1 * nx + beta_2 * ny + convection_3 * nz;
        if (is_boundary_edge)
        {
          if (inflow < 0)
          {
            integrand -= inflow * a_jump[0] * t_average[0];
          }
        }
        else
        {
          integrand += -inflow * a_jump[0] * t_average[0]
            + eta * 0.5 * std::abs(inflow) * a_jump[0] * t_jump[0];
        }

        auto weight = orig_weights[quad_pt_i];
        // Update local matrix
        loc_mat[t_dof][a_dof] += weight * integrand;
      } // endfor ansatz dofs
    } // endfor test dofs
  } //endfor quad points
} // end compute_DG_bilinear_form_local
#ifdef __3D__
template void compute_DG_bilinear_form_local<3>(const ParameterDatabase&
    param_db, const int& base_vec_dim, const std::vector<double>* val_bas_fct, const std::vector<double>*
    val_bas_fct_dx, const std::vector<double>* val_bas_fct_dy, const
    std::vector<double>* val_bas_fct_dz, const std::vector<double>*
    val_bas_fct_neigh, const std::vector<double>* val_bas_fct_neigh_dx, const
    std::vector<double>* val_bas_fct_neigh_dy, const std::vector<double>*
    val_bas_fct_neigh_dz, const std::vector<int>& index_map, const double*
    joint_normal, const double& joint_diam, const TQuadFormula* quad_form_joint,
    const std::vector<double>& orig_weights, double** coefficients, const bool&
    is_boundary_edge, double** loc_mat);
#else
template void compute_DG_bilinear_form_local<2>(const ParameterDatabase&
    param_db, const int& base_vec_dim, const std::vector<double>* val_bas_fct, const std::vector<double>*
    val_bas_fct_dx, const std::vector<double>* val_bas_fct_dy, const
    std::vector<double>* val_bas_fct_dz, const std::vector<double>*
    val_bas_fct_neigh, const std::vector<double>* val_bas_fct_neigh_dx, const
    std::vector<double>* val_bas_fct_neigh_dy, const std::vector<double>*
    val_bas_fct_neigh_dz, const std::vector<int>& index_map, const double*
    joint_normal, const double& joint_diam, const TQuadFormula* quad_form_joint,
    const std::vector<double>& orig_weights, double** coefficients, const bool&
    is_boundary_edge, double** loc_mat);
#endif

template< int d>
void compute_DG_bilinear_form_local(const ParameterDatabase& param_db, const int& base_vec_dim, const
    std::vector<double>* val_bas_fct, const std::vector<double>* val_bas_fct_dx,
    const std::vector<double>* val_bas_fct_dy, const std::vector<double>*
    val_bas_fct_dz, const double* normal, const double& joint_diam, const
    TQuadFormula* quad_form_joint, const std::vector<double>& orig_weights,
    double** coefficients, double** loc_mat)
{
  std::vector<int> dummy_index_map;
  //   val_bas_fct_neigh[0].resize(0);
  compute_DG_bilinear_form_local<d>(param_db, base_vec_dim, val_bas_fct, val_bas_fct_dx,
      val_bas_fct_dy, val_bas_fct_dz, nullptr, nullptr, nullptr, nullptr,
      dummy_index_map, normal, joint_diam, quad_form_joint, orig_weights,
      coefficients, true, loc_mat);
}
#ifdef __3D__
template void compute_DG_bilinear_form_local<3>(const ParameterDatabase&
    param_db, const int& base_vec_dim, const std::vector<double>* val_bas_fct, const std::vector<double>*
    val_bas_fct_dx, const std::vector<double>* val_bas_fct_dy, const
    std::vector<double>* val_bas_fct_dz, const double* normal, const double&
    joint_diam, const TQuadFormula* quad_form_joint, const std::vector<double>&
    orig_weights, double** coefficients, double** loc_mat);
#else
template void compute_DG_bilinear_form_local<2>(const ParameterDatabase&
    param_db, const int& base_vec_dim, const std::vector<double>* val_bas_fct, const std::vector<double>*
    val_bas_fct_dx, const std::vector<double>* val_bas_fct_dy, const
    std::vector<double>* val_bas_fct_dz, const double* normal, const double&
    joint_diam, const TQuadFormula* quad_form_joint, const std::vector<double>&
    orig_weights, double** coefficients, double** loc_mat);
#endif

// Computes the local DG bilinear form for Navier-Stokes problems
template< int d>
void compute_DG_bilinear_form_local_NSE(const ParameterDatabase& param_db, const int& base_vec_dim, const
    std::vector<double>* val_bas_fct, const std::vector<double>* val_bas_fct_dx,
    const std::vector<double>* val_bas_fct_dy, const std::vector<double>*
    /*val_bas_fct_dz*/, const std::vector<double>* val_bas_fct_neigh, const
    std::vector<double>* val_bas_fct_neigh_dx, const std::vector<double>*
    val_bas_fct_neigh_dy, const std::vector<double>* /*val_bas_fct_neigh_dz*/, const
    std::vector<int>& index_map, const double* joint_normal, const double&
    joint_diam, const TQuadFormula* quad_form_joint, const std::vector<double>&
    orig_weights, double** coefficients, const bool& is_boundary_edge, double**
    loc_mat, const int* global_dof, const int* global_dof_neigh, const TBaseCell*
    /*cell*/, int /*cell_no*/, int /*egde_no*/)
{

  // get some information needed for integration
  auto n_quad_pts = quad_form_joint->GetN_QuadPoints();
  auto n_dofs = val_bas_fct[0].size()/base_vec_dim;
  auto n_dofs_neigh = (is_boundary_edge) ? 0 : val_bas_fct_neigh[0].size()/base_vec_dim;
  auto n_dofs_total = n_dofs + n_dofs_neigh;
  auto nx = joint_normal[0];
  auto ny = joint_normal[1];
  //auto nz = (d == 3) ? joint_normal[2] : 0.0;

  // get parameters of DG method
  auto face_sigma = param_db["face_sigma_DG"].get<double>(); // edge penalty
  auto kappa = param_db["symmetry_DG"].get<int>();  // symmetry
  
  for (int quad_pt_i = 0; quad_pt_i < n_quad_pts; ++quad_pt_i)
  {
    //Output::print(" quad point ", quad_pt_i, "   ", quad_form_joint->get_point(quad_pt_i).x);
    // parameter of PDE
    auto viscosity = coefficients[quad_pt_i][0];


    for (unsigned int t_dof = 0; t_dof < n_dofs_total; ++t_dof)
    { // for all test dofs

      // Define jump and averages for test functions index [0] corresponds to
      // the basis functions itself index [1/2/3] corresponds to the derivative
      // wrt. x/y/z for the basis functions.
      // No matter if we are in 2D or 3D we always assume 4 values, i.e. the
      // 4 indices mentioned above, and in 2D we set every entry corresponding
      // to z to 0.
      std::vector<double> t_jump(3 * base_vec_dim);
      std::vector<double> t_average(3 * base_vec_dim);


      // distinguish whether test function are related to the cell or to the
      // neighbor. Compute jump and average  accordingly using the fact that the
      // support of a basis function is exactly the cell it lives on.

      if (t_dof < n_dofs)
      {
        // find if the current t_dof is in the global_dof_neigh and it's index
        // i.e. figure out if t_dof is defined on the shared edge or not
          int global_pos = (n_dofs_neigh > 0) ? global_dof[t_dof] : 0;
          int global_pos_neigh_index = -1;

          for (unsigned int j = 0; j < n_dofs_neigh; ++j)
          {
            if (global_pos == global_dof_neigh[j])
                {
                    global_pos_neigh_index = j;
                    break;
                }
          }

          if (global_pos_neigh_index == -1)
          { // test dof on the boundary edge of the cell
            t_jump[0] = val_bas_fct [quad_pt_i][t_dof];
            t_jump[1] = val_bas_fct_dx [quad_pt_i][t_dof];
            t_jump[2] = val_bas_fct_dy [quad_pt_i][t_dof];
            t_jump[3] = val_bas_fct [quad_pt_i][n_dofs + t_dof];
            t_jump[4] = val_bas_fct_dx [quad_pt_i][n_dofs + t_dof];
            t_jump[5] = val_bas_fct_dy [quad_pt_i][n_dofs + t_dof];

            t_average[0] = val_bas_fct [quad_pt_i][t_dof] / 2;
            t_average[1] = val_bas_fct_dx [quad_pt_i][t_dof] / 2;
            t_average[2] = val_bas_fct_dy [quad_pt_i][t_dof] / 2;
            t_average[3] = val_bas_fct [quad_pt_i][n_dofs + t_dof] / 2;
            t_average[4] = val_bas_fct_dx [quad_pt_i][n_dofs + t_dof] / 2;
            t_average[5] = val_bas_fct_dy [quad_pt_i][n_dofs + t_dof] / 2;

          }
          else
          { // test dof on the shared edge
            t_jump[0] = val_bas_fct[quad_pt_i][t_dof] - val_bas_fct_neigh[index_map[quad_pt_i]][global_pos_neigh_index];
            t_jump[1] = val_bas_fct_dx[quad_pt_i][t_dof] - val_bas_fct_neigh_dx[index_map[quad_pt_i]][global_pos_neigh_index];
            t_jump[2] = val_bas_fct_dy[quad_pt_i][t_dof] - val_bas_fct_neigh_dy[index_map[quad_pt_i]][global_pos_neigh_index];
            t_jump[3] = val_bas_fct[quad_pt_i][n_dofs + t_dof] - val_bas_fct_neigh[index_map[quad_pt_i]][n_dofs_neigh + global_pos_neigh_index];
            t_jump[4] = val_bas_fct_dx[quad_pt_i][n_dofs + t_dof] - val_bas_fct_neigh_dx[index_map[quad_pt_i]][n_dofs_neigh + global_pos_neigh_index];
            t_jump[5] = val_bas_fct_dy[quad_pt_i][n_dofs + t_dof] - val_bas_fct_neigh_dy[index_map[quad_pt_i]][n_dofs_neigh + global_pos_neigh_index];


            t_average[0] = (val_bas_fct[quad_pt_i][t_dof] + val_bas_fct_neigh[index_map[quad_pt_i]][global_pos_neigh_index]) / 2;
            t_average[1] = (val_bas_fct_dx[quad_pt_i][t_dof] + val_bas_fct_neigh_dx[index_map[quad_pt_i]][global_pos_neigh_index]) / 2;
            t_average[2] = (val_bas_fct_dy[quad_pt_i][t_dof] + val_bas_fct_neigh_dy[index_map[quad_pt_i]][global_pos_neigh_index]) / 2;
            t_average[3] = (val_bas_fct[quad_pt_i][n_dofs + t_dof] + val_bas_fct_neigh[index_map[quad_pt_i]][n_dofs_neigh + global_pos_neigh_index]) / 2;
            t_average[4] = (val_bas_fct_dx[quad_pt_i][n_dofs + t_dof] + val_bas_fct_neigh_dx[index_map[quad_pt_i]][n_dofs_neigh + global_pos_neigh_index]) / 2;
            t_average[5] = (val_bas_fct_dy[quad_pt_i][n_dofs + t_dof] + val_bas_fct_neigh_dy[index_map[quad_pt_i]][n_dofs_neigh + global_pos_neigh_index]) / 2;
          }
      }
      else
      {
        // find if the current t_dof is in the global_dof and it's index
        // i.e. figure out if t_dof is defined on the shared edge or not
          int global_pos_neigh = global_dof_neigh[t_dof - n_dofs];
          int global_pos_index = -1;

          for (unsigned int j = 0; j < n_dofs; ++j)
          {
            if (global_pos_neigh == global_dof[j])
                {
                    global_pos_index = j;
                    break;
                }
          }


        // test dof on the boundary edges from the neighbour
        if (global_pos_index == -1)
        {
            t_jump[0] = -val_bas_fct_neigh[index_map[quad_pt_i]][t_dof - n_dofs];
            t_jump[1] = -val_bas_fct_neigh_dx[index_map[quad_pt_i]][t_dof - n_dofs];
            t_jump[2] = -val_bas_fct_neigh_dy[index_map[quad_pt_i]][t_dof - n_dofs];
            t_jump[3] = -val_bas_fct_neigh[index_map[quad_pt_i]][n_dofs_neigh + t_dof - n_dofs];
            t_jump[4] = -val_bas_fct_neigh_dx[index_map[quad_pt_i]][n_dofs_neigh + t_dof - n_dofs];
            t_jump[5] = -val_bas_fct_neigh_dy[index_map[quad_pt_i]][n_dofs_neigh + t_dof - n_dofs];


            t_average[0] = val_bas_fct_neigh[index_map[quad_pt_i]][t_dof - n_dofs] / 2;
            t_average[1] = val_bas_fct_neigh_dx[index_map[quad_pt_i]][t_dof - n_dofs] / 2;
            t_average[2] = val_bas_fct_neigh_dy[index_map[quad_pt_i]][t_dof - n_dofs] / 2;
            t_average[3] = val_bas_fct_neigh[index_map[quad_pt_i]][n_dofs_neigh + t_dof - n_dofs] / 2;
            t_average[4] = val_bas_fct_neigh_dx[index_map[quad_pt_i]][n_dofs_neigh + t_dof - n_dofs] / 2;
            t_average[5] = val_bas_fct_neigh_dy[index_map[quad_pt_i]][n_dofs_neigh + t_dof - n_dofs] / 2;
        }
        else
        {
        // test dof on the shared edge with the neighbour is 0
            t_jump[0] = 0;
            t_jump[1] = 0;
            t_jump[2] = 0;
            t_jump[3] = 0;
            t_jump[4] = 0;
            t_jump[5] = 0;

            t_average[0] = 0;
            t_average[1] = 0;
            t_average[2] = 0;
            t_average[3] = 0;
            t_average[4] = 0;
            t_average[5] = 0;

        }
      }

      if (is_boundary_edge)
        { // Correct values for boundary edges: jump and average on the boundary
          // are defined as values seen from cell
          for (int i = 0; i < 6; ++i)
          {
            t_average[i] *= 2;
          }
        }


      for (unsigned int a_dof = 0; a_dof < n_dofs_total; ++a_dof)
      { // for all ansatz dofs

        // Fix the descriptions at some point
        // Define jump and averages for ansatz functions index [0] corresponds
        // to the basis functions itself index [1/2/3] corresponds to the
        // derivative wrt. x/y/z for the basis functions.  No matter if we are
        // in 2D or 3D we always assume 4 values, i.e. the
        // 4 indices mentioned above, and in 2D we set every entry corresponding
        // to z to 0.
        std::vector<double> a_jump(3 * base_vec_dim);
        std::vector<double> a_average(3 * base_vec_dim);

        // distinguish whether ansatz function are related to the cell or to the
        // neighbor. Compute jump and average  accordingly using the fact that
        // the support of a basis function is exactly the cell it lives on.
        if (a_dof < n_dofs)
        {
            // find if the current a_dof is in the global_dof_neigh and it's index
            // i.e. figure out if a_dof is defined on the shared edge or not
              int global_pos = (n_dofs_neigh > 0) ? global_dof[a_dof] : 0;
              int global_pos_neigh_index = -1;

              for (unsigned int j = 0; j < n_dofs_neigh; ++j)
              {
                if (global_pos == global_dof_neigh[j])
                    {
                        global_pos_neigh_index = j;
                        break;
                    }
              }

              if (global_pos_neigh_index == -1)
              { // test dof on the boundary edge of the cell
                a_jump[0] = val_bas_fct [quad_pt_i][a_dof];
                a_jump[1] = val_bas_fct_dx [quad_pt_i][a_dof];
                a_jump[2] = val_bas_fct_dy [quad_pt_i][a_dof];
                a_jump[3] = val_bas_fct [quad_pt_i][n_dofs + a_dof];
                a_jump[4] = val_bas_fct_dx [quad_pt_i][n_dofs + a_dof];
                a_jump[5] = val_bas_fct_dy [quad_pt_i][n_dofs + a_dof];

                a_average[0] = val_bas_fct [quad_pt_i][a_dof] / 2;
                a_average[1] = val_bas_fct_dx [quad_pt_i][a_dof] / 2;
                a_average[2] = val_bas_fct_dy [quad_pt_i][a_dof] / 2;
                a_average[3] = val_bas_fct [quad_pt_i][n_dofs + a_dof] / 2;
                a_average[4] = val_bas_fct_dx [quad_pt_i][n_dofs + a_dof] / 2;
                a_average[5] = val_bas_fct_dy [quad_pt_i][n_dofs + a_dof] / 2;
               }
              else
              { // ansatz dof on the shared edge
                a_jump[0] = val_bas_fct[quad_pt_i][a_dof] - val_bas_fct_neigh[index_map[quad_pt_i]][global_pos_neigh_index];
                a_jump[1] = val_bas_fct_dx[quad_pt_i][a_dof] - val_bas_fct_neigh_dx[index_map[quad_pt_i]][global_pos_neigh_index];
                a_jump[2] = val_bas_fct_dy[quad_pt_i][a_dof] - val_bas_fct_neigh_dy[index_map[quad_pt_i]][global_pos_neigh_index];
                a_jump[3] = val_bas_fct[quad_pt_i][n_dofs + a_dof] - val_bas_fct_neigh[index_map[quad_pt_i]][n_dofs_neigh + global_pos_neigh_index];
                a_jump[4] = val_bas_fct_dx[quad_pt_i][n_dofs + a_dof] - val_bas_fct_neigh_dx[index_map[quad_pt_i]][n_dofs_neigh + global_pos_neigh_index];
                a_jump[5] = val_bas_fct_dy[quad_pt_i][n_dofs + a_dof] - val_bas_fct_neigh_dy[index_map[quad_pt_i]][n_dofs_neigh + global_pos_neigh_index];


                a_average[0] = (val_bas_fct[quad_pt_i][a_dof] + val_bas_fct_neigh[index_map[quad_pt_i]][global_pos_neigh_index]) / 2;
                a_average[1] = (val_bas_fct_dx[quad_pt_i][a_dof] + val_bas_fct_neigh_dx[index_map[quad_pt_i]][global_pos_neigh_index]) / 2;
                a_average[2] = (val_bas_fct_dy[quad_pt_i][a_dof] + val_bas_fct_neigh_dy[index_map[quad_pt_i]][global_pos_neigh_index]) / 2;
                a_average[3] = (val_bas_fct[quad_pt_i][n_dofs + a_dof] + val_bas_fct_neigh[index_map[quad_pt_i]][n_dofs_neigh + global_pos_neigh_index]) / 2;
                a_average[4] = (val_bas_fct_dx[quad_pt_i][n_dofs + a_dof] + val_bas_fct_neigh_dx[index_map[quad_pt_i]][n_dofs_neigh + global_pos_neigh_index]) / 2;
                a_average[5] = (val_bas_fct_dy[quad_pt_i][n_dofs + a_dof] + val_bas_fct_neigh_dy[index_map[quad_pt_i]][n_dofs_neigh + global_pos_neigh_index]) / 2;
              }
      }
      else
      {
        // find if the current a_dof is in the global_dof and it's index
        // i.e. figure out if a_dof is defined on the shared edge or not
          int global_pos_neigh = global_dof_neigh[a_dof - n_dofs];
          int global_pos_index = -1;

          for (unsigned int j = 0; j < n_dofs; ++j)
          {
            if (global_pos_neigh == global_dof[j])
                {
                    global_pos_index = j;
                    break;
                }
          }


        // ansatz dof on the boundary edges from the neighbour
        if (global_pos_index == -1)
        {
            a_jump[0] = -val_bas_fct_neigh[index_map[quad_pt_i]][a_dof - n_dofs];
            a_jump[1] = -val_bas_fct_neigh_dx[index_map[quad_pt_i]][a_dof - n_dofs];
            a_jump[2] = -val_bas_fct_neigh_dy[index_map[quad_pt_i]][a_dof - n_dofs];
            a_jump[3] = -val_bas_fct_neigh[index_map[quad_pt_i]][n_dofs_neigh + a_dof - n_dofs];
            a_jump[4] = -val_bas_fct_neigh_dx[index_map[quad_pt_i]][n_dofs_neigh + a_dof - n_dofs];
            a_jump[5] = -val_bas_fct_neigh_dy[index_map[quad_pt_i]][n_dofs_neigh + a_dof - n_dofs];


            a_average[0] = val_bas_fct_neigh[index_map[quad_pt_i]][a_dof - n_dofs] / 2;
            a_average[1] = val_bas_fct_neigh_dx[index_map[quad_pt_i]][a_dof - n_dofs] / 2;
            a_average[2] = val_bas_fct_neigh_dy[index_map[quad_pt_i]][a_dof - n_dofs] / 2;
            a_average[3] = val_bas_fct_neigh[index_map[quad_pt_i]][n_dofs_neigh + a_dof - n_dofs] / 2;
            a_average[4] = val_bas_fct_neigh_dx[index_map[quad_pt_i]][n_dofs_neigh + a_dof - n_dofs] / 2;
            a_average[5] = val_bas_fct_neigh_dy[index_map[quad_pt_i]][n_dofs_neigh + a_dof - n_dofs] / 2;
        }
        else
        {
        // test dof on the shared edge with the neighbour is 0
            a_jump[0] = 0;
            a_jump[1] = 0;
            a_jump[2] = 0;
            a_jump[3] = 0;
            a_jump[4] = 0;
            a_jump[5] = 0;

            a_average[0] = 0;
            a_average[1] = 0;
            a_average[2] = 0;
            a_average[3] = 0;
            a_average[4] = 0;
            a_average[5] = 0;
        }

      }

       if (is_boundary_edge)
        { // Correct values for boundary edges: jump and average on the boundary
          // are defined as values seen from cell
          for (int i = 0; i < 6; ++i)
          {
            a_average[i] *= 2;
          }
        }
        
//         if(t_dof == 0)
//         {
//         Output::print("   grad u ", a_average[1], "  ", a_average[2], "  ", a_average[4], "  ", a_average[5]);
//         }
//         if(a_dof == 0)
//         {
//           Output::print("   t_jump ", t_jump[0], "  ", t_jump[3]);
//         }
        double tx = -ny;
        double ty = nx;
        double t_jump_t = t_jump[0]*tx + t_jump[3]*ty; // tangential jump of test function
        double a_jump_t = a_jump[0]*tx + a_jump[3]*ty; // tangential jump of ansatz function
        double integrand = 0;
        integrand -=
               (a_average[1] * nx + a_average[2] * ny) * tx * t_jump_t
             + (a_average[4] * nx + a_average[5] * ny) * ty * t_jump_t
            // symmetry term
             + kappa * (a_jump_t * tx * (t_average[1] * nx + t_average[2] * ny)
             + a_jump_t * ty * (t_average[4] * nx + t_average[5] * ny));
             
             
        // convective term for Navier-Stokes  
        if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE == 5) // Navier-Stokes (3 means Stokes)
      {   
		auto inflow = val_bas_fct[quad_pt_i][t_dof] * nx + val_bas_fct[quad_pt_i][t_dof + n_dofs] * ny;
		//Output::print("boundary_edge:", is_boundary_edge ,"   inflow: ", inflow);
		
		//not sure why the term below is needed 
		/*if (is_boundary_edge)
		{	
		    	
		    if (inflow < 0)
		    {
			integrand -= inflow * a_jump_t  * (t_average[0] * tx + t_average[3] * ty);
		    }
		}
		else
		{*/
		      integrand += -inflow * a_jump_t  * (t_average[0] * tx + t_average[3] * ty) + 0.5 * std::abs(inflow) * a_jump_t * t_jump_t;
		//}
	}

        // penalty term to achieve coercivity
        integrand += face_sigma / joint_diam * a_jump_t * t_jump_t;


        auto weight = orig_weights[quad_pt_i];
        // Update local matrix
        loc_mat[t_dof][a_dof] += viscosity * weight * integrand;
      } // endfor ansatz dofs
    } // endfor test dofs
  } //endfor quad points
} // end compute_DG_bilinear_form_local_NSE

#ifdef __3D__
template void compute_DG_bilinear_form_local_NSE<3>(const ParameterDatabase&
    param_db, const int& base_vec_dim, const std::vector<double>* val_bas_fct, const std::vector<double>*
    val_bas_fct_dx, const std::vector<double>* val_bas_fct_dy, const
    std::vector<double>* val_bas_fct_dz, const std::vector<double>*
    val_bas_fct_neigh, const std::vector<double>* val_bas_fct_neigh_dx, const
    std::vector<double>* val_bas_fct_neigh_dy, const std::vector<double>*
    val_bas_fct_neigh_dz, const std::vector<int>& index_map, const double*
    joint_normal, const double& joint_diam, const TQuadFormula* quad_form_joint,
    const std::vector<double>& orig_weights, double** coefficients, const bool&
    is_boundary_edge, double** loc_mat, const int* global_dof, const int* global_dof_neigh, const TBaseCell* cell, int cell_no, int egde_no);
#else
template void compute_DG_bilinear_form_local_NSE<2>(const ParameterDatabase&
    param_db, const int& base_vec_dim, const std::vector<double>* val_bas_fct, const std::vector<double>*
    val_bas_fct_dx, const std::vector<double>* val_bas_fct_dy, const
    std::vector<double>* val_bas_fct_dz, const std::vector<double>*
    val_bas_fct_neigh, const std::vector<double>* val_bas_fct_neigh_dx, const
    std::vector<double>* val_bas_fct_neigh_dy, const std::vector<double>*
    val_bas_fct_neigh_dz, const std::vector<int>& index_map, const double*
    joint_normal, const double& joint_diam, const TQuadFormula* quad_form_joint,
    const std::vector<double>& orig_weights, double** coefficients, const bool&
    is_boundary_edge, double** loc_mat, const int* global_dof, const int* global_dof_neigh, const TBaseCell* cell, int cell_no, int egde_no);
#endif

template< int d>
void compute_DG_bilinear_form_local_NSE(const ParameterDatabase& param_db, const int& base_vec_dim, const
    std::vector<double>* val_bas_fct, const std::vector<double>* val_bas_fct_dx,
    const std::vector<double>* val_bas_fct_dy, const std::vector<double>*
    val_bas_fct_dz, const double* normal, const double& joint_diam, const
    TQuadFormula* quad_form_joint, const std::vector<double>& orig_weights,
    double** coefficients, double** loc_mat, const int* global_dof)
{
  std::vector<int> dummy_index_map;
  //   val_bas_fct_neigh[0].resize(0);
  compute_DG_bilinear_form_local_NSE<d>(param_db, base_vec_dim, val_bas_fct, val_bas_fct_dx,
      val_bas_fct_dy, val_bas_fct_dz, nullptr, nullptr, nullptr, nullptr,
      dummy_index_map, normal, joint_diam, quad_form_joint, orig_weights,
      coefficients, true, loc_mat, global_dof, nullptr, nullptr, 0, 0);
}
#ifdef __3D__
template void compute_DG_bilinear_form_local_NSE<3>(const ParameterDatabase&
    param_db, const int& base_vec_dim, const std::vector<double>* val_bas_fct, const std::vector<double>*
    val_bas_fct_dx, const std::vector<double>* val_bas_fct_dy, const
    std::vector<double>* val_bas_fct_dz, const double* normal, const double&
    joint_diam, const TQuadFormula* quad_form_joint, const std::vector<double>&
    orig_weights, double** coefficients, double** loc_mat, const int* global_dof);
#else
template void compute_DG_bilinear_form_local_NSE<2>(const ParameterDatabase&
    param_db, const int& base_vec_dim, const std::vector<double>* val_bas_fct, const std::vector<double>*
    val_bas_fct_dx, const std::vector<double>* val_bas_fct_dy, const
    std::vector<double>* val_bas_fct_dz, const double* normal, const double&
    joint_diam, const TQuadFormula* quad_form_joint, const std::vector<double>&
    orig_weights, double** coefficients, double** loc_mat, const int* global_dof);
#endif

// End of changes for the Navier Stokes bilinear form




// Computes the local DG right-hand side for convection-diffusion problems
template< int d >
void compute_DG_rhs_local(const ParameterDatabase& param_db, const TBaseCell*
    cell, const int& joint_index, const FiniteElement& finite_element, const
    std::vector<double>* val_bas_fct, const std::vector<double>* val_bas_fct_dx,
    const std::vector<double>* val_bas_fct_dy, const std::vector<double>*
    val_bas_fct_dz, const double* joint_normal, const double& joint_diam, const
    TQuadFormula* quad_form_facet, const std::vector<double>& orig_weights,
    double** coefficients, const typename
    Template_names<d>::BoundaryValuesFunction* boundary_value, double* loc_rhs)
{
  // get some information needed for integration
  auto n_quad_pts = quad_form_facet->GetN_QuadPoints();
  auto n_dofs = finite_element.GetN_DOF();

  // get parameters of DG method
  auto face_sigma = param_db["face_sigma_DG"].get<double>(); // edge penalty
  auto kappa = param_db["symmetry_DG"].get<int>();  // symmetry

  // geometrical information
  auto nx = joint_normal[0];
  auto ny = joint_normal[1];
  auto nz = (d == 3) ? joint_normal[2] : 0;

  auto boundedge = (const TBoundEdge*) cell->GetJoint(joint_index);
  auto boundary_comp = boundedge->GetBoundComp();
  auto component = boundary_comp->GetID();

  // Prepare evaluation of boundary condition for 3D case
  BFRefElements ref_element;
  ReferenceTransformation_type ref_trans_id;
  std::vector<parmoon::Point> transformed_quad_points;
  if (d == 3)
  {
    ref_element = finite_element.GetBaseFunct()->GetRefElement();
    ref_trans_id = finite_element.GetRefTransID();
    transformed_quad_points = transform_quad_points(cell, ref_element,ref_trans_id, joint_index, quad_form_facet);
  }

  for (int quad_pt_i = 0; quad_pt_i < n_quad_pts; ++quad_pt_i)
  {
    // get value in this quadrature point (in val_bdry_cond)
    double val_bdry_cond = 0;
#ifdef __2D__
    // prepare computation of right-hand side: compute parameters of joint
    double t0, t1;
    boundedge->GetParameters(t0, t1);

    // compute quadrature point located on the boundary
    auto zeta_1D = quad_form_facet->get_point(quad_pt_i).x;
    auto quad_pt_bdry = t0 + 0.5 * (t1 - t0) * (zeta_1D + 1);
    boundary_value(component, quad_pt_bdry, val_bdry_cond);
#else
    boundary_value(component, transformed_quad_points[quad_pt_i].x,
        transformed_quad_points[quad_pt_i].y,
        transformed_quad_points[quad_pt_i].z, val_bdry_cond);
#endif

    // Parameter of PDE
    auto diffusion = coefficients[quad_pt_i][0];
    auto beta_1 = coefficients[quad_pt_i][1];
    auto beta_2 = coefficients[quad_pt_i][2];
    auto convection_3 = (d == 3) ? coefficients[quad_pt_i][3] : 0;

    for (int t_dof = 0; t_dof < n_dofs; ++t_dof)
    { // for all test dofs

      double integrand = 0;
      // contribution of diffusive term
      integrand -= diffusion * kappa * (val_bas_fct_dx[quad_pt_i][t_dof] * nx
          + val_bas_fct_dy[quad_pt_i][t_dof] * ny) * val_bdry_cond;
      integrand -= (d==3) ? diffusion * kappa * val_bdry_cond
        * val_bas_fct_dz[quad_pt_i][t_dof] * nz : 0;

      // penalty term to achieve coercivity (factor 2 for boundary edges)
      integrand += 2 * face_sigma / joint_diam * val_bdry_cond *
        val_bas_fct[quad_pt_i][t_dof];

      // contribution of convective term
      auto inflow = beta_1 * nx + beta_2 * ny + convection_3 * nz;
      if (inflow < 0)
      {
        integrand -= inflow * val_bdry_cond * val_bas_fct[quad_pt_i][t_dof];
      }

      auto weight = orig_weights[quad_pt_i];
      loc_rhs[t_dof] += weight * integrand;
    } // endfor test dofs
  } // endfor quad_pt_i
} // end compute_DG_rhs_local
#ifdef __3D__
template void compute_DG_rhs_local<3>(const ParameterDatabase& param_db, const
    TBaseCell* cell, const int& joint_index, const FiniteElement& finite_element,
    const std::vector<double>* val_bas_fct, const std::vector<double>*
    val_bas_fct_dx, const std::vector<double>* val_bas_fct_dy, const
    std::vector<double>* val_bas_fct_dz, const double* joint_normal, const
    double& joint_diam, const TQuadFormula* quad_form_facet, const
    std::vector<double>& orig_weights, double** coefficients, const typename
    Template_names<3>::BoundaryValuesFunction* boundary_value, double* loc_rhs);
#else
template void compute_DG_rhs_local<2>(const ParameterDatabase& param_db, const
    TBaseCell* cell, const int& joint_index, const FiniteElement& finite_element,
    const std::vector<double>* val_bas_fct, const std::vector<double>*
    val_bas_fct_dx, const std::vector<double>* val_bas_fct_dy, const
    std::vector<double>* val_bas_fct_dz, const double* joint_normal, const
    double& joint_diam, const TQuadFormula* quad_form_facet, const
    std::vector<double>& orig_weights, double** coefficients, const typename
    Template_names<2>::BoundaryValuesFunction* boundary_value, double* loc_rhs);
#endif

// Computes the local DG right-hand side for Navier-Stokes problems
template< int d >
void compute_DG_rhs_local_NSE(const ParameterDatabase& param_db, TBaseCell*
    cell, const int& joint_index, const FiniteElement& finite_element, const
    std::vector<double>* val_bas_fct, const std::vector<double>* val_bas_fct_dx,
    const std::vector<double>* val_bas_fct_dy, const std::vector<double>*
    /*val_bas_fct_dz*/, const double* /*joint_normal*/, const double& joint_diam, const
    TQuadFormula* quad_form_facet, const std::vector<double>& orig_weights,
    double** coefficients, const typename
    Template_names<d>::BoundaryValuesFunction* boundary_value, double* loc_rhs)
{
     // get some information needed for integration
  auto n_quad_pts = quad_form_facet->GetN_QuadPoints();
  auto n_dofs = finite_element.GetN_DOF();

  // get parameters of DG method
  auto face_sigma = param_db["face_sigma_DG"].get<double>(); // edge penalty
  auto kappa = param_db["symmetry_DG"].get<int>();  // symmetry

  auto boundedge = (const TBoundEdge*) cell->GetJoint(joint_index);
  //Output::print("boundedge: ", boundedge);
  auto boundary_comp = boundedge->GetBoundComp();
  //Output::print("boundary_comp: ", boundary_comp);
  auto component = boundary_comp->GetID();
  auto normal = cell->ComputeOuterJointNormal(joint_index);
  auto tangent_x = (-1)*normal[1];
  auto tangent_y = normal[0];

  // Prepare evaluation of boundary condition for 3D case
  BFRefElements ref_element;
  ReferenceTransformation_type ref_trans_id;
  std::vector<parmoon::Point> transformed_quad_points;
  if (d == 3)
  {
    ref_element = finite_element.GetBaseFunct()->GetRefElement();
    ref_trans_id = finite_element.GetRefTransID();
    transformed_quad_points = transform_quad_points(cell, ref_element,ref_trans_id, joint_index, quad_form_facet);
  }

  for (int quad_pt_i = 0; quad_pt_i < n_quad_pts; ++quad_pt_i)
  {
    double viscosity = coefficients[quad_pt_i][0];
    // get value in this quadrature point (in val_bdry_cond)
    double val_bdry_cond = 0;
#ifdef __2D__
    // prepare computation of right-hand side: compute parameters of joint
    double t0, t1;
    boundedge->GetParameters(t0, t1);

    // compute quadrature point located on the boundary
    auto zeta_1D = quad_form_facet->get_point(quad_pt_i).x;
    auto quad_pt_bdry = t0 + 0.5 * (t1 - t0) * (zeta_1D + 1);
    boundary_value(component, quad_pt_bdry, val_bdry_cond);
#else
    boundary_value(component, transformed_quad_points[quad_pt_i].x,
        transformed_quad_points[quad_pt_i].y,
        transformed_quad_points[quad_pt_i].z, val_bdry_cond);
#endif

    for (int t_dof = 0; t_dof < n_dofs; ++t_dof)
    { // for all test dofs
      double integrand = 0;
      double vt = (val_bas_fct[quad_pt_i][t_dof]*tangent_x + val_bas_fct[quad_pt_i][t_dof+n_dofs]*tangent_y);
      //Output::print("  t_dof " , t_dof, "  vt ", vt, "  gt ", val_bdry_cond);
      // penalty term to achieve coercivity (factor 2 for boundary edges)
      integrand +=  face_sigma / joint_diam * val_bdry_cond * vt;
      
      double gradv_n1 = val_bas_fct_dx[quad_pt_i][t_dof] * normal[0]
                      + val_bas_fct_dy[quad_pt_i][t_dof] * normal[1];
      double gradv_n2 = val_bas_fct_dx[quad_pt_i][t_dof+n_dofs] * normal[0]
                      + val_bas_fct_dy[quad_pt_i][t_dof+n_dofs] * normal[1];
      integrand -= kappa * val_bdry_cond * (gradv_n1*tangent_x + gradv_n2*tangent_y);
      
      // contribution of the convective term for Navier Stokes
      /*if (boundedge)
      {
	      if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE == 5) // Navier-Stokes (3 means Stokes)
	      {
		      auto inflow = val_bas_fct [quad_pt_i][t_dof] * normal[0] + val_bas_fct[quad_pt_i][t_dof + n_dofs] * normal[1];
		      Output::print("inflow: ", inflow);
		      if (inflow < 0)
		      {
			 integrand -= inflow * val_bdry_cond * vt;
		       }
	       }  
	}	       */     

      auto weight = orig_weights[quad_pt_i];
      loc_rhs[t_dof] += viscosity * weight * integrand;
    } // endfor test dofs
  } // endfor quad_pt_i


} // end compute_DG_rhs_local_NSE
#ifdef __3D__
template void compute_DG_rhs_local_NSE<3>(const ParameterDatabase& param_db,
    TBaseCell* cell, const int& joint_index, const FiniteElement& finite_element,
    const std::vector<double>* val_bas_fct, const std::vector<double>*
    val_bas_fct_dx, const std::vector<double>* val_bas_fct_dy, const
    std::vector<double>* val_bas_fct_dz, const double* joint_normal, const
    double& joint_diam, const TQuadFormula* quad_form_facet, const
    std::vector<double>& orig_weights, double** coefficients, const typename
    Template_names<3>::BoundaryValuesFunction* boundary_value, double* loc_rhs);
#else
template void compute_DG_rhs_local_NSE<2>(const ParameterDatabase& param_db,
    TBaseCell* cell, const int& joint_index, const FiniteElement& finite_element,
    const std::vector<double>* val_bas_fct, const std::vector<double>*
    val_bas_fct_dx, const std::vector<double>* val_bas_fct_dy, const
    std::vector<double>* val_bas_fct_dz, const double* joint_normal, const
    double& joint_diam, const TQuadFormula* quad_form_facet, const
    std::vector<double>& orig_weights, double** coefficients, const typename
    Template_names<2>::BoundaryValuesFunction* boundary_value, double* loc_rhs);
#endif
// End of changes for the Navier Stokes r.h.s

template< int d >
void copy_local_to_global(const int& n_dofs, const int& n_dofs_neigh, const int*
    global_dof, const int* global_dof_neigh, double** loc_mat, double* loc_rhs,
    typename Template_names<d>::SquareMatrixD* global_matrix, double*
    global_rhs, double** loc_mata22, typename Template_names<d>::SquareMatrixD* global_matrixA22)
{
  // Check consistency of input
  if ( (loc_rhs == nullptr && global_rhs != nullptr ) ||
      (loc_rhs != nullptr && global_rhs == nullptr ) )
  {
    ErrThrow("Don't know what to do if a local right-hand side is given but no",
        " global one or vice versa.");
  }

  // Does a right-hand side appear?
  bool is_rhs = (loc_rhs != nullptr);
  // Does a neighbor cell play a role?
  auto n_dof_total = (global_dof_neigh == nullptr) ? n_dofs :
    n_dofs + n_dofs_neigh;

  for (int t_dof = 0; t_dof < n_dof_total ; ++t_dof)
  {
    auto t_global_number = (t_dof < n_dofs) ? global_dof[t_dof] :
      global_dof_neigh[t_dof - n_dofs];  // global dof number
    if(t_global_number >= global_matrix->get_n_active_rows())
        continue;    
    auto loc_row = loc_mat[t_dof];
    
    double *loc_rowa22;
    if(global_matrixA22 != nullptr)
    {
      loc_rowa22 = loc_mata22[t_dof];
    }
    for (int a_dof = 0; a_dof < n_dof_total; ++a_dof)
    {
      auto a_global_number = (a_dof < n_dofs) ? global_dof[a_dof] :
        global_dof_neigh[a_dof - n_dofs];  // global dof number
      global_matrix->add(t_global_number, a_global_number, loc_row[a_dof]);
      if(global_matrixA22 != nullptr)
      {
        global_matrixA22->add(t_global_number, a_global_number, loc_rowa22[a_dof]);
      }
    } // endfor ansatz dofs
    if (is_rhs && t_dof < n_dofs)
    {
      global_rhs[t_global_number] += loc_rhs[t_dof];
    }
  } // endfor test dofs
} //end copy_local_to_global

#ifdef __3D__
template void copy_local_to_global<3>(const int& n_dofs, const int&
    n_dofs_neigh, const int* global_dof, const int* global_dof_neigh, double**
    loc_mat, double* loc_rhs, typename Template_names<3>::SquareMatrixD*
    global_matrix, double* global_rhs, double** loc_mata22,
    typename Template_names<3>::SquareMatrixD* global_matrixA22);
#else
template void copy_local_to_global<2>(const int& n_dofs, const int&
    n_dofs_neigh, const int* global_dof, const int* global_dof_neigh, double**
    loc_mat, double* loc_rhs, typename Template_names<2>::SquareMatrixD*
    global_matrix, double* global_rhs, double** loc_mata22, 
    typename Template_names<2>::SquareMatrixD* global_matrixA22);
#endif


template< int d >
void copy_local_to_global(const int& n_dofs, const int* global_dof, double**
    loc_mat, double* loc_rhs, typename Template_names<d>::SquareMatrixD*
    global_matrix, double* global_rhs)
{
  copy_local_to_global<d>(n_dofs, 0, global_dof, nullptr, loc_mat, loc_rhs,
      global_matrix, global_rhs);
}

#ifdef __3D__
template void copy_local_to_global<3>(const int& n_dofs, const int* global_dof,
    double** loc_mat, double* loc_rhs, typename
    Template_names<3>::SquareMatrixD* global_matrix, double* global_rhs);
#else
template void copy_local_to_global<2>(const int& n_dofs, const int* global_dof,
    double** loc_mat, double* loc_rhs, typename
    Template_names<2>::SquareMatrixD* global_matrix, double* global_rhs);
#endif


template < int d >
void compute_basis_fct_values_at_joint(const FiniteElement& finite_element,
    const TBaseCell* cell, const int& joint_index, const int& base_vec_dim,
    const TQuadFormula& quad_formula,
    std::vector< std::vector<double> >& original_values,
    std::vector< std::vector<double> >& original_values_dx,
    std::vector< std::vector<double> >& original_values_dy,
    std::vector< std::vector<double> >& original_values_dz,
    const TCollection& coll
    )
{
  // Get some basic information from the finite element on the cell for
  // transformation later
  auto n_base_func = finite_element.GetN_DOF();

  // resize output values to correct size
  auto n_quad_pts = quad_formula.GetN_QuadPoints();
  original_values.resize(n_quad_pts);
  original_values_dx.resize(n_quad_pts);
  original_values_dy.resize(n_quad_pts);
  if (d == 3)
  {
    original_values_dz.resize(n_quad_pts);
  }

  // Prepare cell for transformation of basis functions from reference cell to
  // original cell

  // Isoparametric elements are not supported yet. We need to find the reference
  // transformation and if the element was an isoparametric element then the
  // reference transformation differs from the standard ones.
  if (TDatabase::ParamDB->USE_ISOPARAMETRIC == 1)
  {
    ErrThrow("This function can only handle non-isoparametric elements.");
  }
  // Checking whether an element is a isoparametric element could look like the
  // following lines that are commented out
  // bool IsIsoparametric = false;
  // if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
  // {
  //   auto n_joints = cell->GetN_Joints();
  //   for(int joint_i = 0; joint_i < n_joints; joint_i++)
  //   {
  //     auto joint = cell->GetJoint(joint_i);
  //     auto jointtype = joint->GetType();
  //     if(jointtype == BoundaryEdge)
  //     {
  //       auto bdtype = ((const TBoundEdge *)(joint))->GetBoundComp()->GetType();
  //       if(bdtype != Line)
  //         IsIsoparametric = true;
  //     }
  //     if(jointtype == InterfaceJoint)
  //     {
  //       auto bdtype = ((const TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
  //       if(bdtype != Line)
  //         IsIsoparametric = true;
  //     }
  //     if(jointtype == IsoInterfaceJoint ||
  //         jointtype == IsoBoundEdge)
  //       IsIsoparametric = true;
  //   }
  // }// endif
  auto ref_trans_id = finite_element.GetRefTransID();
  // auto ref_element = finite_element.GetBaseFunct()->GetRefElement();
  // if(IsIsoparametric)
  // {
  //   switch(ref_element)
  //   {
  //     case BFRefElements::BFUnitSquare:
  //       ref_trans_id = ReferenceTransformation_type::QuadIsoparametric;
  //       break;

  //     case BFRefElements::BFUnitTriangle:
  //       ref_trans_id = ReferenceTransformation_type::TriaIsoparametric;
  //       break;

  //     case BFRefElements::BFUnitTetrahedron:
  //       ref_trans_id = ReferenceTransformation_type::TetraIsoparametric;
  //       break;

  //     case BFRefElements::BFUnitHexahedron:
  //       ref_trans_id = ReferenceTransformation_type::HexaIsoparametric;
  //       break;
  //   }
  // } // endif IsIsoparametric
#ifdef __2D__
  auto ref_trans = FEDatabase::GetRefTrans2D(ref_trans_id);
#else
  auto ref_trans = FEDatabase::GetRefTrans3D(ref_trans_id);
#endif
  // auto degree = finite_element.GetBaseFunct()->GetPolynomialDegree();
  // auto qf_id = quad_formula.get_type();
  // switch(ref_trans_id)
  // {
  // #ifdef __2D__
  //   case ReferenceTransformation_type::TriaAffin:
  //   case ReferenceTransformation_type::QuadAffin:
  //   case ReferenceTransformation_type::QuadBilinear:
  //     break;
  //   case ReferenceTransformation_type::TriaIsoparametric:
  //     ((TTriaIsoparametric *)ref_trans)->SetApproximationOrder(degree);
  //     ((TTriaIsoparametric *)ref_trans)->SetQuadFormula(qf_id);
  //     break;
  //   case ReferenceTransformation_type::QuadIsoparametric:
  //     ((TQuadIsoparametric *)ref_trans)->SetApproximationOrder(degree);
  //     ((TQuadIsoparametric *)ref_trans)->SetQuadFormula(qf_id);
  //     break;
  // #else
  //   case ReferenceTransformation_type::TetraAffin:
  //   case ReferenceTransformation_type::HexaAffin:
  //   case ReferenceTransformation_type::HexaTrilinear:
  //     break;
  //   case ReferenceTransformation_type::TetraIsoparametric:
  //     ((TTetraIsoparametric *)ref_trans)->SetApproximationOrder(degree);
  //     ((TTetraIsoparametric *)ref_trans)->SetQuadFormula(qf_id);
  //     break;
  //   case ReferenceTransformation_type::HexaIsoparametric:
  //     ((THexaIsoparametric *)ref_trans)->SetApproximationOrder(degree);
  //     ((THexaIsoparametric *)ref_trans)->SetQuadFormula(qf_id);
  //     break;
  // #endif
  // } // endswitch
  ref_trans->SetCell(cell);

  // Compute the values of the basis functions and their derivatives in
  // quadrature points according to 2D and 3D methods of ParMooN
#ifdef __2D__
  // Get values of basis functions and their derivatives in all quadrature
  // points on REFERENCE cell
  auto ref_values = FEDatabase::GetJointDerivatives2D(
      *finite_element.GetBaseFunct(), quad_formula, joint_index,
      MultiIndex2D::D00);
  auto ref_values_dxi = FEDatabase::GetJointDerivatives2D(
      *finite_element.GetBaseFunct(), quad_formula, joint_index,
      MultiIndex2D::D10);
  auto ref_values_deta = FEDatabase::GetJointDerivatives2D(
      *finite_element.GetBaseFunct(), quad_formula, joint_index,
      MultiIndex2D::D01);

  // Transform values of basis functions and their derivatives for each
  // quadrature point to ORIGINAL cell

  std::vector<double> u_orig(n_base_func * base_vec_dim);
  std::vector<double> u_orig_dx(n_base_func * base_vec_dim);
  std::vector<double> u_orig_dy(n_base_func * base_vec_dim);

  for( int quad_pt_i = 0; quad_pt_i < n_quad_pts; quad_pt_i++ )
  {
    auto quad_pt_x = quad_formula.get_point(quad_pt_i).x;
    ref_trans->GetOrigValues(joint_index, quad_pt_x, n_base_func,
        ref_values[quad_pt_i], ref_values_dxi[quad_pt_i],
        ref_values_deta[quad_pt_i], u_orig.data(), u_orig_dx.data(),
        u_orig_dy.data(), base_vec_dim);

    finite_element.GetBaseFunct()->ChangeBF(&coll, cell, u_orig.data());
    finite_element.GetBaseFunct()->ChangeBF(&coll, cell, u_orig_dx.data());
    finite_element.GetBaseFunct()->ChangeBF(&coll, cell, u_orig_dy.data());

    original_values[quad_pt_i].resize(n_base_func * base_vec_dim);
    original_values_dx[quad_pt_i].resize(n_base_func * base_vec_dim);
    original_values_dy[quad_pt_i].resize(n_base_func * base_vec_dim);

    for (int base_func_j = 0; base_func_j < n_base_func * base_vec_dim; base_func_j++)

    {
      original_values[quad_pt_i][base_func_j] = u_orig[base_func_j];
      original_values_dx[quad_pt_i][base_func_j] = u_orig_dx[base_func_j];
      original_values_dy[quad_pt_i][base_func_j] = u_orig_dy[base_func_j];
    }

    // for (int i = 0; i < base_vec_dim; i++)
    // {
    //   for(int j=0;j < n_base_func; j++)
    //   {
    //     int edge = finite_element.GetFEDesc()->GetJointOfThisDOF(j);
    //     if(edge!=-1) // edge ==-1 means inner dof
    //     {
    //       original_values[quad_pt_i][j+i*n_base_func] *= cell->GetNormalOrientation(edge);
    //       original_values_dx[quad_pt_i][j+i*n_base_func] *= cell->GetNormalOrientation(edge);
    //       original_values_dy[quad_pt_i][j+i*n_base_func] *= cell->GetNormalOrientation(edge);
    //     }
    //   }
    // }


  }
#else
  // Get values of basis functions and their derivatives in all quadrature
  // points on REFERENCE cell
  auto ref_values = FEDatabase::GetJointDerivatives3D(
      *finite_element.GetBaseFunct(), quad_formula, joint_index,
      MultiIndex3D::D000);
  auto ref_values_dxi = FEDatabase::GetJointDerivatives3D(
      *finite_element.GetBaseFunct(), quad_formula, joint_index,
      MultiIndex3D::D100);
  auto ref_values_deta = FEDatabase::GetJointDerivatives3D(
      *finite_element.GetBaseFunct(), quad_formula, joint_index,
      MultiIndex3D::D010);
  auto ref_values_dzeta = FEDatabase::GetJointDerivatives3D(
      *finite_element.GetBaseFunct(), quad_formula, joint_index,
      MultiIndex3D::D001);

  // Transform values of basis functions and their derivatives for each
  // quadrature point to ORIGINAL cell
  auto ref_element = finite_element.GetBaseFunct()->GetRefElement();
  std::vector<double> u_orig(n_base_func * base_vec_dim);
  std::vector<double> u_orig_dx(n_base_func * base_vec_dim);
  std::vector<double> u_orig_dy(n_base_func * base_vec_dim);
  std::vector<double> u_orig_dz(n_base_func * base_vec_dim);
  for( int quad_pt_i = 0; quad_pt_i < n_quad_pts; quad_pt_i++ )
  {
    auto quad_pt_x = quad_formula.get_point(quad_pt_i).x;
    auto quad_pt_y = quad_formula.get_point(quad_pt_i).y;

    auto ref_point_3D = transform(ref_element, joint_index, quad_pt_x,
        quad_pt_y); auto xi = ref_point_3D.x;
    auto eta = ref_point_3D.y;
    auto zeta = ref_point_3D.z;
    ref_trans->GetOrigValues(xi, eta, zeta, n_base_func, ref_values[quad_pt_i],
        ref_values_dxi[quad_pt_i], ref_values_deta[quad_pt_i],
        ref_values_dzeta[quad_pt_i], u_orig.data(), u_orig_dx.data(),
        u_orig_dy.data(), u_orig_dz.data(), base_vec_dim);

    original_values[quad_pt_i].resize(n_base_func * base_vec_dim);
    original_values_dx[quad_pt_i].resize(n_base_func * base_vec_dim);
    original_values_dy[quad_pt_i].resize(n_base_func * base_vec_dim);
    original_values_dz[quad_pt_i].resize(n_base_func * base_vec_dim);
    for (int base_func_j = 0; base_func_j < n_base_func * base_vec_dim; base_func_j++)
    {
      original_values[quad_pt_i][base_func_j] = u_orig[base_func_j];
      original_values_dx[quad_pt_i][base_func_j] = u_orig_dx[base_func_j];
      original_values_dy[quad_pt_i][base_func_j] = u_orig_dy[base_func_j];
      original_values_dz[quad_pt_i][base_func_j] = u_orig_dz[base_func_j];
    }
  }
#endif
}
#ifdef __3D__
template void compute_basis_fct_values_at_joint<3>(const FiniteElement&
    finite_element, const TBaseCell* cell, const int& joint_index, const int&
    base_vec_dim, const TQuadFormula& quad_formula,
    std::vector< std::vector<double> >& original_values,
    std::vector< std::vector<double> >& original_values_dx,
    std::vector< std::vector<double> >& original_values_dy,
    std::vector< std::vector<double> >& original_values_dz,
    const TCollection& coll
    );
#else
template void compute_basis_fct_values_at_joint<2>(const FiniteElement&
    finite_element, const TBaseCell* cell, const int& joint_index, const int&
    base_vec_dim, const TQuadFormula& quad_formula,
    std::vector< std::vector<double> >& original_values,
    std::vector< std::vector<double> >& original_values_dx,
    std::vector< std::vector<double> >& original_values_dy,
    std::vector< std::vector<double> >& original_values_dz,
    const TCollection& coll
    );
#endif

std::vector<parmoon::Point> transform_quad_points( const TBaseCell* cell, const
    BFRefElements& ref_element, const ReferenceTransformation_type&
    ref_trans_id, const int& joint_index, const TQuadFormula* quad_form )
{
  auto n_quad_pts = quad_form->GetN_QuadPoints();

  // xi, eta and zeta are the d dimensional points on the reference cell.
  // Those points are computed below by the transformation of the (d-1)
  // dimensional quadrature points to the d dimensional reference cell.
  std::vector<double> xi(n_quad_pts);
  std::vector<double> eta(n_quad_pts);
#ifdef __3D__
  std::vector<double> zeta(n_quad_pts);
#endif

  for(int quad_pt_i = 0; quad_pt_i < n_quad_pts; quad_pt_i++)
  {
#ifdef __2D__
    // s is the 1 d coordinate of the quadrature point
    auto s = quad_form->get_point(quad_pt_i).x;  // 1D coord of point
    // xi_eta_2D is s transformed to the reference cell, i.e. 2 dimensional
    auto xi_eta_2D = transform(ref_element, joint_index, s);
    xi[quad_pt_i] = xi_eta_2D.x;
    eta[quad_pt_i] = xi_eta_2D.y;
#else
    // s and t are the 2 d coordinates of the quadrature point
    auto s = quad_form->get_point(quad_pt_i).x;  // 2D coord of point
    auto t = quad_form->get_point(quad_pt_i).y;
    // xi_eta__zeta_3D is (s,t) transformed to the reference cell, i.e.
    // 3 dimensional
    auto xi_eta_zeta_3D = transform(ref_element, joint_index, s, t);
    xi[quad_pt_i] = xi_eta_zeta_3D.x;
    eta[quad_pt_i] = xi_eta_zeta_3D.y;
    zeta[quad_pt_i] = xi_eta_zeta_3D.z;
#endif
  } // endfor quad_pt_i

  // x, y and z are the points on the original cell, i.e. the points after
  // transformation of xi, eta and zeta
  std::vector<double> x(n_quad_pts);
  std::vector<double> y(n_quad_pts);
  std::vector<double> z(n_quad_pts);

  FEDatabase::SetCellForRefTrans(cell, ref_trans_id);
#ifdef __2D__
  FEDatabase::GetOrigFromRef(ref_trans_id, n_quad_pts,
      xi.data(), eta.data(), x.data(), y.data());
#else
  FEDatabase::GetOrigFromRef(ref_trans_id, n_quad_pts,
      xi.data(), eta.data(), zeta.data(), x.data(), y.data(), z.data());
#endif

#ifdef __2D__
  int d = 2;
#else
  int d = 3;
#endif

  std::vector<parmoon::Point> transformed_quad_points(n_quad_pts,
      parmoon::Point((unsigned int) d-1));
  for (int quad_pt_i = 0; quad_pt_i < n_quad_pts; ++quad_pt_i)
  {
    if (d == 2)
    {
      transformed_quad_points[quad_pt_i] = parmoon::Point( x[quad_pt_i],
          y[quad_pt_i] ); }
    else
    {
      transformed_quad_points[quad_pt_i] = parmoon::Point( x[quad_pt_i],
          y[quad_pt_i], z[quad_pt_i] );
    }
  }

  return transformed_quad_points;
}

std::vector<double> compute_orig_weights( TBaseCell* cell, const int&
    joint_index, const TQuadFormula* quad_form)
{
  auto n_quad_pts = quad_form->GetN_QuadPoints();
  std::vector<double> orig_weights(n_quad_pts);
  auto shape = cell->GetType();
  // These shape types belong to 2D shapes, see also Shapes.
  if (shape >= 1 && shape <= 4)
  {
    // The facet in this case is a line. The weights have to be transformed to
    // weight[i] = reference_weight[i] * h_E / 2
    // what can be seen from transforming the reference line (-1, 1) to any
    // other line (a, b)
    auto correction_factor = cell->ComputeDiameterOfJoint(joint_index) / 2;
    for (int pt_i = 0; pt_i < n_quad_pts; ++pt_i)
    {
      orig_weights[pt_i] = quad_form->get_weight(pt_i) * correction_factor;
    }
  }
  // These shape types belong to 3D shapes, see also Shapes.
  else if (shape >= 5 && shape <= 7)
  {
    // In 3D a facet is a triangle or a quadrilateral, and hence, the
    // transformation may be non-affine. In this case the correction factor is
    // non constant. In any case the correction factor depends on the shape of
    // the cell and therefore, we need some information about the shape which is
    // encoded in the shape descriptor.

    // Get information about the vertices that build the face with number
    // joint_indext. This information is precisely given by FaceVertex, see also
    // TShapeDesc.
    auto shape_desc = cell->GetShapeDesc();
    const int* face_vertex;
    const int* face_vertex_len;
    int max_face_vertex_len;
    shape_desc->GetFaceVertex(face_vertex, face_vertex_len,
        max_face_vertex_len);

    // Find the position in FaceVertex where the vertices of the joint_index-th
    // face start. The ordering in FaceVertex is such that first the vertices of
    // the 0th face are specified, then of the 1st face and so on. Therefore,
    // the starting point for our current face is the sum of all vertices that
    // belong to faces with a smaller local face number.
    int starting_index = 0;
    for (int facet_nr = 0; facet_nr < joint_index; ++facet_nr)
    {
      starting_index += face_vertex_len[facet_nr];
    }

    // Get the coordinates of three vertices of the face. I try to use the same
    // notation as in the beginning of this if block. To remind you, V_i is not
    // the ith vertex in the cell but the i-th vertex in the face.
    std::vector<std::array<double, 3>> vertices(face_vertex_len[joint_index]);
    for (int vertex_i = 0; vertex_i < face_vertex_len[joint_index]; ++vertex_i)
    {
      cell->GetVertex(face_vertex[starting_index + vertex_i])
        ->GetCoords(vertices[vertex_i][0], vertices[vertex_i][1],
            vertices[vertex_i][2]);
    }

    if (shape == 5)
    {
      // triangular case and therefore an affine transform
      auto xc1 = vertices[1][0] - vertices[0][0];
      auto xc2 = vertices[2][0] - vertices[0][0];

      auto yc1 = vertices[1][1] - vertices[0][1];
      auto yc2 = vertices[2][1] - vertices[0][1];

      auto zc1 = vertices[1][2] - vertices[0][2];
      auto zc2 = vertices[2][2] - vertices[0][2];

      // normal vector
      auto nx = yc1*zc2 - zc1*yc2;
      auto ny = zc1*xc2 - xc1*zc2;
      auto nz = xc1*yc2 - yc1*xc2;
      // determinant of reference trafo
      auto correction_factor = std::sqrt(nx*nx + ny*ny + nz*nz);
      for (int pt_i = 0; pt_i < n_quad_pts; ++pt_i)
      {
        orig_weights[pt_i] = quad_form->get_weight(pt_i) * correction_factor;
      }
    }
    else
    {
      // quadrilateral case and therefore possibly a bilinear transformation
      auto xc1 = (-vertices[0][0] + vertices[1][0] + vertices[2][0]
          - vertices[3][0]) * 0.25;
      auto xc2 = (-vertices[0][0] - vertices[1][0] + vertices[2][0]
          + vertices[3][0]) * 0.25;
      auto xc3 = ( vertices[0][0] - vertices[1][0] + vertices[2][0]
          - vertices[3][0]) * 0.25;

      auto yc1 = (-vertices[0][1] + vertices[1][1] + vertices[2][1]
          - vertices[3][1]) * 0.25;
      auto yc2 = (-vertices[0][1] - vertices[1][1] + vertices[2][1]
          + vertices[3][1]) * 0.25;
      auto yc3 = ( vertices[0][1] - vertices[1][1] + vertices[2][1]
          - vertices[3][1]) * 0.25;

      auto zc1 = (-vertices[0][2] + vertices[1][2] + vertices[2][2]
          - vertices[3][2]) * 0.25;
      auto zc2 = (-vertices[0][2] - vertices[1][2] + vertices[2][2]
          + vertices[3][2]) * 0.25;
      auto zc3 = ( vertices[0][2] - vertices[1][2] + vertices[2][2]
          - vertices[3][2]) * 0.25;

      for (int pt_i = 0; pt_i < n_quad_pts; ++pt_i)
      {
        auto quad_point = quad_form->get_point(pt_i);
        auto nx = (yc1 + quad_point.y * yc3)*(zc2 + quad_point.x * zc3)
          - (zc1 + quad_point.y * zc3)*(yc2 + quad_point.x * yc3);
        auto ny = (zc1 + quad_point.y * zc3)*(xc2 + quad_point.x * xc3)
          - (xc1 + quad_point.y * xc3)*(zc2 + quad_point.x * zc3);
        auto nz = (xc1 + quad_point.y * xc3)*(yc2 + quad_point.x * yc3)
          - (yc1 + quad_point.y * yc3)*(xc2 + quad_point.x * xc3);
        auto correction_factor = std::sqrt(nx*nx+ny*ny+nz*nz);
        orig_weights[pt_i] = quad_form->get_weight(pt_i) * correction_factor;
      }
    }
  }
  else
  {
    ErrThrow("1D shapes are not implemented yet.");
  }
  return orig_weights;;
}


// Computes the local DG bilinear form for convection-diffusion problems
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
    loc_mat, int matrix_flag)
{

  // get some information needed for integration
  auto n_quad_pts = quad_form_joint->GetN_QuadPoints();
  auto n_dofs = val_bas_fct[0].size()/base_vec_dim;
  auto n_dofs_neigh = (is_boundary_edge) ? 0 : val_bas_fct_neigh[0].size()/base_vec_dim;
  auto n_dofs_total = n_dofs + n_dofs_neigh;
  auto nx = joint_normal[0];
  auto ny = joint_normal[1];
  auto nz = (d == 3) ? joint_normal[2] : 0.0;

  // get parameters of DG method
  auto face_sigma = param_db["face_sigma_DG"].get<double>(); // edge penalty
  auto kappa = param_db["symmetry_DG"].get<int>();  // symmetry
  auto eta = param_db["eta_upwind_DG"].get<double>(); // upwind

  if (is_boundary_edge)
  { // adapt boundary penalty, see e.g. Kanschat 2008
    face_sigma *= 2;
  }

  for (int quad_pt_i = 0; quad_pt_i < n_quad_pts; ++quad_pt_i)
  {
    // parameter of PDE
    auto beta_1 = coefficients[quad_pt_i][5];
    auto beta_2 = coefficients[quad_pt_i][6];
    auto convection_3 = (d == 3) ? coefficients[quad_pt_i][3] : 0.0;

    for (unsigned int t_dof = 0; t_dof < n_dofs_total; ++t_dof)
    { // for all test dofs

      // Define jump and averages for test functions index [0] corresponds to
      // the basis functions itself index [1/2/3] corresponds to the derivative
      // wrt. x/y/z for the basis functions.
      // No matter if we are in 2D or 3D we always assume 4 values, i.e. the
      // 4 indices mentioned above, and in 2D we set every entry corresponding
      // to z to 0.
      std::vector<double> t_jump(4 * base_vec_dim);
      std::vector<double> t_average(4 * base_vec_dim);

      // distinguish whether test function are related to the cell or to the
      // neighbor. Compute jump and average  accordingly using the fact that the
      // support of a basis function is exactly the cell it lives on.
      if (t_dof < n_dofs)
      { // test dof from cell
        t_jump[0] = val_bas_fct [quad_pt_i][t_dof];
        t_jump[1] = val_bas_fct_dx [quad_pt_i][t_dof];
        t_jump[2] = val_bas_fct_dy [quad_pt_i][t_dof];
        t_jump[3] = (d==3) ? val_bas_fct_dz [quad_pt_i][t_dof] : 0;        
      }
      else
      { // test dof from neighbor
        t_jump[0] = -val_bas_fct_neigh[index_map[quad_pt_i]][t_dof - n_dofs];
        t_jump[1] = -val_bas_fct_neigh_dx[index_map[quad_pt_i]][t_dof - n_dofs];
        t_jump[2] = -val_bas_fct_neigh_dy[index_map[quad_pt_i]][t_dof - n_dofs];
        t_jump[3] = (d == 3) ? -val_bas_fct_neigh_dz[index_map[quad_pt_i]][t_dof - n_dofs] : 0;
      }
      if (is_boundary_edge)
      { // Correct values for boundary edges: jump and average on the boundary
        // are defined as values seen from cell
        for (int i = 0; i < 4; ++i)
        {
          t_average[i] *= 2;
        }
      }
      for (unsigned int a_dof = 0; a_dof < n_dofs_total; ++a_dof)
      { // for all ansatz dofs
        // Define jump and averages for ansatz functions index [0] corresponds
        // to the basis functions itself index [1/2/3] corresponds to the
        // derivative wrt. x/y/z for the basis functions.  No matter if we are
        // in 2D or 3D we always assume 4 values, i.e. the
        // 4 indices mentioned above, and in 2D we set every entry corresponding
        // to z to 0.
        std::vector<double> a_jump(4 * base_vec_dim);
        std::vector<double> a_average(4 * base_vec_dim);

        // distinguish whether ansatz function are related to the cell or to the
        // neighbor. Compute jump and average  accordingly using the fact that
        // the support of a basis function is exactly the cell it lives on.
        if (a_dof < n_dofs)
        { // ansatz dof from cell
          a_jump[0] = val_bas_fct [quad_pt_i][a_dof];
          a_jump[1] = val_bas_fct_dx [quad_pt_i][a_dof];
          a_jump[2] = val_bas_fct_dy [quad_pt_i][a_dof];
          a_jump[3] = (d == 3) ? val_bas_fct_dz [quad_pt_i][a_dof] : 0;          
        }
        else
        { // ansatz dof from neighbor
          a_jump[0] = -val_bas_fct_neigh [index_map[quad_pt_i]][a_dof - n_dofs];
          a_jump[1] = -val_bas_fct_neigh_dx [index_map[quad_pt_i]][a_dof - n_dofs];
          a_jump[2] = -val_bas_fct_neigh_dy [index_map[quad_pt_i]][a_dof - n_dofs];
          a_jump[3] = (d == 3) ? -val_bas_fct_neigh_dz[index_map[quad_pt_i]][a_dof - n_dofs] : 0;          
        }
        if (is_boundary_edge)
        { // Correct values for boundary edges: jump and average on the boundary
          // are defined as values seen from cell
          for (int i = 0; i < 4; ++i)
          {
            a_average[i] *= 2;
          }
        }

        double integrand = 0;
        
        if (is_boundary_edge)
        {
          ErrThrow("Only the interior is implemented");
        }
        else
        {
          if(matrix_flag == 1)
          {
            integrand += 2 * joint_diam * joint_diam * 
                       (beta_1 * ny * a_jump[1] + beta_2 * ny * a_jump[2])
                     * (beta_1 * ny * t_jump[1] + beta_2 * ny * t_jump[2]);
          }
          if(matrix_flag == 3)
          {
            integrand -= 2 * joint_diam * joint_diam * 
                       (beta_1 * nx * a_jump[1] + beta_2 * nx * a_jump[2])
                     * (beta_1 * nx * t_jump[1] + beta_2 * nx * t_jump[2]);
          }
          
        }

        auto weight = orig_weights[quad_pt_i];
        // Update local matrix
        loc_mat[t_dof][a_dof] += weight * integrand;
      } // endfor ansatz dofs
    } // endfor test dofs
  } //endfor quad points
} // end compute_DG_bilinear_form_local
#ifdef __3D__
template void compute_DG_bilinear_form_local_Oseen<3>(const ParameterDatabase&
    param_db, const int& base_vec_dim, const std::vector<double>* val_bas_fct, const std::vector<double>*
    val_bas_fct_dx, const std::vector<double>* val_bas_fct_dy, const
    std::vector<double>* val_bas_fct_dz, const std::vector<double>*
    val_bas_fct_neigh, const std::vector<double>* val_bas_fct_neigh_dx, const
    std::vector<double>* val_bas_fct_neigh_dy, const std::vector<double>*
    val_bas_fct_neigh_dz, const std::vector<int>& index_map, const double*
    joint_normal, const double& joint_diam, const TQuadFormula* quad_form_joint,
    const std::vector<double>& orig_weights, double** coefficients, const bool&
    is_boundary_edge, double** loc_mat, int matrix_flag);
#else
template void compute_DG_bilinear_form_local_Oseen<2>(const ParameterDatabase&
    param_db, const int& base_vec_dim, const std::vector<double>* val_bas_fct, const std::vector<double>*
    val_bas_fct_dx, const std::vector<double>* val_bas_fct_dy, const
    std::vector<double>* val_bas_fct_dz, const std::vector<double>*
    val_bas_fct_neigh, const std::vector<double>* val_bas_fct_neigh_dx, const
    std::vector<double>* val_bas_fct_neigh_dy, const std::vector<double>*
    val_bas_fct_neigh_dz, const std::vector<int>& index_map, const double*
    joint_normal, const double& joint_diam, const TQuadFormula* quad_form_joint,
    const std::vector<double>& orig_weights, double** coefficients, const bool&
    is_boundary_edge, double** loc_mat, int matrix_flag);
#endif


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
    )
{
  // Get some basic information from the finite element on the cell for
  // transformation later
  auto n_base_func = finite_element.GetN_DOF();

  // resize output values to correct size
  auto n_quad_pts = quad_formula.GetN_QuadPoints();
  u.resize(n_quad_pts);
  ux.resize(n_quad_pts);
  uy.resize(n_quad_pts);
  uxx.resize(n_quad_pts);
  uxy.resize(n_quad_pts);
  uyy.resize(n_quad_pts);
  
  if (d == 3)
  {
    uz.resize(n_quad_pts);
  }

  // Prepare cell for transformation of basis functions from reference cell to
  // original cell

  // Isoparametric elements are not supported yet. We need to find the reference
  // transformation and if the element was an isoparametric element then the
  // reference transformation differs from the standard ones.
  if (TDatabase::ParamDB->USE_ISOPARAMETRIC == 1)
  {
    ErrThrow("This function can only handle non-isoparametric elements.");
  }
  auto ref_trans_id = finite_element.GetRefTransID();
#ifdef __2D__
  auto ref_trans = FEDatabase::GetRefTrans2D(ref_trans_id);
#else
  auto ref_trans = FEDatabase::GetRefTrans3D(ref_trans_id);
#endif

  ref_trans->SetCell(cell);

  // Compute the values of the basis functions and their derivatives in
  // quadrature points according to 2D and 3D methods of ParMooN
#ifdef __2D__
  // Get values of basis functions and their derivatives in all quadrature
  // points on REFERENCE cell
  auto r00 = FEDatabase::GetJointDerivatives2D(
      *finite_element.GetBaseFunct(), quad_formula, joint_index,
      MultiIndex2D::D00);
  auto r10 = FEDatabase::GetJointDerivatives2D(
      *finite_element.GetBaseFunct(), quad_formula, joint_index,
      MultiIndex2D::D10);
  auto r01 = FEDatabase::GetJointDerivatives2D(
      *finite_element.GetBaseFunct(), quad_formula, joint_index,
      MultiIndex2D::D01);
  auto r20 = FEDatabase::GetJointDerivatives2D(
      *finite_element.GetBaseFunct(), quad_formula, joint_index,
      MultiIndex2D::D20);
  auto r11 = FEDatabase::GetJointDerivatives2D(
      *finite_element.GetBaseFunct(), quad_formula, joint_index,
      MultiIndex2D::D11);
  auto r02 = FEDatabase::GetJointDerivatives2D(
      *finite_element.GetBaseFunct(), quad_formula, joint_index,
      MultiIndex2D::D02);
  // Transform values of basis functions and their derivatives for each
  // quadrature point to ORIGINAL cell

  std::vector<double> o00(n_base_func * base_vec_dim);
  std::vector<double> o10(n_base_func * base_vec_dim);
  std::vector<double> o01(n_base_func * base_vec_dim);
  std::vector<double> o20(n_base_func * base_vec_dim);
  std::vector<double> o11(n_base_func * base_vec_dim);
  std::vector<double> o02(n_base_func * base_vec_dim);

  for( int quad_pt_i = 0; quad_pt_i < n_quad_pts; quad_pt_i++ )
  {
    auto quad_pt_x = quad_formula.get_point(quad_pt_i).x;
    auto quad_pt_y = quad_formula.get_point(quad_pt_i).y;
    ref_trans->GetOrigAllDerivatives(quad_pt_x, quad_pt_y, n_base_func,
      r00[quad_pt_i], r10[quad_pt_i], r01[quad_pt_i],
      r20[quad_pt_i], r11[quad_pt_i], r02[quad_pt_i],
      o00.data(), o10.data(), o01.data(), o20.data(), 
      o11.data(), o02.data(), base_vec_dim);

    finite_element.GetBaseFunct()->ChangeBF(&coll, cell, o00.data());
    finite_element.GetBaseFunct()->ChangeBF(&coll, cell, o10.data());
    finite_element.GetBaseFunct()->ChangeBF(&coll, cell, o01.data());
    finite_element.GetBaseFunct()->ChangeBF(&coll, cell, o20.data());
    finite_element.GetBaseFunct()->ChangeBF(&coll, cell, o11.data());
    finite_element.GetBaseFunct()->ChangeBF(&coll, cell, o20.data());

    u[quad_pt_i].resize(n_base_func * base_vec_dim);
    ux[quad_pt_i].resize(n_base_func * base_vec_dim);
    uy[quad_pt_i].resize(n_base_func * base_vec_dim);
    uxx[quad_pt_i].resize(n_base_func * base_vec_dim);
    uxy[quad_pt_i].resize(n_base_func * base_vec_dim);
    uyy[quad_pt_i].resize(n_base_func * base_vec_dim);

    for (int base_func_j = 0; base_func_j < n_base_func * base_vec_dim; base_func_j++)

    {
      u[quad_pt_i][base_func_j] = o00[base_func_j];
      ux[quad_pt_i][base_func_j] = o10[base_func_j];
      uy[quad_pt_i][base_func_j] = o01[base_func_j];
      uxx[quad_pt_i][base_func_j] = o20[base_func_j];
      uxy[quad_pt_i][base_func_j] = o11[base_func_j];
      uyy[quad_pt_i][base_func_j] = o02[base_func_j];
    }    
  }
#endif
}
#ifdef __3D__
template void compute_basis_fct_values_at_joint<3>(const FiniteElement& finite_element,
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
#else
template void compute_basis_fct_values_at_joint<2>(const FiniteElement& finite_element,
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
#endif


#ifdef __2D__
void compute_DG_bilinear_form_local_Oseen(const ParameterDatabase& param_db, const int& base_vec_dim, const std::vector<double>* val_bas_fct, const std::vector<double>* val_bas_fct_dx, const std::vector<double>* val_bas_fct_dy, const std::vector<double>* val_bas_fct_dxx, const std::vector<double>* val_bas_fct_dxy, const std::vector<double>* val_bas_fct_dyy, const std::vector<double>* val_bas_fct_dz, const std::vector<double>* val_bas_fct_neigh, const std::vector<double>* val_bas_fct_neigh_dx, const std::vector<double>* val_bas_fct_neigh_dy, const std::vector<double>* val_bas_fct_neigh_dxx, const std::vector<double>* val_bas_fct_neigh_dxy, const std::vector<double>* val_bas_fct_neigh_dyy, const std::vector<double>* val_bas_fct_neigh_dz, const std::vector<int>& index_map, const double* joint_normal, const double& joint_diam, const TQuadFormula* quad_form_joint, const std::vector<double>& orig_weights, double ** coefficients, const bool& is_boundary_edge, double ** loc_mat, int matrix_flag)
{
  // get some information needed for integration
  auto n_quad_pts = quad_form_joint->GetN_QuadPoints();
  auto n_dofs = val_bas_fct[0].size()/base_vec_dim;
  auto n_dofs_neigh = (is_boundary_edge) ? 0 : val_bas_fct_neigh[0].size()/base_vec_dim;
  auto n_dofs_total = n_dofs + n_dofs_neigh;
  auto nx = joint_normal[0];
  auto ny = joint_normal[1];  

  // get parameters of DG method
  auto face_sigma = param_db["face_sigma_DG"].get<double>(); // edge penalty
  auto kappa = param_db["symmetry_DG"].get<int>();  // symmetry
  auto eta = param_db["eta_upwind_DG"].get<double>(); // upwind
  for (int quad_pt_i = 0; quad_pt_i < n_quad_pts; ++quad_pt_i)
  {
    // parameter of PDE
    auto beta_1 = coefficients[quad_pt_i][5];
    auto beta_2 = coefficients[quad_pt_i][6];
    auto beta_1x = coefficients[quad_pt_i][8];
    auto beta_2x = coefficients[quad_pt_i][9];
    auto beta_1y = coefficients[quad_pt_i][10];
    auto beta_2y = coefficients[quad_pt_i][11];
    
    auto beta_1xy = coefficients[quad_pt_i][12];
    auto beta_1yx = coefficients[quad_pt_i][13];
    auto beta_1xx = coefficients[quad_pt_i][14];
    auto beta_1yy = coefficients[quad_pt_i][15];
    
    auto beta_2xy = coefficients[quad_pt_i][16];
    auto beta_2yx = coefficients[quad_pt_i][17];
    auto beta_2xx = coefficients[quad_pt_i][18];
    auto beta_2yy = coefficients[quad_pt_i][19];

    for (unsigned int t_dof = 0; t_dof < n_dofs_total; ++t_dof)
    { // for all test dofs

      // Define jump and averages for test functions index [0] corresponds to
      // the basis functions itself index [1/2/3] corresponds to the derivative
      // wrt. x/y/z for the basis functions.
      // No matter if we are in 2D or 3D we always assume 4 values, i.e. the
      // 4 indices mentioned above, and in 2D we set every entry corresponding
      // to z to 0.
      std::vector<double> t_jump(6 * base_vec_dim);
      std::vector<double> t_average(4 * base_vec_dim);

      // distinguish whether test function are related to the cell or to the
      // neighbor. Compute jump and average  accordingly using the fact that the
      // support of a basis function is exactly the cell it lives on.
      if (t_dof < n_dofs)
      { // test dof from cell
        t_jump[0] = val_bas_fct [quad_pt_i][t_dof];
        t_jump[1] = val_bas_fct_dx [quad_pt_i][t_dof];
        t_jump[2] = val_bas_fct_dy [quad_pt_i][t_dof];
        t_jump[3] = val_bas_fct_dxx[quad_pt_i][t_dof];
        t_jump[4] = val_bas_fct_dxy[quad_pt_i][t_dof];
        t_jump[5] = val_bas_fct_dyy[quad_pt_i][t_dof];
      }
      else
      { // test dof from neighbor
        t_jump[0] = -val_bas_fct_neigh[index_map[quad_pt_i]][t_dof - n_dofs];
        t_jump[1] = -val_bas_fct_neigh_dx[index_map[quad_pt_i]][t_dof - n_dofs];
        t_jump[2] = -val_bas_fct_neigh_dy[index_map[quad_pt_i]][t_dof - n_dofs];
        t_jump[3] = -val_bas_fct_neigh_dxx[index_map[quad_pt_i]][t_dof - n_dofs];
        t_jump[4] = -val_bas_fct_neigh_dxy[index_map[quad_pt_i]][t_dof - n_dofs];
        t_jump[5] = -val_bas_fct_neigh_dyy[index_map[quad_pt_i]][t_dof - n_dofs];
      }
      if (is_boundary_edge)
      { // Correct values for boundary edges: jump and average on the boundary
        // are defined as values seen from cell
        for (int i = 0; i < 4; ++i)
        {
          t_average[i] *= 2;
        }
      }
      for (unsigned int a_dof = 0; a_dof < n_dofs_total; ++a_dof)
      { // for all ansatz dofs
        // Define jump and averages for ansatz functions index [0] corresponds
        // to the basis functions itself index [1/2/3] corresponds to the
        // derivative wrt. x/y/z for the basis functions.  No matter if we are
        // in 2D or 3D we always assume 4 values, i.e. the
        // 4 indices mentioned above, and in 2D we set every entry corresponding
        // to z to 0.
        std::vector<double> a_jump(6 * base_vec_dim);
        std::vector<double> a_average(4 * base_vec_dim);

        // distinguish whether ansatz function are related to the cell or to the
        // neighbor. Compute jump and average  accordingly using the fact that
        // the support of a basis function is exactly the cell it lives on.
        if (a_dof < n_dofs)
        { // ansatz dof from cell
          a_jump[0] = val_bas_fct [quad_pt_i][a_dof];
          a_jump[1] = val_bas_fct_dx [quad_pt_i][a_dof];
          a_jump[2] = val_bas_fct_dy [quad_pt_i][a_dof];
          a_jump[3] = val_bas_fct_dxx [quad_pt_i][a_dof];
          a_jump[4] = val_bas_fct_dxy [quad_pt_i][a_dof];
          a_jump[5] = val_bas_fct_dyy [quad_pt_i][a_dof];
        }
        else
        { // ansatz dof from neighbor
          a_jump[0] = -val_bas_fct_neigh [index_map[quad_pt_i]][a_dof - n_dofs];
          a_jump[1] = -val_bas_fct_neigh_dx [index_map[quad_pt_i]][a_dof - n_dofs];
          a_jump[2] = -val_bas_fct_neigh_dy [index_map[quad_pt_i]][a_dof - n_dofs];
          a_jump[3] = -val_bas_fct_neigh_dxx[index_map[quad_pt_i]][a_dof - n_dofs];
          a_jump[4] = -val_bas_fct_neigh_dxy[index_map[quad_pt_i]][a_dof - n_dofs];
          a_jump[5] = -val_bas_fct_neigh_dyy[index_map[quad_pt_i]][a_dof - n_dofs];
        }
        if (is_boundary_edge)
        { // Correct values for boundary edges: jump and average on the boundary
          // are defined as values seen from cell
          for (int i = 0; i < 4; ++i)
          {
            a_average[i] *= 2;
          }
        }

        double integrand = 0;
        
        double factor1 = TDatabase::ParamDB->P1 * std::pow(joint_diam, 2);
        double factor2 = TDatabase::ParamDB->P1 * std::pow(joint_diam, 4);
        double factor3 = TDatabase::ParamDB->P1 * std::pow(joint_diam, 6);
        
        if (is_boundary_edge)
        {
          ErrThrow("Only the interior is implemented");
        }
        else
        {
          if(matrix_flag == 1)
          {
            integrand += factor1 * 
                       (beta_1 * ny * a_jump[1] + beta_2 * ny * a_jump[2])
                     * (beta_1 * ny * t_jump[1] + beta_2 * ny * t_jump[2]);
            integrand -= factor2 * 
                        (beta_1y * a_jump[1] + beta_1 * a_jump[4] 
                       + beta_2y * a_jump[2] + beta_2 * a_jump[5])
                       *(beta_1y * t_jump[1] + beta_1 * t_jump[4] 
                       + beta_2y * t_jump[2] + beta_2 * t_jump[5]);
            double grad_curl_ansatz 
                       = ( beta_1xy * a_jump[1] + beta_1y * a_jump[3]
                         + beta_1x  * a_jump[4] + beta_1  * 0.0
                         + beta_2xy * a_jump[2] + beta_2y * a_jump[4]
                         + beta_2x  * a_jump[5] + beta_2  * 0.0) +
                         ( beta_1yy * a_jump[1] + beta_1y * a_jump[4]
                         + beta_1y  * a_jump[4] + beta_1  * 0.0
                         + beta_2yy * a_jump[1] + beta_2y * a_jump[5]
                         + beta_2y  * a_jump[5] + beta_2  * 0.0);
            double grad_curl_test = 
                         ( beta_1xy * t_jump[1] + beta_1y * t_jump[3]
                         + beta_1x  * t_jump[4] + beta_1  * 0.0
                         + beta_2xy * t_jump[2] + beta_2y * t_jump[4]
                         + beta_2x  * t_jump[5] + beta_2  * 0.0) +
                         ( beta_1yy * t_jump[1] + beta_1y * t_jump[4]
                         + beta_1y  * t_jump[4] + beta_1  * 0.0
                         + beta_2yy * t_jump[1] + beta_2y * t_jump[5]
                         + beta_2y  * t_jump[5] + beta_2  * 0.0);
                         
            integrand -= factor3 * grad_curl_ansatz * grad_curl_test;
          }
          if(matrix_flag == 3)
          {
            integrand -= factor1 * 
                       (beta_1 * nx * a_jump[1] + beta_2 * nx * a_jump[2])
                     * (beta_1 * nx * t_jump[1] + beta_2 * nx * t_jump[2]);
            integrand += factor2 * 
                        (beta_1x * a_jump[1] + beta_1 * a_jump[3] 
                       + beta_2x * a_jump[2] + beta_2 * a_jump[4])
                     * ( beta_1x * t_jump[1] + beta_1 * t_jump[3] 
                       + beta_2x * t_jump[2] + beta_2 * t_jump[4]);
                     
            double grad_curl_ansatz
                   = ( beta_1xx * a_jump[1] + beta_1x * a_jump[5]
                     + beta_1x  * a_jump[3] + beta_1  * 0.0
                     + beta_2xx * a_jump[2] + beta_2x * a_jump[4]
                     + beta_2x  * a_jump[4] + beta_2  * 0.0 ) + 
                     ( beta_1xy * a_jump[2] + beta_1x * a_jump[4]
                     + beta_1y  * a_jump[3] + beta_1  * 0.0 
                     + beta_2xy * a_jump[2] + beta_2x * a_jump[5]
                     + beta_2y  * a_jump[4] + beta_2 * 0.0 );
            double grad_curl_test
                   = ( beta_1xx * t_jump[1] + beta_1x * t_jump[5]
                     + beta_1x  * t_jump[3] + beta_1  * 0.0
                     + beta_2xx * t_jump[2] + beta_2x * t_jump[4]
                     + beta_2x  * t_jump[4] + beta_2  * 0.0 ) + 
                     ( beta_1xy * t_jump[2] + beta_1x * t_jump[4]
                     + beta_1y  * t_jump[3] + beta_1  * 0.0 
                     + beta_2xy * t_jump[2] + beta_2x * t_jump[5]
                     + beta_2y  * t_jump[4] + beta_2 * 0.0 );
            integrand += factor3 * grad_curl_ansatz * grad_curl_test;
          }
          
        }

        auto weight = orig_weights[quad_pt_i];
        // Update local matrix
        loc_mat[t_dof][a_dof] += weight * integrand;
      } // endfor ansatz dofs
    } // endfor test dofs
  } //endfor quad points
}

void compute_CIP_Stab_Oseen(const ParameterDatabase& param_db, 
    const int& base_vec_dim, 
    const std::vector<double>* val_bas_fct, 
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
    const std::vector<int>& index_map, const double* joint_normal, const double& joint_diam, const TQuadFormula* quad_form_joint, const std::vector<double>& orig_weights, double ** coefficients, 
    double ** parameters, 
    const bool& is_boundary_edge, const int* global_dof, const int* global_dof_neigh, const TBaseCell* cell, int cell_no, int matrix_flag, double ** loc_matA11, double **loc_matA12, 
    double **loc_matA21, double **loc_matA22)
{
  // get some information needed for integration
  auto n_quad_pts = quad_form_joint->GetN_QuadPoints();
  auto n_dofs = val_bas_fct[0].size()/base_vec_dim;
  auto n_dofs_neigh = (is_boundary_edge) ? 0 : val_bas_fct_neigh[0].size()/base_vec_dim;
  auto n_dofs_total = n_dofs + n_dofs_neigh;
  auto nx = joint_normal[0];
  auto ny = joint_normal[1];
  //auto nz = (d == 3) ? joint_normal[2] : 0.0;

  // get parameters of DG method
  auto face_sigma = param_db["face_sigma_DG"].get<double>(); // edge penalty
  auto kappa = param_db["symmetry_DG"].get<int>();  // symmetry
  
  for (int quad_pt_i = 0; quad_pt_i < n_quad_pts; ++quad_pt_i)
  {
    //Output::print(" quad point ", quad_pt_i, "   ", quad_form_joint->get_point(quad_pt_i).x);
    // parameter of PDE
    auto beta_1  = coefficients[quad_pt_i][5];
    auto beta_2  = coefficients[quad_pt_i][6];
    auto beta_1x = coefficients[quad_pt_i][8];
    auto beta_2x = coefficients[quad_pt_i][9];
    auto beta_1y = coefficients[quad_pt_i][10];
    auto beta_2y = coefficients[quad_pt_i][11];
    
    auto beta_1xy = coefficients[quad_pt_i][12];
    // auto beta_1yx = coefficients[quad_pt_i][13];
    auto beta_1xx = coefficients[quad_pt_i][14];
    auto beta_1yy = coefficients[quad_pt_i][15];
    
    auto beta_2xy = coefficients[quad_pt_i][16];
    // auto beta_2yx = coefficients[quad_pt_i][17];
    auto beta_2xx = coefficients[quad_pt_i][18];
    auto beta_2yy = coefficients[quad_pt_i][19];
    
    if(parameters && TDatabase::ParamDB->FLOW_PROBLEM_TYPE == 6)
    {
      beta_1   = parameters[quad_pt_i][0];
      beta_2   = parameters[quad_pt_i][1];
      beta_1x  = parameters[quad_pt_i][2];
      beta_2x  = parameters[quad_pt_i][3];
      beta_1y  = parameters[quad_pt_i][4];
      beta_2y  = parameters[quad_pt_i][5];
      beta_1xy = parameters[quad_pt_i][6];
      beta_1xx = parameters[quad_pt_i][7];
      beta_1yy = parameters[quad_pt_i][8];
      beta_2xy = parameters[quad_pt_i][9];
      beta_2xx = parameters[quad_pt_i][10];
      beta_2yy = parameters[quad_pt_i][11];
      // Output::print(parameters[quad_pt_i][0], " ", parameters[quad_pt_i][1], 
      //              " ", parameters[quad_pt_i][5], " ", parameters[quad_pt_i][8]);
    }


    for (unsigned int t_dof = 0; t_dof < n_dofs_total; ++t_dof)
    { // for all test dofs
      std::vector<double> t_jump(6 * base_vec_dim, 0.);
      std::vector<double> t_average(3 * base_vec_dim);
      if (t_dof < n_dofs)
      {
        // find if the current t_dof is in the global_dof_neigh and it's index
        // i.e. figure out if t_dof is defined on the shared edge or not
          int global_pos = (n_dofs_neigh > 0) ? global_dof[t_dof] : 0;
          int global_pos_neigh_index = -1;

          for (unsigned int j = 0; j < n_dofs_neigh; ++j)
          {
            if (global_pos == global_dof_neigh[j])
            {
              global_pos_neigh_index = j;
              break;
            }
          }

          if (global_pos_neigh_index == -1)
          { // test dof on the boundary edge of the cell
            t_jump[0] = val_bas_fct [quad_pt_i][t_dof];
            t_jump[1] = val_bas_fct_dx [quad_pt_i][t_dof];
            t_jump[2] = val_bas_fct_dy [quad_pt_i][t_dof];
            
            t_jump[3] = val_bas_fct_dxx[quad_pt_i][t_dof];
            t_jump[4] = val_bas_fct_dxy[quad_pt_i][t_dof];
            t_jump[5] = val_bas_fct_dyy[quad_pt_i][t_dof];
          }
          else
          { // test dof on the shared edge
            t_jump[0] =  val_bas_fct[quad_pt_i][t_dof] 
                       - val_bas_fct_neigh[index_map[quad_pt_i]][global_pos_neigh_index];
            t_jump[1] =  val_bas_fct_dx[quad_pt_i][t_dof] 
                       - val_bas_fct_neigh_dx[index_map[quad_pt_i]][global_pos_neigh_index];
            t_jump[2] =  val_bas_fct_dy[quad_pt_i][t_dof] 
                       - val_bas_fct_neigh_dy[index_map[quad_pt_i]][global_pos_neigh_index];
            t_jump[3] =  val_bas_fct_dxx[quad_pt_i][t_dof] 
                       - val_bas_fct_neigh_dxx[index_map[quad_pt_i]][global_pos_neigh_index]; 
            t_jump[4] =  val_bas_fct_dxy[quad_pt_i][t_dof] 
                       - val_bas_fct_neigh_dxy[index_map[quad_pt_i]][global_pos_neigh_index]; 
            t_jump[5] =  val_bas_fct_dyy[quad_pt_i][t_dof] 
                       - val_bas_fct_neigh_dyy[index_map[quad_pt_i]][global_pos_neigh_index];             
          }
      }
      else
      {
        // find if the current t_dof is in the global_dof and it's index
        // i.e. figure out if t_dof is defined on the shared edge or not
          int global_pos_neigh = global_dof_neigh[t_dof - n_dofs];
          int global_pos_index = -1;

          for (unsigned int j = 0; j < n_dofs; ++j)
          {
            if (global_pos_neigh == global_dof[j])
            {
              global_pos_index = j;
              break;
            }
          }
          if (global_pos_index == -1)
          {
            t_jump[0] = -val_bas_fct_neigh[index_map[quad_pt_i]][t_dof - n_dofs];
            t_jump[1] = -val_bas_fct_neigh_dx[index_map[quad_pt_i]][t_dof - n_dofs];
            t_jump[2] = -val_bas_fct_neigh_dy[index_map[quad_pt_i]][t_dof - n_dofs];
            t_jump[3] = -val_bas_fct_neigh_dxx[index_map[quad_pt_i]][t_dof - n_dofs];
            t_jump[4] = -val_bas_fct_neigh_dxy[index_map[quad_pt_i]][t_dof - n_dofs];
            t_jump[5] = -val_bas_fct_neigh_dyy[index_map[quad_pt_i]][t_dof - n_dofs];
        }
        else
        {
        // test dof on the shared edge with the neighbour is 0
            t_jump[0] = 0;
            t_jump[1] = 0;
            t_jump[2] = 0;
            t_jump[3] = 0;
            t_jump[4] = 0;
            t_jump[5] = 0;
        }
      }//else
      if (is_boundary_edge)
      { // Correct values for boundary edges: jump and average on the boundary
          // are defined as values seen from cell
        ErrThrow("We are not considering the boundary edges");
      }
      for (unsigned int a_dof = 0; a_dof < n_dofs_total; ++a_dof)
      {
        std::vector<double> a_jump(6 * base_vec_dim, 0.);
        if (a_dof < n_dofs)
        {
           int global_pos = (n_dofs_neigh > 0) ? global_dof[a_dof] : 0;
           int global_pos_neigh_index = -1;
           for (unsigned int j = 0; j < n_dofs_neigh; ++j)
           {
             if (global_pos == global_dof_neigh[j])
             {
               global_pos_neigh_index = j;
               break;
             }
           }
           if (global_pos_neigh_index == -1)
           { // test dof on the boundary edge of the cell
             a_jump[0] = val_bas_fct [quad_pt_i][a_dof];
             a_jump[1] = val_bas_fct_dx [quad_pt_i][a_dof];
             a_jump[2] = val_bas_fct_dy [quad_pt_i][a_dof];
             a_jump[3] = val_bas_fct_dxx [quad_pt_i][a_dof];
             a_jump[4] = val_bas_fct_dxy [quad_pt_i][a_dof];
             a_jump[5] = val_bas_fct_dyy [quad_pt_i][a_dof];
           }
           else
           { // ansatz dof on the shared edge
             a_jump[0] =  val_bas_fct[quad_pt_i][a_dof] 
                        - val_bas_fct_neigh[index_map[quad_pt_i]][global_pos_neigh_index];
             a_jump[1] =  val_bas_fct_dx[quad_pt_i][a_dof] 
                        - val_bas_fct_neigh_dx[index_map[quad_pt_i]][global_pos_neigh_index];
             a_jump[2] =  val_bas_fct_dy[quad_pt_i][a_dof] 
                        - val_bas_fct_neigh_dy[index_map[quad_pt_i]][global_pos_neigh_index];
             a_jump[3] =  val_bas_fct_dxx[quad_pt_i][a_dof] 
                        - val_bas_fct_neigh_dxx[index_map[quad_pt_i]][global_pos_neigh_index];
             a_jump[4] =  val_bas_fct_dxy[quad_pt_i][a_dof] 
                        - val_bas_fct_neigh_dxy[index_map[quad_pt_i]][global_pos_neigh_index];
             a_jump[5] =  val_bas_fct_dyy[quad_pt_i][a_dof] 
                        - val_bas_fct_neigh_dyy[index_map[quad_pt_i]][global_pos_neigh_index];
           }
        }
        else
        {
          int global_pos_neigh = global_dof_neigh[a_dof - n_dofs];
          int global_pos_index = -1;

          for (unsigned int j = 0; j < n_dofs; ++j)
          {
            if (global_pos_neigh == global_dof[j])
            {
              global_pos_index = j;
              break;
            }
          }
          if (global_pos_index == -1)
          {
            a_jump[0] = -val_bas_fct_neigh[index_map[quad_pt_i]][a_dof - n_dofs];
            a_jump[1] = -val_bas_fct_neigh_dx[index_map[quad_pt_i]][a_dof - n_dofs];
            a_jump[2] = -val_bas_fct_neigh_dy[index_map[quad_pt_i]][a_dof - n_dofs];
            a_jump[3] = -val_bas_fct_neigh_dxx[index_map[quad_pt_i]][a_dof - n_dofs];
            a_jump[4] = -val_bas_fct_neigh_dxy[index_map[quad_pt_i]][a_dof - n_dofs];
            a_jump[5] = -val_bas_fct_neigh_dyy[index_map[quad_pt_i]][a_dof - n_dofs];
          }
          else
          {       
            a_jump[0] = 0;
            a_jump[1] = 0;
            a_jump[2] = 0;
            a_jump[3] = 0;
            a_jump[4] = 0;
            a_jump[5] = 0;
          }
        }
        if (is_boundary_edge)
        {
          ErrThrow("Only the interior is implemented");
        }
        double integranda11 = 0.;
        double integranda12 = 0.;
        double integranda21 = 0.;
        double integranda22 = 0.;
        //double tx = -ny;
        //double ty = nx;
        double factor1 = TDatabase::ParamDB->P1 * std::pow(joint_diam, 2);
        double factor2 = TDatabase::ParamDB->P2 * std::pow(joint_diam, 4);
        double factor3 = TDatabase::ParamDB->P3 * std::pow(joint_diam, 6);
        switch(TDatabase::ParamDB->CIP_TYPE)
        {
          case 0:
          {
            integranda11 += factor1 * std::fabs(beta_1*nx + beta_2 * ny)*
                      (nx * a_jump[1] + ny * a_jump[2])
                    * (nx * t_jump[1] + ny * t_jump[2]);
                    
            integranda22 += integranda11;
            break;
          }
          case 1:
          {
            double anzatz = (beta_1 * a_jump[1] + beta_2 * a_jump[2]);
            double test   = (beta_1 * t_jump[1] + beta_2 * t_jump[2]);
            integranda11 += factor1 * anzatz *ny * test * ny;
            integranda12 -= factor1 * test * ny * anzatz * nx;
            integranda21 -= factor1 * anzatz * nx * test * ny;
            integranda22 += factor1 * anzatz * nx * test * nx;
            
            double ansatz_x = beta_1x * a_jump[1] + beta_1 * a_jump[3]
                            + beta_2x * a_jump[2] + beta_2 * a_jump[4];
            double ansatz_y = beta_1y * a_jump[1] + beta_1 * a_jump[4]
                            + beta_2y * a_jump[2] + beta_2 * a_jump[5];
            double test_x   = beta_1x * t_jump[1] + beta_1 * t_jump[3]
                            + beta_2x * t_jump[2] + beta_2 * t_jump[4];
            double test_y   = beta_1y * t_jump[1] + beta_1 * t_jump[4]
                            + beta_2y * t_jump[2] + beta_2 * t_jump[5];
            integranda11 += factor2 * ansatz_y * test_y;
            integranda12 -= factor2 * ansatz_y * test_x;
            integranda21 -= factor2 * ansatz_x * test_y;
            integranda22 += factor2 * ansatz_x * test_x;
            // S3 A11 block
            double ansatz_p1 = (beta_1x + beta_2y) * a_jump[4] + beta_1xy * a_jump[1]
                              + beta_2xy * a_jump[2] + beta_1y * a_jump[3] 
                              + beta_2x * a_jump[5];
            double test_p1 = (beta_1x + beta_2y) * t_jump[4] + beta_1xy * t_jump[1]
                            + beta_2xy * t_jump[2] + beta_1y * t_jump[3] 
                            + beta_2x * t_jump[5];
            double ansatz_p2 = beta_1yy * a_jump[1] + beta_2yy * a_jump[2]
                             + 2.0* beta_1y * a_jump[4] + 2.0* beta_2y * a_jump[5];
            double test_p2   = beta_1yy * t_jump[1] + beta_2yy * t_jump[2]
                             + 2.0* beta_1y * t_jump[4] + 2.0* beta_2y * t_jump[5];
            integranda11 += factor3 * (ansatz_p1 * test_p1 + ansatz_p2 * test_p2);
            // S3 A12 block 
            test_p1 = (beta_1x + beta_2y) * t_jump[4] + beta_1xy * t_jump[1] 
                     + beta_2xy * t_jump[2] + beta_1y * t_jump[3] 
                     + beta_2x * t_jump[5];
            ansatz_p1 = beta_1xx * a_jump[1] + 2.0*beta_1x * a_jump[3] 
                      + beta_2xx * a_jump[2] + 2.0 * beta_2x * a_jump[4];
            ansatz_p2 = (beta_1x + beta_2y) * a_jump[4] + beta_1xy * a_jump[1]
                       + beta_2xy * a_jump[2] + beta_1y * a_jump[3] 
                       + beta_2x * a_jump[5];
            test_p2 = beta_1yy * t_jump[1] + beta_2yy * t_jump[2] 
                    + 2.0*beta_1y * t_jump[4] + 2.0 * beta_2y * t_jump[5];
            integranda12 -= factor3 * (test_p1 * ansatz_p1 + test_p2 * ansatz_p2);
            // S3 A21 block
            ansatz_p1 = (beta_1x + beta_2y) * a_jump[4] + beta_1xy * a_jump[1] 
                       + beta_2xy * a_jump[2] + beta_1y * a_jump[3] 
                       + beta_2x * a_jump[5];
            test_p1 = beta_1xx * t_jump[1] + 2.0*beta_1x * t_jump[3] 
                    + beta_2xx * t_jump[2] + 2.0 * beta_2x * t_jump[4];
            ansatz_p2 = beta_1yy * a_jump[1] + beta_2yy * a_jump[2] 
                      + 2.*beta_1y * a_jump[4] + 2.*beta_2y * a_jump[5];
            test_p2 = (beta_1x + beta_2y) * t_jump[4] + beta_1xy * t_jump[1] 
                     + beta_2xy * t_jump[2] + beta_1y * t_jump[3] 
                     + beta_2x * t_jump[5];
            integranda21 -= factor3 * (ansatz_p1 * test_p1 + ansatz_p2 * test_p2);
            //S3 A22 block
            ansatz_p1 = beta_1xx * a_jump[1] + 2.*beta_1x * a_jump[3]
                      + beta_2xx * a_jump[2] + 2.*beta_2x * a_jump[4];
            test_p1 = beta_1xx * t_jump[1] + 2.*beta_1x * t_jump[3]
                    + beta_2xx * t_jump[2] + 2.*beta_2x * t_jump[4];
            ansatz_p2 = (beta_1x + beta_2y) * a_jump[4] + beta_1xy * a_jump[1]
                      + beta_2xy * a_jump[2] + beta_1y * a_jump[3] 
                      + beta_2x * a_jump[5];
            test_p2 = (beta_1x + beta_2y) * t_jump[4] + beta_1xy * t_jump[1]
                      + beta_2xy * t_jump[2] + beta_1y * t_jump[3] 
                      + beta_2x * t_jump[5];
            integranda22 += factor3 * (ansatz_p1 * test_p1 + ansatz_p2 * test_p2);
            break;
          }
          case 2:
          {
            break;
          }
        }        
        auto weight = orig_weights[quad_pt_i];
        loc_matA11[t_dof][a_dof] += weight * integranda11;
        loc_matA12[t_dof][a_dof] += weight * integranda12;
        loc_matA21[t_dof][a_dof] += weight * integranda21;
        loc_matA22[t_dof][a_dof] += weight * integranda22;
      }//endfor ansatz dofs
    } //endfor t_dof    
  }// endfor quad_pt_i
}
#endif

#ifdef __2D__

void AddJumpStabilizationCIP(int n_fespaces, const TFESpace2D ** fespaces, int n_sqmatrices, TSquareMatrix2D ** sqmatrices, BoundCondFunct2D ** BoundaryConditions, BoundValueFunct2D ** BoundaryValues, const ParameterDatabase& param_db, LocalAssembling2D& la)
{
  auto db = default_dg_database();
  db.merge(param_db, true);
  // Prepare loop over cells -> mark cells to access cell numbers
  auto collection = fespaces[0]->GetCollection();
  collection->mark_all_cells(); // set cell_index
  auto n_cells = collection->GetN_Cells();  // number of cells

  // Set up quadrature rule on facets
  // Find out maximal degree needed for integration, i.e. the largest degree
  // over all cells. Note that on single cells this degree can be too large if
  // the degree differs between cell, and hence a computational overhead exists.
  // Fix it, if you think it is a substantial problem.
  int max_degree = 0;
  for (int cell_i = 0; cell_i < n_cells; ++cell_i)
  {
    for(int fe_i = 0; fe_i < n_fespaces; fe_i++)
    {
      auto& element = fespaces[fe_i]->get_fe(cell_i);
      auto degree = element.GetBaseFunct()->GetPolynomialDegree();
      max_degree = (degree > max_degree) ? degree : max_degree;
    }
  }
  // There will be integrands of the type base_fct * base_fct. Hence, the
  // quadrature rule has to be exact for at least 2 * degree.
  // This in principle allows for different polynomial degrees in two
  // adjacent cells, since max_degree is large enough, but at the price that the
  // quadrature rule may be too complicated in possibly many cells.
  BFRefElements joint_type;
  
  // Currently in 3D mixed grids are not possible, i.e. the mesh
  // consists of either hexahedra or tetrahedra. Therefore, only one
  // quadrature rule is needed depending on the type of cells. A criterion for
  // distinguishing between tetrahedra and hexahedra is e.g. the number of
  // vertices in the first cell ( 4 for tetrahedra and 6 for hexahedra ).
  joint_type = BFRefElements::BFUnitLine;
  
  auto quad_formula_joint = QuadratureFormulaDatabase::qf_from_degree(
      2 * max_degree, joint_type);
  auto n_quad_pts = quad_formula_joint->GetN_QuadPoints();
  
  auto coeff_fct = la.GetCoeffFct();
  std::vector<double*> parameters(n_quad_pts);
  auto n_params = la.GetN_Parameters();
  std::vector<double> aux(n_quad_pts * n_params);
  for(int quad_pt_i = 0; quad_pt_i < n_quad_pts; quad_pt_i++)
  {
    parameters[quad_pt_i] = aux.data() + quad_pt_i * n_params;
  }
  // Perform loop over cells to loop over local edges
  for (int cell_i = 0; cell_i < n_cells; ++cell_i)
  { // loop over cells
    //Output::print("cell ", cell_i);
    auto finite_element = fespaces[0]->get_fe(cell_i);
    auto cell = collection->GetCell(cell_i);
    auto n_joints = cell->GetN_Joints();
    
    la.GetParameters(*quad_formula_joint, cell_i, parameters.data());
    
    cell-> SetNormalOrientation();

    // for all edges
    for (int joint_i = 0; joint_i < n_joints; ++joint_i)
    {
      //Output::print(" joint ", joint_i);
      // Calculate the values of the basis function and their derivatives in
      // quad pts for this edge
      std::vector<std::vector<double>> val_bas_fct;
      std::vector<std::vector<double>> val_bas_fct_dx;
      std::vector<std::vector<double>> val_bas_fct_dy;
      std::vector<std::vector<double>> val_bas_fct_dz;
      std::vector<std::vector<double>> val_bas_fct_uxx;
      std::vector<std::vector<double>> val_bas_fct_uxy;
      std::vector<std::vector<double>> val_bas_fct_uyy;
      // val_bas_fct[j] contains the values of the basis functions on the
      // original cell evaluated at quadrature point j.
      // val_bas_fct_dx (_dy, _dz) are the derivatives wrt. x (y, z)
      // This is filled in the next lines
      auto base_vec_dim = finite_element.GetBaseFunct()->GetBaseVectDim();
      compute_basis_fct_values_at_joint<2>(finite_element, cell, joint_i, 
                base_vec_dim, *quad_formula_joint, val_bas_fct, val_bas_fct_dx,
               val_bas_fct_dy, 
               val_bas_fct_uxx, val_bas_fct_uxy, val_bas_fct_uyy,
               val_bas_fct_dz, *collection);
      

      auto ref_element = finite_element.GetBaseFunct()->GetRefElement();
      auto ref_trans_id = finite_element.GetRefTransID();
      auto transformed_quad_points = transform_quad_points(cell, ref_element,
      ref_trans_id, joint_i, quad_formula_joint);



      // Check for neighbours to distinguish between boundary and inner edges
      auto joint = cell->GetJoint(joint_i);
      auto neigh = joint->GetNeighbour(cell);

      // Not every facet has to be considered. Only Dirichlet boundary facets
      // and interior facet have to be considered. Neumann boundary facets are
      // already treated by previously called assemble routines. Each interior
      // facets only needs to be considered once. This is checked by the
      // following boolean.
      bool consider_joint = false;

      // Storage for values of neighbor basis function evaluated at the quad
      // pts. This is filled only in the case of interior edges.
      std::vector<std::vector<double>> val_bas_fct_neigh;
      std::vector<std::vector<double>> val_bas_fct_neigh_dx;
      std::vector<std::vector<double>> val_bas_fct_neigh_dy;
      std::vector<std::vector<double>> val_bas_fct_neigh_uxx;
      std::vector<std::vector<double>> val_bas_fct_neigh_uxy;
      std::vector<std::vector<double>> val_bas_fct_neigh_uyy;
      std::vector<std::vector<double>> val_bas_fct_neigh_dz;
      std::vector<int> index_map;
      int cell_nr_neigh = -1; // dummy value to remove warning
      if (neigh)
      { // interior edge
        // In MPI case unfortunately the unused cells are not deleted correctly.
        // This results in the fact that cell think they have neighbours which
        // are not in the collection. The only way I found to exclude those
        // cells is to check the clipboard even if I try to avoid this. The
        // clipboard of these "ghost cells" is set to -1 as for cells with
        // hanging nodes. But those cells have to be excluded
        if (neigh->GetClipBoard() >= 0 )
        {
          cell_nr_neigh = neigh->GetCellIndex(); // cell number of neighbor

          // Every interior edge has to be considered only once. We consider it
          // only at their first appearance as local edge, i.e. if cell_i is
          // smaller than cell_nr_neigh. In this case also the normal is defined
          // as the outer unit normal of cell i.
          if (cell_i < cell_nr_neigh)
          { // interior edges was not considered yet
            consider_joint = true;

            // find out the local joint number corresponding to joint_i
            auto facet_nr_neigh = joint->get_joint_nr_in_cell(neigh);

            // val_bas_fct_neigh[j] contains the values of the basis functions
            // on the original cell evaluated at quadrature point j.
            // val_bas_fct_neigh_dx (_y, _z) are the derivatives wrt. x (y, z)
            // This is filled in the next line
            auto fe_neigh = fespaces[0]->get_fe(cell_nr_neigh);
            int base_vec_dim_ne = fe_neigh.GetBaseFunct()->GetBaseVectDim();            
            compute_basis_fct_values_at_joint<2>(fe_neigh, neigh, facet_nr_neigh, 
                base_vec_dim_ne, *quad_formula_joint, val_bas_fct_neigh, val_bas_fct_neigh_dx,
               val_bas_fct_neigh_dy, 
               val_bas_fct_neigh_uxx, val_bas_fct_neigh_uxy, val_bas_fct_neigh_uyy,
               val_bas_fct_neigh_dz, *collection);
            // Compute an index map that maps the indices of the quadrature
            // points seen from the cell to the quadrature points seen form the
            // neighboring cell
            index_map = compute_index_map_for_facet_quad_points_of_neigh<2>(
                cell, neigh, joint_i, facet_nr_neigh, finite_element, fe_neigh,
                quad_formula_joint );
          } // endif cell_i < cell_nr_neigh
        }
      } // endif neigh
      else
      {
        
      }

      if (consider_joint)
      {
        // Compute coefficients of PDE in quadrature points on edge_i
        std::vector<double*> coefficients(n_quad_pts);
        std::vector<double> aux(n_quad_pts * n_local_coeff);
        for(int quad_pt_i = 0; quad_pt_i < n_quad_pts; quad_pt_i++)
        {
          coefficients[quad_pt_i] = aux.data() + quad_pt_i * n_local_coeff;
        }
        compute_coefficients<2>(cell, joint_i, finite_element,
            quad_formula_joint, coeff_fct, coefficients.data());

        // compute geometrical information needed for integration
        auto facet_diam = cell->ComputeDiameterOfJoint(joint_i);
        auto face_normal = cell->ComputeOuterJointNormal(joint_i);

        // allocate local matrix. Upper-left square of size n_dofs^2
        // correspond to dofs from cell.
        // In case of interior edges, lower-right square of size
        // n_dofs_neigh^2 corresponds to dofs of neighbor cell, upper-right
        // and lower-left part of size n_dofs * n_dofs_neigh correspond to
        // coupling between test and ansatz dofs of cell and neighbour and
        // other way around.
        auto n_dofs = finite_element.GetN_DOF();
        int n_dofs_neigh = 0;
        if (neigh)
        {
          auto element_neigh = fespaces[0]->get_fe(cell_nr_neigh);
          n_dofs_neigh = element_neigh.GetN_DOF();
        }
        auto n_dof_total = n_dofs + n_dofs_neigh;
        std::vector<double*> loc_matA11(n_dof_total);
        std::vector<double*> loc_matA12(n_dof_total);
        std::vector<double*> loc_matA21(n_dof_total);
        std::vector<double*> loc_matA22(n_dof_total);
        std::vector<double> pointer_vec11(n_dof_total * n_dof_total, 0);
        std::vector<double> pointer_vec12(n_dof_total * n_dof_total, 0);
        std::vector<double> pointer_vec21(n_dof_total * n_dof_total, 0);
        std::vector<double> pointer_vec22(n_dof_total * n_dof_total, 0);
        for (int i = 0; i < n_dof_total; ++i)
        {
          loc_matA11[i] = &pointer_vec11[i * n_dof_total]; // local matrix
          loc_matA12[i] = &pointer_vec12[i * n_dof_total]; // local matrix
          loc_matA21[i] = &pointer_vec21[i * n_dof_total]; // local matrix
          loc_matA22[i] = &pointer_vec22[i * n_dof_total]; // local matrix
        }

        auto orig_weights = compute_orig_weights(cell, joint_i,
            quad_formula_joint);
        auto base_vec_dim = finite_element.GetBaseFunct()->GetBaseVectDim();

        auto global_dof = fespaces[0]->GetGlobalDOF(cell_i);
        if (neigh)
        {
            
#ifdef __2D__              
              auto global_dof_neigh = fespaces[0]->GetGlobalDOF(cell_nr_neigh);
                compute_CIP_Stab_Oseen(db, base_vec_dim, 
                        val_bas_fct.data(),
                        val_bas_fct_dx.data(), 
                        val_bas_fct_dy.data(),
                        val_bas_fct_uxx.data(), 
                        val_bas_fct_uxy.data(),
                        val_bas_fct_uyy.data(),
                        val_bas_fct_dz.data(), 
                        val_bas_fct_neigh.data(),
                        val_bas_fct_neigh_dx.data(), 
                        val_bas_fct_neigh_dy.data(),
                        val_bas_fct_neigh_uxx.data(),
                        val_bas_fct_neigh_uxy.data(),
                        val_bas_fct_neigh_uyy.data(),
                        val_bas_fct_neigh_dz.data(), 
                        index_map, face_normal.data(),
                        facet_diam, quad_formula_joint, orig_weights,
                        coefficients.data(), 
                        parameters.data(),
                        false, 
                        global_dof, 
                        global_dof_neigh, cell, cell_i, 0, 
                        loc_matA11.data(), loc_matA12.data(), 
                        loc_matA21.data(), loc_matA22.data());
#endif
        }
        else
        {
          
        }

        // compute local right-hand side
        std::vector<double> loc_rhs(n_dofs, 0);
        if (!neigh)
        {
          //NOTE: nothing to do on the boundary because we are only interested in the interior
        }

        // add local to global
        if (neigh && neigh->GetClipBoard() >= 0)
        {          
          auto global_dof_neigh = fespaces[0]->GetGlobalDOF(cell_nr_neigh);          
          auto n_dof_total = (global_dof_neigh == nullptr) ? n_dofs : n_dofs + n_dofs_neigh;
          
          for(int t_dof = 0; t_dof < n_dof_total; t_dof++)
          {
            auto t_global_number = (t_dof < n_dofs) ? global_dof[t_dof] :
            global_dof_neigh[t_dof - n_dofs];  // global dof number
            
            if(t_global_number >= sqmatrices[0]->get_n_active_rows())
              continue;   
            
            auto loc_rowa11 = loc_matA11.data()[t_dof];
            auto loc_rowa22 = loc_matA22.data()[t_dof];
            auto loc_rowa12 = loc_matA12.data()[t_dof];
            auto loc_rowa21 = loc_matA21.data()[t_dof];
            
            for(int a_dof = 0; a_dof < n_dof_total; a_dof++)
            {
              auto a_global_number = (a_dof < n_dofs) ? global_dof[a_dof] : 
                                          global_dof_neigh[a_dof - n_dofs];
              sqmatrices[0]->add(t_global_number, a_global_number, loc_rowa11[a_dof]);
              sqmatrices[1]->add(t_global_number, a_global_number, loc_rowa12[a_dof]);
              sqmatrices[2]->add(t_global_number, a_global_number, loc_rowa21[a_dof]);
              sqmatrices[3]->add(t_global_number, a_global_number, loc_rowa22[a_dof]);
            }//endfor ansatz dof
          }//endfor test dofs
        }        
      }
    } // endfor loop over edges
  } // endfor loop over cells
}
#endif


