// ======================================================================
// @(#)BoundaryAssembling2D.C        05/18/16
//
// Functions for (external and internal) boundary integral
//
// authors: Alfonso Caiazzo and Laura Blank
// ======================================================================
/*
   List of possible improvements
 * link each cell with its FE (in order not to use a global function to get FE properties)
 * create an FE class, in order not to recover all properties via FEDatabase::FunctionName(FeID)
 * rewrite the GetFormulaData using (e.g.) vector<> class
 * avoid to define double* with MaxN... (use vector instead?)
 * alternative: pass a list of joints (each with its own FE space on it)?
 * the edge list could be generated outside this class
 */

#include "FEDatabase.h"
#include <BoundaryAssembling2D.h>
#include <BoundEdge.h>
#include <Collection.h>
#include "BaseCell.h"
#include "Database.h"
#include "InterfaceJoint.h"
#include "TriaIsoparametric.h"
#include "QuadIsoparametric.h"
#include <QuadratureFormulaDatabase.h>


//========================================================================
void BoundaryAssembling2D::rhs_g_v_n(BlockVector &rhs,
    const TFESpace2D *U_Space,
    BoundValueFunct2D *given_boundary_data,
    int boundary_component_id,
    double mult)
{
  // Create a list of those boundary edges that are on the boundary component with given ID
  std::vector<TBoundEdge*> boundaryEdgeList;
  const TCollection *coll = U_Space->GetCollection();
  coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);

  rhs_g_v_n(rhs, U_Space, given_boundary_data, boundaryEdgeList, mult);
}

void BoundaryAssembling2D::rhs_g_v_n(BlockVector &rhs,
    const TFESpace2D *U_Space,
    BoundValueFunct2D *given_boundary_data,
    std::vector<TBoundEdge*> &boundaryEdgeList,
    double mult)
{
 // Go through all boundary edges on the current boundary component
  for (size_t m = 0; m < boundaryEdgeList.size(); m++)
  {
    // current edge
    TBoundEdge *boundedge = boundaryEdgeList[m];
    TBaseCell *cell = boundedge->GetNeighbour(0);
    // mapping from local (cell) DOF to global DOF
    const int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex());

    // --------------------------------------------------------------
    ///@todo put the following part into FESpace2D::getEdgeQuadratureData
    // compute values of all basis functions
    // and their first partial derivatives at all quadrature points
    // get basis dimension and FE space data of the current cell
    auto element = U_Space->get_fe(cell->GetCellIndex());
    int BaseVectDim = 1; // we assume only scalar FE
    int joint_id = boundedge->get_index_in_neighbour(cell);
    // get a quadrature formula good enough for the argument of the integral
    int fe_degree = element.GetBaseFunct()->GetPolynomialDegree();
    auto LineQuadFormula = QuadratureFormulaDatabase::qf_from_degree(
        2 * fe_degree, BFRefElements::BFUnitLine);
    std::vector<double> quadWeights, quadPoints;
    get_quadrature_formula_data(quadPoints, quadWeights, *LineQuadFormula);

    std::vector< std::vector<double> > uorig, u_dx_orig, u_dy_orig;
    get_original_values(element, joint_id, cell, quadPoints, BaseVectDim,
        uorig, u_dx_orig, u_dy_orig,  *LineQuadFormula);

    // get normal and edge length
    double x_0, x_1, y_0, y_1;
    boundedge->get_vertices(x_0, y_0, x_1, y_1);
    double joint_length = boundedge->get_length();
    // normal vector to this boundary (normalized)
    double n1, n2;
    boundedge->get_normal(n1, n2);
    double reference_joint_length = 2;
    double transformationDeterminant = joint_length/reference_joint_length;

    // ***** COMPUTE INTEGRAL *****
    for (unsigned int k = 0; k < quadPoints.size(); k++)
    {
      // get the value of (\nabla u - pI) on the boundary component (here denoted by g) 
      double value;

      if (given_boundary_data != nullptr)
      {
        double x = x_0 + (quadPoints[k]+1.)/2. * (x_1-x_0);
        double y = y_0 + (quadPoints[k]+1.)/2. * (y_1-y_0);
        double T;
        boundedge->GetBoundComp()->GetTofXY(x, y, T);
        given_boundary_data(boundedge->GetBoundComp()->GetID(), T, value);
      }
      else
      {
        value = 1;
      }

      double commonFactor = mult * quadWeights[k] * transformationDeterminant;

      for (unsigned int l = 0; l < uorig[k].size(); l++)
      {
        int global_dof_from_local = DOF[l];

        // if the DOF is Dirichlet, continue
        if (global_dof_from_local < U_Space->get_n_active())
        {
          // updating rhs: int_gamma rhsval v \cdot n
          double v1 = uorig[k][l]; // value of test function (vtest = v1 = v2)
          double v2 = v1;

          // add for both components
          rhs.block(0)[global_dof_from_local] += commonFactor * value * v1 * n1;

          if (rhs.n_blocks() == 2)
          {
            rhs.block(1)[global_dof_from_local] += commonFactor * value * v2 * n2;
          }
        }
      }
    }
  }
}

//=================================================================================
void BoundaryAssembling2D::rhs_uD_v(BlockVector &rhs,
    const TFESpace2D *U_Space,
    BoundValueFunct2D *given_boundary_data1,
    BoundValueFunct2D *given_boundary_data2,
    int boundary_component_id,
    double mult,
    bool rescale_by_h)
{
  std::vector<TBoundEdge*> boundaryEdgeList;
  const TCollection *coll = U_Space->GetCollection();
  coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);
  rhs_uD_v(rhs, U_Space, given_boundary_data1, given_boundary_data2,
      boundaryEdgeList, mult, rescale_by_h);
}

void BoundaryAssembling2D::rhs_uD_v(BlockVector &rhs,
    const TFESpace2D *U_Space,
    BoundValueFunct2D *given_boundary_data1,
    BoundValueFunct2D *given_boundary_data2,
    std::vector<TBoundEdge*> &boundaryEdgeList,
    double mult,
    bool rescale_by_h)
{
  int ActiveBound = U_Space->get_n_active();

  for (size_t m = 0; m < boundaryEdgeList.size(); m++)
  {
    TBoundEdge *boundedge = boundaryEdgeList[m];
    TBaseCell *cell = boundedge->GetNeighbour(0);

    // get basis dimension and FE space data of cell i
    auto element = U_Space->get_fe(cell->GetCellIndex());

    int BDComponent = boundedge->GetBoundComp()->GetID();

    int BaseVectDim = 1; // we assume only scalar FE // Only for BDM and RT elements \neq 1
    int joint_id = boundedge->get_index_in_neighbour(cell);

    // get a quadrature formula good enough for the velocity FE space (here exact to 2*fe_degree)
    int fe_degree = element.GetBaseFunct()->GetPolynomialDegree();
    auto LineQuadFormula = QuadratureFormulaDatabase::qf_from_degree(
        2*fe_degree, BFRefElements::BFUnitLine);
    std::vector<double> quadWeights, quadPoints;
    get_quadrature_formula_data(quadPoints, quadWeights, *LineQuadFormula);

    // compute values of all basis functions at all quadrature points
    std::vector< std::vector<double> > uorig, u_dx_orig, u_dy_orig;
    get_original_values(element, joint_id, cell, quadPoints, BaseVectDim, uorig,
                        u_dx_orig, u_dy_orig, *LineQuadFormula);

    double x_0, x_1, y_0, y_1;
    boundedge->get_vertices(x_0, y_0, x_1, y_1);
    // compute length of the edge
    double joint_length = boundedge->get_length();

    // quadrature
    for (unsigned int k = 0; k < quadPoints.size(); k++)
    {
      ///@attention in 1D the reference joint is [-1,1] => length = 2
      double reference_joint_length = 2;
      double x = x_0+(quadPoints[k]+1.)/2.*(x_1-x_0);
      double y = y_0+(quadPoints[k]+1.)/2.*(y_1-y_0);

      double T;
      boundedge->GetBoundComp()->GetTofXY(x, y, T);

      // get the boundary values of rhs
      double value1, value2;
      if (given_boundary_data1 != nullptr)
      {
        given_boundary_data1(BDComponent, T, value1);
      }
      else
      {
        Output::print<1>("WARNING: Due to missing input (nullptr) in "
            "BoundaryAssembling2D::rhs_uD_v(), a default value for "
            "given_boundary_data1 is used, which might not fit your example.");
        value1 = 1.;
      }

      if(given_boundary_data2 != nullptr)
      {
        given_boundary_data2(BDComponent, T, value2);
      }
      else
      {
        Output::print<1>("WARNING: Due to missing input (nullptr) in "
            "BoundaryAssembling2D::rhs_uD_v(), a default value for "
            "given_boundary_data2 is used, which might not fit your example.");
        value2 = 0.0;
      }

      // mapping from local (cell) DOF to global DOF
      const int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex()); //BeginIndex[i];

      for (unsigned int l = 0; l < uorig[k].size(); l++)
      {
        int global_dof_from_local = DOF[l];

        // if the DOF is Dirichlet, continue
        if (global_dof_from_local >= ActiveBound)
          continue;

        // updating rhs: int_gamma rhsval[2] v
        double v1 = uorig[k][l]; // value of test function (vtest = vx = vy)
        double v2 = v1;

        // add for both components
        if (!rescale_by_h)
        {
          rhs.block(0)[global_dof_from_local] += mult * quadWeights[k] * value1 * v1 *
            (joint_length/reference_joint_length);

          rhs.block(1)[global_dof_from_local] += mult * quadWeights[k] * value2 * v2 *
            (joint_length/reference_joint_length);
        }
        else
        {
          rhs.block(0)[global_dof_from_local] += (mult * quadWeights[k] * value1 * v1 *
              (joint_length/reference_joint_length)) / joint_length;

          rhs.block(1)[global_dof_from_local] += (mult * quadWeights[k] * value2 * v2 *
              (joint_length/reference_joint_length)) /joint_length;
        }
      }
    }
  }
}

//========================================================================
void BoundaryAssembling2D::rhs_gradv_n_uD(BlockVector &rhs,
    const TFESpace2D *U_Space,
    BoundValueFunct2D *given_boundary_data1,
    BoundValueFunct2D *given_boundary_data2,
    int boundary_component_id,
    double mult)
{
  std::vector<TBoundEdge*> boundaryEdgeList;
  const TCollection *coll = U_Space->GetCollection();
  coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);
  rhs_gradv_n_uD(rhs, U_Space, given_boundary_data1, given_boundary_data2, boundaryEdgeList, mult);
}

void BoundaryAssembling2D::rhs_gradv_n_uD(BlockVector &rhs,
    const TFESpace2D *U_Space,
    BoundValueFunct2D *given_boundary_data1,
    BoundValueFunct2D *given_boundary_data2,
    std::vector<TBoundEdge*> &boundaryEdgeList,
    double mult)
{
  int ActiveBound = U_Space->get_n_active();

  for (size_t m = 0; m < boundaryEdgeList.size(); m++)
  {
    TBoundEdge *boundedge = boundaryEdgeList[m];
    TBaseCell *cell = boundedge->GetNeighbour(0);

    // get basis dimension and FE space data of cell i
    auto element = U_Space->get_fe(cell->GetCellIndex());

    int BDComponent = boundedge->GetBoundComp()->GetID();

    int BaseVectDim = 1; // we assume only scalar FE // Only for BDM and RT elements \neq 1
    int joint_id = boundedge->get_index_in_neighbour(cell);

    // get a quadrature formula good enough for the velocity FE space (here exact to 2*fe_degree)
    int fe_degree = element.GetBaseFunct()->GetPolynomialDegree();
    auto LineQuadFormula = QuadratureFormulaDatabase::qf_from_degree(
        2*fe_degree, BFRefElements::BFUnitLine);
    std::vector<double> quadWeights, quadPoints;
    get_quadrature_formula_data(quadPoints, quadWeights, *LineQuadFormula);

    // compute values of all basis functions at all quadrature points
    std::vector< std::vector<double> > uorig, u_dx_orig, u_dy_orig;
    get_original_values(element, joint_id, cell, quadPoints, BaseVectDim, uorig,
                        u_dx_orig, u_dy_orig, *LineQuadFormula);

    double x_0, x_1, y_0, y_1;
    boundedge->get_vertices(x_0, y_0, x_1, y_1);
    // compute length of the edge
    double joint_length = boundedge->get_length();
    // normal vector to this boundary (normalized)
    double n1, n2;
    boundedge->get_normal(n1, n2);

    // quadrature
    for (unsigned int k = 0; k < quadPoints.size(); k++)
    {
      ///@attention in 1D the reference joint is [-1,1] => length = 2
      double reference_joint_length = 2;
      double x = x_0 + (quadPoints[k]+1.)/2. * (x_1-x_0);
      double y = y_0 + (quadPoints[k]+1.)/2. * (y_1-y_0);

      double T;
      boundedge->GetBoundComp()->GetTofXY(x, y, T);

      // get the boundary values of rhs
      double value1, value2;
      if (given_boundary_data1 != nullptr)
      {
        given_boundary_data1(BDComponent, T, value1);
      }
      else
      {
        Output::print<1>("WARNING: Due to missing input (nullptr) in "
            "BoundaryAssembling2D::rhs_gradv_n_uD(), a default value for "
            "given_boundary_data1 is used, which might not fit your example.");
        value1 = 1;
      }

      if(given_boundary_data2 != nullptr)
      {
        given_boundary_data2(BDComponent, T, value2);
      }
      else
      {
        Output::print<1>("WARNING: Due to missing input (nullptr) in "
            "BoundaryAssembling2D::rhs_gradv_n_uD(), a default value for "
            "given_boundary_data2 is used, which might not fit your example.");
        value2 = 0.0;
      }

      // mapping from local (cell) DOF to global DOF
      const int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex()); //BeginIndex[i];

      for (unsigned int l = 0; l < uorig[k].size(); l++)
      {
        int global_dof_from_local = DOF[l];

        // if the DOF is Dirichlet, continue
        if (global_dof_from_local >= ActiveBound)
          continue;

        // updating rhs: int_gamma rhsval[2] v
        double v1_dx = u_dx_orig[k][l]; // value of test function (vtest = (vx,vy) )
        double v1_dy = u_dy_orig[k][l];
        double v2_dy = v1_dy;
        double v2_dx = v1_dx;

        rhs.block(0)[global_dof_from_local] += mult * quadWeights[k] * value1 *
            (v1_dx * n1 + v1_dy * n2) *  (joint_length/reference_joint_length);

        rhs.block(1)[global_dof_from_local] += mult * quadWeights[k] * value2 *
            (v2_dx * n1 + v2_dy * n2) * (joint_length/reference_joint_length);
      }
    }
  }
}

//======================================================================
void BoundaryAssembling2D::matrix_u_n_v_n(BlockFEMatrix &M,
    const TFESpace2D *U_Space,
    int boundary_component_id,
    double mult,
    bool rescale_by_h)
{
  std::vector<TBoundEdge*> boundaryEdgeList;
  const TCollection *coll = U_Space->GetCollection();
  coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);

  matrix_u_n_v_n(M, U_Space, boundaryEdgeList, mult, rescale_by_h);
}

/**
 * @attention this functions assumes implicitely Matrix Type 14.
 * This means that the blocks are ordered like: A11,A12,B1t,A21,A22,B2t,B1,B2,C)
 * hence we need: blocks[0], blocks[1], blocks[3], blocks[4]
 * for A11, A12, A21, A22
 **/
void BoundaryAssembling2D::matrix_u_n_v_n(BlockFEMatrix &M,
    const TFESpace2D *U_Space,
    std::vector<TBoundEdge*> &boundaryEdgeList,
    double mult,
    bool rescale_by_h)
{
  int ActiveBound = U_Space->get_n_active();

  std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();

  /**
   * @todo: check if the matrix structure is correct:
   * we need 4 square matrices with the same FE spaces
   */

  for (size_t m = 0; m < boundaryEdgeList.size(); m++)
  {
    TBoundEdge *boundedge = boundaryEdgeList[m];
    TBaseCell *cell = boundedge->GetNeighbour(0);
    // get basis dimension and FE space data of cell i
    auto element = U_Space->get_fe(cell->GetCellIndex());

    int BaseVectDim = 1; // we assume only scalar FE
    int joint_id = boundedge->get_index_in_neighbour(cell);

    // get a quadrature formula good enough for the velocity FE space
    int fe_degree = element.GetBaseFunct()->GetPolynomialDegree();
    auto LineQuadFormula = QuadratureFormulaDatabase::qf_from_degree(
        2*fe_degree, BFRefElements::BFUnitLine);
    std::vector<double> quadWeights, quadPoints;
    get_quadrature_formula_data(quadPoints, quadWeights, *LineQuadFormula);

    //FEDatabase::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(
    //  this->LineQuadFormula);

    // compute values of all basis functions at all quadrature points
    std::vector< std::vector<double> > uorig, uxorig, uyorig;
    get_original_values(element, joint_id, cell, quadPoints, BaseVectDim, uorig,
                        uxorig, uyorig, *LineQuadFormula);

    double x0, x1, y0, y1;
    boundedge->get_vertices(x0, y0, x1, y1);
    // compute length of the edge
    double joint_length = boundedge->get_length();
    // normal vector to this boundary (normalized)
    double n1, n2;
    boundedge->get_normal(n1, n2);

    // quadrature
    for (unsigned int k = 0; k < quadPoints.size(); k++)
    {
      ///@attention in 1D the reference joint is [-1,1] => length = 2
      double reference_joint_length = 2;

      // mapping from local(cell) DOF to global DOF
      const int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex()); //BeginIndex[i];

      double scale_factor = mult * quadWeights[k] * (joint_length/reference_joint_length);
      if (rescale_by_h)
      {
        scale_factor = scale_factor/joint_length;
      }

      // loop on test functions
      for (unsigned int l1 = 0; l1 < uorig[k].size(); l1++)
      {
        int test_DOF = DOF[l1];

        // if the DOF is Dirichlet, continue
        if(test_DOF >= ActiveBound)
          continue;

        double v1 = uorig[k][l1];
        double v2 = v1; // x and y component have the same FE space

        // loop on ansatz functions
        for (unsigned int l2 = 0; l2 < uorig[k].size(); l2++)
        {
          int ansatz_DOF = DOF[l2];
          double u1 = uorig[k][l2];
          double u2 = u1; // x and y component have the same FE space

          // (see the note about blocks at the beginning of the function)
          blocks[0]->add( test_DOF, ansatz_DOF, scale_factor * (v1 * n1) * (u1 * n1) ); // A11
          blocks[1]->add( test_DOF, ansatz_DOF, scale_factor * (v1 * n1) * (u2 * n2) ); // A12
          blocks[3]->add( test_DOF, ansatz_DOF, scale_factor * (v2 * n2) * (u1 * n1) ); // A21
          blocks[4]->add( test_DOF, ansatz_DOF, scale_factor * (v2 * n2) * (u2 * n2) ); // A22
        }
      }
    }
  }
}

//======================================================================
void BoundaryAssembling2D::rhs_uD_n_v_n(BlockVector &rhs,
        const TFESpace2D *U_Space,
        BoundValueFunct2D *given_boundary_data1,
        BoundValueFunct2D *given_boundary_data2,
    int boundary_component_id,
    double mult,
    bool rescale_by_h)
{
  std::vector<TBoundEdge*> boundaryEdgeList;
  const TCollection *coll = U_Space->GetCollection();
  coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);

  rhs_uD_n_v_n(rhs, U_Space, given_boundary_data1, given_boundary_data2,
          boundaryEdgeList, mult, rescale_by_h);
}

/**
 * @attention this functions assumes implicitly Matrix Type 14.
 * This means that the blocks are ordered like: A11,A12,B1t,A21,A22,B2t,B1,B2,C)
 * hence we need: blocks[0], blocks[1], blocks[3], blocks[4]
 * for A11, A12, A21, A22
 **/
void BoundaryAssembling2D::rhs_uD_n_v_n(BlockVector &rhs,
        const TFESpace2D *U_Space,
        BoundValueFunct2D *given_boundary_data1,
        BoundValueFunct2D *given_boundary_data2,
    std::vector<TBoundEdge*> &boundaryEdgeList,
    double mult,
    bool rescale_by_h)
{
  int ActiveBound = U_Space->get_n_active();

  /**
   * @todo: check if the matrix structure is correct:
   * we need 4 square matrices with the same FE spaces
   */

  for (size_t m = 0; m < boundaryEdgeList.size(); m++)
  {
    TBoundEdge *boundedge = boundaryEdgeList[m];
    TBaseCell *cell = boundedge->GetNeighbour(0);
    // get basis dimension and FE space data of cell i
    auto element = U_Space->get_fe(cell->GetCellIndex());

    int BaseVectDim = 1; // we assume only scalar FE
    int joint_id = boundedge->get_index_in_neighbour(cell);

    // get a quadrature formula good enough for the velocity FE space
    int fe_degree = element.GetBaseFunct()->GetPolynomialDegree();
    auto LineQuadFormula = QuadratureFormulaDatabase::qf_from_degree(
        2*fe_degree, BFRefElements::BFUnitLine);
    std::vector<double> quadWeights, quadPoints;
    get_quadrature_formula_data(quadPoints, quadWeights, *LineQuadFormula);

    //FEDatabase::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(
    //  this->LineQuadFormula);

    // compute values of all basis functions at all quadrature points
    std::vector< std::vector<double> > uorig, uxorig, uyorig;
    get_original_values(element, joint_id, cell, quadPoints, BaseVectDim, uorig,
                        uxorig, uyorig, *LineQuadFormula);

    double x_0, x_1, y_0, y_1;
    boundedge->get_vertices(x_0, y_0, x_1, y_1);
    // compute length of the edge
    double joint_length = boundedge->get_length();
    // normal vector to this boundary (normalized)
    double n1, n2;
    boundedge->get_normal(n1, n2);

    // quadrature
    for (unsigned int k = 0; k < quadPoints.size(); k++)
    {
      ///@attention in 1D the reference joint is [-1,1] => length = 2
      double reference_joint_length = 2;
      double x = x_0 + (quadPoints[k]+1.)/2. * (x_1-x_0);
      double y = y_0 + (quadPoints[k]+1.)/2. * (y_1-y_0);

      double T;
      boundedge->GetBoundComp()->GetTofXY(x, y, T);

      int BDComponent = boundedge->GetBoundComp()->GetID();

      // get the boundary values of rhs
      double value1, value2;
      if (given_boundary_data1 != nullptr)
      {
        given_boundary_data1(BDComponent, T, value1);
	
      }
      else
      {
        Output::print<1>("WARNING: Due to missing input (nullptr) in "
                "BoundaryAssembling2D::rhs_gradv_n_uD(), a default value for "
                "given_boundary_data1 is used, which might not fit your example.");
        value1 = 1;
      }

      if(given_boundary_data2 != nullptr)
      {
        given_boundary_data2(BDComponent, T, value2);
      }
      else
      {
        Output::print<1>("WARNING: Due to missing input (nullptr) in "
                "BoundaryAssembling2D::rhs_gradv_n_uD(), a default value for "
                "given_boundary_data2 is used, which might not fit your example.");
        value2 = 0.0;
      }

      // mapping from local(cell) DOF to global DOF
      const int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex()); //BeginIndex[i];

        for (unsigned int l = 0; l < uorig[k].size(); l++)
        {
          int global_dof_from_local = DOF[l];

          // if the DOF is Dirichlet, continue
          if (global_dof_from_local >= ActiveBound)
            continue;

          // updating rhs: int_gamma rhsval[2] v
          double v1 = uorig[k][l]; // value of test function (vtest = (vx,vy) )
          double v2 = v1;
          double u_n = value1 * n1 + value2 * n2;
          if (!rescale_by_h)
          {
            rhs.block(0)[global_dof_from_local] += mult * quadWeights[k] * u_n  * v1 * n1 *
                          (joint_length/reference_joint_length);
            rhs.block(1)[global_dof_from_local] += mult * quadWeights[k] * u_n  * v2 * n2 *
                                      (joint_length/reference_joint_length);
          }
          else
          {
            rhs.block(0)[global_dof_from_local] += mult * quadWeights[k] * u_n  * v1 * n1 *
                                      (joint_length/reference_joint_length) /joint_length;
           rhs.block(1)[global_dof_from_local] += mult * quadWeights[k] * u_n * v2 * n2 *
                                                  (joint_length/reference_joint_length) /joint_length;
          }
        }
    }
  }
}

//===============================================================================
void BoundaryAssembling2D::matrix_and_rhs_corner_stabilization(BlockFEMatrix &M,
        BlockVector &rhs,
        const TFESpace2D *U_Space,
        BoundValueFunct2D *given_boundary_data1,
        BoundValueFunct2D *given_boundary_data2,
        const std::vector<size_t>& nitsche_id,
        double mult)
{ 

  std::vector<TBoundEdge*> boundaryEdgeList;
  const TCollection *coll = U_Space->GetCollection();
  coll->get_boundary_edge_list(boundaryEdgeList);
  matrix_and_rhs_corner_stabilization(M,rhs, U_Space, given_boundary_data1,
          given_boundary_data2, boundaryEdgeList, nitsche_id, mult);

}

/**
 * @attention this functions assumes implicitely Matrix Type 14.
 * This means that the blocks are ordered like: A11,A12,B1t,A21,A22,B2t,B1,B2,C)
 * hence we need: blocks[0], blocks[1], blocks[3], blocks[4]
 * for A11, A12, A21, A22
 **/

void BoundaryAssembling2D::matrix_and_rhs_corner_stabilization(BlockFEMatrix &M,
        BlockVector &rhs,
        const TFESpace2D *U_Space,
        BoundValueFunct2D *given_boundary_data1,
        BoundValueFunct2D *given_boundary_data2,
        std::vector<TBoundEdge*> &boundaryEdgeList,
        const std::vector<size_t>& nitsche_id,
        double mult)
{
  //int ActiveBound = U_Space->GetN_ActiveDegrees();

  std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
  /**
   * @todo: check if the matrix structure is correct:
   * we need 4 square matrices with the same FE spaces
   */

  // find a 'pair' of corner edges
  for (size_t m = 0; m < boundaryEdgeList.size(); m++)
  {
    TBoundEdge *boundedge_1 = boundaryEdgeList[m];

    // Check if boundary edge is on boundary components of Nitsche type
    for (size_t km = 0; km < nitsche_id.size(); km++)
    {
      if (boundedge_1->GetBoundComp()->GetID() == (int)nitsche_id[km])
      {
        TBoundEdge *boundedge_2;

        for (size_t k = m+1; k < boundaryEdgeList.size(); k++)
        {
          boundedge_2 = boundaryEdgeList[k];

          for (size_t kl = 0; kl < nitsche_id.size(); kl++)
          {
            // Check if boundary edge is on boundary components of Nitsche type
            if (boundedge_2->GetBoundComp()->GetID() == (int)nitsche_id[kl])
            {
              // get normals
              double nE1_x, nE2_x, nE1_y, nE2_y;
              boundedge_1->get_normal(nE1_x, nE1_y);
              boundedge_2->get_normal(nE2_x, nE2_y);
              double abs_nE1_cdot_nE2 = std::abs(nE1_x * nE2_x + nE1_y * nE2_y);

              if ( std::abs( abs_nE1_cdot_nE2 - 1.) > 1e-12 )  // a.b=|a||b|std::cos(alpha) (not 0, 180, 360 degrees)
              {
                // get coordinates of relevant (Nitsche) corners xc,yc
                double x0_E1, y0_E1, x1_E1, y1_E1, x0_E2, y0_E2, x1_E2, y1_E2;
                boundedge_1->get_vertices(x0_E1, y0_E1, x1_E1, y1_E1);
                boundedge_2->get_vertices(x0_E2, y0_E2, x1_E2, y1_E2);

                //Check if the edges E1 and E2 share a vertex
                std::vector<double> xc, yc;
                if ( ( (std::abs(x0_E1 - x0_E2) <= 1.e-12) & (std::abs(y0_E1 - y0_E2) <= 1.e-12) ) ||
                     ( (std::abs(x0_E1 - x1_E2) <= 1.e-12) & (std::abs(y0_E1 - y1_E2) <= 1.e-12) ) )
                {
                  xc.push_back(x0_E1);
                  yc.push_back(y0_E1);
                  Output::print<4>("Detected a pair of corner edges with corner (xc,yc) = (",
                          xc[0], ", ", yc[0],").");

                  int locdof_corner_1, locdof_corner_2;
                  find_cornerDofs_in_boundarycells(xc, yc, U_Space, boundedge_1, boundedge_2,
                          locdof_corner_1, locdof_corner_2);

                  // mapping from local(cell) DOF to global DOF
                  TBaseCell *cell_1 = boundedge_1->GetNeighbour(0);
                  const int *DOF = U_Space->GetGlobalDOF(cell_1->GetCellIndex());
                  int test_DOF = DOF[locdof_corner_1];
                  int ansatz_DOF = test_DOF;

                  //see the note about blocks at the beginning of the function)
                  // In each corner point, there is exactly one basis function which is nonzero at (xc,yc).
                  // In fact it has the value (1,1) in (xc,yc).
                  blocks[0]->add( test_DOF, ansatz_DOF, mult * (nE1_x - nE2_x) * (nE1_x - nE2_x) ); // (1,0)^T * (n_E1-n_E2) --> A11
                  blocks[1]->add( test_DOF, ansatz_DOF, mult * (nE1_x - nE2_x) * (nE1_y - nE2_y) ); // A12
                  blocks[3]->add( test_DOF, ansatz_DOF, mult * (nE1_x - nE2_x) * (nE1_y - nE2_y) ); // A21
                  blocks[4]->add( test_DOF, ansatz_DOF, mult * (nE1_y - nE2_y) * (nE1_y - nE2_y)  ); // (0,1)^T * (n_E1-n_E2) --> A22

                  double value1 = 0.0, value2 = 0.0, T = 0.0;
                  if (given_boundary_data1 != nullptr)
                  {
                    boundedge_1->GetBoundComp()->GetTofXY(xc[0], yc[0], T);
                    given_boundary_data1(boundedge_1->GetBoundComp()->GetID(), T, value1);
		    
                  }
                  else
                  {
                    Output::print<1>("WARNING: Due to missing input (nullptr) in "
                            "BoundaryAssembling2D::matrix_and_rhs_corner_stabilization(), a default value"
                            " for given_boundary_data1 is used, which might not fit your example.");
                    value1 = 0.;//1;
                  }

                  if (given_boundary_data2 != nullptr)
                  {
                    boundedge_1->GetBoundComp()->GetTofXY(xc[0], yc[0], T);
                    given_boundary_data2(boundedge_1->GetBoundComp()->GetID(), T, value2);
                  }
                  else
                  {
                    Output::print<1>("WARNING: Due to missing input (nullptr) in "
                            "BoundaryAssembling2D::matrix_and_rhs_corner_stabilization(), a default value"
                            " for given_boundary_data2 is used, which might not fit your example.");
                    value2 = 0.;//1;
                  }

                  // add to both velocity components of the rhs
                  rhs.block(0)[test_DOF] += mult * ( value1*(nE1_x - nE2_x) + value2* (nE1_y - nE2_y) ) * (nE1_x - nE2_x);
                  rhs.block(1)[test_DOF] += mult * ( value1*(nE1_x - nE2_x) + value2* (nE1_y - nE2_y) ) * (nE1_y - nE2_y);
                }
                else if ( ( (std::abs(x1_E1 - x1_E2) <= 1.e-12) & (std::abs(y1_E1 - y1_E2) <= 1.e-12)) ||
                          ( (std::abs(x1_E1 - x0_E2) <= 1.e-12) & (std::abs(y1_E1 - y0_E2) <= 1.e-12) ))
                {
                  xc.push_back(x1_E1);
                  yc.push_back(y1_E1);
                  Output::print<4>("Detected a pair of corner edges with corner (xc,yc) = (", xc[0], ", ", yc[0],").");

                  int locdof_corner_1, locdof_corner_2;
                  find_cornerDofs_in_boundarycells(xc, yc, U_Space, boundedge_1, boundedge_2,
                          locdof_corner_1, locdof_corner_2);

                  // mapping from local(cell) DOF to global DOF
                  TBaseCell *cell_1 = boundedge_1->GetNeighbour(0);
                  const int *DOF = U_Space->GetGlobalDOF(cell_1->GetCellIndex());
                  int test_DOF = DOF[locdof_corner_1];
                  int ansatz_DOF = test_DOF;

                  //see the note about blocks at the beginning of the function)
                  // In each corner point, there is exactly one basis function which is nonzero at (xc,yc).
                  // In fact it has the value (1,1) in (xc,yc).
                  blocks[0]->add( test_DOF, ansatz_DOF, mult * (nE1_x - nE2_x) * (nE1_x - nE2_x) ); // (1,0)^T * (n_E1-n_E2) --> A11
                  blocks[1]->add( test_DOF, ansatz_DOF, mult * (nE1_x - nE2_x) * (nE1_y - nE2_y) ); // A12
                  blocks[3]->add( test_DOF, ansatz_DOF, mult * (nE1_x - nE2_x) * (nE1_y - nE2_y) ); // A21
                  blocks[4]->add( test_DOF, ansatz_DOF, mult * (nE1_y - nE2_y) * (nE1_y - nE2_y)  ); // (0,1)^T * (n_E1-n_E2) --> A22


                                    /*blocks[0]->add( test_DOF, ansatz_DOF, mult * (nE1_1 - nE2_1) * (nE1_1 - nE2_1) ); // (1,0)^T * (n_E1-n_E2) --> A11
                  blocks[1]->add( test_DOF, ansatz_DOF, mult * (nE1_2 - nE2_2) * (nE1_1 - nE2_1) ); // A12
                  blocks[3]->add( test_DOF, ansatz_DOF, mult * (nE1_1 - nE2_1) * (nE1_2 - nE2_2) ); // A21
                  blocks[4]->add( test_DOF, ansatz_DOF, mult * (nE1_2 - nE2_2) * (nE1_2 - nE2_2) ); // (0,1)^T * (n_E1-n_E2) --> A22
*/
                  double value1 = 0.0, value2 = 0.0, T = 0.0;

                  if (given_boundary_data1 != nullptr)
                  {
                    boundedge_1->GetBoundComp()->GetTofXY(xc[0], yc[0], T);
                    given_boundary_data1(boundedge_1->GetBoundComp()->GetID(), T, value1);
                  }
                  else
                  {
                    Output::print<1>("WARNING: Due to missing input (nullptr) in "
                            "BoundaryAssembling2D::matrix_and_rhs_corner_stabilization(), a default value"
                            " for given_boundary_data1 is used, which might not fit your example.");
                    value1 = 0.;//1;
                  }

                  if (given_boundary_data2 != nullptr)
                  {
                    boundedge_1->GetBoundComp()->GetTofXY(xc[0], yc[0], T);
                    given_boundary_data2(boundedge_1->GetBoundComp()->GetID(), T, value2);
                  }
                  else
                  {
                    Output::print<1>("WARNING: Due to missing input (nullptr) in "
                            "BoundaryAssembling2D::matrix_and_rhs_corner_stabilization(), a default value"
                            " for given_boundary_data2 is used, which might not fit your example.");
                    value2 = 0.;//1;
                  }

                  // add to both velocity components of the rhs
                  rhs.block(0)[test_DOF] += mult * ( value1*(nE1_x - nE2_x) + value2*(nE1_y - nE2_y) ) * (nE1_x - nE2_x);
                  rhs.block(1)[test_DOF] += mult * ( value1*(nE1_x - nE2_x) + value2*(nE1_y - nE2_y) ) * (nE1_y - nE2_y);
                }
              }
            }
          }
        }
      }
    }
  }
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void BoundaryAssembling2D::find_cornerDofs_in_boundarycells(
    const std::vector<double>& xc, const std::vector<double>& yc,
    const TFESpace2D *U_Space,
    TBoundEdge *boundedge_1, TBoundEdge *boundedge_2, 
    int &locdof_corner_1, int &locdof_corner_2)
{
  TBaseCell *cell_1 = boundedge_1->GetNeighbour(0);
  TBaseCell *cell_2 = boundedge_2->GetNeighbour(0);

  int N_BaseFunct = U_Space->get_n_local_dof(cell_1->GetCellIndex());

  for (int locdof = 0; locdof < N_BaseFunct; locdof++)
  {
    for (unsigned int i = 0; i < xc.size(); i++)
    {
      const int *DOF = U_Space->GetGlobalDOF(cell_1->GetCellIndex());
      int globaldof = DOF[locdof];
      double x, y; 
      U_Space->GetDOFPosition(globaldof, x, y);
      if ( (x == xc[i]) & (y == yc[i]) )
      { 
        locdof_corner_1 = locdof;
        break;
      }
    }
  }
  //Output::print("number of cell_1: ", cell_1->GetCellIndex());
  //Output::print("number of cell_2: ", cell_2->GetCellIndex());

  if (cell_1->GetCellIndex()  == cell_2->GetCellIndex() )
  {
    locdof_corner_2 = locdof_corner_1;
    Output::print("WARNING in BoundaryAssembling2D::find_cornerDofs_in_boundarycells(): "
            "Detected a single corner cell. The theory for the corner stabilization does not "
            "allow corner nodes that are solely contained in one cell. Decompose them!");
  }
  else
  {
    int N_BaseFunct = U_Space->get_n_local_dof(cell_2->GetCellIndex());

    for (int locdof = 0; locdof < N_BaseFunct; locdof++)
    {
      for (unsigned int i = 0; i < xc.size(); i++)
      {
        const int *DOF = U_Space->GetGlobalDOF(cell_2->GetCellIndex());
        int globaldof = DOF[locdof];
        double x, y;
        U_Space->GetDOFPosition(globaldof, x, y);
        if ( (x == xc[i]) & (y == yc[i]))
        { 
          locdof_corner_2 = locdof;
          break;
        }
      }
    }
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// <d(u.n)/dtau, d(v.n)/dtau> at corners of the domain
void BoundaryAssembling2D::matrix_gradu_n_t_gradv_n_t(BlockFEMatrix &M,
    const TFESpace2D *U_Space,
    int boundary_component_id,
    double mult)
{
  std::vector<TBoundEdge*> boundaryEdgeList;
  const TCollection *coll = U_Space->GetCollection();
  coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);

  matrix_gradu_n_t_gradv_n_t(M, U_Space, boundaryEdgeList, mult);
}
/**
 * @attention this functions assumes implicitely Matrix Type 14.
 * This means that the blocks are ordered like: A11,A12,B1t,A21,A22,B2t,B1,B2,C)
 * hence we need: blocks[0], blocks[1], blocks[3], blocks[4]
 * for A11, A12, A21, A22
 **/
void BoundaryAssembling2D::matrix_gradu_n_t_gradv_n_t(BlockFEMatrix &M,
    const TFESpace2D *U_Space,
    std::vector<TBoundEdge*> &boundaryEdgeList,
    double mult)
{
  int ActiveBound = U_Space->get_n_active();

  std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
  /**
   * @todo: check if the matrix structure is correct:
   * we need 4 square matrices with the same FE spaces
   */

  for (size_t m = 0; m < boundaryEdgeList.size(); m++)
  {
    TBoundEdge *boundedge = boundaryEdgeList[m];
    TBaseCell *cell = boundedge->GetNeighbour(0);
    // get basis dimension and FE space data of cell i
    auto element = U_Space->get_fe(cell->GetCellIndex());

    int BaseVectDim = 1; // we assume only scalar FE
    int joint_id = boundedge->get_index_in_neighbour(cell);

    // get a quadrature formula good enough for the velocity FE space
    int fe_degree = element.GetBaseFunct()->GetPolynomialDegree();
    auto LineQuadFormula = QuadratureFormulaDatabase::qf_from_degree(
        2*fe_degree, BFRefElements::BFUnitLine);
    std::vector<double> quadWeights, quadPoints;
    get_quadrature_formula_data(quadPoints, quadWeights, *LineQuadFormula);

    //FEDatabase::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(
    //  this->LineQuadFormula);

    // compute values of all basis functions at all quadrature points
    std::vector< std::vector<double> > uorig, uxorig, uyorig;
    get_original_values(element, joint_id, cell, quadPoints, BaseVectDim, uorig,
                        uxorig, uyorig, *LineQuadFormula);

    double x0, x1, y0, y1;
    boundedge->get_vertices(x0, y0, x1, y1);
    // compute length of the edge
    double joint_length = boundedge->get_length();
    // normal vector to this boundary (normalized)
    double t1, t2, n1, n2;
    boundedge->get_tangent(t1, t2);
    boundedge->get_normal(n1, n2);

    // quadrature
    for (unsigned int k = 0; k < quadPoints.size(); k++)
    {
      ///@attention in 1D the reference joint is [-1,1] => length = 2
      double reference_joint_length = 2;

      // mapping from local (cell) DOF to global DOF
      const int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex()); //BeginIndex[i];

      // loop on test functions
      double scale_factor = mult * quadWeights[k] * (joint_length/reference_joint_length);
      for (unsigned int l1 = 0; l1 < uorig[k].size(); l1++)
      {
        int test_DOF = DOF[l1];

        // if the DOF is Dirichlet, continue
        if(test_DOF >= ActiveBound)
          continue;

        double v1_dx = uxorig[k][l1];
        double v1_dy = uyorig[k][l1];
        double v2_dx = v1_dx; // v1 and v2 component have the same FE space
        double v2_dy = v1_dy; // v1 and v2 component have the same FE space

        // loop over ansatz functions
        for (unsigned int l2 = 0; l2 < uorig[k].size(); l2++)
        {
          int ansatz_DOF = DOF[l2];
          double u1_dx = uxorig[k][l2];
          double u1_dy = uyorig[k][l2];
          double u2_dx = u1_dx; // u1 and u2 component have the same FE space
          double u2_dy = u1_dy; // u1 and u2 component have the same FE space

          // (see the note about blocks at the beginning of the function)
          blocks[0]->add( test_DOF, ansatz_DOF, scale_factor * ( (u1_dx * n1 * t1) * (v1_dx * n1 * t1) 
                + (u1_dx * n1 * t1) * (v1_dy * n1 * t2)
                + (u1_dy * n1 * t2) * (v1_dx * n1 * t1)
                + (u1_dy * n1 * t2) * (v1_dy * n1 * t2)) ); // A11
          blocks[1]->add( test_DOF, ansatz_DOF, scale_factor * ( (u2_dx * n2 * t1) * (v1_dx * n1 * t1) 
                + (u2_dx * n2 * t1) * (v1_dy * n1 * t2)
                + (u2_dy * n2 * t2) * (v1_dx * n1 * t1)
                + (u2_dy * n2 * t2) * (v1_dy * n1 * t2)) ); // A12
          blocks[3]->add( test_DOF, ansatz_DOF, scale_factor * ( (u1_dx * n1 * t1) * (v2_dx * n2 * t1) 
                + (u1_dx * n1 * t1) * (v2_dy * n2 * t2)
                + (u1_dy * n1 * t2) * (v2_dx * n2 * t1)
                + (u1_dy * n1 * t2) * (v2_dy * n2 * t2)) );// A21
          blocks[4]->add( test_DOF, ansatz_DOF, scale_factor * (  (u2_dx * n2 * t1) * (v2_dx * n2 * t1) 
                + (u2_dx * n2 * t1) * (v2_dy * n2 * t2)
                + (u2_dy * n2 * t2) * (v2_dx * n2 * t1)
                + (u2_dy * n2 * t2) * (v2_dy * n2 * t2)) ); // A22
        }
      } 
    }
  } 
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void BoundaryAssembling2D::matrix_gradu_n_v(BlockFEMatrix &M,
    const TFESpace2D *U_Space,
    int boundary_component_id,
    double mult)
{
  std::vector<TBoundEdge*> boundaryEdgeList;
  const TCollection *coll = U_Space->GetCollection();
  coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);
  matrix_gradu_n_v(M, U_Space, boundaryEdgeList, mult, boundary_component_id);
}
/**
 * @attention this functions assumes implicitely Matrix Type 14.
 * This means that the blocks are ordered like: A11,A12,B1t,A21,A22,B2t,B1,B2,C)
 * hence we need: blocks[0], blocks[1], blocks[3], blocks[4]
 * for A11, A12, A21, A22
 **/
void BoundaryAssembling2D::matrix_gradu_n_v(BlockFEMatrix &M,
    const TFESpace2D *U_Space,
    std::vector<TBoundEdge*> &boundaryEdgeList,
    double mult,
    int)
{
  int ActiveBound = U_Space->get_n_active();

  std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
  /**
   * @todo: check if the matrix structure is correct:
   * we need 4 square matrices with the same FE spaces
   */

  for (size_t m = 0; m < boundaryEdgeList.size(); m++)
  {
    TBoundEdge *boundedge = boundaryEdgeList[m];
    TBaseCell *cell = boundedge->GetNeighbour(0);
    // get basis dimension and FE space data of cell i
    auto element = U_Space->get_fe(cell->GetCellIndex());

    int BaseVectDim = 1; // we assume only scalar FE
    int joint_id = boundedge->get_index_in_neighbour(cell);

    // get a quadrature formula good enough for the velocity FE space
    int fe_degree = element.GetBaseFunct()->GetPolynomialDegree();
    auto LineQuadFormula = QuadratureFormulaDatabase::qf_from_degree(
        2*fe_degree, BFRefElements::BFUnitLine);

    std::vector<double> quadWeights, quadPoints;
    get_quadrature_formula_data(quadPoints, quadWeights, *LineQuadFormula);

    //FEDatabase::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(
    //  this->LineQuadFormula);

    // compute values of all basis functions at all quadrature points
    std::vector< std::vector<double> > uorig, uxorig, uyorig;
    get_original_values(element, joint_id, cell, quadPoints, BaseVectDim, uorig,
                        uxorig, uyorig, *LineQuadFormula);

    double x0, x1, y0, y1;
    boundedge->get_vertices(x0,  y0, x1, y1);
    // compute length of the edge
    double joint_length = boundedge->get_length();
    // normal vector to this boundary (normalized)
    double n1, n2;
    boundedge->get_normal(n1, n2);

    // quadrature
    for (unsigned int k = 0; k < quadPoints.size(); k++)
    {
      ///@attention in 1D the reference joint is [-1,1] => length = 2
      double reference_joint_length = 2;

      // mapping from local(cell) DOF to global DOF
      const int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex()); //BeginIndex[i];

      // loop on test functions
      double scale_factor = mult * quadWeights[k] * (joint_length/reference_joint_length);
      for(unsigned int l1 = 0; l1 < uorig[k].size(); l1++)
      {
        int test_DOF = DOF[l1];

        // if the DOF is Dirichlet, continue
        if(test_DOF >= ActiveBound)
          continue;

        double v1 = uorig[k][l1];
        double v2 = v1; // x and y component have the same FE space

        // loop on ansatz functions
        for(unsigned int l2 = 0; l2 < uorig[k].size(); l2++)
        {
          int ansatz_DOF = DOF[l2];

          double u1_dx = uxorig[k][l2];
          double u1_dy = uyorig[k][l2];
          double u2_dx = u1_dx; // u1 and u2 component have the same FE space
          double u2_dy = u1_dy; // u1 and u2 component have the same FE space

          // (see the note about blocks at the beginning of the function)
          blocks[0]->add( test_DOF, ansatz_DOF, scale_factor * v1 * u1_dx * n1 ); // A11
          blocks[0]->add( test_DOF, ansatz_DOF, scale_factor * v1 * u1_dy * n2 ); // A11
          blocks[4]->add(test_DOF, ansatz_DOF, scale_factor * v2 * u2_dy * n2 ); // A22
          blocks[4]->add(test_DOF, ansatz_DOF, scale_factor * v2 * u2_dx * n1 ); // A22
        }
      }
    }
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void BoundaryAssembling2D::matrix_gradv_n_u(BlockFEMatrix &M,
    const TFESpace2D *U_Space,
    int boundary_component_id,
    double mult)
{
  std::vector<TBoundEdge*> boundaryEdgeList;
  const TCollection *coll = U_Space->GetCollection();
  coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);
  matrix_gradv_n_u(M, U_Space, boundaryEdgeList, mult);
}
/**
 * @attention this functions assumes implicitely Matrix Type 14.
 * This means that the blocks are ordered like: A11,A12,B1t,A21,A22,B2t,B1,B2,C)
 * hence we need: blocks[0], blocks[1], blocks[3], blocks[4]
 * for A11, A12, A21, A22
 **/
void BoundaryAssembling2D::matrix_gradv_n_u(BlockFEMatrix &M,
    const TFESpace2D *U_Space,
    std::vector<TBoundEdge*> &boundaryEdgeList,
    double mult)
{
  int ActiveBound = U_Space->get_n_active();

  std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
  /**
   * @todo: check if the matrix structure is correct:
   * we need 4 square matrices with the same FE spaces
   */

  for(size_t m = 0; m < boundaryEdgeList.size(); m++)
  {
    TBoundEdge *boundedge = boundaryEdgeList[m];
    TBaseCell *cell = boundedge->GetNeighbour(0);
    // get basis dimension and FE space data of cell i
    auto element = U_Space->get_fe(cell->GetCellIndex());

    int BaseVectDim = 1; // we assume only scalar FE
    int joint_id = boundedge->get_index_in_neighbour(cell);

    // get a quadrature formula good enough for the velocity FE space
    int fe_degree = element.GetBaseFunct()->GetPolynomialDegree();
    auto LineQuadFormula = QuadratureFormulaDatabase::qf_from_degree(
        2*fe_degree, BFRefElements::BFUnitLine);
    std::vector<double> quadWeights, quadPoints;
    get_quadrature_formula_data(quadPoints, quadWeights, *LineQuadFormula);

    //FEDatabase::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(
    //  this->LineQuadFormula);

    // compute values of all basis functions at all quadrature points
    std::vector< std::vector<double> > uorig, uxorig, uyorig;
    get_original_values(element, joint_id, cell, quadPoints, BaseVectDim, uorig,
                        uxorig, uyorig, *LineQuadFormula);

    double x0, x1, y0, y1;
    boundedge->get_vertices(x0, y0, x1, y1);
    // compute length of the edge
    double joint_length = boundedge->get_length();
    // normal vector to this boundary (normalized)
    double n1, n2;
    boundedge->get_normal(n1, n2);

    // quadrature
    for(unsigned int k = 0; k < quadPoints.size(); k++)
    {
      ///@attention in 1D the reference joint is [-1,1] => length = 2
      double reference_joint_length = 2;

      // mapping from local(cell) DOF to global DOF
      const int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex()); //BeginIndex[i];

      // loop on test functions
      double scale_factor = mult * quadWeights[k] * (joint_length/reference_joint_length);
      for (unsigned int l1 = 0; l1 < uorig[k].size(); l1++)
      {
        int test_DOF = DOF[l1];

        // if the DOF is Dirichlet, continue
        if(test_DOF >= ActiveBound)
          continue;

        double v1_dx = uxorig[k][l1];
        double v1_dy = uyorig[k][l1];
        double v2_dx = v1_dx;
        double v2_dy = v1_dy; // x and y component have the same FE space

        // loop on ansatz functions
        for(unsigned int l2 = 0; l2 < uorig[k].size(); l2++)
        {
          int ansatz_DOF = DOF[l2];

          double u1 = uorig[k][l2];
          double u2 = u1; // x and y component have the same FE space

          // (see the note about blocks at the beginning of the function)
          blocks[0]->add( test_DOF, ansatz_DOF, scale_factor * u1 * v1_dx * n1 ); // A11
          blocks[0]->add( test_DOF, ansatz_DOF, scale_factor * u1 * v1_dy * n2 ); // A11
          blocks[4]->add( test_DOF, ansatz_DOF, scale_factor * u2 * v2_dx * n1 ); // A22
          blocks[4]->add( test_DOF, ansatz_DOF, scale_factor * u2 * v2_dy * n2 ); // A22
        }
      }
    }
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void BoundaryAssembling2D::matrix_u_v(BlockFEMatrix &M,
    const TFESpace2D *U_Space,
    int boundary_component_id,
    double mult,
    bool rescale_by_h)
{
  std::vector<TBoundEdge*> boundaryEdgeList;
  const TCollection *coll = U_Space->GetCollection();
  coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);
  matrix_u_v(M, U_Space, boundaryEdgeList, mult, rescale_by_h);
}
/**
 * @attention this functions assumes implicitely Matrix Type 14.
 * This means that the blocks are ordered like: A11, A12, B1T, A21, A22, B2T, B1, B2, C)
 * hence we need: blocks[0], blocks[1], blocks[3], blocks[4]
 * for A11, A12, A21, A22
 **/
void BoundaryAssembling2D::matrix_u_v(BlockFEMatrix &M,
    const TFESpace2D *U_Space,
    std::vector<TBoundEdge*> &boundaryEdgeList,
    double mult,
    bool rescale_by_h)
{
  int ActiveBound = U_Space->get_n_active();

  std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
  /**
   * @todo: check if the matrix structure is correct:
   * we need 4 square matrices with the same FE spaces
   */

  for(size_t m = 0; m < boundaryEdgeList.size(); m++)
  {
    TBoundEdge *boundedge = boundaryEdgeList[m];
    TBaseCell *cell = boundedge->GetNeighbour(0);
    // get basis dimension and FE space data of cell i
    auto element = U_Space->get_fe(cell->GetCellIndex());

    int BaseVectDim = 1; // we assume only scalar FE
    int joint_id = boundedge->get_index_in_neighbour(cell);

    // get a quadrature formula good enough for the velocity FE space
    int fe_degree = element.GetBaseFunct()->GetPolynomialDegree();
    auto LineQuadFormula = QuadratureFormulaDatabase::qf_from_degree(
        2*fe_degree, BFRefElements::BFUnitLine);
    std::vector<double> quadWeights, quadPoints;
    get_quadrature_formula_data(quadPoints, quadWeights, *LineQuadFormula);

    //FEDatabase::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(
    //  this->LineQuadFormula);

    // compute values of all basis functions at all quadrature points
    std::vector< std::vector<double> > uorig, uxorig, uyorig;
    get_original_values(element, joint_id, cell, quadPoints, BaseVectDim, uorig,
                        uxorig, uyorig, *LineQuadFormula);

    double x0, x1, y0, y1;
    boundedge->get_vertices(x0, y0, x1, y1);
    // compute length of the edge
    double joint_length = boundedge->get_length();

    // quadrature
    for(unsigned int k = 0; k < quadPoints.size(); k++)
    {
      ///@attention in 1D the reference joint is [-1,1] => length = 2
      double reference_joint_length = 2;

      // mapping from local(cell) DOF to global DOF
      const int *DOF = U_Space->GetGlobalDOF(cell->GetCellIndex()); //BeginIndex[i];

      // rescale local integral (Nitsche)
      double scale_factor;
      if (rescale_by_h)
      {
        scale_factor = ( mult * quadWeights[k] * (joint_length/reference_joint_length)) /joint_length;
      }
      else
      {
        scale_factor = mult * quadWeights[k] * (joint_length/reference_joint_length);
      }

      // loop on test functions
      for(unsigned int l1 = 0; l1 < uorig[k].size(); l1++)
      {
        int test_DOF = DOF[l1];

        // if the DOF is Dirichlet, continue
        if(test_DOF >= ActiveBound)
          continue;

        double v1 = uorig[k][l1];
        double v2 = v1; // x and y component have the same FE space

        // loop on ansatz functions
        for(unsigned int l2 = 0; l2 < uorig[k].size(); l2++)
        {
          int ansatz_DOF = DOF[l2];
          double u1 = uorig[k][l2];
          double u2 = u1; // x and y component have the same FE space

          // (see the note about blocks at the beginning of the function)
          blocks[0]->add( test_DOF, ansatz_DOF, scale_factor * u1 * v1 ); // A11
          blocks[4]->add( test_DOF, ansatz_DOF, scale_factor * u2 * v2 ); // A22
        }
      }
    }
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// (the edge list can be generated outside this class)
void BoundaryAssembling2D::matrix_q_u_n(BlockFEMatrix &M,
    const TFESpace2D *U_Space,
    const TFESpace2D *P_Space,
    int boundary_component_id,
    double mult)
{
  std::vector<TBoundEdge*> boundaryEdgeList;
  const TCollection *coll = U_Space->GetCollection();
  coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);
  matrix_q_u_n(M, U_Space, P_Space, boundaryEdgeList, mult);
}

void BoundaryAssembling2D::matrix_q_u_n(BlockFEMatrix &M,
    const TFESpace2D *U_Space,
    const TFESpace2D *P_Space,
    std::vector<TBoundEdge*> &boundaryEdgeList,
    double mult)
{
  int ActiveBound = U_Space->get_n_active();

  std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
  /**
   * @todo: check if the matrix structure is correct:
   * we need 4 square matrices with the same FE spaces
   */

  for(size_t m = 0; m < boundaryEdgeList.size(); m++)
  {
    TBoundEdge *boundedge = boundaryEdgeList[m];
    TBaseCell *cell = boundedge->GetNeighbour(0);
    // get basis dimension and FE space data of cell i
    auto element_u = U_Space->get_fe(cell->GetCellIndex());
    // get basis dimension and FE space data of cell i
    auto element_p = P_Space->get_fe(cell->GetCellIndex());

    int BaseVectDim = 1; // we assume only scalar FE
    int joint_id = boundedge->get_index_in_neighbour(cell);

    // get a quadrature formula good enough for the velocity FE space
    int fe_degree_U = element_u.GetBaseFunct()->GetPolynomialDegree();
    int fe_degree_P = element_p.GetBaseFunct()->GetPolynomialDegree();
    auto LineQuadFormula = QuadratureFormulaDatabase::qf_from_degree(
        fe_degree_P+fe_degree_U, BFRefElements::BFUnitLine);
    std::vector<double> quadWeights, quadPoints;
    get_quadrature_formula_data(quadPoints, quadWeights, *LineQuadFormula);

    //FEDatabase::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(
    //   this->LineQuadFormula);

    // compute values of all basis functions at all quadrature points
    std::vector<std::vector<double>> uorig, uxorig, uyorig;
    get_original_values(element_u, joint_id, cell, quadPoints, BaseVectDim,
                        uorig, uxorig, uyorig, *LineQuadFormula);

    int BaseVectDim_P = 1; // we assume only scalar FE; nur bei Raviart-Thomas & BDM \neq 1

    // compute values of all basis functions at all quadrature points
    std::vector<std::vector<double>> porig, pxorig, pyorig;
    get_original_values(element_p, joint_id, cell, quadPoints, BaseVectDim_P,
                        porig, pxorig, pyorig, *LineQuadFormula);

    double x0, x1, y0, y1;
    boundedge->get_vertices(x0, y0, x1, y1);
    // compute length of the edge
    double joint_length = boundedge->get_length();
    // normal vector to this boundary (normalized)
    double n1, n2;
    boundedge->get_normal(n1, n2);

    // mapping from local(cell) DOF to global DOF
    // int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];
    const int *DOF_P = P_Space->GetGlobalDOF(cell->GetCellIndex());
    const int *DOF_U = U_Space->GetGlobalDOF(cell->GetCellIndex());

    // quadrature
    for( unsigned int k = 0; k < quadPoints.size(); k++ )
    {
      ///@attention in 1D the reference joint is [-1,1] => length = 2
      double reference_joint_length = 2;

      // loop on test functions
      double scale_factor = mult * quadWeights[k] * (joint_length/reference_joint_length);

      for( unsigned int l1 = 0; l1 < uorig[k].size(); l1++ )
      {
        int ansatz_DOF = DOF_U[l1];

        double u1 = uorig[k][l1];
        double u2 = u1; // x and y component have the same FE space

        // loop on ansatz functions
        for(unsigned int l2 = 0; l2 < porig[k].size(); l2++)
        {
          int test_DOF = DOF_P[l2];
          // if the DOF is Dirichlet, continue
          if(test_DOF >= ActiveBound)
            continue;
          double q = porig[k][l2];
          blocks[6]->add( test_DOF, ansatz_DOF, scale_factor * q * u1 * n1 ); // B1
          blocks[7]->add( test_DOF, ansatz_DOF, scale_factor * q * u2 * n2 ); // B2
        }
      }
    } 
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void BoundaryAssembling2D::rhs_q_uD_n(BlockVector &rhs,
    const TFESpace2D *U_Space,
    const TFESpace2D *P_Space,
    BoundValueFunct2D *given_boundary_data1,
    BoundValueFunct2D *given_boundary_data2,
    int boundary_component_id,
    double mult)
{
  std::vector<TBoundEdge*> boundaryEdgeList;
  const TCollection *coll = P_Space->GetCollection();
  coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);

  rhs_q_uD_n(rhs, U_Space, P_Space, given_boundary_data1, given_boundary_data2, boundaryEdgeList, mult);
}

void BoundaryAssembling2D::rhs_q_uD_n(BlockVector &rhs,
    const TFESpace2D *U_Space,
    const TFESpace2D *P_Space,
    BoundValueFunct2D *given_boundary_data1,
    BoundValueFunct2D *given_boundary_data2,
    std::vector<TBoundEdge*> &boundaryEdgeList,
    double mult)
{
  // =========================================
  int ActiveBound = P_Space->get_n_active();

  for(size_t m = 0; m < boundaryEdgeList.size(); m++)
  {
    TBoundEdge *boundedge = boundaryEdgeList[m];
    TBaseCell *cell = boundedge->GetNeighbour(0);

    int BaseVectDim_P = 1; // we assume only scalar FE
    int joint_id = boundedge->get_index_in_neighbour(cell);
    auto element_u = U_Space->get_fe(cell->GetCellIndex());
    // get basis dimension and FE space data of cell i
    auto element_p = P_Space->get_fe(cell->GetCellIndex());

    // get a quadrature formula good enough for the velocity FE space
    int fe_degree_U = element_u.GetBaseFunct()->GetPolynomialDegree();
    int fe_degree_P = element_p.GetBaseFunct()->GetPolynomialDegree();

    auto LineQuadFormula = QuadratureFormulaDatabase::qf_from_degree(
        fe_degree_U + fe_degree_P, BFRefElements::BFUnitLine);
    std::vector<double> quadWeights, quadPoints;
    get_quadrature_formula_data(quadPoints, quadWeights, *LineQuadFormula);

    // compute values of all basis functions at all quadrature points
    std::vector< std::vector<double> > porig, pxorig, pyorig;
    get_original_values(element_p, joint_id, cell, quadPoints, BaseVectDim_P,
                        porig, pxorig, pyorig, *LineQuadFormula);

    double x_0, x_1, y_0, y_1;
    boundedge->get_vertices(x_0, y_0, x_1, y_1);
    // compute length of the edge
    double joint_length = boundedge->get_length();
    // normal vector to this boundary (normalized)
    double n1,n2; // n=(n_1, n_2)^T
    boundedge->get_normal(n1, n2);

    // quadrature
    for(unsigned int k = 0; k < quadPoints.size(); k++)
    {
      ///@attention in 1D the reference joint is [-1,1] => length = 2
      double reference_joint_length = 2;
      double x = x_0 + (quadPoints[k]+1.)/2.*(x_1-x_0);
      double y = y_0 + (quadPoints[k]+1.)/2.*(y_1-y_0);

      double T;
      boundedge->GetBoundComp()->GetTofXY(x, y, T);

      int BDComponent = boundedge->GetBoundComp()->GetID();

      // get the boundary values of rhs
      double value1, value2;
      if(given_boundary_data1 != nullptr)
      {
        given_boundary_data1(BDComponent, T, value1);
      }
      else
      {
        Output::print<1>("WARNING: Due to missing input (nullptr) in "
            "BoundaryAssembling2D::rhs_q_uD_n(), a default value for "
            "given_boundary_data1 is used, which might not fit your example.");
        value1 = 1.;
      }

      if(given_boundary_data2 != nullptr)
      {
        given_boundary_data2(BDComponent, T, value2);
      }
      else
      {
        Output::print<1>("WARNING: Due to missing input (nullptr) in "
            "BoundaryAssembling2D::rhs_q_uD_n(), a default value for "
            "given_boundary_data2 is used, which might not fit your example.");
        value2 = 1.;//0.0;
      }

      // mapping from local (cell) DOF to global DOF
      const int *DOF = P_Space->GetGlobalDOF(cell->GetCellIndex()); //BeginIndex[i];

      for(unsigned int l = 0; l < porig[k].size(); l++)
      {
        int global_dof_from_local = DOF[l];

        // if the DOF is Dirichlet, continue
        if(global_dof_from_local >= ActiveBound)
          continue;

        // updating rhs: int_gamma rhsval v \cdot n
        double q = porig[k][l]; // value of test function (vtest = vx = vy)
        // add for both components
        rhs.block(2)[global_dof_from_local] += mult * quadWeights[k] * q * (value1 * n1) *
          (joint_length/reference_joint_length);
        rhs.block(2)[global_dof_from_local] += mult * quadWeights[k] * q * (value2 * n2) *
          (joint_length/reference_joint_length);
      }   
    }
  } 
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void BoundaryAssembling2D::matrix_p_v_n(BlockFEMatrix &M,
    const TFESpace2D *U_Space,
    const TFESpace2D *P_Space,
    int boundary_component_id,
    double mult)
{
  std::vector<TBoundEdge*> boundaryEdgeList;
  const TCollection *coll = U_Space->GetCollection();
  coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);
  matrix_p_v_n(M,U_Space, P_Space, boundaryEdgeList, mult);
}

void BoundaryAssembling2D::matrix_p_v_n(BlockFEMatrix &M,
    const TFESpace2D *U_Space,
    const TFESpace2D *P_Space,
    std::vector<TBoundEdge*> &boundaryEdgeList,
    double mult)
{
  int ActiveBound = U_Space->get_n_active();

  std::vector<std::shared_ptr<FEMatrix>> blocks = M.get_blocks_uniquely();
  /**
   * @todo: check if the matrix structure is correct:
   * we need 4 square matrices with the same FE spaces
   */

  for(size_t m = 0; m < boundaryEdgeList.size(); m++)
  {
    TBoundEdge *boundedge = boundaryEdgeList[m];
    TBaseCell *cell = boundedge->GetNeighbour(0);
    // get basis dimension and FE space data of cell i
    auto element_u = U_Space->get_fe(cell->GetCellIndex());
    // get basis dimension and FE space data of cell i
    auto element_p = P_Space->get_fe(cell->GetCellIndex());

    int BaseVectDim = 1; // we assume only scalar FE
    int joint_id = boundedge->get_index_in_neighbour(cell);

    // get a quadrature formula good enough for the velocity FE space
    int fe_degree_U = element_u.GetBaseFunct()->GetPolynomialDegree();
    int fe_degree_P = element_p.GetBaseFunct()->GetPolynomialDegree();
    auto LineQuadFormula = QuadratureFormulaDatabase::qf_from_degree(
        fe_degree_P*fe_degree_U, BFRefElements::BFUnitLine);
    std::vector<double> quadWeights, quadPoints;
    get_quadrature_formula_data(quadPoints, quadWeights, *LineQuadFormula);

    //FEDatabase::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(
    //  this->LineQuadFormula);

    // compute values of all basis functions at all quadrature points
    std::vector< std::vector<double> > uorig, uxorig, uyorig;
    get_original_values(element_u, joint_id, cell, quadPoints, BaseVectDim,
                        uorig, uxorig, uyorig, *LineQuadFormula);

    int BaseVectDim_P = 1; // we assume only scalar FE; nur bei Raviart-Thomas & BDM \neq 1

    // compute values of all basis functions at all quadrature points
    std::vector< std::vector<double>> porig, pxorig, pyorig;
    get_original_values(element_p, joint_id, cell, quadPoints, BaseVectDim_P,
                        porig, pxorig, pyorig, *LineQuadFormula);

    double x0, x1, y0, y1;
    boundedge->get_vertices(x0, y0, x1, y1);
    // compute length of the edge
    double joint_length = boundedge->get_length();
    // normal vector to this boundary (normalized)
    double n1, n2;
    boundedge->get_normal(n1, n2);

    // mapping from local(cell) DOF to global DOF
    // int *DOF = GlobalNumbers + BeginIndex[cell->GetCellIndex()]; //BeginIndex[i];
    const int *DOF_P=P_Space->GetGlobalDOF(cell->GetCellIndex());
    const int *DOF_U=U_Space->GetGlobalDOF(cell->GetCellIndex());

    // quadrature
    for(unsigned int k = 0; k < quadPoints.size(); k++)
    {
      ///@attention in 1D the reference joint is [-1,1] => length = 2
      double reference_joint_length = 2;

      // loop on test functions
      double scale_factor = mult * quadWeights[k] * (joint_length/reference_joint_length);
      for(unsigned int l1 = 0 ; l1 < uorig[k].size() ; l1++)
      {
        int test_DOF = DOF_U[l1];

        // if the DOF is Dirichlet, continue
        if(test_DOF >= ActiveBound)
          continue;

        double v1 = uorig[k][l1];
        double v2 = v1; // x and y component have the same FE space

        // loop on ansatz functions
        for(unsigned int l2 = 0; l2 < porig[k].size(); l2++)
        {
          int ansatz_DOF = DOF_P[l2];
          double p = porig[k][l2];
          // ToDo: Check if block reference is correct:
          // (see the note about blocks at the beginning of the function)
          blocks[2]->add( test_DOF, ansatz_DOF, scale_factor * p * v1 * n1 ); // B1
          blocks[5]->add( test_DOF, ansatz_DOF, scale_factor * p * v2 * n2 ); // B2
        }  
      }
    }
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void BoundaryAssembling2D::get_quadrature_formula_data(std::vector<double> &P,
    std::vector<double> &W, const TQuadFormula& qf1)
{
  // get the type of required quadrature (include/FE/Enumerations.h)
  // initialize points and weights of quadrature
  unsigned int nQuadPoints = qf1.GetN_QuadPoints();
  P.resize(nQuadPoints);
  W.resize(nQuadPoints);

  for (unsigned int i = 0; i < nQuadPoints; i++)
  {
    P[i] = qf1.get_point(i).x;
    W[i] = qf1.get_weight(i);
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void BoundaryAssembling2D::get_original_values(const FiniteElement& FE, int joint_id,
    TBaseCell *cell,
    const std::vector<double>& quadPoints,
    int BaseVectDim,
    std::vector< std::vector<double> > &originalValues,
    std::vector< std::vector<double> > &originalValues_x,
    std::vector< std::vector<double> > &originalValues_y,
    const TQuadFormula & qf)
{
  BFRefElements RefElement = FE.GetBaseFunct()->GetRefElement();
  int N_BaseFunct = FE.GetN_DOF();

  originalValues.resize(quadPoints.size());
  originalValues_x.resize(quadPoints.size());
  originalValues_y.resize(quadPoints.size());
  
  bool IsIsoparametric = false;
  if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
  {
    auto n_edges = cell->GetN_Joints();
    for(int i=0;i<n_edges;i++)
    {
      auto joint = cell->GetJoint(i);
      JointType jointtype = joint->GetType();
      if(jointtype == BoundaryEdge)
      {
        BoundTypes bdtype = ((const TBoundEdge *)(joint))->GetBoundComp()->GetType();
        if(bdtype != Line)
          IsIsoparametric = true;
      }
      if(jointtype == InterfaceJoint)
      {
        BoundTypes bdtype = ((const TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
        if(bdtype != Line)
          IsIsoparametric = true;
      }
      if(jointtype == IsoInterfaceJoint ||
          jointtype == IsoBoundEdge)
        IsIsoparametric = true;
    }
  }// endif
  auto ref_trans_id = FE.GetRefTransID();
  if(IsIsoparametric)
  {
    switch(RefElement)
    {
      case BFRefElements::BFUnitSquare:
        ref_trans_id = ReferenceTransformation_type::QuadIsoparametric;
        break;

      case BFRefElements::BFUnitTriangle:
        ref_trans_id = ReferenceTransformation_type::TriaIsoparametric;
        break;
      default:
        ErrThrow("unexpected Reference element ", RefElement);
        break;
    }
  } // endif IsIsoparametric
  auto ref_trans2D = FEDatabase::GetRefTrans2D(ref_trans_id);
  auto degree = FE.GetBaseFunct()->GetPolynomialDegree();
  auto qf_id = qf.get_type();
  switch(ref_trans_id)
  {
    case ReferenceTransformation_type::TriaAffin:
    case ReferenceTransformation_type::QuadAffin:
    case ReferenceTransformation_type::QuadBilinear:
      break;
    case ReferenceTransformation_type::TriaIsoparametric:
      ((TTriaIsoparametric *)ref_trans2D)->SetApproximationOrder(degree);
      ((TTriaIsoparametric *)ref_trans2D)->SetQuadFormula(qf_id);
      break;
    case ReferenceTransformation_type::QuadIsoparametric:
      ((TQuadIsoparametric *)ref_trans2D)->SetApproximationOrder(degree);
      ((TQuadIsoparametric *)ref_trans2D)->SetQuadFormula(qf_id);
      break;
    default:
      ErrThrow("unexpected Reference transformation ",
               static_cast<int>(ref_trans_id));
      break;
  } // endswitch
  ref_trans2D->SetCell(cell);
  

  for( unsigned int k = 0; k < quadPoints.size(); k++ )
  {
    // ----------------------------
    // for each quadrature point, get values and derivatives of basis functions
    // according to selected quadrature formula
    ///@todo write something like:
    // GetQuadratureValues(..., uorig, uxorig, uyorig)?
    // needs:
    // FE.GetBaseFunct2D(),qf,joint_id, cell
    // RefElement, quadPoints[k], N_BaseFunct
    double **reference_values = FEDatabase::GetJointDerivatives2D(
      *FE.GetBaseFunct(), qf, joint_id, MultiIndex2D::D00);
    double **derivative_xi = FEDatabase::GetJointDerivatives2D(
      *FE.GetBaseFunct(), qf, joint_id, MultiIndex2D::D10);
    double **derivative_eta = FEDatabase::GetJointDerivatives2D(
      *FE.GetBaseFunct(), qf, joint_id, MultiIndex2D::D01);

    std::vector<double> uorig(N_BaseFunct), uxorig(N_BaseFunct), uyorig(N_BaseFunct);
    ref_trans2D->GetOrigValues(joint_id, quadPoints[k],N_BaseFunct,
                               reference_values[k], derivative_xi[k],
                               derivative_eta[k], uorig.data(), uxorig.data(),
                               uyorig.data(), BaseVectDim);

    originalValues[k].resize(N_BaseFunct);
    originalValues_x[k].resize(N_BaseFunct);
    originalValues_y[k].resize(N_BaseFunct);
    for (int ib=0; ib<N_BaseFunct; ib++)
    {
      originalValues[k][ib] = uorig[ib];
      originalValues_x[k][ib] = uxorig[ib];
      originalValues_y[k][ib] = uyorig[ib];
    }
  }
}


// =================================================================================
void BoundaryAssembling2D::nitsche_bc(BlockFEMatrix &s_matrix,BlockVector &s_rhs,
              const TFESpace2D * v_space, const TFESpace2D *p_space,
              BoundValueFunct2D * U1, BoundValueFunct2D *U2,
              int bd_comp, double gamma,
              double mu, double sigma, double L_0,
              int sym_u, int sym_p)
{
  if (TDatabase::ParamDB->NSTYPE != 14)
  {
    Output::print<1>("WARNING BoundaryAssembling2D::nitsche_bc(..): "
            "The NSTYPE is not equal to 14. This might result in errors related to bad access.");
  }

  //============================== PENALTY TERMS ===================================
  // mueff * gamma/h * (u,v)
  matrix_u_v(s_matrix, v_space, bd_comp, gamma * mu, true);  // rescale local integral by edge values

  // mueff * gamma/h * (uD,v) [rhs]
  rhs_uD_v(s_rhs, v_space, U1, U2, bd_comp, gamma * mu, true);   // rescale local integral by edge values

  // sigma * L_0^2 * gamma/h (u.n,v.n)
  matrix_u_n_v_n(s_matrix, v_space,  bd_comp, gamma * sigma * L_0 * L_0, true); // rescale local integral by edge values

  // sigma * L_0^2 * gamma/h (uD.n,v.n)
  rhs_uD_n_v_n(s_rhs, v_space, U1, U2, bd_comp, gamma * sigma * L_0 * L_0, true);   // rescale local integral by edge values

  //=========================== PENALTY-FREE TERMS =================================
  // - (mu grad(u)n,v)
  matrix_gradu_n_v(s_matrix, v_space, bd_comp, -1. * mu); 

  // - sign_u * (u,mu grad(v)n) [sign_u=1: symmetrix, -1: skew-symmetric]
  matrix_gradv_n_u(s_matrix, v_space, bd_comp, (-1) * sym_u * mu);

  // - sign_u * (uD,mu grad(v)n) [rhs]
  rhs_gradv_n_uD(s_rhs, v_space, U1, U2, bd_comp, (-1) * sym_u * mu );

  // (pn,v)
  matrix_p_v_n(s_matrix, v_space, p_space, bd_comp, 1.);

  // sign_p * (u,qn)
  matrix_q_u_n(s_matrix, v_space, p_space, bd_comp, sym_p);

  // sign_p * (uD,qn) [rhs]
  rhs_q_uD_n(s_rhs, v_space, p_space, U1, U2, bd_comp, sym_p);
}

// =================================================================================
void BoundaryAssembling2D::nitsche_bc_nonlinear_iteration(BlockFEMatrix &s_matrix,
							  const TFESpace2D * v_space,
							  int bd_comp, double gamma,
							  double mu, double sigma, double L_0,
							  int sym_u)
{
  if (TDatabase::ParamDB->NSTYPE != 14)
  {
    Output::print<1>("WARNING BoundaryAssembling2D::nitsche_bc(..): "
            "The NSTYPE is not equal to 14. This might result in errors related to bad access.");
  }

  //============================== PENALTY TERMS ===================================
  // mueff * gamma/h * (u,v)
  matrix_u_v(s_matrix, v_space, bd_comp, gamma * mu, true);  // rescale local integral by edge values


  // sigma * L_0^2 * gamma/h (u.n,v.n)
  matrix_u_n_v_n(s_matrix, v_space,  bd_comp, gamma * sigma * L_0 * L_0, true); // rescale local integral by edge values

  //=========================== PENALTY-FREE TERMS =================================
  // - (mu grad(u)n,v)
  matrix_gradu_n_v(s_matrix, v_space, bd_comp, -1. * mu); 

  // - sign_u * (u,mu grad(v)n) [sign_u=1: symmetrix, -1: skew-symmetric]
  matrix_gradv_n_u(s_matrix, v_space, bd_comp, (-1) * sym_u * mu);


}
