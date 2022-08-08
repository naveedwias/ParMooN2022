#include <Database.h>
#include <Joint.h>
#include <BoundEdge.h>
#include "FEDatabase.h"
#include <FEFunction2D.h>
//#include <bits/c++config.h>
#include <string.h>
#include "QuadIsoparametric.h"
#include "TriaIsoparametric.h"
#include <NodalFunctional.h>
#include <MainUtilities.h>
#include <MooNMD_Io.h>
#include <AuxParam2D.h>
#include <GridCell.h>
#include <ConvDiff.h>
#include <Assemble_DG.h>

#include <InterfaceJoint.h>
#include <IsoBoundEdge.h>
#include <BdLine.h>
#include <LinAlg.h>
#include "ErrorEstimator.h"
#include "BlockVector.h"

#include <fstream>
#include <sstream>
#include <dirent.h> 
#include <type_traits>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>
#include <cstring> // memset

#include <BaseCell.h>
#include <BoundaryAssembling2D.h>
#include "QuadratureFormulaDatabase.h"

void OnlyDirichlet(int, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

TFEFunction2D::TFEFunction2D() :
  Name("dummy_fe_fct_2d")
{
  FESpace2D=nullptr;
  Values=nullptr;
}

/** constructor with vector initialization */
TFEFunction2D::TFEFunction2D(std::shared_ptr<const TFESpace2D> fespace2D,
                             const std::string& name, double *values)
: Name(name)
{
  Output::print<5>("Constructor of TFEFunction2D");
  FESpace2D=fespace2D;
  Values=values;
}

TFEFunction2D::TFEFunction2D(const ErrorEstimator<2>& estimator,
                             BlockVector& values)
 : Name("estimator")
{
  Output::print<5>("Constructor of TFEFunction2D using an ErrorEstimator");
  auto collection = estimator.get_collection();
  int n_cells = collection->GetN_Cells();
  FESpace2D = std::shared_ptr<const TFESpace2D>(
    new TFESpace2D(collection, "estimator_space",
                   BoundConditionNoBoundCondition, 0));
  if(FESpace2D->get_n_dof() != n_cells)
  {
    ErrThrow("picewise constant 2d FE space not properly constructed");
  }
  values.copy_structure(BlockVector(n_cells));
  Values = values.get_entries();
  
  
  for(int i = 0; i < n_cells; ++i)
  {
    auto DOF = FESpace2D->GetGlobalDOF(i);
    Values[DOF[0]] = estimator.get_eta_K(i);
  }
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/** Calculate L2-errors to a given function at the boundary */
/** written by Laura Blank 27.02.18 */
void TFEFunction2D::GetL2BoundaryError(BoundValueFunct2D *Exact,
                                       TAuxParam2D *, int,
                                       const TFESpace2D **fespaces,
                                       double *final_boundary_error_l2,
                                       bool rescale_by_h_E)
{
  // ########################################################################
  // loop over all Nitsche edges
  // ########################################################################
  // all spaces use the same Coll
  const TCollection *Coll = fespaces[0]->GetCollection();

  final_boundary_error_l2[0] = 0;
  for(int k = 0; k < TDatabase::ParamDB->n_nitsche_boundary; k++)
  {
    int boundary_component_id = TDatabase::ParamDB->nitsche_boundary_id[k];
    // Create a list of those boundary edges that are on the boundary component with given ID
    std::vector<TBoundEdge*> boundaryEdgeList;
    Coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);

    double boundary_error_l2_on_Component = 0;

    for(size_t m = 0; m < boundaryEdgeList.size(); m++)
    {
      TBoundEdge *boundedge = boundaryEdgeList[m];
      TBaseCell *cell = boundedge->GetNeighbour(0);
      // get basis dimension and FE space data of the current cell
      // Here it is assumed that only one local FE space is used for the hole mesh!
      auto element = FESpace2D->get_fe(cell->GetCellIndex());

      // calculate all needed values and derivatives of the basis functions
      if(element.GetBaseFunct()->GetBaseVectDim() > 1)
      {
        Output::print("The method TFEFunction2D::GetL2BoundaryErrors"
                      "is not designed for vector valued functions");
        Output::print("No errors were computed");
        return;
      }

      int BaseVectDim = 1; 
      int joint_id = boundedge->get_index_in_neighbour(cell);
      // get a quadrature formula good enough for the argument of the integral
      int fe_degree = element.GetBaseFunct()->GetPolynomialDegree();
      auto LineQuadFormula = QuadratureFormulaDatabase::qf_from_degree(
          std::max(TDatabase::ParamDB->INPUT_QUAD_RULE, 2*fe_degree),
          BFRefElements::BFUnitLine);
      std::vector<double> quadWeights, quadPoints;
      BoundaryAssembling2D::get_quadrature_formula_data(quadPoints, quadWeights,
                                                        *LineQuadFormula);
      // compute values of all basis functions and their first partial derivatives at all quadrature points
      std::vector< std::vector<double> > uorig, u_dx_orig ,u_dy_orig;
      BoundaryAssembling2D::get_original_values(element, joint_id, cell,
                                                quadPoints, BaseVectDim, uorig,
                                                u_dx_orig, u_dy_orig,
                                                *LineQuadFormula);

      int N_ = element.GetN_DOF();

      double FEFunctValues[MaxN_BaseFunctions2D];
      auto DOF = FESpace2D->GetGlobalDOF(cell->GetCellIndex());
      for(int l = 0; l < N_; l++)
      {
        FEFunctValues[l] = Values[DOF[l]];
      }

      double summing_boundary_error_l2_on_edge = 0;
      double edge_length = boundedge->get_length();
      double reference_edge_length = 2; // [-1,1] is the reference edge here

      for(size_t j = 0; j < quadPoints.size(); j++)
      {
        double value = 0;
        for(int l = 0; l < N_; l++)
        { 
          value += FEFunctValues[l] * uorig[j][l]; // compute u_h|T = \sum \alpha_i \Phi_i
        }
        auto comp = boundedge->GetBoundComp();
        int comp_ID = comp->GetID();
        double t0, t1;
        boundedge->GetParameters(t0, t1);
        double t = t0 + 0.5 * (t1-t0) * (quadPoints[j]+1);
        double ExactVal, Diff;
        Exact(comp_ID, t, ExactVal);
        Diff = ExactVal - value;

        if(rescale_by_h_E)
        {
          summing_boundary_error_l2_on_edge += quadWeights[j] * edge_length/reference_edge_length * 1/edge_length * Diff * Diff; 
        }
        else
        {
          summing_boundary_error_l2_on_edge += quadWeights[j] * edge_length/reference_edge_length * Diff * Diff; 
        }
      }
      boundary_error_l2_on_Component += summing_boundary_error_l2_on_edge;
    }
    final_boundary_error_l2[0] +=  boundary_error_l2_on_Component;
  }
 }

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/** Calculate L2-errors to a given function at the boundary (without the global database)*/
void TFEFunction2D::GetL2BoundaryError(BoundValueFunct2D *Exact,
                                       TAuxParam2D *, int,
                                       const TFESpace2D **fespaces,
                                       double *final_boundary_error_l2,
                                       int boundary_component_id,
                                       bool rescale_by_h_E)
{
  // ########################################################################
  // loop over all Nitsche edges
  // ########################################################################
  // all spaces use the same Coll
  const TCollection *Coll = fespaces[0]->GetCollection();

  final_boundary_error_l2[0] = 0;

    // Create a list of those boundary edges that are on the boundary component with given ID
    std::vector<TBoundEdge*> boundaryEdgeList;
    Coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);

    double boundary_error_l2_on_Component = 0;

    for(size_t m = 0; m < boundaryEdgeList.size(); m++)
    {
      TBoundEdge *boundedge = boundaryEdgeList[m];
      TBaseCell *cell = boundedge->GetNeighbour(0);
      // get basis dimension and FE space data of the current cell
      // Here it is assumed that only one local FE space is used for the hole mesh!
      auto element = FESpace2D->get_fe(cell->GetCellIndex());

      if(element.GetBaseFunct()->GetBaseVectDim() > 1)
      {
        Output::print("The method TFEFunction2D::GetL2BoundaryErrors",
                      "is not designed for vector valued functions");
        Output::print("No errors were computed");
        return;
      }

      int BaseVectDim = 1;
      int joint_id = boundedge->get_index_in_neighbour(cell);
      // get a quadrature formula good enough for the argument of the integral
      int fe_degree = element.GetBaseFunct()->GetPolynomialDegree();
      auto LineQuadFormula = QuadratureFormulaDatabase::qf_from_degree(
          std::max(TDatabase::ParamDB->INPUT_QUAD_RULE, 2*fe_degree),
          BFRefElements::BFUnitLine);
      std::vector<double> quadWeights, quadPoints;
      BoundaryAssembling2D::get_quadrature_formula_data(quadPoints, quadWeights,
                                                        *LineQuadFormula);
      // compute values of all basis functions and their first partial derivatives at all quadrature points
      std::vector< std::vector<double> > uorig, u_dx_orig ,u_dy_orig;
      BoundaryAssembling2D::get_original_values(element, joint_id, cell,
                                                quadPoints, BaseVectDim, uorig,
                                                u_dx_orig, u_dy_orig,
                                                *LineQuadFormula);

      int N_ = element.GetN_DOF();

      double FEFunctValues[MaxN_BaseFunctions2D];
      const int *DOF = FESpace2D->GetGlobalDOF(cell->GetCellIndex());
      for(int l = 0; l < N_; l++)
      {
        FEFunctValues[l] = Values[DOF[l]];
      }

      double summing_boundary_error_l2_on_edge = 0;
      double edge_length = boundedge->get_length();
      double reference_edge_length = 2; // [-1,1] is the reference edge here

      for(size_t j = 0; j < quadPoints.size(); j++)
      {
        double value = 0;
        for(int l = 0; l < N_; l++)
        {
          value += FEFunctValues[l] * uorig[j][l]; // compute u_h|T = \sum \alpha_i \Phi_i
        }
        auto comp = boundedge->GetBoundComp();
        int comp_ID = comp->GetID();
        double t0, t1;
        boundedge->GetParameters(t0, t1);
        double t = t0 + 0.5 * (t1-t0) * (quadPoints[j]+1);
        double ExactVal, Diff;
        Exact(comp_ID, t, ExactVal);
        Diff = ExactVal - value;

        if(rescale_by_h_E)
        {
          summing_boundary_error_l2_on_edge += quadWeights[j] * edge_length/reference_edge_length * 1/edge_length * Diff * Diff;
        }
        else
        {
          summing_boundary_error_l2_on_edge += quadWeights[j] * edge_length/reference_edge_length * Diff * Diff;
        }
      }
      boundary_error_l2_on_Component += summing_boundary_error_l2_on_edge;
    }
    final_boundary_error_l2[0] +=  boundary_error_l2_on_Component;
 }

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/** calculate errors to given function */
void TFEFunction2D::GetErrors(DoubleFunct2D *Exact, int N_Derivatives,
    MultiIndex2D *NeededDerivatives,
    int N_Errors, ErrorMethod *ErrorMeth, 
    CoeffFct2D Coeff, 
    TAuxParam2D *Aux,
    int n_fespaces, const TFESpace2D **fespaces,
    double *errors, bool,
    std::function<bool(const TBaseCell*, int)>funct,
    const ParameterDatabase& param_db) const
{
  
  bool* SecondDer = new bool[n_fespaces];
//   bool SecondDer[n_fespaces];
  // initialize the array SecondDer as either all true or all false.
  {
    // find out if one of the indices in NeededDerivatives refers to a second
    // order derivative
    bool need_second_derivatives = false;
    for(int i = 0; i < N_Derivatives; ++i)
    {
      if(NeededDerivatives[i] == MultiIndex2D::D20
         || NeededDerivatives[i] == MultiIndex2D::D11
         || NeededDerivatives[i] == MultiIndex2D::D02)
        need_second_derivatives = true;
    }
    for(int i = 0; i < n_fespaces; i++)
      SecondDer[i] = need_second_derivatives;
  }

  double *aux;
  double *Param[MaxN_QuadPoints_2D];
  int N_Parameters = Aux->GetN_Parameters();
  aux = new double [MaxN_QuadPoints_2D*N_Parameters];
  for(int j = 0; j < MaxN_QuadPoints_2D; j++)
    Param[j] = aux + j*N_Parameters;

  double *Derivatives[MaxN_QuadPoints_2D];
  aux = new double [MaxN_QuadPoints_2D*N_Derivatives];
  for(int j = 0; j < MaxN_QuadPoints_2D; j++)
    Derivatives[j] = aux + j*N_Derivatives;

  double *ExactVal[MaxN_QuadPoints_2D];
  aux = new double [MaxN_QuadPoints_2D * 4];
  for(int j = 0; j < MaxN_QuadPoints_2D; j++)
    ExactVal[j] = aux + j*4;

  // 20 <= number of term
  double *AuxArray[MaxN_QuadPoints_2D];
  aux = new double [MaxN_QuadPoints_2D*20];
  for(int j = 0; j < MaxN_QuadPoints_2D; j++)
    AuxArray[j] = aux + j*20;

  double max_error = 0.;
  memset(errors, 0, N_Errors * sizeof(double));

  // Generate databases needed for DG error
  auto db = LocalAssembling<2>::default_local_assembling_database();
  db.merge(default_dg_database(), true);
  db.merge(param_db, false); // update this database with given values
  // ########################################################################
  // loop over all cells
  // ########################################################################
  // all spaces use same Coll
  const TCollection *Coll = fespaces[0]->GetCollection();
  int N_Cells = Coll->GetN_Cells();
  double x_coord_max_error, y_coord_max_error;

  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);
  double AbsDetjk[MaxN_QuadPoints_2D];
  for(int k = 0; k < MaxN_QuadPoints_2D; ++k) AbsDetjk[k] = 1.;

  if (Coll->includes_hanging_vertices())
  { // hanging node
    Output::warn("Computation of DG error does not work for adaptive grids yet."
        );
  }
  Coll->mark_all_cells(); // needed for DG error


  for(int i = 0; i < N_Cells; i++)
  {
    TBaseCell *cell = Coll->GetCell(i);
    if ( funct(cell,-4711) )
      continue;

    // get basis dimension and FE space data of cell i
    auto element = FESpace2D->get_fe(i);
    if(element.GetBaseFunct()->GetBaseVectDim() > 1)
    {
      Output::print("For vector valued basis functions, you should use ",
                    "TFEFunction2D::GetErrorsForVectorValuedFunction instead ",
                    "of TFEFunction2D::GetErrors");
      Output::print("No errors were computed");
      return;
    }

    ReferenceTransformation_type RefTrans = element.GetRefTransID();
    // get quadrature formula id, be carefull: TDatabase::ParamDB->INPUT_QUAD_RULE is initialized with 0. It might be the case that the resulting Quadrature rule is not accurate enough for the exact solution

    FEDatabase::SetCellForRefTrans(cell, RefTrans);

    // get quadrature coordinates on original cell (also AbsDetjk is filled)
    double hK = cell->Get_hK(TDatabase::ParamDB->CELL_MEASURE);

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    std::vector< const FiniteElement*> used_elements(n_fespaces, nullptr);
    for(int j = 0; j < n_fespaces; j++)
    {
      auto& element = fespaces[j]->get_fe(i);
      used_elements[j] = &element;
    }

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    // Compute transformation of basis functions (and their derivatives) to original cell and the evaluation for all quadrature points 
    int N_QuadPoints;
    const double *weights, *X, *Y;
    FEDatabase::GetOrig(used_elements, Coll, cell, SecondDer, qf_ref, qf_orig);
    qf_orig.GetFormulaData(N_QuadPoints, weights, X, Y);
    
    if(N_Parameters>0)
      Aux->GetParameters(qf_orig, i, Param);

    // calculate all needed derivatives of this FE function
    int N_ = element.GetN_DOF();

    double FEFunctValues[MaxN_BaseFunctions2D];
    const int *DOF = FESpace2D->GetGlobalDOF(i);
    for(int l = 0; l < N_; l++)
    {
      FEFunctValues[l] = Values[DOF[l]];
    }

    for(int k = 0; k < N_Derivatives; k++)
    { // create a pointer to the QuadPoint-values of the basis on original element and its derivatives
      double **OrigFEValues = FEDatabase::GetOrigElementValues(
        *element.GetBaseFunct(), NeededDerivatives[k]);
      for(int j = 0; j < N_QuadPoints; j++)
      {
        double *Orig = OrigFEValues[j];
        double value = 0;

        for(int l = 0; l < N_; l++)
        {
          value += FEFunctValues[l] * Orig[l];  // compute u_h|T = \sum \alpha_i \Phi_i
        }                                         // endfor l
        // Derivatives[j][0] contains values of the FEFunction at QuadPoints
        // Derivatives[j][1] contains values of x derivatives of the FEFunction at QuadPoints
        // Derivatives[j][2] contains values of y derivatives of the FEFunction at QuadPoints↲
        Derivatives[j][k] = value;
      }                                           // endfor j
    }                                             // endfor k

    for(int j = 0; j < N_QuadPoints; j++)
    {
      auto p = qf_orig.get_point(j);
      Exact(p.x, p.y, ExactVal[j]);

      // Compute L^\inf-errors of the function in all quadrature points
      // and find the coordinates of the maximum
      double loc_error = std::abs(*ExactVal[j] - Derivatives[j][0]);
      //Output::print("Pointwise_Error in: ", p.x, " ", p.y, " is ", loc_error);
      if(loc_error > max_error)
      {
        max_error = loc_error;
        x_coord_max_error = p.x;
        y_coord_max_error = p.y;
      }
    }

    if(Coeff)
      Coeff(N_QuadPoints, X, Y, Param, AuxArray);

    double* LocError = new double[N_Errors];
    ErrorMeth(N_QuadPoints, {{X, Y}}, AbsDetjk, weights, hK, Derivatives, 
        ExactVal, AuxArray, LocError);


    // ########################################################################
    // DG-error
    // ########################################################################
    if(db["space_discretization_type"].is("dg"))
    {
      // IDEA: The  volume part of the DG norm is computed during the error
      // method, hence only the edge contributions have to be considered.
      // Therefore, do a loop over the edges and compute the respective
      // quantities.

      // Prepare quadrature for integration along edges
      auto deg = element.GetBaseFunct()->GetPolynomialDegree();  // degree
      // As degree for the quadrature rule at least 15 is chosen since the
      // exact solution can be of arbitrary degree / not even a polynomial at
      // all. Why 15? Just a choice. What would you suggest?
      // On the other hand, at least the square of the discrete solution has
      // to be integrated exactly -> 2*deg. This could be too less, if either
      // the diffusion or the convection or the reaction terms are
      // non-constant polynomials, since b*\nabla(u_h) + c u_h u_h has to be
      // integrated
      auto quad_form_1D = QuadratureFormulaDatabase::qf_from_degree(
          std::max(15,2*deg), BFRefElements::BFUnitLine);
      auto n_quad_pts_1D = quad_form_1D->GetN_QuadPoints();  // # quadrature pts

      // Compute dg solution uh and nabla(uh) in quadrature points on CELL
      auto n_loc_dof = element.GetN_DOF();  // # local dof
      auto DOF2 = FESpace2D->GetGlobalDOF(i); // global numbers of local dofs
      std::vector<double> FEFunctValues2(n_loc_dof);  // values corresp. to dofs
      for(int l = 0; l < n_loc_dof; l++)
      {
        FEFunctValues2[l] = Values[DOF2[l]];
      }

      std::vector<std::vector<double>> val_u_dg(n_quad_pts_1D,
          std::vector<double> (N_Derivatives));
      // val_u_dg[j][0] contains values of the FEFunction at QuadPoint j
      // val_u_dg[j][1/2] contains values of x/y derivative of the FEFunction at
      // QuadPoint j

      std::vector<double> quad_points;
      quad_points.resize(n_quad_pts_1D);
      for (int quad_pt_i = 0; quad_pt_i < n_quad_pts_1D; ++quad_pt_i)
      {
        quad_points[quad_pt_i] = quad_form_1D->get_point(quad_pt_i).x;
      }
      auto n_edges = cell->GetN_Edges();  // # edges
      for (int edge_i = 0; edge_i < n_edges; ++edge_i)
      {
        std::vector<std::vector<double>> orig_values;
        std::vector<std::vector<double>> orig_values_x_deriv;
        std::vector<std::vector<double>> orig_values_y_deriv;
        // orig_values[j] contains the values of the basis functions on the
        // original cell evaluated at quadrature point j.
        // orig_values_x_deriv (_y_) are the derivatives wrt. x (y)
        // This is filled in the next line
        BoundaryAssembling2D::get_original_values(*used_elements[0], edge_i,
            cell, quad_points, element.GetBaseFunct()->GetBaseVectDim(),
            orig_values, orig_values_x_deriv, orig_values_y_deriv,
            *quad_form_1D);

        for(int derivative_i = 0; derivative_i < N_Derivatives; derivative_i++)
        {
          for(int j = 0; j < n_quad_pts_1D; j++)
          {
            std::vector<double> aux_values;
            if (derivative_i == 0)
            {
              aux_values = orig_values[j];
            }
            else if (derivative_i == 1)
            {
              aux_values = orig_values_x_deriv[j];
            }
            else if (derivative_i == 2)
            {
              aux_values = orig_values_y_deriv[j];
            }
            else
            {
                ErrThrow("More derivatives needed than implemented.");
            }

            double value = 0;
            // compute u_h|T = \sum coeff_i x val_basis_fct_i
            for(int l = 0; l < n_loc_dof; l++)
            {
              value += aux_values[l] * FEFunctValues2[l];
            }                                         // endfor l
            val_u_dg[j][derivative_i] = value;
          }                                           // endfor j
        }                                             // endfor derivative_i

        // Compute exact solution u and nabla(u) in quadrature points
        double** val_u_ex = new double*[n_quad_pts_1D];
        aux = new double [n_quad_pts_1D * 4];
        for(int j = 0; j < n_quad_pts_1D; j++)
        {
          val_u_ex[j] = aux + j*4;
        }

        auto RefTrans = FEDatabase::GetOrig(used_elements, Coll, cell, 
                                            SecondDer, qf_ref, qf_orig);

        double* xi = new double [n_quad_pts_1D];
        double* eta = new double  [n_quad_pts_1D];
        auto bf2Drefelements = element.GetBaseFunct()->GetRefElement();
        for(int j = 0; j < n_quad_pts_1D; j++)
        {
          auto zeta = quad_form_1D->get_point(j).x;
          switch (bf2Drefelements)
          { // transformation from 1D to 2D reference elements
            case BFRefElements::BFUnitSquare:
              if (edge_i == 0)
              {
                xi[j] = zeta;
                eta[j] = -1;
              }
              else if (edge_i == 1)
              {
                xi[j] = 1;
                eta[j] = zeta;
              }
              else if (edge_i == 2)
              {
                xi[j] = -zeta;
                eta[j] = 1;
              }
              else
              {
                xi[j] = -1;
                eta[j] = -zeta;
              }
              break;
            case BFRefElements::BFUnitTriangle:
              if (edge_i == 0)
              {
                xi[j] = (1+zeta)/2;
                eta[j] = 0;
              }
              else if (edge_i == 1)
              {
                xi[j] = (-zeta+1)/2;
                eta[j] = (zeta+1)/2;
              }
              else
              {
                xi[j] = 0;
                eta[j] = (-zeta +1)/2;
              }
              break;
            default:
              ErrThrow("Unknown reference shape on cell ", i, '.');
          }

        } // endfor j
        double* X = new double [n_quad_pts_1D]; // x coord on original cell
        double* Y = new double [n_quad_pts_1D]; // y coord on original cell
        FEDatabase::GetOrigFromRef(RefTrans, n_quad_pts_1D, xi, eta, X, Y);
        for(int j = 0; j < n_quad_pts_1D; j++)
        { // compute exact values
          Exact(X[j], Y[j], val_u_ex[j]);
        }                                           // endfor j

        // Allocate space for jump and average of eh:=uh-u along edges. This
        // will be filled acording to whether the edge is an inner or boundary
        // edge. The exact solution is considered to be continuous along the
        // edge, and hence its jump across interior interfaces is zero. Along
        // the boundary the jump and the average of eh are defined as eh.
        // What to do with non-continuous exact solutions?
        std::vector<std::vector<double>> val_jump_eh(n_quad_pts_1D,
          std::vector<double> (N_Derivatives));
        std::vector<std::vector<double>> val_average_eh(n_quad_pts_1D,
          std::vector<double> (N_Derivatives));

        // we have to distinguish between Dirichlet and inner edges on the one
        // hand, and non-dirichlet edges on the other hand. The latter does not
        // play a role in the dg bilinearform
        bool is_Dirichlet_or_inner = true;

        // compute geometrical information about edge, needed for integration
        // and setting face_sigma_adaption
        double x0, y0, x1, y1, z0, z1;
        cell->GetVertex(edge_i)->GetCoords(x0, y0, z0);
        cell->GetVertex((edge_i+1) % n_edges)->GetCoords(x1, y1, z1);
        auto hE = std::sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
        auto nx = (y1-y0) / hE;
        auto ny = (x0-x1) / hE;

        // adaption parameter according to bilinearform. This will be modified
        // depending on the type of edge (inner, Dirichlet)
        double face_sigma_adaption = 1;

        // Check for neighbours to distinguish boundary and inner edges.
        auto joint = cell->GetJoint(edge_i);
        auto neigh = joint->GetNeighbour(cell);
        if (neigh)
        { // inner edge
          auto cell_nr_neigh = neigh->GetClipBoard(); // cell number of neighbor
          if (cell_nr_neigh < 0)
          { // hanging node
            continue;
          }

          // Every interior edge has to be considered only once. We consider it
          // only at their first apperance as local edge, i.e. if i is smaller
          // than cell_nr_neigh. In this case also the normal is defined as the
          // outer unit normal of cell i that was calculated before.
          if (i < cell_nr_neigh)
          { // interior edges was not considered yet

            auto cell_neigh = Coll->GetCell(cell_nr_neigh);
            int edge_nr_neigh=0;
            while(neigh->GetJoint(edge_nr_neigh)->GetNeighbour(cell_neigh)!=cell)
            {
              edge_nr_neigh ++;
            }

            // Compute dg solution uh and nabla(uh) in quadrature points on
            // NEIGHBOR.
            // ATTENTION: here it is assumed that on each cell there are the
            // same number of dofs or at least the quadrature rule is good
            // enough
            auto element_neigh = FESpace2D->get_fe(cell_nr_neigh);
            auto n_loc_dof_neigh = element_neigh.GetN_DOF();  // # local dof
            auto DOF3 = FESpace2D->GetGlobalDOF(cell_nr_neigh); // global
              // numbers of local dofs
            std::vector<double> FEFunctValues3(n_loc_dof_neigh);  // values
              // corresp. to dofs
            for(int dof_l = 0; dof_l < n_loc_dof_neigh; dof_l++)
            {
              FEFunctValues3[dof_l] = Values[DOF3[dof_l]];
            }

            std::vector<std::vector<double>> val_u_dg_neigh(n_quad_pts_1D,
                std::vector<double> (N_Derivatives));
            std::vector<std::vector<double>> orig_values_neigh;
            std::vector<std::vector<double>> orig_values_neigh_x_deriv;
            std::vector<std::vector<double>> orig_values_neigh_y_deriv;
            BoundaryAssembling2D::get_original_values(
                fespaces[0]->get_fe(cell_nr_neigh), edge_nr_neigh,
                cell_neigh, quad_points,
                element.GetBaseFunct()->GetBaseVectDim(), orig_values_neigh,
                orig_values_neigh_x_deriv, orig_values_neigh_y_deriv,
                *quad_form_1D);

            for(int derivative_i2 = 0; derivative_i2 < N_Derivatives;
                derivative_i2++)
            {
              for(int quad_pt_j = 0; quad_pt_j < n_quad_pts_1D; quad_pt_j++)
              {
                std::vector<double> values_bf;  // values of basis functions
                if (derivative_i2 == 0)
                {
                  values_bf = orig_values_neigh[quad_pt_j];
                }
                else if (derivative_i2 == 1)
                {
                  values_bf = orig_values_neigh_x_deriv[quad_pt_j];
                }
                else if (derivative_i2 == 2)
                {
                  values_bf = orig_values_neigh_y_deriv[quad_pt_j];
                }
                else
                {
                  ErrThrow("More derivatives needed than implemented.");
                }

                // compute u_h|T = \sum coeff_i x val_basis_fct_i
                double value2 = 0;
                for(int loc_dof_l = 0; loc_dof_l < n_loc_dof_neigh; loc_dof_l++)
                {
                  value2 += values_bf[loc_dof_l] * FEFunctValues3[loc_dof_l];
                } // endfor loc_dof_l
                val_u_dg_neigh[quad_pt_j][derivative_i2] = value2;
              } // endfor quad_pt_j
            } // endfor derivative_i2

            // Check wether the quadrature points seen from each cell are in the
            // same order and change order of val_u_dg_neigh if necessary
            std::vector< const FiniteElement* > used_elements_neigh(n_fespaces,
                nullptr);
            for(int fe_j = 0; fe_j < n_fespaces; fe_j++)
            {
              auto& element = fespaces[fe_j]->get_fe(cell_nr_neigh);
              used_elements_neigh[fe_j] = &element;
            }
            TQuadFormula qf_ref_neigh(qf_ref);
            TQuadFormula qf_orig_neigh(qf_ref);
            auto RefTrans_neigh = FEDatabase::GetOrig(used_elements_neigh,
                                                      Coll, cell_neigh,
                                                      SecondDer, qf_ref_neigh,
                                                      qf_orig_neigh);
            double* xi_neigh = new double;
            double* eta_neigh = new double;
            auto bf2Drefelements = element_neigh.GetBaseFunct()->
              GetRefElement();
            auto zeta_neigh = quad_form_1D->get_point(0).x;
            switch (bf2Drefelements)
            { // transformation from 1D to 2D reference elements
              case BFRefElements::BFUnitSquare:
                if (edge_nr_neigh == 0)
                {
                  xi_neigh[0] = zeta_neigh;
                  eta_neigh[0] = -1;
                }
                else if (edge_nr_neigh == 1)
                {
                  xi_neigh[0] = 1;
                  eta_neigh[0] = zeta_neigh;
                }
                else if (edge_nr_neigh == 2)
                {
                  xi_neigh[0] = -zeta_neigh;
                  eta_neigh[0] = 1;
                }
                else
                {
                  xi_neigh[0] = -1;
                  eta_neigh[0] = -zeta_neigh;
                }
                break;
              case BFRefElements::BFUnitTriangle:
                if (edge_nr_neigh == 0)
                {
                  xi_neigh[0] = (1+zeta_neigh)/2;
                  eta_neigh[0] = 0;
                }
                else if (edge_nr_neigh == 1)
                {
                  xi_neigh[0] = (-zeta_neigh+1)/2;
                  eta_neigh[0] = (zeta_neigh+1)/2;
                }
                else
                {
                  xi_neigh[0] = 0;
                  eta_neigh[0] = (-zeta_neigh +1)/2;
                }
                break;
              default:
                ErrThrow("Unknown reference shape on cell ", cell_nr_neigh, '.'
                    );
            }
            double* X_neigh = new double; // x coord on original cell
            double* Y_neigh = new double; // y coord on original cell
            FEDatabase::GetOrigFromRef(RefTrans_neigh, 1, xi_neigh, eta_neigh,
                                       X_neigh, Y_neigh);

            bool is_correct_order = false;
            if (X_neigh[0] == X[0] && Y_neigh[0] == Y[0])
            {
              is_correct_order = true;
            }
            if (is_correct_order == false)
            { // swap first-last, second-second last, ...
              for (int quad_pt_i = 0; quad_pt_i < n_quad_pts_1D / 2;
                  ++quad_pt_i)
              { // only perform half of the swaps since else you would swap back
                for (int der_j = 0; der_j < N_Derivatives; ++der_j)
                {
                  auto swap_val = val_u_dg_neigh[n_quad_pts_1D-1-quad_pt_i]
                    [der_j];
                  val_u_dg_neigh[n_quad_pts_1D-1-quad_pt_i][der_j] =
                    val_u_dg_neigh[quad_pt_i][der_j];
                  val_u_dg_neigh[quad_pt_i][der_j] = swap_val;
                } // endfor der_j
              } //endfor quad_pt_i
            } //endif is_correct_order

            // Compute jumps and averages of eh:=uh-u, nabla(eh)
            for (int quad_pt_i = 0; quad_pt_i < n_quad_pts_1D; ++quad_pt_i)
            {
              for (int der_j = 0; der_j < N_Derivatives; ++der_j)
              {
                val_jump_eh[quad_pt_i][der_j] = val_u_dg[quad_pt_i][der_j]-
                  val_u_dg_neigh[quad_pt_i][der_j];
                val_average_eh[quad_pt_i][der_j] =
                  (val_u_dg[quad_pt_i][der_j] +
                   val_u_dg_neigh[quad_pt_i][der_j]) / 2 - val_u_ex[quad_pt_i][der_j];
              } // endfor der_j
            } //endfor quad_pt_i
            delete xi_neigh;
            delete eta_neigh;
            delete X_neigh;
            delete Y_neigh;

          } // endif i < cell_nr_neigh
          else
          { // edge was already considered
            for (int quad_pt_i = 0; quad_pt_i < n_quad_pts_1D; ++quad_pt_i)
            {
              for (int der_j = 0; der_j < N_Derivatives; ++der_j)
              {
                val_jump_eh[quad_pt_i][der_j] = 0;
                val_average_eh[quad_pt_i][der_j] = 0;
              } // endfor der_j
            } //endfor quad_pt_i
          }
        } //endif inner edge
        else
        { // boundary edge
          if(joint->GetType() == BoundaryEdge ||
              joint->GetType() == IsoBoundEdge ||
              joint->GetType() == InterfaceJoint)
          {
            double t0, t1;
            const TBoundComp* BoundComp;
            if(joint->GetType() == BoundaryEdge||
                joint->GetType() == InterfaceJoint)
            {
              auto boundedge = (const TBoundEdge*) joint;
              BoundComp = boundedge->GetBoundComp();
              boundedge->GetParameters(t0, t1);
            }
            else
            {
              auto isoboundedge = (const TIsoBoundEdge*) joint;
              BoundComp = isoboundedge->GetBoundComp();
              isoboundedge->GetParameters(t0, t1);
            }
            // get id of the boundary component
            auto comp = BoundComp->GetID();
            // get type of the boundary condition at the beginning
            // and at the end of the current edge
            auto bound_cond_fct = FESpace2D->get_boundary_condition();
            BoundCond Cond0, Cond1;
            if (t0 < t1)
            {
              double eps = 1e-12;
              bound_cond_fct(comp, t0+eps, Cond0);
              bound_cond_fct(comp, t1-eps, Cond1);
            }
            else
            {
              double eps = 1e-12;
              bound_cond_fct(comp, t0-eps, Cond0);
              bound_cond_fct(comp, t1+eps, Cond1);
            }
            if(Cond0 == Cond1)
            {
              // Dirichlet
              if (Cond0== DIRICHLET)
              {
                // Compute jumps and averages of eh:=uh-u, nabla(eh)
                for (int quad_pt_i = 0; quad_pt_i < n_quad_pts_1D; ++quad_pt_i)
                {
                  for (int der_j = 0; der_j < N_Derivatives; ++der_j)
                  {
                    val_jump_eh[quad_pt_i][der_j] = val_u_dg[quad_pt_i][der_j]-
                      val_u_ex[quad_pt_i][der_j];
                    val_average_eh[quad_pt_i][der_j] =
                      val_u_dg[quad_pt_i][der_j]-val_u_ex[quad_pt_i][der_j];
                  } // endfor der_j
                } //endfor quad_pt_i
                face_sigma_adaption = 2;
              } // endif Dirichlet
              else
              { // Other edges than Dirichlet or inner ones are not considered
                // in bilinearform
                is_Dirichlet_or_inner = false;
              } // endelse non-Dirichlet
            } // endif Cond0 == Cond1
            else
            {
              ErrThrow("Two different boundary conditions on edge ", edge_i);
            }
          } // endif joint->GetType()
        } // endif boundary edge

        // compute the integrals and add them accordingly to bilinearform
        if (is_Dirichlet_or_inner)
        { // Only Dirichlet and innter edges are used in dg bilinearform
          auto face_sigma = db["face_sigma_DG"].get<double>();
          face_sigma *= face_sigma_adaption; // adapted stabilizing factor
          auto eta = db["eta_upwind_DG"].get<double>();

          double** parameters = new double*[n_quad_pts_1D];
          int n_parameters = Aux->GetN_Parameters();
          double* aux_dg = new double [n_quad_pts_1D * n_parameters];
          for(int j = 0; j < n_quad_pts_1D; j++)
          {
            parameters[j] = aux_dg + j * n_parameters;
          }

          // 20 <= number of term
          double** coefficients = new double*[n_quad_pts_1D];
          double* aux_dg2 = new double [n_quad_pts_1D * 20];
          for(int j = 0; j < n_quad_pts_1D; j++)
          {
            coefficients[j] = aux_dg2 + j*20;
          }

          if(N_Parameters>0)
          {
            Aux->GetParameters(*quad_form_1D, i, parameters);
          }
          if(Coeff)
          {
            Coeff(n_quad_pts_1D, X, Y, parameters, coefficients);
          }

          double integral = 0;
          
          if(Coeff)
          {
            for (int quad_pt_i = 0; quad_pt_i < n_quad_pts_1D; ++quad_pt_i)
            {
              double integrand = 0;
              // add terms to integrand according to DG bilinearform. The norm is
              // chosen to be ||e_h||^2 := A(e_h,e_h) where A is the NON-SYMMETRIC
              // bilinearform, i.e. only the volume terms, the penalty terms and
              // convective terms play a role.

              // penalty term to achieve coercitivity
              integrand += face_sigma/hE * val_jump_eh[quad_pt_i][0] *
                val_jump_eh[quad_pt_i][0];
              // contribution of convective term
              // This is exactly the upwind part of the bilinearform (see e.g.
              // - Kanschat: "Discontinuous Galerkin Methods for Viscous
              //   Incompressible Flow" (2008), chapter 1.4 or
              // - Houston, Schwab, Süli: SIAM J. Numer. Anal., 39(6), 2133–2163,
              //   doi = 10.1137/S0036142900374111)
              // even though it uses the point of view of stabilisation (see e.g.
              // - Brezzi, Marini, Süli: Mathematical Models and Methods in Applied
              //   Sciences, 14(12), 1893-1903, doi = 10.1142/S0218202504003866, or
              // - Di Pietro, Ern: "Mathematical Aspects of Discontinuous Galerkin
              //   Methods", doi = 10.1007/978-3-642-22980-0, chapter 2.3)
              auto inflow = coefficients[quad_pt_i][1] * nx +
                coefficients[quad_pt_i][2] * ny;
              if (inflow < 0 || face_sigma == db["face_sigma_DG"].get<double>())
              {
                // inflow boundary or interior boundary. Interi
                integrand -= (inflow) *
                  val_jump_eh[quad_pt_i][0] * val_average_eh[quad_pt_i][0];
              }
              if (neigh)
              {
                integrand += eta * 0.5 * std::abs(coefficients[quad_pt_i][1] * nx
                    + coefficients[quad_pt_i][2] * ny) * val_jump_eh[quad_pt_i][0]
                  * val_jump_eh[quad_pt_i][0];
              }
              auto weight = quad_form_1D->get_weight(quad_pt_i) *  hE / 2;
              integral += weight * integrand;
            }
          }
          
          if(N_Errors > 3)
            LocError[3] += integral;
          
          delete [] aux_dg;
          delete [] aux_dg2;
          
          delete [] coefficients;
          delete [] parameters;

        } //endif is_Dirichlet_or_inner
        delete [] aux;
        delete [] xi;
        delete [] eta;
        delete [] X;
        delete [] Y;
        delete [] val_u_ex;
      } //endfor edge_i
    } // end DG-stuff
    // ########################################################################

#ifdef __2D__
    if (!(ErrorMeth == &conv_diff_l2_h1_linf_error<2>))
    {
      for(int j = 0; j < N_Errors; j++)
        errors[j] += LocError[j];
    }
    else
    {
      for(int j = 0; j < N_Errors-1; j++)
        errors[j] += LocError[j];
      // L_infty error
      if(errors[N_Errors-1] <  LocError[N_Errors-1])
        errors[N_Errors-1] = LocError[N_Errors-1];
    }
#endif
    
    delete[] LocError;
  } // end for i, loop over cells

#ifdef __2D__
  if (!(ErrorMeth == &L1Error))
  {
    if (!(ErrorMeth == &conv_diff_l2_h1_linf_error<2>))
    {
      for(int j = 0; j < N_Errors; j++)
        errors[j] = std::sqrt(errors[j]);
    }
    else
    {
      for(int j = 0; j < N_Errors-1; j++)
        errors[j] = std::sqrt(errors[j]);
    }
  }
#endif

  if(max_error != 0.)
    Output::print<1>("Maximal pointwise error appeared in: ", x_coord_max_error,
                     " ", y_coord_max_error, " value: ", max_error);
  else
    Output::print<1>("Maximal pointwise error is zero");

  delete [] AuxArray[0];
  delete [] ExactVal[0];
  delete [] Derivatives[0];
  delete [] Param[0];

  delete[] SecondDer;
}                                                 // TFEFunction2D::GetErrors


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double TFEFunction2D::get_L2_norm_on_boundary(int boundary_component) const
{
  auto collection = FESpace2D->GetCollection();
  auto n_cells = collection->GetN_Cells();
  double l2_norm_on_boundary = 0.;
  for(int i_cell = 0; i_cell < n_cells; ++i_cell)
  {
    auto cell = collection->GetCell(i_cell);
    auto n_joints = cell->GetN_Joints();
    for(int i_joint = 0; i_joint < n_joints; ++i_joint)
    {
      auto joint = cell->GetJoint(i_joint);
      if(joint->GetType() == BoundaryEdge)
      {
        auto *boundary_joint = dynamic_cast<TBoundEdge*>(joint);
        auto component = boundary_joint->GetBoundComp();
        if(boundary_component < 0 || component->GetID() == boundary_component)
        {
          // found a joint which is on the desired boundary
          auto fe = FESpace2D->get_fe(i_cell);
          auto basis_functions = fe.GetBaseFunct();
          auto fe_degree = basis_functions->GetPolynomialDegree();
          auto n_local_basis_functions = basis_functions->GetDimension();
          auto LineQuadFormula = QuadratureFormulaDatabase::qf_from_degree(
              fe_degree, BFRefElements::BFUnitLine);
          std::vector<double> quadWeights, quadPoints;
          BoundaryAssembling2D::get_quadrature_formula_data(
            quadPoints, quadWeights, *LineQuadFormula);
          std::vector<std::vector<double>> uorig, u_dx_orig ,u_dy_orig;
          BoundaryAssembling2D::get_original_values(
            fe, i_joint, cell, quadPoints, basis_functions->GetBaseVectDim(), 
            uorig, u_dx_orig, u_dy_orig, *LineQuadFormula);
          
          auto local_to_global_dof = FESpace2D->GetGlobalDOF(i_cell);
          double* FEFunctValues = new double[n_local_basis_functions];
          for(int l = 0; l < n_local_basis_functions; l++)
          {
            FEFunctValues[l] = Values[local_to_global_dof[l]];
          }
          double edge_length = boundary_joint->get_length();
          double reference_edge_length = 2; // [-1,1] is the reference edge
          
          for(auto j = 0u; j < quadPoints.size(); j++)
          {
            double value = 0;
            for(int l = 0; l < n_local_basis_functions; l++)
            { 
              value += FEFunctValues[l] * uorig[j][l]; // compute u_h|T = \sum \alpha_i \Phi_i
            }
            l2_norm_on_boundary += quadWeights[j] * (value * value) 
                                   * edge_length/reference_edge_length;
          }
          
          delete[] FEFunctValues;
        }
      }
    }
  }
  return std::sqrt(l2_norm_on_boundary);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double TFEFunction2D::get_L2_norm() const
{
  auto collection = FESpace2D->GetCollection();
  auto n_cells = collection->GetN_Cells();
  double l2_norm = 0.;
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);
  for(int i_cell = 0; i_cell < n_cells; ++i_cell)
  {
    auto cell = collection->GetCell(i_cell);
    auto fe = FESpace2D->get_fe(i_cell);
    auto basis_functions = fe.GetBaseFunct();
    auto n_local_basis_functions = basis_functions->GetDimension();
    bool SecondDer = false;

    // Compute transformation of basis functions (and their derivatives) to 
    // original cell and the evaluation for all quadrature points 
    std::vector<const FiniteElement*> used_fe(1, &fe);
    FEDatabase::GetOrig(used_fe, collection, cell, &SecondDer, qf_ref, qf_orig);
    
    int N_QuadPoints = qf_orig.GetN_QuadPoints();

    // calculate all needed derivatives of this FE function
    double* FEFunctValues = new double[n_local_basis_functions];
    const int *DOF = FESpace2D->GetGlobalDOF(i_cell);
    for(int l = 0; l < n_local_basis_functions; l++)
    {
      FEFunctValues[l] = Values[DOF[l]];
    }

    // create a pointer to the QuadPoint-values of the basis on original element and its derivatives
    double **OrigFEValues = FEDatabase::GetOrigElementValues(
        *basis_functions, MultiIndex2D::D00);
    for(int j = 0; j < N_QuadPoints; j++)
    {
      double *Orig = OrigFEValues[j];
      double value = 0;

      for(int l = 0; l < n_local_basis_functions; l++)
      {
        value += FEFunctValues[l] * Orig[l];
      }
      l2_norm += value * value * qf_orig.get_weight(j);
    }
    
    delete[] FEFunctValues;
  }
  return std::sqrt(l2_norm);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std::pair<double,double> TFEFunction2D::get_L2_and_H1_norm() const
{
  auto collection = FESpace2D->GetCollection();
  auto n_cells = collection->GetN_Cells();
  double l2_norm = 0.;
  double h1_norm = 0.;
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);
  for(int i_cell = 0; i_cell < n_cells; ++i_cell)
  {
    auto cell = collection->GetCell(i_cell);
    auto fe = FESpace2D->get_fe(i_cell);
    auto basis_functions = fe.GetBaseFunct();
    auto n_local_basis_functions = basis_functions->GetDimension();
    bool SecondDer = false;

    // Compute transformation of basis functions (and their derivatives) to 
    // original cell and the evaluation for all quadrature points 
    std::vector<const FiniteElement*> used_fe(1, &fe);
    FEDatabase::GetOrig(used_fe, collection, cell, &SecondDer, qf_ref, qf_orig);
    
    int N_QuadPoints = qf_orig.GetN_QuadPoints();

    // calculate all needed derivatives of this FE function
    double* FEFunctValues = new double[n_local_basis_functions];
    const int *DOF = FESpace2D->GetGlobalDOF(i_cell);
    for(int l = 0; l < n_local_basis_functions; l++)
    {
      FEFunctValues[l] = Values[DOF[l]];
    }

    // create a pointer to the QuadPoint-values of the basis on original element and its derivatives
    double **OrigFEValues = FEDatabase::GetOrigElementValues(
        *basis_functions, MultiIndex2D::D00);
    double **OrigFED10 = FEDatabase::GetOrigElementValues(*basis_functions,
                                                          MultiIndex2D::D10);
    double **OrigFED01 = FEDatabase::GetOrigElementValues(*basis_functions,
                                                          MultiIndex2D::D01);
    for(int j = 0; j < N_QuadPoints; j++)
    {
      double *Orig = OrigFEValues[j];
      double *Orig10 = OrigFED10[j];
      double *Orig01 = OrigFED01[j];
      double value = 0;
      double d10 = 0;
      double d01 = 0;

      for(int l = 0; l < n_local_basis_functions; l++)
      {
        value += FEFunctValues[l] * Orig[l];
        d10 += FEFunctValues[l] * Orig10[l];
        d01 += FEFunctValues[l] * Orig01[l];
      }
      l2_norm += value * value * qf_orig.get_weight(j);
      h1_norm += (d10*d10 + d01*d01) * qf_orig.get_weight(j);
    }
    
    delete[] FEFunctValues;
  }
  return {std::sqrt(l2_norm), std::sqrt(h1_norm)};
}

//==========================================================================

void TFEFunction2D::getValuesFromOriginalOnes(double* retValues, int valuesOffset, int N_BaseFunct, int cellNumber, double* uorig, double* uxorig, double* uyorig) const
{
  retValues[0] = 0.0;
  retValues[1] = 0.0;
  retValues[2] = 0.0;
  
  const int* Numbers = FESpace2D->GetGlobalDOF(cellNumber);
  
  for(int j = 0; j < N_BaseFunct; j++)
  {
    double val = Values[valuesOffset + Numbers[j]];
    retValues[0] += uorig[j]*val;
    retValues[1] += uxorig[j]*val;
    retValues[2] += uyorig[j]*val;
  }
}

//==========================================================================

void TFEFunction2D::FindGradient(double x, double y, double *values) const
{
  int N_Found = 0;
  
  double valuesFromCells[3];
  
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;

  const TCollection* Coll = FESpace2D->GetCollection();
  int N_Cells = Coll->GetN_Cells();
  for(int i = 0; i < N_Cells; i++)
  {
    TBaseCell* cell = Coll->GetCell(i);
    if(cell->PointInCell(parmoon::Point(x,y)))
    {
      N_Found++;
      FindGradientLocal(cell, i, x, y, valuesFromCells);
      
      for(int j = 0; j < 3; j++)
        values[j] += valuesFromCells[j];
    }
  }

  if(N_Found > 0)
  {
    for(int i = 0; i < 3; i++)
      values[i] /= N_Found;
  }
  else
  {
    ErrThrow(" Point not found ", "x ", x, " y ", y);
  }
}


void TFEFunction2D::FindGradientLocal(const TBaseCell *cell, int cell_no, 
                                      double x, double y, double *values) const
{
  const FiniteElement& FE_Obj = FESpace2D->get_fe(cell_no);
  const BaseFunctions* bf = FE_Obj.GetBaseFunct();
  
  int N_BaseFunct = bf->GetDimension();
  int BaseVectDim = bf->GetBaseVectDim();
  
  double *uorig = new double[N_BaseFunct*BaseVectDim];
  double *uxorig = new double[N_BaseFunct*BaseVectDim];
  double *uyorig = new double[N_BaseFunct*BaseVectDim];
  
  
  if(BaseVectDim == 1)
  {
    getOrigValuesForCell(cell_no, x, y, uorig, uxorig, uyorig);
    
    getValuesFromOriginalOnes(values, 0, N_BaseFunct, cell_no, uorig, uxorig, uyorig);
  }
  else
  {
    // for vector valued basis functions (Raviart-Thomas elements), some signs
    // must be changed
    double xi, eta;
    
    auto Coll = FESpace2D->GetCollection();
    ReferenceTransformation_type RefTrans = FE_Obj.GetRefTransID();
    
    // set cell for reference transformation
    FEDatabase::SetCellForRefTrans(cell, RefTrans);
    
    // find local coordinates of the given point
    FEDatabase::GetRefFromOrig(RefTrans, x, y, xi, eta);

    double *uref = new double[N_BaseFunct*BaseVectDim];
    double *uxiref = new double[N_BaseFunct*BaseVectDim];
    double *uetaref = new double[N_BaseFunct*BaseVectDim];

    // get values and derivatives of basis functions on the
    // reference mesh cell
    bf->GetDerivatives(MultiIndex2D::D00, xi, eta, uref);
    bf->GetDerivatives(MultiIndex2D::D10, xi, eta, uxiref);
    bf->GetDerivatives(MultiIndex2D::D01, xi, eta, uetaref);

    // compute values on the original mesh cell
    FEDatabase::GetOrigValues(RefTrans, xi, eta, bf, Coll, (TGridCell *)cell,
                              uref, uxiref, uetaref, uorig, uxorig, uyorig);
    const double eps = 1e-20;
    
    auto Numbers = FESpace2D->GetGlobalDOF(cell_no);
    for(int i = 0; i < BaseVectDim; i++)
    {
      double u = 0;
      double ux = 0;
      double uy = 0;
      for(int j=0;j<N_BaseFunct;j++)
      {
        int k = Numbers[j];
        double val = Values[k];
        u += uorig[j]*val;
        ux += uxorig[j]*val;
        uy += uyorig[j]*val;
      }

      if (std::abs(ux)<eps)
        ux=0.0;
      if (std::abs(uy)<eps)
        uy=0.0;

      values[    3*i] = u;
      values[1 + 3*i] = ux;
      values[2 + 3*i] = uy;
    }
    
    delete [] uref;
    delete [] uxiref;
    delete [] uetaref;
  }
  

  delete [] uorig;
  delete [] uxorig;
  delete [] uyorig;
}


/** determine the value of function at the given point which lies within the
  cell *cell. This also works for vector valued basis functions as are used 
  for Raviart-Thomas elements.
 */
void TFEFunction2D::FindValueLocal(const TBaseCell *cell, int cell_no, double x,
                                   double y, double *values) const
{
  double xi, eta;
  auto Coll = FESpace2D->GetCollection();

  auto FE_Obj = FESpace2D->get_fe(cell_no);
  ReferenceTransformation_type RefTrans = FE_Obj.GetRefTransID();

  // set cell for reference transformation
  FEDatabase::SetCellForRefTrans(cell, RefTrans);

  // find local coordinates of the given point
  FEDatabase::GetRefFromOrig(RefTrans, x, y, xi, eta);
  // cout << " xi: " << xi << endl;
  // cout << "eta: " << eta << endl;

  // get base function object
  auto bf = FE_Obj.GetBaseFunct();
  int N_BaseFunct = bf->GetDimension();
  int BaseVectDim = bf->GetBaseVectDim();

  // allocate arrays
  double *uorig = new double[N_BaseFunct*BaseVectDim];
  double *uxorig = new double[N_BaseFunct*BaseVectDim];
  double *uyorig = new double[N_BaseFunct*BaseVectDim];

  double *uref = new double[N_BaseFunct*BaseVectDim];
  double *uxiref = new double[N_BaseFunct*BaseVectDim];
  double *uetaref = new double[N_BaseFunct*BaseVectDim];

  // get values and derivatives of basis functions on the
  // reference mesh cell
  bf->GetDerivatives(MultiIndex2D::D00, xi, eta, uref);
  bf->GetDerivatives(MultiIndex2D::D10, xi, eta, uxiref);
  bf->GetDerivatives(MultiIndex2D::D01, xi, eta, uetaref);

  // compute values on the original mesh cell
  FEDatabase::GetOrigValues(RefTrans, xi, eta, bf, Coll, (TGridCell *)cell,
                            uref, uxiref, uetaref, uorig, uxorig, uyorig);

  auto Numbers = FESpace2D->GetGlobalDOF(cell_no);
  for (int i=0; i<BaseVectDim; i++)
  {
    double u = 0;
    for(int j=0;j<N_BaseFunct;j++)
    {
      int k = Numbers[j];
      double val = Values[k];
      u += uorig[j+i*N_BaseFunct]*val;
    }

    values[i] = u;
  }

  // set memory free
  delete [] uorig;
  delete [] uxorig;
  delete [] uyorig;
  delete [] uref;
  delete [] uxiref;
  delete [] uetaref;
}


/** calculate the interpolation of an exact function */
void TFEFunction2D::Interpolate(DoubleFunct2D *Exact, DoubleFunct2D *Exact2)
{
  int N_Points;
  const double *xi, *eta;
  TRefTrans2D *rt;
  int BaseVectDim = FESpace2D->get_fe(0).GetBaseFunct()->GetBaseVectDim();
  if((BaseVectDim > 1 && !Exact2) || (BaseVectDim == 1 && Exact2))
  {
    ErrThrow("for vector valued basis functions you need to provide two "
             "functions to Interpolate; for scalar functions only one.");
  }
  double X[MaxN_PointsForNodal2D], Y[MaxN_PointsForNodal2D];
  double PointValues[MaxN_PointsForNodal2D*BaseVectDim];
  double FunctionalValues[MaxN_PointsForNodal2D*BaseVectDim];
  double FctVal[4];
  // begin code

  auto Coll = FESpace2D->GetCollection();
  int N_Cells = Coll->GetN_Cells();
  int N_DOFs = FESpace2D->get_n_dof();

  memset(Values, 0, sizeof(double)*N_DOFs);

  for(int i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);
    auto Element = FESpace2D->get_fe(i);
    auto nf = Element.GetNodalFunctional();
    nf->GetPointsForAll(N_Points, xi, eta);
    int N_LocalDOFs = Element.GetN_DOF();

    int PolynomialDegree = Element.GetBaseFunct()->GetPolynomialDegree();
    int ApproxOrder = Element.GetBaseFunct()->GetAccuracy();

    BFRefElements RefElement = Element.GetBaseFunct()->GetRefElement();

    auto QuadFormula = QuadratureFormulaDatabase::qf_from_degree(
        3*PolynomialDegree, RefElement);
    ReferenceTransformation_type RefTrans = Element.GetRefTransID();

    bool IsIsoparametric = TDatabase::ParamDB->USE_ISOPARAMETRIC
                           && cell->has_isoparametric_joint();
    if(IsIsoparametric)
    {
      switch(RefElement)
      {
        case BFRefElements::BFUnitSquare:
          RefTrans = ReferenceTransformation_type::QuadIsoparametric;
          break;

        case BFRefElements::BFUnitTriangle:
          RefTrans = ReferenceTransformation_type::TriaIsoparametric;
          break;
          
        default:
          ;
      }
    }                                             // endif IsIsoparametric

    rt = FEDatabase::GetRefTrans2D(RefTrans);
    switch(RefTrans)
    {
      case ReferenceTransformation_type::TriaAffin:
      case ReferenceTransformation_type::QuadAffin:
      case ReferenceTransformation_type::QuadBilinear:
        break;
      case ReferenceTransformation_type::QuadIsoparametric:
        ((TQuadIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((TQuadIsoparametric *)rt)->SetQuadFormula(QuadFormula->get_type());
        break;
      case ReferenceTransformation_type::TriaIsoparametric:
        ((TTriaIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((TTriaIsoparametric *)rt)->SetQuadFormula(QuadFormula->get_type());
        break;
      default:
        ErrThrow("unknown reftrans id: ", static_cast<int>(RefTrans));
        break;
    }
    rt->SetCell(cell);
    rt->GetOrigFromRef(N_Points, xi, eta, X, Y);

    for(int j=0;j<N_Points;j++)
    {
      Exact(X[j], Y[j], FctVal);
      PointValues[j] = FctVal[0];
      if(Exact2)
      {
        Exact2(X[j], Y[j], FctVal);
        PointValues[j+N_Points] = FctVal[0];
      }
    }

    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValues,
        FunctionalValues);

    auto DOF = FESpace2D->GetGlobalDOF(i);

    for(int j=0;j<N_LocalDOFs;j++)
      Values[DOF[j]] = FunctionalValues[j];
  }
}

/** calculate parameters which are connected to a mesh cell */
void TFEFunction2D::GetMeshCellParams(DoubleFunct2D* Exact, int N_Derivatives, 
    MultiIndex2D* NeededDerivatives,
    int N_Errors, ErrorMethod* ErrorMeth,
    CoeffFct2D Coeff, TAuxParam2D* Aux,
    int n_fespaces,
    const TFESpace2D** fespaces,
    double* errors, double* parameters)
{
  int i,j,k,l;
  int N_Cells, N_Points, N_Parameters, N_;
  const double *weights;
  const double *X, *Y;
  std::array<double, MaxN_QuadPoints_2D> AbsDetjk;
  std::fill(AbsDetjk.begin(), AbsDetjk.end(), 1.);
  double *Param[MaxN_QuadPoints_2D], *aux;
  double *Derivatives[MaxN_QuadPoints_2D];
  double *ExactVal[MaxN_QuadPoints_2D];
  double *AuxArray[MaxN_QuadPoints_2D];
  double **OrigFEValues, *Orig, value;
  double FEFunctValues[MaxN_BaseFunctions2D];
  double LocError[4];
  double hK;
  bool *SecondDer;

  SecondDer = new bool[n_fespaces];
  for(i=0;i<n_fespaces;i++)
    SecondDer[i] = false;

  N_Parameters = Aux->GetN_Parameters();
  aux = new double [MaxN_QuadPoints_2D*N_Parameters];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Param[j] = aux + j*N_Parameters;

  aux = new double [MaxN_QuadPoints_2D*N_Derivatives];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Derivatives[j] = aux + j*N_Derivatives;

  aux = new double [MaxN_QuadPoints_2D * 4];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    ExactVal[j] = aux + j*4;

  // 20 <= number of term
  aux = new double [MaxN_QuadPoints_2D*20];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    AuxArray[j] = aux + j*20;

  for(i=0;i<N_Errors;i++)
    errors[i] = 0.0;

  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);
  // ########################################################################
  // loop over all cells
  // ########################################################################
  auto Coll = fespaces[0]->GetCollection();         // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);

    hK = cell->GetDiameter();

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    std::vector<const FiniteElement*> used_elements(n_fespaces, nullptr);
    for(j=0;j<n_fespaces;j++)
    {
      auto& element = fespaces[j]->get_fe(i);
      used_elements[j] = &element;
    }

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    FEDatabase::GetOrig(used_elements, Coll, cell, SecondDer, qf_ref, qf_orig);
    qf_orig.GetFormulaData(N_Points, weights, X, Y);

    if(N_Parameters>0)
      Aux->GetParameters(qf_orig, i, Param);

    // calculate all needed derivatives of this FE function
    auto element = FESpace2D->get_fe(i);
    N_ = element.GetN_DOF();

    auto DOF = FESpace2D->GetGlobalDOF(i);
    for(l=0;l<N_;l++)
      FEFunctValues[l] = Values[DOF[l]];

    for(k=0;k<N_Derivatives;k++)
    {
      OrigFEValues = FEDatabase::GetOrigElementValues(*element.GetBaseFunct(),
                                                      NeededDerivatives[k]);
      for(j=0;j<N_Points;j++)
      {
        Orig = OrigFEValues[j];
        value = 0;
        for(l=0;l<N_;l++)
        {
          value += FEFunctValues[l] * Orig[l];
        }                                         // endfor l
        Derivatives[j][k] = value;
      }                                           // endfor j
    }                                             // endfor k

    for(j=0;j<N_Points;j++)
      Exact(X[j], Y[j], ExactVal[j]);

    if(Coeff)
      Coeff(N_Points, X, Y, Param, AuxArray);

    ErrorMeth(N_Points, {{X, Y}}, AbsDetjk.data(), weights, hK, Derivatives,
        ExactVal, AuxArray, LocError);

    for(j=0;j<N_Errors;j++)
    {
      errors[j] += LocError[j];
      parameters[i + j *N_Cells] = LocError[j];
    }

  }                                               // endfor i

  for(j=0;j<N_Errors;j++)
    errors[j] = std::sqrt(errors[j]);

  delete Param[0];
  delete AuxArray[0];
  delete SecondDer;
  delete ExactVal[0];
  delete Derivatives[0];
}                                                 // TFEFunction2D::GetErrors


/*
   Compute errors for one vector valued function. This is needed if you use for 
   example Raviart-Thomas elements and can't compute the error separetely for each
   component. This happens if you want to compute the L2-error of the divergence.

   Output is double *errors with length 3:
   errors[0] -- L2-Error for this function
   errors[1] -- L2-Error of divergence of this function 
   errors[2] -- H1-semi-Norm-Error of this function
 */
void TFEFunction2D::GetErrorsForVectorValuedFunction(
    DoubleFunct2D * const * const Exact, 
    ErrorMethod * const ErrorMeth, 
    double * const errors) const
{
  // set all errors to zero at first
  memset(errors,0,3*sizeof(double));
  auto coll = FESpace2D->GetCollection(); 
  const int n_cells = coll->GetN_Cells();
  TQuadFormula qf_orig(QuadratureFormula_type::BaryCenterTria); // dummy type

  for(int i = 0; i < n_cells; i++)
  {
    TBaseCell *cell = coll->GetCell(i);
    auto element = FESpace2D->get_fe(i);
    // number of basis functions
    int n_base_functs = element.GetN_DOF();
    const BaseFunctions* bf = element.GetBaseFunct();
    const int * DOF = FESpace2D->GetGlobalDOF(i);
    ReferenceTransformation_type ref_trans = element.GetRefTransID();
    FEDatabase::SetCellForRefTrans(cell, ref_trans);
    BFRefElements ref_element = element.GetBaseFunct()->GetRefElement();
    // get quadrature formula id
    const TQuadFormula *qf = QuadratureFormulaDatabase::qf_from_degree(
          TDatabase::ParamDB->INPUT_QUAD_RULE, ref_element);

    int n_points = qf->GetN_QuadPoints(); // number of quadrature points in one cell
    // get quadrature coordinates on original cell
    FEDatabase::GetOrigFromRef(ref_trans, *qf, qf_orig);
    // will store values of exact solution at all quadrature points in one cell
    std::vector<double*> ExactVal(n_points, nullptr);
    for(int j = 0; j < n_points; j++)
    {
      auto p = qf_orig.get_point(j);
      ExactVal[j] = new double[8]; // 4 values for each of the two components 
      Exact[0](p.x, p.y, ExactVal[j]); // x-component
      Exact[1](p.x, p.y, ExactVal[j]+4); // y-component
    }
    // will store values of this FEFunction and its first derivatives at all 
    // quadrature points
    // the six means values and 2 first derivatives for both components
    double** Derivatives = new double*[n_points];
    for(int j = 0; j < n_points; j++)
    {
      Derivatives[j] = new double[6];  
    }
    // Get the function values of this FE-function at the local dofs.
    double* FEFunctValues = new double[n_base_functs];
    for(int l = 0; l < n_base_functs; l++)
    {
      FEFunctValues[l] = Values[DOF[l]];
    }

    double* uref = new double[n_base_functs*2]; // 2 because vector valued basis functions
    double* uxiref = new double[n_base_functs*2]; // 2 because vector valued basis functions
    double* uetaref = new double[n_base_functs*2]; // 2 because vector valued basis functions

    double*** AllOrigValues = new double**[3];
    
    for(int iOrig = 0; iOrig < 3; iOrig++)
    {
      AllOrigValues[iOrig] = new double*[n_points];
      
      for(int jOrig = 0; jOrig < n_points; jOrig++)
        AllOrigValues[iOrig][jOrig] = new double[2*n_base_functs];
    }
    
    
    std::vector<double> weights(n_points); // quadrature weights
    std::vector<double> AbsDetjk(n_points, 1.);
    std::vector<double> XYZ(2*n_points);
    for(int k = 0; k < n_points; k++)
    {
      auto p = qf->get_point(k);
      bf->GetDerivatives(MultiIndex2D::D00, p.x, p.y, uref);
      bf->GetDerivatives(MultiIndex2D::D10, p.x, p.y, uxiref);
      bf->GetDerivatives(MultiIndex2D::D01, p.x, p.y, uetaref);
      // Piola transform
      FEDatabase::GetOrigValues(ref_trans, p.x, p.y, bf, coll,
                                (TGridCell *) cell, uref, uxiref, uetaref,
                                AllOrigValues[0][k], AllOrigValues[1][k],
                                AllOrigValues[2][k]);
      p = qf_orig.get_point(k);
      weights[k] = qf_orig.get_weight(k);
      XYZ[k] = p.x;
      XYZ[k+n_points] = p.y;
    }

    // loop over all needed derivatives
    for(int k = 0; k < 3; k++)
    {
      // loop over all quadrature points
      for(int j = 0; j < n_points; j++) 
      {
        double value_x = 0;
        double value_y = 0;
        // loop over all basis functions
        for(int l = 0; l < n_base_functs; l++)
        {
          value_x += FEFunctValues[l] * AllOrigValues[k][j][l];
          value_y += FEFunctValues[l] * AllOrigValues[k][j][l+n_base_functs];
        }
        Derivatives[j][k] = value_x;
        Derivatives[j][k+3] = value_y;
      }
    }
    double hK=1; // we set it one here since it is not needed in 
    // ErrorMeth=L2H1Errors
    double LocError[3]; // L^2 error in value, divergence and first derivative
    ErrorMeth(n_points, {{XYZ.data(), XYZ.data()+n_points}}, AbsDetjk.data(),
              weights.data(), hK, Derivatives, ExactVal.data(), nullptr,
              LocError);
    for(int j = 0; j < 3; j++) 
    {
      errors[j] += LocError[j];
    }
    // delete everything which was created with "new" within this loop
    // otherwise one would get (many) memory leaks
    for (int j = 0; j < n_points; j++)
    {
      delete [] ExactVal[j]; ExactVal[j] = nullptr;
      delete [] Derivatives[j]; Derivatives[j] = nullptr;
    }
    
    delete[] Derivatives;
    delete[] FEFunctValues;
    delete[] uref;
    delete[] uxiref;
    delete[] uetaref;
    
    for(int iOrig = 0; iOrig < 3; iOrig++)
    {
      for(int jOrig = 0; jOrig < n_points; jOrig++)
        delete[] AllOrigValues[iOrig][jOrig];
      
      delete[] AllOrigValues[iOrig];
    }
    delete[] AllOrigValues;
    
  } // end loop over all cells

  for(int j = 0; j < 3; j++)  
    errors[j] = std::sqrt(errors[j]);
}


void TFEFunction2D::Interpolate(TFEFunction2D *OldFeFunction)
{
  int i,j, N_Cells;
  int N_LocalDOFs;
  int N_Points;

  const double *xi, *eta;
  double X[MaxN_PointsForNodal2D], Y[MaxN_PointsForNodal2D];
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_PointsForNodal2D];
  double values[4];


  const TCollection* Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();


  for(i=0;i<N_Cells;i++)
  {
    TBaseCell* cell = Coll->GetCell(i);
    const FiniteElement& Element = FESpace2D->get_fe(i);
    const NodalFunctional* nf = Element.GetNodalFunctional();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element.GetN_DOF();

    ReferenceTransformation_type RefTrans = Element.GetRefTransID();

    
    bool IsIsoparametric = TDatabase::ParamDB->USE_ISOPARAMETRIC
                           && cell->has_isoparametric_joint();
    if(IsIsoparametric)
    {
      BFRefElements RefElement = Element.GetBaseFunct()->GetRefElement();
      
      switch(RefElement)
      {
        case BFRefElements::BFUnitSquare:
          RefTrans = ReferenceTransformation_type::QuadIsoparametric;
          break;

        case BFRefElements::BFUnitTriangle:
          RefTrans = ReferenceTransformation_type::TriaIsoparametric;
          break;
          
        default:
          ;
      }
    }
    

    FEDatabase::SetCellForRefTrans(cell, RefTrans);
    FEDatabase::GetOrigFromRef(RefTrans, N_Points, xi, eta, X, Y);

    for(j=0;j<N_Points;j++)
    {
      OldFeFunction->FindGradient(X[j], Y[j], values);
      PointValues[j] = values[0]; 
    }

    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValues,
        FunctionalValues);

    const int* DOF = FESpace2D->GetGlobalDOF(i);

    for(j=0;j<N_LocalDOFs;j++)
    {
      Values[DOF[j]] = FunctionalValues[j];
    }
  } //for i
}

/** ************************************************************************* */
void TFEFunction2D::add(AnalyticFunction f)
{
  int N_Points;
  const double *xi, *eta;
  // begin code

  auto Coll = FESpace2D->GetCollection();
  int N_Cells = Coll->GetN_Cells();
  int N_DOFs = FESpace2D->get_n_dof();
  std::vector<int> IntIndex(N_DOFs, 0);

  for(int i = 0; i < N_Cells; i++)
  {
    const TBaseCell * cell = Coll->GetCell(i);
    const FiniteElement& Element = FESpace2D->get_fe(i);
    auto nf = Element.GetNodalFunctional();
    nf->GetPointsForAll(N_Points, xi, eta);
    int N_LocalDOFs = Element.GetN_DOF();
    
    if(Element.GetBaseFunct()->GetBaseVectDim() != 1)
      ErrThrow("TFEFunction2D::add not implemented for vector-valued basis "
               "functions");

    BFRefElements RefElement = Element.GetBaseFunct()->GetRefElement();
    ReferenceTransformation_type RefTrans = Element.GetRefTransID();
    bool IsIsoparametric = TDatabase::ParamDB->USE_ISOPARAMETRIC
                           && cell->has_isoparametric_joint();
    if(IsIsoparametric)
    {
      RefTrans = (RefElement == BFRefElements::BFUnitSquare) 
                 ? ReferenceTransformation_type::QuadIsoparametric
                 : ReferenceTransformation_type::TriaIsoparametric;
    }

    FEDatabase::SetCellForRefTrans(cell, RefTrans);
    std::vector<double> X(N_Points, 0.);
    std::vector<double> Y(N_Points, 0.);
    FEDatabase::GetOrigFromRef(RefTrans, N_Points, xi, eta, X.data(), Y.data());
    std::vector<double> PointValues(N_Points, 0.);

    for(int j = 0; j < N_Points; j++)
    {
      PointValues[j] = f(cell, i, {{X[j], Y[j]}});
    }
    std::vector<double> FunctionalValues(N_LocalDOFs, 0.);
    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValues.data(),
        FunctionalValues.data());

    const int * DOF = FESpace2D->GetGlobalDOF(i);

    for(int j = 0; j < N_LocalDOFs; j++)
    {
      if(IntIndex[DOF[j]] == 0)
        Values[DOF[j]] += FunctionalValues[j];
      IntIndex[DOF[j]] ++;
    }
  }

  for(int i = 0; i < N_DOFs; i++)
  {
    if(IntIndex[i] == 0)
    {     
      ErrThrow("TFEFunction2D::add: unable to set dof ", i);
    }
  }
}

/** ************************************************************************* */
void TFEFunction2D::add_constant(double b)
{
  int length = this->GetLength();
  // vector of the same length as this TFEFunction2D. It represents a function
  // which has the constant value 'mean' for all nodal functionals. The last 
  // step in this projection will be to substract this vector from the vector of
  // this TFEFunction2D
  // for standard P_k or Q_k finite elements this is a constant function
  std::vector<double> interpol(length);

  auto coll = FESpace2D->GetCollection();
  const int n_cells = coll->GetN_Cells();
  for(int i = 0; i < n_cells; i++)
  {
    TBaseCell *cell = coll->GetCell(i); // current cell
    auto element = FESpace2D->get_fe(i);
    // finite element on the current cell
    const int n_loc_dof = element.GetN_DOF(); // number of local dofs
    const int * DOF = FESpace2D->GetGlobalDOF(i);
    const NodalFunctional *nf = element.GetNodalFunctional();
    int n_points; // number of evaluation points to compute nodal functionals
    const double *xi, *eta; //coordinates of evaluation points in reference cell
    nf->GetPointsForAll(n_points, xi, eta);
    double *point_values = new double[n_points];
    for(int j = 0; j < n_points; j++)
      point_values[j] = b;
    // evaluate nodal functionals
    double *functional_values = new double[n_loc_dof];
    nf->GetAllFunctionals(coll, cell, point_values, functional_values);
    for(int j = 0; j < n_loc_dof; j++)
      interpol[DOF[j]] = functional_values[j];

    delete [] point_values;
    delete [] functional_values;
  }

  // the vector 'interpol' is now complete
  // substract it from the vector of this TFEFunction2D
  for(int i = 0; i < length; i++)
    Values[i] += interpol[i];
}

/** ************************************************************************* */
void TFEFunction2D::compute_integral_and_measure(double& integral,
    double& measure) const
{
  integral = 0.0; // variable to store integral value of this TFEFunction2D
  measure = 0.0; // variable to store the measure of the domain

  // loop over all cells, find out integral value of this FEFunction2D and the`
  // measure of its domain
  for(int i = 0; i < FESpace2D->GetCollection()->GetN_Cells(); i++)
  {
    double integral_local;
    double measure_local;
    this->compute_integral_and_measure(i, integral_local, measure_local);
    integral += integral_local;
    measure += measure_local;
  }
}

/** ************************************************************************* */
void TFEFunction2D::compute_integral_and_measure(int cell_i, double& integral,
    double& measure) const
{
  
  // find out integral value of this FEFunction2D and the measure of cell_i
  integral = 0.0; // variable to store integral value of this TFEFunction2D
  measure = 0.0; // variable to store the measure of the domain

  // calculate values on original element (i.e. prepare reference
  // transformation)
  auto coll = FESpace2D->GetCollection();
  auto cell = coll->GetCell(cell_i); // current cell
  auto element = FESpace2D->get_fe(cell_i);

  bool IsIsoparametric = TDatabase::ParamDB->USE_ISOPARAMETRIC
                         && cell->has_isoparametric_joint();
  
  auto degree = element.GetBaseFunct()->GetPolynomialDegree();
  auto ref_ele = element.GetBaseFunct()->GetRefElement();
  auto qf_ref = QuadratureFormulaDatabase::qf_from_degree(degree, ref_ele);
  TQuadFormula qf_orig(*qf_ref);

  auto ref_trans_id = FESpace2D->get_fe(cell_i).GetRefTransID();
  if(IsIsoparametric)
  {
    switch(ref_ele)
    {
      case BFRefElements::BFUnitSquare:
        ref_trans_id = ReferenceTransformation_type::QuadIsoparametric;
        break;

      case BFRefElements::BFUnitTriangle:
        ref_trans_id = ReferenceTransformation_type::TriaIsoparametric;
        break;
        
      default:
        ;
    }
  } // endif IsIsoparametric
  auto ref_trans2D = FEDatabase::GetRefTrans2D(ref_trans_id);
  switch(ref_trans_id)
  {
    case ReferenceTransformation_type::TriaAffin:
    case ReferenceTransformation_type::QuadAffin:
    case ReferenceTransformation_type::QuadBilinear:
      break;
    case ReferenceTransformation_type::TriaIsoparametric:
      ((TTriaIsoparametric *)ref_trans2D)->SetApproximationOrder(degree);
      ((TTriaIsoparametric *)ref_trans2D)->SetQuadFormula(qf_ref->get_type());
      break;
    case ReferenceTransformation_type::QuadIsoparametric:
      ((TQuadIsoparametric *)ref_trans2D)->SetApproximationOrder(degree);
      ((TQuadIsoparametric *)ref_trans2D)->SetQuadFormula(qf_ref->get_type());
      break;
      
    default:
      ;
  } // endswitch
  ref_trans2D->SetCell(cell);
  ref_trans2D->GetOrigFromRef(*qf_ref, qf_orig);

  // local integration (loop over all quadrature points)
  for(int point_j = 0; point_j < qf_orig.GetN_QuadPoints(); point_j++)
  {
    // local transformed values on this quadrature point
    double value; // value of this TFEFunction2D at this quadrature point
    auto pt = qf_orig.get_point(point_j);
    this->FindValueLocal(cell, cell_i, pt.x, pt.y, &value);

    auto weight = qf_orig.get_weight(point_j);
    integral += weight * value;
    measure += weight;
  } // endfor point_j
}


/** ************************************************************************* */
void TFEFunction2D::compute_edge_integral_and_length(int cell_i, int edge_j,
    double& integral, double& length) const
{
  // find local used elements on this cell
  auto finite_element = this->GetFESpace()->get_fe(cell_i);
  auto base_functions = finite_element.GetBaseFunct();
  auto degree = base_functions->GetPolynomialDegree();

  auto quad_form_1D = QuadratureFormulaDatabase::qf_from_degree(
      degree, BFRefElements::BFUnitLine);

  // Find the coordinates of the quadrature points on original cell and store
  // in x, y
  // Transform 1D quadrature rule to 2D on reference cell and compute x, y
  auto n_quad_pts = quad_form_1D->GetN_QuadPoints();
  std::vector<double> xi(n_quad_pts);
  std::vector<double> eta(n_quad_pts);
  auto ref_element = base_functions->GetRefElement();
  for(int quad_pt_i = 0; quad_pt_i < n_quad_pts; quad_pt_i++)
  {
    auto zeta_1D = quad_form_1D->get_point(quad_pt_i).x;
    // transform 1D quad point to 2D
    auto xi_eta_2D = transform(ref_element, edge_j, zeta_1D);
    xi[quad_pt_i] = xi_eta_2D.x;
    eta[quad_pt_i] = xi_eta_2D.y;
  } // endfor quad_pt_i
  // transform quad points to original cell
  std::vector<double> x(n_quad_pts);
  std::vector<double> y(n_quad_pts);
  auto coll = FESpace2D->GetCollection();
  auto cell = coll->GetCell(cell_i);
  FEDatabase::GetRefTrans2D(finite_element.GetRefTransID())->SetCell(cell);
  FEDatabase::GetOrigFromRef(finite_element.GetRefTransID(), n_quad_pts,
                             xi.data(), eta.data(), x.data(), y.data());

  // find out integral value of this FEFunction2D and the measure of the edge
  integral = 0.0;
  length = 0.0;

  // local integration (loop over all quadrature points)
  for(int point_j = 0; point_j < n_quad_pts; point_j++)
  {
    // local transformed values on this quadrature point
    double value; // value of this TFEFunction2D at this quadrature point
    this->FindValueLocal(cell, cell_i, x[point_j], y[point_j], &value);

    auto weight = quad_form_1D->get_weight(point_j);
    integral += weight * value;
  } // endfor point_j
  double x0, y0, x1, y1, z0;
  cell->GetVertex(edge_j)->GetCoords(x0, y0, z0);
  cell->GetVertex((edge_j+1) % cell->GetN_Edges())->GetCoords(x1, y1, z0);
  length = std::sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
  integral *= length / 2;
}


/** ************************************************************************* */
double TFEFunction2D::compute_mean() const
{
  double integral, measure;
  this->compute_integral_and_measure(integral, measure);

  return integral/measure;
}

/** ************************************************************************* */
void TFEFunction2D::project_into_L20(double a)
{
  // compute current integral and measure of the domain:
  double integral, measure;
  this->compute_integral_and_measure(integral, measure);
  double new_mean = (integral - a)/measure;

  // substract the new mean to this FEFunction2D
  this->add_constant(-new_mean);
}


TFEFunction2D& TFEFunction2D::operator*=(double alpha)
{
  int N_Active = FESpace2D->get_n_active();
  for (int i=0; i<N_Active; i++)
  {
    Values[i] *= alpha;
  }
  return *this;
}

TFEFunction2D & TFEFunction2D::operator+=(const TFEFunction2D & rhs)
{
  if(FESpace2D != rhs.FESpace2D)
  {
    ErrThrow("ERROR: TFEFunction2D::operator+=() The two arguments ",
             "have different fe spaces. Exiting");
  }
  if(Values == rhs.Values)
  {
    ErrThrow("ERROR: TFEFunction2D::operator+=() The two arguments ",
             "have the same solution vector. This operation would be ",
             "equivalent to a multiplication by 2! Exiting");
  }
  int N_Active = FESpace2D->get_n_active();
  for (int i=0; i<N_Active; i++)
  {
    Values[i] += rhs.Values[i];
  }
  return *this;
}

TFEFunction2D & TFEFunction2D::operator=(const TFEFunction2D & other)
{
  this->Name        = other.Name;
  this->FESpace2D   = other.FESpace2D;
  this->Values      = other.Values;

  return *this;
}

// this function computes the values of the basis functions on a grid on the 
// reference cell. When using these values it is important to go through them in
// exactly the same order as here.
// nppd = number of points per dim
void update_data(std::vector<double>& bf_values, int nppd,
                 BaseFunction_type& current_type, const BaseFunctions* bf)
{
  if(current_type == bf->GetID())
    return;
  current_type = bf->GetID();
  int tria = (bf->GetRefElement() == BFRefElements::BFUnitTriangle);
  auto reference_element = bf->GetRefElement();
  int n_base_functs = bf->GetDimension();
  int n_points = tria ? ((nppd+1) * (nppd+2)/2) : (nppd+1) * (nppd+1);
  bf_values.resize(n_points * n_base_functs);
  if(reference_element == BFRefElements::BFUnitSquare)
  {
    int index = 0;
    for(int ix = 0; ix <= nppd; ++ix)
    {
      double x = -1. + ix*2./nppd;
      for(int iy = 0; iy <= nppd; ++iy)
      {
        double y = -1. + iy*2./nppd;
        bf->GetDerivatives(MultiIndex2D::D00, x, y,
                           &bf_values[index*n_base_functs]);
        index++;
      }
    }
  }
  else if(reference_element == BFRefElements::BFUnitTriangle)
  {
    int index = 0;
    for(int ix = 0; ix <= nppd; ++ix)
    {
      double x = (double)ix / nppd;
      for(int iy = 0; iy+ix <= nppd; ++iy)
      {
        double y = (double)iy / nppd;
        bf->GetDerivatives(MultiIndex2D::D00, x, y,
                           &bf_values[index*n_base_functs]);
        index++;
      }
    }
  }
  else
  {
    ErrThrow("Reference element ", reference_element,
             " not supported here.");
  }
}


// Compute minimal and maximal value of nodal functionals of a Pk element in
// a given cell evaluated at an FE-function, where the Pk element are of the
// same order as the given FE-function. Attention: This works only correctly if
// the reference transformation preserves point values.
std::pair<double, double> compute_cell_min_max_nodal_fctn(const
    TFEFunction2D& fe_function, int cell_i)
{
  auto fespace = fe_function.GetFESpace();
  std::pair<double, double> minmax(1e100, -1e100);
  // Here the nodal functionals of the Pk elements are used that
  // are of the same order to compute the minimum and maximum value in each
  // cell. Due to the Kronecker-Delta property of the Pk elements finding the
  // minimum and maximum corresponds to finding the minimal and maximal value
  // of the FE-vector in the Pk basis. This is in particular easy if the
  // FESpace itself is Pk space. Otherwise, the local FE-vector has to be
  // transformed.
  auto is_pk_space = !fespace->is_discontinuous();
  auto values = fe_function.GetValues();
  auto global_dofs = fespace->GetGlobalDOF(cell_i);
  if (is_pk_space)
  {
    auto n_loc_dof = fespace->get_n_local_dof(cell_i);
    for (unsigned int dof_j = 0; dof_j < n_loc_dof; ++dof_j)
    {
      auto val = values[global_dofs[dof_j]];
      minmax.first = std::min(minmax.first, val);
      minmax.second = std::max(minmax.second, val);
    }
  }
  else
  {
    // The idea here is to locally map the polynomial basis functions of the
    // discontinuous elements onto the basis functions of a sufficiently large
    // continuous finite element and determine the min (max) value in this basis
    // by evaluating the min (max) of the vector.
    // For this a transformation matrix on the reference cell is computed. The
    // FE-vector is then computed by multiplying the FE-vector to the matrix.
    // Actually the values have to be transformed from the reference cell to the
    // original cell, but for value-preserving transformation this is not
    // needed. Therefore, this method works not correctly for instance for Hdiv
    // elements.

    // Get the matrix. This information can be gathered from
    // FEDatabase::GetProlongationMatrix2D
    auto fe_original = fespace->get_fe(cell_i);
    auto fe_type_original = fe_original.GetID();
    auto degree = fe_original.GetBaseFunct()->GetPolynomialDegree();
    auto ele_for_shape = TFESpace2D::get_element_for_shape(degree, 2);
    auto shape = fespace->GetCollection()->GetCell(cell_i)->GetType();
    auto element = ele_for_shape[shape];

    auto matrix = FEDatabase::GetProlongationMatrix2D(fe_type_original, NoRef,
                                                      element, 0);

    // Apply matrix vector multiplication on reference cell. Transformation back
    // to the original cell is not used here, see also above.
    auto n_loc_dof = fespace->get_n_local_dof(cell_i);
    for (unsigned int dof_i = 0; dof_i < n_loc_dof; ++dof_i)
    {
      double val = 0;
      for (unsigned int dof_j = 0; dof_j < n_loc_dof; ++dof_j)
      {
        val += matrix[dof_j + dof_i * MaxN_BaseFunctions2D]
          * values[global_dofs[dof_j]];
      }
      minmax.first = std::min(val, minmax.first);
      minmax.second = std::max(val, minmax.second);
    }
  }
  if (minmax.first > minmax.second)
  {
    Output::warn<1>("TFEFunction2D::compute_cell_wise_min_max_function_eval",
        "Method was not successful on cell ", cell_i);
  }
  return minmax;
}

// Compute minimal and maximal value of uh using point evaluations
std::pair<double, double> compute_cell_min_max_function_eval(const
    TFEFunction2D& fe_function, std::vector<double>& bf_values, int
    n_points_per_dim, BaseFunction_type& current_type, const BaseFunctions* bf,
    int cell_i)
{
  int n_base_functs = bf->GetDimension();
  update_data(bf_values, n_points_per_dim, current_type, bf);
  int tria = (bf->GetRefElement() == BFRefElements::BFUnitTriangle);
  int index = 0;
  auto global_dofs = fe_function.GetFESpace()->GetGlobalDOF(cell_i);
  auto values = fe_function.GetValues();
  std::pair<double, double> minmax(1e100, -1e100);
  for(int ix = 0; ix <= n_points_per_dim; ++ix)
  {
    for(int iy = 0; iy+ix*tria <= n_points_per_dim; ++iy)
    {
      double val = 0.;
      for(int i = 0; i < n_base_functs; ++i)
      {
        int global_dof = global_dofs[i];
        val += bf_values[index*n_base_functs + i] * values[global_dof];
      }
      index++;
      minmax.first = std::min(minmax.first, val);
      minmax.second = std::max(minmax.second, val);
    }
  }
  if (minmax.first > minmax.second)
  {
    Output::warn<1>("TFEFunction2D::compute_cell_wise_min_max_function_eval",
        "Method was not successful on cell ", cell_i);
  }
  return minmax;
}


std::pair<double, double> TFEFunction2D::compute_cell_min_max(int cell_nr,
    std::vector<double>& bf_values, BaseFunction_type& current_type,  bool
    use_pk_nodal_fctn) const
{
  // Gather some information used in the case use_nodal_fctn == false
  constexpr int n_points_per_dim = 41;

  // Find local minimum and maximum over cell
  std::pair<double, double> minmax;
  if (use_pk_nodal_fctn)
  {
    minmax = compute_cell_min_max_nodal_fctn(*this, cell_nr);
  }
  else
  {
    auto bf = FESpace2D->get_fe(cell_nr).GetBaseFunct();
    minmax = compute_cell_min_max_function_eval(*this, bf_values,
        n_points_per_dim, current_type, bf, cell_nr);
  }
  // Update min and max accordingly
  return minmax;
}

/** ************************************************************************* */
void TFEFunction2D::MinMax(double & min, double & max, const bool&
    use_pk_nodal_fctn) const
{
  //double t = GetTime();
  if(this->GetFESpace()->GetBaseVectDim() != 1)
  {
    ErrThrow("computing the minimum and maximum of a vector valued function is "
        "not implemented.");
  }

  min = 1e100;
  max = -1e100;

  // Find local minimum and maximum over all cells and update the global minimum
  // and maximum accordingly
  auto n_cells = GetFESpace()->GetCollection()->GetN_Cells();
  std::vector<double> bf_values; // values of basis functions on reference cell
  BaseFunction_type current_type = BF_C_L_P0_1D; // dummy type
  for (int cell_i = 0; cell_i < n_cells; ++cell_i)
  {
    auto minmax = this->compute_cell_min_max(cell_i, bf_values, current_type,
        use_pk_nodal_fctn);
    // Update min and max accordingly
    min = std::min(min, minmax.first);
    max = std::max(max, minmax.second);
  }

  ////t = GetTime() - t;
  ////Output::print<3>("time for min/max computation: ", t, " seconds");
}

/** ************************************************************************* */
void TFEFunction2D::PrintMinMax(const std::string& name, const bool&
    use_pk_nodal_fctn) const
{
  double min, max;
  this->MinMax(min, max, use_pk_nodal_fctn);
  if( min <= max )
  {
    if(name.empty())
      Output::print<1>(this->Name, " min ", std::setprecision(12), min, " max ", max);
    else
      Output::print<1>(name, " min ", std::setprecision(12), min, "  max ", max);
  }
  else
  {
    Output::print<1>("WARNING: TFEFunction2D::MinMax was not successful!");
  }
}

/** ************************************************************************* */
double TFEFunction2D::compute_mean_oscillation(double u_min, double u_max, const
    bool& use_pk_nodal_fctn) const
{
  if(this->GetFESpace()->GetBaseVectDim() != 1)
  {
    ErrThrow("computing the minimum and maximum of a vector valued function is "
        "not implemented.");
  }

  // Gather some information used in the case use_nodal_fctn == false
  constexpr int n_points_per_dim = 41;
  BaseFunction_type current_type = BF_C_L_P0_1D; // dummy type
  std::vector<double> bf_values; // values of basis functions on reference cell

  // Find local minimum and maximum over all cells and update the mean
  // oscillation accordingly
  double mean_osc = 0;
  auto n_cells = GetFESpace()->GetCollection()->GetN_Cells();
  for (int cell_i = 0; cell_i < n_cells; ++cell_i)
  {
    std::pair<double, double> minmax;
    if (use_pk_nodal_fctn)
    {
      minmax = compute_cell_min_max_nodal_fctn(*this, cell_i);
    }
    else
    {
      auto bf = FESpace2D->get_fe(cell_i).GetBaseFunct();
       minmax = compute_cell_min_max_function_eval(*this, bf_values,
          n_points_per_dim, current_type, bf, cell_i);
    }
    // Update min and max accordingly
    mean_osc += std::max(0., minmax.second - u_max) + std::max(0., u_min
        - minmax.first);
  }
  mean_osc /= n_cells;
  return mean_osc;
}

/** ************************************************************************* */
/** set Dirichlet values according to boundary conditions */
void TFEFunction2D::SetDirichletBC(BoundCondFunct2D *BoundaryCondition,
    BoundValueFunct2D *BoundaryValue)
{
  const TBoundComp2D* BoundComp;
  BoundCond Cond0, Cond1;
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_BaseFunctions2D];
  double t0, t1, s, t;
  int comp;
  int i, l, m;
  int N_Cells, N_Joints, N_EdgeDOF;
  int *EdgeDOF, N_EdgePoints;
  const double *EdgePoints;
  double eps=1e-4;

  auto Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    auto cell  = Coll->GetCell(i);
    auto Element = FESpace2D->get_fe(i);
    auto nf = Element.GetNodalFunctional();
    auto FEDesc_Obj = Element.GetFEDesc();

    nf->GetPointsForEdge(N_EdgePoints, EdgePoints);

    auto DOF = FESpace2D->GetGlobalDOF(i);

    N_Joints = cell->GetN_Edges();
    for(m=0;m<N_Joints;m++)
    {
      auto joint = cell->GetJoint(m);
      if(joint->GetType() == BoundaryEdge ||
          joint->GetType() == IsoBoundEdge)
      {
        if(joint->GetType() == BoundaryEdge)
        {
          auto boundedge = (const TBoundEdge *)joint;
          BoundComp = boundedge->GetBoundComp();
          boundedge->GetParameters(t0, t1);
        }
        else
        {
          auto isoboundedge = (const TIsoBoundEdge *)joint;
          BoundComp = isoboundedge->GetBoundComp();
          isoboundedge->GetParameters(t0, t1);
        }
        // get id of the boundary component
        comp=BoundComp->GetID();
        // get type of the boundary condition at the beginning
        // and at the end of the current edge
        if (t0 < t1)
        {
          BoundaryCondition(comp, t0+eps, Cond0);
          BoundaryCondition(comp, t1-eps, Cond1);
        }
        else
        {
          BoundaryCondition(comp, t0-eps, Cond0);
          BoundaryCondition(comp, t1+eps, Cond1);
        }
        // only one boundary condition per edge allowed
        if(Cond0 == Cond1)
        {
          if(Cond0 == DIRICHLET)
          {
            // read boundary values for each quadrature point
            for(l=0;l<N_EdgePoints;l++)
            {
              s = EdgePoints[l];
              t = 0.5*(t0*(1-s) + t1*(1+s));
              BoundaryValue(comp, t, PointValues[l]);
            } // endfor l
            // compute boundary values for each dof on the 
            // boundary edge with the nodal functionals
            nf->GetEdgeFunctionals(Coll, cell, m, PointValues,
                FunctionalValues);
            EdgeDOF = FEDesc_Obj->GetJointDOF(m);
            N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();
            // save boundary values of each dof on the boundary
            // edge in the rhs
            for(l=0;l<N_EdgeDOF;l++)
            {
              Values[DOF[EdgeDOF[l]]] = FunctionalValues[l];
              // cout << i <<  setw(25) << Values[DOF[EdgeDOF[l]]]<< endl;
            }
          } // endif Cond0==DIRICHLET
        } // endif (Cond0==Cond1)
        else
        {
          ErrThrow("different boundary condition on one edge are not allowed!");
        }
      } // endif (boundary joint)
    } // endfor m<N_Joints
  } // endfor i<N_Cells
}


double TFEFunction2D::compute_cell_average(const int& cell_nr) const
{
  double integral;
  double measure;
  this->compute_integral_and_measure(cell_nr, integral, measure);
  double cell_av = integral / measure;
  return cell_av;
}

void TFEFunction2D::compute_cell_averages(double* cell_averages)
{
  auto n_cells = this->GetFESpace()->GetCollection()->GetN_Cells();
  for (int cell_i = 0; cell_i < n_cells; ++cell_i)
  {
    cell_averages[cell_i] = this->compute_cell_average(cell_i);
  }
}

void TFEFunction2D::getOrigValuesForCell(int cellNumber, double x, double y, double* uorig, double* uxorig, double* uyorig) const
{
  const FiniteElement FE_Obj = FESpace2D->get_fe(cellNumber);
  ReferenceTransformation_type refTrans = FE_Obj.GetRefTransID();
  TBaseCell* cellPtr = FESpace2D->GetCollection()->GetCell(cellNumber);
  
  bool isIsoparametric = TDatabase::ParamDB->USE_ISOPARAMETRIC
                          && cellPtr->has_isoparametric_joint();
  if(isIsoparametric)
  {
    BFRefElements refElement = FE_Obj.GetBaseFunct()->GetRefElement();
    switch(refElement)
    {
      case BFRefElements::BFUnitSquare:
        refTrans = ReferenceTransformation_type::QuadIsoparametric;
        break;

      case BFRefElements::BFUnitTriangle:
        refTrans = ReferenceTransformation_type::TriaIsoparametric;
        break;
        
      default:
        ;
    }
  }     
  
  FEDatabase::SetCellForRefTrans(cellPtr, refTrans);
  
  double xi, eta;
  FEDatabase::GetRefFromOrig(refTrans, x, y, xi, eta);
  
  const BaseFunctions* bf = FE_Obj.GetBaseFunct();
  int N_BaseFunct = bf->GetDimension();
  
  double *uref = new double[N_BaseFunct];
  double *uxiref = new double[N_BaseFunct];
  double *uetaref = new double[N_BaseFunct];

  bf->GetDerivatives(MultiIndex2D::D00, xi, eta, uref);
  bf->GetDerivatives(MultiIndex2D::D10, xi, eta, uxiref);
  bf->GetDerivatives(MultiIndex2D::D01, xi, eta, uetaref);
  
  FEDatabase::GetOrigValues(refTrans, xi, eta, bf, FESpace2D->GetCollection(),
                            (TGridCell *)cellPtr, uref, uxiref, uetaref,
                            uorig, uxorig, uyorig);
  
  delete[] uref;
  delete[] uxiref;
  delete[] uetaref;
}

