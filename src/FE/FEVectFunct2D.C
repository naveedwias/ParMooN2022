#ifdef _MPI
# include "mpi.h"
#endif

#include <FEVectFunct2D.h>
#include <FEDatabase.h>
#include <NodalFunctional.h>
#include <QuadIsoparametric.h>
#include <TriaIsoparametric.h>
#include <Database.h>
#include <AuxParam2D.h>
#include <GridCell.h>

#include <Joint.h>
#include <BoundEdge.h>
#include <InterfaceJoint.h>

#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <MooNMD_Io.h>
// #include <malloc.h>
#include <dirent.h> 
#include <cstring> // memset
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <BoundaryAssembling2D.h>
#include "QuadratureFormulaDatabase.h"

/// Default constructor. Constructs an empty object.
TFEVectFunct2D::TFEVectFunct2D()
{
  N_Components = 0;
}


/** constructor with vector initialization */
TFEVectFunct2D::TFEVectFunct2D(std::shared_ptr<const TFESpace2D> fespace2D,
                               const std::string& name, double *values,
                               int n_components)
  : TFEFunction2D(fespace2D, name, values)
{
  N_Components = n_components;
}

TFEVectFunct2D& TFEVectFunct2D::operator=( const TFEVectFunct2D & other)
{
  //call base class copy assignment
  TFEFunction2D::operator=(other);

  this->N_Components  = other.N_Components;

  return *this;
}

//====================================================================
/** calculate errors to given vector function */
void TFEVectFunct2D::GetDeformationTensorErrors( 
  DoubleFunct2D *Exact, DoubleFunct2D *Exact1,
  int N_Derivatives,
  MultiIndex2D *NeededDerivatives,
  int N_Errors, TFEFunction2D::ErrorMethod *ErrorMeth, 
  const CoeffFct2D& Coeff, 
  TAuxParam2D *Aux,
  int n_fespaces, TFESpace2D **fespaces,
  double *errors)
{
  int i,j,k,l;
  int N_Cells, N_Points, N_Parameters, N_;  //N_U;
  int Used[N_FEs2D];
  TFESpace2D *fespace;
  const double *weights;
  const double *X, *Y;
  double *Param[MaxN_QuadPoints_2D], *aux, *aux1, *aux2, *aux3;
  std::vector<double> AbsDetjk(MaxN_QuadPoints_2D, 1.);
  double *Derivatives[2*MaxN_QuadPoints_2D];
  double *ExactVal[2*MaxN_QuadPoints_2D];
  double *AuxArray[MaxN_QuadPoints_2D];
  double **OrigFEValues, *Orig, value, value1;
  double FEFunctValues[MaxN_BaseFunctions2D];
  double FEFunctValues1[MaxN_BaseFunctions2D];
  double LocError[4], *Values0,*Values1;
  double hK;
  bool *SecondDer;

  SecondDer = new bool[n_fespaces];
  for(i=0;i<n_fespaces;i++)
    SecondDer[i] = false;

  N_Parameters = Aux->GetN_Parameters();
  aux1 = new double [MaxN_QuadPoints_2D*N_Parameters];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Param[j] = aux1 + j*N_Parameters;

  aux2 = new double [2*MaxN_QuadPoints_2D*N_Derivatives];
  for(j=0;j<2*MaxN_QuadPoints_2D;j++)
    Derivatives[j] = aux2 + j*N_Derivatives;
  
  aux3 = new double [2*MaxN_QuadPoints_2D * 4];
  for(j=0;j<2*MaxN_QuadPoints_2D;j++)
    ExactVal[j] = aux3 + j*4;

  // 20 <= number of term
  aux = new double [MaxN_QuadPoints_2D*20]; 
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    AuxArray[j] = aux + j*20;

  fespace = fespaces[0];
  Values0 = Values;
  Values1 = Values+FESpace2D->get_n_dof();

  for(i=0;i<N_Errors;i++)
    errors[i] = 0.0;

  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);
// ########################################################################
// loop over all cells
// ########################################################################
  auto Coll = fespaces[0]->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);

    hK = cell->GetDiameter();

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    memset(Used, 0, N_FEs2D*sizeof(int));
    std::vector<const FiniteElement*> used_fe(n_fespaces, nullptr);
    for(j=0;j<n_fespaces;j++)
    {
      used_fe[j] = &fespaces[j]->get_fe(i);
    }

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    FEDatabase::GetOrig(used_fe, Coll, cell, SecondDer, qf_ref, qf_orig);
    qf_orig.GetFormulaData(N_Points, weights, X, Y);

    if(N_Parameters>0)
      Aux->GetParameters(qf_orig, i, Param); 

    // calculate all needed derivatives of this FE function
    auto element = fespace->get_fe(i);
    N_ = element.GetN_DOF();

    auto DOF = fespace->GetGlobalDOF(i);
    for(l=0;l<N_;l++)
    {
      FEFunctValues[l] = Values0[DOF[l]];
      FEFunctValues1[l] = Values1[DOF[l]];
    }

    // for all needed derivatives
    for(k=0;k<N_Derivatives;k++)
    {
      OrigFEValues = FEDatabase::GetOrigElementValues(*element.GetBaseFunct(),
                                                      NeededDerivatives[k]);
      // for all quadrature points
      for(j=0;j<N_Points;j++)
      {
        Orig = OrigFEValues[j];
        value = 0;
        value1 = 0;
        for(l=0;l<N_;l++)
        {
          value += FEFunctValues[l] * Orig[l];
          value1 += FEFunctValues1[l] * Orig[l];
        } // endfor l
        Derivatives[j][k] = value;
        Derivatives[j+N_Points][k] = value1;
      } // endfor j
    } // endfor k

    // exact value for first component
    for(j=0;j<N_Points;j++)
      Exact(X[j], Y[j], ExactVal[j]);

    // exact value for second component
    for(j=0;j<N_Points;j++)
      Exact1(X[j], Y[j], ExactVal[j+N_Points]);

    if(Coeff)
      Coeff(N_Points, X, Y, Param, AuxArray);      

    ErrorMeth(N_Points, {{X, Y}}, AbsDetjk.data(), weights, hK, Derivatives, 
              ExactVal, AuxArray, LocError);

    for(j=0;j<N_Errors;j++)
      errors[j] += LocError[j];

  } // endfor i

  for(j=0;j<N_Errors;j++)
  {
    if (errors[j] > 0)
    errors[j] = std::sqrt(errors[j]);
  }

  delete aux;
  delete aux1;
  delete aux2;
  delete aux3;
  delete SecondDer;

} // TFEVectFunct2D::GetDeformationTensorErrors

//================================================================
std::pair<double, double> TFEVectFunct2D::get_L2_norm_divergence_curl() const
{
  std::vector<double> values(2);
  auto f = [](std::vector<double>& v, std::array<double, 8> e)
           {
             v[0] += e[4] + e[7]; // divergence
             v[1] += e[5] - e[6]; // curl
           };
  this->get_functional_value(values, f);
  return {values[0], values[1]};

} // TFEVectFunct2D::get_L2_norm_divergence_curl

//================================================================
void TFEVectFunct2D::get_functional_value(std::vector<double>& values,
  const std::function<void(std::vector<double>&, std::array<double, 8>)>&
    functional)
const
{
  bool SecondDer[1] = { false };
  int N_Points;
  double FEValues0[MaxN_BaseFunctions2D];
  double FEValues1[MaxN_BaseFunctions2D];
  
  auto n_values = values.size();
  std::fill(values.begin(), values.end(), 0.0);

  double *Values0 = Values;
  double *Values1 = Values+FESpace2D->get_n_dof();

  auto Coll = FESpace2D->GetCollection();
  auto N_Cells = Coll->GetN_Cells();
  
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);

  for(int i = 0; i < N_Cells; i++)
  {
    auto cell = Coll->GetCell(i);
    auto element = FESpace2D->get_fe(i);
    std::vector<const FiniteElement*> used_fe(1, &element);
    FEDatabase::GetOrig(used_fe, Coll, cell, SecondDer, qf_ref, qf_orig);
    N_Points = qf_orig.GetN_QuadPoints();

    // calculate all needed derivatives of this FE function
    auto N_Bf = element.GetN_DOF();

    auto DOF = FESpace2D->GetGlobalDOF(i);
    
    for(int j = 0; j < N_Bf; j++)
    {
      int k = DOF[j];
      FEValues0[j] = Values0[k];
      FEValues1[j] = Values1[k];
    }

    auto& bf = *element.GetBaseFunct();
    auto OrigFEValues = FEDatabase::GetOrigElementValues(bf, MultiIndex2D::D00);
    auto OrigFEValuesX = FEDatabase::GetOrigElementValues(bf,MultiIndex2D::D10);
    auto OrigFEValuesY = FEDatabase::GetOrigElementValues(bf,MultiIndex2D::D01);

    // for all quadrature points
    for(int j = 0; j < N_Points; j++)
    {
      double * Orig = OrigFEValues[j];
      double * OrigX = OrigFEValuesX[j];
      double * OrigY = OrigFEValuesY[j];
      std::vector<double> local_values(n_values, 0.0);
      auto p = qf_orig.get_point(j);
      std::array<double, 8> evaluations{{p.x, p.y, 0., 0., 0., 0., 0., 0.}};
      for(int l = 0; l < N_Bf; l++)
      {
        evaluations[2] += Orig[l] * FEValues0[l];
        evaluations[3] += Orig[l] * FEValues1[l];
        evaluations[4] += OrigX[l] * FEValues0[l];
        evaluations[5] += OrigX[l] * FEValues1[l];
        evaluations[6] += OrigY[l] * FEValues0[l];
        evaluations[7] += OrigY[l] * FEValues1[l];
      } // endfor l
      functional(local_values, evaluations);
      double local_weight = qf_orig.get_weight(j);
      for(auto k = 0ul; k < n_values; ++k)
      {
        values[k] += local_weight * (local_values[k] * local_values[k]);
      }
    } // endfor j
  } // endfor i
  for(auto k = 0ul; k < n_values; ++k)
  {
    values[k] = std::sqrt(values[k]);
  }
} // TFEVectFunct2D::get_functional_value


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/** calculate L2-norm of (u-u_h).n-error - written by Laura Blank 01.03.18*/
double TFEVectFunct2D::GetL2NormNormalComponentError(BoundValueFunct2D *Exact_u1, BoundValueFunct2D *Exact_u2, bool rescale_by_h_E)
{
  double *Values0 = Values;
  double *Values1 = Values+FESpace2D->get_n_dof();

  auto Coll = FESpace2D->GetCollection();

  // Initialize Pointer (necessary: Resource allocation is initialization Paradigm)
  double *ExactVal_u1[MaxN_QuadPoints_2D], *ExactVal_u2[MaxN_QuadPoints_2D];
  double *aux1 = new double [MaxN_QuadPoints_2D * 4]{0.};
  double *aux2 = new double [MaxN_QuadPoints_2D * 4]{0.};
  for(int ii = 0; ii < MaxN_QuadPoints_2D; ii++)
  {
    ExactVal_u1[ii] = aux1 + ii*4;
    ExactVal_u2[ii] = aux2 + ii*4;
  }

  double FEFunctValues0[MaxN_BaseFunctions2D], FEFunctValues1[MaxN_BaseFunctions2D]; 
  double final_boundary_error_l2[1];
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
      auto element = FESpace2D->get_fe(cell->GetCellIndex());
 
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
      int BaseVectDim = 1;
      BoundaryAssembling2D::get_original_values(
          element, joint_id, cell, quadPoints, BaseVectDim, uorig,
          u_dx_orig, u_dy_orig, *LineQuadFormula);

      // calculate all needed derivatives of this FE function                 
      int N_Bf = element.GetN_DOF();

      const int *DOF = FESpace2D->GetGlobalDOF(cell->GetCellIndex());
      for(int l = 0; l < N_Bf; l++)
      {
        FEFunctValues0[l] = Values0[DOF[l]];
        FEFunctValues1[l] = Values1[DOF[l]];
      }

      double summing_boundary_error_l2_on_edge = 0;
      double edge_length = boundedge->get_length();
      double reference_edge_length = 2; // [-1,1] is the reference edge here↲   
      // normal vector to this boundary (normalized)
      double n1, n2;
      boundedge->get_normal(n1, n2);

      for(size_t j = 0; j < quadPoints.size(); j++)
      {
        double value = 0;
        double value_u1 = 0;
        double value_u2 = 0;
        double exact_val = 0;
        for(int l = 0; l < N_Bf; l++)
        { // compute u_h|T = \sum \alpha_i \Phi_i 
          value_u1 += FEFunctValues0[l] * uorig[j][l];
          value_u2 += FEFunctValues1[l] * uorig[j][l];
        } 
        value += value_u1 * n1 + value_u2 * n2;

        auto comp = boundedge->GetBoundComp();
        int comp_ID = comp->GetID();
        double t0, t1;
        boundedge->GetParameters(t0, t1);
        double t = t0 + 0.5 * (t1-t0) * (quadPoints[j]+1); 
        Exact_u1(comp_ID, t, ExactVal_u1[j][0]);
        Exact_u2(comp_ID, t, ExactVal_u2[j][0]);
        exact_val += ExactVal_u1[j][0] * n1 + ExactVal_u2[j][0] * n2;

        if(rescale_by_h_E)
        { 
          summing_boundary_error_l2_on_edge += quadWeights[j] * edge_length/reference_edge_length * 1/edge_length  * (exact_val - value) * (exact_val - value);
        }
        else
        {
          summing_boundary_error_l2_on_edge += quadWeights[j] * edge_length/reference_edge_length * (exact_val - value) * (exact_val - value);
        }
      }
      boundary_error_l2_on_Component += summing_boundary_error_l2_on_edge;
    } 
    final_boundary_error_l2[0] +=  boundary_error_l2_on_Component;
  }
  final_boundary_error_l2[0] = std::sqrt(final_boundary_error_l2[0]);
  return final_boundary_error_l2[0];
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/** calculate L2-norm of (u-u_h).n-error without the global database */
double TFEVectFunct2D::GetL2NormNormalComponentError(BoundValueFunct2D *Exact_u1, BoundValueFunct2D *Exact_u2, int boundary_component_id, bool rescale_by_h_E)
{
  double *Values0 = Values;
  double *Values1 = Values+FESpace2D->get_n_dof();
  auto Coll = FESpace2D->GetCollection();

  // Initialize Pointer (necessary: Resource allocation is initialization Paradigm)
  double *ExactVal_u1[MaxN_QuadPoints_2D], *ExactVal_u2[MaxN_QuadPoints_2D];
  double *aux1 = new double [MaxN_QuadPoints_2D * 4]{0.};
  double *aux2 = new double [MaxN_QuadPoints_2D * 4]{0.};
  for(int ii = 0; ii < MaxN_QuadPoints_2D; ii++)
  {
    ExactVal_u1[ii] = aux1 + ii*4;
    ExactVal_u2[ii] = aux2 + ii*4;
  }

  double FEFunctValues0[MaxN_BaseFunctions2D], FEFunctValues1[MaxN_BaseFunctions2D];
  double final_boundary_error_l2[1];
  final_boundary_error_l2[0] = 0;

    // Create a list of those boundary edges that are on the boundary component with given ID
    std::vector<TBoundEdge*> boundaryEdgeList;
    Coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);

    double boundary_error_l2_on_Component = 0;

    for(size_t m = 0; m < boundaryEdgeList.size(); m++)
    {
      TBoundEdge *boundedge = boundaryEdgeList[m];
      TBaseCell *cell = boundedge->GetNeighbour(0);
      auto element = FESpace2D->get_fe(cell->GetCellIndex());

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
      int BaseVectDim = 1;
      BoundaryAssembling2D::get_original_values(
          element, joint_id, cell, quadPoints, BaseVectDim, uorig, u_dx_orig,
          u_dy_orig, *LineQuadFormula);

      // calculate all needed derivatives of this FE function
      int N_Bf = element.GetN_DOF();

      const int *DOF = FESpace2D->GetGlobalDOF(cell->GetCellIndex());
      for(int l = 0; l < N_Bf; l++)
      {
        FEFunctValues0[l] = Values0[DOF[l]];
        FEFunctValues1[l] = Values1[DOF[l]];
      }

      double summing_boundary_error_l2_on_edge = 0;
      double edge_length = boundedge->get_length();
      double reference_edge_length = 2; // [-1,1] is the reference edge here↲
      // normal vector to this boundary (normalized)
      double n1, n2;
      boundedge->get_normal(n1, n2);

      for(size_t j = 0; j < quadPoints.size(); j++)
      {
        double value = 0;
        double value_u1 = 0;
        double value_u2 = 0;
        double exact_val = 0;
        for(int l = 0; l < N_Bf; l++)
        { // compute u_h|T = \sum \alpha_i \Phi_i
          value_u1 += FEFunctValues0[l] * uorig[j][l];
          value_u2 += FEFunctValues1[l] * uorig[j][l];
        }
        value += value_u1 * n1 + value_u2 * n2;

        auto comp = boundedge->GetBoundComp();
        int comp_ID = comp->GetID();
        double t0, t1;
        boundedge->GetParameters(t0, t1);
        double t = t0 + 0.5 * (t1-t0) * (quadPoints[j]+1);
        Exact_u1(comp_ID, t, ExactVal_u1[j][0]);
        Exact_u2(comp_ID, t, ExactVal_u2[j][0]);
        exact_val += ExactVal_u1[j][0] * n1 + ExactVal_u2[j][0] * n2;

        if(rescale_by_h_E)
        {
          summing_boundary_error_l2_on_edge += quadWeights[j] * edge_length/reference_edge_length * 1/edge_length  * (exact_val - value) * (exact_val - value);
        }
        else
        {
          summing_boundary_error_l2_on_edge += quadWeights[j] * edge_length/reference_edge_length * (exact_val - value) * (exact_val - value);
        }
      }
      boundary_error_l2_on_Component += summing_boundary_error_l2_on_edge;
    }
    final_boundary_error_l2[0] =  boundary_error_l2_on_Component;

  return final_boundary_error_l2[0];
}

//==========================================================================
/** calculate L2-norm of divergence error - written by Laura Blank 03.01.18*/
double TFEVectFunct2D::GetL2NormDivergenceError(DoubleFunct2D *Exact_u1,DoubleFunct2D *Exact_u2)
{
  double *Values0 = Values;
  double *Values1 = Values+FESpace2D->get_n_dof();

  auto Coll = FESpace2D->GetCollection();
  int N_Cells = Coll->GetN_Cells();

  // Initialize Pointer (necessary: Resource allocation is initialization Paradigm)
  double *ExactVal_u1[MaxN_QuadPoints_2D];
  double *ExactVal_u2[MaxN_QuadPoints_2D];
  double *aux1 = new double [MaxN_QuadPoints_2D * 4]{0.};
  double *aux2 = new double [MaxN_QuadPoints_2D * 4]{0.};
  for(int ii = 0; ii < MaxN_QuadPoints_2D; ii++)
  {
    ExactVal_u1[ii] = aux1 + ii*4;
    ExactVal_u2[ii] = aux2 + ii*4;
  }

  double loc_div_values_exact_solution = 0.;
  double loc_div_values_FEsolution = 0.;
  double FEFunctValues0[MaxN_BaseFunctions2D];
  double FEFunctValues1[MaxN_BaseFunctions2D];
  double **OrigFEValuesX, *OrigX, value;
  double **OrigFEValuesY, *OrigY;
  bool SecondDer[1] = { false };

  double div_error = 0.;
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);
  // loop over all cells
  for(int i = 0; i < N_Cells; i++)
  {
    auto cell = Coll->GetCell(i);
    auto element = FESpace2D->get_fe(i);

    // compute transformation to reference cell
    std::vector<const FiniteElement*> used_fe(1, &element);
    FEDatabase::GetOrig(used_fe, Coll, cell, SecondDer, qf_ref, qf_orig);
    int N_Points = qf_orig.GetN_QuadPoints();

    // calculate all needed derivatives of this FE function
    int N_Bf = element.GetN_DOF();

    const int *DOF = FESpace2D->GetGlobalDOF(i);
    for(int jj = 0; jj < N_Bf; jj++)
    {
      int k = DOF[jj];
      FEFunctValues0[jj] = Values0[k];
      FEFunctValues1[jj] = Values1[k];
    }
    auto& bf = *element.GetBaseFunct();
    OrigFEValuesX = FEDatabase::GetOrigElementValues(bf, MultiIndex2D::D10);
    OrigFEValuesY = FEDatabase::GetOrigElementValues(bf, MultiIndex2D::D01);

    // loop over all quadrature points
    for(int j = 0; j < N_Points; j++)
    {
      OrigX = OrigFEValuesX[j];
      OrigY = OrigFEValuesY[j];
      value = 0;
      for(int l = 0; l < N_Bf; l++)
      {
        value += FEFunctValues0[l] * OrigX[l] + FEFunctValues1[l] * OrigY[l];
      }
      loc_div_values_FEsolution = value;

      auto p = qf_orig.get_point(j);
      Exact_u1(p.x, p.y, ExactVal_u1[j]);
      Exact_u2(p.x, p.y, ExactVal_u2[j]);
      loc_div_values_exact_solution = ExactVal_u1[j][1] + ExactVal_u2[j][2];

      div_error += qf_orig.get_weight(j) * std::abs(loc_div_values_FEsolution - loc_div_values_exact_solution) * std::abs(loc_div_values_FEsolution - loc_div_values_exact_solution);
    } // endfor j
  } // endfor i
  delete [] aux1;
  delete [] aux2;
  div_error = std::sqrt(div_error);
  return div_error;
} // TFEVectFunct2D::GetL2NormDivergenceError

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/** calculate int_{Gamma_i} u.n on a given boundary component*/
double TFEVectFunct2D::compute_flux(int boundary_component_id) const
{
  // get the FE collection
  auto Coll = FESpace2D->GetCollection();

  // Create a list of those boundary edges that are on the boundary component with given ID
  std::vector<TBoundEdge*> boundaryEdgeList;
  Coll->get_edge_list_on_component(boundary_component_id, boundaryEdgeList);

  
  double flux = 0.;
  for(size_t m = 0; m < boundaryEdgeList.size(); m++)
  {
    
    TBoundEdge *boundedge = boundaryEdgeList[m];
    TBaseCell *cell = boundedge->GetNeighbour(0);
    auto element = FESpace2D->get_fe(cell->GetCellIndex());

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
    int BaseVectDim = 1;
    BoundaryAssembling2D::get_original_values(element, joint_id, cell,
                                              quadPoints, BaseVectDim, uorig,
                                              u_dx_orig, u_dy_orig,
                                              *LineQuadFormula);
    
    // calculate all needed derivatives of this FE function
    double flux_on_edge = 0;
    double edge_length = boundedge->get_length();
    //double reference_edge_length = 2; // [-1,1] is the reference edge here↲
    // normal vector to this boundary (normalized)
    double n1, n2;
    boundedge->get_normal(n1, n2);

    double x0, x1, y0, y1;
    boundedge->get_vertices(x0, y0, x1, y1);
    //Output::print(" v0 = (",x0," , ",y0,"), v1 = (",x1," , ",y1,")");
    double val1[3], val2[3];
    FindVectGradient(x0, y0, val1, val2);
    //Output::print(" u(v0) = (",val1[0]," , ",val2[0],")");
    flux_on_edge = val1[0];
    FindVectGradient(x1, y1, val1, val2);
    //Output::print(" u(v1) = (",val1[0]," , ",val2[0],")");

    flux_on_edge = (flux_on_edge + val1[0])/2.*edge_length;
    
    flux += flux_on_edge;
  }

  return flux;
}

//==========================================================================
void TFEVectFunct2D::Interpolate(TFEVectFunct2D *OldVectFunct)
{
 int N_Points;
 int PolynomialDegree, ApproxOrder;

 const double *xi, *eta;
 double X[MaxN_PointsForNodal2D], Y[MaxN_PointsForNodal2D];
 double PointValues1[MaxN_PointsForNodal2D];
 double PointValues2[MaxN_PointsForNodal2D];
 double FunctionalValues[MaxN_PointsForNodal2D];
 double val1[3], val2[3];
  
  const TCollection* Coll = FESpace2D->GetCollection();
  int N_Cells = Coll->GetN_Cells();
  int N_DOFs = FESpace2D->get_n_dof();
  
  memset(Values, 0, sizeof(double)*N_Components*N_DOFs);

    
  for(int i = 0; i < N_Cells; i++)
  {
    TBaseCell* cell = Coll->GetCell(i);
    const FiniteElement& Element = FESpace2D->get_fe(i);
    const NodalFunctional* nf = Element.GetNodalFunctional();
    nf->GetPointsForAll(N_Points, xi, eta);
    int N_LocalDOFs = Element.GetN_DOF();

    PolynomialDegree = Element.GetBaseFunct()->GetPolynomialDegree();
    ApproxOrder = Element.GetBaseFunct()->GetAccuracy();

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
    }
 

    TRefTrans2D* rt = FEDatabase::GetRefTrans2D(RefTrans);
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
    }
    rt->SetCell(cell);
    rt->GetOrigFromRef(N_Points, xi, eta, X, Y);

    for(int j = 0; j < N_Points; j++)
    {
     OldVectFunct->FindVectGradient(X[j], Y[j], val1, val2);
     PointValues1[j] = val1[0];
     PointValues2[j] = val2[0];
    }

    const int* DOF = FESpace2D->GetGlobalDOF(i);
    
    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValues1, FunctionalValues); 
    for(int j = 0; j < N_LocalDOFs; j++)
      Values[DOF[j]] = FunctionalValues[j];
    
    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValues2, FunctionalValues);     
    for(int j = 0; j < N_LocalDOFs; j++)
      Values[N_DOFs + DOF[j]] = FunctionalValues[j];

      
  } // for i 
} // Interpolate
   
   
/** determine the value of a vect function and its first derivatives at
 the given point */
void TFEVectFunct2D::FindVectGradient(double x, double y, double *val1, double *val2) const
{
  
  double valuesFromCells1[3];
  double valuesFromCells2[3];
  
  int N_Found = 0;
  
  val1[0] = 0;
  val1[1] = 0;
  val1[2] = 0;

  val2[0] = 0;
  val2[1] = 0;
  val2[2] = 0;

  auto Coll = FESpace2D->GetCollection();
  int N_Cells = Coll->GetN_Cells();  

  for(int i = 0; i < N_Cells; i++)
  {
    auto cell = Coll->GetCell(i);
    
    if(cell->PointInCell(parmoon::Point(x,y)))
    {
      N_Found++;
      FindVectGradientLocal(i, x, y, valuesFromCells1, valuesFromCells2);
      
      for(int j = 0; j < 3; j++)
      {
        val1[j] += valuesFromCells1[j];
        val2[j] += valuesFromCells2[j];
      }
    }
  }

  
  if(N_Found>0)
  {
    for(int i = 0; i < 3; i++)
    {
      val1[i] /= N_Found;
      val2[i] /= N_Found;
    }
  }
  else 
  {
    ErrThrow("Point (", x, " , ", y, ") not found !!!!!");
  }
}

void TFEVectFunct2D::FindVectGradientLocal(int cell_no,
                                           double x, double y, double *val1,
                                           double *val2) const
{
  FiniteElement FE_Obj = FESpace2D->get_fe(cell_no);
  
  const BaseFunctions* bf = FE_Obj.GetBaseFunct();
  int N_BaseFunct = bf->GetDimension();
  
  double *uorig = new double[N_BaseFunct];
  double *uxorig = new double[N_BaseFunct];
  double *uyorig = new double[N_BaseFunct];
  
  getOrigValuesForCell(cell_no, x, y, uorig, uxorig, uyorig);
  
  getValuesFromOriginalOnes(val1, 0, N_BaseFunct, cell_no, uorig, uxorig, uyorig);
  getValuesFromOriginalOnes(val2, FESpace2D->get_n_dof(), N_BaseFunct,
                            cell_no, uorig, uxorig, uyorig);
  
  delete [] uorig;
  delete [] uxorig;
  delete [] uyorig;
}

void TFEVectFunct2D::FindValueLocal(const TBaseCell* cell, int cell_no, 
				    double x, double y, double* values) const
{
 this->TFEFunction2D::FindValueLocal(cell, cell_no, x, y, values);
 auto u2 = this->GetComponent(1);
 u2->FindValueLocal(cell, cell_no, x, y, values+1);
}


TFEVectFunct2D& TFEVectFunct2D::operator*=(double alpha)
{
  int N_Active = FESpace2D->get_n_active();
  int Length = FESpace2D->get_n_dof();
  for(int n=0; n<N_Components; n++)
  {
    for (int i=0; i<N_Active; i++)
    {
      Values[i+n*Length] *= alpha;
    }
  }
  return *this;
}

TFEVectFunct2D & TFEVectFunct2D::operator+=(const TFEVectFunct2D & rhs)
{
  if(FESpace2D != rhs.FESpace2D)
  {
    ErrThrow("ERROR: TFEVectFunct2D::operator+=() The two arguments "
             "have different fe spaces. Exiting");
  }
  if(Values == rhs.Values)
  {
    ErrThrow("ERROR: TFEVectFunct2D::operator+=() The two arguments "
             "have the same solution vector. This operation would be "
             "equivalent to a multiplication by 2! Exiting");
  }
  if(N_Components != rhs.N_Components)
  {
    ErrThrow("ERROR: TFEVectFunct2D::operator+=() The two arguments "
             "have different number of components. You are trying to add a ",
             rhs.N_Components, "-dimensional vector to a ", N_Components,
             "-dimensional vector! Exiting");
  }
  int N_Active = FESpace2D->get_n_active();
  int Length = FESpace2D->get_n_dof();
  for(int n=0; n<N_Components; n++)
  {
    for (int i=0; i<N_Active; i++)
    {
      Values[i+n*Length] += rhs.Values[i+n*Length];
    }
  }
  return *this;
}
