#include <NSEErrorEstimator.h>
#include "BaseCell.h"
#include "Enumerations_fe.h"
#include <FEFunction2D.h>
#include <FEVectFunct2D.h>
#include <Database.h>
#include <array>
#include "FEDatabase.h"
#include <BoundEdge.h>
#include "AuxParam2D.h"
#include "QuadratureFormulaDatabase.h"

#include <algorithm>

constexpr long loc_maxn_quadpoints_2d = MaxN_QuadPoints_2D * 40; 
// needed derivatives
const std::vector<MultiIndex2D> derivatives_u {MultiIndex2D::D00,
                                               MultiIndex2D::D01,
                                               MultiIndex2D::D10,
                                               MultiIndex2D::D02,
                                               MultiIndex2D::D20};
const std::vector<MultiIndex2D> derivatives_p {MultiIndex2D::D00,
                                               MultiIndex2D::D01,
                                               MultiIndex2D::D10};

/**
 * Data structure holding relevant data of the edges.
 *
 * XEdge1D, YEdge1D: coordinates of the edges
 * xi1D, eta1D: arrays holding the line quadrature points on reference cell of the form arr[BaseFunction][edgeIdx][quadraturePoint] = coordinate
 * xietaval, xideriv, etaderiv: arrays holding values of the derivatives at the quadrature points, same structure as xi1D and eta1D
 * AbsDetjk1D: determinant of the affine mapping for each edge of the form arr[edgeIdx][quadPoint] = value
 * xyval_ref1D: values in reference cell
 * xderiv_ref1D, yderiv_ref1D: derivative values in reference cell
 * xyval_1D: values in original cell
 * xderiv_1D, yderiv_1D: values of derivatives in original cell
 * 
 * note: it seems much of this data could as well be taken from the FEDatabase.
 */
template <int d>
struct NSEErrorEstimator<d>::EdgeData
{
  public:
    // pointer to backing array for xietaval_ref1D, xideriv_ref1D, etaderiv_ref1D
    //std::shared_ptr<double> xieta_ref1D_data {};
    std::vector<double> xieta_ref1D_data;
    // pointer to backing array for xi1D, eta1D
    //std::shared_ptr<double> xi_eta_1D_data {};
    std::vector<double> xi_eta_1D_data;
    // pointer to backing array for gradient and value of function on edge
    //std::shared_ptr<double> xyval_xyderiv_1D_data {};
    std::vector<double> xyval_xyderiv_1D_data;
    // pointer to backing array for gradient and value of function on reference edge
    //std::shared_ptr<double> xyval_xyderiv_ref_1D_data {};
    std::vector<double> xyval_xyderiv_ref_1D_data;

    // coordinates of the edges
    std::vector<std::vector<double>> XEdge1D, YEdge1D;

    // determinant of the affine mapping for each edge
    std::vector<std::vector<double>> AbsDetjk1D{4, std::vector<double>(loc_maxn_quadpoints_2d)};
    // values and derivative values in reference cell
    std::vector<std::array<double *, MaxN_QuadPoints_1D>> xyval_ref1D{4};
    std::vector<std::array<double *, MaxN_QuadPoints_1D>> xderiv_ref1D{4};
    std::vector<std::array<double *, MaxN_QuadPoints_1D>> yderiv_ref1D{4};
    // values and derivative values in original cell
    std::vector<double *> xyval_1D{4};
    std::vector<double *> xderiv_1D{4};
    std::vector<double *> yderiv_1D{4};

    // mapping (base_function, edge, quadPoint) -> value
    std::vector<std::array<double *, 4>> xi1D{N_BaseFuncts2D}, eta1D{N_BaseFuncts2D};

    // mapping (base_function2D, edge, quadPoint, BaseFunction) -> value
    std::vector<std::array<std::array<double *, MaxN_QuadPoints_1D>, 4>> xietaval_ref1D{N_BaseFuncts2D};
    std::vector<std::array<std::array<double *, MaxN_QuadPoints_1D>, 4>> xideriv_ref1D{N_BaseFuncts2D};
    std::vector<std::array<std::array<double *, MaxN_QuadPoints_1D>, 4>> etaderiv_ref1D{N_BaseFuncts2D};

    const double *weights1D;
    const double *zeta;
    unsigned int n_points_1D;
    // maximal number of base functions 2d
    unsigned int max_n_base_funct_2d;

    EdgeData(unsigned int max_n_base_functions_2d)
    {
      this->max_n_base_funct_2d = max_n_base_functions_2d;

      // initialize structures holding function values and derivatives on reference edge
      xyval_xyderiv_ref_1D_data.resize(3 * 4 * MaxN_QuadPoints_1D * max_n_base_funct_2d);
      {
        //double *ptr = xyval_xyderiv_ref_1D_data.get();
        auto *ptr = xyval_xyderiv_ref_1D_data.data();
        xyval_ref1D[0][0] = &ptr[0];
        xderiv_ref1D[0][0] = &ptr[4 * MaxN_QuadPoints_1D * max_n_base_funct_2d];
        yderiv_ref1D[0][0] = &ptr[2 * 4 * MaxN_QuadPoints_1D * max_n_base_funct_2d];
        for (size_t ii = 0; ii < 4; ii++)
        {
          for (size_t jj = 0; jj < MaxN_QuadPoints_1D; jj++)
          {
            xyval_ref1D[ii][jj] = xyval_ref1D[0][0] + (ii * MaxN_QuadPoints_1D * max_n_base_funct_2d + jj * max_n_base_funct_2d);
            xderiv_ref1D[ii][jj] = xderiv_ref1D[0][0] + (ii * MaxN_QuadPoints_1D * max_n_base_funct_2d + jj * max_n_base_funct_2d);
            yderiv_ref1D[ii][jj] = yderiv_ref1D[0][0] + (ii * MaxN_QuadPoints_1D * max_n_base_funct_2d + jj * max_n_base_funct_2d);
          }
        }
      }

      // initialize structures holding (xi, eta) values on edge
      xi_eta_1D_data.resize(N_BaseFuncts2D * 4 * MaxN_QuadPoints_1D * 2);
      {
        // back xi1D and eta1D by xi_eta_1D_data
        //double *ptr = xi_eta_1D_data.get();
        auto *ptr = xi_eta_1D_data.data();
        xi1D[0][0] = &ptr[0];
        eta1D[0][0] = &ptr[N_BaseFuncts2D * 4 * MaxN_QuadPoints_1D];
        for (auto ii = 0; ii < N_BaseFuncts2D; ii++)
        {
          for (auto jj = 0; jj < 4; jj++)
          {
            xi1D[ii][jj] = xi1D[0][0] + (ii * 4 * MaxN_QuadPoints_1D + jj * MaxN_QuadPoints_1D);
            eta1D[ii][jj] = eta1D[0][0] + (ii * 4 * MaxN_QuadPoints_1D + jj * MaxN_QuadPoints_1D);
          }
        }
      }
      // initialize structures holding values and derivatives on reference edge
      xieta_ref1D_data.resize(N_BaseFuncts2D * 4 * MaxN_QuadPoints_1D * max_n_base_funct_2d * 3);
      {
        // back xietaval_ref1D, xideriv_ref1D, etaderiv_ref1D by xieta_ref1D_data
        //double *ptr = xieta_ref1D_data.get();
        auto *ptr = xieta_ref1D_data.data();
        xietaval_ref1D[0][0][0] = &ptr[0];
        xideriv_ref1D[0][0][0] = &ptr[N_BaseFuncts2D * 4 * MaxN_QuadPoints_1D * max_n_base_funct_2d];
        etaderiv_ref1D[0][0][0] = &ptr[N_BaseFuncts2D * 4 * MaxN_QuadPoints_1D * max_n_base_funct_2d * 2];
        for (size_t ii = 0; ii < N_BaseFuncts2D; ii++)
        {
          for (size_t jj = 0; jj < 4; jj++)
          {
            for (size_t kk = 0; kk < MaxN_QuadPoints_1D; kk++)
            {
              xietaval_ref1D[ii][jj][kk] = xietaval_ref1D[0][0][0] + (kk * max_n_base_funct_2d + jj * MaxN_QuadPoints_1D * max_n_base_funct_2d + ii * 4 * MaxN_QuadPoints_1D * max_n_base_funct_2d);
              xideriv_ref1D[ii][jj][kk] = xideriv_ref1D[0][0][0] + (kk * max_n_base_funct_2d + jj * MaxN_QuadPoints_1D * max_n_base_funct_2d + ii * 4 * MaxN_QuadPoints_1D * max_n_base_funct_2d);
              etaderiv_ref1D[ii][jj][kk] = etaderiv_ref1D[0][0][0] + (kk * max_n_base_funct_2d + jj * MaxN_QuadPoints_1D * max_n_base_funct_2d + ii * 4 * MaxN_QuadPoints_1D * max_n_base_funct_2d);
            }
          }
        }
      }
    }

    void setQuadratureFormula(const TQuadFormula &quadFormula)
    {
      if(quadFormula.get_dimension() != 1)
        ErrThrow("only 1D quadrature formulas can be set here");
      quadFormula.GetFormulaData((int &) n_points_1D, weights1D, zeta);

      // initialize structures holding function values and derivatives on edge
      //xyval_xyderiv_1D_data.reset(new double[3 * 4 * 3 * n_points_1D]);
      xyval_xyderiv_1D_data.resize(3 * 4 * 3 * n_points_1D);
      {
        //double *ptr = xyval_xyderiv_1D_data.get();
        auto *ptr = xyval_xyderiv_1D_data.data();
        for (size_t i = 0; i < 4; i++)
        {
          xyval_1D[i] = &ptr[i * 3 * n_points_1D];
          xderiv_1D[i] = &ptr[1 * 4 * 3 * n_points_1D + i * 3 * n_points_1D];
          yderiv_1D[i] = &ptr[2 * 4 * 3 * n_points_1D + i * 3 * n_points_1D];
        }
      }
      XEdge1D = std::vector<std::vector<double>>(4, std::vector<double>(n_points_1D));
      YEdge1D = std::vector<std::vector<double>>(4, std::vector<double>(n_points_1D));
    }
};

template <int d>
struct NSEErrorEstimator<d>::EdgeRefData
{
  protected:
    unsigned int refDataQuadPoints1D;
    std::vector<double> neigh_data{};
    std::vector<double> ref_data{};
    std::vector<double> refNeigh1D_data{};

  public:
    unsigned int max_n_base_functions_2d;
    std::vector<double> FEFunctValuesNeigh;
    double *xi1DNeigh, *eta1DNeigh, *X1DNeigh, *Y1DNeigh, *X1DCell, *Y1DCell;
    double *xderiv_Neigh1D, *yderiv_Neigh1D, *xyval_Neigh1D, *xderiv_Cell1D, *yderiv_Cell1D, *xyval_Cell1D;

    std::array<std::array<double *, MaxN_QuadPoints_1D>, MaxN_BaseFunctions2D> xietaval_refNeigh1D;
    std::array<std::array<double *, MaxN_QuadPoints_1D>, MaxN_BaseFunctions2D> xideriv_refNeigh1D;
    std::array<std::array<double *, MaxN_QuadPoints_1D>, MaxN_BaseFunctions2D> etaderiv_refNeigh1D;
    std::vector<double> absdet1D;

    std::vector<double *> xyval_refNeigh1D{}, xderiv_refNeigh1D{}, yderiv_refNeigh1D{};

    EdgeRefData(unsigned int max_n_base_functions_2d) : refDataQuadPoints1D(0), max_n_base_functions_2d(max_n_base_functions_2d)
    {
      FEFunctValuesNeigh = std::vector<double>(3 * max_n_base_functions_2d);
      absdet1D = std::vector<double>(loc_maxn_quadpoints_2d);
      unsigned int k = MaxN_BaseFunctions2D * MaxN_QuadPoints_1D * max_n_base_functions_2d;
      //refNeigh1D_data.reset(new double[3*k], std::default_delete<double[]>());
      refNeigh1D_data.resize(3 * k);
      //double *refNeigh1D_ptr = refNeigh1D_data.get();
      double *refNeigh1D_ptr = &refNeigh1D_data[0];
      xietaval_refNeigh1D[0][0] = &refNeigh1D_ptr[0];
      xideriv_refNeigh1D[0][0] = &refNeigh1D_ptr[0] + k;
      etaderiv_refNeigh1D[0][0] = &refNeigh1D_ptr[0] + 2 * k;
      for (unsigned int i = 0; i < MaxN_BaseFunctions2D; i++)
      {
        unsigned int n = i * MaxN_BaseFunctions2D * max_n_base_functions_2d;
        for (unsigned int j = 0; j < MaxN_QuadPoints_1D; j++)
        {
          unsigned int l = n + j * max_n_base_functions_2d;
          xietaval_refNeigh1D[i][j] = xietaval_refNeigh1D[0][0] + l;
          xideriv_refNeigh1D[i][j] = xideriv_refNeigh1D[0][0] + l;
          etaderiv_refNeigh1D[i][j] = etaderiv_refNeigh1D[0][0] + l;
        }
      }
    }

    void updateNQuadPoints1D(unsigned int newNQuadPoints1D)
    {
      if (this->refDataQuadPoints1D != newNQuadPoints1D)
      {
        refDataQuadPoints1D = newNQuadPoints1D;
        // assign neigh_data
        {
          //neigh_data.reset(new double[24 * newNQuadPoints1D], std::default_delete<double[]>());
          neigh_data.resize(24 * newNQuadPoints1D);
          //double *ptr = neigh_data.get();
          auto *ptr = &neigh_data[0];
          xi1DNeigh = &ptr[0];
          eta1DNeigh = &ptr[newNQuadPoints1D];
          X1DNeigh = &ptr[2 * newNQuadPoints1D];
          Y1DNeigh = &ptr[3 * newNQuadPoints1D];
          X1DCell = &ptr[4 * newNQuadPoints1D];
          Y1DCell = &ptr[5 * newNQuadPoints1D];

          xderiv_Neigh1D = &ptr[6 * newNQuadPoints1D];
          yderiv_Neigh1D = &ptr[9 * newNQuadPoints1D];
          xyval_Neigh1D = &ptr[12 * newNQuadPoints1D];
          xderiv_Cell1D = &ptr[15 * newNQuadPoints1D];
          yderiv_Cell1D = &ptr[18 * newNQuadPoints1D];
          xyval_Cell1D = &ptr[21 * newNQuadPoints1D];
        }
        // assign ref_data
        {
          //ref_data.reset(new double[3 * newNQuadPoints1D * MaxN_BaseFunctions2D], std::default_delete<double[]>());
          ref_data.resize(3 * newNQuadPoints1D * MaxN_BaseFunctions2D);
          //double *ptr = ref_data.get();
          auto *ptr = &ref_data[0];
          //ref_data.resize(3 * newNQuadPoints1D * MaxN_BaseFunctions2D);
          xyval_refNeigh1D.resize(newNQuadPoints1D);
          xderiv_refNeigh1D.resize(newNQuadPoints1D);
          yderiv_refNeigh1D.resize(newNQuadPoints1D);
          for (size_t i = 0; i < newNQuadPoints1D; i++)
          {
            xyval_refNeigh1D[i] = &ptr[i * MaxN_BaseFunctions2D];
            xderiv_refNeigh1D[i] = &ptr[newNQuadPoints1D * MaxN_BaseFunctions2D + i * MaxN_BaseFunctions2D];
            yderiv_refNeigh1D[i] = &ptr[2 * newNQuadPoints1D * MaxN_BaseFunctions2D + i * MaxN_BaseFunctions2D];
          }
        }
      }
    }
};

//
// Helper functions
//
  std::vector<double> getWeights(const double hK, const double /*delta*/,
                                 const double *coeff)
  {
    std::vector<double> alpha(N_NSE2D_ESTIMATOR_TYPES);
    alpha[0] = hK * hK;              // weight for H^1 estimator
    alpha[1] = hK * hK * hK * hK;    // weight for L^2 estimator
    alpha[2] = 1;                    // weight for energy norm estimator
    if (hK * hK / coeff[0] < 1)
      alpha[2] = hK * hK / coeff[0]; // update weight for energy norm estimator
    return alpha;
  }

template<int d>
NSEErrorEstimator<d>::NSEErrorEstimator(const ParameterDatabase& param_db)
: ErrorEstimator<d>(param_db)
{
  if((int)ErrorEstimator<d>::db["estimator_type"] < 0 || 
    (int)ErrorEstimator<d>::db["estimator_type"] >= N_NSE2D_ESTIMATOR_TYPES)
    ErrThrow("unknown CDErrorEstimator type ", 
             (int)ErrorEstimator<d>::db["estimator_type"], 
             ". Please choose a value between 0 and ", N_NSE2D_ESTIMATOR_TYPES);
  estimatorType = (NSE2DErrorEstimatorType)(int)ErrorEstimator<d>::db["estimator_type"];
  conform_grid = true;
  is_nse = true;
  if(param_db.contains("problem_type"))
  {
    is_nse = (param_db["problem_type"].is(5) || param_db["problem_type"].is(6));
  }
  
//   if(estimatorType != NSE2DErrorEstimatorType::gradient_indicator)
//   {
//     // the code to use other estimator types is included in this file, but it is
//     // not working. Please fix this.
//     ErrThrow("estimator types other than the simple gradient indicator do not "
//              "work currently. This needs to be fixed.");
//   }
}


template <int d>
void NSEErrorEstimator<d>::estimate(const Example_NSE &ex,
                                    const FEVectFunct &fe_function_u,
                                    const FEFunction &fe_function_p,
                                    TAuxParam2D &Aux)
{
  //taking from the base class ErrorEstimator
  this->eta_K.clear();
  //Initializing the global error to zero
  this->estimated_global_error = 0;

  // this call is obligatory such that
  // the refinement strategy knows the current collection
  this->currentCollection = fe_function_u.GetFESpace()->GetCollection();

  this->eta_K.resize(this->currentCollection->GetN_Cells());
  this->maximal_local_error = 0.;

  // spaces u and p
  constexpr int n_fespaces = 2;

  FESpace const * fe_spaces[2] = {fe_function_u.GetFESpace2D().get(),
                                  fe_function_p.GetFESpace2D().get()};

  /**
  * Maximal number of base functions for u and p, respectively
  */
  const unsigned int max_n_base_functions_u = fe_spaces[0]->get_max_n_local_dof();
  const unsigned int max_n_base_functions_p = fe_spaces[1]->get_max_n_local_dof();
  const unsigned int max_n_base_functions = max_n_base_functions_u + max_n_base_functions_p;

  int N_Points, N_Edges, N_;
  std::vector<int> Used(N_FEs2D);

  double *Param[loc_maxn_quadpoints_2d];
  double *Derivatives[3 * loc_maxn_quadpoints_2d];
  double *AuxArray[loc_maxn_quadpoints_2d];
  std::vector<double> FEFunctValues(3 * max_n_base_functions);
  double estimated_global_errors[N_NSE2D_ESTIMATOR_TYPES];
  double estimated_local_errors[N_NSE2D_ESTIMATOR_TYPES];

  // ########################################################################
  // store information in local arrays
  // ########################################################################
  // note that char is used instead of bool to prevent specialization of the vector
  std::vector<char> NeedsSecondDer(n_fespaces);

  auto used_fe_u = fe_spaces[0]->GetUsedElements();
  auto used_fe_p = fe_spaces[1]->GetUsedElements();
  std::vector<FE_type> UsedElements;
  UsedElements.reserve(used_fe_u.size() + used_fe_p.size());//preallocate memory
  UsedElements.insert(UsedElements.end(), used_fe_u.begin(), used_fe_u.end());
  UsedElements.insert(UsedElements.end(), used_fe_p.begin(), used_fe_p.end());
  std::sort(UsedElements.begin(), UsedElements.end());
  auto last = std::unique(UsedElements.begin(), UsedElements.end());
  UsedElements.erase(last, UsedElements.end());
  int NUsedElements = UsedElements.size();

  // ########################################################################
  // prepare error estimates
  // ########################################################################

  EdgeData edgeData(max_n_base_functions);
  EdgeRefData edgeRefData(max_n_base_functions);

  auto *values_u = fe_function_u.GetValues();                      // values of fe function
  auto *values_p = fe_function_p.GetValues();

  this->currentCollection->mark_all_cells();

  // initialize some quantities
  for (size_t i = 0; i < N_NSE2D_ESTIMATOR_TYPES; i++)
  {
    estimated_global_errors[i] = 0.0;
  }
  double max_loc_err = 0;

  // ########################################################################
  // calculate values of base functions and derivatives on ref element
  // ########################################################################

  {
    // fe for velocity
    auto CurrentElement = (fe_function_u.GetFESpace2D()->GetUsedElements())[0];
    const FiniteElement element{CurrentElement};
    // get 1d quad fromula
    auto l = element.GetBaseFunct()->GetPolynomialDegree();
    auto LineQuadFormula = QuadratureFormulaDatabase::qf_from_degree(
        2 * l, BFRefElements::BFUnitLine);
    edgeData.setQuadratureFormula(*LineQuadFormula);
  }
  for (int i = 0; i < NUsedElements; i++)
  {
    auto CurrentElement = UsedElements[i];
    const FiniteElement element{CurrentElement};
    const int N_Points1D = edgeData.n_points_1D;
    auto bf = element.GetBaseFunct();
    auto BaseFunct = bf->GetID();
    auto bf2Drefelements = element.GetBaseFunct()->GetRefElement();

    switch (bf2Drefelements)                      // compute coordinates of line quadrature
    {                                            // points in reference cell
      // quadrilateral cell
      case BFRefElements::BFUnitSquare :                          // edge 0
        for (int j = 0; j < N_Points1D; j++)                 // for all quadrature points
        {
          edgeData.xi1D[BaseFunct][0][j] = edgeData.zeta[j];
          edgeData.eta1D[BaseFunct][0][j] = -1;
          bf->GetDerivatives(MultiIndex2D::D00, edgeData.zeta[j], -1, edgeData.xietaval_ref1D[BaseFunct][0][j]);
          bf->GetDerivatives(MultiIndex2D::D10, edgeData.zeta[j], -1, edgeData.xideriv_ref1D[BaseFunct][0][j]);
          bf->GetDerivatives(MultiIndex2D::D01, edgeData.zeta[j], -1, edgeData.etaderiv_ref1D[BaseFunct][0][j]);
        }                                        // edge 1
        for (int j = 0; j < N_Points1D; j++)               // for all quadrature points
        {
          edgeData.xi1D[BaseFunct][1][j] = 1;
          edgeData.eta1D[BaseFunct][1][j] = edgeData.zeta[j];
          bf->GetDerivatives(MultiIndex2D::D00, 1, edgeData.zeta[j], edgeData.xietaval_ref1D[BaseFunct][1][j]);
          bf->GetDerivatives(MultiIndex2D::D10, 1, edgeData.zeta[j], edgeData.xideriv_ref1D[BaseFunct][1][j]);
          bf->GetDerivatives(MultiIndex2D::D01, 1, edgeData.zeta[j], edgeData.etaderiv_ref1D[BaseFunct][1][j]);
        }                                        // edge 2
        for (int j = 0; j < N_Points1D; j++)                 // for all quadrature points
        {
          edgeData.xi1D[BaseFunct][2][j] = -edgeData.zeta[j];
          edgeData.eta1D[BaseFunct][2][j] = 1;
          bf->GetDerivatives(MultiIndex2D::D00, -edgeData.zeta[j], 1, edgeData.xietaval_ref1D[BaseFunct][2][j]);
          bf->GetDerivatives(MultiIndex2D::D10, -edgeData.zeta[j], 1, edgeData.xideriv_ref1D[BaseFunct][2][j]);
          bf->GetDerivatives(MultiIndex2D::D01, -edgeData.zeta[j], 1, edgeData.etaderiv_ref1D[BaseFunct][2][j]);
        }                                         // edge 3
        for (int j = 0; j < N_Points1D; j++)                  // for all quadrature points
        {
          edgeData.xi1D[BaseFunct][3][j] = -1;
          edgeData.eta1D[BaseFunct][3][j] = -edgeData.zeta[j];
          bf->GetDerivatives(MultiIndex2D::D00, -1, -edgeData.zeta[j], edgeData.xietaval_ref1D[BaseFunct][3][j]);
          bf->GetDerivatives(MultiIndex2D::D10, -1, -edgeData.zeta[j], edgeData.xideriv_ref1D[BaseFunct][3][j]);
          bf->GetDerivatives(MultiIndex2D::D01, -1, -edgeData.zeta[j], edgeData.etaderiv_ref1D[BaseFunct][3][j]);
        }
        break;

      case BFRefElements::BFUnitTriangle :                        // triangular cell
        for (int j = 0; j < N_Points1D; j++)                 // for all quadrature points
        {
          edgeData.xi1D[BaseFunct][0][j] = (edgeData.zeta[j] + 1) / 2;
          edgeData.eta1D[BaseFunct][0][j] = 0;
          bf->GetDerivatives(MultiIndex2D::D00, (edgeData.zeta[j] + 1) / 2, 0, edgeData.xietaval_ref1D[BaseFunct][0][j]);
          bf->GetDerivatives(MultiIndex2D::D10, (edgeData.zeta[j] + 1) / 2, 0, edgeData.xideriv_ref1D[BaseFunct][0][j]);
          bf->GetDerivatives(MultiIndex2D::D01, (edgeData.zeta[j] + 1) / 2, 0, edgeData.etaderiv_ref1D[BaseFunct][0][j]);
        }                                       // edge 1
        for (int j = 0; j < N_Points1D; j++)                // for all quadrature points
        {
          edgeData.xi1D[BaseFunct][1][j] = (-edgeData.zeta[j] + 1) / 2;
          edgeData.eta1D[BaseFunct][1][j] = (edgeData.zeta[j] + 1) / 2;
          bf->GetDerivatives(MultiIndex2D::D00, (-edgeData.zeta[j] + 1) / 2, (edgeData.zeta[j] + 1) / 2, edgeData.xietaval_ref1D[BaseFunct][1][j]);
          bf->GetDerivatives(MultiIndex2D::D10, (-edgeData.zeta[j] + 1) / 2, (edgeData.zeta[j] + 1) / 2, edgeData.xideriv_ref1D[BaseFunct][1][j]);
          bf->GetDerivatives(MultiIndex2D::D01, (-edgeData.zeta[j] + 1) / 2, (edgeData.zeta[j] + 1) / 2, edgeData.etaderiv_ref1D[BaseFunct][1][j]);
        }                                       // edge 2
        for (int j = 0; j < N_Points1D; j++)                // for all quadrature points
        {
          edgeData.xi1D[BaseFunct][2][j] = 0;
          edgeData.eta1D[BaseFunct][2][j] = (-edgeData.zeta[j] + 1) / 2;
          bf->GetDerivatives(MultiIndex2D::D00, 0, (-edgeData.zeta[j] + 1) / 2, edgeData.xietaval_ref1D[BaseFunct][2][j]);
          bf->GetDerivatives(MultiIndex2D::D10, 0, (-edgeData.zeta[j] + 1) / 2, edgeData.xideriv_ref1D[BaseFunct][2][j]);
          bf->GetDerivatives(MultiIndex2D::D01, 0, (-edgeData.zeta[j] + 1) / 2, edgeData.etaderiv_ref1D[BaseFunct][2][j]);
        }
        break;
        
      default:
        ;
    }
  }

  // aux data
  // get number of parameters of equation
  auto N_Parameters = Aux.GetN_Parameters();
  std::vector<double> ParamData(loc_maxn_quadpoints_2d * N_Parameters + 3 * loc_maxn_quadpoints_2d * derivatives_u.size() + loc_maxn_quadpoints_2d * 20);
  {

    auto *ParamDataPtr = &ParamData[0];
    // set pointers
    for (auto j = 0; j < loc_maxn_quadpoints_2d; j++)
    {
      Param[j] = &ParamDataPtr[0] + j * N_Parameters;
      AuxArray[j] = &ParamDataPtr[loc_maxn_quadpoints_2d * N_Parameters + 3 * loc_maxn_quadpoints_2d * derivatives_u.size()] + j * 20;
    }
    for (auto j = 0; j < 3 * loc_maxn_quadpoints_2d; j++)
    {
      Derivatives[j] = &ParamDataPtr[loc_maxn_quadpoints_2d * N_Parameters] + j * derivatives_u.size();
    }
  }
  std::vector< const FiniteElement*> LocalUsedElements(n_fespaces, nullptr);
  std::vector<const BaseFunctions*> LocBF(n_fespaces);
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);
  
  // ########################################################################
  // loop over all cells
  // ########################################################################
  // for all cells in the collection
  for(int cellIdx = 0, n_cells = this->currentCollection->GetN_Cells();
      cellIdx < n_cells; cellIdx++)
  {
    // next cell
    TBaseCell *cell = this->currentCollection->GetCell(cellIdx);

    const TFESpace2D *fe_space_u = fe_function_u.GetFESpace2D().get();
    const TFESpace2D *fe_space_p = fe_function_p.GetFESpace2D().get();

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    for (int j = 0; j < n_fespaces; j++)
    {
      // j == 0: velocity
      // j == 1: pressure
      auto& element = fe_spaces[j]->get_fe(cellIdx);
      LocalUsedElements[j] = &element;
      LocBF[j] = element.GetBaseFunct();
      // velocity needs 2nd derivative
      NeedsSecondDer[j] = (j == 0);
    }

    // ####################################################################
    // calculate values on original element
    // ####################################################################

    // get reference transformation
    ReferenceTransformation_type RefTrans = FEDatabase::GetOrig(
      LocalUsedElements, this->currentCollection, cell,
      (bool *) NeedsSecondDer.data(), qf_ref, qf_orig);
    const double *xi, *eta, *weights, *X, *Y;
    qf_ref.GetFormulaData(N_Points, weights, xi, eta);
    qf_orig.GetFormulaData(N_Points, weights, X, Y);
    // get parameters of equ.
    if (N_Parameters > 0)
    {
      Aux.GetParameters(qf_orig, cellIdx, Param);
    }

    // velocity first
    {
      // calculate all needed derivatives of this FE function
      auto& element = fe_space_u->get_fe(cellIdx);
      auto bf = element.GetBaseFunct();
      N_ = element.GetN_DOF();
      // dof of current mesh cell
      auto DOF_u = fe_space_u->GetGlobalDOF(cellIdx);
      for (int l = 0; l < N_; l++)
      {
        // u fe values of dofs
        FEFunctValues[l] = values_u[DOF_u[l]];
        // v fe values of dofs
        FEFunctValues[l + max_n_base_functions] = values_u[DOF_u[l] + fe_function_u.GetLength()];
      }

      // compute values for all derivatives
      // in all quadrature points
      // in original mesh cell
      {
        double **OrigFEValues;
        double value[2];
        // for all derivatives
        for (size_t k = 0; k < derivatives_u.size(); k++)
        {
          // get values in original cell
          OrigFEValues = FEDatabase::GetOrigElementValues(*bf,
                                                          derivatives_u[k]);
          // for all quadrature points
          for (int j = 0; j < N_Points; j++)
          {
            // value in original cell
            double *Orig = OrigFEValues[j];
            value[0] = value[1] = 0;
            // for all basis functions
            for (int l = 0; l < N_; l++)
            {
              // accumulate u value of derivative in point j
              value[0] += FEFunctValues[l] * Orig[l];
              // accumulate v value of derivative in point j
              value[1] += FEFunctValues[l + max_n_base_functions] * Orig[l];
            }
            // k-th u derivative
            Derivatives[j][k] = value[0];
            // k-th v derivative
            Derivatives[j + loc_maxn_quadpoints_2d][k] = value[1];
          }
        }
      }

      // get coefficients of pde
      if (ex.get_coeffs())
      {
        (ex.get_coeffs())(N_Points, X, Y, Param, AuxArray);
      }
    }
    // prepare 1D quadrature formula
    // pressure first
    // finite element on cell
    auto& element_p = fe_space_p->get_fe(cellIdx);
    // basis functions
    auto bf_p = element_p.GetBaseFunct();
    auto BaseFunctP = bf_p->GetID();
    // # basis functions
    int N_P = element_p.GetN_DOF();
    {

      auto DOF_p = fe_space_p->GetGlobalDOF(cellIdx);    // dof of current mesh cell
      //cout << "cell no " << i << " begin " << *DOF_p <<
      // " " << *GlobalNumbersP << " " << BeginIndexP[i] << endl;
      for (int l = 0; l < N_P; l++)
      {
        FEFunctValues[l + 2 * max_n_base_functions] = values_p[DOF_p[l]];    // p fe values of dofs
        //cout <<  FEFunctValues[l+2*max_n_base_functions] << endl;
      }

      {
        double value;
        double **OrigFEValues;
        for (size_t k = 0; k < derivatives_p.size(); k++)            // for all derivatives
        {                                       // get values in original cell
          OrigFEValues = FEDatabase::GetOrigElementValues(*bf_p,
                                                          derivatives_p[k]);
          for (int j = 0; j < N_Points; j++)               // for all quadrature points
          {
            double *Orig = OrigFEValues[j];             // value in original cell
            value = 0;
            for (int l = 0; l < N_P; l++)                   // for all basis functions
            {
              value += FEFunctValues[l + 2 * max_n_base_functions] * Orig[l]; // accumulate p value of derivative in point j
            }
            Derivatives[j + 2 * loc_maxn_quadpoints_2d][k] = value;// for k-th v derivative
          }                                     // endfor j
        } // endfor k
      }
    }

    unsigned int n_quadrature_points_1d;
    {
      unsigned int l = bf_p->GetPolynomialDegree();
      auto qf1D = QuadratureFormulaDatabase::qf_from_degree(
          2 * l, BFRefElements::BFUnitLine);
      {
        int n_points1d_val;
        qf1D->GetFormulaData(n_points1d_val, edgeData.weights1D, edgeData.zeta);
        n_quadrature_points_1d = (unsigned int) n_points1d_val;
      }
    }

    // # edges
    N_Edges = cell->GetN_Edges();
    // pressure second
    {
      // loop over all edges of cell
      for (int edgeIndex = 0; edgeIndex < N_Edges; edgeIndex++)
      {
        // get original coordinates of edge quad. points
        FEDatabase::GetOrigFromRef(RefTrans, (int) n_quadrature_points_1d,
                                   edgeData.xi1D[BaseFunctP][edgeIndex],
                                   edgeData.eta1D[BaseFunctP][edgeIndex],
                                   edgeData.XEdge1D[edgeIndex].data(),
                                   edgeData.YEdge1D[edgeIndex].data());
        // get values and derivatives in original cell
        for (unsigned int k = 0; k < n_quadrature_points_1d; k++)
        {
          FEDatabase::GetOrigValues(
            RefTrans, edgeData.xi1D[BaseFunctP][edgeIndex][k],
            edgeData.eta1D[BaseFunctP][edgeIndex][k],
            bf_p, this->currentCollection, (TGridCell *) cell,
            edgeData.xietaval_ref1D[BaseFunctP][edgeIndex][k],
            edgeData.xideriv_ref1D[BaseFunctP][edgeIndex][k],
            edgeData.etaderiv_ref1D[BaseFunctP][edgeIndex][k],
            edgeData.xyval_ref1D[edgeIndex][k],
            edgeData.xderiv_ref1D[edgeIndex][k],
            edgeData.yderiv_ref1D[edgeIndex][k]);
        }

        {
          double value;
          size_t m;
          // for all quadrature points
          for (size_t k = 0; k < n_quadrature_points_1d; k++)
          {
            value = 0;
            // for all basis functions
            for (int l = 0; l < N_P; l++)
            {
              m = (size_t) (l + 2 * max_n_base_functions);
              value += FEFunctValues[m] * edgeData.xyval_ref1D[edgeIndex][k][l]; // accumulate value of derivative
            } // endfor l
            m = k + 2 * n_quadrature_points_1d;
            edgeData.xyval_1D[edgeIndex][m] = value; // for k-th
          } // endfor k
        }
      }     // endfor j
    }
    // velocity second
    {
      auto& element = fe_space_u->get_fe(cellIdx);
      auto bf = element.GetBaseFunct();
      auto BaseFunct = bf->GetID();
      unsigned int m;
      double val[6];
      for (int edgeIndex = 0; edgeIndex < N_Edges; edgeIndex++)                  // loop over all edges of cell
      {                                     // get original coordinates of edge quad. points
        FEDatabase::GetOrigFromRef(RefTrans, (int) n_quadrature_points_1d,
                                   edgeData.xi1D[BaseFunct][edgeIndex],
                                   edgeData.eta1D[BaseFunct][edgeIndex],
                                   edgeData.XEdge1D[edgeIndex].data(), 
                                   edgeData.YEdge1D[edgeIndex].data());
        for (size_t k = 0; k < n_quadrature_points_1d; k++)  // get values and derivatives in original cell
        {
          FEDatabase::GetOrigValues(
            RefTrans, edgeData.xi1D[BaseFunct][edgeIndex][k], 
            edgeData.eta1D[BaseFunct][edgeIndex][k],
            bf, this->currentCollection, (TGridCell *) cell,
            edgeData.xietaval_ref1D[BaseFunct][edgeIndex][k], 
            edgeData.xideriv_ref1D[BaseFunct][edgeIndex][k],
            edgeData.etaderiv_ref1D[BaseFunct][edgeIndex][k], 
            edgeData.xyval_ref1D[edgeIndex][k],
            edgeData.xderiv_ref1D[edgeIndex][k], 
            edgeData.yderiv_ref1D[edgeIndex][k]);
        }

        for (size_t k = 0; k < n_quadrature_points_1d; k++)     // for all quadrature points
        {
          val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = 0;

          for (int l = 0; l < N_; l++)       // for all basis functions
          {
            m = (int) (l + max_n_base_functions);
            val[0] += FEFunctValues[l] * edgeData.xyval_ref1D[edgeIndex][k][l]; // accumulate value of derivative
            val[1] += FEFunctValues[l] * edgeData.xderiv_ref1D[edgeIndex][k][l]; // accumulate value of derivative
            val[2] += FEFunctValues[l] * edgeData.yderiv_ref1D[edgeIndex][k][l]; // accumulate value of derivative
            val[3] += FEFunctValues[m] * edgeData.xyval_ref1D[edgeIndex][k][l]; // accumulate value of derivative
            val[4] += FEFunctValues[m] * edgeData.xderiv_ref1D[edgeIndex][k][l]; // accumulate value of derivative
            val[5] += FEFunctValues[m] * edgeData.yderiv_ref1D[edgeIndex][k][l]; // accumulate value of derivative
          } // endfor l
          m = (int) (k + n_quadrature_points_1d);
          edgeData.xyval_1D[edgeIndex][k] = val[0]; // for k-th
          edgeData.xderiv_1D[edgeIndex][k] = val[1]; // for k-th
          edgeData.yderiv_1D[edgeIndex][k] = val[2]; // for k-th
          edgeData.xyval_1D[edgeIndex][m] = val[3]; // for k-th
          edgeData.xderiv_1D[edgeIndex][m] = val[4]; // for k-th
          edgeData.yderiv_1D[edgeIndex][m] = val[5]; // for k-th
        } // endfor k
      }     // endfor j
    }

    {
      // estimate local errors
      calculateEtaK(fe_function_u, fe_function_p, cell, (unsigned int) N_Points,
                    n_quadrature_points_1d, weights, Derivatives, AuxArray, ex,
                    edgeData, edgeRefData, &estimated_local_errors[0]);
    }


    for (size_t k = 0; k < N_NSE2D_ESTIMATOR_TYPES; k++)
    {
      // update global error estimates
      estimated_global_errors[k] += estimated_local_errors[k];
    }
    // update maximal local error estimate
    if (estimated_local_errors[int(estimatorType)] > max_loc_err)
    {
      max_loc_err = estimated_local_errors[int(estimatorType)];
    }

    this->eta_K[cellIdx] = estimated_local_errors[int(estimatorType)];
  }

  // compute global error estimates
  for (size_t i = 1; i < N_NSE2D_ESTIMATOR_TYPES; i++)
  {
    estimated_global_errors[i] = std::sqrt(estimated_global_errors[i]);
  }
  this->estimated_global_error = estimated_global_errors[0];
  // compute maximal local error estimate
  this->maximal_local_error = std::sqrt(max_loc_err);
}

template <int d>
void NSEErrorEstimator<d>::calculateEtaK(const TFEVectFunct2D &fe_function_u,
                                         const TFEFunction2D &fe_function_p,
                                         const TBaseCell *cell,
                                         unsigned int N_Points,
                                         unsigned int N_Points1D,
                                         const double *weights,
                                         double **Derivatives,
                                         double **coeffs,
                                         const Example_NSE &example,
                                         EdgeData &edgeData,
                                         EdgeRefData &edgeRefData,
                                         double *estimated_local_error)
{

  double estimated_error[N_NSE2D_ESTIMATOR_TYPES];
  double delta = TDatabase::ParamDB->FILTER_WIDTH_CONSTANT;
  int ee_verbose = 0;
  auto fe_space_u = fe_function_u.GetFESpace2D();
  auto fe_space_p = fe_function_p.GetFESpace2D();
  const double * values_u = fe_function_u.GetValues();
  const double * values_p = fe_function_p.GetValues();
  
  auto collection_u = fe_space_u->GetCollection();

  // initalize local error estimates
  for (size_t i = 0; i < N_NSE2D_ESTIMATOR_TYPES; i++)
    estimated_error[i] = 0.0;

  if (estimatorType == NSE2DErrorEstimatorType::gradient_indicator)
  {
    double w;
    double e1, e2, e3, e4;
    // for all quadrature points
    for (unsigned int i = 0; i < N_Points; i++)
    {
      w = weights[i];
      // all u derivatives in quadrature points
      auto *deriv = Derivatives[i];
      // x derivative
      e1 = deriv[0];
      // y derivative
      e2 = deriv[1];
      deriv = Derivatives[i + loc_maxn_quadpoints_2d];
      // all v derivatives in quadrature points
      e3 = deriv[0];
      e4 = deriv[1];
      estimated_error[int(estimatorType)] += w * (e1 * e1 + e2 * e2 + e3 * e3 + e4 * e4);
    }
  }
  else if (int(estimatorType) > 0)
  {

    edgeRefData.updateNQuadPoints1D(N_Points1D);

    int N_U = fe_function_u.GetLength();


    bool check_cont_u = false;
    bool check_cont_p;
    if (TDatabase::ParamDB->ANSATZ_ORDER > 0)
      check_cont_u = 1;
    if (TDatabase::ParamDB->ANSATZ_ORDER < -1)
      check_cont_u = 1;
    // this has to be set for all continuous pressure spaces
    if ((TDatabase::ParamDB->VELOCITY_SPACE > 0) && (TDatabase::ParamDB->VELOCITY_SPACE < 10))
        check_cont_p = 1;
    else
        check_cont_p = 0;

/*************************************************************************/
/*  strong residual                                                      */
/*************************************************************************/
    auto hK = cell->GetDiameter();
    double strong_residual = 0;
    // for all quadrature points
    double *coeff = coeffs[0];
    for (size_t i = 0; i < N_Points; i++)
    {
      coeff = coeffs[i];

      // all u derivatives in quadrature points
      auto *derivatives_u = Derivatives[i];
      auto *derivatives_v = Derivatives[i + loc_maxn_quadpoints_2d];
      auto *derivatives_p = Derivatives[i + 2 * loc_maxn_quadpoints_2d];
      auto w = weights[i];
      double e1, e2;
      if (is_nse)
      {
        // strong residual, u - comp
        e1 = -coeff[0] * (derivatives_u[3] + derivatives_u[4]) + derivatives_u[2] * derivatives_u[0] + derivatives_v[2] * derivatives_u[1] + derivatives_p[0] - coeff[1];
        // strong residual, v - comp
        e2 = -coeff[0] * (derivatives_v[3] + derivatives_v[4]) + derivatives_u[2] * derivatives_v[0] + derivatives_v[2] * derivatives_v[1] + derivatives_p[1] - coeff[2];
      }
      else
      {
        // strong residual, u - comp
        e1 = -coeff[0] * (derivatives_u[3] + derivatives_u[4]) + derivatives_p[0] - coeff[1];
        // strong residual, v - comp
        e2 = -coeff[0] * (derivatives_v[3] + derivatives_v[4]) + derivatives_p[1] - coeff[2];
      }
      // L^2 norm
      strong_residual += w * (e1 * e1 + e2 * e2);
    }
    {
      auto alphas = getWeights(hK, delta, coeffs[N_Points - 1]);
      for (size_t i = 1; i < 4; i++)
      {
        estimated_error[i] = alphas[i - 1] * strong_residual;
      }
    }
/*************************************************************************/
/*  compute jumps across the edges                                       */
/*************************************************************************/
    const unsigned int N_Edges = (unsigned int) cell->GetN_Edges();
    // loop over all edges of cell
    for (size_t edgeIdx = 0; edgeIdx < N_Edges; edgeIdx++)
    {
      const TJoint *joint = cell->GetJoint(int(edgeIdx));
      // boundary edge
      if ((joint->GetType() == BoundaryEdge) || (joint->GetType() == IsoBoundEdge))
      {
        double t0, t1;
        double x0, x1, y0, y1;
        TBoundEdge *bdryEdge = (TBoundEdge *) joint;
        BoundCond Cond0;
        // get boundary component
        const TBoundComp2D *BoundComp = bdryEdge->GetBoundComp();
        bdryEdge->GetParameters(t0, t1);     // parameter interval
        int comp = BoundComp->GetID();              // boundary id
        example.get_bc()[0](comp, (t0 + t1) / 2.0, Cond0); // type of boundary condition
        // at midpoint of boundary
        switch (Cond0)
        {
          // boundary is Dirichlet
          case DIRICHLET:
          {
              // no error
              break;
          }
          // boundary is Neumann
          case NEUMANN:
          {
            // coordinates at begin of parameter interval
            bdryEdge->GetXYofT(t0, x0, y0);
            // coordinates at end of parameter interval
            bdryEdge->GetXYofT(t1, x1, y1);
            // outer normal vector
            double nx = y1 - y0;
            double ny = x0 - x1;
            // length of edge
            double hE = std::sqrt(nx * nx + ny * ny);
            // normalized normal vector
            nx /= hE;
            ny /= hE;
            double jump = 0;
            // compute difference to Neumann condition
            for (size_t i = 0; i < N_Points1D; i++)
            {
              size_t m = i + N_Points1D;
              x0 = edgeData.XEdge1D[edgeIdx][i];
              y0 = edgeData.YEdge1D[edgeIdx][i];
              // coordinates at quadrature points
              bdryEdge->GetXYofT(t0, x0, y0);
              double neumann_data;
              // Neumann data
              example.boundary_data[0](comp, t0, neumann_data);
              double e1 = coeff[0] * (edgeData.xderiv_1D[edgeIdx][i] * nx + edgeData.yderiv_1D[edgeIdx][i] * ny) - neumann_data;
              // Neumann data
              example.get_bd()[1](comp, t0, neumann_data);
              double e2 = coeff[0] * (edgeData.xderiv_1D[edgeIdx][m] * nx + edgeData.yderiv_1D[edgeIdx][m] * ny) - neumann_data;
              double w = edgeData.weights1D[i] * hE / 2.0;
              // integral on the edge
              jump += w * (e1 * e1 + e2 * e2);
            }
            double beta[N_NSE2D_ESTIMATOR_TYPES];
            // weight for H^1 estimator
            beta[0] = hE;
            // weight for L^2 estimator
            beta[1] = hE * hE * hE;
            // weight for energy norm estimator
            beta[2] = hE / coeff[0];
            if (1.0 / std::sqrt(coeff[0]) < beta[2])
              beta[2] = 1.0 / std::sqrt(coeff[0]);
            for (size_t i = 1; i < 4; i++)
              estimated_error[i] += beta[i - 1] * jump;
            break;
          }
          default:
          {
            ErrThrow("Only few BC implementations done ");
          }
        }
      }
      else
      {
        // inner edge.
        const int *TmpEdVer;
        // get refinement descriptor
        const TRefDesc *refdesc = cell->GetRefDesc();
        refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVer);
        // get vertices of face j
        const TVertex *ver0 = cell->GetVertex(TmpEdVer[2 * edgeIdx]);
        const TVertex *ver1 = cell->GetVertex(TmpEdVer[2 * edgeIdx + 1]);
        // coordinates of face j
        auto cell_x0 = cell->GetVertex(TmpEdVer[2 * edgeIdx])->GetX();
        auto cell_y0 = cell->GetVertex(TmpEdVer[2 * edgeIdx])->GetY();
        auto cell_x1 = cell->GetVertex(TmpEdVer[2 * edgeIdx + 1])->GetX();
        auto cell_y1 = cell->GetVertex(TmpEdVer[2 * edgeIdx + 1])->GetY();
        double jump = 0;
        const TBaseCell *neigh = joint->GetNeighbour(cell);
        auto nx = cell_y1 - cell_y0;              // compute normal
        auto ny = cell_x0 - cell_x1;
        auto hE = std::sqrt(nx * nx + ny * ny);               // length of edge
        nx /= hE;                             // normalized normal vector
        ny /= hE;
        if (ee_verbose > 1)
        {
          cout << " A " << cell_x0 << " " << cell_y0;
          cout << " B " << cell_x1 << " " << cell_y1;
          cout << " n " << nx << " " << ny;
        }

/*************************************************************************/
/*  no neighbour, find neighbour of parent                               */
/*************************************************************************/
        if (!neigh)
        {
          const int *TmpEdVerParent, *TmpoEnlE;
          const int *TmpCE;
          int MaxLen1;
          // there is no neighbour on the same level
          //  => finer cell in 1 regularity
          const TBaseCell *parent = cell->GetParent();                      // parent cell
          refdesc = parent->GetRefDesc();
          refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVerParent);
          refdesc->GetChildEdge(TmpCE, MaxLen1);
          refdesc->GetNewEdgeOldEdge(TmpoEnlE);
          size_t l = 0;
          while (parent->GetChild((int) l) != cell) l++;           // local child number
          int parent_edge = TmpCE[l * MaxLen1 + edgeIdx];               // number of father edge
          parent_edge = TmpoEnlE[parent_edge];            // number of father edge

          const TJoint *parent_joint = parent->GetJoint(parent_edge);
          neigh = parent_joint->GetNeighbour(parent);        // neighbour to parent
          if (!neigh)
          {
            cout << "Hier sollte man aber nicht hinkommen 2 !" << endl;
          }

          int neigh_edge = 0;
          while (neigh->GetJoint(neigh_edge)->GetNeighbour(neigh) != parent) neigh_edge++;
          int part = 0;
          {
            const TVertex *ver2 = neigh->GetVertex(TmpEdVerParent[2 * neigh_edge]);          // vertices of edge
            const TVertex *ver3 = neigh->GetVertex(TmpEdVerParent[2 * neigh_edge + 1]);
            if (ver1 == ver2)
            {
              // first part of long edge
              part = -1;
            }
            else if (ver0 == ver3)
            {
              // second part of long edge
              part = 1;
            }
            else
            {
              cout << "Hier sollte man aber nicht hinkommen 4 !" << endl;
            }
          }

          int neigh_N_ = neigh->GetClipBoard();                  // number of neighbour in iterator
          if (neigh_N_ == -1)
          {
            cout << "Hier sollte man aber nicht hinkommen 3 !" << endl;
          }
          auto& element_neighbor = fe_space_u->get_fe(neigh_N_); // finite element on neighbour
          auto bf_neighbor = element_neighbor.GetBaseFunct();
          BaseFunction_type BaseFunctNeigh = bf_neighbor->GetID(); // basis functions on neighbour
          int N_Neigh = element_neighbor.GetN_DOF(); // number of basis functions

          BFRefElements bf2DrefelementsNeigh = bf_neighbor->GetRefElement();   // reference cell of neighbour

          if (conform_grid)
          {
            switch (bf2DrefelementsNeigh)                // compute coordinates of line quadrature
            {                                    // points in reference cell
              case BFRefElements::BFUnitSquare :                  // edge 0
                if (neigh_edge == 0)
                {
                  for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                  {
                    edgeRefData.xi1DNeigh[i] = -edgeData.zeta[i];
                    edgeRefData.eta1DNeigh[i] = -1;
                  }
                }
                if (neigh_edge == 1)
                {                               // edge 1
                  for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                  {
                    edgeRefData.xi1DNeigh[i] = 1;
                    edgeRefData.eta1DNeigh[i] = -edgeData.zeta[i];
                  }
                }
                if (neigh_edge == 2)
                {                               // edge 2
                  for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                  {
                    edgeRefData.xi1DNeigh[i] = edgeData.zeta[i];
                    edgeRefData.eta1DNeigh[i] = 1;
                  }
                }

                if (neigh_edge == 3)
                {                               // edge 3
                  for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                  {
                    edgeRefData.xi1DNeigh[i] = -1;
                    edgeRefData.eta1DNeigh[i] = edgeData.zeta[i];
                  }
                }
                break;

              case BFRefElements::BFUnitTriangle :
                if (neigh_edge == 0)
                {
                  for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                  {
                    edgeRefData.xi1DNeigh[i] = (-edgeData.zeta[i] + 1) / 2;
                    edgeRefData.eta1DNeigh[i] = 0;
                  }
                }
                if (neigh_edge == 1)
                {
                  for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                  {
                    edgeRefData.xi1DNeigh[i] = (edgeData.zeta[i] + 1) / 2;
                    edgeRefData.eta1DNeigh[i] = (-edgeData.zeta[i] + 1) / 2;
                  }
                }
                if (neigh_edge == 2)
                {
                  for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                  {
                    edgeRefData.xi1DNeigh[i] = 0;
                    edgeRefData.eta1DNeigh[i] = (edgeData.zeta[i] + 1) / 2;
                  }
                }
                break;
                
              default:
                ;
            }
          }
          else
          {
            // this is only for 1-regular triangulations
            // compute coordinates of line quadrature
            // points in reference cell
            switch (bf2DrefelementsNeigh)
            {
              case BFRefElements::BFUnitSquare :                  // edge 0
                if (neigh_edge == 0)
                {
                  for (size_t i = 0; i < N_Points1D; i++) // for all quadrature points
                  {
                    edgeRefData.xi1DNeigh[i] = (-edgeData.zeta[i] + part) / 2;
                    edgeRefData.eta1DNeigh[i] = -1;
                  }
                }
                if (neigh_edge == 1)
                {                               // edge 1
                  for (size_t i = 0; i < N_Points1D; i++) // for all quadrature points
                  {
                    edgeRefData.xi1DNeigh[i] = 1;
                    edgeRefData.eta1DNeigh[i] = (-edgeData.zeta[i] + part) / 2;
                  }
                }
                if (neigh_edge == 2)
                {                               // edge 2
                  for (size_t i = 0; i < N_Points1D; i++) // for all quadrature points
                  {
                    edgeRefData.xi1DNeigh[i] = (edgeData.zeta[i] - part) / 2;
                    edgeRefData.eta1DNeigh[i] = 1;
                  }
                }

                if (neigh_edge == 3)
                {                               // edge 3
                  for (size_t i = 0; i < N_Points1D; i++) // for all quadrature points
                  {
                    edgeRefData.xi1DNeigh[i] = -1;
                    edgeRefData.eta1DNeigh[i] = (edgeData.zeta[i] - part) / 2;
                  }
                }
                break;

              case BFRefElements::BFUnitTriangle :
                if (neigh_edge == 0)
                {
                  for (size_t i = 0; i < N_Points1D; i++) // for all quadrature points
                  {
                    if (part == -1) part = 0;
                    edgeRefData.xi1DNeigh[i] = ((-edgeData.zeta[i] + 1) / 2 + part) / 2;
                    edgeRefData.eta1DNeigh[i] = 0;
                  }
                }
                if (neigh_edge == 1)
                {
                  for (size_t i = 0; i < N_Points1D; i++) // for all quadrature points
                  {
                    if (part == 1) part = 0;
                    edgeRefData.xi1DNeigh[i] = ((edgeData.zeta[i] + 1) / 2 - part) / 2;
                    if (part == 0) part = 1;
                    if (part == -1) part = 0;
                    edgeRefData.eta1DNeigh[i] = ((-edgeData.zeta[i] + 1) / 2 + part) / 2;
                    if (part == 0) part = -1;
                  }
                }
                if (neigh_edge == 2)
                {
                  for (size_t i = 0; i < N_Points1D; i++) // for all quadrature points
                  {
                    if (part == 1) part = 0;
                    edgeRefData.xi1DNeigh[i] = 0;
                    edgeRefData.eta1DNeigh[i] = ((edgeData.zeta[i] + 1) / 2 - part) / 2;
                  }
                }
                break;
                
              default:
                ;
            }
          }

          if (ee_verbose > 1)
            for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
              cout << "xiN " << edgeRefData.xi1DNeigh[i] << " etaN " << edgeRefData.eta1DNeigh[i] << endl;

          // compute gradients in reference cell of the neighbour
          for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
          {
            bf_neighbor->GetDerivatives(MultiIndex2D::D00, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i]);
            bf_neighbor->GetDerivatives(MultiIndex2D::D10, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.xideriv_refNeigh1D[BaseFunctNeigh][i]);
            bf_neighbor->GetDerivatives(MultiIndex2D::D01, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.etaderiv_refNeigh1D[BaseFunctNeigh][i]);
          }
          auto RefTransNeigh = element_neighbor.GetRefTransID();          // reftrafo of neighbour
          FEDatabase::SetCellForRefTrans(neigh, RefTransNeigh);

          auto DOF = fe_space_u->GetGlobalDOF(neigh_N_);
          for (int i = 0; i < N_Neigh; i++)
          {
            edgeRefData.FEFunctValuesNeigh[i] = values_u[DOF[i]];
            edgeRefData.FEFunctValuesNeigh[i + edgeRefData.max_n_base_functions_2d] = values_u[DOF[i] + N_U]; // v values
            if (ee_verbose > 1)
              cout << " value " << edgeRefData.FEFunctValuesNeigh[i] << endl;
          }
          for (size_t i = 0; i < N_Points1D; i++)  // get values and derivatives in original cell
          {
            FEDatabase::GetOrigValues(
              RefTransNeigh, edgeRefData.xi1DNeigh[i],
              edgeRefData.eta1DNeigh[i],
              bf_neighbor, collection_u, (TGridCell *) neigh,
              edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i],
              edgeRefData.xideriv_refNeigh1D[BaseFunctNeigh][i],
              edgeRefData.etaderiv_refNeigh1D[BaseFunctNeigh][i],
              edgeRefData.xyval_refNeigh1D[i],
              edgeRefData.xderiv_refNeigh1D[i],
              edgeRefData.yderiv_refNeigh1D[i]);
          }

          double val[6];
          for (size_t i = 0; i < N_Points1D; i++)     // for all quadrature points
          {
            val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = 0;
            for (int l = 0; l < N_Neigh; l++)       // for all basis functions
            {
              size_t m = l + edgeRefData.max_n_base_functions_2d;
              val[0] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.xderiv_refNeigh1D[i][l]; // accumulate value of derivative
              val[1] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.yderiv_refNeigh1D[i][l]; // accumulate value of derivative
              val[2] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.xyval_refNeigh1D[i][l]; // accumulate value of derivative
              val[3] += edgeRefData.FEFunctValuesNeigh[m] * edgeRefData.xderiv_refNeigh1D[i][l]; // accumulate value of derivative
              val[4] += edgeRefData.FEFunctValuesNeigh[m] * edgeRefData.yderiv_refNeigh1D[i][l]; // accumulate value of derivative
              val[5] += edgeRefData.FEFunctValuesNeigh[m] * edgeRefData.xyval_refNeigh1D[i][l]; // accumulate value of derivative
              if (ee_verbose > 1)
                cout << l << "  " << edgeRefData.xderiv_refNeigh1D[i][l] << "  " <<
                  edgeRefData.yderiv_refNeigh1D[i][l] << "  " << edgeRefData.FEFunctValuesNeigh[l] << endl;
            } // endfor l
            size_t m = i + N_Points1D;
            edgeRefData.xderiv_Neigh1D[i] = val[0]; // for k-th
            edgeRefData.yderiv_Neigh1D[i] = val[1]; // for k-th
            edgeRefData.xyval_Neigh1D[i] = val[2]; // for k-th
            edgeRefData.xderiv_Neigh1D[m] = val[3]; // for k-th
            edgeRefData.yderiv_Neigh1D[m] = val[4]; // for k-th
            edgeRefData.xyval_Neigh1D[m] = val[5]; // for k-th
          } // endfor i


          FEDatabase::GetOrigFromRef(RefTransNeigh, N_Points1D,
                                     edgeRefData.xi1DNeigh,
                                     edgeRefData.eta1DNeigh,
                                     edgeRefData.X1DNeigh,
                                     edgeRefData.Y1DNeigh);
          // pressure second
          // compute values at the quadrature points on the edge of
          // the neighbour element

          auto& eleNeigh = fe_space_p->get_fe(neigh_N_);

          BaseFunctNeigh = eleNeigh.GetBaseFunct_ID();    // basis functions on neighbout
          auto bfNeigh = eleNeigh.GetBaseFunct();
          N_Neigh = fe_space_p->get_n_local_dof(neigh_N_);                    // number of basis functions

          // compute gradients in reference cell of the neighbour
          for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
          {
            bfNeigh->GetDerivatives(MultiIndex2D::D00, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i]);
          }

          RefTransNeigh = eleNeigh.GetRefTransID();          // reftrafo of neighbour
          FEDatabase::SetCellForRefTrans(neigh, RefTransNeigh);

          auto DOFP = fe_space_p->GetGlobalDOF(neigh_N_);

          for (int i = 0; i < N_Neigh; i++)
          {
            edgeRefData.FEFunctValuesNeigh[i + 2 * edgeRefData.max_n_base_functions_2d] = values_p[DOFP[i]]; // p values
            if (ee_verbose > 1)
              cout << " p value " << edgeRefData.FEFunctValuesNeigh[i + 2 * edgeRefData.max_n_base_functions_2d] << endl;
          }
          for (size_t i = 0; i < N_Points1D; i++)  // get values and derivatives in original cell
          {
            FEDatabase::GetOrigValues(
              RefTransNeigh, edgeRefData.xi1DNeigh[i],
              edgeRefData.eta1DNeigh[i],
              bfNeigh, collection_u, (TGridCell *) neigh,
              edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i],
              edgeRefData.xideriv_refNeigh1D[BaseFunctNeigh][i],
              edgeRefData.etaderiv_refNeigh1D[BaseFunctNeigh][i],
              edgeRefData.xyval_refNeigh1D[i],
              edgeRefData.xderiv_refNeigh1D[i],
              edgeRefData.yderiv_refNeigh1D[i]);
          }

          size_t m;
          for (size_t i = 0; i < N_Points1D; i++)     // for all quadrature points on edge
          {
            val[0] = 0;
            for (int l = 0; l < N_Neigh; l++)       // for all basis functions
            {
              m = l + 2 * edgeRefData.max_n_base_functions_2d;
              val[0] += edgeRefData.FEFunctValuesNeigh[m] * edgeRefData.xyval_refNeigh1D[i][l]; // accumulate value of derivative
            } // endfor l
            m = i + 2 * N_Points1D;
            edgeRefData.xyval_Neigh1D[m] = val[0]; // for k-th
            //cout << "p value " << xyval_Neigh1D[m] << " " << xyval_1D[j][m] << endl;
          } // endfor i

          jump = 0.0;
          double absdetjk1D = hE / 2;
          for (size_t i = 0; i < N_Points1D; i++)           // compute jump
          {
            m = i + N_Points1D;
            l = m + N_Points1D;
            if ((std::abs(edgeData.XEdge1D[edgeIdx][i] - edgeRefData.X1DNeigh[i]) + std::abs(edgeData.YEdge1D[edgeIdx][i] - edgeRefData.Y1DNeigh[i])) > 1e-8)
              cout << " wrong quad points_a " << edgeData.XEdge1D[edgeIdx][i] << " , " << edgeData.YEdge1D[edgeIdx][i]
                << "   " << edgeRefData.X1DNeigh[i] << " , " << edgeRefData.Y1DNeigh[i] << endl;
            if (check_cont_u)
            {
              if (std::abs(edgeRefData.xyval_Neigh1D[i] - edgeData.xyval_1D[edgeIdx][i]) > 1e-8) {
                cout << " cell " << cell->GetCellIndex() << " i " << i << " uval_a " << edgeData.xyval_1D[edgeIdx][i] << " uneigh_a " << edgeRefData.xyval_Neigh1D[i] << endl;
              }
              if (std::abs(edgeRefData.xyval_Neigh1D[m] - edgeData.xyval_1D[edgeIdx][m]) > 1e-8) {
                cout << " cell " << cell->GetCellIndex() << " i " << i << " vval_a " << edgeData.xyval_1D[edgeIdx][m] << " vneigh_a " << edgeRefData.xyval_Neigh1D[m] << endl;
              }
            }
            if (check_cont_p)
            {
              if (std::abs(edgeRefData.xyval_Neigh1D[l] - edgeData.xyval_1D[edgeIdx][l]) > 1e-8)
                  cout << " i " << i << " pval_a " << edgeData.xyval_1D[edgeIdx][l] << " pneigh_a " << edgeRefData.xyval_Neigh1D[l] << endl;
            }
            auto e1 = coeff[0] * ((edgeData.xderiv_1D[edgeIdx][i] - edgeRefData.xderiv_Neigh1D[i]) * nx
                                  + (edgeData.yderiv_1D[edgeIdx][i] - edgeRefData.yderiv_Neigh1D[i]) * ny)
                      - (edgeData.xyval_1D[edgeIdx][l] - edgeRefData.xyval_Neigh1D[l]) * nx;
            auto e2 = coeff[0] * ((edgeData.xderiv_1D[edgeIdx][m] - edgeRefData.xderiv_Neigh1D[m]) * nx
                                  + (edgeData.yderiv_1D[edgeIdx][m] - edgeRefData.yderiv_Neigh1D[m]) * ny)
                      - (edgeData.xyval_1D[edgeIdx][l] - edgeRefData.xyval_Neigh1D[l]) * ny;
            if (ee_verbose > 1)
              cout << i << " jumpx " << edgeData.xderiv_1D[edgeIdx][i] << " " << edgeRefData.xderiv_Neigh1D[i] << endl;
            auto w = edgeData.weights1D[i] * absdetjk1D;
            jump += w * (e1 * e1 + e2 * e2);                       // integral on the edge
          }
          if (ee_verbose > 1)
            cout << "jump " << jump << endl;

          double beta[N_NSE2D_ESTIMATOR_TYPES];
          beta[0] = hE;                          // weight for H^1 estimator
          beta[1] = hE * hE * hE;                    // weight for L^2 estimator
          beta[2] = hE / coeff[0];                 // weight for energy norm estimator
          if (1.0 / std::sqrt(coeff[0]) < beta[2])
            beta[2] = 1.0 / std::sqrt(coeff[0]);
          for (size_t i = 1; i < 4; i++)
            estimated_error[i] += beta[i - 1] * jump / 2.0;

        }
        else
        {
          /*************************************************************************/
          /*  neighbour is not on the finest level, find children of neighbour     */
          /*************************************************************************/
          int n = neigh->GetClipBoard();
          if (n == -1)
          {
            // the neighbour is no member of the collection
            // check whether the children of neigh are in collection
            // find the local edge of neigh on which cell is -> l

            unsigned int edge2neigh = 0;
            while (neigh->GetJoint(edge2neigh)->GetNeighbour(neigh) != cell)
              edge2neigh++; // find connections between cells

            const int *TmpEdVerNeigh, *TmpoEnE, *TmpLen1, *TmpEC, *TmpLen2, *TmpoEnlE;
            int MaxLen1, MaxLen2;
            refdesc = neigh->GetRefDesc();                          // ref desc of neighbour
            refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVerNeigh);// get edges
            refdesc->GetOldEdgeNewEdge(TmpoEnE, TmpLen1, MaxLen1);// get connection to child edges
            refdesc->GetEdgeChild(TmpEC, TmpLen2, MaxLen2);       // get cell belonging to child edge (TmpEC)
            refdesc->GetOldEdgeNewLocEdge(TmpoEnlE);              // get local no.s of child edge

            // not general  !!!
            const unsigned int N_child = conform_grid ? 1 : 2;

            // find children of neigh on face l -> child
            for (size_t k = 0; k < N_child && TmpoEnE; k++)
            {
              // edge child, not general !!!
              auto edge1 = TmpoEnE[edge2neigh * MaxLen1 + k];
              // local number of child cell
              auto chnum1 = TmpEC[edge1 * MaxLen2];
              // child cell
              auto child = neigh->GetChild(chnum1);

              int MaxLen3;
              // get local indices of child edge
              const int *TmpECI, *TmpLen3;
              refdesc->GetEdgeChildIndex(TmpECI, TmpLen3, MaxLen3);
              // local index of child edge
              auto l_child = TmpECI[edge1 * MaxLen3];

              auto refdesc_child = child->GetRefDesc();          // ref desc of child
              refdesc_child->GetShapeDesc()->GetEdgeVertex(TmpEdVer);// conn. edge -> vertices
              auto ver2 = child->GetVertex(TmpEdVer[2 * l_child]);          // vertices of edge
              auto ver3 = child->GetVertex(TmpEdVer[2 * l_child + 1]);

              if (ee_verbose > 1)
              {
                cout << "ver 0 " << ver0->GetX() << "  " << ver0->GetY() << endl;
                cout << "ver 1 " << ver1->GetX() << "  " << ver1->GetY() << endl;
                cout << "ver 2 " << ver2->GetX() << "  " << ver2->GetY() << endl;
                cout << "ver 3 " << ver3->GetX() << "  " << ver3->GetY() << endl;
              }
              int part = 0;
              if (ver1 == ver2)
              {
                part = 1;
              }
              else if (ver0 == ver3)
              {
                part = -1;
              }
              else
              {
                cout << " something wrong 5 " << endl;
              }

              // now from point of view of child cell -> cell becomes the neighbour
              // prepare intergration for the half part of edge j

              int neigh_N_ = cell->GetClipBoard();           // number of original cell  in iterator
              if (neigh_N_ == -1)
              {
                cout << "Hier sollte man aber nicht hinkommen 33 !" << endl;
              }
              auto eleNeigh = fe_space_u->get_fe(neigh_N_);
              auto bfNeigh = eleNeigh.GetBaseFunct();
              auto BaseFunctNeigh = bfNeigh->GetID();    // basis functions on neighbour
              auto N_Neigh = eleNeigh.GetN_DOF();                    // number of basis functions

              auto bf2DrefelementsNeigh = bfNeigh->GetRefElement();   // reference cell of neighbour

              size_t neigh_edge = edgeIdx;
              if (conform_grid)
              {
                switch (bf2DrefelementsNeigh)                // compute coordinates of line quadrature
                {                                    // points in reference cell
                  case BFRefElements::BFUnitSquare :                  // edge 0
                    if (neigh_edge == 0)
                    {
                      for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                      {
                        edgeRefData.xi1DNeigh[i] = -edgeData.zeta[i];
                        edgeRefData.eta1DNeigh[i] = -1;
                      }
                    }
                    if (neigh_edge == 1)
                    {                               // edge 1
                      for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                      {
                        edgeRefData.xi1DNeigh[i] = 1;
                        edgeRefData.eta1DNeigh[i] = -edgeData.zeta[i];
                      }
                    }
                    if (neigh_edge == 2)
                    {                               // edge 2
                      for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                      {
                        edgeRefData.xi1DNeigh[i] = edgeData.zeta[i];
                        edgeRefData.eta1DNeigh[i] = 1;
                      }
                    }

                    if (neigh_edge == 3)
                    {                               // edge 3
                      for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                      {
                        edgeRefData.xi1DNeigh[i] = -1;
                        edgeRefData.eta1DNeigh[i] = edgeData.zeta[i];
                      }
                    }
                    break;

                  case BFRefElements::BFUnitTriangle :
                    if (neigh_edge == 0)
                    {
                      for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                      {
                        edgeRefData.xi1DNeigh[i] = (-edgeData.zeta[i] + 1) / 2;
                        edgeRefData.eta1DNeigh[i] = 0;
                      }
                    }
                    if (neigh_edge == 1)
                    {
                      for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                      {
                        edgeRefData.xi1DNeigh[i] = (edgeData.zeta[i] + 1) / 2;
                        edgeRefData.eta1DNeigh[i] = (-edgeData.zeta[i] + 1) / 2;
                      }
                    }
                    if (neigh_edge == 2)
                    {
                      for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                      {
                        edgeRefData.xi1DNeigh[i] = 0;
                        edgeRefData.eta1DNeigh[i] = (edgeData.zeta[i] + 1) / 2;
                      }
                    }
                    break;
                    
                  default:
                    ;
                }
              }
              else
              {
                // compute coordinates of line quadrature
                // this is only for 1-regular triangulations
                switch (bf2DrefelementsNeigh)
                {
                  // points in reference cell
                  case BFRefElements::BFUnitSquare :
                    // edge 0
                    if (neigh_edge == 0)
                    {
                      for (size_t i = 0; i < N_Points1D; i++)
                      {
                        edgeRefData.xi1DNeigh[i] = (-edgeData.zeta[i] + part) / 2;
                        edgeRefData.eta1DNeigh[i] = -1;
                      }
                    }
                    // edge 1
                    if (neigh_edge == 1)
                    {
                      // for all quadrature points
                      for (size_t i = 0; i < N_Points1D; i++)
                      {
                        edgeRefData.xi1DNeigh[i] = 1;
                        edgeRefData.eta1DNeigh[i] = (-edgeData.zeta[i] + part) / 2;
                      }
                    }
                    // edge 2
                    if (neigh_edge == 2)
                    {
                      for (size_t i = 0; i < N_Points1D; i++)
                      {
                        edgeRefData.xi1DNeigh[i] = (edgeData.zeta[i] - part) / 2;
                        edgeRefData.eta1DNeigh[i] = 1;
                      }
                    }
                    // edge 3
                    if (neigh_edge == 3)
                    {
                      for (size_t i = 0; i < N_Points1D; i++)
                      {
                        edgeRefData.xi1DNeigh[i] = -1;
                        edgeRefData.eta1DNeigh[i] = (edgeData.zeta[i] - part) / 2;
                      }
                    }
                    break;

                  case BFRefElements::BFUnitTriangle :
                    if (neigh_edge == 0)
                    {
                      for (size_t i = 0; i < N_Points1D; i++)
                      {
                        if (part == -1) part = 0;
                        edgeRefData.xi1DNeigh[i] = ((-edgeData.zeta[i] + 1) / 2 + part) / 2;
                        edgeRefData.eta1DNeigh[i] = 0;
                      }
                    }
                    if (neigh_edge == 1)
                    {
                      for (size_t i = 0; i < N_Points1D; i++)
                      {
                        if (part == 1) part = 0;
                        edgeRefData.xi1DNeigh[i] = ((edgeData.zeta[i] + 1) / 2 - part) / 2;
                        if (part == 0) part = 1;
                        if (part == -1) part = 0;
                        edgeRefData.eta1DNeigh[i] = ((-edgeData.zeta[i] + 1) / 2 + part) / 2;
                        if (part == 0) part = -1;
                      }
                    }
                    if (neigh_edge == 2)
                    {
                      for (size_t i = 0; i < N_Points1D; i++)
                      {
                        if (part == 1) part = 0;
                        edgeRefData.xi1DNeigh[i] = 0;
                        edgeRefData.eta1DNeigh[i] = ((edgeData.zeta[i] + 1) / 2 - part) / 2;
                      }
                    }
                    break;
                    
                  default:
                    ;
                }
              }

              if (ee_verbose > 1)
                for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                  cout << "xiN " << edgeRefData.xi1DNeigh[i] << " etaN " << edgeRefData.eta1DNeigh[i] << endl;

              // compute gradients in reference cell of the neighbour
              for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
              {
                bfNeigh->GetDerivatives(MultiIndex2D::D00, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i]);
                bfNeigh->GetDerivatives(MultiIndex2D::D10, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.xideriv_refNeigh1D[BaseFunctNeigh][i]);
                bfNeigh->GetDerivatives(MultiIndex2D::D01, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.etaderiv_refNeigh1D[BaseFunctNeigh][i]);
              }
              auto RefTransNeigh = eleNeigh.GetRefTransID();          // reftrafo of neighbour
              FEDatabase::SetCellForRefTrans(cell, RefTransNeigh);

              auto DOF = fe_space_u->GetGlobalDOF(neigh_N_);
              for (int i = 0; i < N_Neigh; i++)
              {
                edgeRefData.FEFunctValuesNeigh[i] = values_u[DOF[i]];
                edgeRefData.FEFunctValuesNeigh[i + edgeRefData.max_n_base_functions_2d] = values_u[DOF[i] + N_U]; // v values
                if (ee_verbose > 1)
                  cout << " value " << edgeRefData.FEFunctValuesNeigh[i] <<
                    " " << edgeRefData.FEFunctValuesNeigh[i + edgeRefData.max_n_base_functions_2d] << endl;
              }

              for (size_t i = 0; i < N_Points1D; i++)  // get values and derivatives in original cell
              {
                FEDatabase::GetOrigValues(
                  RefTransNeigh, edgeRefData.xi1DNeigh[i],
                  edgeRefData.eta1DNeigh[i],
                  bfNeigh, collection_u, (TGridCell *) neigh,
                  edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i],
                  edgeRefData.xideriv_refNeigh1D[BaseFunctNeigh][i],
                  edgeRefData.etaderiv_refNeigh1D[BaseFunctNeigh][i],
                  edgeRefData.xyval_refNeigh1D[i],
                  edgeRefData.xderiv_refNeigh1D[i],
                  edgeRefData.yderiv_refNeigh1D[i]);
              }

              double val[6];
              for (size_t i = 0; i < N_Points1D; i++)     // for all quadrature points
              {
                val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = 0;
                for (int l = 0; l < N_Neigh; l++)       // for all basis functions
                {
                  size_t m = l + edgeRefData.max_n_base_functions_2d;
                  val[0] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.xderiv_refNeigh1D[i][l]; // accumulate value of derivative
                  val[1] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.yderiv_refNeigh1D[i][l]; // accumulate value of derivative
                  val[2] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.xyval_refNeigh1D[i][l]; // accumulate value of derivative
                  val[3] += edgeRefData.FEFunctValuesNeigh[m] * edgeRefData.xderiv_refNeigh1D[i][l]; // accumulate value of derivative
                  val[4] += edgeRefData.FEFunctValuesNeigh[m] * edgeRefData.yderiv_refNeigh1D[i][l]; // accumulate value of derivative
                  val[5] += edgeRefData.FEFunctValuesNeigh[m] * edgeRefData.xyval_refNeigh1D[i][l]; // accumulate value of derivative
                  if (ee_verbose > 1)
                    cout << l << "  " << edgeRefData.xderiv_refNeigh1D[i][l] << "  " <<
                      edgeRefData.yderiv_refNeigh1D[i][l] << "  " << edgeRefData.FEFunctValuesNeigh[l] << endl;
                } // endfor l
                size_t m = i + N_Points1D;
                edgeRefData.xderiv_Cell1D[i] = val[0]; // for k-th
                edgeRefData.yderiv_Cell1D[i] = val[1]; // for k-th
                edgeRefData.xyval_Cell1D[i] = val[2]; // for k-th
                edgeRefData.xderiv_Cell1D[m] = val[3]; // for k-th
                edgeRefData.yderiv_Cell1D[m] = val[4]; // for k-th
                edgeRefData.xyval_Cell1D[m] = val[5]; // for k-th
              } // endfor i

              FEDatabase::GetOrigFromRef(RefTransNeigh, N_Points1D,
                                         edgeRefData.xi1DNeigh,
                                         edgeRefData.eta1DNeigh,
                                         edgeRefData.X1DCell,
                                         edgeRefData.Y1DCell);

              // pressure second
              // compute values at the quadrature points on the edge of
              // the cell

              auto eleNeigh_p = fe_space_p->get_fe(neigh_N_);
              auto bfNeigh_p = eleNeigh_p.GetBaseFunct();

              BaseFunctNeigh = bfNeigh_p->GetID();    // basis functions on neighbout
              N_Neigh = eleNeigh_p.GetN_DOF();                    // number of basis functions


              // compute gradients in reference cell of the neighbour
              for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
              {
                bfNeigh_p->GetDerivatives(MultiIndex2D::D00, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i]);
              }

              RefTransNeigh = eleNeigh_p.GetRefTransID();          // reftrafo of neighbour
              FEDatabase::SetCellForRefTrans(cell, RefTransNeigh);

              auto DOFP = fe_space_p->GetGlobalDOF(neigh_N_);

              for (int i = 0; i < N_Neigh; i++)
              {
                edgeRefData.FEFunctValuesNeigh[i + 2 * edgeRefData.max_n_base_functions_2d] = values_p[DOFP[i]]; // p values
                if (ee_verbose > 1)
                  cout << " p value " << edgeRefData.FEFunctValuesNeigh[i + 2 * edgeRefData.max_n_base_functions_2d] << endl;
              }
              for (size_t i = 0; i < N_Points1D; i++)  // get values and derivatives in original cell
              {
                FEDatabase::GetOrigValues(
                  RefTransNeigh, edgeRefData.xi1DNeigh[i],
                  edgeRefData.eta1DNeigh[i],
                  bfNeigh_p, collection_u, (TGridCell *) neigh,
                  edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i],
                  edgeRefData.xideriv_refNeigh1D[BaseFunctNeigh][i],
                  edgeRefData.etaderiv_refNeigh1D[BaseFunctNeigh][i],
                  edgeRefData.xyval_refNeigh1D[i],
                  edgeRefData.xderiv_refNeigh1D[i],
                  edgeRefData.yderiv_refNeigh1D[i]);
              }

              for (size_t i = 0; i < N_Points1D; i++)     // for all quadrature points on edge
              {
                val[0] = 0;
                for (int l = 0; l < N_Neigh; l++)       // for all basis functions
                {
                  size_t m = l + 2 * edgeRefData.max_n_base_functions_2d;
                  val[0] += edgeRefData.FEFunctValuesNeigh[m] * edgeRefData.xyval_refNeigh1D[i][l]; // accumulate value of derivative
                } // endfor l
                size_t m = i + 2 * N_Points1D;
                edgeRefData.xyval_Cell1D[m] = val[0]; // for k-th
                //cout << "p value " << xyval_Cell1D[m] << endl;
              } // endfor i

              // prepare integration for the child of the neighbour belong to the half part
              // of edge j

              neigh_N_ = child->GetClipBoard(); // number of neighbour in iterator
              eleNeigh = fe_space_u->get_fe(neigh_N_);
              bfNeigh = eleNeigh.GetBaseFunct();

              BaseFunctNeigh = bfNeigh->GetID();    // basis functions on neighbout
              N_Neigh = eleNeigh.GetN_DOF();                    // number of basis functions

              bf2DrefelementsNeigh = bfNeigh->GetRefElement();   // reference cell of neighbour

              neigh_edge = (size_t) l_child;
              switch (bf2DrefelementsNeigh)                // compute coordinates of line quadrature
              {                                    // points in reference cell
                case BFRefElements::BFUnitSquare :                  // edge 0
                  if (neigh_edge == 0)
                  {
                    for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                    {
                      edgeRefData.xi1DNeigh[i] = edgeData.zeta[i];
                      edgeRefData.eta1DNeigh[i] = -1;
                    }
                  }
                  if (neigh_edge == 1)
                  {                               // edge 1
                    for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                    {
                      edgeRefData.xi1DNeigh[i] = 1;
                      edgeRefData.eta1DNeigh[i] = edgeData.zeta[i];
                    }
                  }
                  if (neigh_edge == 2)
                  {                               // edge 2
                    for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                    {
                      edgeRefData.xi1DNeigh[i] = -edgeData.zeta[i];
                      edgeRefData.eta1DNeigh[i] = 1;
                    }
                  }

                  if (neigh_edge == 3)
                  {                               // edge 3
                    for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                    {
                      edgeRefData.xi1DNeigh[i] = -1;
                      edgeRefData.eta1DNeigh[i] = -edgeData.zeta[i];
                    }
                  }
                  break;

                case BFRefElements::BFUnitTriangle :
                  if (neigh_edge == 0)
                  {
                    for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                    {
                      edgeRefData.xi1DNeigh[i] = (edgeData.zeta[i] + 1) / 2;
                      edgeRefData.eta1DNeigh[i] = 0;
                    }
                  }
                  if (neigh_edge == 1)
                  {
                    for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                    {
                      edgeRefData.xi1DNeigh[i] = (-edgeData.zeta[i] + 1) / 2;
                      edgeRefData.eta1DNeigh[i] = (edgeData.zeta[i] + 1) / 2;
                    }
                  }
                  if (neigh_edge == 2)
                  {
                    for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                    {
                      edgeRefData.xi1DNeigh[i] = 0;
                      edgeRefData.eta1DNeigh[i] = (-edgeData.zeta[i] + 1) / 2;
                    }
                  }
                  break;
                  
                default:
                  ;
              }

              if (ee_verbose > 1)
                for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                  cout << "xiN " << edgeRefData.xi1DNeigh[i] << " etaN " << edgeRefData.eta1DNeigh[i] << endl;

              // compute gradients in reference cell of the neighbour
              for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
              {
                bfNeigh->GetDerivatives(MultiIndex2D::D00, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i]);
                bfNeigh->GetDerivatives(MultiIndex2D::D10, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.xideriv_refNeigh1D[BaseFunctNeigh][i]);
                bfNeigh->GetDerivatives(MultiIndex2D::D01, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.etaderiv_refNeigh1D[BaseFunctNeigh][i]);
              }
              RefTransNeigh = eleNeigh.GetRefTransID();          // reftrafo of neighbour
              FEDatabase::SetCellForRefTrans(child, RefTransNeigh);

              {
                DOF = fe_space_u->GetGlobalDOF(neigh_N_);
                for (int i = 0; i < N_Neigh; i++) {
                  edgeRefData.FEFunctValuesNeigh[i] = values_u[DOF[i]];       // u values
                  edgeRefData.FEFunctValuesNeigh[i + edgeRefData.max_n_base_functions_2d] = values_u[DOF[i] + N_U]; // v values
                  if (ee_verbose > 1)
                    cout << " value " << edgeRefData.FEFunctValuesNeigh[i] <<
                      " " << edgeRefData.FEFunctValuesNeigh[i + edgeRefData.max_n_base_functions_2d] << endl;
                  }
              }

              // get values and derivatives in original cell
              for (size_t i = 0; i < N_Points1D; i++)
              {
                FEDatabase::GetOrigValues(
                  RefTransNeigh, edgeRefData.xi1DNeigh[i],
                  edgeRefData.eta1DNeigh[i],
                  bfNeigh, collection_u, (TGridCell *) neigh,
                  edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i],
                  edgeRefData.xideriv_refNeigh1D[BaseFunctNeigh][i],
                  edgeRefData.etaderiv_refNeigh1D[BaseFunctNeigh][i],
                  edgeRefData.xyval_refNeigh1D[i],
                  edgeRefData.xderiv_refNeigh1D[i],
                  edgeRefData.yderiv_refNeigh1D[i]);
              }

              // for all quadrature points
              for (size_t i = 0; i < N_Points1D; i++)
              {
                val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = 0;
                // for all basis functions
                for (int l = 0; l < N_Neigh; l++)
                {
                  size_t m = l + edgeRefData.max_n_base_functions_2d;
                  val[0] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.xderiv_refNeigh1D[i][l];
                  val[1] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.yderiv_refNeigh1D[i][l];
                  val[2] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.xyval_refNeigh1D[i][l];
                  val[3] += edgeRefData.FEFunctValuesNeigh[m] * edgeRefData.xderiv_refNeigh1D[i][l];
                  val[4] += edgeRefData.FEFunctValuesNeigh[m] * edgeRefData.yderiv_refNeigh1D[i][l];
                  val[5] += edgeRefData.FEFunctValuesNeigh[m] * edgeRefData.xyval_refNeigh1D[i][l];
                  if (ee_verbose > 1)
                  {
                    cout << l << "  " << edgeRefData.xderiv_refNeigh1D[i][l] << "  " <<
                    edgeRefData.yderiv_refNeigh1D[i][l] << "  " << edgeRefData.FEFunctValuesNeigh[l] << endl;
                  }
                }
                auto m = i + N_Points1D;
                edgeRefData.xderiv_Neigh1D[i] = val[0];
                edgeRefData.yderiv_Neigh1D[i] = val[1];
                edgeRefData.xyval_Neigh1D[i] = val[2];
                edgeRefData.xderiv_Neigh1D[m] = val[3];
                edgeRefData.yderiv_Neigh1D[m] = val[4];
                edgeRefData.xyval_Neigh1D[m] = val[5];
              }

              FEDatabase::GetOrigFromRef(RefTransNeigh, N_Points1D,
                                         edgeRefData.xi1DNeigh,
                                         edgeRefData.eta1DNeigh,
                                         edgeRefData.X1DNeigh,
                                         edgeRefData.Y1DNeigh);

              // pressure second
              // compute values at the quadrature points on the edge of
              // the neighbour element

              eleNeigh_p = fe_space_p->get_fe(neigh_N_);
              bfNeigh_p = eleNeigh_p.GetBaseFunct();

              BaseFunctNeigh = bfNeigh_p->GetID();    // basis functions on neighbout
              N_Neigh = eleNeigh_p.GetN_DOF();                    // number of basis functions


              // compute gradients in reference cell of the neighbour
              for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
              {
                bfNeigh_p->GetDerivatives(MultiIndex2D::D00, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i]);
              }

              RefTransNeigh = eleNeigh_p.GetRefTransID();          // reftrafo of neighbour
              FEDatabase::SetCellForRefTrans(child, RefTransNeigh);

              DOFP = fe_space_p->GetGlobalDOF(neigh_N_);

              for (int i = 0; i < N_Neigh; i++)
              {
                edgeRefData.FEFunctValuesNeigh[i + 2 * edgeRefData.max_n_base_functions_2d] = values_p[DOFP[i]]; // p values
                if (ee_verbose > 1)
                  cout << " p value " << edgeRefData.FEFunctValuesNeigh[i + 2 * edgeRefData.max_n_base_functions_2d] << endl;
              }
              for (size_t i = 0; i < N_Points1D; i++)  // get values and derivatives in original cell
              {
                FEDatabase::GetOrigValues(
                  RefTransNeigh, edgeRefData.xi1DNeigh[i],
                  edgeRefData.eta1DNeigh[i],
                  bfNeigh_p, collection_u, (TGridCell *) neigh,
                  edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i],
                  edgeRefData.xideriv_refNeigh1D[BaseFunctNeigh][i],
                  edgeRefData.etaderiv_refNeigh1D[BaseFunctNeigh][i],
                  edgeRefData.xyval_refNeigh1D[i],
                  edgeRefData.xderiv_refNeigh1D[i],
                  edgeRefData.yderiv_refNeigh1D[i]);
              }

              for (size_t i = 0; i < N_Points1D; i++)     // for all quadrature points on edge
              {
                val[0] = 0;
                for (int l = 0; l < N_Neigh; l++)       // for all basis functions
                {
                  auto m = l + 2 * edgeRefData.max_n_base_functions_2d;
                  val[0] += edgeRefData.FEFunctValuesNeigh[m] * edgeRefData.xyval_refNeigh1D[i][l]; // accumulate value of derivative
                } // endfor l
                auto m = i + 2 * N_Points1D;
                edgeRefData.xyval_Neigh1D[m] = val[0]; // for k-th
                //cout << "p value " << xyval_Neigh1D[m] << " " << xyval_1D[j][m] << endl;
              } // endfor i

              jump = 0.0;
              auto absdetjk1D = hE / (2.0 * N_child);       // only half edge is considered
              for (size_t i = 0; i < N_Points1D; i++)           // compute jump
              {
                auto m = i + N_Points1D;
                auto l = m + N_Points1D;
                if ((std::abs(edgeRefData.X1DCell[i] - edgeRefData.X1DNeigh[i]) + std::abs(edgeRefData.Y1DCell[i] - edgeRefData.Y1DNeigh[i])) > 1e-8)
                    cout << " wrong quad points_c " << edgeData.XEdge1D[edgeIdx][i] << " , " << edgeData.YEdge1D[edgeIdx][i]
                    << "   " << edgeRefData.X1DNeigh[i] << " , " << edgeRefData.Y1DNeigh[i] << endl;
                if (check_cont_u)
                {
                  if (std::abs(edgeRefData.xyval_Neigh1D[i] - edgeRefData.xyval_Cell1D[i]) > 1e-8)
                  {
                    cout << " cell " << cell->GetCellIndex() << " i " << i << " uval_c " << edgeRefData.xyval_Cell1D[i] << " uneigh_c " << edgeRefData.xyval_Neigh1D[i] << endl;
                  }
                  if (std::abs(edgeRefData.xyval_Neigh1D[m] - edgeRefData.xyval_Cell1D[m]) > 1e-8)
                  {
                    cout << " cell " << cell->GetCellIndex() << " i " << i << " vval_c " << edgeRefData.xyval_Cell1D[m] << " vneigh_c " << edgeRefData.xyval_Neigh1D[m] << endl;
                    }
                }
                if (check_cont_p)
                {
                  if (std::abs(edgeRefData.xyval_Neigh1D[l] - edgeRefData.xyval_Cell1D[l]) > 1e-8)
                    cout << " i " << i << " pval_c " << edgeRefData.xyval_Cell1D[l] << " pneigh_c " << edgeRefData.xyval_Neigh1D[l] << endl;
                }
                /**
                 * TODO: This part was missing, intentionally?
                 */
                auto e1 = coeff[0] * ((edgeRefData.xderiv_Cell1D[i] - edgeRefData.xderiv_Neigh1D[i]) * nx
                                      + (edgeRefData.yderiv_Cell1D[i] - edgeRefData.yderiv_Neigh1D[i]) * ny)
                          - (edgeRefData.xyval_Cell1D[l] - edgeRefData.xyval_Neigh1D[l]) * nx;
                auto e2 = coeff[0] * ((edgeRefData.xderiv_Cell1D[m] - edgeRefData.xderiv_Neigh1D[m]) * nx
                                      + (edgeRefData.yderiv_Cell1D[m] - edgeRefData.yderiv_Neigh1D[m]) * ny)
                          - (edgeRefData.xyval_Cell1D[l] - edgeRefData.xyval_Neigh1D[l]) * ny;
                auto w = edgeData.weights1D[i] * absdetjk1D;
                jump += w * (e1 * e1 + e2 * e2);                       // integral on the edge
              }

              if (ee_verbose > 1)
                cout << "jump_c " << jump << endl;
              auto hE2 = hE / N_child;
              double beta[N_NSE2D_ESTIMATOR_TYPES];
              beta[0] = hE2;                          // weight for H^1 estimator
              beta[1] = hE2 * hE2 * hE2;                    // weight for L^2 estimator
              beta[2] = hE2 / coeff[0];                 // weight for energy norm estimator
              if (1.0 / std::sqrt(coeff[0]) < beta[2])
                  beta[2] = 1.0 / std::sqrt(coeff[0]);
              for (size_t i = 1; i < 4; i++)
                  estimated_error[i] += beta[i - 1] * jump / 2.0;
            }
          }
          else
          {
            /*************************************************************************/
            /*  neighbour is on the finest level                                     */
            /*************************************************************************/
            // find the finite element on the other side
            // find the local edge of neigh on which cell is -> l
            unsigned int neigh_edge = 0;
            const int *TmpEdVerNeigh;
            while (neigh->GetJoint(neigh_edge)->GetNeighbour(neigh) != cell) neigh_edge++;
            refdesc = neigh->GetRefDesc();
            refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVerNeigh);
            ver0 = cell->GetVertex(TmpEdVer[2 * edgeIdx]);
            ver1 = cell->GetVertex(TmpEdVer[2 * edgeIdx + 1]);
            auto ver2 = neigh->GetVertex(TmpEdVerNeigh[2 * neigh_edge]);          // vertices of edge
            auto ver3 = neigh->GetVertex(TmpEdVerNeigh[2 * neigh_edge + 1]);
            if (((ver0 == ver2) && (ver1 == ver3)) || ((ver0 == ver3) && (ver1 == ver2)));
            else
            {
              cout << "wrong edge " << endl;
            }
            // velocity first
            // compute gradient at the quadrature points on the edge of
            // the neighbour element

            auto neigh_N_ = neigh->GetClipBoard();                  // number of neighbour in iterator
            auto eleNeigh = fe_space_u->get_fe(neigh_N_);
            auto bfNeigh = eleNeigh.GetBaseFunct();
            auto BaseFunctNeigh = bfNeigh->GetID();    // basis functions on neighbout
            auto N_Neigh = eleNeigh.GetN_DOF();                    // number of basis functions

            auto bf2DrefelementsNeigh = bfNeigh->GetRefElement();   // referenz cell of neighbour

            switch (bf2DrefelementsNeigh)                // compute coordinates of line quadrature
            {                                    // points in reference cell
              case BFRefElements::BFUnitSquare :                  // edge 0
                if (neigh_edge == 0)
                {
                  for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                  {
                    edgeRefData.xi1DNeigh[i] = -edgeData.zeta[i];
                    edgeRefData.eta1DNeigh[i] = -1;
                  }
                }
                if (neigh_edge == 1)
                {                               // edge 1
                  for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                  {
                    edgeRefData.xi1DNeigh[i] = 1;
                    edgeRefData.eta1DNeigh[i] = -edgeData.zeta[i];
                  }
                }
                if (neigh_edge == 2)
                {                               // edge 2
                  for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                  {
                    edgeRefData.xi1DNeigh[i] = edgeData.zeta[i];
                    edgeRefData.eta1DNeigh[i] = 1;
                  }
                }

                if (neigh_edge == 3)
                {                               // edge 3
                  for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                  {
                    edgeRefData.xi1DNeigh[i] = -1;
                    edgeRefData.eta1DNeigh[i] = edgeData.zeta[i];
                  }
                }
                break;

              case BFRefElements::BFUnitTriangle :
                if (neigh_edge == 0)
                {
                  for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                  {
                    edgeRefData.xi1DNeigh[i] = (-edgeData.zeta[i] + 1) / 2;
                    edgeRefData.eta1DNeigh[i] = 0;
                  }
                }
                if (neigh_edge == 1)
                {
                  for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                  {
                    edgeRefData.xi1DNeigh[i] = (edgeData.zeta[i] + 1) / 2;
                    edgeRefData.eta1DNeigh[i] = (-edgeData.zeta[i] + 1) / 2;
                  }
                }
                if (neigh_edge == 2)
                {
                  for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                  {
                    edgeRefData.xi1DNeigh[i] = 0;
                    edgeRefData.eta1DNeigh[i] = (edgeData.zeta[i] + 1) / 2;
                  }
                }
                break;
              
              default:
                ;
            }

            if (ee_verbose > 1)
              for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
                cout << "xiN " << edgeRefData.xi1DNeigh[i] << " etaN " << edgeRefData.eta1DNeigh[i] << endl;

            // compute gradients in reference cell of the neighbour
            for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
            {
              bfNeigh->GetDerivatives(MultiIndex2D::D00, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i]);
              bfNeigh->GetDerivatives(MultiIndex2D::D10, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.xideriv_refNeigh1D[BaseFunctNeigh][i]);
              bfNeigh->GetDerivatives(MultiIndex2D::D01, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.etaderiv_refNeigh1D[BaseFunctNeigh][i]);
            }
            auto RefTransNeigh = eleNeigh.GetRefTransID();          // reftrafo of neighbour
            FEDatabase::SetCellForRefTrans(neigh, RefTransNeigh);

            auto DOF = fe_space_u->GetGlobalDOF(neigh_N_);
            for (int i = 0; i < N_Neigh; i++)
            {
              edgeRefData.FEFunctValuesNeigh[i] = values_u[DOF[i]];       // u values
              edgeRefData.FEFunctValuesNeigh[i + edgeRefData.max_n_base_functions_2d] = values_u[DOF[i] + N_U]; // v values
              if (ee_verbose > 1)
                cout << " value " << edgeRefData.FEFunctValuesNeigh[i] <<
                  " " << edgeRefData.FEFunctValuesNeigh[i + edgeRefData.max_n_base_functions_2d] << endl;
            }
            for (size_t i = 0; i < N_Points1D; i++)  // get values and derivatives in original cell
            {
              FEDatabase::GetOrigValues(
                RefTransNeigh, edgeRefData.xi1DNeigh[i],
                edgeRefData.eta1DNeigh[i],
                bfNeigh, collection_u, (TGridCell *) neigh,
                edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i],
                edgeRefData.xideriv_refNeigh1D[BaseFunctNeigh][i],
                edgeRefData.etaderiv_refNeigh1D[BaseFunctNeigh][i],
                edgeRefData.xyval_refNeigh1D[i],
                edgeRefData.xderiv_refNeigh1D[i],
                edgeRefData.yderiv_refNeigh1D[i]);
            }

            double val[6];
            for (size_t i = 0; i < N_Points1D; i++)     // for all quadrature points on edge
            {
              val[0] = val[1] = val[2] = val[3] = val[4] = val[5] = 0;
              for (int l = 0; l < N_Neigh; l++)       // for all basis functions
              {
                auto m = l + edgeRefData.max_n_base_functions_2d;
                val[0] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.xderiv_refNeigh1D[i][l]; // accumulate value of derivative
                val[1] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.yderiv_refNeigh1D[i][l]; // accumulate value of derivative
                val[2] += edgeRefData.FEFunctValuesNeigh[l] * edgeRefData.xyval_refNeigh1D[i][l]; // accumulate value of derivative
                val[3] += edgeRefData.FEFunctValuesNeigh[m] * edgeRefData.xderiv_refNeigh1D[i][l]; // accumulate value of derivative
                val[4] += edgeRefData.FEFunctValuesNeigh[m] * edgeRefData.yderiv_refNeigh1D[i][l]; // accumulate value of derivative
                val[5] += edgeRefData.FEFunctValuesNeigh[m] * edgeRefData.xyval_refNeigh1D[i][l]; // accumulate value of derivative
                if (ee_verbose > 1)
                    cout << l << "  " << edgeRefData.xderiv_refNeigh1D[i][l] << "  " <<
                    edgeRefData.yderiv_refNeigh1D[i][l] << "  " << edgeRefData.FEFunctValuesNeigh[l] << endl;
              } // endfor l
              auto m = i + N_Points1D;
              edgeRefData.xderiv_Neigh1D[i] = val[0]; // for k-th
              edgeRefData.yderiv_Neigh1D[i] = val[1]; // for k-th
              edgeRefData.xyval_Neigh1D[i] = val[2]; // for k-th
              edgeRefData.xderiv_Neigh1D[m] = val[3]; // for k-th
              edgeRefData.yderiv_Neigh1D[m] = val[4]; // for k-th
              edgeRefData.xyval_Neigh1D[m] = val[5]; // for k-th
            } // endfor i

            // just for testing quad points, may be deleted later
            FEDatabase::GetOrigFromRef(RefTransNeigh, N_Points1D,
                                       edgeRefData.xi1DNeigh,
                                       edgeRefData.eta1DNeigh,
                                       edgeRefData.X1DNeigh,
                                       edgeRefData.Y1DNeigh);
            // pressure second
            // compute values at the quadrature points on the edge of
            // the neighbour element

            auto eleNeigh_p = fe_space_p->get_fe(neigh_N_);
            auto bfNeigh_p = eleNeigh_p.GetBaseFunct();
            BaseFunctNeigh = bfNeigh_p->GetID();    // basis functions on neighbout
            N_Neigh = eleNeigh_p.GetN_DOF();                    // number of basis functions

            // compute gradients in reference cell of the neighbour
            for (size_t i = 0; i < N_Points1D; i++)         // for all quadrature points
            {
              bfNeigh_p->GetDerivatives(MultiIndex2D::D00, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i]);
            }

            RefTransNeigh = eleNeigh_p.GetRefTransID();          // reftrafo of neighbour
            FEDatabase::SetCellForRefTrans(neigh, RefTransNeigh);

            auto DOFP = fe_space_p->GetGlobalDOF(neigh_N_);

            for (int i = 0; i < N_Neigh; i++)
            {
              edgeRefData.FEFunctValuesNeigh[i + 2 * edgeRefData.max_n_base_functions_2d] = values_p[DOFP[i]]; // p values
              if (ee_verbose > 1)
                cout << " p value " << edgeRefData.FEFunctValuesNeigh[i + 2 * edgeRefData.max_n_base_functions_2d] << endl;
            }
            for (size_t i = 0; i < N_Points1D; i++)  // get values and derivatives in original cell
            {
              FEDatabase::GetOrigValues(
                RefTransNeigh, edgeRefData.xi1DNeigh[i],
                edgeRefData.eta1DNeigh[i],
                bfNeigh_p, collection_u, (TGridCell *) neigh,
                edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i],
                edgeRefData.xideriv_refNeigh1D[BaseFunctNeigh][i],
                edgeRefData.etaderiv_refNeigh1D[BaseFunctNeigh][i],
                edgeRefData.xyval_refNeigh1D[i],
                edgeRefData.xderiv_refNeigh1D[i],
                edgeRefData.yderiv_refNeigh1D[i]);
            }

            for (size_t i = 0; i < N_Points1D; i++)     // for all quadrature points on edge
            {
              val[0] = 0;
              for (int l = 0; l < N_Neigh; l++)       // for all basis functions
              {
                auto m = l + 2 * edgeRefData.max_n_base_functions_2d;
                val[0] += edgeRefData.FEFunctValuesNeigh[m] * edgeRefData.xyval_refNeigh1D[i][l]; // accumulate value of derivative
              } // endfor l
              auto m = i + 2 * N_Points1D;
              edgeRefData.xyval_Neigh1D[m] = val[0]; // for k-th
              //cout << "p value " << xyval_Neigh1D[m] << " " << xyval_1D[j][m] << endl;
            } // endfor i

            jump = 0.0;
            auto absdetjk1D = hE / 2;
            for (size_t i = 0; i < N_Points1D; i++)           // compute jump
            {
              auto m = i + N_Points1D;
              auto l = m + N_Points1D;

              if ((std::abs(edgeData.XEdge1D[edgeIdx][i] - edgeRefData.X1DNeigh[i]) + std::abs(edgeData.YEdge1D[edgeIdx][i] - edgeRefData.Y1DNeigh[i])) > 2e-8)
              {
                cout << " wrong quad points_b " << edgeData.XEdge1D[edgeIdx][i] << " , " << edgeData.YEdge1D[edgeIdx][i] << "   " << edgeRefData.X1DNeigh[i] << " , " << edgeRefData.Y1DNeigh[i] << endl;
                cout << "\tby : " << (std::abs(edgeData.XEdge1D[edgeIdx][i] - edgeRefData.X1DNeigh[i]) + std::abs(edgeData.YEdge1D[edgeIdx][i] - edgeRefData.Y1DNeigh[i])) << endl;
              }
              if (check_cont_u)
              {
                if (std::abs(edgeRefData.xyval_Neigh1D[i] - edgeData.xyval_1D[edgeIdx][i]) > 1e-8)
                  cout << " i " << i << " uval_b " << edgeData.xyval_1D[edgeIdx][i] << " uneigh_b " << edgeRefData.xyval_Neigh1D[i] << endl;
                if (std::abs(edgeRefData.xyval_Neigh1D[m] - edgeData.xyval_1D[edgeIdx][m]) > 1e-8)
                  cout << " i " << i << " vval_b " << edgeData.xyval_1D[edgeIdx][m] << " vneigh_b " << edgeRefData.xyval_Neigh1D[m] << endl;
              }
              if (check_cont_p)
              {
                if (std::abs(edgeRefData.xyval_Neigh1D[l] - edgeData.xyval_1D[edgeIdx][l]) > 1e-8)
                  cout << " i " << i << " pval_b " << edgeData.xyval_1D[edgeIdx][l] << " pneigh_b " << edgeRefData.xyval_Neigh1D[l] << endl;
              }
              auto e1 = coeff[0] * ((edgeData.xderiv_1D[edgeIdx][i] - edgeRefData.xderiv_Neigh1D[i]) * nx
                                    + (edgeData.yderiv_1D[edgeIdx][i] - edgeRefData.yderiv_Neigh1D[i]) * ny)
                        - (edgeData.xyval_1D[edgeIdx][l] - edgeRefData.xyval_Neigh1D[l]) * nx;
              auto e2 = coeff[0] * ((edgeData.xderiv_1D[edgeIdx][m] - edgeRefData.xderiv_Neigh1D[m]) * nx
                                    + (edgeData.yderiv_1D[edgeIdx][m] - edgeRefData.yderiv_Neigh1D[m]) * ny)
                        - (edgeData.xyval_1D[edgeIdx][l] - edgeRefData.xyval_Neigh1D[l]) * ny;
              if (ee_verbose > 1)
                cout << i << " jumpx " << edgeData.xderiv_1D[edgeIdx][i] << " " << edgeRefData.xderiv_Neigh1D[i] << endl;
              auto w = edgeData.weights1D[i] * absdetjk1D;
              jump += w * (e1 * e1 + e2 * e2);                       // integral on the edge
            }
            // stelle 1
            if (ee_verbose > 1)
              cout << "jump " << jump << endl;
            double beta[N_NSE2D_ESTIMATOR_TYPES];
            beta[0] = hE;                          // weight for H^1 estimator
            beta[1] = hE * hE * hE;                    // weight for L^2 estimator
            beta[2] = hE / coeff[0];                 // weight for energy norm estimator
            if (1.0 / std::sqrt(coeff[0]) < beta[2])
                beta[2] = 1.0 / std::sqrt(coeff[0]);
            for (size_t i = 1; i < 4; i++)
                estimated_error[i] += beta[i - 1] * jump / 2.0;


          }  // end neighbour is member of the collection
        } // end neighbour on the finer level
      }     // end inner edge
    }         // end for j
  }

  for (int i = 0; i < N_NSE2D_ESTIMATOR_TYPES; i++)
  {
    estimated_local_error[i] = estimated_error[i];
  }
}

template<int d>
void NSEErrorEstimator<d>::info()
{
  if(this->currentCollection == nullptr)
  {
    Output::print("NSEErrorEstimator<", d, ">: the method 'estimate' has not "
                  "been called yet, therefore there is no further information "
                  "available on this object.");
    return;
  }
  Output::print("NSEErrorEstimator<", d, "> for ",
                this->currentCollection->GetN_Cells(), " cells");
  Output::dash("maximal_local_error: ", this->maximal_local_error);
  Output::dash("estimated global error: ", this->estimated_global_error);
}


std::ostream &operator<<(std::ostream &os, NSE2DErrorEstimatorType &type)
{
  switch (type)
  {
    case NSE2DErrorEstimatorType::gradient_indicator:
      os << "Gradient indicator";
      break;
    case NSE2DErrorEstimatorType::residual_estimator_h1:
      os << "H1 residual estimator";
      break;
    case NSE2DErrorEstimatorType::residual_estimator_l2:
      os << "L2 residual estimator";
      break;
    case NSE2DErrorEstimatorType::residual_estimator_energy_quasi_robust:
      os << "Energy residual estimator quasi robust";
      break;
    case NSE2DErrorEstimatorType::gradient_recovery:
      os << "Gradient recovery";
      break;
    case NSE2DErrorEstimatorType::implicit_estimator_neumann:
      os << "Implicit neumann estimator";
      break;
  }
  os << " (database value = " << int(type) << ")";
  return os;
}

#ifdef __3D__
//template class NSEErrorEstimator<3>;
#else
template class NSEErrorEstimator<2>;
#endif
