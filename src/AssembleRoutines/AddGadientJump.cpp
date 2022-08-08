#include "AddGadientJump.hpp"
#include "array"
#include "vector"
#include "Database.h"
#include "FEDatabase.h"
#include "QuadratureFormulaDatabase.h"
#include "memory"
#include <BoundEdge.h>
#include <cmath>
#include "BaseCell.h"

struct JointData
{
  protected:
    std::vector<double> xieta_ref1D_data;
  public:
    // constants
    constexpr static int MaxN_QuadPoints_1D_loc = MaxN_QuadPoints_1D;
    constexpr static int N_BaseFuncts2D_loc = N_BaseFuncts2D;
    
    std::array<std::vector<double>, 4> XEdge1D{}, YEdge1D{};//, X1DNeigh{}, Y1DNeigh{}; 
    
    // array holding the line quadrature points on reference cell
    double xi1D[N_BaseFuncts2D_loc][4][MaxN_QuadPoints_1D_loc]{};
    double eta1D[N_BaseFuncts2D_loc][4][MaxN_QuadPoints_1D_loc]{};
    
    std::vector<std::array<std::array<double *, MaxN_QuadPoints_1D_loc>, 4>> xietaval_ref1D{N_BaseFuncts2D_loc};
    std::vector<std::array<std::array<double *, MaxN_QuadPoints_1D_loc>, 4>> xideriv_ref1D{N_BaseFuncts2D_loc};
    std::vector<std::array<std::array<double *, MaxN_QuadPoints_1D_loc>, 4>> etaderiv_ref1D{N_BaseFuncts2D_loc};
    //values and derivative in reference cell
    std::vector<std::array<std::array<double, MaxN_QuadPoints_1D_loc>, 4>> xyval_ref1D{};
    std::vector<std::array<std::array<double, MaxN_QuadPoints_1D_loc>, 4>> xderiv_ref1D{};
    std::vector<std::array<std::array<double, MaxN_QuadPoints_1D_loc>, 4>> yderiv_ref1D{};
    
    // values and derivative values in original cell
    std::array<std::vector<double>, 4> xyval_1D{};
    std::array<std::vector<double>, 4> xderiv_1D{};
    std::array<std::vector<double>, 4> yderiv_1D{};
    
    JointData()
    {
      xieta_ref1D_data.resize(N_BaseFuncts2D_loc * 4 * MaxN_QuadPoints_1D_loc * N_BaseFuncts2D_loc * 3);
      
      auto *ptr = xieta_ref1D_data.data();
      xietaval_ref1D[0][0][0] = &ptr[0];
      xideriv_ref1D[0][0][0] = &ptr[N_BaseFuncts2D_loc * 4 * MaxN_QuadPoints_1D_loc * N_BaseFuncts2D_loc];
      etaderiv_ref1D[0][0][0] = &ptr[N_BaseFuncts2D_loc * 4 * MaxN_QuadPoints_1D_loc * N_BaseFuncts2D_loc * 2];
      
      for(size_t ii = 0; ii < N_BaseFuncts2D_loc; ii++)
      {
        for(size_t jj = 0; jj < 4; jj++)
        {
          for(size_t kk = 0; kk < MaxN_QuadPoints_1D_loc; kk++)
          {
            xietaval_ref1D[ii][jj][kk] = xietaval_ref1D[0][0][0] + (kk * N_BaseFuncts2D_loc + jj * MaxN_QuadPoints_1D_loc * N_BaseFuncts2D_loc + ii * 4 * MaxN_QuadPoints_1D_loc * N_BaseFuncts2D_loc);
              xideriv_ref1D[ii][jj][kk] = xideriv_ref1D[0][0][0] + (kk * N_BaseFuncts2D_loc + jj * MaxN_QuadPoints_1D_loc * N_BaseFuncts2D_loc + ii * 4 * MaxN_QuadPoints_1D_loc * N_BaseFuncts2D_loc);
              etaderiv_ref1D[ii][jj][kk] = etaderiv_ref1D[0][0][0] + (kk * N_BaseFuncts2D_loc + jj * MaxN_QuadPoints_1D_loc * N_BaseFuncts2D_loc + ii * 4 * MaxN_QuadPoints_1D_loc * N_BaseFuncts2D_loc);
          }
        }
      }
    }
    
    void updateEdgeData(unsigned int max_n_base_funct_2d)
    {
      xietaval_ref1D.resize(max_n_base_funct_2d);
      xideriv_ref1D.resize(max_n_base_funct_2d);
      etaderiv_ref1D.resize(max_n_base_funct_2d);
      xyval_ref1D.resize(max_n_base_funct_2d);
      xderiv_ref1D.resize(max_n_base_funct_2d);
      yderiv_ref1D.resize(max_n_base_funct_2d);
      
      for(unsigned int i = 0; i < 4; i++)
      {
        XEdge1D[i] = std::vector<double>(max_n_base_funct_2d);
        YEdge1D[i] = std::vector<double>(max_n_base_funct_2d);
        xyval_1D[i] = std::vector<double>{};
        xderiv_1D[i] = std::vector<double>{};
        yderiv_1D[i] = std::vector<double>{};
      }
    }
};

struct JointRefData
{
  protected:
    std::vector<double> refNeigh1D_data;
    unsigned int refDataQuadPoints1D;
  public:
    constexpr static int MaxN_QuadPoints_1D_loc = MaxN_QuadPoints_1D;
    
    unsigned int max_n_base_funct_2d;
    std::array<std::array<double *, MaxN_QuadPoints_1D_loc>, MaxN_BaseFunctions2D> xietaval_refNeigh1D;
    std::array<std::array<double *, MaxN_QuadPoints_1D_loc>, MaxN_BaseFunctions2D> xideriv_refNeigh1D;
    std::array<std::array<double *, MaxN_QuadPoints_1D_loc>, MaxN_BaseFunctions2D> etaderiv_refNeigh1D;
    
    std::vector<double> xi1DNeigh;
    std::vector<double> eta1DNeigh;
    std::vector<double> X1DNeigh;
    std::vector<double> Y1DNeigh;
    std::vector<double> X1DCell;
    std::vector<double> Y1DCell;
    std::vector<double> xderiv_Neigh1D;
    std::vector<double> yderiv_Neigh1D;
    std::vector<double> xyval_Neigh1D;
    std::vector<double> xderiv_Cell1D;
    std::vector<double> yderiv_Cell1D;
    std::vector<double> xyval_Cell1D;
    std::array<double *, MaxN_QuadPoints_1D_loc> xyval_refNeigh1D{};
    std::array<double *, MaxN_QuadPoints_1D_loc> xderiv_refNeigh1D{};
    std::array<double *, MaxN_QuadPoints_1D_loc> yderiv_refNeigh1D{};
    std::vector<double> FEFunctValuesNeigh;
    
    JointRefData (unsigned int max_n_base_funct_2d)
     : refDataQuadPoints1D{0}, max_n_base_funct_2d{max_n_base_funct_2d}
     {
       unsigned int k = MaxN_BaseFunctions2D * MaxN_QuadPoints_1D_loc * max_n_base_funct_2d;
       refNeigh1D_data = std::vector<double>(3 * k);
       xietaval_refNeigh1D[0][0] = &refNeigh1D_data[0];
       xideriv_refNeigh1D[0][0] = &refNeigh1D_data[0] + k;
       etaderiv_refNeigh1D[0][0] = &refNeigh1D_data[0] + 2 * k;
       for(unsigned int i = 0; i < MaxN_BaseFunctions2D; i++)
       {
         unsigned int n = i * MaxN_BaseFunctions2D * max_n_base_funct_2d;
         for(unsigned int j = 0; j < MaxN_QuadPoints_1D_loc; j++)
         {
           unsigned int l = n + j * max_n_base_funct_2d;
           xietaval_refNeigh1D[i][j] = xietaval_refNeigh1D[0][0] + l;
           xideriv_refNeigh1D[i][j] = xideriv_refNeigh1D[0][0] + l;
           etaderiv_refNeigh1D[i][j] = etaderiv_refNeigh1D[0][0] + l;
         }
       }
     }
     
     void updateQuadPointData(unsigned int refDataQuadPoints1D)
    {
      if (this->refDataQuadPoints1D != refDataQuadPoints1D)
      {
        this->refDataQuadPoints1D = refDataQuadPoints1D;
        xi1DNeigh.resize(refDataQuadPoints1D);
        eta1DNeigh.resize(refDataQuadPoints1D);
        X1DNeigh.resize(refDataQuadPoints1D);
        Y1DNeigh.resize(refDataQuadPoints1D);
        X1DCell.resize(refDataQuadPoints1D);
        Y1DCell.resize(refDataQuadPoints1D);
        xderiv_Neigh1D.resize(refDataQuadPoints1D);
        yderiv_Neigh1D.resize(refDataQuadPoints1D);
        xyval_Neigh1D.resize(refDataQuadPoints1D);
        xderiv_Cell1D.resize(refDataQuadPoints1D);
        yderiv_Cell1D.resize(refDataQuadPoints1D);
        xyval_Cell1D.resize(refDataQuadPoints1D);     
        
        FEFunctValuesNeigh.resize((3 * refDataQuadPoints1D + 1) * MaxN_BaseFunctions2D);
        unsigned int k = MaxN_BaseFunctions2D;
        for(size_t i = 0; i < refDataQuadPoints1D; i++)
        {
          xyval_refNeigh1D[i] = &FEFunctValuesNeigh.data()[0] + k;
          k += MaxN_BaseFunctions2D;
          xderiv_refNeigh1D[i] = &FEFunctValuesNeigh.data()[0] + k;
          k += MaxN_BaseFunctions2D;
          yderiv_refNeigh1D[i] = &FEFunctValuesNeigh.data()[0] + k;
          k += MaxN_BaseFunctions2D;
        }
      }
    }
};

void get_values_ref(BFRefElements bfref, const TQuadFormula& qf, const BaseFunctions* bf,
  JointData& edgeData, int bf_loc)
{
  unsigned int N_Points1D = qf.GetN_QuadPoints();
  switch(bfref)
  {
    case BFRefElements::BFUnitSquare:
    {
      //edge 0
      for(auto j = 0u; j < N_Points1D; j++)
      {
        double zeta = qf.get_point(j).x;
        edgeData.xi1D[bf_loc][0][j] = zeta;
        edgeData.eta1D[bf_loc][0][j] = -1;
        
        bf->GetDerivatives(MultiIndex2D::D00, zeta, -1, edgeData.xietaval_ref1D[bf_loc][0][j]);
        bf->GetDerivatives(MultiIndex2D::D10, zeta, -1, edgeData.xideriv_ref1D[bf_loc][0][j]);
        bf->GetDerivatives(MultiIndex2D::D01, zeta, -1, edgeData.etaderiv_ref1D[bf_loc][0][j]);
      }
      //edge 1
      for(auto j = 0u; j < N_Points1D; j++)
      {
        double zeta = qf.get_point(j).x;
        edgeData.xi1D[bf_loc][1][j] = 1;
        edgeData.eta1D[bf_loc][1][j] = zeta;
        bf->GetDerivatives(MultiIndex2D::D00, 1, zeta, edgeData.xietaval_ref1D[bf_loc][1][j]);
        bf->GetDerivatives(MultiIndex2D::D10, 1, zeta, edgeData.xideriv_ref1D[bf_loc][1][j]);
        bf->GetDerivatives(MultiIndex2D::D01, 1, zeta, edgeData.etaderiv_ref1D[bf_loc][1][j]);
      }
      // edge 2
      for(auto j = 0u; j < N_Points1D; j++)
      {
        double zeta = qf.get_point(j).x;
        edgeData.xi1D[bf_loc][2][j] = -zeta;
        edgeData.eta1D[bf_loc][2][j] = 1;
        bf->GetDerivatives(MultiIndex2D::D00, -zeta, 1, edgeData.xietaval_ref1D[bf_loc][2][j]);
        bf->GetDerivatives(MultiIndex2D::D10, -zeta, 1, edgeData.xideriv_ref1D[bf_loc][2][j]);
        bf->GetDerivatives(MultiIndex2D::D01, -zeta, 1, edgeData.etaderiv_ref1D[bf_loc][2][j]);
      }
      // edge 3
      for(auto j = 0u; j < N_Points1D; j++)
      {
        double zeta = qf.get_point(j).x;
        edgeData.xi1D[bf_loc][3][j] = -1;
        edgeData.eta1D[bf_loc][3][j] = -zeta;
        bf->GetDerivatives(MultiIndex2D::D00, -1, -zeta, edgeData.xietaval_ref1D[bf_loc][3][j]);
        bf->GetDerivatives(MultiIndex2D::D10, -1, -zeta, edgeData.xideriv_ref1D[bf_loc][3][j]);
        bf->GetDerivatives(MultiIndex2D::D01, -1, -zeta, edgeData.etaderiv_ref1D[bf_loc][3][j]);
      }
      break;
    }
    case BFRefElements::BFUnitTriangle:
    {
      // edge 0
      for(auto j = 0u; j < N_Points1D; j++)
      {
        double zeta = qf.get_point(j).x;
        edgeData.xi1D[bf_loc][0][j] = (zeta + 1) / 2;
        edgeData.eta1D[bf_loc][0][j] = 0;
        bf->GetDerivatives(MultiIndex2D::D00, (zeta + 1) / 2, 0, edgeData.xietaval_ref1D[bf_loc][0][j]);
        bf->GetDerivatives(MultiIndex2D::D10, (zeta + 1) / 2, 0, edgeData.xideriv_ref1D[bf_loc][0][j]);
        bf->GetDerivatives(MultiIndex2D::D01, (zeta + 1) / 2, 0, edgeData.etaderiv_ref1D[bf_loc][0][j]);
      }
      // edge 1
      for(auto j = 0u; j < N_Points1D; j++)
      {
        double zeta = qf.get_point(j).x;        
        edgeData.xi1D[bf_loc][1][j] = (-zeta + 1) / 2;
        edgeData.eta1D[bf_loc][1][j] = (zeta + 1) / 2;
        bf->GetDerivatives(MultiIndex2D::D00,(-zeta + 1)/2, (zeta + 1)/2,edgeData.xietaval_ref1D[bf_loc][1][j]);
        bf->GetDerivatives(MultiIndex2D::D10,(-zeta + 1)/2, (zeta + 1)/2,edgeData.xideriv_ref1D[bf_loc][1][j]);
        bf->GetDerivatives(MultiIndex2D::D01,(-zeta + 1)/2, (zeta + 1)/2,edgeData.etaderiv_ref1D[bf_loc][1][j]);
      }
      // edge 2
      for(auto j = 0u; j < N_Points1D; j++)
      {
        double zeta = qf.get_point(j).x;
        edgeData.xi1D[bf_loc][2][j] = 0;
        edgeData.eta1D[bf_loc][2][j] = (-zeta + 1) / 2;
        bf->GetDerivatives(MultiIndex2D::D00, 0, (-zeta + 1) / 2, edgeData.xietaval_ref1D[bf_loc][2][j]);
        bf->GetDerivatives(MultiIndex2D::D10, 0, (-zeta + 1) / 2, edgeData.xideriv_ref1D[bf_loc][2][j]);
        bf->GetDerivatives(MultiIndex2D::D01, 0, (-zeta + 1) / 2, edgeData.etaderiv_ref1D[bf_loc][2][j]);
      }
      break;
    }
    case BFRefElements::BFUnitHexahedron:
    case BFRefElements::BFUnitTetrahedron:
    case BFRefElements::BFUnitLine:
      ErrThrow("ComputeAddStab: 3D case is not supported ye");
      break;
  }
  Output::print("get_values_ref done");
}

void eta_to_edge_ref(BFRefElements bf2DrefelementsNeigh, int neigh_edge,
                     const TQuadFormula& qf, JointRefData& edgeData)
{
  unsigned int N_Points1D = qf.GetN_QuadPoints();
  switch(bf2DrefelementsNeigh)
  {
    case BFRefElements::BFUnitSquare:
    {
      // edge 0
      if(neigh_edge == 0)
      {
        for(unsigned int i = 0; i < N_Points1D; i++)
        {
          edgeData.xi1DNeigh[i] = -qf.get_point(i).x;
          edgeData.eta1DNeigh[i] = -1;
        }
      }
      // edge 1
      if(neigh_edge == 1)
      {
        for(unsigned int i = 0; i < N_Points1D; i++)
        {
          edgeData.xi1DNeigh[i] = 1;
          edgeData.eta1DNeigh[i] =-qf.get_point(i).x;
        }
      }
      // edge 2
      if(neigh_edge == 2)
      {
        for(unsigned int i = 0; i < N_Points1D; i++)
        {
          edgeData.xi1DNeigh[i] = qf.get_point(i).x;
          edgeData.eta1DNeigh[i] = 1;
        }
      }
      // edge 3
      if(neigh_edge == 3)
      {
        for(unsigned int i = 0; i < N_Points1D; i++)
        {
          edgeData.xi1DNeigh[i] = -1;
          edgeData.eta1DNeigh[i] = qf.get_point(i).x;
        }
      }
      break;
    }
    case BFRefElements::BFUnitTriangle:
    {
      //edge 0
      if(neigh_edge == 0)
      {
        for(unsigned int i = 0; i < N_Points1D; i++)
        {
          // for all quadrature points
          edgeData.xi1DNeigh[i] = (-qf.get_point(i).x + 1) / 2;
          edgeData.eta1DNeigh[i] = 0;
        }
      }
      // edge 1
      if(neigh_edge == 1)
      {
        for(unsigned int i = 0; i < N_Points1D; i++)
        {
          // for all quadrature points
          edgeData.xi1DNeigh[i] = (qf.get_point(i).x + 1) / 2;
          edgeData.eta1DNeigh[i] = (-qf.get_point(i).x + 1) / 2;
        }
      }
      //edge 2
      if(neigh_edge == 2)
      {
        for(unsigned int i = 0; i < N_Points1D; i++)
        {
          // for all quadrature points
          edgeData.xi1DNeigh[i] = 0;
          edgeData.eta1DNeigh[i] = (qf.get_point(i).x + 1) / 2;
        }
      }
      break;
    }
    case BFRefElements::BFUnitHexahedron:
    case BFRefElements::BFUnitTetrahedron:
    case BFRefElements::BFUnitLine:
      ErrThrow("ComputeAddStab: case is not supported yet");
      break;
  }
}

void ComputeAddStab(const TFESpace2D** fe_spaces, LocalAssembling2D& la, 
                    TSquareMatrix2D **sqmatrices, TSquareMatrix2D **A22)
{
  const unsigned int stride = 40;
  std::vector<double> aux(MaxN_QuadPoints_2D * stride);
  
  std::vector<double *> coeffPerQuadPoint;
  for(unsigned int i=0; i<MaxN_QuadPoints_2D; i++)
    coeffPerQuadPoint.push_back(&aux[i * stride]);
  
  auto fespace = fe_spaces[0];
  
  JointData edgeData{};
  
  const unsigned int max_n_base_funct_2d = fespace->get_max_n_local_dof();  
  JointRefData edgeRefData{max_n_base_funct_2d};

  edgeData.updateEdgeData(max_n_base_funct_2d);
  
  double jump_xyval[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];
  double jump_xderiv[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];
  double jump_yderiv[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];
            
  auto coll = fespace->GetCollection();
  coll->mark_all_cells();
  
  auto n_used_elements = fespace->GetN_UsedElements();
  // store used elements 
  std::vector<FE_type> UsedElements(n_used_elements);
  {
    std::vector<int> Used = std::vector<int>(N_FEs2D);
    auto n = fespace->GetN_UsedElements();
    auto use = fespace->GetUsedElements();
    for(auto j = 0; j < n; j++)
    {
      FE_type currElement = use[j];
      Used[currElement] = 1;
    }
    int N_used = 0;
    for(auto i = 0; i<N_FEs2D; i++)
      if(Used[i]) 
        N_used++;
    
    size_t j = 0;
    for(auto i = 0; i < N_FEs2D; i++)
    {
      if(Used[i])
      {
        UsedElements[j] = (FE_type) i;
        j++;
      }
    }
  }
  
  int N_Points1D;
  int N_Points;
  const TQuadFormula *quadFormula;
  bool secondDerivatives[] = {true};
  for(auto i = 0; i<n_used_elements; i++)
  {
    FE_type usedElement = UsedElements.at(i);
    const FiniteElement element{usedElement};
    quadFormula = QuadratureFormulaDatabase::qf_from_degree(
      2*element.GetBaseFunct()->GetPolynomialDegree(), BFRefElements::BFUnitLine);
    N_Points1D = quadFormula->GetN_QuadPoints();
    if(N_Points1D > edgeData.MaxN_QuadPoints_1D_loc)
    {
      ErrThrow("AddGadientJump: too many quad points", N_Points1D);
    }
    
    for(auto edge = 0; edge < 4; edge++)
    {
      edgeData.xyval_1D[edge].resize((unsigned long) N_Points1D);
      edgeData.xderiv_1D[edge].resize((unsigned long) N_Points1D);
      edgeData.yderiv_1D[edge].resize((unsigned long) N_Points1D);
    }
    
    BaseFunction_type baseFunct2D = element.GetBaseFunct_ID();
    const BaseFunctions *bf = element.GetBaseFunct();
    BFRefElements bf2Drefelements = bf->GetRefElement();
    
    int baseFunct_loc = 0;
    {
      int j;
      for(j = 0; j < n_used_elements; j++)
      {
        if ((int) baseFunct2D == (int) UsedElements[j])
          break;
      }
      baseFunct_loc = j;
    }
    //baseFunct_loc = element.GetBaseFunct_ID();
    get_values_ref(bf2Drefelements, *quadFormula, bf, edgeData, baseFunct_loc);
  }
  
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);
  
  auto N_Cells = coll->GetN_Cells();
  for(auto cellIdx = 0; cellIdx < N_Cells; cellIdx++)
  {
    auto cell = coll->GetCell(cellIdx);
    cell->SetClipBoard(cellIdx);
  }
  double integral = 0.; double integrand = 0.;
  double integral2 = 0.; double integrant2 = 0.;
  for (auto cellIdx = 0; cellIdx < N_Cells; cellIdx++)
  {
    auto cell = coll->GetCell(cellIdx);
    auto element = fespace->get_fe(cellIdx);    
    const BaseFunctions *bf = element.GetBaseFunct();
    int n_base_functs = element.GetN_DOF(); // N_
    auto DOF = fespace->GetGlobalDOF(cellIdx);
    
    std::vector<const FiniteElement*> used_fe(1, &element);
    ReferenceTransformation_type refTrans = FEDatabase::GetOrig(
      used_fe, coll, cell, secondDerivatives, qf_ref, qf_orig);
    N_Points = qf_orig.GetN_QuadPoints();
    
    // do we need fefunctions??????
    auto Coeffs = la.GetCoeffFct();
    
    const double *X = qf_orig.get_xi();
    const double *Y = qf_orig.get_eta();
    Coeffs(N_Points, X, Y, nullptr, &coeffPerQuadPoint[0]);
    
    const TQuadFormula *quadFormula = QuadratureFormulaDatabase::qf_from_degree(
      2 * bf->GetPolynomialDegree(), BFRefElements::BFUnitLine);
    N_Points1D = quadFormula->GetN_QuadPoints();
    
    int BaseFunct_loc = 0;
    {
      int j = 0;
      for(j = 0; j < n_used_elements; j++)
      {
        if( (int) bf->GetID() == (int) fespace->GetUsedElements()[j])
        {
          break;
        }
      }
      BaseFunct_loc = j;
    }
    
//     BaseFunct_loc = element.GetBaseFunct_ID();

    for(auto edgeIdx = 0; edgeIdx < cell->GetN_Edges(); edgeIdx++)
    {
      FEDatabase::GetOrigFromRef(refTrans, N_Points1D, 
        edgeData.xi1D[BaseFunct_loc][edgeIdx],
        edgeData.eta1D[BaseFunct_loc][edgeIdx],
        edgeData.XEdge1D[edgeIdx].data(),
        edgeData.YEdge1D[edgeIdx].data()
      );
      
      // values and derivatives on original cell
      for(int qp = 0; qp < N_Points1D; qp++)
      {
        FEDatabase::GetOrigValues(
            refTrans, edgeData.xi1D[BaseFunct_loc][edgeIdx][qp],
            edgeData.eta1D[BaseFunct_loc][edgeIdx][qp],
            bf, coll, (TGridCell *) cell,
            edgeData.xietaval_ref1D[BaseFunct_loc][edgeIdx][qp],
            edgeData.xideriv_ref1D[BaseFunct_loc][edgeIdx][qp],
            edgeData.etaderiv_ref1D[BaseFunct_loc][edgeIdx][qp],
            edgeData.xyval_ref1D[edgeIdx][qp].data(),
            edgeData.xderiv_ref1D[edgeIdx][qp].data(),
            edgeData.yderiv_ref1D[edgeIdx][qp].data()); 
      }
    }
    
    edgeRefData.updateQuadPointData(N_Points1D);
    auto N_Edges = cell->GetN_Edges();
    for(int edgeIdx = 0; edgeIdx < N_Edges; edgeIdx++)
    {
      auto neigh = cell->GetJoint(edgeIdx)->GetNeighbour(cell);
      
      Coeffs(N_Points, edgeData.XEdge1D[edgeIdx].data(), edgeData.YEdge1D[edgeIdx].data(), nullptr, &coeffPerQuadPoint[0]);
      
      if(neigh)
      {
        auto q = neigh->GetClipBoard();
        if(cellIdx < q)
        {
          auto elementNeigh = fespace->get_fe(q);
          auto bf_neigh = elementNeigh.GetBaseFunct();
          auto n_neigh = elementNeigh.GetN_DOF();
          auto DOF_neigh = fespace->GetGlobalDOF(q);
          auto BaseFunctNeigh = elementNeigh.GetBaseFunct_ID();
          std::vector<const FiniteElement*> used_elements(1, &elementNeigh);
          auto refTransNeigh = FEDatabase::GetOrig(used_elements, coll, neigh, 
                                                   secondDerivatives, qf_ref, qf_orig);
          
          int neigh_edge = 0;
          {
            while(neigh->GetJoint(neigh_edge)->GetNeighbour(neigh) != cell)
              neigh_edge++;
          }
          
          // reference of neighbour
          BFRefElements bf2DrefelementsNeigh = bf_neigh->GetRefElement();
          eta_to_edge_ref(bf2DrefelementsNeigh, neigh_edge, *quadFormula, edgeRefData);
          // compute gradients in reference cell of the neighbour
          for(int i = 0; i < N_Points1D; i++)
          {
            bf_neigh->GetDerivatives(MultiIndex2D::D00, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][i]);
            bf_neigh->GetDerivatives(MultiIndex2D::D10, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.xideriv_refNeigh1D[BaseFunctNeigh][i]);
            bf_neigh->GetDerivatives(MultiIndex2D::D01, edgeRefData.xi1DNeigh[i], edgeRefData.eta1DNeigh[i], edgeRefData.etaderiv_refNeigh1D[BaseFunctNeigh][i]);
          }
          
          FEDatabase::GetOrigFromRef(refTransNeigh, N_Points1D, 
              edgeRefData.xi1DNeigh.data(),
              edgeRefData.eta1DNeigh.data(),
              edgeRefData.X1DNeigh.data(),
              edgeRefData.Y1DNeigh.data() );
          
          if(edgeRefData.X1DNeigh[neigh_edge] == edgeData.XEdge1D[neigh_edge][0]
              && edgeRefData.Y1DNeigh[neigh_edge] == edgeData.YEdge1D[neigh_edge][0])
          {
            for(int qp = 0; qp < N_Points1D; qp++)
              FEDatabase::GetOrigValues(refTransNeigh, edgeRefData.xi1DNeigh[qp],
                        edgeRefData.eta1DNeigh[qp], bf_neigh, coll, (TGridCell *) neigh,
                        edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][qp],
                        edgeRefData.xideriv_refNeigh1D.data()[BaseFunctNeigh][qp],
                        edgeRefData.etaderiv_refNeigh1D.data()[BaseFunctNeigh][qp],
                        edgeRefData.xyval_refNeigh1D.data()[qp],
                        edgeRefData.xderiv_refNeigh1D.data()[qp],
                        edgeRefData.yderiv_refNeigh1D.data()[qp]);
          }
          else
          {
            for(int qp = 0; qp < N_Points1D; qp++)
              FEDatabase::GetOrigValues(refTransNeigh, edgeRefData.xi1DNeigh[qp],
                        edgeRefData.eta1DNeigh[qp], bf_neigh, coll, (TGridCell *) neigh,
                        edgeRefData.xietaval_refNeigh1D[BaseFunctNeigh][qp],
                        edgeRefData.xideriv_refNeigh1D.data()[BaseFunctNeigh][qp],
                        edgeRefData.etaderiv_refNeigh1D.data()[BaseFunctNeigh][qp],
                        edgeRefData.xyval_refNeigh1D.data()[qp],
                        edgeRefData.xderiv_refNeigh1D.data()[qp],
                        edgeRefData.yderiv_refNeigh1D.data()[qp]);
            // Output::print("Inverse the order of the quadrature oints on neighbour edge !");
          }
          
          for(int k = 0; k < n_base_functs; k++)
          {
            int index = 0;
            int l = 0;
            while( l < n_neigh && index == 0)
            {
              if(DOF[k] == DOF_neigh[l])
                index = 1;
              l++;
            }
            l = l-1;
            if(index ==1)
            {
              for(int qp = 0; qp < N_Points1D; qp++)
              {
                jump_xyval[qp][k] = edgeData.xyval_ref1D[edgeIdx][qp][k]-edgeRefData.xyval_refNeigh1D[qp][l];                                    
                jump_xderiv[qp][k]=edgeData.xderiv_ref1D[edgeIdx][qp][k]-edgeRefData.xderiv_refNeigh1D[qp][l];
                jump_yderiv[qp][k]=edgeData.yderiv_ref1D[edgeIdx][qp][k]-edgeRefData.yderiv_refNeigh1D[qp][l];
               /* {
                  Output::print("x: ", edgeData.xyval_ref1D[edgeIdx][qp][k],
                                  " xn: ", edgeRefData.xyval_refNeigh1D[qp][l],
                                  " j= ", jump_xyval[qp][k],
                                  " bf: ", DOF[k], " c: ", cellIdx,
                                  " e: ", edgeIdx,
                                  " qp: ", qp);
                    Output::print("xd: ", edgeData.xderiv_ref1D[edgeIdx][qp][k],
                                  " xdn: ", edgeRefData.xderiv_refNeigh1D[qp][l],
                                  " jump= ", jump_xderiv[qp][k],
                                  " bf: ", DOF[k], " c: ", cellIdx,
                                  " e: ", edgeIdx, " qp: ", qp);
                    Output::print("yd: ", edgeData.yderiv_ref1D[edgeIdx][qp][k],
                                  " ydn: ", edgeRefData.yderiv_refNeigh1D[qp][l],
                                  " j= ", jump_yderiv[qp][k],
                                  " bf: ", DOF[k], " c: ", cellIdx,
                                  " e: ", edgeIdx, " qp: ", qp);
                }*/
              }
            }//endif index == 1
            if(index == 0)
            {
              for(int qp = 0; qp < N_Points1D; qp++)
              {
                jump_xyval[qp][k]  =  edgeData.xyval_ref1D[edgeIdx][qp][k];
                jump_xderiv[qp][k] =  edgeData.xderiv_ref1D[edgeIdx][qp][k];
                jump_yderiv[qp][k] =  edgeData.yderiv_ref1D[edgeIdx][qp][k];
                
                /*{
                  Output::print("x: ", edgeData.xyval_ref1D[edgeIdx][qp][k],
                                  " j= ", jump_xyval[qp][k],
                                  " bf: ", DOF[k], " c: ", cellIdx,
                                  " e: ", edgeIdx,
                                  " qp: ", qp);
                    Output::print("xd: ", edgeData.xderiv_ref1D[edgeIdx][qp][k],
                                  " jump= ", jump_xderiv[qp][k],
                                  " bf: ", DOF[k], " c: ", cellIdx,
                                  " e: ", edgeIdx, " qp: ", qp);
                    Output::print("yd: ", edgeData.yderiv_ref1D[edgeIdx][qp][k],
                                  " j= ", jump_yderiv[qp][k],
                                  " bf: ", DOF[k], " c: ", cellIdx,
                                  " e: ", edgeIdx, " qp: ", qp);
                }*/
              }
            }
          }//endfor k < n_base_functs
          // basis functions of neighbour cell q which are not in local FE-space of cellIdx
          for(int l = 0; l < n_neigh; l++)
          {
            int index = 0;
            int k = 0;
            while (k < n_base_functs && index == 0)
            {
              if(DOF_neigh[l] == DOF[k])
                index++;
              k++;
            }
            k = k-1;
            // if l of neighbour q is not in local fe-space of cellIdx
            if(index == 0)
            {
              for(int qp = 0; qp < N_Points1D; qp++)
              {
                jump_xyval[qp][l+n_base_functs] = - edgeRefData.xyval_refNeigh1D[qp][l];
                jump_xderiv[qp][l+n_base_functs] = - edgeRefData.xderiv_refNeigh1D[qp][l];
                jump_yderiv[qp][l+n_base_functs] = - edgeRefData.yderiv_refNeigh1D[qp][l];
                
                /*{
                  Output::print(
                                  " xn: ", edgeRefData.xyval_refNeigh1D[qp][l],
                                  " j= ", jump_xyval[qp][k],
                                  " bf: ", DOF[k], " c: ", cellIdx,
                                  " e: ", edgeIdx,
                                  " qp: ", qp);
                    Output::print(
                                  " xdn: ", edgeRefData.xderiv_refNeigh1D[qp][l],
                                  " jump= ", jump_xderiv[qp][k],
                                  " bf: ", DOF[k], " c: ", cellIdx,
                                  " e: ", edgeIdx, " qp: ", qp);
                    Output::print(
                                  " ydn: ", edgeRefData.yderiv_refNeigh1D[qp][l],
                                  " j= ", jump_yderiv[qp][k],
                                  " bf: ", DOF[k], " c: ", cellIdx,
                                  " e: ", edgeIdx, " qp: ", qp);
                }*/
              }
            }
          }//endfor l < n_neigh
          
          // compute integrals with jump derivatives
          double x0, y0, x1, y1;
          cell->GetVertex(edgeIdx)->GetCoords(x0, y0);
          cell->GetVertex((edgeIdx+1) % N_Edges)->GetCoords(x1,y1);
          double nx = (y1-y0);
          double ny = (x0-x1);
          auto hE = std::sqrt(nx * nx + ny * ny);
          nx /= hE;
          ny /= hE;
          auto delta1 = TDatabase::ParamDB->FACE_SIGMA * hE * hE;
          auto sqm1Entries = sqmatrices[0]->GetEntries();
          auto rowPtr1 = sqmatrices[0]->get_row_ptr();
          auto colInd1 = sqmatrices[0]->get_vector_columns();
          auto fespace1 = sqmatrices[0]->GetFESpace2D();
          auto active = fespace1->get_n_active();
          
          auto a22entries = A22[0]->GetEntries();
          
          /*for(auto j=0;j<N_Points1D; j++)
          {
            Output::print(" b1: ", coeffPerQuadPoint[j][5], " b2: ",
                          coeffPerQuadPoint[j][6]);
          }*/
          //NOTE: test function from cell and ansatz from the cell
          for(int bfIdx = 0; bfIdx < n_base_functs; bfIdx++)
          {
            auto doff_bfIdx = DOF[bfIdx];
            if(doff_bfIdx  >= active)
              continue;
            
            for(int ansatz_dof = 0; ansatz_dof < n_base_functs; ansatz_dof++)
            {
              auto dof_ansatz = DOF[ansatz_dof];
              integral = 0.;
              integral2 = 0.;
              for(auto qp = 0; qp < N_Points1D; qp++)
              {
                switch (TDatabase::ParamDB->CIP_TYPE)
                {
                  case 0:
                    integrand = delta1 * std::abs(nx*coeffPerQuadPoint[qp][5] + ny*coeffPerQuadPoint[qp][6]) *
                                    (nx * jump_xderiv[qp][ansatz_dof] + ny * jump_yderiv[qp][ansatz_dof] ) *
                                    (nx * jump_xderiv[qp][bfIdx] + ny * jump_yderiv[qp][bfIdx]);
                    integrant2 = integrand;
                    break;
                  case 1:
                    integrand = delta1 * std::abs(nx*coeffPerQuadPoint[qp][5] + ny*coeffPerQuadPoint[qp][6]) *
                                    (coeffPerQuadPoint[qp][5] * ny * jump_xderiv[qp][ansatz_dof] 
                                   + coeffPerQuadPoint[qp][6] * ny * jump_yderiv[qp][ansatz_dof] ) *
                                    (coeffPerQuadPoint[qp][5] * ny * jump_xderiv[qp][bfIdx] 
                                   + coeffPerQuadPoint[qp][6] * ny * jump_yderiv[qp][bfIdx]);
                                  
                    integrant2 = -delta1 * std::abs(nx*coeffPerQuadPoint[qp][5] + ny*coeffPerQuadPoint[qp][6]) *
                                (coeffPerQuadPoint[qp][5] * nx * jump_xderiv[qp][ansatz_dof] 
                               + coeffPerQuadPoint[qp][6] * nx * jump_yderiv[qp][ansatz_dof] ) *
                                (coeffPerQuadPoint[qp][5] * nx * jump_xderiv[qp][bfIdx] 
                               + coeffPerQuadPoint[qp][6] * nx * jump_yderiv[qp][bfIdx]);
                    break;
                }
                double weights = quadFormula->get_weight(qp) * hE/2.;
                integral += weights * integrand;
                integral2 += weights * integrant2;
              }
              // update the matrix 
              int found = 0;
              for(auto row = rowPtr1[doff_bfIdx]; row < rowPtr1[doff_bfIdx+1]; row++)
              {
                if(colInd1[row] == dof_ansatz)
                {
                  sqm1Entries[row] += integral;
                  a22entries[row]  += integral2;
                  found = 1;
                  break;
                }
              }
              if(!found)
              {
                ErrThrow("error in the test and ansatz function from the cell");
              }
            }//endfor ansatz_dof
          }//endfor doff_bfIdx test
          
          //NOTE: test function from the cell and ansatz from the neighbour
          for(int testIdx = 0; testIdx < n_base_functs; testIdx++)
          {
            auto doff_test = DOF[testIdx];
            if(doff_test >= active)
              continue;
            for(int anzIdx = 0; anzIdx < n_neigh; anzIdx++)
            {
              auto doff_ansatz = DOF_neigh[anzIdx];
              int index = 0;
              int l = 0;
              while(l < n_base_functs && index == 0)
              {
                if(doff_ansatz == DOF[l])
                  index = 1;
                l++;
              }
              if(index == 1)
                continue;
              integral = 0.;
              integral2 = 0.;
              for(auto qp = 0; qp < N_Points1D; qp++)
              {
                switch(TDatabase::ParamDB->CIP_TYPE)
                {
                  case 0:
                    integrand = delta1 * std::abs(nx*coeffPerQuadPoint[qp][5] + ny*coeffPerQuadPoint[qp][6]) *
                            (jump_xderiv[qp][anzIdx + n_base_functs] * nx
                            + jump_yderiv[qp][anzIdx + n_base_functs] * ny) *
                            (jump_xderiv[qp][testIdx]*nx + jump_yderiv[qp][testIdx]*ny);                    
                    integrant2 = integrand;
                    break;
                  case 1:
                    integrand = -delta1 * std::abs(nx*coeffPerQuadPoint[qp][5] + ny*coeffPerQuadPoint[qp][6]) *
                            (coeffPerQuadPoint[qp][5] * jump_xderiv[qp][anzIdx + n_base_functs] * ny
                            + coeffPerQuadPoint[qp][6] * jump_yderiv[qp][anzIdx + n_base_functs] * ny) *
                            (coeffPerQuadPoint[qp][5] * jump_xderiv[qp][testIdx]*ny 
                            + coeffPerQuadPoint[qp][6] * jump_yderiv[qp][testIdx]*ny);
                    
                    integrant2 = delta1 * std::abs(nx*coeffPerQuadPoint[qp][5] + ny*coeffPerQuadPoint[qp][6]) *
                            (coeffPerQuadPoint[qp][5] * jump_xderiv[qp][anzIdx + n_base_functs] * nx
                            + coeffPerQuadPoint[qp][6] * jump_yderiv[qp][anzIdx + n_base_functs] * nx) *
                            (coeffPerQuadPoint[qp][5] * jump_xderiv[qp][testIdx]*nx 
                            + coeffPerQuadPoint[qp][6] * jump_yderiv[qp][testIdx]*nx);
                    break;
                }
                double weights = quadFormula->get_weight(qp) * hE/2.0;
                integral += weights * integrand;
                integral2 += weights * integrant2;
              }
              // update matrix
              int found = 0;
              for(auto row = rowPtr1[doff_test]; row < rowPtr1[doff_test+1]; row++)
              {
                if(colInd1[row] == doff_ansatz)
                {
                  sqm1Entries[row] += integral;
                  a22entries[row] += integral2;
                  found = 1;
                  break;
                }
              }//endfor row 
              if(!found)
              {
                ErrThrow("error in the test function cell, ansatz function neighbour");
              }
            }//endfor anzIdx < n_neigh
          }//endfor testIdx < n_base_functs
          
          //NOTE: test neighbour anzats cell 
          for(int testIdx = 0; testIdx < n_neigh; testIdx++)
          {
            auto doff_test = DOF_neigh[testIdx];
            if(doff_test >= active)
              continue;
            int index = 0;
            int l = 0;
            while(l < n_base_functs && index == 0)
            {
              if(doff_test == DOF[l])
                index = 1;
              l++;
            }
            if(index == 1)
              continue;            
            for(int anzIdx = 0; anzIdx < n_base_functs; anzIdx++)
            {
              auto doff_ansatz = DOF[anzIdx];
              integral = 0;
              integral2 = 0.;
              for(auto qp = 0; qp < N_Points1D; qp++)
              {
                switch (TDatabase::ParamDB->CIP_TYPE)
                {
                  case 0:
                    integrand = delta1 * std::abs(nx*coeffPerQuadPoint[qp][5] + ny*coeffPerQuadPoint[qp][6]) * 
                     (jump_xderiv[qp][anzIdx] * nx + jump_yderiv[qp][anzIdx] * ny) *
                     (jump_xderiv[qp][testIdx + n_base_functs]*nx + jump_yderiv[qp][testIdx + n_base_functs]*ny);
                    integrant2 = integrand;
                    break;
                  case 1:
                    integrand = delta1 * std::abs(nx*coeffPerQuadPoint[qp][5] + ny*coeffPerQuadPoint[qp][6]) * 
                            (coeffPerQuadPoint[qp][5] * jump_xderiv[qp][anzIdx] * ny
                            + coeffPerQuadPoint[qp][6] * jump_yderiv[qp][anzIdx] * ny) *
                            (coeffPerQuadPoint[qp][5] * jump_xderiv[qp][testIdx + n_base_functs]*ny 
                            + coeffPerQuadPoint[qp][6] * jump_yderiv[qp][testIdx + n_base_functs]*ny);
                     integrant2 = -delta1 * std::abs(nx*coeffPerQuadPoint[qp][5] + ny*coeffPerQuadPoint[qp][6]) *
                            (coeffPerQuadPoint[qp][5] * jump_xderiv[qp][anzIdx] * nx
                            + coeffPerQuadPoint[qp][6] * jump_yderiv[qp][anzIdx] * nx) *
                            (coeffPerQuadPoint[qp][5] * jump_xderiv[qp][testIdx + n_base_functs]*nx 
                            + coeffPerQuadPoint[qp][6] * jump_yderiv[qp][testIdx + n_base_functs]*nx);
                    break;
                }
                double weights = quadFormula->get_weight(qp) * hE/2.0;
                integral += weights * integrand;
                integral2 += weights * integrant2;
              }
              int found = 0;
              for(auto row = rowPtr1[doff_test]; row < rowPtr1[doff_test+1]; row++)
              {
                if(colInd1[row] == doff_ansatz)
                {
                  sqm1Entries[row] += integral;
                  a22entries[row] += integral2;
                  found = 1;
                  break;
                }
              }
              if(!found)
              {
                ErrThrow("error in the test function neighbour, ansatz function cell");
              }
            }//endfor anzats
          }//endfor test
          
          //NOTE: test and anzats from neighbour
          for(auto testIdx = 0; testIdx < n_neigh; testIdx++)
          {
            auto doff_test = DOF_neigh[testIdx];
            if(doff_test >= active)
              continue;
            int index = 0;
            int l = 0;
            while(l < n_base_functs && index == 0)
            {
              if(doff_test == DOF[l])
                index = 1;
              l++;
            }
            if(index == 1)
              continue;
            for(auto anzIdx = 0; anzIdx < n_neigh; anzIdx++)
            {
              auto doff_ansatz = DOF_neigh[anzIdx];
              index = 0;
              l = 0;
              while(l < n_base_functs && index == 0)
              {
                if(doff_ansatz == DOF[l])
                  index = 1;
                l++;
              }
              if(index==1)
                continue;
              
              integral = 0.;
              integral2 = 0.;
              for(auto qp = 0; qp < N_Points1D; qp++)
              {
                switch (TDatabase::ParamDB->CIP_TYPE)
                {
                  case 0:
                    integrand = delta1 * std::abs(nx*coeffPerQuadPoint[qp][5] + ny*coeffPerQuadPoint[qp][6]) *
                        (jump_xderiv[qp][anzIdx + n_base_functs] * nx + jump_yderiv[qp][anzIdx + n_base_functs] * ny) *
                        (jump_xderiv[qp][testIdx + n_base_functs]*nx + jump_yderiv[qp][testIdx + n_base_functs]*ny);
                    
                    integrant2 = integrand;
                    break;
                  case 1:
                    integrand = delta1 * std::abs(nx*coeffPerQuadPoint[qp][5] + ny*coeffPerQuadPoint[qp][6]) *
                        (coeffPerQuadPoint[qp][5] * jump_xderiv[qp][anzIdx + n_base_functs] * ny
                        + coeffPerQuadPoint[qp][6] * jump_yderiv[qp][anzIdx + n_base_functs] * ny) *
                        (coeffPerQuadPoint[qp][5] * jump_xderiv[qp][testIdx + n_base_functs]*ny 
                        + coeffPerQuadPoint[qp][6] * jump_yderiv[qp][testIdx + n_base_functs]*ny);
                   
                   integrant2 = -delta1 * std::abs(nx*coeffPerQuadPoint[qp][5] + ny*coeffPerQuadPoint[qp][6]) *
                        (coeffPerQuadPoint[qp][5] * jump_xderiv[qp][anzIdx + n_base_functs] * nx
                      + coeffPerQuadPoint[qp][6] * jump_yderiv[qp][anzIdx + n_base_functs] * nx) *
                        (coeffPerQuadPoint[qp][5] * jump_xderiv[qp][testIdx + n_base_functs]*nx 
                      + coeffPerQuadPoint[qp][6] * jump_yderiv[qp][testIdx + n_base_functs]*nx);
                    break;
                }
                double weights = quadFormula->get_weight(qp) * hE/2.0;
                integral += weights * integrand;                
                integral2 += weights * integrant2;
              }
              int found = 0;
              for(auto row = rowPtr1[doff_test]; row < rowPtr1[doff_test+1]; row++)
              {
                if(colInd1[row] == doff_ansatz)
                {
                  sqm1Entries[row] += integral;
                  a22entries[row] += integral2;
                  found = 1;
                  break;
                }
              }
              if(!found)
              {
                ErrThrow("error in the test-ansatz from neighbour");
              }
            }
          }
        }
        //Output::print("cell: ", cellIdx, " there is a neighbour ");
      }//endif neigh
      else
      {
       // Output::print("cell: ", cellIdx," there is No Neighbour");
      }      
    }//endfor edgeIdx < N_Edges    
  }//endfor cellIdx < N_Cells
  
//   auto rowPtr1 = sqmatrices[0]->get_row_ptr();
//   auto sqm1Entries = sqmatrices[0]->GetEntries();
//   auto colInd1 = sqmatrices[0]->get_vector_columns();
//   auto row = sqmatrices[0]->get_n_rows();
//   
//   for(int i=0;i<row;i++)
//   {
//     int end=rowPtr1[i+1];
//     for(int j=rowPtr1[i];j<end;j++)
//     {
//       // cout << j << endl;
//       cout << "Matrix: " << setw(5) << i << setw(5) << colInd1[j] << "   ";
//       cout << setw(10) << sqm1Entries[j] << endl;
//     }
//   }
//       cout << endl;  
}
