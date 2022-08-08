#include "FEDatabase.h"
#include "AllRefTrans.h"
#include "AllRefTrans3D.h"
#include "Database.h"
#include "FiniteElement.h"
#include "GridCell.h"
#include "JointEqN.h"
#include "LinAlg.h"
#include "MooNMD_Io.h"
#include "NodalFunctional.h"
#include "QuadratureFormulaDatabase.h"
#include <algorithm>
#include <map>
#include <vector>

namespace FEDatabase
{
  /** values of FE functions and their derivatives on the
   * corresponding reference element
   */
  std::map<std::tuple<BaseFunction_type, QuadratureFormula_type, MultiIndex2D>,
           double **> RefElementValues2D;
  std::map<std::tuple<BaseFunction_type, QuadratureFormula_type, MultiIndex3D>,
           double **> RefElementValues3D;

  /** derivatives of FE base functions on the edges of the corresponding
   * reference element
   */
  std::map<std::tuple<BaseFunction_type, QuadratureFormula_type, int,
                      MultiIndex2D>, double **>
      JointDerivatives2D;
  std::map<std::tuple<BaseFunction_type, QuadratureFormula_type, int,
                      MultiIndex3D>, double **>
      JointDerivatives3D;

  /** values of FE functions and their derivatives on the
   * current element
   */
  double **OrigElementValues2D[N_BaseFuncts2D][N_MultiIndices2D] = {{nullptr}};
  double **OrigElementValues3D[N_BaseFuncts3D][N_MultiIndices3D] = {{nullptr}};

  /** reference transformations */
  std::array<TRefTrans2D*, N_RefTrans2D+N_RefTrans1D> ReferenceTrans2D = {};
  
  std::array<TRefTrans3D*, N_RefTrans3D+N_RefTrans2D+N_RefTrans1D>
    ReferenceTrans3D = {};

  /** prolongation matrix storage */
  double *ProlongationMatrix2D[N_BaseFuncts2D][N_REFDESC][N_BaseFuncts2D]
                              [MAXN_CHILDREN] = {{{{ nullptr }}}};
  double *ProlongationMatrix3D[N_BaseFuncts3D][N_REFDESC][N_BaseFuncts3D]
                              [MAXN_CHILDREN] = {{{{ nullptr }}}};

  /** function restriction matrix storage */
  double *RestrictionMatrix2D[N_BaseFuncts2D][N_REFDESC][N_BaseFuncts2D]
                             [MAXN_CHILDREN] = {{{{ nullptr }}}};
  double *RestrictionMatrix3D[N_BaseFuncts3D][N_REFDESC][N_BaseFuncts3D]
                             [MAXN_CHILDREN] = {{{{ nullptr }}}};

//==============================================================================
void destroy()
{
  for(auto it : RefElementValues2D)
  {
    delete [] it.second[0];
    delete [] it.second;
  }
  RefElementValues2D.clear();
  for(auto it : JointDerivatives2D)
  {
    delete [] it.second[0];
    delete [] it.second;
  }
  JointDerivatives2D.clear();
  for(int i = 0; i < N_BaseFuncts2D; ++i)
  {
    for(int j = 0; j < N_MultiIndices2D; ++j)
    {
      if(OrigElementValues2D[i][j] != nullptr)
      {
        delete [] OrigElementValues2D[i][j][0];
        delete [] OrigElementValues2D[i][j];
        OrigElementValues2D[i][j] = nullptr;
      }
    }
  }
  for(int i = 0; i < N_RefTrans2D+N_RefTrans1D; ++i)
  {
    delete ReferenceTrans2D[i];
    ReferenceTrans2D[i] = nullptr;
  }
  for(int i = 0; i < N_BaseFuncts2D; ++i)
  {
    for(int j = 0; j < N_REFDESC; ++j)
    {
      for(int k = 0; k < N_BaseFuncts2D; ++k)
      {
        for(int l = 0; l < MAXN_CHILDREN; ++l)
        {
          delete [] ProlongationMatrix2D[i][j][k][l];
          ProlongationMatrix2D[i][j][k][l] = nullptr;
          delete [] RestrictionMatrix2D[i][j][k][l];
          RestrictionMatrix2D[i][j][k][l] = nullptr;
        }
      }
    }
  }
  // 3D 
  for(auto it : RefElementValues3D)
  {
    delete [] it.second[0];
    delete [] it.second;
  }
  RefElementValues3D.clear();
  for(auto it : JointDerivatives3D)
  {
    delete [] it.second[0];
    delete [] it.second;
  }
  JointDerivatives3D.clear();
  for(int i = 0; i < N_BaseFuncts3D; ++i)
  {
    for(int j = 0; j < N_MultiIndices3D; ++j)
    {
      if(OrigElementValues3D[i][j] != nullptr)
      {
        delete [] OrigElementValues3D[i][j][0];
        delete [] OrigElementValues3D[i][j];
        OrigElementValues3D[i][j] = nullptr;
      }
    }
  }
  for(int i = 0; i < N_RefTrans3D+N_RefTrans2D+N_RefTrans1D; ++i)
  {
    delete ReferenceTrans3D[i];
    ReferenceTrans3D[i] = nullptr;
  }
  for(int i = 0; i < N_BaseFuncts3D; ++i)
  {
    for(int j = 0; j < N_REFDESC; ++j)
    {
      for(int k = 0; k < N_BaseFuncts3D; ++k)
      {
        for(int l = 0; l < MAXN_CHILDREN; ++l)
        {
          delete [] ProlongationMatrix3D[i][j][k][l];
          ProlongationMatrix3D[i][j][k][l] = nullptr;
          delete [] RestrictionMatrix3D[i][j][k][l];
          RestrictionMatrix3D[i][j][k][l] = nullptr;
        }
      }
    }
  }
}

//==============================================================================
TGridCell * GenerateRefElement2D(BFRefElements ref_element_id)
{
  TGridCell *Cell = nullptr;
  TVertex *v[4];
  TJointEqN *b[4];
  switch(ref_element_id)
  {
    case BFRefElements::BFUnitSquare:
      Cell = new TGridCell(TDatabase::RefDescDB[Rectangle], 0);
      v[0] = new TVertex(-1.0, -1.0);
      v[1] = new TVertex( 1.0, -1.0);
      v[2] = new TVertex( 1.0,  1.0);
      v[3] = new TVertex(-1.0,  1.0);
    
      Cell->SetClipBoard(Refinement);

      for(int i=0;i<4;i++)
      {
        b[i] = new TJointEqN(Cell);
        Cell->SetVertex(i, v[i]);
        Cell->SetJoint(i, b[i]);
      } // endfor i
      break;
    case BFRefElements::BFUnitTriangle:
      Cell = new TGridCell(TDatabase::RefDescDB[Triangle], 0);
      v[0] = new TVertex(0.0, 0.0);
      v[1] = new TVertex(1.0, 0.0);
      v[2] = new TVertex(0.0, 1.0);
    
      Cell->SetClipBoard(Refinement);

      for(int i=0;i<3;i++)
      {
        b[i] = new TJointEqN(Cell);
        Cell->SetVertex(i, v[i]);
        Cell->SetJoint(i, b[i]);
      } // endfor i
      break;
    default:
      cout << "unknown reference element!" << endl;
  }
  return Cell;
}

//==============================================================================
TGridCell* GenerateRefElement3D(BFRefElements ref_element_id)
{
  TGridCell *Cell =nullptr;
  TVertex *v[8];
  TJointEqN *b[6];
  int i;

  switch(ref_element_id)
  {
    case BFRefElements::BFUnitHexahedron:
      Cell = new TGridCell(TDatabase::RefDescDB[Hexahedron], 0);
      v[0] = new TVertex(-1.0, -1.0, -1.0);
      v[1] = new TVertex( 1.0, -1.0, -1.0);
      v[2] = new TVertex( 1.0,  1.0, -1.0);
      v[3] = new TVertex(-1.0,  1.0, -1.0);
      v[4] = new TVertex(-1.0, -1.0,  1.0);
      v[5] = new TVertex( 1.0, -1.0,  1.0);
      v[6] = new TVertex( 1.0,  1.0,  1.0);
      v[7] = new TVertex(-1.0,  1.0,  1.0);
    
      Cell->SetClipBoard(Refinement);

      for(i=0;i<8;i++)
        Cell->SetVertex(i, v[i]);

      for(i=0;i<6;i++)
      {
        b[i] = new TJointEqN(Cell);
        Cell->SetJoint(i, b[i]);
      }

    break;

    case BFRefElements::BFUnitTetrahedron:
      Cell = new TGridCell(TDatabase::RefDescDB[Tetrahedron], 0);
      v[0] = new TVertex(0.0, 0.0, 0.0);
      v[1] = new TVertex(1.0, 0.0, 0.0);
      v[2] = new TVertex(0.0, 1.0, 0.0);
      v[3] = new TVertex(0.0, 0.0, 1.0);
    
      Cell->SetClipBoard(Refinement);

      for(i=0;i<4;i++)
      {
        b[i] = new TJointEqN(Cell);
        Cell->SetVertex(i, v[i]);
        Cell->SetJoint(i, b[i]);
      } // endfor i
    break;

    default:
      ErrThrow("unknown Reference element in 3D ", ref_element_id);
      break;
  }

  return Cell;
}

//==============================================================================
double *GetProlongationMatrix2D(FE_type parent, Refinements refine,
                                FE_type child, int childnumber)
{
  double *ret, *ret2;
  int j,k,l;
  int N_Coarse, N_Points, N_Children;
  const double *xi, *eta;
  double X[MaxN_PointsForNodal2D], Y[MaxN_PointsForNodal2D];
  double AllPointValues[MaxN_PointsForNodal2D][MaxN_BaseFunctions2D*2];
  double PointValues[MaxN_PointsForNodal2D*2];
  const TRefDesc *RefDesc;
  BaseFunction_type Coarse, Fine;
  const TGridCell *cell;
  BFRefElements RefElement;
  ReferenceTransformation_type F_K = ReferenceTransformation_type::TriaAffin; 
  //avoid uninit warning
  TRefTrans2D *rt;

  const FiniteElement CoarseElement(parent);
  Coarse = CoarseElement.GetBaseFunct_ID();
  const FiniteElement FineElement(child);
  Fine = FineElement.GetBaseFunct_ID();
  int basevect_dim = FineElement.GetBaseFunct()->GetBaseVectDim();

  ret = ProlongationMatrix2D[Coarse][refine][Fine][childnumber]; 

  if(ret == nullptr)
  {
    // cerr << "ret == nullptr" << endl;
    // prolongation matrix was not generated yet

    auto BaseFunctions = CoarseElement.GetBaseFunct();
    N_Coarse = BaseFunctions->GetDimension();

    ret2 = new double[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];

    TGridCell *RefCell = GenerateRefElement2D(BaseFunctions->GetRefElement());
    
    if(refine == NoRef)
    {
      ret = new double[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
      auto nf = FineElement.GetNodalFunctional();
      nf->GetPointsForAll(N_Points, xi, eta);

      RefElement = BaseFunctions->GetRefElement();

      switch(RefElement)
      {
        case BFRefElements::BFUnitSquare:
          rt = GetRefTrans2D(ReferenceTransformation_type::QuadAffin);
          ((TQuadAffin *)rt)->SetCell(RefCell);
          F_K = ReferenceTransformation_type::QuadAffin;
          break;
        case BFRefElements::BFUnitTriangle:
          rt = GetRefTrans2D(ReferenceTransformation_type::TriaAffin);
          ((TTriaAffin *)rt)->SetCell(RefCell);
          F_K = ReferenceTransformation_type::TriaAffin;
          break;
        default:
          F_K = ReferenceTransformation_type::TriaAffin;
          break;
      }
      GetOrigFromRef(F_K, N_Points, xi, eta, X, Y);

      for(k=0;k<N_Points;k++)
        BaseFunctions->GetDerivatives(MultiIndex2D::D00, X[k], Y[k], 
                                      AllPointValues[k]);

      for(k=0;k<N_Coarse;k++)
      {
        for(l=0;l<N_Points;l++)
          PointValues[l] = AllPointValues[l][k];
        if(basevect_dim != 1)
          PointValues[l+N_Points] = AllPointValues[l][k+N_Points];

        nf->GetAllFunctionals(nullptr, nullptr, PointValues,
                              ret2+k*MaxN_BaseFunctions2D);
      }

      for(k=0;k<MaxN_BaseFunctions2D;k++)
        for(l=0;l<MaxN_BaseFunctions2D;l++)
          ret[k*MaxN_BaseFunctions2D+l] = 
                        ret2[l*MaxN_BaseFunctions2D+k];

      ProlongationMatrix2D[Coarse][refine][Fine][childnumber] = ret;
    }
    else
    {
      RefDesc = TDatabase::RefDescDB[refine+N_SHAPES];
      N_Children = RefDesc->GetN_Children();
      
      RefCell->SetRefDesc(RefDesc);
      RefCell->Refine(1);

      for(j=0;j<N_Children;j++)
      {
        ret = new double[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  
        cell = (const TGridCell *)((const TGridCell *)RefCell)->GetChild(j);
        //int N_Fine = FineElement->GetBaseFunct2D()->GetDimension();
  
        auto nf = FineElement.GetNodalFunctional();
        nf->GetPointsForAll(N_Points, xi, eta);

        switch(cell->GetType())
        {
          case Quadrangle:
          case Parallelogram:
          case Rectangle:
            rt = GetRefTrans2D(ReferenceTransformation_type::QuadAffin);
            ((TQuadAffin *)rt)->SetCell(cell);
            F_K = ReferenceTransformation_type::QuadAffin;
          break;
          case Triangle:
            rt = GetRefTrans2D(ReferenceTransformation_type::TriaAffin);
            ((TTriaAffin *)rt)->SetCell(cell);
            F_K = ReferenceTransformation_type::TriaAffin;
          break;
         default:
           Output::print("Not a 2D cell type");
         break;
        }
        GetOrigFromRef(F_K ,N_Points, xi, eta, X, Y);
  
        for(k=0;k<N_Points;k++)
          BaseFunctions->GetDerivatives(MultiIndex2D::D00, X[k], Y[k], 
                        AllPointValues[k]);
  
        for(k=0;k<N_Coarse;k++)
        {
          for(l=0;l<N_Points;l++)
            PointValues[l] = AllPointValues[l][k];
          if(basevect_dim != 1)
            for(l=0;l<N_Points;l++)
              PointValues[l+N_Points] = AllPointValues[l][k+N_Points];
  
          nf->GetAllFunctionals(nullptr, nullptr, PointValues,
                                ret2+k*MaxN_BaseFunctions2D);
        }

        for(k=0;k<MaxN_BaseFunctions2D;k++)
          for(l=0;l<MaxN_BaseFunctions2D;l++)
            ret[k*MaxN_BaseFunctions2D+l] = ret2[l*MaxN_BaseFunctions2D+k];

        ProlongationMatrix2D[Coarse][refine][Fine][j] = ret;
      } // endfor j
    }

    ret = ProlongationMatrix2D[Coarse][refine][Fine][childnumber]; 

    RefCell->Derefine();
    for(int c = 0; c < RefCell->GetN_Vertices(); ++c)
      delete RefCell->GetVertex(c);
    delete (const TGridCell *)RefCell;
    delete [] ret2;
  }

  return ret;
}

//==============================================================================
#ifdef __3D__
double *GetProlongationMatrix3D(FE_type parent, Refinements refine,
                                FE_type child, int childnumber)
{ 
  double *ret, *ret2;
  int j,k,l;
  int N_Coarse, N_Points, N_Children; // int N_Fine;
  const double *xi, *eta, *zeta;
  double X[MaxN_PointsForNodal3D], Y[MaxN_PointsForNodal3D];
  double Z[MaxN_PointsForNodal3D];
  double AllPointValues[MaxN_PointsForNodal3D][MaxN_BaseFunctions3D];
  double PointValues[MaxN_PointsForNodal3D];
  const TRefDesc *RefDesc;
  TGridCell *RefCell;
  const TGridCell * cell;
  BFRefElements RefElement;
  ReferenceTransformation_type F_K = ReferenceTransformation_type::TetraAffin; 
//avoid uninit warning
  TRefTrans3D *rt;

  const FiniteElement CoarseElement(parent);
  int Coarse = static_cast<int>(CoarseElement.GetBaseFunct_ID())
               - N_BaseFuncts1D - N_BaseFuncts2D;
  const FiniteElement FineElement(child);
  int Fine = static_cast<int>(FineElement.GetBaseFunct_ID())
             - N_BaseFuncts1D - N_BaseFuncts2D;

  ret = ProlongationMatrix3D[Coarse][refine][Fine][childnumber]; 

  if(ret == nullptr)
  {
    // cerr << "ret == nullptr" << endl;
    // prolongation matrix was not generated yet

    auto BaseFunctions = CoarseElement.GetBaseFunct();
    N_Coarse = BaseFunctions->GetDimension();

    if(refine == NoRef)
    {
      ret = new double[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
      ret2 = new double[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];

      RefCell = GenerateRefElement3D(BaseFunctions->GetRefElement());
      auto nf = FineElement.GetNodalFunctional();
      nf->GetPointsForAll(N_Points, xi, eta, zeta);

      RefElement = BaseFunctions->GetRefElement();

      switch(RefElement)
      {
        case BFRefElements::BFUnitHexahedron:
          rt = GetRefTrans3D(ReferenceTransformation_type::HexaAffin);
          ((THexaAffin *)rt)->SetCell(RefCell);
          F_K = ReferenceTransformation_type::HexaAffin;
          break;
        case BFRefElements::BFUnitTetrahedron:
          rt = GetRefTrans3D(ReferenceTransformation_type::TetraAffin);
          ((TTetraAffin *)rt)->SetCell(RefCell);
          F_K = ReferenceTransformation_type::TetraAffin;
          break;
        default:
          F_K = ReferenceTransformation_type::TetraAffin;
      }
      GetOrigFromRef(F_K, N_Points, xi, eta, zeta, X, Y, Z);

      for(k=0;k<N_Points;k++)
        BaseFunctions->GetDerivatives(MultiIndex3D::D000, X[k], Y[k], Z[k],
                                      AllPointValues[k]);

      for(k=0;k<N_Coarse;k++)
      {
        for(l=0;l<N_Points;l++)
          PointValues[l] = AllPointValues[l][k];

        nf->GetAllFunctionals(nullptr, nullptr, PointValues,
                              ret2+k*MaxN_BaseFunctions3D);
      }

      for(k=0;k<MaxN_BaseFunctions3D;k++)
        for(l=0;l<MaxN_BaseFunctions3D;l++)
          ret[k*MaxN_BaseFunctions3D+l] = 
                        ret2[l*MaxN_BaseFunctions3D+k];

      ProlongationMatrix3D[Coarse][refine][Fine][childnumber] = ret;
    }
    else
    {
      RefDesc = TDatabase::RefDescDB[refine+N_SHAPES];
      N_Children = RefDesc->GetN_Children();
  
      RefCell = GenerateRefElement3D(BaseFunctions->GetRefElement());
      RefCell->SetRefDesc(RefDesc);
      RefCell->Refine(1);
  
      for(j=0;j<N_Children;j++)
      {
        ret = new double[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
        ret2 = new double[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  
        cell = (const TGridCell *)((const TGridCell *)RefCell)->GetChild(j);
        Fine = static_cast<int>(FineElement.GetBaseFunct_ID())
               - N_BaseFuncts1D - N_BaseFuncts2D;
        auto nf = FineElement.GetNodalFunctional();
        nf->GetPointsForAll(N_Points, xi, eta, zeta);
  
        switch(cell->GetType())
        {
          case Hexahedron:
          case Brick:
            rt = GetRefTrans3D(ReferenceTransformation_type::HexaAffin);
            ((THexaAffin *)rt)->SetCell(cell);
            F_K = ReferenceTransformation_type::HexaAffin;
          break;
          case Tetrahedron:
            rt = GetRefTrans3D(ReferenceTransformation_type::TetraAffin);
            ((TTetraAffin *)rt)->SetCell(cell);
            F_K = ReferenceTransformation_type::TetraAffin;
          break;
          default:
            break;
        }
        GetOrigFromRef(F_K ,N_Points, xi, eta, zeta, X, Y, Z);
  
        for(k=0;k<N_Points;k++)
          BaseFunctions->GetDerivatives(MultiIndex3D::D000, X[k], Y[k], Z[k],
                        AllPointValues[k]);
  
        for(k=0;k<N_Coarse;k++)
        {
          for(l=0;l<N_Points;l++)
            PointValues[l] = AllPointValues[l][k];
  
          nf->GetAllFunctionals(nullptr, nullptr, PointValues,
                                ret2+k*MaxN_BaseFunctions3D);
        }

        for(k=0;k<MaxN_BaseFunctions3D;k++)
          for(l=0;l<MaxN_BaseFunctions3D;l++)
            ret[k*MaxN_BaseFunctions3D+l] = ret2[l*MaxN_BaseFunctions3D+k];
        delete[] ret2;
  
        ProlongationMatrix3D[Coarse][refine][Fine][j] = ret;
      } // endfor j
    }

    ret = ProlongationMatrix3D[Coarse][refine][Fine][childnumber]; 

    // RefCell->Derefine();
    // delete (TGridCell *)RefCell;
  }

  return ret;
}
#endif // 3D

//==============================================================================
double *GetRestrictionMatrix2D(FE_type parent, Refinements refine,
                               FE_type child, int childnumber)
{
  double *ret;
  int i,j,k,l1, l2;
  int N_Coarse, N_Fine, N_Children;
  double AllPointValues[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];
  const TRefDesc *RefDesc;
  BaseFunction_type Coarse, Fine;
  TBaseCell *RefCell;
  ReferenceTransformation_type F_K;
  TRefTrans2D *rt;

  double G[MaxN_BaseFunctions2D][MaxN_BaseFunctions2D];
  double Gret[MaxN_BaseFunctions2D][MaxN_BaseFunctions2D];
  double R[MaxN_BaseFunctions2D][MaxN_BaseFunctions2D];
  const TQuadFormula *Formula;
  double **CoarseBFData, **FineBFData, *PointData;
  double *FinePointData;
  int N_QuadPoints;
  int LDA=MaxN_BaseFunctions2D;
  
  auto FindLocalQuadFormula2D = [](const FiniteElement & Element)
    -> const TQuadFormula*
  {
    // find adequate quadrature formula for all elements
    int PolynomialDegree = Element.GetBaseFunct()->GetPolynomialDegree();
    auto reference_element = Element.GetBaseFunct()->GetRefElement();
    auto qf = QuadratureFormulaDatabase::qf_from_degree(
        2*PolynomialDegree, reference_element);
    return qf;
  };

  const FiniteElement CoarseElement(parent);
  Coarse = CoarseElement.GetBaseFunct_ID();
  const FiniteElement FineElement(child);
  Fine = FineElement.GetBaseFunct_ID();

  ret = RestrictionMatrix2D[Coarse][refine][Fine][childnumber]; 

  if(ret == nullptr)
  {
    // restriction matrix was not generated yet

    auto BaseFunctions = CoarseElement.GetBaseFunct();
    N_Coarse = BaseFunctions->GetDimension();

    std::fill((double *)G,
              (double *)G + MaxN_BaseFunctions2D*MaxN_BaseFunctions2D, 0.0);

    // build matrix G, gij = (uiH, ujH)
    Formula = FindLocalQuadFormula2D(CoarseElement);
    N_QuadPoints = Formula->GetN_QuadPoints();

    CoarseBFData = GetRefElementValues(*CoarseElement.GetBaseFunct(), *Formula,
                                       MultiIndex2D::D00);
    for(k=0;k<N_QuadPoints;k++)
    {
      PointData = CoarseBFData[k];
      double w = Formula->get_weight(k);
      for(i=0;i<N_Coarse;i++)
      {
        for(j=0;j<N_Coarse;j++)
        {
          // G is symmetric
          G[j][i] += w*PointData[i]*PointData[j];
        }
      }
    } // endfor k

    std::copy_n((double *)G, MaxN_BaseFunctions2D*MaxN_BaseFunctions2D,
                (double *)Gret);

    if(refine == NoRef)
    {
      std::fill((double *)R,
                (double *)R + MaxN_BaseFunctions2D*MaxN_BaseFunctions2D, 0.0);

      BaseFunctions = FineElement.GetBaseFunct();
      N_Fine = BaseFunctions->GetDimension();

      // build matrix R, rij = (uiH, ujh)
      Formula = FindLocalQuadFormula2D(FineElement);
      N_QuadPoints = Formula->GetN_QuadPoints();
      
      FineBFData = GetRefElementValues(*FineElement.GetBaseFunct(), *Formula,
                                       MultiIndex2D::D00);
      CoarseBFData = GetRefElementValues(*CoarseElement.GetBaseFunct(),
                                         *Formula, MultiIndex2D::D00);
      for(k=0;k<N_QuadPoints;k++)
      {
        FinePointData = FineBFData[k];
        PointData = CoarseBFData[k];
        double w = Formula->get_weight(k);
        for(l1=0;l1<N_Coarse;l1++)
        {
          for(l2=0;l2<N_Fine;l2++)
          {
            // = r(l1,l2) if row stored
            R[l2][l1] += w*PointData[l1]*FinePointData[l2];
          }
        }
      } // endfor k

      std::copy_n((double *)Gret, MaxN_BaseFunctions2D*MaxN_BaseFunctions2D,
                  (double *)G);

      // determine inv(G)*R
      SolveMultipleSystems((double *)G, (double *)R, 
                           N_Coarse, LDA, LDA, N_Fine);

      ret = new double[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
      for(l1=0;l1<N_Coarse;l1++)
        for(l2=0;l2<N_Fine;l2++)
          ret[l1*MaxN_BaseFunctions2D+l2] = R[l2][l1]; 
  
      RestrictionMatrix2D[Coarse][refine][Fine][childnumber] = ret;
    }
    else
    {
      RefDesc = TDatabase::RefDescDB[refine+N_SHAPES];
      N_Children = RefDesc->GetN_Children();
  
      RefCell = GenerateRefElement2D(BaseFunctions->GetRefElement());
      RefCell->SetRefDesc(RefDesc);
      RefCell->Refine(1);
  
      for(j=0;j<N_Children;j++)
      {
        // cout << "child: " << j << endl;
        std::fill((double *)R,
                  (double *)R + MaxN_BaseFunctions2D*MaxN_BaseFunctions2D, 0.0);
  
        N_Fine = FineElement.GetN_DOF();
  
        Formula = FindLocalQuadFormula2D(FineElement);
        N_QuadPoints = Formula->GetN_QuadPoints();
        FineBFData = GetRefElementValues(*FineElement.GetBaseFunct(), *Formula,
                                         MultiIndex2D::D00);
  
        F_K = FineElement.GetRefTransID();
  
        const TBaseCell* child_cell = ((const TGridCell *)RefCell)->GetChild(j);
        switch(F_K)
        {
          case ReferenceTransformation_type::QuadAffin:
          case ReferenceTransformation_type::QuadBilinear:
          case ReferenceTransformation_type::QuadIsoparametric:
            rt = GetRefTrans2D(ReferenceTransformation_type::QuadAffin);
            ((TQuadAffin *)rt)->SetCell(child_cell);
            F_K = ReferenceTransformation_type::QuadAffin;
            break;
          case ReferenceTransformation_type::TriaAffin:
          case ReferenceTransformation_type::TriaIsoparametric:
            rt = GetRefTrans2D(ReferenceTransformation_type::TriaAffin);
            ((TTriaAffin *)rt)->SetCell(child_cell);
            break;
          default:
            ErrThrow("unexpected reference transformation ",
                     static_cast<int>(F_K));
            break;
        }
        TQuadFormula qf_orig(*Formula);
        GetOrigFromRef(F_K, *Formula, qf_orig);
  
        for(k=0;k<N_QuadPoints;k++)
        {
          auto p = Formula->get_point(k);
          BaseFunctions->GetDerivatives(MultiIndex2D::D00, p.x, p.y,
                                        AllPointValues[k]);
        }

        for(k=0;k<N_QuadPoints;k++)
        {
          FinePointData = FineBFData[k];
          PointData = AllPointValues[k];
          double w = qf_orig.get_weight(k);
          for(l1=0;l1<N_Coarse;l1++)
          {
            for(l2=0;l2<N_Fine;l2++)
            {
              R[l2][l1] += w*PointData[l1]*FinePointData[l2];
            }
          }
        } // endfor k
       
        std::copy_n((double *)Gret, MaxN_BaseFunctions2D*MaxN_BaseFunctions2D,
                    (double *)G);

        // determine inv(G)*R
        SolveMultipleSystems((double *)G, (double *)R, N_Coarse, 
                             LDA, LDA, N_Fine);
  
        ret = new double[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
        // change from column to row storage
        for(l1=0;l1<N_Coarse;l1++)
          for(l2=0;l2<N_Fine;l2++)
            ret[l1*MaxN_BaseFunctions2D+l2] = R[l2][l1]; 
        RestrictionMatrix2D[Coarse][refine][Fine][j] = ret;
      } // endfor j

      RefCell->Derefine();
      delete (TGridCell *)RefCell;

    }

    ret = RestrictionMatrix2D[Coarse][refine][Fine][childnumber]; 

  }

  return ret;
}

//==============================================================================
#ifdef __3D__
double *GetRestrictionMatrix3D(FE_type parent, Refinements refine,
                               FE_type child, int childnumber)
{ 
  double *ret; // double *ret2;
  int i,j,k, l1, l2;
  int N_Coarse, N_Fine, N_Children; // int N_Points;
  double AllPointValues[MaxN_QuadPoints_3D][MaxN_BaseFunctions3D];
  //double PointValues[MaxN_PointsForNodal3D];
  const TRefDesc *RefDesc;
  TBaseCell *RefCell; // *cell;
  //TNodalFunctional3D *nf;
  ReferenceTransformation_type F_K;
  TRefTrans3D *rt;

  double G[MaxN_BaseFunctions3D][MaxN_BaseFunctions3D];
  double Gret[MaxN_BaseFunctions3D][MaxN_BaseFunctions3D];
  double R[MaxN_BaseFunctions3D][MaxN_BaseFunctions3D];
  const TQuadFormula *Formula;
  TQuadFormula qf_orig(QuadratureFormula_type::BaryCenterTetra); // dummy type
  double **CoarseBFData, **FineBFData, *PointData;
  double *FinePointData;
  double AbsDetjk[MaxN_QuadPoints_3D];
  int LDA=MaxN_BaseFunctions3D;
  
  auto FindLocalQuadFormula3D = [](const FiniteElement & local_used_element)
    -> const TQuadFormula *
  {
    // find adequate quadrature formula for all elements
    // and find needed reference transformation
    const BaseFunctions *bf = local_used_element.GetBaseFunct();
    int PolynomialDegree = bf->GetPolynomialDegree();
    const TQuadFormula * qf = QuadratureFormulaDatabase::qf_from_degree(
        2*PolynomialDegree, bf->GetRefElement());
    return qf;
  };

  const FiniteElement CoarseElement(parent);
  int Coarse = static_cast<int>(CoarseElement.GetBaseFunct_ID())
               - N_BaseFuncts1D - N_BaseFuncts2D;
  const FiniteElement FineElement(child);
  int Fine = static_cast<int>(FineElement.GetBaseFunct_ID())
             - N_BaseFuncts1D - N_BaseFuncts2D;

  ret = RestrictionMatrix3D[Coarse][refine][Fine][childnumber]; 

  if(ret == nullptr)
  {
    // restriction matrix was not generated yet

    auto BaseFunctions = CoarseElement.GetBaseFunct();
    N_Coarse = BaseFunctions->GetDimension();

    std::fill((double *)G,
              (double *)G + MaxN_BaseFunctions3D*MaxN_BaseFunctions3D, 0.0);

    // build matrix G, gij = (uiH, ujH)
    Formula = FindLocalQuadFormula3D(CoarseElement);
    int N_QuadPoints = Formula->GetN_QuadPoints();

    CoarseBFData = GetRefElementValues(*BaseFunctions, *Formula,
                                       MultiIndex3D::D000);
    for(k=0;k<N_QuadPoints;k++)
    {
      PointData = CoarseBFData[k];
      double w = Formula->get_weight(k);
      for(i=0;i<N_Coarse;i++)
      {
        for(j=0;j<N_Coarse;j++)
        {
          // G is symmetric
          G[j][i] += w*PointData[i]*PointData[j];
        }
      }
    } // endfor k

    std::copy_n((double *)G, MaxN_BaseFunctions3D*MaxN_BaseFunctions3D,
                (double *)Gret);

    if(refine == NoRef)
    {
      std::fill((double *)R,
                (double *)R + MaxN_BaseFunctions3D*MaxN_BaseFunctions3D, 0.0);

      BaseFunctions = FineElement.GetBaseFunct();
      N_Fine = BaseFunctions->GetDimension();

      // build matrix R, rij = (uiH, ujh)
      Formula = FindLocalQuadFormula3D(FineElement);
      int N_QuadPoints = Formula->GetN_QuadPoints();
      
      FineBFData = GetRefElementValues(*FineElement.GetBaseFunct(), *Formula,
                                       MultiIndex3D::D000);
      CoarseBFData = GetRefElementValues(*CoarseElement.GetBaseFunct(),
                                         *Formula, MultiIndex3D::D000);
      for(k=0;k<N_QuadPoints;k++)
      {
        FinePointData = FineBFData[k];
        PointData = CoarseBFData[k];
        double w = Formula->get_weight(k);
        for(l1=0;l1<N_Coarse;l1++)
        {
          for(l2=0;l2<N_Fine;l2++)
          {
            // = r(l1,l2) if row stored
            R[l2][l1] += w*PointData[l1]*FinePointData[l2];
          }
        }
      } // endfor k

      std::copy_n((double *)Gret, MaxN_BaseFunctions3D*MaxN_BaseFunctions3D,
                  (double *)G);

      // determine inv(G)*R
      SolveMultipleSystems((double *)G, (double *)R, 
                           N_Coarse, LDA, LDA, N_Fine);

      ret = new double[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
      for(l1=0;l1<N_Coarse;l1++)
        for(l2=0;l2<N_Fine;l2++)
          ret[l1*MaxN_BaseFunctions3D+l2] = R[l2][l1]; 
  
      RestrictionMatrix3D[Coarse][refine][Fine][childnumber] = ret;
    }
    else
    {
      RefDesc = TDatabase::RefDescDB[refine+N_SHAPES];
      N_Children = RefDesc->GetN_Children();
  
      RefCell = GenerateRefElement3D(BaseFunctions->GetRefElement());
      RefCell->SetRefDesc(RefDesc);
      RefCell->Refine(1);
  
      for(j=0;j<N_Children;j++)
      {
        // cout << "child: " << j << endl;
        std::fill((double *)R,
                  (double *)R + MaxN_BaseFunctions3D*MaxN_BaseFunctions3D, 0.0);
  
//        cell = RefCell->GetChild(j);
        Fine = static_cast<int>(FineElement.GetBaseFunct_ID())
               - N_BaseFuncts1D - N_BaseFuncts2D;
        auto FineBF = FineElement.GetBaseFunct();
        N_Fine = FineBF->GetDimension();
  
        Formula = FindLocalQuadFormula3D(FineElement);
        FineBFData = GetRefElementValues(*FineBF, *Formula, MultiIndex3D::D000);
  
        F_K = FineElement.GetRefTransID();
  
        const TBaseCell * child_cell = ((const TGridCell*)RefCell)->GetChild(j);
        switch(F_K)
        {
          case ReferenceTransformation_type::HexaAffin:
          case ReferenceTransformation_type::HexaTrilinear:
          case ReferenceTransformation_type::HexaIsoparametric:
            rt = GetRefTrans3D(ReferenceTransformation_type::HexaAffin);
            ((THexaAffin *)rt)->SetCell(child_cell);
            F_K = ReferenceTransformation_type::HexaAffin;
            break;
          case ReferenceTransformation_type::TetraAffin:
          case ReferenceTransformation_type::TetraIsoparametric:
            rt = GetRefTrans3D(ReferenceTransformation_type::TetraAffin);
            ((TTetraAffin *)rt)->SetCell(child_cell);
            break;
          default:
            ErrThrow("unexpected reference transformation ",
                     static_cast<int>(F_K));
            break;
        }
        GetOrigFromRef(F_K ,*Formula, qf_orig);
  
        unsigned int N_QuadPoints = Formula->GetN_QuadPoints();
        for(unsigned int k=0;k<N_QuadPoints;k++)
        {
          auto p = Formula->get_point(k);
          BaseFunctions->GetDerivatives(MultiIndex3D::D000, p.x, p.y, p.z,
                                        AllPointValues[k]);
        }

        for(unsigned int k=0;k<N_QuadPoints;k++)
        {
          FinePointData = FineBFData[k];
          PointData = AllPointValues[k];
          double w = Formula->get_weight(k)*AbsDetjk[k];
          for(l1=0;l1<N_Coarse;l1++)
          {
            for(l2=0;l2<N_Fine;l2++)
            {
              R[l2][l1] += w*PointData[l1]*FinePointData[l2];
            }
          }
        } // endfor k
       
        std::copy_n((double *)Gret, MaxN_BaseFunctions3D*MaxN_BaseFunctions3D,
                    (double *)G);

        // determine inv(G)*R
        SolveMultipleSystems((double *)G, (double *)R, N_Coarse, 
                             LDA, LDA, N_Fine);
  
        ret = new double[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
        // change from column to row storage
        for(l1=0;l1<N_Coarse;l1++)
          for(l2=0;l2<N_Fine;l2++)
            ret[l1*MaxN_BaseFunctions3D+l2] = R[l2][l1]; 
        RestrictionMatrix3D[Coarse][refine][Fine][j] = ret;
      } // endfor j

      // RefCell->Derefine();
      // delete (TGridCell *)RefCell;

    }

    ret = RestrictionMatrix3D[Coarse][refine][Fine][childnumber]; 

  }

  return ret;
}

#endif // 3D

//==============================================================================
double **GetRefElementValues(const BaseFunctions& bf, const TQuadFormula& qf,
                             MultiIndex2D MultiIndex)
{
  auto index = std::make_tuple(bf.GetID(), qf.get_type(), MultiIndex);
  auto it = RefElementValues2D.find(index);
  if(it == RefElementValues2D.end())
  {
    if(qf.get_dimension() != 2)
      ErrThrow("you want to evaluate 2D basis functions on a cell at quadrature"
               " points, but the quadrature formula is not a 2D formula");
    int n_quad_points = qf.GetN_QuadPoints();
    int n_functions = bf.GetDimension();
    int BaseVectDim = bf.GetBaseVectDim();
    auto Values = new double* [n_quad_points];
    auto AllValues = new double [n_quad_points*n_functions*BaseVectDim];
    for(int j=0;j<n_quad_points;j++)
      Values[j] = AllValues+j*n_functions*BaseVectDim;
    bf.GetDerivatives(MultiIndex, &qf, Values);
    RefElementValues2D[std::make_tuple(bf.GetID(), qf.get_type(), MultiIndex)]
        = Values;
  }
  return RefElementValues2D[index];
}

//==============================================================================
double **GetRefElementValues(const BaseFunctions& bf, const TQuadFormula& qf,
                             MultiIndex3D MultiIndex)
{
  auto index = std::make_tuple(bf.GetID(), qf.get_type(), MultiIndex);
  auto it = RefElementValues3D.find(index);
  if(it == RefElementValues3D.end())
  {
    if(qf.get_dimension() != 3)
      ErrThrow("you want to evaluate 3D basis functions on a cell at quadrature"
               " points, but the quadrature formula is not a 3D formula");
    int n_quad_points = qf.GetN_QuadPoints();
    int n_functions = bf.GetDimension();
    int BaseVectDim = bf.GetBaseVectDim();
    auto Values = new double* [n_quad_points];
    auto AllValues = new double [n_quad_points*n_functions*BaseVectDim];
    for(int j=0;j<n_quad_points;j++)
      Values[j] = AllValues+j*n_functions*BaseVectDim;
    bf.GetDerivatives(MultiIndex, &qf, Values);
    RefElementValues3D[std::make_tuple(bf.GetID(), qf.get_type(), MultiIndex)]
        = Values;
  }
  return RefElementValues3D[index];
}

//==============================================================================
double **GetJointDerivatives2D(const BaseFunctions& bf, const TQuadFormula& qf,
                               int joint, MultiIndex2D MultiIndex)
{
  auto index = std::make_tuple(bf.GetID(), qf.get_type(), joint, MultiIndex);
  auto it = JointDerivatives2D.find(index);
  if(it == JointDerivatives2D.end())
  {
    if(qf.get_dimension() != 1)
      ErrThrow("you want to evaluate 2D basis functions on a joint at "
               "quadrature points, but the quadrature formula is not a 1D "
               "formula");
    // data not generated yet
    int n_quad_points = qf.GetN_QuadPoints();
    int n_functions = bf.GetDimension();
    int BaseVectDim = bf.GetBaseVectDim();
    auto Values = new double* [n_quad_points];
    auto AllValues = new double [n_quad_points*n_functions*BaseVectDim];
    for(int j=0;j<n_quad_points;j++)
      Values[j] = AllValues+j*n_functions*BaseVectDim;
    bf.GetDerivatives(MultiIndex, &qf, joint, Values);
    JointDerivatives2D[std::make_tuple(bf.GetID(), qf.get_type(), joint,
                                       MultiIndex)] = Values;
  }
  return JointDerivatives2D[index];
}

//==============================================================================
double **GetJointDerivatives3D(const BaseFunctions& bf, const TQuadFormula& qf,
                               int joint, MultiIndex3D MultiIndex)
{
  auto index = std::make_tuple(bf.GetID(), qf.get_type(), joint, MultiIndex);
  auto it = JointDerivatives3D.find(index);
  if(it == JointDerivatives3D.end())
  {
    if(qf.get_dimension() != 2)
      ErrThrow("you want to evaluate 3D basis functions on a joint at "
               "quadrature points, but the quadrature formula is not a 2D "
               "formula");
    // data not generated yet
    int n_quad_points = qf.GetN_QuadPoints();
    int n_functions = bf.GetDimension();
    int BaseVectDim = bf.GetBaseVectDim();
    auto Values = new double* [n_quad_points];
    auto AllValues = new double [n_quad_points*n_functions*BaseVectDim];
    for(int j=0;j<n_quad_points;j++)
      Values[j] = AllValues+j*n_functions*BaseVectDim;
    bf.GetDerivatives(MultiIndex, &qf, joint, Values);
    JointDerivatives3D[std::make_tuple(bf.GetID(), qf.get_type(), joint,
                                       MultiIndex)] = Values;
  }
  return JointDerivatives3D[index];
}

//==============================================================================
double **GetOrigElementValues(const BaseFunctions& bf, MultiIndex2D MultiIndex)
{
  auto values = OrigElementValues2D[bf.GetID()-N_BaseFuncts1D]
                                   [static_cast<int>(MultiIndex)];
  if(values==nullptr)
  {
    int n_basis_functions = bf.GetDimension();
    int basevect_dim = bf.GetBaseVectDim();
    values = new double* [MaxN_QuadPoints_2D];
    auto aux = new double [MaxN_QuadPoints_2D*n_basis_functions*basevect_dim];
    for(int j=0;j<MaxN_QuadPoints_2D;j++)
      values[j] = aux+j*n_basis_functions*basevect_dim;
    OrigElementValues2D[bf.GetID()-N_BaseFuncts1D][static_cast<int>(MultiIndex)]
      = values;
  }
  return values;
}

//==============================================================================
double **GetOrigElementValues(const BaseFunctions& bf, MultiIndex3D MultiIndex)
{
  auto values = OrigElementValues3D[bf.GetID()-N_BaseFuncts1D-N_BaseFuncts2D]
                                   [static_cast<int>(MultiIndex)];
  if(values==nullptr)
  {
    int n_basis_functions = bf.GetDimension();
    int basevect_dim = bf.GetBaseVectDim();
    values = new double* [MaxN_QuadPoints_3D];
    auto aux = new double [MaxN_QuadPoints_3D*n_basis_functions*basevect_dim];
    for(int j=0;j<MaxN_QuadPoints_3D;j++)
      values[j] = aux+j*n_basis_functions*basevect_dim;
    OrigElementValues3D[bf.GetID()-N_BaseFuncts1D-N_BaseFuncts2D]
                       [static_cast<int>(MultiIndex)] = values;
  }
  return values;
}

//==============================================================================
void GetOrigValues(ReferenceTransformation_type RefTrans, double xi, double eta,
                   const BaseFunctions *bf, const TCollection *Coll,
                   const TGridCell *cell,
                   double *uref, double *uxiref, double *uetaref,
                   double *uorig, double *uxorig, double *uyorig)
{
  TRefTrans2D *rt = GetRefTrans2D(RefTrans);
  int N_BaseFunct = bf->GetDimension();
  int BaseVectDim = bf->GetBaseVectDim();
  rt->GetOrigValues(xi, eta, N_BaseFunct, uref, uxiref, uetaref, uorig, uxorig,
                    uyorig, BaseVectDim);
  bf->ChangeBF(Coll, cell, uorig);
  bf->ChangeBF(Coll, cell, uxorig);
  bf->ChangeBF(Coll, cell, uyorig);
}

//==============================================================================
#ifdef __3D__
void GetOrigValues(ReferenceTransformation_type RefTrans,
                   double xi, double eta, double zeta, const BaseFunctions *bf,
                   const TCollection *Coll, const TBaseCell *cell,
                   double *uref, double *uxiref, double *uetaref,
                   double *uzetaref,
                   double *uorig, double *uxorig, double *uyorig,
                   double *uzorig)
{
  const TRefTrans3D *rt = GetRefTrans3D(RefTrans);
  int N_BaseFunct = bf->GetDimension();
  int BaseVectDim = bf->GetBaseVectDim();
  rt->GetOrigValues(xi, eta, zeta, N_BaseFunct, uref, uxiref, uetaref, uzetaref,
                    uorig, uxorig, uyorig, uzorig, BaseVectDim);
  bf->ChangeBF(Coll, cell, uorig);
  bf->ChangeBF(Coll, cell, uxorig);
  bf->ChangeBF(Coll, cell, uyorig);
  bf->ChangeBF(Coll, cell, uzorig);
}
#endif // 3D

//==============================================================================
void GetOrigValues(ReferenceTransformation_type RefTrans, double zeta,
                   const BaseFunctions *bf, int edgeNumber,
                   const TCollection *Coll, const TBaseCell *cell,
                   double *uref, double *uxiref, double *uetaref,
                   double *uorig, double *uxorig, double *uyorig)
{
  TRefTrans2D *rt = GetRefTrans2D(RefTrans);
  int N_BaseFunct = bf->GetDimension();
  int BaseVectDim = bf->GetBaseVectDim();
  rt->GetOrigValues(edgeNumber, zeta, N_BaseFunct, uref, uxiref, uetaref,
                    uorig, uxorig, uyorig, BaseVectDim);
  bf->ChangeBF(Coll, cell, uorig);
  bf->ChangeBF(Coll, cell, uxorig);
  bf->ChangeBF(Coll, cell, uyorig);
}

//==============================================================================
void GetOrigValues(ReferenceTransformation_type RefTrans,
                   const std::vector<const BaseFunctions *>& BaseFuncts,
                   const TCollection *Coll, const TBaseCell *cell,
                   const TQuadFormula& qf, bool* Needs2ndDer)
{
#ifdef __2D__
  auto rt = GetRefTrans2D(RefTrans);
  int n_quad_points = qf.GetN_QuadPoints();
  unsigned int N_Sets = BaseFuncts.size();
  for(unsigned int i=0;i<N_Sets;i++)
  {
    const BaseFunctions& bf = *BaseFuncts[i];
    int BaseVectDim = bf.GetBaseVectDim();
    int N_Functs = bf.GetDimension();

    auto refvaluesD00 = GetRefElementValues(bf, qf, MultiIndex2D::D00);
    auto refvaluesD10 = GetRefElementValues(bf, qf, MultiIndex2D::D10);
    auto refvaluesD01 = GetRefElementValues(bf, qf, MultiIndex2D::D01);

    auto origvaluesD00 = GetOrigElementValues(bf, MultiIndex2D::D00);
    auto origvaluesD10 = GetOrigElementValues(bf, MultiIndex2D::D10);
    auto origvaluesD01 = GetOrigElementValues(bf, MultiIndex2D::D01);

    if(!Needs2ndDer[i])
    {
      for(int j=0;j < n_quad_points;j++)
      {
        auto p = qf.get_point(j);
        auto refD00 = refvaluesD00[j];
        auto refD10 = refvaluesD10[j];
        auto refD01 = refvaluesD01[j];

        auto origD10 = origvaluesD10[j];
        auto origD01 = origvaluesD01[j];
        auto origD00 = origvaluesD00[j];
        rt->GetOrigValues(p.x, p.y, N_Functs, refD00, refD10, refD01, origD00,
                          origD10, origD01, BaseVectDim);
      }
    }
    else
    {
      auto refvaluesD20 = GetRefElementValues(bf, qf, MultiIndex2D::D20);
      auto origvaluesD20 = GetOrigElementValues(bf, MultiIndex2D::D20);

      auto refvaluesD11 = GetRefElementValues(bf, qf, MultiIndex2D::D11);
      auto origvaluesD11 = GetOrigElementValues(bf, MultiIndex2D::D11);

      auto refvaluesD02 = GetRefElementValues(bf, qf, MultiIndex2D::D02);
      auto origvaluesD02 = GetOrigElementValues(bf, MultiIndex2D::D02);

      for(int j=0;j < n_quad_points;j++)
      {
        auto p = qf.get_point(j);
        auto refD00 = refvaluesD00[j];
        auto refD10 = refvaluesD10[j];
        auto refD01 = refvaluesD01[j];
        auto refD20 = refvaluesD20[j];
        auto refD11 = refvaluesD11[j];
        auto refD02 = refvaluesD02[j];
      
        auto origD10 = origvaluesD10[j];
        auto origD01 = origvaluesD01[j];
        auto origD00 = origvaluesD00[j];
        auto origD20 = origvaluesD20[j];
        auto origD11 = origvaluesD11[j];
        auto origD02 = origvaluesD02[j];
        rt->GetOrigAllDerivatives(
          p.x, p.y, N_Functs, refD00, refD10, refD01, refD20, refD11, refD02,
          origD00, origD10, origD01, origD20, origD11, origD02, BaseVectDim);
      }
      bf.ChangeBF(Coll, cell, n_quad_points, origvaluesD20);
      bf.ChangeBF(Coll, cell, n_quad_points, origvaluesD11);
      bf.ChangeBF(Coll, cell, n_quad_points, origvaluesD02);
    }
    bf.ChangeBF(Coll, cell, n_quad_points, origvaluesD00);
    bf.ChangeBF(Coll, cell, n_quad_points, origvaluesD10);
    bf.ChangeBF(Coll, cell, n_quad_points, origvaluesD01);
  }
#else // 2D -> 3D
  auto rt = GetRefTrans3D(RefTrans);
  int n_quad_points = qf.GetN_QuadPoints();
  unsigned int N_Sets = BaseFuncts.size();

  for(unsigned int i=0;i<N_Sets;i++)
  {
    const BaseFunctions& bf = *BaseFuncts[i];
    int N_Functs = bf.GetDimension();
    int BaseVectDim = bf.GetBaseVectDim();

    auto refvaluesD000 = GetRefElementValues(bf, qf, MultiIndex3D::D000);
    auto origvaluesD000 = GetOrigElementValues(bf, MultiIndex3D::D000);

    auto refvaluesD100 = GetRefElementValues(bf, qf, MultiIndex3D::D100);
    auto origvaluesD100 = GetOrigElementValues(bf, MultiIndex3D::D100);

    auto refvaluesD010 = GetRefElementValues(bf, qf, MultiIndex3D::D010);
    auto origvaluesD010 = GetOrigElementValues(bf, MultiIndex3D::D010);

    auto refvaluesD001 = GetRefElementValues(bf, qf, MultiIndex3D::D001);
    auto origvaluesD001 = GetOrigElementValues(bf, MultiIndex3D::D001);

    if(!Needs2ndDer[i])
    {
      for(int j=0;j<n_quad_points;j++)
      {
        auto p = qf.get_point(j);
        auto refD000 = refvaluesD000[j];
        auto refD100 = refvaluesD100[j];
        auto refD010 = refvaluesD010[j];
        auto refD001 = refvaluesD001[j];
      
        auto origD000 = origvaluesD000[j];
        auto origD100 = origvaluesD100[j];
        auto origD010 = origvaluesD010[j];
        auto origD001 = origvaluesD001[j];
        rt->GetOrigValues(p.x, p.y, p.z, N_Functs, refD000, refD100, refD010,
                          refD001, origD000, origD100, origD010, origD001,
                          BaseVectDim);
      }
    }
    else
    {
      auto refvaluesD200 = GetRefElementValues(bf, qf, MultiIndex3D::D200);
      auto origvaluesD200 = GetOrigElementValues(bf, MultiIndex3D::D200);

      auto refvaluesD110 = GetRefElementValues(bf, qf, MultiIndex3D::D110);
      auto origvaluesD110 = GetOrigElementValues(bf, MultiIndex3D::D110);

      auto refvaluesD101 = GetRefElementValues(bf, qf, MultiIndex3D::D101);
      auto origvaluesD101 = GetOrigElementValues(bf, MultiIndex3D::D101);

      auto refvaluesD020 = GetRefElementValues(bf, qf, MultiIndex3D::D020);
      auto origvaluesD020 = GetOrigElementValues(bf, MultiIndex3D::D020);

      auto refvaluesD011 = GetRefElementValues(bf, qf, MultiIndex3D::D011);
      auto origvaluesD011 = GetOrigElementValues(bf, MultiIndex3D::D011);

      auto refvaluesD002 = GetRefElementValues(bf, qf, MultiIndex3D::D002);
      auto origvaluesD002 = GetOrigElementValues(bf, MultiIndex3D::D002);

      for(int j=0;j < n_quad_points;j++)
      {
        auto p = qf.get_point(j);
        auto refD000 = refvaluesD000[j];
        auto refD100 = refvaluesD100[j];
        auto refD010 = refvaluesD010[j];
        auto refD001 = refvaluesD001[j];
        auto refD200 = refvaluesD200[j];
        auto refD110 = refvaluesD110[j];
        auto refD101 = refvaluesD101[j];
        auto refD020 = refvaluesD020[j];
        auto refD011 = refvaluesD011[j];
        auto refD002 = refvaluesD002[j];
      
        auto origD000 = origvaluesD000[j];
        auto origD100 = origvaluesD100[j];
        auto origD010 = origvaluesD010[j];
        auto origD001 = origvaluesD001[j];
        auto origD200 = origvaluesD200[j];
        auto origD110 = origvaluesD110[j];
        auto origD101 = origvaluesD101[j];
        auto origD020 = origvaluesD020[j];
        auto origD011 = origvaluesD011[j];
        auto origD002 = origvaluesD002[j];
        rt->GetOrigAllDerivatives(
          p.x, p.y, p.z, N_Functs, refD000, refD100, refD010, refD001,
          refD200, refD110, refD101, refD020, refD011, refD002,
          origD000, origD100, origD010, origD001,
          origD200, origD110, origD101, origD020, origD011, origD002,
          BaseVectDim);
      }
      bf.ChangeBF(Coll, cell, n_quad_points, origvaluesD200);
      bf.ChangeBF(Coll, cell, n_quad_points, origvaluesD110);
      bf.ChangeBF(Coll, cell, n_quad_points, origvaluesD101);
      bf.ChangeBF(Coll, cell, n_quad_points, origvaluesD020);
      bf.ChangeBF(Coll, cell, n_quad_points, origvaluesD011);
      bf.ChangeBF(Coll, cell, n_quad_points, origvaluesD002);
    } // endfor Needs2ndDer[i]
    bf.ChangeBF(Coll, cell, n_quad_points, origvaluesD000);
    bf.ChangeBF(Coll, cell, n_quad_points, origvaluesD100);
    bf.ChangeBF(Coll, cell, n_quad_points, origvaluesD010);
    bf.ChangeBF(Coll, cell, n_quad_points, origvaluesD001);
  } // endfor i
#endif
}

//==============================================================================
TRefTrans2D *GetRefTrans2D(ReferenceTransformation_type reftrans)
{
  int index = static_cast<int>(reftrans);
  if(ReferenceTrans2D[index] == nullptr)
  {
    switch(reftrans)
    {
      case ReferenceTransformation_type::TriaAffin:
        ReferenceTrans2D[index] = new TTriaAffin();
        break;
      case ReferenceTransformation_type::QuadAffin:
        ReferenceTrans2D[index] = new TQuadAffin();
        break;
      case ReferenceTransformation_type::QuadBilinear:
        ReferenceTrans2D[index] = new TQuadBilinear();
        break;
      case ReferenceTransformation_type::TriaIsoparametric:
        ReferenceTrans2D[index] = new TTriaIsoparametric();
        break;
      case ReferenceTransformation_type::QuadIsoparametric:
        ReferenceTrans2D[index] = new TQuadIsoparametric();
        break;
      default:
        ErrThrow("the index ", index, 
                 " does not map to a valid 2D reference transformation.");
        break;
    }
  }
  return ReferenceTrans2D[index];
}

//==============================================================================
#ifdef __3D__
TRefTrans3D *GetRefTrans3D(ReferenceTransformation_type reftrans)
{
  int index = static_cast<int>(reftrans);
  if(ReferenceTrans3D[index] == nullptr)
  {
    switch(reftrans)
    {
      case ReferenceTransformation_type::TetraAffin:
        ReferenceTrans3D[index] = new TTetraAffin();
        break;
      case ReferenceTransformation_type::TetraIsoparametric:
        ReferenceTrans3D[index] = new TTetraIsoparametric();
        break;
      case ReferenceTransformation_type::HexaAffin:
        ReferenceTrans3D[index] = new THexaAffin();
        break;
      case ReferenceTransformation_type::HexaTrilinear:
        ReferenceTrans3D[index] = new THexaTrilinear();
        break;
      case ReferenceTransformation_type::HexaIsoparametric:
        ReferenceTrans3D[index] = new THexaIsoparametric();
        break;
      default:
        ErrThrow("the index ", index, 
                 " does not map to a valid 3D reference transformation.");
        break;
    }
  }
  return ReferenceTrans3D[index];
}
#endif // 3D

//==============================================================================
void GetOrigFromRef(ReferenceTransformation_type RefTrans, int n_points,
                    const double *xi, const double *eta, double *X, double *Y)
{
  TRefTrans2D *rt = GetRefTrans2D(RefTrans);
  rt->GetOrigFromRef(n_points, xi, eta, X, Y);
}

//==============================================================================
#ifdef __3D__
void GetOrigFromRef(ReferenceTransformation_type RefTrans, int n_points,
                    const double *xi, const double *eta, const double *zeta,
                    double *X, double *Y, double *Z)
{
  const TRefTrans3D *rt = GetRefTrans3D(RefTrans);
  rt->GetOrigFromRef(n_points, xi, eta, zeta, X, Y, Z);
}
#endif // 3D

//==============================================================================
void GetOrigFromRef(ReferenceTransformation_type RefTrans,
                    const TQuadFormula& qf_ref, TQuadFormula& qf_orig)
{
#ifdef __2D__
  const TRefTrans2D *rt = GetRefTrans2D(RefTrans);
#else // 2D -> 3D
  const TRefTrans3D *rt = GetRefTrans3D(RefTrans);
#endif
  rt->GetOrigFromRef(qf_ref, qf_orig);
}

//==============================================================================
ReferenceTransformation_type GetOrig(std::vector<const FiniteElement*> used_fe,
                                     const TCollection *Coll,
                                     const TBaseCell *cell, bool *Needs2ndDer,
                                     TQuadFormula& qf_ref,
                                     TQuadFormula& qf_orig)
{
  // remove duplicates from used_fe, to avoid doing computations multiple times,
  // but we need to keep the order, because otherwise `Needs2ndDer` is no longer
  // compatible. To do so we use a std::set.
  {
    struct FEPointerCompare
    {
      bool operator()(const FiniteElement * fe1, const FiniteElement * fe2) const
      {
        return fe1->GetID() < fe2->GetID();
      }
    };
    std::set<const FiniteElement*, FEPointerCompare> fe_set;
    std::vector<const FiniteElement*> used_fe2;
    for(auto fe : used_fe)
    {
      if(fe_set.insert(fe).second)
        used_fe2.push_back(fe);
    }
    used_fe = used_fe2;
  }
  unsigned int N_LocalUsedElements = used_fe.size();
  std::vector<const BaseFunctions*> BaseFuncts(N_LocalUsedElements, nullptr);
  
#ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif


  // find adequate quadrature formula for all elements
  // and find needed reference transformation
  auto RefTrans = ReferenceTransformation_type::TriaAffin; // smallest value
  int MaxPolynomialDegree = 0;
  int MaxApproxOrder = 0;
  BFRefElements RefElement = BFRefElements::BFUnitTriangle;
  
  for(unsigned int i=0;i<N_LocalUsedElements;i++)
  {
    auto& element = *used_fe[i];
    BaseFuncts[i] = element.GetBaseFunct();
    MaxPolynomialDegree = std::max(
      element.GetBaseFunct()->GetPolynomialDegree(), MaxPolynomialDegree);
    MaxApproxOrder = std::max(element.GetBaseFunct()->GetAccuracy(),
                              MaxApproxOrder);

    RefTrans = std::max(element.GetRefTransID(), RefTrans);
    RefElement = element.GetBaseFunct()->GetRefElement();
  }

  // number of terms in products that need to be assembled
  int N_terms = (TDatabase::ParamDB->INTERNAL_PROBLEM_LINEAR ? 2 : 3);
  qf_ref = *QuadratureFormulaDatabase::qf_from_degree(
      N_terms*MaxPolynomialDegree, RefElement);
  switch(RefElement)
  {
    case BFRefElements::BFUnitSquare:
      if(TDatabase::ParamDB->INTERNAL_QUAD_RULE == 95)
      {
        qf_ref = *QuadratureFormulaDatabase::qf_from_degree(3, RefElement);
      }
      if(TDatabase::ParamDB->INTERNAL_QUAD_RULE == 97)
      {
        qf_ref = *QuadratureFormulaDatabase::qf_from_degree(18, RefElement);
      }
      if(TDatabase::ParamDB->INTERNAL_QUAD_RULE == 96)
      {
        qf_ref = *QuadratureFormulaDatabase::qf_from_degree(9, RefElement);
      }
      if(TDatabase::ParamDB->INTERNAL_QUAD_RULE == 99)
      {
        qf_ref = *QuadratureFormulaDatabase::qf_from_degree(18, RefElement);
      }
      if(TDatabase::ParamDB->INTERNAL_QUAD_QUAD<N_terms*MaxPolynomialDegree)
      {
        Output::print<2>("Quadrature formula for quadrilateral is ",
                         qf_ref.get_type());
        TDatabase::ParamDB->INTERNAL_QUAD_QUAD = N_terms*MaxPolynomialDegree;
        TDatabase::ParamDB->INTERNAL_MESH_CELL_TYPE = 4;
      }
      break;
    case BFRefElements::BFUnitTriangle:
      // more accurate quad rule for JohnMaubachTobiska1997.h
      if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == 97)
      {
        qf_ref = *QuadratureFormulaDatabase::qf_from_degree(11, RefElement);
      }
      if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == 98)
      {
        qf_ref = *QuadratureFormulaDatabase::qf_from_degree(5, RefElement);
      }
      if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == 99)
      {
        qf_ref = *QuadratureFormulaDatabase::qf_from_degree(20, RefElement);
      }
      if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == 96)
      {
        qf_ref = *QuadratureFormulaDatabase::qf_from_degree(20, RefElement);
      }
      if (TDatabase::ParamDB->INTERNAL_QUAD_RULE == -22)
      {
        qf_ref = *QuadratureFormulaDatabase::qf_from_degree(21, RefElement);
      }
      if(TDatabase::ParamDB->INTERNAL_QUAD_TRIA<N_terms*MaxPolynomialDegree)
      {
        Output::print<2>("Quadrature formula for triangles is ",
                         qf_ref.get_type());
        TDatabase::ParamDB->INTERNAL_QUAD_TRIA = N_terms*MaxPolynomialDegree;
        TDatabase::ParamDB->INTERNAL_MESH_CELL_TYPE = 3;
      }
      break;
    case BFRefElements::BFUnitHexahedron:
      if (TDatabase::ParamDB->INTERNAL_QUAD_HEXA<N_terms*MaxPolynomialDegree)
      {
#ifdef _MPI
        if(rank==0)
#endif
        {
          Output::info<2>("Quadrature formula for hexahedra is ",
                          qf_ref.get_type());
        }
        TDatabase::ParamDB->INTERNAL_QUAD_HEXA = N_terms*MaxPolynomialDegree;
      }
      break;
    case BFRefElements::BFUnitTetrahedron:
      if(TDatabase::ParamDB->INTERNAL_QUAD_TETRA < N_terms*MaxPolynomialDegree)
      {
#ifdef _MPI
        if(rank==0)
#endif
        {
          Output::info<2>("Quadrature formula for tetrahedra is ",
                          qf_ref.get_type());
        }
        TDatabase::ParamDB->INTERNAL_QUAD_TETRA = N_terms*MaxPolynomialDegree;
      }
      break;
    default:
      ErrThrow("wrong reference element ", RefElement);
      break;
  } // endswitch

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
      
      case BFRefElements::BFUnitHexahedron:
        RefTrans = ReferenceTransformation_type::HexaIsoparametric;
      break;

      case BFRefElements::BFUnitTetrahedron:
        RefTrans = ReferenceTransformation_type::TetraIsoparametric;
      break;
      
      default:
        ErrThrow("unexpected reference element ", RefElement);
        break;
    }
  } // endif IsIsoparametric

  // calculate the values of base functions and their derivatives
  // on the original element
#ifdef __2D__
  TRefTrans2D *rt = GetRefTrans2D(RefTrans);
#else
  TRefTrans3D *rt = GetRefTrans3D(RefTrans);
#endif
  switch(RefTrans)
  {
#ifdef __2D__
    case ReferenceTransformation_type::TriaAffin:
    case ReferenceTransformation_type::QuadAffin:
    case ReferenceTransformation_type::QuadBilinear:
      break;
    case ReferenceTransformation_type::TriaIsoparametric:
      ((TTriaIsoparametric *)rt)->SetApproximationOrder(MaxApproxOrder);
      ((TTriaIsoparametric *)rt)->SetQuadFormula(qf_ref.get_type());
    break;
    case ReferenceTransformation_type::QuadIsoparametric:
      ((TQuadIsoparametric *)rt)->SetApproximationOrder(MaxApproxOrder);
      ((TQuadIsoparametric *)rt)->SetQuadFormula(qf_ref.get_type());
    break;
#else
    case ReferenceTransformation_type::TetraAffin:
    case ReferenceTransformation_type::HexaAffin:
    case ReferenceTransformation_type::HexaTrilinear:
    break;
    case ReferenceTransformation_type::TetraIsoparametric:
      ((TTetraIsoparametric *)rt)->SetApproximationOrder(MaxApproxOrder);
      ((TTetraIsoparametric *)rt)->SetQuadFormula(qf_ref.get_type());
    break;
    case ReferenceTransformation_type::HexaIsoparametric:
      ((THexaIsoparametric *)rt)->SetApproximationOrder(MaxApproxOrder);
      ((THexaIsoparametric *)rt)->SetQuadFormula(qf_ref.get_type());
    break;
#endif
    default:
      ErrThrow("unexpected reference transformation ",
               static_cast<int>(RefTrans));
      break;
  } // endswitch
  rt->SetCell(cell);
  GetOrigValues(RefTrans, BaseFuncts, Coll, cell, qf_ref, Needs2ndDer);
  rt->GetOrigFromRef(qf_ref, qf_orig);

  return RefTrans;
}

//==============================================================================
void GetRefFromOrig(ReferenceTransformation_type RefTrans, double X, double Y,
                    double &xi, double &eta)
{
  TRefTrans2D *rt = GetRefTrans2D(RefTrans);
  rt->GetRefFromOrig(X, Y, xi, eta);
}

//==============================================================================
#ifdef __3D__
void GetRefFromOrig(ReferenceTransformation_type RefTrans,
                    double X, double Y, double Z,
                    double &xi, double &eta, double &zeta)
{
  const TRefTrans3D *rt = GetRefTrans3D(RefTrans);
  rt->GetRefFromOrig(X, Y, Z,  xi, eta, zeta);
}
#endif // 3D

//==============================================================================
void SetCellForRefTrans(const TBaseCell *cell,
                        ReferenceTransformation_type reftrans)
{
#ifdef __2D__
  TRefTrans2D *rt = GetRefTrans2D(reftrans);
#else // 2D -> 3D
  TRefTrans3D *rt = GetRefTrans3D(reftrans);
#endif
  rt->SetCell(cell);
}


} // end of namespace FEDatabase
