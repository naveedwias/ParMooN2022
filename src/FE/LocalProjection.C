// =======================================================================
// LocalProjection.C
//
// Purpose:   routines for local projection stabilization
//
// Author:    Gunar Matthies  2007/03/06
//
// =======================================================================

#include <Database.h>
#include <MooNMD_Io.h>

#include "FEDatabase.h"
#ifdef __2D__
  #include <FEFunction2D.h>
  #include "Matrix2D.h"
  #include "SquareMatrix2D.h"
#else  
  #include <FEFunction3D.h>
  #include "Matrix3D.h"
  #include "SquareMatrix3D.h"
#endif
#include "BaseCell.h"
#include <ConvDiff.h>
#include "LinAlg.h"

#include <cmath>
#include <string.h>
#include <stdlib.h>

#ifdef __2D__
FE_type GetElement2D(const TBaseCell *cell, int CoarseOrder)
{
  FE_type ele = ( FE_type )0;
  Shapes shapetype;

  shapetype = cell->GetType();
  switch(shapetype)
  {
    // regularly refined quadrilateral
    case Quadrangle:
    case Parallelogram:
    case Rectangle:
      switch(CoarseOrder)
      {
        case 0:
          ele = C_Q0_2D_Q_M;
        break;

        case 1:
          ele = D_P1_2D_Q_M;
        break;

        case 2:
          ele = D_P2_2D_Q_M;
        break;

        case 3:
          ele = D_P3_2D_Q_M;
        break;

        case 4:
          ele = D_P4_2D_Q_M;
        break;

        default:
          if(CoarseOrder<0)
          {
            ele = C_Q00_2D_Q_M;
          }
          else
          {
            ErrThrow("CoarseOrder: ", CoarseOrder,
                     " Projection space is defined up to order 4");
          }
      } // end switch CoarseOrder
    break; // end regularly refined quadrilateral

    case Triangle:
      switch(CoarseOrder)
      {
        case 0:
          ele = C_P0_2D_T_A;
        break;

        case 1:
          ele = D_P1_2D_T_A;
        break;

        case 2:
          ele = D_P2_2D_T_A;
        break;

        case 3:
          ele = D_P3_2D_T_A;
        break;

        case 4:
          ele = D_P4_2D_T_A;
        break;

        default:
          if(CoarseOrder<0)
          {
            ele = C_P00_2D_T_A;
          }
          else
          {
            ErrThrow("CoarseOrder: ", CoarseOrder,
                     " Projection space is defined up to order 4");
          }
      } // end switch CoarseOrder
    break;
    default:
      ErrThrow("Invalid shape");
  } // end switch reftype
  return ele;
}

/** Navier--Stokes type 2 (NSTYPE==2) with C*/
/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        TMatrix *B1T, TMatrix *B2T, TMatrix *C,
        double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, value, value1, value2;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  double *r1, *r2, *r3;
  const int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  const int *BTRowPtr, *BTKCol;
  double *AEntries, *B1Entries, *B2Entries;
  double *B1TEntries, *B2TEntries;
  int N_Active;
  double *CEntries;
  const int *CRowPtr, *CKCol;

  ARowPtr = A->get_row_ptr();
  AKCol = A->get_vector_columns();
  AEntries = A->GetEntries();

  BRowPtr = B1->get_row_ptr();
  BKCol = B1->get_vector_columns();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

  BTRowPtr = B1T->get_row_ptr();
  BTKCol = B1T->get_vector_columns();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();

  CRowPtr = C->get_row_ptr();
  CKCol = C->get_vector_columns();
  CEntries = C->GetEntries();

  N_UDOF = A->get_n_rows();
  N_PDOF = B1->get_n_rows();

  u1 = x;
  u2 = u1+N_UDOF;
  p  = u2+N_UDOF;

  v1 = b;
  v2 = v1+N_UDOF;
  q  = v2+N_UDOF;

  r1 = r;
  r2 = r1+N_UDOF;
  r3 = r2+N_UDOF;

  N_Active = A->get_n_active_rows();

  j = ARowPtr[0];
  for(i=0;i<N_UDOF;i++)
  {
    s = v1[i];
    t = v2[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s -= value * u1[index];
      t -= value * u2[index];
    }
    r1[i] = s;
    r2[i] = t;
  } // endfor i

  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = q[i];
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value1 = B1Entries[j];
      value2 = B2Entries[j];
      s -= value1 * u1[index] + value2 * u2[index];
    } // endfor j
    r3[i] = s;
  } // endfor i

  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
    }
    r1[i] -= s;
    r2[i] -= t;
  } // endfor i

  j = CRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = CRowPtr[i+1];
    for(;j<k;j++)
    {
      s += CEntries[j] * p[CKCol[j]]; // plus is right sign
    }
    r3[i] -= s;
  } // endfor i
}

// TODO Disabled due to use of deprecated parameter.
//void Defect_NSE2C(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
//{
//  int N_UDOF,N_PDOF;
//
//  CoupledDefect(A[0], B[0], B[1], B[2], B[3], B[4], x, b, r);
//  N_UDOF = A[0]->get_n_rows();
//  N_PDOF = B[0]->get_n_rows();
//  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) //deprecated!
//    IntoL20Vector2D(r+2*N_UDOF, N_PDOF, TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
//  return;
//}

/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A, TMatrix *B1, TMatrix *B2,
        TMatrix *B1T, TMatrix *B2T, TMatrix *C,
        double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, value, value1, value2;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  const int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  const int *BTRowPtr, *BTKCol;
  double *AEntries, *B1Entries, *B2Entries;
  double *B1TEntries, *B2TEntries;
  int N_Active;
  double *CEntries;
  const int *CRowPtr, *CKCol;

  ARowPtr = A->get_row_ptr();
  AKCol = A->get_vector_columns();
  AEntries = A->GetEntries();

  BRowPtr = B1->get_row_ptr();
  BKCol = B1->get_vector_columns();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

  BTRowPtr = B1T->get_row_ptr();
  BTKCol = B1T->get_vector_columns();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();

  CRowPtr = C->get_row_ptr();
  CKCol = C->get_vector_columns();
  CEntries = C->GetEntries();

  N_UDOF = A->get_n_rows();
  N_PDOF = B1->get_n_rows();

  u1 = x;
  u2 = u1+N_UDOF;
  p  = u2+N_UDOF;

  v1 = y;
  v2 = v1+N_UDOF;
  q  = v2+N_UDOF;

  N_Active = A->get_n_active_rows();
  j = ARowPtr[0];

  for(i=0;i<N_UDOF;i++)
  {
    s = 0;
    t = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = AEntries[j];
      s += value * u1[index];
      t += value * u2[index];
    }
    v1[i] = s;
    v2[i] = t;
  } // endfor i

  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value1 = B1Entries[j];
      value2 = B2Entries[j];
      s += value1 * u1[index] + value2 * u2[index];
    } // endfor j
    q[i] = s;
  } // endfor i

  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
    }
    v1[i] += s;
    v2[i] += t;
  } // endfor i

  j = CRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = CRowPtr[i+1];
    for(;j<k;j++)
    {
      s += CEntries[j] * p[CKCol[j]]; // plus is right sign
    } // endfor j
    q[i] += s;
  } // endfor i

  return;
}

void MatVect_NSE2C(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], B[0], B[1], B[2], B[3], B[4], x, y);
  return;
}

// stabilisation of full gradient (velocity or pressure)
void UltraLocalProjection(void* A, bool ForPressure)
{
  int i,j,k,l,m;
  int N_Cells;
  int CoarseOrder;
  int N_CoarseDOF, N_DOF;
  std::shared_ptr<const TFESpace2D> fespace;
  TBaseCell *cell;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  const int *RowPtr, *KCol;
  double *Entries;
  int OrderDiff;
  double lpcoeff, lpexponent;

  if(!(TDatabase::ParamDB->LP_FULL_GRADIENT) && !(ForPressure))
  {
    ErrThrow("Local projection stabilization is implemented only for full gradient!");
  }

  if(ForPressure)
  {
    lpcoeff = -(TDatabase::ParamDB->LP_PRESSURE_COEFF);
    lpexponent = TDatabase::ParamDB->LP_PRESSURE_EXPONENT;
    OrderDiff = TDatabase::ParamDB->LP_PRESSURE_ORDER_DIFFERENCE;
  }
  else
  {
    lpcoeff = TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF;
    lpexponent = TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT;
    OrderDiff = TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE;
  }

  if(ForPressure)
  {
    fespace = ((TMatrix2D *)A)->GetTestSpace2D();
    ActiveBound = -1;
    RowPtr = ((TMatrix2D *)A)->get_row_ptr();
    KCol = ((TMatrix2D *)A)->get_vector_columns();
    Entries = ((TMatrix2D *)A)->GetEntries();
    // cout << "for pressure" << endl;
  }
  else
  {
    fespace = ((TSquareMatrix2D *)A)->GetFESpace2D();
    ActiveBound = fespace->get_n_active();
    RowPtr = ((TSquareMatrix2D *)A)->get_row_ptr();
    KCol = ((TSquareMatrix2D *)A)->get_vector_columns();
    Entries = ((TSquareMatrix2D *)A)->GetEntries();
//     cout << "not for pressure" << endl;
  }

  auto Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();
    auto CurrentElement = fespace->get_fe(i);

    auto BF = CurrentElement.GetBaseFunct();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    const FiniteElement CoarseElement(GetElement2D(cell, CoarseOrder));
    auto CoarseBF = CoarseElement.GetBaseFunct();

    // quadrature formula on cell
    std::vector<const FiniteElement*> UsedElements{{&CoarseElement, &CurrentElement}};
    FEDatabase::GetOrig(UsedElements, Coll, cell, SecondDer, qf_ref, qf_orig);
    N_Points = qf_orig.GetN_QuadPoints();

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = FEDatabase::GetOrigElementValues(*CoarseBF,
                                                    MultiIndex2D::D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = qf_orig.get_weight(j);
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && std::abs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*sizeof(double));
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    PCValues = FEDatabase::GetOrigElementValues(*CoarseBF, MultiIndex2D::D00);

    ChildValuesX = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D10);
    ChildValuesY = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D01);

    memset(H, 0, N_CoarseDOF*2*N_DOF*sizeof(double));

    memset(LocMatrix, 0, N_DOF*N_DOF*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      w = qf_orig.get_weight(j);
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*2*N_DOF+l      ] += val*ChildValueX[l];
          H[k*2*N_DOF+l+N_DOF] += val*ChildValueY[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*( ChildValueX[k]*ChildValueX[l]
                                     +ChildValueY[k]*ChildValueY[l]);
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*2*N_DOF*sizeof(double));

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 2*N_DOF, 2*N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2]*( H[i1*2*N_DOF+l      ]*H[i2*2*N_DOF+m      ]
                                           +H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF]);
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m      ];
          s += P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l      ];
          s += P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    auto DOF = fespace->GetGlobalDOF(i);

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if((dof<ActiveBound) || (ForPressure))
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              Entries[p] += lpcoeff*std::pow(hK,lpexponent)*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // 

double UltraLocalError(TFEFunction2D *uh, DoubleFunct2D *ExactU,
        double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l;
  int N_Cells;
  int CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TBaseCell *cell;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double s;
  int i1, i2;
  double hK;
  double *Values;
  double error, locerror;
  double exactval[4];

  error = 0.0;

  auto fespace = uh->GetFESpace2D();
  Values = uh->GetValues();

  auto Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);
  
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);

  for(i=0;i<N_Cells;i++)
  {
    locerror = 0.0;
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    auto DOF = fespace->GetGlobalDOF(i);

    auto CurrentElement = fespace->get_fe(i);

    auto BF = CurrentElement.GetBaseFunct();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    const FiniteElement CoarseElement(GetElement2D(cell, CoarseOrder));
    auto CoarseBF = CoarseElement.GetBaseFunct();

    // quadrature formula on cell
    std::vector<const FiniteElement*> UsedElements{{&CoarseElement, &CurrentElement}};
    FEDatabase::GetOrig(UsedElements, Coll, cell, SecondDer, qf_ref, qf_orig);
    N_Points = qf_orig.GetN_QuadPoints();

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = FEDatabase::GetOrigElementValues(*CoarseBF, 
                                                    MultiIndex2D::D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = qf_orig.get_weight(j);
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && std::abs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*sizeof(double));
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    PCValues = FEDatabase::GetOrigElementValues(*CoarseBF, MultiIndex2D::D00);

    ChildValuesX = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D10);
    ChildValuesY = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D01);

    // only two right-hand sides (x and y derivative)
    memset(H, 0, N_CoarseDOF*2*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];

      // calculate gradient of discrete uh in this quadrature point
      valx = 0.0;
      valy = 0.0;
      for(k=0;k<N_DOF;k++)
      {
        val = Values[DOF[k]];
        valx += ChildValueX[k]*val;
        valy += ChildValueY[k]*val;
      }

      auto p = qf_orig.get_point(j);
      // get gradient of exact u
      ExactU(p.x, p.y, exactval);

      valx -= exactval[1];
      valy -= exactval[2];

      w = qf_orig.get_weight(j);
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        H[k*2  ] += val*valx;
        H[k*2+1] += val*valy;
      } // end for k

      // grad-grad term
      locerror += w*(valx*valx + valy*valy);

    } // end for j
    memcpy(P, H, N_CoarseDOF*2*sizeof(double));

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 2, 2);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    s = 0;
    for(i1=0;i1<N_CoarseDOF;i1++)
      for(i2=0;i2<N_CoarseDOF;i2++)
        s += Gsave[i1*N_CoarseDOF+i2]*( H[i1*2  ]*H[i2*2  ]
                                       +H[i1*2+1]*H[i2*2+1]);
    locerror += s;

    // grad-proj coupling
    s = 0;
    for(i2=0;i2<N_CoarseDOF;i2++)
    {
      s += P[i2*2  ] * H[i2*2  ];
      s += P[i2*2+1] * H[i2*2+1];
    }
    for(i1=0;i1<N_CoarseDOF;i1++)
    {
      s += P[i1*2  ] * H[i1*2  ];
      s += P[i1*2+1] * H[i1*2+1];
    }
    locerror -= s;

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    error += lpcoeff*std::pow(hK,lpexponent)*locerror;
  } // endfor i

  return error;
} // UltraLocalError

void AddStreamlineTerm(TSquareMatrix2D* A, TFEFunction2D *uh1,
                       TFEFunction2D *uh2,
                       double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m;
  int N_Cells;
  int CoarseOrder;
  int N_CoarseDOF, N_DOF;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValues, *ChildValue;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  const int *RowPtr, *KCol;
  double *Entries;
  double *Values1, *Values2;
  double BValue[MaxN_BaseFunctions2D];

  auto fespace = A->GetFESpace2D();
  ActiveBound = fespace->get_n_active();
  RowPtr = A->get_row_ptr();
  KCol = A->get_vector_columns();
  Entries = A->GetEntries();
  // cout << "" << endl;

  auto Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);
  
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);

  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    auto CurrentElement = fespace->get_fe(i);

    auto BF = CurrentElement.GetBaseFunct();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    const FiniteElement CoarseElement(GetElement2D(cell, CoarseOrder));
    auto CoarseBF = CoarseElement.GetBaseFunct();

    // quadrature formula on cell
    std::vector<const FiniteElement*> UsedElements{{&CoarseElement, &CurrentElement}};
    FEDatabase::GetOrig(UsedElements, Coll, cell, SecondDer, qf_ref, qf_orig);
    N_Points = qf_orig.GetN_QuadPoints();

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = FEDatabase::GetOrigElementValues(*CoarseBF,
                                                    MultiIndex2D::D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = qf_orig.get_weight(j);
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && std::abs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*sizeof(double));
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    auto DOF =  fespace->GetGlobalDOF(i);

    Values1 = uh1->GetValues();
    Values2 = uh2->GetValues();

    PCValues = FEDatabase::GetOrigElementValues(*CoarseBF, MultiIndex2D::D00);

    ChildValuesX = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D10);
    ChildValuesY = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D01);
    ChildValues  = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D00);

    memset(H, 0, N_CoarseDOF*N_DOF*sizeof(double));

    memset(LocMatrix, 0, N_DOF*N_DOF*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValue  = ChildValues[j];
      w = qf_orig.get_weight(j);
      valx = 0.0;
      valy = 0.0;
      // compute components of uh in j
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        valx += ChildValue[k]*Values1[l];
        valy += ChildValue[k]*Values2[l];
      }
      for(k=0;k<N_DOF;k++)
      {
        BValue[k] = valx*ChildValueX[k] + valy*ChildValueY[k];
      }
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*N_DOF+l      ] += val*BValue[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*BValue[k]*BValue[l];
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*N_DOF*sizeof(double));

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, N_DOF, N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2] * H[i1*N_DOF+l] * H[i2*N_DOF+m];
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*N_DOF+l      ] * H[i2*N_DOF+m      ];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*N_DOF+m      ] * H[i1*N_DOF+l      ];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if(dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              Entries[p] += lpcoeff*std::pow(hK,lpexponent)*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // AddStreamlineTerm

void AddStreamlineTermPWConst(TSquareMatrix2D* A, TFEFunction2D *uh1,
                              TFEFunction2D *uh2,
                              double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m;
  int N_Cells;
  int CoarseOrder;
  int N_CoarseDOF, N_DOF;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValues, *ChildValue;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  const int *RowPtr, *KCol;
  double *Entries;
  double *Values1, *Values2;
  double BValue[MaxN_BaseFunctions2D];

  auto fespace = A->GetFESpace2D();
  ActiveBound = fespace->get_n_active();
  RowPtr = A->get_row_ptr();
  KCol = A->get_vector_columns();
  Entries = A->GetEntries();
  // cout << "" << endl;

  auto Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);
  
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);

  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    auto CurrentElement = fespace->get_fe(i);

    auto BF = CurrentElement.GetBaseFunct();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    const FiniteElement CoarseElement(GetElement2D(cell, CoarseOrder));
    auto CoarseBF = CoarseElement.GetBaseFunct();

    // quadrature formula on cell
    std::vector<const FiniteElement*> UsedElements{{&CoarseElement, &CurrentElement}};
    FEDatabase::GetOrig(UsedElements, Coll, cell, SecondDer, qf_ref,
                           qf_orig);
    N_Points = qf_orig.GetN_QuadPoints();

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = FEDatabase::GetOrigElementValues(*CoarseBF, 
MultiIndex2D::D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = qf_orig.get_weight(j);
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && std::abs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*sizeof(double));
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    auto DOF = fespace->GetGlobalDOF(i);

    Values1 = uh1->GetValues();
    Values2 = uh2->GetValues();

    PCValues = FEDatabase::GetOrigElementValues(*CoarseBF, MultiIndex2D::D00);

    ChildValuesX = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D10);
    ChildValuesY = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D01);
    ChildValues  = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D00);

    memset(H, 0, N_CoarseDOF*N_DOF*sizeof(double));

    memset(LocMatrix, 0, N_DOF*N_DOF*sizeof(double));

    // calculate pw constant approximation of velocity field (uh1, uh2)
    val  = 0.0;
    valx = 0.0;
    valy = 0.0;
    for(j=0;j<N_Points;j++)
    {
      ChildValue  = ChildValues[j];
      w = qf_orig.get_weight(j);
      // compute components of uh in j
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        val  += w*ChildValue[k]*1;
        valx += w*ChildValue[k]*Values1[l];
        valy += w*ChildValue[k]*Values2[l];
      }
    }
    valx /= val;
    valy /= val;

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValue  = ChildValues[j];
      w = qf_orig.get_weight(j);
      /*
      valx = 0.0;
      valy = 0.0;
      // compute components of uh in j
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        valx += ChildValue[k]*Values1[l];
        valy += ChildValue[k]*Values2[l];
      }
      */
      for(k=0;k<N_DOF;k++)
      {
        BValue[k] = valx*ChildValueX[k] + valy*ChildValueY[k];
      }
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*N_DOF+l      ] += val*BValue[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*BValue[k]*BValue[l];
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*N_DOF*sizeof(double));

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, N_DOF, N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2] * H[i1*N_DOF+l] * H[i2*N_DOF+m];
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*N_DOF+l      ] * H[i2*N_DOF+m      ];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*N_DOF+m      ] * H[i1*N_DOF+l      ];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if(dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              Entries[p] += lpcoeff*std::pow(hK,lpexponent)*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // AddStreamlineTermPWConst

void AddDivergenceTerm(TSquareMatrix2D *A11,TSquareMatrix2D *A12,
                       TSquareMatrix2D *A21,TSquareMatrix2D *A22,
                       double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m;
  int N_Cells;
  int CoarseOrder;
  int N_CoarseDOF, N_DOF;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val;
  double LocMatrixA11[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixA12[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double LocMatrixA22[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s11, s12, s22;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  const int *RowPtr, *KCol;
  double *EntriesA11, *EntriesA12, *EntriesA21, *EntriesA22;

  auto fespace = A11->GetFESpace2D();
  ActiveBound = fespace->get_n_active();
  RowPtr = A11->get_row_ptr();
  KCol = A11->get_vector_columns();
  EntriesA11 = A11->GetEntries();
  EntriesA12 = A12->GetEntries();
  EntriesA21 = A21->GetEntries();
  EntriesA22 = A22->GetEntries();
  // cout << "" << endl;


  auto Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);
  
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);

  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    auto CurrentElement = fespace->get_fe(i);

    auto BF = CurrentElement.GetBaseFunct();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    const FiniteElement CoarseElement(GetElement2D(cell, CoarseOrder));
    auto CoarseBF = CoarseElement.GetBaseFunct();

    // quadrature formula on cell
    std::vector<const FiniteElement*> UsedElements{{&CoarseElement, &CurrentElement}};
    FEDatabase::GetOrig(UsedElements, Coll, cell, SecondDer, qf_ref, qf_orig);
    N_Points = qf_orig.GetN_QuadPoints();

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = FEDatabase::GetOrigElementValues(*CoarseBF,
                                                    MultiIndex2D::D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = qf_orig.get_weight(j);
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && std::abs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*sizeof(double));
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    auto DOF = fespace->GetGlobalDOF(i);

    PCValues = FEDatabase::GetOrigElementValues(*CoarseBF, MultiIndex2D::D00);

    ChildValuesX = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D10);
    ChildValuesY = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D01);

    memset(H, 0, N_CoarseDOF*2*N_DOF*sizeof(double));

    memset(LocMatrixA11, 0, N_DOF*N_DOF*sizeof(double));
    memset(LocMatrixA12, 0, N_DOF*N_DOF*sizeof(double));
    memset(LocMatrixA22, 0, N_DOF*N_DOF*sizeof(double));
    // since A21=transpose(A12) A21 is not explicitly needed

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      w = qf_orig.get_weight(j);
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*2*N_DOF+l      ] += val*ChildValueX[l];
          H[k*2*N_DOF+l+N_DOF] += val*ChildValueY[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrixA11[k*N_DOF+l] += w*ChildValueX[k]*ChildValueX[l];
          LocMatrixA12[k*N_DOF+l] += w*ChildValueX[k]*ChildValueY[l];
          LocMatrixA22[k*N_DOF+l] += w*ChildValueY[k]*ChildValueY[l];
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*2*N_DOF*sizeof(double));

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 2*N_DOF, 2*N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
        s12 = 0;
        s22 = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
          {
            s11 += Gsave[i1*N_CoarseDOF+i2]
                  *H[i1*2*N_DOF+l      ]*H[i2*2*N_DOF+m      ];
            s12 += Gsave[i1*N_CoarseDOF+i2]
                  *H[i1*2*N_DOF+l      ]*H[i2*2*N_DOF+m+N_DOF];
            s22 += Gsave[i1*N_CoarseDOF+i2]
                  *H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF];
          }
        LocMatrixA11[l*N_DOF+m] += s11;
        LocMatrixA12[l*N_DOF+m] += s12;
        LocMatrixA22[l*N_DOF+m] += s22;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s11 = 0;
        s12 = 0;
        s22 = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s11 += P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m      ];
          s12 += P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m+N_DOF];
          s22 += P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s11 += P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l      ];
          s12 += P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l      ];
          s22 += P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF];
        }
        LocMatrixA11[l*N_DOF+m] -= s11;
        LocMatrixA12[l*N_DOF+m] -= s12;
        LocMatrixA22[l*N_DOF+m] -= s22;
//         LocMatrixA21[m*N_DOF+l] = LocMatrixA12[l*N_DOF+m];
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if (dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              EntriesA11[p] += lpcoeff*std::pow(hK,lpexponent)*LocMatrixA11[l*N_DOF+m];
              EntriesA12[p] += lpcoeff*std::pow(hK,lpexponent)*LocMatrixA12[l*N_DOF+m];
              // since A21=transpose(A12)
              EntriesA21[p] += lpcoeff*std::pow(hK,lpexponent)*LocMatrixA12[m*N_DOF+l];
              EntriesA22[p] += lpcoeff*std::pow(hK,lpexponent)*LocMatrixA22[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // AddDivergenceTerm

/** Navier--Stokes type 4 (NSTYPE==4) with C*/
/** matrix * vector for coupled Stokes / Navier-Stokes system */
void CoupledMatVect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A21,
                    TSquareMatrix *A22, TMatrix *B1, TMatrix *B2,
                    TMatrix *B1T, TMatrix *B2T,
                    TMatrix *C,
                    double *x, double *y)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, value, value1, value2,value3;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  const int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  const int *BTRowPtr, *BTKCol;
  double *A11Entries, *B1Entries, *B2Entries;
  double *B1TEntries, *B2TEntries;
  double *A12Entries, *A21Entries, *A22Entries;
  double *CEntries;
  const int *CRowPtr, *CKCol;
  int N_Active;

  ARowPtr = A11->get_row_ptr();
  AKCol = A11->get_vector_columns();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();

  BRowPtr = B1->get_row_ptr();
  BKCol = B1->get_vector_columns();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

  BTRowPtr = B1T->get_row_ptr();
  BTKCol = B1T->get_vector_columns();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();

  CRowPtr = C->get_row_ptr();
  CKCol = C->get_vector_columns();
  CEntries = C->GetEntries();

  N_UDOF = A11->get_n_rows();
  N_PDOF = B1->get_n_rows();

  u1 = x;
  u2 = u1+N_UDOF;
  p  = u2+N_UDOF;

  v1 = y;
  v2 = v1+N_UDOF;
  q  = v2+N_UDOF;

  N_Active = A11->get_n_active_rows();
  j = ARowPtr[0];
  // real dof
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A12Entries[j];
      value2 = A21Entries[j];
      value3 = A22Entries[j];
      s += value * u1[index] + value1 * u2[index];
      t += value2* u1[index] + value3 * u2[index];
    }
    v1[i] = s;
    v2[i] = t;
  } // endfor i
  // Dirichlet and hanging nodes
  j = ARowPtr[N_Active];
  for(i=N_Active;i<N_UDOF;i++)
  {
    s = 0;
    t = 0;
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A22Entries[j];
      s += value * u1[index];
      t += value1 * u2[index];
    }
    v1[i] = s;
    v2[i] = t;
  } // endfor i

  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value1 = B1Entries[j];
      value2 = B2Entries[j];
      s += value1 * u1[index] + value2 * u2[index];
    } // endfor j
    q[i] = s;
  } // endfor i

  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
    }
    v1[i] += s;
    v2[i] += t;
  } // endfor i

  j = CRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = CRowPtr[i+1];
    for(;j<k;j++)
    {
      s += CEntries[j] * p[CKCol[j]]; // plus is right sign
    } // endfor j
    q[i] += s;
  } // endfor i

  return;
}

void MatVect_NSE4C(TSquareMatrix **A, TMatrix **B, double *x, double *y)
{
  CoupledMatVect(A[0], A[1], A[2], A[3], B[0], B[1], B[2], B[3], B[4], x, y);
  return;
}

/** r := b - A * x */
void CoupledDefect(TSquareMatrix *A11, TSquareMatrix *A12, TSquareMatrix *A21,
                   TSquareMatrix *A22, TMatrix *B1, TMatrix *B2,
                   TMatrix *B1T, TMatrix *B2T,
                   TMatrix *C,
                   double *x, double *b, double *r)
{
  int N_UDOF, N_PDOF;
  int i,j,k,index;
  double s, t, value, value1, value2, value3;
  double *u1, *u2, *p;
  double *v1, *v2, *q;
  double *r1, *r2, *r3;
  const int *ARowPtr, *BRowPtr, *AKCol, *BKCol;
  const int *BTRowPtr, *BTKCol;
  double *A11Entries, *B1Entries, *B2Entries;
  double *B1TEntries, *B2TEntries;
  double *A12Entries, *A21Entries, *A22Entries;
  int N_Active;
  double *CEntries;
  const int *CRowPtr, *CKCol;

  ARowPtr = A11->get_row_ptr();
  AKCol = A11->get_vector_columns();
  A11Entries = A11->GetEntries();
  A12Entries = A12->GetEntries();
  A21Entries = A21->GetEntries();
  A22Entries = A22->GetEntries();

  BRowPtr = B1->get_row_ptr();
  BKCol = B1->get_vector_columns();

  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

  BTRowPtr = B1T->get_row_ptr();
  BTKCol = B1T->get_vector_columns();

  B1TEntries = B1T->GetEntries();
  B2TEntries = B2T->GetEntries();

  CRowPtr = C->get_row_ptr();
  CKCol = C->get_vector_columns();
  CEntries = C->GetEntries();

  N_UDOF = A11->get_n_rows();
  N_PDOF = B1->get_n_rows();

  u1 = x;
  u2 = u1+N_UDOF;
  p  = u2+N_UDOF;

  v1 = b;
  v2 = v1+N_UDOF;
  q  = v2+N_UDOF;

  r1 = r;
  r2 = r1+N_UDOF;
  r3 = r2+N_UDOF;

  N_Active = A11->get_n_active_rows();

  j = ARowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = v1[i];
    t = v2[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A12Entries[j];
      value2 = A21Entries[j];
      value3 = A22Entries[j];
      s -= value * u1[index] + value1 * u2[index];
      t -= value2* u1[index] + value3 * u2[index];
    }
    r1[i] = s;
    r2[i] = t;
  } // endfor i

  j = ARowPtr[N_Active];
  for(i=N_Active;i<N_UDOF;i++)
  {
    s = v1[i];
    t = v2[i];
    k = ARowPtr[i+1];
    for(;j<k;j++)
    {
      index = AKCol[j];
      value = A11Entries[j];
      value1 = A22Entries[j];
      s -= value * u1[index];
      t -= value1 * u2[index];
    }
    r1[i] = s;
    r2[i] = t;
  } // endfor i

  j = BRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = q[i];
    k = BRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BKCol[j];
      value1 = B1Entries[j];
      value2 = B2Entries[j];
      s -= value1 * u1[index] + value2 * u2[index];
    } // endfor j
    r3[i] = s;
  } // endfor i

  j = BTRowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = 0;
    t = 0;
    k = BTRowPtr[i+1];
    for(;j<k;j++)
    {
      index = BTKCol[j];
      value1 = B1TEntries[j];
      value2 = B2TEntries[j];
      value = p[index];
      s += value1 * value;
      t += value2 * value;
    }
    r1[i] -= s;
    r2[i] -= t;
  } // endfor i

  j = CRowPtr[0];
  for(i=0;i<N_PDOF;i++)
  {
    s = 0;
    k = CRowPtr[i+1];
    for(;j<k;j++)
    {
      s += CEntries[j] * p[CKCol[j]]; // plus is right sign
    }
    r3[i] -= s;
  } // endfor i
}

// Disabled due to use of deprecated parameter.
//void Defect_NSE4C(TSquareMatrix **A, TMatrix **B, double *x, double *b, double *r)
//{
//  int N_UDOF,N_PDOF;
//
//  CoupledDefect(A[0], A[1], A[2], A[3], B[0], B[1], B[2], B[3], B[4], x, b, r);
//  N_UDOF = A[0]->get_n_rows();
//  N_PDOF = B[0]->get_n_rows();
//  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) //deprecated!
//    IntoL20Vector2D(r+2*N_UDOF, N_PDOF,TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE);
//  return;
//}

double UltraLocalErrorDivergence(TFEFunction2D *uh1, TFEFunction2D *uh2,
                       DoubleFunct2D *ExactU1, DoubleFunct2D *ExactU2,
                       double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l;
  int N_Cells;
  int CoarseOrder;
  int N_CoarseDOF, N_DOF;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double s;
  int i1, i2;
  double hK;
  double *Values1, *Values2;
  double error, locerror;
  double exactval1[4], exactval2[4];
  double div;

  error = 0.0;

  auto fespace = uh1->GetFESpace2D();

  auto Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);

  for(i=0;i<N_Cells;i++)
  {
    locerror = 0.0;
    auto cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    auto DOF = fespace->GetGlobalDOF(i);

    auto CurrentElement = fespace->get_fe(i);

    auto BF = CurrentElement.GetBaseFunct();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    const FiniteElement CoarseElement(GetElement2D(cell, CoarseOrder));
    auto CoarseBF = CoarseElement.GetBaseFunct();

    // quadrature formula on cell
    std::vector<const FiniteElement*> UsedElements{{&CoarseElement, &CurrentElement}};
    FEDatabase::GetOrig(UsedElements, Coll, cell, SecondDer, qf_ref, qf_orig);
    N_Points = qf_orig.GetN_QuadPoints();

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = FEDatabase::GetOrigElementValues(*CoarseBF,
                                                    MultiIndex2D::D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = qf_orig.get_weight(j);
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && std::abs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*sizeof(double));
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    Values1 = uh1->GetValues();
    Values2 = uh2->GetValues();

    PCValues = FEDatabase::GetOrigElementValues(*CoarseBF, MultiIndex2D::D00);

    ChildValuesX = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D10);
    ChildValuesY = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D01);

    // only one right-hand side (div (u-uh))
    memset(H, 0, N_CoarseDOF*1*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];

      // calculate x derivative of uh1 and y derivative of uh2 in this quadrature point
      valx = 0.0;
      valy = 0.0;
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        valx += ChildValueX[k]*Values1[l];
        valy += ChildValueY[k]*Values2[l];
      }

      // get x derivative of u1 and y derivative of u2
      auto p = qf_orig.get_point(j);
      ExactU1(p.x, p.y, exactval1);
      ExactU2(p.x, p.y, exactval2);

      valx -= exactval1[1];
      valy -= exactval2[2];

      div = valx+valy;

      w = qf_orig.get_weight(j);
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        H[k  ] += val*div;
      } // end for k

      // id-id term
      locerror += w*div*div;

    } // end for j
    memcpy(P, H, N_CoarseDOF*1*sizeof(double));

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 1, 1);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    s = 0;
    for(i1=0;i1<N_CoarseDOF;i1++)
      for(i2=0;i2<N_CoarseDOF;i2++)
        s += Gsave[i1*N_CoarseDOF+i2]*H[i1  ]*H[i2  ];
    locerror += s;

    // id-proj coupling
    s = 0;
    for(i2=0;i2<N_CoarseDOF;i2++)
    {
      s += P[i2  ] * H[i2  ];
    }
    for(i1=0;i1<N_CoarseDOF;i1++)
    {
      s += P[i1  ] * H[i1  ];
    }
    locerror -= s;

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    error += lpcoeff*std::pow(hK,lpexponent)*locerror;
  } // endfor i

  return error;
} // UltraLocalError

double UltraLocalErrorStreamline(TFEFunction2D *uh, DoubleFunct2D *ExactU,
                       TFEFunction2D *b1, TFEFunction2D *b2,
                       double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l;
  int N_Cells;
  int CoarseOrder;
  int N_CoarseDOF, N_DOF;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValues, *ChildValue;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double s;
  int i1, i2;
  double hK;
  double *Values;
  double *BValues1, *BValues2;
  double error, locerror;
  double exactval[4];
  double valb1, valb2, StreamlineDerivative;

  error = 0.0;

  auto fespace = uh->GetFESpace2D();

  auto Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);
  
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);

  for(i=0;i<N_Cells;i++)
  {
    locerror = 0.0;
    auto cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    auto DOF = fespace->GetGlobalDOF(i);

    auto CurrentElement = fespace->get_fe(i);

    auto BF = CurrentElement.GetBaseFunct();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    const FiniteElement CoarseElement(GetElement2D(cell, CoarseOrder));
    auto CoarseBF = CoarseElement.GetBaseFunct();

    // quadrature formula on cell
    std::vector<const FiniteElement*> UsedElements{{&CoarseElement, &CurrentElement}};
    FEDatabase::GetOrig(UsedElements, Coll, cell, SecondDer, qf_ref, qf_orig);
    N_Points = qf_orig.GetN_QuadPoints();

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = FEDatabase::GetOrigElementValues(*CoarseBF,
                                                    MultiIndex2D::D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = qf_orig.get_weight(j);
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && std::abs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*sizeof(double));
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    Values = uh->GetValues();

    BValues1 = b1->GetValues();
    BValues2 = b2->GetValues();
    
    PCValues = FEDatabase::GetOrigElementValues(*CoarseBF, MultiIndex2D::D00);

    ChildValuesX = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D10);
    ChildValuesY = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D01);
    ChildValues  = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D00);

    // only one right-hand side (b.grad(ui))
    memset(H, 0, N_CoarseDOF*1*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValue  = ChildValues[j];

      // calculate derivative of uh and components of b in this quadrature point
      valx = 0.0;
      valy = 0.0;
      valb1 = 0.0;
      valb2 = 0.0;
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        valb1 += ChildValue[k]*BValues1[l];
        valb2 += ChildValue[k]*BValues2[l];
        valx += ChildValueX[k]*Values[l];
        valy += ChildValueY[k]*Values[l];
      }

      // get gradient of exact u
      auto p = qf_orig.get_point(j);
      ExactU(p.x, p.y, exactval);

      valx -= exactval[1];
      valy -= exactval[2];

      StreamlineDerivative = valb1*valx + valb2*valy;

      w = qf_orig.get_weight(j);
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        H[k  ] += val*StreamlineDerivative;
      } // end for k

      // id-id term
      locerror += w*StreamlineDerivative*StreamlineDerivative;

    } // end for j
    memcpy(P, H, N_CoarseDOF*1*sizeof(double));

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 1, 1);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    s = 0;
    for(i1=0;i1<N_CoarseDOF;i1++)
      for(i2=0;i2<N_CoarseDOF;i2++)
        s += Gsave[i1*N_CoarseDOF+i2]*H[i1  ]*H[i2  ];
    locerror += s;

    // id-proj coupling
    s = 0;
    for(i2=0;i2<N_CoarseDOF;i2++)
    {
      s += P[i2  ] * H[i2  ];
    }
    for(i1=0;i1<N_CoarseDOF;i1++)
    {
      s += P[i1  ] * H[i1  ];
    }
    locerror -= s;

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    error += lpcoeff*std::pow(hK,lpexponent)*locerror;
  } // endfor i

  return error;
} // UltraLocalError

double UltraLocalErrorStreamlinePWConst(TFEFunction2D *uh, DoubleFunct2D *ExactU,
                       TFEFunction2D *b1, TFEFunction2D *b2,
                       double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l;
  int N_Cells;
  int CoarseOrder;
  int N_CoarseDOF, N_DOF;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValues, *ChildValue;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double s;
  int i1, i2;
  double hK;
  double *Values;
  double *BValues1, *BValues2;
  double error, locerror;
  double exactval[4];
  double valb1, valb2, StreamlineDerivative;

  error = 0.0;

  auto fespace = uh->GetFESpace2D();

  auto Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);
  
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);

  for(i=0;i<N_Cells;i++)
  {
    locerror = 0.0;
    auto cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    auto DOF = fespace->GetGlobalDOF(i);

    auto CurrentElement = fespace->get_fe(i);

    auto BF = CurrentElement.GetBaseFunct();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    const FiniteElement CoarseElement(GetElement2D(cell, CoarseOrder));
    auto CoarseBF = CoarseElement.GetBaseFunct();

    // quadrature formula on cell
    std::vector<const FiniteElement*> UsedElements{{&CoarseElement, &CurrentElement}};
    FEDatabase::GetOrig(UsedElements, Coll, cell, SecondDer, qf_ref,
                           qf_orig);
    N_Points = qf_orig.GetN_QuadPoints();

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = FEDatabase::GetOrigElementValues(*CoarseBF, 
MultiIndex2D::D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = qf_orig.get_weight(j);
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && std::abs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*sizeof(double));
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    Values = uh->GetValues();

    BValues1 = b1->GetValues();
    BValues2 = b2->GetValues();
    
    PCValues = FEDatabase::GetOrigElementValues(*CoarseBF, MultiIndex2D::D00);

    ChildValuesX = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D10);
    ChildValuesY = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D01);
    ChildValues  = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D00);

    // only one right-hand side (b.grad(ui))
    memset(H, 0, N_CoarseDOF*1*sizeof(double));

    // calculate pw constant approximation of velocity field (uh1, uh2)
    val  = 0.0;
    valx = 0.0;
    valy = 0.0;
    for(j=0;j<N_Points;j++)
    {
      ChildValue  = ChildValues[j];
      w = qf_orig.get_weight(j);
      // compute components of uh in j
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        val  += w*ChildValue[k]*1;
        valx += w*ChildValue[k]*BValues1[l];
        valy += w*ChildValue[k]*BValues2[l];
      }
    }
    valb1 = valx/val;
    valb2 = valy/val;

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValue  = ChildValues[j];

      // calculate derivative of uh and components of b in this quadrature point
      valx = 0.0;
      valy = 0.0;
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        valx += ChildValueX[k]*Values[l];
        valy += ChildValueY[k]*Values[l];
      }

      // get gradient of exact u
      auto p = qf_orig.get_point(j);
      ExactU(p.x, p.y, exactval);

      valx -= exactval[1];
      valy -= exactval[2];

      StreamlineDerivative = valb1*valx + valb2*valy;

      w = qf_orig.get_weight(j);
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        H[k  ] += val*StreamlineDerivative;
      } // end for k

      // id-id term
      locerror += w*StreamlineDerivative*StreamlineDerivative;

    } // end for j
    memcpy(P, H, N_CoarseDOF*1*sizeof(double));

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 1, 1);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    s = 0;
    for(i1=0;i1<N_CoarseDOF;i1++)
      for(i2=0;i2<N_CoarseDOF;i2++)
        s += Gsave[i1*N_CoarseDOF+i2]*H[i1  ]*H[i2  ];
    locerror += s;

    // id-proj coupling
    s = 0;
    for(i2=0;i2<N_CoarseDOF;i2++)
    {
      s += P[i2  ] * H[i2  ];
    }
    for(i1=0;i1<N_CoarseDOF;i1++)
    {
      s += P[i1  ] * H[i1  ];
    }
    locerror -= s;

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    error += lpcoeff*std::pow(hK,lpexponent)*locerror;
  } // endfor i

  return error;
} // UltraLocalErrorPWConst


//**************************************************************
//  UltraLocalProjectionStreamlinePLaplacian
//  ultra LPS for scalar equation
//  projection of streamline derivative
//  p-Laplacian term can be added
//**************************************************************

void UltraLocalProjectionStreamlinePLaplacian(TSquareMatrix2D* A, 
                TFEFunction2D *uh,
                const CoeffFct2D& Coeff)
{
  int i,j,k,l,m;
  int N_Cells;
  int CoarseOrder;
  int N_CoarseDOF, N_DOF;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double rhs_l2_proj[MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s, loc_proj_L2;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end, ij, N_Edges;
  const int *RowPtr, *KCol;
  double *Entries;
  double *Values;
  double *coeffs, *params, sx, sy, crosswind_uh[MaxN_QuadPoints_2D], loc_proj;
  double BValue[MaxN_BaseFunctions2D];
  int OrderDiff;
  double lpcoeff, lpexponent, stab_coeff, norm_b, lpcoeff_crosswind, lpexponent_crosswind;
  
  Output::print<2>("LPS streamline");

  coeffs = new double[20];
  params = new double[10];
  memset(params, 0, 10 * sizeof(double));

  lpcoeff = TDatabase::ParamDB->LP_STREAMLINE_COEFF;
  lpexponent = TDatabase::ParamDB->LP_STREAMLINE_EXPONENT;
  OrderDiff = TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE;
  lpcoeff_crosswind = TDatabase::ParamDB->LP_CROSSWIND_COEFF;
  lpexponent_crosswind = TDatabase::ParamDB->LP_CROSSWIND_EXPONENT;
  
  // get fespace and matrices
  auto fespace = A->GetFESpace2D();
  ActiveBound = fespace->get_n_active();
  RowPtr = A->get_row_ptr();
  KCol = A->get_vector_columns();
  Entries = A->GetEntries();

  // get values of fe function if available
  if (uh != nullptr)
  {
      Values = uh->GetValues();     
  }
  else
  {
      if (TDatabase::ParamDB->LP_CROSSWIND)
      {
    ErrThrow("for crosswind LPS the current finite element function has to be given to the routine");
      }
  }

  // get collection
  auto Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);
  
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);
  
  // loop over the cells
  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);
    // number of edges
    N_Edges = cell->GetN_Edges();
    // get diameter of the cell
    switch (TDatabase::ParamDB->CELL_MEASURE)
    {
  case 0: // diameter
      hK = cell->GetDiameter();
      break;
  case 1: // with reference map
      hK = cell->GetLengthWithReferenceMap();
      break;
  case 2: // shortest edge
      hK = cell->GetShortestEdge();
      break;
  case 3: // measure
      hK = cell->GetMeasure();
      hK = std::pow(hK,1.0/3.0);
      break;
  case 4: // mesh cell in convection direction
            TDatabase::ParamDB->INTERNAL_HK_CONVECTION = -1;
            sx = sy = 0;
      for (ij=0;ij<N_Edges;ij++)
      {
    TDatabase::ParamDB->INTERNAL_VERTEX_X[ij] = cell->GetVertex(ij)->GetX();
    sx += TDatabase::ParamDB->INTERNAL_VERTEX_X[ij];
    TDatabase::ParamDB->INTERNAL_VERTEX_Y[ij] = cell->GetVertex(ij)->GetY();
    sy += TDatabase::ParamDB->INTERNAL_VERTEX_Y[ij];
      }
      if (N_Edges==3)
    TDatabase::ParamDB->INTERNAL_VERTEX_X[3] = -4711;
      // center of mesh cell
      sx /= N_Edges;
      sy /= N_Edges;
      hK = cell->GetDiameter(); 
      // get coefficients in center of mesh cell 
      Coeff(1, &sx ,&sy, &params, &coeffs);
      hK = Mesh_size_in_convection_direction<2>(hK, {{coeffs[1], coeffs[2]}});
      break;
  default: // diameter
      hK = cell->GetDiameter();
      break;
    }
    
    // get finite element in the cell
    auto CurrentElement = fespace->get_fe(i);
    // get basis functions of the finite element
    auto BF = CurrentElement.GetBaseFunct();
    N_DOF = BF->GetDimension();

    // compute order of the local coarse space
    CoarseOrder = BF->GetAccuracy() - OrderDiff;

    const FiniteElement CoarseElement(GetElement2D(cell, CoarseOrder));
    auto CoarseBF = CoarseElement.GetBaseFunct();

    // quadrature formula on cell
    std::vector<const FiniteElement*> UsedElements{{&CoarseElement, &CurrentElement}};
    FEDatabase::GetOrig(UsedElements, Coll, cell, SecondDer, qf_ref, qf_orig);
    N_Points = qf_orig.GetN_QuadPoints();
    std::vector<double> X(N_Points);
    std::vector<double> Y(N_Points);
    for(int j=0;j<N_Points;j++)
    {
      auto p = qf_orig.get_point(j);
      X[i] = p.x;
      Y[i] = p.y;
    }

    // number of dof in the local coarse space
    N_CoarseDOF = CoarseBF->GetDimension();
    // get function values for the basis functions of the coarse space
    CoarseValues = FEDatabase::GetOrigElementValues(*CoarseBF,
                                                    MultiIndex2D::D00);
    // initialize array G, stores mass matrix of coarse space
    memset(G, 0, N_CoarseDOF*N_CoarseDOF*sizeof(double));
    // loop over the quadrature points
    for(j=0;j<N_Points;j++)
    {
      // values of the coarse basis functions in the quad points
      CoarseValue = CoarseValues[j];
      // factor for numerical quadrature
      w = qf_orig.get_weight(j);
      // loop over the basis functions of the coarse space
      for(k=0;k<N_CoarseDOF;k++)
      {
    // first factor of integrand
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
      // update integral 
           G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && std::abs(G[0])<1e-10)
      G[0] = 1;
    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*sizeof(double));

    // degrees of freedom of the fine space 
    auto DOF = fespace->GetGlobalDOF(i);
    // this should be the same as CoarseValues ???    
    //PCValues = FEDatabase::GetOrigElementValues(*CoarseBF, D00);
    PCValues = CoarseValues;
    // get derivatives of fe functions of approximation space
    ChildValuesX = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D10);
    ChildValuesY = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D01);
    //double **ChildValues  = FEDatabase::GetOrigElementValues(*BF, D00);
    // initialize array H, holds mixed products of coarse basis 
    // functions and streamline derivatives of fine basis functions
    memset(H, 0, N_CoarseDOF*N_DOF*sizeof(double));
    // initialize array LocMatrix, holds products of derivatives
    // of fine basis functions
    memset(LocMatrix, 0, N_DOF*N_DOF*sizeof(double));
    // loop over the quadrature points
    for(j=0;j<N_Points;j++)
    {
      // get all values in the quad point j 
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      // factor in the integral
      w = qf_orig.get_weight(j);
      // values of convection in this quadrature point
      Coeff(1, &X[j] , &Y[j], &params, &coeffs);
      // compute streamline derivative
      for(k=0;k<N_DOF;k++)
      {
    BValue[k] = coeffs[1]*ChildValueX[k] + coeffs[2]*ChildValueY[k];
    //BValue[k] = -coeffs[2]*ChildValueX[k] + coeffs[1]*ChildValueY[k];
      }
      // loop over the basis functions of the coarse space
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
  // compute products of coarse basis function and
  // streamline derivatives of fine basis functions
        for(l=0;l<N_DOF;l++)
        {
          H[k*N_DOF+l      ] += val*BValue[l];
        } // end for l
      } // end for kcoeffs[2]

      // fine-fine couplings
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
      LocMatrix[k*N_DOF+l] += w*BValue[k]*BValue[l];
        }
      }
    } // end for j
    // save H for later use
    memcpy(P, H, N_CoarseDOF*N_DOF*sizeof(double));
     // solve G * X = H, solution X stored on H
    // the right hand side and the solution are stored column wise
    // right hand side: a finecoeffs[2] fct. tested with all coarse fcts. 
    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, N_DOF, N_DOF);

    // update LocMatrix
    // proj-proj coupling (coarse-coarse coupling)
    // l - test function index
    for(l=0;l<N_DOF;l++)
    {
      // m - ansatz function index
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2] * H[i1*N_DOF+l] * H[i2*N_DOF+m];
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling (fine-coarse couplings)
    // l - test function index
    for(l=0;l<N_DOF;l++)
    {
      // m - ansatz function index
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*N_DOF+l      ] * H[i2*N_DOF+m      ];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*N_DOF+m      ] * H[i1*N_DOF+l      ];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l
    // compute stabilzation coefficient
    switch(TDatabase::ParamDB->LP_COEFF_TYPE)
    {
  case 0:
      stab_coeff = lpcoeff*std::pow(hK,lpexponent);
      break;
  case 1:
      // h^2/epsilon
      stab_coeff = hK*hK/coeffs[0];
      norm_b = std::sqrt(coeffs[1]*coeffs[1]+coeffs[2]*coeffs[2]);
      if (norm_b > 1e-10)
      {
    // h/|b|
    val = hK/norm_b;
    if (val < stab_coeff)
        stab_coeff = val;
      }
      stab_coeff *= lpcoeff;
      break;
  default:
      stab_coeff = lpcoeff*std::pow(hK,lpexponent);
      break;
    }
 
    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if(dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
    // parameter is lpcoeff*std::pow(hK,lpexponent)
               Entries[p] += stab_coeff*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l

    //***************************************************
    // LPS with SOLD term
    //***************************************************
    if (TDatabase::ParamDB->LP_CROSSWIND)
    {
  // recover mass matrix of coarse space
  memcpy(G, Gsave, N_CoarseDOF*N_CoarseDOF*sizeof(double));
  memset(rhs_l2_proj, 0, N_CoarseDOF*sizeof(double));
  // compute L2 projection of crosswind derivative
  // this gives just an array (one value for each quad point)
  // mass matrix of coarse space already computed: G
  // compute right hand side
  // loop over the quadrature points
  for(j=0;j<N_Points;j++)
  {
      // get all values in the quad point j 
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      // values of the coarse basis functions in the quad points
      CoarseValue = CoarseValues[j];
      // factor for numerical quadrature
      w = qf_orig.get_weight(j);
      // compute crosswind derivative in this point
      // calculate gradient of discrete uh in this quadrature point
      valx = 0.0;
      valy = 0.0;
      for(k=0;k<N_DOF;k++)
      {
    l = DOF[k];
    val = Values[l];
    valx += ChildValueX[k]*val;
    valy += ChildValueY[k]*val;
      }
      // values of convection in this quadrature point
      Coeff(1, &X[j] , &Y[j], &params, &coeffs);
      
      // calculate crosswind derivative in this quadrature point
      crosswind_uh[j] = -coeffs[2] * valx + coeffs[1] * valy;
      //Output::print(crosswind_uh[j]);
      w *= crosswind_uh[j];
      // loop over the test functions of the coarse space
      // right hand side of linear system for local projection
      for(k=0;k<N_CoarseDOF;k++)
      {
    rhs_l2_proj[k] += w*CoarseValue[k];
      } // end for k
      //Output::print(rhs_l2_proj[0] << " :: ");
  } // end for j  
  // solution is on rhs_l2_proj
  // these are coefficients wrt to coarse space basis
  //Output::print("G ", G[0], " ", rhs_l2_proj[0]);
  SolveLinearSystemLapack(G, rhs_l2_proj, N_CoarseDOF, N_CoarseDOF);
  //Output::print("G ", G[0], " ", rhs_l2_proj[0]);
  // recover mass matrix of coarse space
  memcpy(G, Gsave, N_CoarseDOF*N_CoarseDOF*sizeof(double));
 
  // initialize array H, holds mixed products of coarse basis 
  // functions and crosswind derivatives of fine basis functions
  memset(H, 0, N_CoarseDOF*N_DOF*sizeof(double));
  // initialize array LocMatrix, holds products of derivatives
  // of fine basis functions
  memset(LocMatrix, 0, N_DOF*N_DOF*sizeof(double));

  // loop over the quadrature points
  // this computes the projections of the crosswind derivatives 
  // for each basis function of the fine space
  for(j=0;j<N_Points;j++)
  {
      // get all values in the quad point j 
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      // factor in the integral
      w = qf_orig.get_weight(j);
      // values of convection in this quadrature point
      Coeff(1, &X[j] , &Y[j], &params, &coeffs);
      // compute crosswind derivative
      for(k=0;k<N_DOF;k++)
      {
    BValue[k] = -coeffs[2]*ChildValueX[k] + coeffs[1]*ChildValueY[k];
        }
      // loop over the basis functions of the coarse space
      for(k=0;k<N_CoarseDOF;k++)
      {
    val = w*PCValue[k];
    // compute products of coarse basis function and
    // crosswind derivatives of fine basis functions
    for(l=0;l<N_DOF;l++)
    {
        H[k*N_DOF+l] += val*BValue[l];
    } // end for l
      } // end for k
  } // end for j

  // solve G * X = H, solution X stored on H
  // the right hand side and the solution are stored column wise
  // right hand side: a fine fct. tested with all coarse fcts. 
  SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, N_DOF, N_DOF);

  // compute L^2 norm of fluctuations 
  if (TDatabase::ParamDB->LP_CROSSWIND_COEFF_TYPE >= 200)
  {
      loc_proj_L2 = 0.0;
      for(j=0;j<N_Points;j++)
      {
    // get all values in the quad point j 
    PCValue = PCValues[j];
    ChildValueX = ChildValuesX[j];
    ChildValueY = ChildValuesY[j];
    // factor in the integral
    w = qf_orig.get_weight(j);
    // values of convection in this quadrature point
    Coeff(1, &X[j] , &Y[j], &params, &coeffs);
    // compute crosswind derivative
    for(k=0;k<N_DOF;k++)
    {
        BValue[k] = -coeffs[2]*ChildValueX[k] + coeffs[1]*ChildValueY[k];
    }
    // compute projection of crosswind derivative for basis functions
    for(k=0;k<N_DOF;k++)
    {
        P[k] = 0;
        for (l=0; l < N_CoarseDOF; l++)
      P[k] += H[l*N_DOF+k] * PCValue[l];
    }
    // calculate projection of crosswind derivative discrete uh in this quadrature point
    val = 0.0;
    for(k=0;k<N_CoarseDOF;k++)
    {
        val += rhs_l2_proj[k] * PCValue[k]; 
    }
    // local projection of crosswind derivative of current solution in this quadrature point
    // this is the coefficient
    loc_proj = val - crosswind_uh[j];
    // update norm
    loc_proj_L2 += w * loc_proj * loc_proj;
      } // end for j
      // L2 norm and scale such that independent of h
      loc_proj_L2 = std::sqrt(loc_proj_L2/cell->GetMeasure());
  }

  // loop over the quadrature points
  // this computes the LPS term
  // for each basis function of the fine space
  for(j=0;j<N_Points;j++)
  {
      // get all values in the quad point j 
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      // factor in the integral
      w = qf_orig.get_weight(j);
      // values of convection in this quadrature point
      Coeff(1, &X[j] , &Y[j], &params, &coeffs);
      // compute crosswind derivative
      for(k=0;k<N_DOF;k++)
      {
    BValue[k] = -coeffs[2]*ChildValueX[k] + coeffs[1]*ChildValueY[k];
      }
      // compute projection of crosswind derivative for basis functions
      for(k=0;k<N_DOF;k++)
      {
    P[k] = 0;
    for (l=0; l < N_CoarseDOF; l++)
        P[k] += H[l*N_DOF+k] * PCValue[l];
      }
      // calculate projection of crosswind derivative discrete uh in this quadrature point
      if (TDatabase::ParamDB->LP_CROSSWIND_COEFF_TYPE < 100)
      {
    val = 0.0;
    for(k=0;k<N_CoarseDOF;k++)
    {
        val += rhs_l2_proj[k] * PCValue[k]; 
    }
    // local projection of crosswind derivative of current solution in this quadrature point
    // this is the coefficient
    loc_proj = std::abs(val - crosswind_uh[j]);
      }
      // linear case
      if ((TDatabase::ParamDB->LP_CROSSWIND_COEFF_TYPE >= 100) && (TDatabase::ParamDB->LP_CROSSWIND_COEFF_TYPE <200))
    loc_proj = 1.0;
      // local average
      if (TDatabase::ParamDB->LP_CROSSWIND_COEFF_TYPE >= 200)
      {
    loc_proj = loc_proj_L2;
      }

      // update factor in integral
      w *= loc_proj;
      // update local matrix
      for(k=0;k<N_DOF;k++)
      {
    for(l=0;l<N_DOF;l++)
    {
        // fine-fine
        val = BValue[k]*BValue[l];
        // coarse-fine
        val -= BValue[k]*P[l];
        val -= BValue[l]*P[k];
        // coarse-coarse
        val += P[l]*P[k];       
        LocMatrix[k*N_DOF+l] += w*val;
    }
      }
  } // end for j


  // compute stabilzation coefficient
  switch(TDatabase::ParamDB->LP_CROSSWIND_COEFF_TYPE)
  {
      case 0:
      case 100:
      case 200:
    stab_coeff = lpcoeff_crosswind*std::pow(hK,lpexponent_crosswind);
    break;
      case 1:
      case 101:
      case 201:
    // norm squared of convection
    norm_b = coeffs[1]*coeffs[1]+coeffs[2]*coeffs[2];
    if (norm_b > 1e-20)
    {
        // h/|b|
        stab_coeff = lpcoeff_crosswind*hK/norm_b;
    }
    else
        stab_coeff = 0.0;
    break;
      default:
    stab_coeff = lpcoeff_crosswind*std::pow(hK,lpexponent_crosswind);
    break;
  }
    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if(dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
    // parameter is lpcoeff*std::pow(hK,lpexponent)
               Entries[p] += stab_coeff*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
    } // end LP_CROSSWIND

  } // endfor i
  delete params;
  delete coeffs;
  Output::print<2>("LPS streamline done");
} // UltraLocalProjectionStreamlinePLaplacian


//**************************************************************
// computes the local projection to Q0 on a coarse grid
// uh      -  current fe function
// uh_proj -  fe function for projection, pw constant
//**************************************************************

void LocalProjectionCoarseGridQ0(TFEFunction2D *uh,
         TFEFunction2D *uh_proj,
               const CoeffFct2D& Coeff,
               int convection_flag)
{
  int i, j, iq, index;
  int N_Cells, N_Edges;
  double sx, sy, b1, b2, area, detJK, weight_det, value;
  double *coeffs, *params, *Values_coarse,  x[4], y[4], val[4];
  double gauss2_x[4]=
  {
    -0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746,
    -0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746
  };
  double gauss2_y[4]=
  {
    -0.57735026918962576450914878050195746,
    -0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746
  };
  const TBaseCell *cell, *parent_cell, *child_cell;

  Output::print("compute local projection to Q0 on coarse grid");
  coeffs = new double[20];
  params = new double[10];
  memset(params, 0, 10 * sizeof(double));

  auto fespace_fine = uh->GetFESpace2D();
  //const TFESpace2D *fespace_coarse = uh_proj->GetFESpace2D();
  // get collection, same for fine and coarse fe space
  auto Coll = fespace_fine->GetCollection();
  N_Cells = Coll->GetN_Cells();
  
  // get values of fe function 
  //double *Values_fine = uh->GetValues();     
  Values_coarse = uh_proj->GetValues();     
  // pw constant space
  memset(Values_coarse, 0.0, N_Cells * sizeof(double));
  
  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  // loop over the cells
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);   
    // number of edges
    N_Edges = cell->GetN_Edges();
    // get barycenter
    sx = sy = 0;
    for (j=0;j<N_Edges;j++)
    {
      x[j] = cell->GetVertex(j)->GetX(); 
      sx += x[j];
      y[j] = cell->GetVertex(j)->GetY();
      sy += y[j];
    }
    sx /= N_Edges;
    sy /= N_Edges;
    // area of parallelogramm with vector product
    area = std::abs((x[1]-x[0])*(y[3]-y[0]) - (x[3]-x[0])*(y[1]-y[0]));
    // functional determinant
    detJK = area/4.0;
    // get coefficients in center of mesh cell 
    Coeff(1, &sx ,&sy, &params, &coeffs);
    // take the desired direction
    switch(convection_flag)
    {
      case 0:
      // convection
        b1 = coeffs[1];
        b2 = coeffs[2];
        break;
      case 1:
      // normalized crosswind convection
        b2 = std::sqrt(coeffs[1]*coeffs[1] + coeffs[2] * coeffs[2]);
  b1 = -coeffs[2]/b2;
  b2 = coeffs[1]/b2;
  break;
    }
    // compute integral of b times grad_uh
    // loop over the quadrature points
     for (iq = 0;iq < 4; iq++)
    {
       sx = (x[1]-x[0])* gauss2_x[iq] + (x[3]-x[0])*gauss2_y[iq] + x[1] + x[3];
       sx /= 2.0;
       sy = (y[1]-y[0])* gauss2_x[iq] + (y[3]-y[0])*gauss2_y[iq] + y[1] + y[3];
       sy /= 2.0;
       uh->FindGradientLocal(cell,i,sx,sy,val);
       value = b1 * val[1] + b2 * val[2];
       weight_det = detJK;
       Values_coarse[i] += value * weight_det;
    }   
  }
   
  // compute the mean values on the macros
    // loop over the cells
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    // cell not yet child of a macro cell
    if (cell->GetClipBoard()>=0)
    {
      parent_cell = cell->GetParent();
      if (parent_cell == nullptr)
      {
        Output::print("no coarse grid to project !!!");
        delete params;
        delete coeffs;
  return;
      }
      area = 0.0;
      value = 0.0;
      for (j=0;j<4;j++)
      {
  child_cell = parent_cell->GetChild(j);
  // get index of pw constant space
  index = child_cell->GetClipBoard();
  // integral of projected quantity
  value += Values_coarse[index];
  // area of macro cell
  area += child_cell->GetMeasure();
      }
      value /= area;
      for (j=0;j<4;j++)
      {
  child_cell = parent_cell->GetChild(j);
  // get index of pw constant space
  index = child_cell->GetClipBoard();
        // projection
  Values_coarse[index] = value;
  child_cell->SetClipBoard(-1);
      }
    }
  }

  delete[] params;
  delete[] coeffs;
}


void LocalProjectionCrossWindCoarseGridQ0(TDomain *, int,
            TFEFunction2D *uh,
            TFEFunction2D *uh_proj,
            const CoeffFct2D& Coeff,
            double *rhs,
                        int convection_flag) 
{

  ErrThrow("This method was dependent on two old global Database parameters ",
           "which were removed. If you want to use it, take care of the FIXMEs.");

  int i, j, k, iq, N_Cells, N_Edges, N_Unknowns, index, dof_index;
  double  sx, sy, b1, b2, area, detJK, weight_det, value, value1, value_proj, hK;
  double *coeffs, *params, *Values_coarse,  x[4], y[4], val[4];
  double *sol, *sol_copy, *proj_test, area_coarse, norm_b, stab_coeff_cw;
  double lpcoeff_crosswind = TDatabase::ParamDB->LP_CROSSWIND_COEFF;
  double gauss2_x[4]=
  {
    -0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746,
    -0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746
  };
  double gauss2_y[4]=
  {
    -0.57735026918962576450914878050195746,
    -0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746
  };
  TCollection *coll_coarse, *coll;
  const TBaseCell *cell, *child_cell;
 
   Output::print("update rhs of crosswind local projection to Q0 on coarse grid");
  // get coarse grid
  // FIXME that parameter is gone! coll_coarse=Domain->GetCollection(It_EQ, mg_level+TDatabase::ParamDB->SC_COARSEST_LEVEL_SCALAR-1);
//  if (coll_coarse == nullptr)
//  {
//    Output::print("No coarse grid !!!");
//    return;
//  }
// END FIXME
   coeffs = new double[20];
  params = new double[10];
  memset(params, 0, 10 * sizeof(double));

  // get fine grid
  // FIXME that parameter is gone! coll = Domain->GetCollection(It_EQ, mg_level+TDatabase::ParamDB->SC_COARSEST_LEVEL_SCALAR);
  //N_Cells = coll->GetN_Cells();
  N_Cells = 1000; //A dummy number.
  //END FIXME
  // initialise ClipBoard for fine grid
  for(i=0;i<N_Cells;i++)
    coll->GetCell(i)->SetClipBoard(i);
  
  // fespace on the fine grid
  auto fespace = uh->GetFESpace2D();
  // copy solution
  N_Unknowns = uh->GetLength();
  sol = uh->GetValues();
  sol_copy = new double[N_Unknowns];
  memcpy(sol_copy,sol, N_Unknowns*sizeof(double));

  proj_test = new double[N_Unknowns];
  // the already computed projections of the current solution 
  Values_coarse = uh_proj->GetValues();  
  // loop over coarse grid
  N_Cells = coll_coarse->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = coll_coarse->GetCell(i);  
    area_coarse = cell->GetMeasure();
    hK = cell->GetDiameter();
    memset(proj_test,0,N_Unknowns*sizeof(double));
    for (j=0;j<4;j++)
    {
        // get child cell
  child_cell = cell->GetChild(j);
  // get degrees of freedom in the cell on the fine grid
  index = child_cell->GetClipBoard();
  auto dof = fespace->GetGlobalDOF(index);
        // number of edges
        N_Edges = child_cell->GetN_Edges();
  // for each degree of freedom
       for (k=0;k<N_Edges;k++)
       {
         x[k] = child_cell->GetVertex(k)->GetX(); 
         sx += x[k];
         y[k] = child_cell->GetVertex(k)->GetY();
         sy += y[k];
       }
       sx /= N_Edges;
       sy /= N_Edges;
       // area of parallelogramm with vector product
       area = std::abs((x[1]-x[0])*(y[3]-y[0]) - (x[3]-x[0])*(y[1]-y[0]));
       // functional determinant
       detJK = area/4.0;
       // get coefficients in center of mesh cell 
       Coeff(1, &sx ,&sy, &params, &coeffs);
       // take the desired direction
       switch(convection_flag)
       {
         case 0:
         // convection
           b1 = coeffs[1];
           b2 = coeffs[2];
           break;
        case 1:
        // normalized crosswind convection
           b2 = std::sqrt(coeffs[1]*coeffs[1] + coeffs[2] * coeffs[2]);
           b1 = -coeffs[2]/b2;
     b2 = coeffs[1]/b2;
     break;
        }
        // loop over the local dofs
  for (k=0;k<N_Edges;k++)
  {
    memset(sol,0,N_Unknowns*sizeof(double));
    dof_index = dof[k];
    sol[dof_index] = 1.0;
          // compute integral of b times grad_vh
          // loop over the quadrature points
          for (iq = 0;iq < 4; iq++)
          {
             sx = (x[1]-x[0])* gauss2_x[iq] + (x[3]-x[0])*gauss2_y[iq] + x[1] + x[3];
             sx /= 2.0;
             sy = (y[1]-y[0])* gauss2_x[iq] + (y[3]-y[0])*gauss2_y[iq] + y[1] + y[3];
             sy /= 2.0;
             uh->FindGradientLocal(child_cell,index,sx,sy,val);
             value = b1 * val[1] + b2 * val[2];
             weight_det = detJK;
             proj_test[dof_index] += value * weight_det/area_coarse;
           }   
           Output::print(dof_index, " ", proj_test[dof_index], ":");
  } 
    } // end j
    // projection of test function for macro cell i computed
    // now compute contribution to rhs
    for (j=0;j<4;j++)
    {
        // get child cell
  child_cell = cell->GetChild(j);
  // get degrees of freedom in the cell on the fine grid
  index = child_cell->GetClipBoard();
  auto dof = fespace->GetGlobalDOF(index);
        // number of edges
        N_Edges = child_cell->GetN_Edges();
  // for each degree of freedom
       for (k=0;k<N_Edges;k++)
       {
         x[k] = child_cell->GetVertex(k)->GetX(); 
         sx += x[k];
         y[k] = child_cell->GetVertex(k)->GetY();
         sy += y[k];
       }
       sx /= N_Edges;
       sy /= N_Edges;
       // area of parallelogramm with vector product
       area = std::abs((x[1]-x[0])*(y[3]-y[0]) - (x[3]-x[0])*(y[1]-y[0]));
       // functional determinant
       detJK = area/4.0;
       // get coefficients in center of mesh cell 
       Coeff(1, &sx ,&sy, &params, &coeffs);
       // take the desired direction
       switch(convection_flag)
       {
         case 0:
         // convection
           b1 = coeffs[1];
           b2 = coeffs[2];
           break;
        case 1:
        // normalized crosswind convection
           b2 = std::sqrt(coeffs[1]*coeffs[1] + coeffs[2] * coeffs[2]);
           b1 = -coeffs[2]/b2;
     b2 = coeffs[1]/b2;
     break;
        }
        // value of projection
  value_proj = Values_coarse[index];
        // loop over the local dofs
  for (k=0;k<N_Edges;k++)
  {
    dof_index = dof[k];
          // compute integral of b times grad_vh
          // loop over the quadrature points
    value1 = 0;
          for (iq = 0;iq < 4; iq++)
          {
             sx = (x[1]-x[0])* gauss2_x[iq] + (x[3]-x[0])*gauss2_y[iq] + x[1] + x[3];
             sx /= 2.0;
             sy = (y[1]-y[0])* gauss2_x[iq] + (y[3]-y[0])*gauss2_y[iq] + y[1] + y[3];
             sy /= 2.0;
             uh->FindGradientLocal(child_cell,index,sx,sy,val);
             value = b1 * val[1] + b2 * val[2] - value_proj;
       value = std::abs(value) * value * proj_test[dof_index];
             weight_det = detJK;
             value1 += value * weight_det;
           }
           norm_b = std::sqrt(b1*b1+b2*b2);
           if (norm_b > 1e-20)
           {
           // h/|b|
             stab_coeff_cw = 2*lpcoeff_crosswind*hK/norm_b;
           }
           else
             stab_coeff_cw = 0.0;
     
           rhs[dof_index] += stab_coeff_cw * value1;
     //Output::print( -hK*value1);
  } 
    } // end j
   
  }
  
  // recover solution
  memcpy(sol,sol_copy, N_Unknowns*sizeof(double));

  delete[] proj_test;
  delete[] sol_copy; 
  delete[] params;
  delete[] coeffs;

}

// Adaptive post processing basd on Friedhelm Schiweck talk at MAFELAP 09 - Sashi
void AdaptivePostProcess(TFEFunction2D *FeFunction, double *PostSol, bool)
{
  int i, j, k, l, CoarseOrder, N_Points, N_U;
  int N_Cells, N_CoarseDOF, N_DOF;
  int OrderDiff;

  double **CoarseValues, *CoarseValue, *LPS_sol;
  double val, w;
  const double *xi, *eta;
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double rhs[MaxN_BaseFunctions2D];
  double Values[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];
  double PCValues[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];
  double *Value, sol;
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_PointsForNodal2D];
  double *W, maxbubble=-1E8, minbubble=1E8;


  bool SecondDer[1] = { false };


  auto fespace = FeFunction->GetFESpace2D();
  auto coll = fespace->GetCollection();
  N_Cells = coll->GetN_Cells();
  LPS_sol = FeFunction->GetValues();
  N_U = FeFunction->GetLength();

  OrderDiff = TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE;
  W = new double[N_U];
  memset(W, 0, N_U*sizeof(double));

  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);

  for(i=0;i<N_Cells;i++)
   {
    auto cell = coll->GetCell(i);

    auto CurrentElement = fespace->get_fe(i);
    auto BF = CurrentElement.GetBaseFunct();
    //BaseFunct2D BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    const FiniteElement CoarseElement(GetElement2D(cell, CoarseOrder));
    auto CoarseBF = CoarseElement.GetBaseFunct();
    N_CoarseDOF = CoarseBF->GetDimension();
    std::vector<const FiniteElement*> UsedElements(1, &CoarseElement);
    // quadrature formula on cell
    FEDatabase::GetOrig(UsedElements, coll, cell, SecondDer, qf_ref, qf_orig);
    N_Points = qf_orig.GetN_QuadPoints();

    CoarseValues = FEDatabase::GetOrigElementValues(*CoarseBF,
                                                    MultiIndex2D::D00);

    auto DOF = fespace->GetGlobalDOF(i);

    for(j=0;j<N_Points;j++)
    {
      auto p = qf_ref.get_point(j);
      BF->GetDerivatives(MultiIndex2D::D00, p.x, p.y, Values[j]);
    }


    memset(G, 0, N_CoarseDOF*N_CoarseDOF*sizeof(double));
    memset(rhs, 0, N_CoarseDOF*sizeof(double));

    for(j=0;j<N_Points;j++)
     {
      CoarseValue = CoarseValues[j];
      Value = Values[j];
      w = qf_orig.get_weight(j);

//    find the lps solution at this quadrature point
      sol=0.;
      for(k=0;k<N_DOF;k++)
       sol+= LPS_sol[DOF[k]]*Value[k];

      for(k=0;k<N_CoarseDOF;k++)
       {
        val = w*CoarseValue[k];
        rhs[k] += sol*val;

        for(l=0;l<N_CoarseDOF;l++)
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];

       } // end for k
     } // for(j=0;j<N_Points;j++)

//     for(j=0;j<N_CoarseDOF;j++)
//       for(k=0;k<N_CoarseDOF;k++)
//         cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
//     cout  << endl;
//     for(j=0;j<N_CoarseDOF;j++)
//         cout << j << "  " << rhs[j] << endl;

    SolveLinearSystemLapack(G, rhs, N_CoarseDOF, N_CoarseDOF);

// interpolate the discont solution to the original LPS space

    auto nf = CurrentElement.GetNodalFunctional();
    nf->GetPointsForAll(N_Points, xi, eta);

    for(j=0;j<N_Points;j++)
     CoarseBF->GetDerivatives(MultiIndex2D::D00, xi[j], eta[j], PCValues[j]);

    memset(PointValues, 0, N_Points*sizeof(double));

    for(j=0;j<N_Points;j++)
     for(k=0;k<N_CoarseDOF;k++)
      PointValues[j] +=  rhs[k]*PCValues[j][k];

    nf->GetAllFunctionals(coll, cell, PointValues,
                          FunctionalValues);

    for(j=0;j<N_DOF;j++)
     {
      PostSol[DOF[j]] += FunctionalValues[j];
      W[DOF[j]] += 1.;
      if(j==N_DOF-1)
        {
//          cout<< j<< " PostSol " << PostSol[DOF[j]] << endl;
          if (maxbubble< PostSol[DOF[j]]) maxbubble = PostSol[DOF[j]];
          if (minbubble > PostSol[DOF[j]]) minbubble = PostSol[DOF[j]];
        }
     }
   } // for(i=0;i<N_Cells;i++)

  for(i=0;i<N_U;i++)
   PostSol[i] /=W[i];

//   cout<<" maxbubble " << maxbubble  <<" minbubble " << minbubble  << endl;
//     exit(0);
}


//Fefunction can be a different FEspace - Sashi
void AddALEStreamlineLPS(TSquareMatrix2D* A, int N_FeFunct, TFEFunction2D **FeFunct,
                         double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m;
  int N_Cells;
  int CoarseOrder;
  int N_CoarseDOF, N_DOF, N_UDOF, N_WDOF;
  std::shared_ptr<const TFESpace2D> W_fespace;
  TFEFunction2D *uh1, *uh2, *wh1, *wh2;  
  
  bool SecondDer[2] = { false, false };
  int N_Points;
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
//   double **ChildValues, *ChildValue;
  double **U_Values, *U_Value; 
  double **W_Values, *W_Value;
  
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  const int *RowPtr, *KCol;
  double *Entries;
  double *Values1, *Values2, *Values3, *Values4;
  double BValue[MaxN_BaseFunctions2D];

  auto fespace = A->GetFESpace2D();
  ActiveBound = fespace->get_n_active();
  RowPtr = A->get_row_ptr();
  KCol = A->get_vector_columns();
  Entries = A->GetEntries();
  // cout << "" << endl;

  if(N_FeFunct==2)
   {
    uh1 = FeFunct[0];
    uh2 = FeFunct[1];      
   }
  else if(N_FeFunct==4)
   {
    uh1 = FeFunct[0];
    uh2 = FeFunct[1];   
    wh1 = FeFunct[2];
    wh2 = FeFunct[3]; 
   }
  else
   {
    cout << "N_FeFunct must be 2 or 4 " <<endl;
    exit(0);    
   }
   
  auto U_fespace = uh1->GetFESpace2D();
  Values1 = uh1->GetValues();
  Values2 = uh2->GetValues();
    
  if(N_FeFunct==4) 
   {
    W_fespace = wh1->GetFESpace2D();  
    Values3 = wh1->GetValues();
    Values4 = wh2->GetValues();       
   }
  
  auto Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTria); // dummy type
  TQuadFormula qf_orig(qf_ref);

  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    auto CurrentElement = fespace->get_fe(i);
    auto BF = CurrentElement.GetBaseFunct();
    N_DOF = BF->GetDimension();

    auto U_CurrentElement = U_fespace->get_fe(i);
    auto U_BF = U_CurrentElement.GetBaseFunct();
    N_UDOF = U_BF->GetDimension();
    
    auto W_CurrentElement = W_fespace->get_fe(i);
    auto W_BF = W_CurrentElement.GetBaseFunct();
    N_WDOF = W_BF->GetDimension();
    
    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    const FiniteElement CoarseElement(GetElement2D(cell, CoarseOrder));
    auto CoarseBF = CoarseElement.GetBaseFunct();

    // quadrature formula on cell
    std::vector<const FiniteElement*> UsedElements{{&CoarseElement, &CurrentElement}};
    FEDatabase::GetOrig(UsedElements, Coll, cell, SecondDer, qf_ref, qf_orig);
    N_Points = qf_orig.GetN_QuadPoints();

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = FEDatabase::GetOrigElementValues(*CoarseBF,
                                                    MultiIndex2D::D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = qf_orig.get_weight(j);
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && std::abs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*sizeof(double));
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    auto DOF = fespace->GetGlobalDOF(i);
    auto U_DOF = U_fespace->GetGlobalDOF(i);
    auto W_DOF = W_fespace->GetGlobalDOF(i);

    PCValues = FEDatabase::GetOrigElementValues(*CoarseBF, MultiIndex2D::D00);

    ChildValuesX = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D10);
    ChildValuesY = FEDatabase::GetOrigElementValues(*BF, MultiIndex2D::D01);

    U_Values  = FEDatabase::GetOrigElementValues(*U_BF, MultiIndex2D::D00);
    W_Values  = FEDatabase::GetOrigElementValues(*W_BF, MultiIndex2D::D00);

    memset(H, 0, N_CoarseDOF*N_DOF*sizeof(double));
    memset(LocMatrix, 0, N_DOF*N_DOF*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      U_Value  = U_Values[j];
      W_Value  = W_Values[j];      
      
      w = qf_orig.get_weight(j);
      valx = 0.0;
      valy = 0.0;     
      // compute components of uh in j
      for(k=0;k<N_UDOF;k++)
      {
        l = U_DOF[k];
        valx += U_Value[k]*Values1[l];
        valy += U_Value[k]*Values2[l];
      }
      // sub mesh velo (uh-wh)
      for(k=0;k<N_WDOF;k++)
      {
        l = W_DOF[k];
        valx -= W_Value[k]*Values3[l];
        valy -= W_Value[k]*Values4[l];
      }

      for(k=0;k<N_DOF;k++)
      {
        BValue[k] = valx*ChildValueX[k] + valy*ChildValueY[k];
      }
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*N_DOF+l      ] += val*BValue[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*BValue[k]*BValue[l];
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*N_DOF*sizeof(double));

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, N_DOF, N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2] * H[i1*N_DOF+l] * H[i2*N_DOF+m];
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*N_DOF+l      ] * H[i2*N_DOF+m      ];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*N_DOF+m      ] * H[i1*N_DOF+l      ];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if(dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              Entries[p] += lpcoeff*std::pow(hK,lpexponent)*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // AddStreamlineTerm
#else


FE_type GetElement3D(TBaseCell *cell, int CoarseOrder)
{
  FE_type ele = (FE_type)0;
  Shapes shapetype;

  shapetype = cell->GetType();
  switch(shapetype)
  {
    // regularly refined Hexahedron
    case Hexahedron:
    case Brick:
      switch(CoarseOrder)
      {
        case 0:
          ele = C_Q0_3D_H_M;
        break;

        case 1:
          ele = D_P1_3D_H_M;
        break;

        case 2:
          ele = D_P2_3D_H_M;
        break;

        case 3:
          ele = D_P3_3D_H_M;
        break;

        default:
          if(CoarseOrder<0)
          {
            ele = C_Q00_3D_H_M;
          }
          else
          {
            ErrThrow("CoarseOrder: ", CoarseOrder,
                     " Projection space is defined up to order 3");
          }
      } // end switch CoarseOrder
    break; // end regularly refined quadrilateral

    case Tetrahedron:
      switch(CoarseOrder)
      {
        case 0:
          ele = C_P0_3D_T_A;
        break;

        case 1:
          ele = D_P1_3D_T_A;
        break;

        default:
          if(CoarseOrder<0)
          {
            ele = C_P00_3D_T_A;
          }
          else
          {
            ErrThrow("CoarseOrder: ", CoarseOrder,
                     " Projection space is defined up to order 1");
          }
      } // end switch CoarseOrder
    break;
    default:
      ErrThrow("Invalid shape");
  } // end switch reftype
  return ele;
}

// ADDED ON 17.06.2011 BY SASHI
void AddStreamlineTerm(TSquareMatrix3D* A, TFEFunction3D *uh1,
                       TFEFunction3D *uh2, TFEFunction3D *uh3,
                       double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m;
  int N_Cells;
  int CoarseOrder;
  int N_CoarseDOF, N_DOF;
  int N_Points;
  int i1, i2;
  int ActiveBound, dof;
  int p, end;
  const int *RowPtr, *KCol;

  bool SecondDer[2] = { false, false };

  double G[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double Gsave[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double H[2*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double P[2*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValuesZ, *ChildValueZ;  
  double **ChildValues, *ChildValue;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy, valz;
  double LocMatrix[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double s, hK;
  double *Entries;
  double *Values1, *Values2, *Values3;
  double BValue[MaxN_BaseFunctions3D];

  auto fespace = A->GetFESpace3D();
  ActiveBound = fespace->get_n_active();
  RowPtr = A->get_row_ptr();
  KCol = A->get_vector_columns();
  Entries = A->GetEntries();
  // cout << "" << endl;

  auto Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);
  
  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTetra); // dummy type
  TQuadFormula qf_orig(qf_ref);

  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    auto CurrentElement = fespace->get_fe(i);

    auto BF = CurrentElement.GetBaseFunct();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    const FiniteElement CoarseElement(GetElement3D(cell, CoarseOrder));
    auto CoarseBF = CoarseElement.GetBaseFunct();

    // quadrature formula on cell
    std::vector<const FiniteElement*> UsedElements{{&CoarseElement, &CurrentElement}};
    FEDatabase::GetOrig(UsedElements, Coll, cell, SecondDer, qf_ref, qf_orig);
    N_Points = qf_orig.GetN_QuadPoints();

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = FEDatabase::GetOrigElementValues(*CoarseBF,
                                                    MultiIndex3D::D000);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = qf_orig.get_weight(j);
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && std::abs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*sizeof(double));
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    auto DOF = fespace->GetGlobalDOF(i);

    Values1 = uh1->GetValues();
    Values2 = uh2->GetValues();
    Values3 = uh3->GetValues();
    
    PCValues = FEDatabase::GetOrigElementValues(*CoarseBF,
                                                MultiIndex3D::D000);

    ChildValuesX = FEDatabase::GetOrigElementValues(*BF, MultiIndex3D::D100);
    ChildValuesY = FEDatabase::GetOrigElementValues(*BF, MultiIndex3D::D010);
    ChildValuesZ = FEDatabase::GetOrigElementValues(*BF, MultiIndex3D::D001);
    ChildValues  = FEDatabase::GetOrigElementValues(*BF, MultiIndex3D::D000);

    memset(H, 0, N_CoarseDOF*N_DOF*sizeof(double));

    memset(LocMatrix, 0, N_DOF*N_DOF*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValueZ = ChildValuesZ[j];     
      ChildValue  = ChildValues[j];
      w = qf_orig.get_weight(j);
      valx = 0.0;
      valy = 0.0;
      valz = 0.0;
      
      // compute components of uh in j
      for(k=0;k<N_DOF;k++)
      {
        l = DOF[k];
        valx += ChildValue[k]*Values1[l];
        valy += ChildValue[k]*Values2[l];
        valz += ChildValue[k]*Values3[l];        
      }
      for(k=0;k<N_DOF;k++)
      {
        BValue[k] = valx*ChildValueX[k] + valy*ChildValueY[k] + valz*ChildValueZ[k];
      }
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*N_DOF+l] += val*BValue[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*BValue[k]*BValue[l];
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*N_DOF*sizeof(double));

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, N_DOF, N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2] * H[i1*N_DOF+l] * H[i2*N_DOF+m];
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*N_DOF+l      ] * H[i2*N_DOF+m      ];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*N_DOF+m      ] * H[i1*N_DOF+l      ];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if(dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              Entries[p] += lpcoeff*std::pow(hK,lpexponent)*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // AddStreamlineTerm


// stabilisation of full gradient (scalar)
void UltraLocalProjection(TSquareMatrix3D* A, 
                          double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m;
  int N_Cells;
  int CoarseOrder;
  int N_CoarseDOF, N_DOF;
  int N_Points;
  int i1, i2;
  int ActiveBound, dof;
  int p, end;
  const int *RowPtr, *KCol;

//  Shapes shapetype;
  bool SecondDer[2] = { false, false };

  double G[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double Gsave[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double H[3*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double P[3*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double **CoarseValues, *CoarseValue; 
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValuesZ, *ChildValueZ;  
  double **PCValues;
  double *PCValue;
  double w, val;  // double valx, valy, valz;
  double LocMatrix[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double s, hK;
  double *Entries;
//  double BValue[MaxN_BaseFunctions3D];

  auto fespace = A->GetFESpace3D();
  ActiveBound = fespace->get_n_active();
  RowPtr = A->get_row_ptr();
  KCol = A->get_vector_columns();
  Entries = A->GetEntries();
  // cout << "" << endl;

  auto Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTetra); // dummy type
  TQuadFormula qf_orig(qf_ref);

  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    auto CurrentElement = fespace->get_fe(i);

    auto BF = CurrentElement.GetBaseFunct();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    const FiniteElement CoarseElement(GetElement3D(cell, CoarseOrder));
    auto CoarseBF = CoarseElement.GetBaseFunct();

    // quadrature formula on cell
    std::vector<const FiniteElement*> UsedElements{{&CoarseElement, &CurrentElement}};
    FEDatabase::GetOrig(UsedElements, Coll, cell, SecondDer, qf_ref, qf_orig);
    N_Points = qf_orig.GetN_QuadPoints();
    
    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = FEDatabase::GetOrigElementValues(*CoarseBF,
                                                    MultiIndex3D::D000);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = qf_orig.get_weight(j);
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && std::abs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*sizeof(double));
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    auto DOF = fespace->GetGlobalDOF(i);
    
    PCValues = FEDatabase::GetOrigElementValues(*CoarseBF,
                                                MultiIndex3D::D000);

    ChildValuesX = FEDatabase::GetOrigElementValues(*BF, MultiIndex3D::D100);
    ChildValuesY = FEDatabase::GetOrigElementValues(*BF, MultiIndex3D::D010);
    ChildValuesZ = FEDatabase::GetOrigElementValues(*BF, MultiIndex3D::D001);
 
    memset(H, 0, N_CoarseDOF*3*N_DOF*sizeof(double));

    memset(LocMatrix, 0, N_DOF*N_DOF*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValueZ = ChildValuesZ[j];     
      w = qf_orig.get_weight(j);

      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*3*N_DOF+l      ] += val*ChildValueX[l];
          H[k*3*N_DOF+l+N_DOF] += val*ChildValueY[l];
          H[k*3*N_DOF+l+2*N_DOF] += val*ChildValueZ[l];  
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*(  ChildValueX[k]*ChildValueX[l]
                                     + ChildValueY[k]*ChildValueY[l]
                                     + ChildValueZ[k]*ChildValueZ[l]);
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*3*N_DOF*sizeof(double));

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 3*N_DOF, 3*N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*N_DOF+k] << endl;
    */
    
    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2]*( H[i1*3*N_DOF+l      ]*H[i2*3*N_DOF+m      ]
                                           +H[i1*3*N_DOF+l+N_DOF]*H[i2*3*N_DOF+m+N_DOF]
                                           +H[i1*3*N_DOF+l+2*N_DOF]*H[i2*3*N_DOF+m+2*N_DOF]);
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*3*N_DOF+l      ] * H[i2*3*N_DOF+m      ];
          s += P[i2*3*N_DOF+l+N_DOF] * H[i2*3*N_DOF+m+N_DOF];
          s += P[i2*3*N_DOF+l+2*N_DOF] * H[i2*3*N_DOF+m+2*N_DOF];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*3*N_DOF+m      ] * H[i1*3*N_DOF+l      ];
          s += P[i1*3*N_DOF+m+N_DOF] * H[i1*3*N_DOF+l+N_DOF];
          s += P[i1*3*N_DOF+m+2*N_DOF] * H[i1*3*N_DOF+l+2*N_DOF];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if(dof<ActiveBound)
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              Entries[p] += lpcoeff*std::pow(hK,lpexponent)*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // AddStreamlineTerm


//**************************************************************
//  UltraLocalProjection
//  checked: Volker John 08/02/19
//**************************************************************
#ifdef __2D__
void UltraLocalProjection(void* A, bool ForPressure, CoeffFct2D Coeff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *coll;
  TFESpace2D *fespace;
  FE2D UsedElements[2];
  int N_UsedElements;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  const double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s, sx, sy;
  int i1, i2, ij, N_Edges;
  double hK, *coeffs, *params;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *Entries;
  int OrderDiff;
  double lpcoeff, lpexponent;

  Output::print("LPS full gradient");

  coeffs = new double[20];
  params = new double[10];
  memset(params, 0, 10 * sizeof(double));

  if(!(TDatabase::ParamDB->LP_FULL_GRADIENT) && !(ForPressure))
  {
    Output::print("Local projection stabilization is implemented only for full gradient!");
  }

  if(ForPressure)
  {
    lpcoeff = TDatabase::ParamDB->LP_PRESSURE_COEFF;
    lpexponent = TDatabase::ParamDB->LP_PRESSURE_EXPONENT;
    OrderDiff = TDatabase::ParamDB->LP_PRESSURE_ORDER_DIFFERENCE;
  }
  else
  {
    lpcoeff = TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF;
    lpexponent = TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT;
    OrderDiff = TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE;
  }

  // get fespace and matrices
  if(ForPressure)
  {
    fespace = (TFESpace2D*)(((TMatrix2D *)A)->GetStructure().GetTestSpace());
    ActiveBound = -1;
    RowPtr = ((TMatrix2D *)A)->get_row_ptr();
    KCol = ((TMatrix2D *)A)->get_vector_columns();
    Entries = ((TMatrix2D *)A)->GetEntries();
    // cout << "for pressure" << endl;
  }
  else
  {
    fespace = ((TSquareMatrix2D *)A)->GetFESpace();
    ActiveBound = fespace->get_n_active();
    RowPtr = ((TSquareMatrix2D *)A)->get_row_ptr();
    KCol = ((TSquareMatrix2D *)A)->get_vector_columns();
    Entries = ((TSquareMatrix2D *)A)->GetEntries();
    // cout << "not for pressure" << endl;
  }
  // get collection
  coll = fespace->GetCollection();
  N_Cells = coll->GetN_Cells();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    coll->GetCell(i)->SetClipBoard(i);
  // loop over the cells
  for(i=0;i<N_Cells;i++)
  {
      // get cell
    cell = coll->GetCell(i);
    // get diameter of the cell
    switch (TDatabase::ParamDB->CELL_MEASURE)
    {
  case 0: // diameter
      hK = cell->GetDiameter();
      break;
  case 1: // with reference map
      hK = cell->GetLengthWithReferenceMap();
      break;
  case 2: // shortest edge
      hK = cell->GetShortestEdge();
      break;
  case 3: // measure
      hK = cell->GetMeasure();
      hK = std::pow(hK,1.0/3.0);
      break;
  case 4: // mesh cell in convection direction
      N_Edges = cell->GetN_Edges();
      sx = sy = 0;
      for (ij=0;ij<N_Edges;ij++)
      {
    TDatabase::ParamDB->INTERNAL_VERTEX_X[ij] = cell->GetVertex(ij)->GetX();
    sx += TDatabase::ParamDB->INTERNAL_VERTEX_X[ij];
    TDatabase::ParamDB->INTERNAL_VERTEX_Y[ij] = cell->GetVertex(ij)->GetY();
    sy += TDatabase::ParamDB->INTERNAL_VERTEX_Y[ij];
      }
      if (N_Edges==3)
    TDatabase::ParamDB->INTERNAL_VERTEX_X[3] = -4711;
      // center of mesh cell
      sx /= N_Edges;
      sy /= N_Edges;
      hK = cell->GetDiameter(); 
      // get coefficients in center of mesh cell 
      Coeff(1, &sx ,&sy, &params, &coeffs);
      hK = Mesh_size_in_convection_direction(hK, coeffs[1], coeffs[2]); 
      break;
  default: // diameter
      hK = cell->GetDiameter();
      break;
    }
    // get finite element in the cell
    auto CurrentElement = fespace->get_fe(i);
    // get basis functions of the finite element
    auto BF = CurrentElement.GetBaseFunct();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();
    // compute order of the local coarse space
    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // get type of the mesh cell
    shapetype = cell->GetType();

    // determine finite element of the local coarse space
    switch(shapetype)
    {
      // regularly refined quadrilateral
      // discontinuous finite elements  
      case Quadrangle:
      case Parallelogram:
      case Rectangle:
        switch(CoarseOrder)
        {
          case 0:
            UsedElements[0] = C_Q0_2D_Q_M;
          break;

          case 1:
            UsedElements[0] = D_P1_2D_Q_M;
          break;

          case 2:
            UsedElements[0] = D_P2_2D_Q_M;
          break;

          case 3:
            UsedElements[0] = D_P3_2D_Q_M;
          break;

          case 4:
            UsedElements[0] = D_P4_2D_Q_M;
          break;

          default:
            ErrThrow("Projection space is defined upto order 4");
            break;
        } // end switch CoarseOrder
      break; // end regularly refined quadrilateral

      case Triangle:
        switch(CoarseOrder)
        {
          case 0:
            UsedElements[0] = C_P0_2D_T_A;
          break;

          case 1:
            UsedElements[0] = D_P1_2D_T_A;
          break;

          case 2:
            UsedElements[0] = D_P2_2D_T_A;
          break;

          case 3:
            UsedElements[0] = D_P3_2D_T_A;
          break;

          case 4:
            UsedElements[0] = D_P4_2D_T_A;
          break;

          default:
            ErrThrow("Projection space is defined upto order 4");
            break;
        }
      break;
      default: 
        ErrThrow("No coarse finite elements for mesh cell type ", shapetype,
                 " implemented !!!");
        break;
    } // end switch reftype

    // approximation space (index 1) and projection space (index 0)
    N_UsedElements = 2;
    UsedElements[1] = CurrentElement.GetID();

    const FiniteElement CoarseElement(UsedElements[0]);
    auto CoarseBF = CoarseElement.GetBaseFunct();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    FEDatabase::GetOrig(N_UsedElements, UsedElements, coll, cell, SecondDer,
                        N_Points, xi, eta, weights, X, Y, AbsDetjk);
    // number of dof in the local coarse space
    N_CoarseDOF = CoarseBF->GetDimension();
    // get function values for the basis functions of the coarse space
    CoarseValues = FEDatabase::GetOrigElementValues(CoarseBF_ID, D00);
    // initialize array G, stores mass matrix of coarse space
    memset(G, 0, N_CoarseDOF*N_CoarseDOF*sizeof(double));
    // loop over the quadrature points
    for(j=0;j<N_Points;j++)
    {
      // values of the coarse basis functions in the quad points
      CoarseValue = CoarseValues[j];
      // factor for numerical quadrature
      w = AbsDetjk[j]*weights[j];
      // loop over the basis functions of the coarse space
      for(k=0;k<N_CoarseDOF;k++)
      {
    // first factor of integrand
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
      // update integral 
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && std::abs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*sizeof(double));
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    // this should be the same as CoarseValues ???
    PCValues = FEDatabase::GetOrigElementValues(CoarseBF_ID, D00);
    // get derivatives of fe functions of approximation space
    ChildValuesX = FEDatabase::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = FEDatabase::GetOrigElementValues(BF_ID, D01);
    // initialize array H, holds mixed products of coarse basis 
    // functions and derivatives of fine basis functions
    memset(H, 0, N_CoarseDOF*2*N_DOF*sizeof(double));
    // initialize array LocMatrix, holds products of derivatives
    // of fine basis functions
    memset(LocMatrix, 0, N_DOF*N_DOF*sizeof(double));
    // loop over the quadrature points
    for(j=0;j<N_Points;j++)
    {
  // get all values in the quad point j 
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      // factor in the integral
      w = AbsDetjk[j]*weights[j];
      // loop over the basis functions of the coarse space
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
  // compute products of coarse basis function and
  // derivatives of fine basis functions
  // values for the same quad point are stored after
  // each other (x deriv of all fine fcts., y deriv
  // of all fine fcts.) tested with the k-th coarse fct. 
        for(l=0;l<N_DOF;l++)
        {
          H[k*2*N_DOF+l      ] += val*ChildValueX[l];
          H[k*2*N_DOF+l+N_DOF] += val*ChildValueY[l];
        } // end for l
      } // end for k

      // grad-grad matrix (fine-fine coupling)
      // LocMatrix will store the local update of the global
      // system matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*( ChildValueX[k]*ChildValueX[l]
                                     +ChildValueY[k]*ChildValueY[l]);
        }
      }
    } // end for j
    // copy array H to array P
    memcpy(P, H, N_CoarseDOF*2*N_DOF*sizeof(double));

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // solve G * X = H, solution X stored on H
    // the right hand side and the solution are stored column wise
    // right hand side: a fine fct. tested with all coarse fcts. 
    // X - coefficients of local L^2 projection of grad of fine function
    //     ((first fct.)_x, first coarse fct.), ((second fct.)_x, first coarse fct.), 
    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 2*N_DOF, 2*N_DOF);
    //SolveMultipleSystems(G, H, N_CoarseDOF, N_CoarseDOF, 2*N_DOF, 2*N_DOF);

    /*// checked 08/02/19
    double val;
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
    {
  for(l=0;l<2*N_DOF;l++)
  {
      val = 0;
      
      for(k=0;k<N_CoarseDOF;k++)
      {
    val += G[j*N_CoarseDOF+k] * H[l+k*2*N_DOF];
    Output::print(l+2*k*N_DOF);
      }
      Output::print(j, " ", l, " ", P[j*2*N_DOF+l], " ", val);
  }
    }
    exit(1);
    */

    // update LocMatrix
    // proj-proj coupling (coarse-coarse coupling)
    // l - test function index
    for(l=0;l<N_DOF;l++)
    {
  // m - ansatz function index
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2]*( H[i1*2*N_DOF+l      ]*H[i2*2*N_DOF+m      ]
                                           +H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF]);
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling (fine-coarse couplings)
    // l - test function index
    for(l=0;l<N_DOF;l++)
    {
  // m - ansatz function index
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m      ];
          s += P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l      ];
          s += P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    // get numbers of global dof belonging to this mesh cell
    DOF = fespace->GetGlobalDOF(i);

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if((dof<ActiveBound) || (ForPressure))
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
    // parameter is lpcoeff*std::pow(hK,lpexponent)
              Entries[p] += lpcoeff*std::pow(hK,lpexponent)*LocMatrix[l*N_DOF+m];
        break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
  delete params;
  delete coeffs;
} // UltraLocalProjection



bool TestCell(TBaseCell *cell)
{
  int i, N_;
  double x, y;
  bool val = true;
  
  N_ = cell->GetN_Vertices();
  
  for(i=0;i<N_;i++)
  {
    cell->GetVertex(i)->GetCoords(x,y);
    if(x<TDatabase::ParamDB->P6 && y<TDatabase::ParamDB->P6)
      val = val && true;
    else 
      val = val && false;    
  }
  return val;
}

// stabilisation of function velocity
void UltraLocalProjectionFunction(void* A, bool ForPressure)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D UsedElements[2];
  int N_UsedElements;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  const double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  int *RowPtr, *KCol;
  double *Entries;
  int OrderDiff;
  double lpcoeff, lpexponent;

  if(!(TDatabase::ParamDB->LP_FULL_GRADIENT) && !(ForPressure))
  {
    ErrThrow("Local projection stabilization is implemented only for full gradient!");
  }

  if(ForPressure)
  {
    lpcoeff = -(TDatabase::ParamDB->LP_PRESSURE_COEFF);
    lpexponent = TDatabase::ParamDB->LP_PRESSURE_EXPONENT;
    OrderDiff = TDatabase::ParamDB->LP_PRESSURE_ORDER_DIFFERENCE;
  }
  else
  {
    lpcoeff = TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF;
    lpexponent = TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT;
    OrderDiff = TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE;
  }

  if(ForPressure)
  {
    fespace = (TFESpace2D*)(((TMatrix2D *)A)->GetStructure().GetTestSpace());
    ActiveBound = -1;
    RowPtr = ((TMatrix2D *)A)->get_row_ptr();
    KCol = ((TMatrix2D *)A)->get_vector_columns();
    Entries = ((TMatrix2D *)A)->GetEntries();
    // cout << "for pressure" << endl;
  }
  else
  {
    fespace = ((TSquareMatrix2D *)A)->GetFESpace();
    ActiveBound = fespace->get_n_active();
    RowPtr = ((TSquareMatrix2D *)A)->get_row_ptr();
    KCol = ((TSquareMatrix2D *)A)->get_vector_columns();
    Entries = ((TSquareMatrix2D *)A)->GetEntries();
    // cout << "not for pressure" << endl;
  }

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    auto CurrentElement = fespace->get_fe(i);

    auto BF = CurrentElement.GetBaseFunct();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approximation space (index 1) and projection space (index 0)
    N_UsedElements = 2;
    UsedElements[1] = CurrentElement.GetID();

    const FiniteElement CoarseElement(UsedElements[0]);
    auto CoarseBF = CoarseElement.GetBaseFunct();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    FEDatabase::GetOrig(N_UsedElements, UsedElements, Coll, cell, SecondDer,
                        N_Points, xi, eta, weights, X, Y, AbsDetjk);

    // cout << "N_Points: " << N_Points << endl;

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = FEDatabase::GetOrigElementValues(CoarseBF_ID, D00);

    // cout << "N_CoarseDOF: " << N_CoarseDOF << endl;

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && std::abs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*sizeof(double));
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    PCValues = FEDatabase::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = FEDatabase::GetOrigElementValues(BF_ID, D00);
    ChildValuesY = FEDatabase::GetOrigElementValues(BF_ID, D00);

    memset(H, 0, N_CoarseDOF*2*N_DOF*sizeof(double));

    memset(LocMatrix, 0, N_DOF*N_DOF*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*2*N_DOF+l      ] += val*ChildValueX[l];
          H[k*2*N_DOF+l+N_DOF] += val*ChildValueY[l];
        } // end for l
      } // end for k

      // (u,v)--matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*( ChildValueX[k]*ChildValueX[l]
                                     +ChildValueY[k]*ChildValueY[l]);
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*2*N_DOF*sizeof(double));

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 2*N_DOF, 2*N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // (\pi u, \pi v) proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2]*( H[i1*2*N_DOF+l      ]*H[i2*2*N_DOF+m      ]
                                           +H[i1*2*N_DOF+l+N_DOF]*H[i2*2*N_DOF+m+N_DOF]);
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // function-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*2*N_DOF+l      ] * H[i2*2*N_DOF+m      ];
          s += P[i2*2*N_DOF+l+N_DOF] * H[i2*2*N_DOF+m+N_DOF];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*2*N_DOF+m      ] * H[i1*2*N_DOF+l      ];
          s += P[i1*2*N_DOF+m+N_DOF] * H[i1*2*N_DOF+l+N_DOF];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    DOF = fespace->GetGlobalDOF(i);

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if((dof<ActiveBound) || (ForPressure))
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              Entries[p] += 1000.*hK*hK*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
} // UltraLocalProjection

double UltraLocalErrorSmooth(TFEFunction2D *uh, DoubleFunct2D *ExactU,
        double lpcoeff, double lpexponent, int OrderDiff)
{
  int i,j,k,l,m,n;
  int N_Cells;
  int *DOF;
  int CellOrder, CoarseOrder;
  int N_CoarseDOF, N_DOF;
  TCollection *Coll;
  TFESpace2D *fespace;
  FE2D CurrEleID, UsedElements[2];
  int N_UsedElements;
  FiniteElement *CurrentElement, *CoarseElement;
  TBaseFunct2D *BF, *CoarseBF;
  BaseFunct2D BF_ID, CoarseBF_ID;
  TBaseCell *cell;
  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  const double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double Gsave[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double H[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double P[2*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy;
  double LocMatrix[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  double *Values;
  double error, locerror;
  double exactval[4];

  error = 0.0;

  fespace = uh->GetFESpace2D();
  Values = uh->GetValues();

  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    locerror = 0.0;
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();
    
    if((TestCell(cell)) == false) continue;
    
    DOF = fespace->GetGlobalDOF(i);

    auto CurrentElement = fespace->get_fe(i);
    CurrEleID = CurrentElement.GetID();

    BF = CurrentElement.GetBaseFunct();
    BF_ID = BF->GetID();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    UsedElements[0] = GetElement2D(cell, CoarseOrder);

    // approximation space (index 1) and projection space (index 0)
    N_UsedElements = 2;
    UsedElements[1] = CurrEleID;

    const FiniteElement CoarseElement(UsedElements[0]);
    CoarseBF = CoarseElement.GetBaseFunct();
    CoarseBF_ID = CoarseBF->GetID();

    // quadrature formula on cell
    FEDatabase::GetOrig(N_UsedElements, UsedElements, Coll, cell, SecondDer,
                        N_Points, xi, eta, weights, X, Y, AbsDetjk);

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = FEDatabase::GetOrigElementValues(CoarseBF_ID, D00);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && std::abs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*sizeof(double));
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    PCValues = FEDatabase::GetOrigElementValues(CoarseBF_ID, D00);

    ChildValuesX = FEDatabase::GetOrigElementValues(BF_ID, D10);
    ChildValuesY = FEDatabase::GetOrigElementValues(BF_ID, D01);

    // only two right-hand sides (x and y derivative)
    memset(H, 0, N_CoarseDOF*2*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];

      // calculate gradient of discrete uh in this quadrature point
      valx = 0.0;
      valy = 0.0;
      for(k=0;k<N_DOF;k++)
      {
        val = Values[DOF[k]];
        valx += ChildValueX[k]*val;
        valy += ChildValueY[k]*val;
      }

      // get gradient of exact u
      ExactU(X[j], Y[j], exactval);

      valx -= exactval[1];
      valy -= exactval[2];

      w = AbsDetjk[j]*weights[j];
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        H[k*2  ] += val*valx;
        H[k*2+1] += val*valy;
      } // end for k

      // grad-grad term
      locerror += w*(valx*valx + valy*valy);

    } // end for j
    memcpy(P, H, N_CoarseDOF*2*sizeof(double));

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 2, 2);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    s = 0;
    for(i1=0;i1<N_CoarseDOF;i1++)
      for(i2=0;i2<N_CoarseDOF;i2++)
        s += Gsave[i1*N_CoarseDOF+i2]*( H[i1*2  ]*H[i2*2  ]
                                       +H[i1*2+1]*H[i2*2+1]);
    locerror += s;

    // grad-proj coupling
    s = 0;
    for(i2=0;i2<N_CoarseDOF;i2++)
    {
      s += P[i2*2  ] * H[i2*2  ];
      s += P[i2*2+1] * H[i2*2+1];
    }
    for(i1=0;i1<N_CoarseDOF;i1++)
    {
      s += P[i1*2  ] * H[i1*2  ];
      s += P[i1*2+1] * H[i1*2+1];
    }
    locerror -= s;

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    error += lpcoeff*std::pow(hK,lpexponent)*locerror;
  } // endfor i

  return error;
} // UltraLocalProjection

#endif // __2D__

// stabilisation of full gradient (velocity or pressure)
void UltraLocalProjection3D(void* A, bool ForPressure)
{
  int i,j,k,l,m;
  int N_Cells;
  int CoarseOrder;
  int N_CoarseDOF, N_DOF;
  std::shared_ptr<const TFESpace3D> fespace;
//  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double G[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double Gsave[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double H[3*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double P[3*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValuesZ, *ChildValueZ;
  double **PCValues;
  double *PCValue;
  double w, val;
  double LocMatrix[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double s;
  int i1, i2;
  double hK;
  int ActiveBound, dof;
  int p, end;
  const int *RowPtr, *KCol;
  double *Entries;
  int OrderDiff;
  double lpcoeff, lpexponent;

  if(!(TDatabase::ParamDB->LP_FULL_GRADIENT) && !(ForPressure))
  {
    ErrThrow("Local projection stabilization is implemented only for full gradient!");
  }

  if(ForPressure)
  {
    lpcoeff = -(TDatabase::ParamDB->LP_PRESSURE_COEFF);
    lpexponent = TDatabase::ParamDB->LP_PRESSURE_EXPONENT;
    OrderDiff = TDatabase::ParamDB->LP_PRESSURE_ORDER_DIFFERENCE;
  }
  else
  {
    lpcoeff = TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF;
    lpexponent = TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT;
    OrderDiff = TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE;
  }

  if(ForPressure)
  {
    fespace = ((TMatrix3D *)A)->GetTestSpace3D();
    ActiveBound = -1;
    RowPtr = ((TMatrix3D *)A)->get_row_ptr();
    KCol = ((TMatrix3D *)A)->get_vector_columns();
    Entries = ((TMatrix3D *)A)->GetEntries();
    // cout << "for pressure" << endl;
  }
  else
  {
    fespace = ((TSquareMatrix3D *)A)->GetFESpace3D();
    ActiveBound = fespace->get_n_active();
    RowPtr = ((TSquareMatrix3D *)A)->get_row_ptr();
    KCol = ((TSquareMatrix3D *)A)->get_vector_columns();
    Entries = ((TSquareMatrix3D *)A)->GetEntries();
    // cout << "not for pressure" << endl;
  }

  auto Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTetra); // dummy type
  TQuadFormula qf_orig(qf_ref);

  for(i=0;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    auto CurrentElement = fespace->get_fe(i);

    auto BF = CurrentElement.GetBaseFunct();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    const FiniteElement CoarseElement(GetElement3D(cell, CoarseOrder));
    auto CoarseBF = CoarseElement.GetBaseFunct();

    // quadrature formula on cell
    std::vector<const FiniteElement*> UsedElements{{&CoarseElement, &CurrentElement}};
    FEDatabase::GetOrig(UsedElements, Coll, cell, SecondDer, qf_ref, qf_orig);
    N_Points = qf_orig.GetN_QuadPoints();
    
    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = FEDatabase::GetOrigElementValues(*CoarseBF,
                                                    MultiIndex3D::D000);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = qf_orig.get_weight(j);
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && std::abs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*sizeof(double));

    PCValues = FEDatabase::GetOrigElementValues(*CoarseBF,
                                                MultiIndex3D::D000);

    ChildValuesX = FEDatabase::GetOrigElementValues(*BF, MultiIndex3D::D100);
    ChildValuesY = FEDatabase::GetOrigElementValues(*BF, MultiIndex3D::D010);
    ChildValuesZ = FEDatabase::GetOrigElementValues(*BF, MultiIndex3D::D001);

    memset(H, 0, N_CoarseDOF*3*N_DOF*sizeof(double));

    memset(LocMatrix, 0, N_DOF*N_DOF*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValueZ = ChildValuesZ[j];
      w = qf_orig.get_weight(j);
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        for(l=0;l<N_DOF;l++)
        {
          H[k*3*N_DOF+l        ] += val*ChildValueX[l];
          H[k*3*N_DOF+l+N_DOF  ] += val*ChildValueY[l];
          H[k*3*N_DOF+l+2*N_DOF] += val*ChildValueZ[l];
        } // end for l
      } // end for k

      // grad-grad matrix
      for(k=0;k<N_DOF;k++)
      {
        for(l=0;l<N_DOF;l++)
        {
          LocMatrix[k*N_DOF+l] += w*( ChildValueX[k]*ChildValueX[l]
                                     +ChildValueY[k]*ChildValueY[l]
                                     +ChildValueZ[k]*ChildValueZ[l]);
        }
      }
    } // end for j
    memcpy(P, H, N_CoarseDOF*3*N_DOF*sizeof(double));

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 3*N_DOF, 3*N_DOF);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
        for(i1=0;i1<N_CoarseDOF;i1++)
          for(i2=0;i2<N_CoarseDOF;i2++)
            s += Gsave[i1*N_CoarseDOF+i2]*( H[i1*3*N_DOF+l        ]*H[i2*3*N_DOF+m        ]
                                           +H[i1*3*N_DOF+l+N_DOF  ]*H[i2*3*N_DOF+m+N_DOF  ]
                                           +H[i1*3*N_DOF+l+2*N_DOF]*H[i2*3*N_DOF+m+2*N_DOF]);
        LocMatrix[l*N_DOF+m] += s;
      } // endfor m
    } // endfor l

    // grad-proj coupling
    for(l=0;l<N_DOF;l++)
    {
      for(m=0;m<N_DOF;m++)
      {
        s = 0;
//        for(i3=0;i3<N_CoarseDOF;i3++)
//        {
//          s += P[i3*3*N_DOF+l        ] * H[i3*3*N_DOF+m        ];
//          s += P[i3*3*N_DOF+l+N_DOF  ] * H[i3*3*N_DOF+m+N_DOF  ];
//          s += P[i3*3*N_DOF+l+2*N_DOF] * H[i3*3*N_DOF+m+2*N_DOF];
//        }
        for(i2=0;i2<N_CoarseDOF;i2++)
        {
          s += P[i2*3*N_DOF+l        ] * H[i2*3*N_DOF+m        ];
          s += P[i2*3*N_DOF+l+N_DOF  ] * H[i2*3*N_DOF+m+N_DOF  ];
          s += P[i2*3*N_DOF+l+2*N_DOF] * H[i2*3*N_DOF+m+2*N_DOF];
        }
        for(i1=0;i1<N_CoarseDOF;i1++)
        {
          s += P[i1*3*N_DOF+m        ] * H[i1*3*N_DOF+l        ];
          s += P[i1*3*N_DOF+m+N_DOF  ] * H[i1*3*N_DOF+l+N_DOF  ];
          s += P[i1*3*N_DOF+m+2*N_DOF] * H[i1*3*N_DOF+l+2*N_DOF];
        }
        LocMatrix[l*N_DOF+m] -= s;
      } // end for m
    } // end for l

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    auto DOF = fespace->GetGlobalDOF(i);

    // add to global matrix
    for(l=0;l<N_DOF;l++)
    {
      dof = DOF[l];
      if((dof<ActiveBound) || (ForPressure))
      {
        for(m=0;m<N_DOF;m++)
        {
          end = RowPtr[dof+1];
          for(p=RowPtr[dof];p<end;p++)
            if(KCol[p] == DOF[m])
            {
              Entries[p] += lpcoeff*std::pow(hK,lpexponent)*LocMatrix[l*N_DOF+m];
              break;
            }
        } // endfor m
      } // endif dof<ActiveBound
    } // endfor l
  } // endfor i
  cout << "\n end of lps \n" << endl;
} // UltraLocalProjection

double UltraLocalError3D(TFEFunction3D *uh, DoubleFunct3D *ExactU,
        double lpcoeff, double lpexponent, int OrderDiff)
{
  // cout << "\n start of lps error \n" << endl;
  
  int i,j,k,l;
  int N_Cells;
  int CoarseOrder;
  int N_CoarseDOF, N_DOF;
//  Shapes shapetype;
  bool SecondDer[2] = { false, false };
  int N_Points;
  double G[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double Gsave[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double H[3*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double P[3*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double **CoarseValues, *CoarseValue;
  double **ChildValuesX, *ChildValueX;
  double **ChildValuesY, *ChildValueY;
  double **ChildValuesZ, *ChildValueZ;
  double **PCValues;
  double *PCValue;
  double w, val, valx, valy, valz;
//  double LocMatrix[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  double s;
  int i1, i2;
  double hK;
//  int ActiveBound, dof;
//  int p, end;
  double *Values;
  double error, locerror;
  double exactval[5];

  error = 0.0;

  auto fespace = uh->GetFESpace3D();
  Values = uh->GetValues();

  auto Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // initialise ClipBoard
  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  TQuadFormula qf_ref(QuadratureFormula_type::BaryCenterTetra); // dummy type
  TQuadFormula qf_orig(qf_ref);

  for(i=0;i<N_Cells;i++)
  {
    locerror = 0.0;
    auto cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    auto DOF = fespace->GetGlobalDOF(i);

    auto CurrentElement = fespace->get_fe(i);

    auto BF = CurrentElement.GetBaseFunct();
    N_DOF = BF->GetDimension();

    CoarseOrder = BF->GetAccuracy() - OrderDiff;
    // cout << "CoarseOrder: " << CoarseOrder << endl;

    const FiniteElement CoarseElement(GetElement3D(cell, CoarseOrder));
    auto CoarseBF = CoarseElement.GetBaseFunct();

    // quadrature formula on cell
    std::vector<const FiniteElement*> UsedElements{{&CoarseElement, &CurrentElement}};
    FEDatabase::GetOrig(UsedElements, Coll, cell, SecondDer, qf_ref, qf_orig);
    N_Points = qf_orig.GetN_QuadPoints();

    N_CoarseDOF = CoarseBF->GetDimension();
    CoarseValues = FEDatabase::GetOrigElementValues(*CoarseBF,
                                                    MultiIndex3D::D000);

    memset(G, 0, N_CoarseDOF*N_CoarseDOF*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      CoarseValue = CoarseValues[j];
      w = qf_orig.get_weight(j);
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*CoarseValue[k];
        for(l=0;l<N_CoarseDOF;l++)
        {
          G[k*N_CoarseDOF+l] += val*CoarseValue[l];
        } // end for l
      } // end for k
    } // end for j

    // Hack: to deal with D_h = {0}
    if(N_CoarseDOF == 1 && std::abs(G[0])<1e-10)
      G[0] = 1;

    // save G for later use
    memcpy(Gsave, G, N_CoarseDOF*N_CoarseDOF*sizeof(double));
    /*
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<N_CoarseDOF;k++)
        cout << j << " " << k << " " << G[N_CoarseDOF*j+k] << endl;
    */

    PCValues = FEDatabase::GetOrigElementValues(*CoarseBF,
                                                MultiIndex3D::D000);

    ChildValuesX = FEDatabase::GetOrigElementValues(*BF, MultiIndex3D::D100);
    ChildValuesY = FEDatabase::GetOrigElementValues(*BF, MultiIndex3D::D010);
    ChildValuesZ = FEDatabase::GetOrigElementValues(*BF, MultiIndex3D::D001);

    // only two right-hand sides (x and y derivative)
    memset(H, 0, N_CoarseDOF*3*sizeof(double));

    for(j=0;j<N_Points;j++)
    {
      PCValue = PCValues[j];
      ChildValueX = ChildValuesX[j];
      ChildValueY = ChildValuesY[j];
      ChildValueZ = ChildValuesZ[j];

      // calculate gradient of discrete uh in this quadrature point
      valx = 0.0;
      valy = 0.0;
      valz = 0.0;
      for(k=0;k<N_DOF;k++)
      {
        val = Values[DOF[k]];
        valx += ChildValueX[k]*val;
        valy += ChildValueY[k]*val;
        valz += ChildValueZ[k]*val;
      }

      // get gradient of exact u
      auto p = qf_orig.get_point(j);
      ExactU(p.x, p.y, p.z, exactval);

      valx -= exactval[1];
      valy -= exactval[2];
      valz -= exactval[3];

      w = qf_orig.get_weight(j);
      for(k=0;k<N_CoarseDOF;k++)
      {
        val = w*PCValue[k];
        H[k*3  ] += val*valx;
        H[k*3+1] += val*valy;
        H[k*3+2] += val*valz;
      } // end for k

      // grad-grad term
      locerror += w*(valx*valx + valy*valy + valz*valz);

    } // end for j
    memcpy(P, H, N_CoarseDOF*3*sizeof(double));

    /*
    cout << "vor" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    SolveMultipleSystemsNew(G, H, N_CoarseDOF, N_CoarseDOF, 3, 3);

    /*
    cout << "nach" << endl;
    for(j=0;j<N_CoarseDOF;j++)
      for(k=0;k<2*N_DOF;k++)
        cout << j << " " << k << " " << H[j*2*N_DOF+k] << endl;
    */

    // proj-proj coupling
    s = 0;
    for(i1=0;i1<N_CoarseDOF;i1++)
      for(i2=0;i2<N_CoarseDOF;i2++)
        s += Gsave[i1*N_CoarseDOF+i2]*( H[i1*3  ]*H[i2*3  ]
                                       +H[i1*3+1]*H[i2*3+1]
                                       +H[i1*3+2]*H[i2*3+2]);
    locerror += s;

    // grad-proj coupling
    s = 0;
//    for(i3=0;i3<N_CoarseDOF;i3++)
//    {
//      s += P[i3*3  ] * H[i3*3  ];
//      s += P[i3*3+1] * H[i3*3+1];
//      s += P[i3*3+2] * H[i3*3+2];
//    }
    for(i2=0;i2<N_CoarseDOF;i2++)
    {
      s += P[i2*3  ] * H[i2*3  ];
      s += P[i2*3+1] * H[i2*3+1];
      s += P[i2*3+2] * H[i2*3+2];
    }
    for(i1=0;i1<N_CoarseDOF;i1++)
    {
      s += P[i1*3  ] * H[i1*3  ];
      s += P[i1*3+1] * H[i1*3+1];
      s += P[i1*3+2] * H[i1*3+2];
    }
    locerror -= s;

    /*
    cout << "Matrix: " << endl;
     for(l=0;l<N_DOF;l++)
        for(m=0;m<N_DOF;m++)
          cout << l << " " << m << " " << LocMatrix[l*N_DOF+m] << endl;
    */

    error += lpcoeff*std::pow(hK,lpexponent)*locerror;
  } // endfor i

  // error = 0;
  cout << " end of local error " << endl;
  return error;
} // UltraLocalProjection


#endif
