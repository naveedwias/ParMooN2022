/**
 * @file GridTransfer.C
 *
 * Implements grid transfer operations declared in GridTransfer.h.
 *
 * TODO Enable transfering multiple functions at once.
 * TODO Re-enable MPI.
 *
 * @date 2016/05/10
 * @author Clemens Bartsch
 */

#include <GridTransfer.h>
#include <MooNMD_Io.h>

#include <FESpace2D.h>
#include <FESpace3D.h>
#include "FEDatabase.h"
#include "BaseCell.h"

#include <cstring>

#ifdef _MPI
#include<ParFECommunicator3D.h>
#endif

#ifdef __2D__

void GridTransfer::Prolongate(
    const TFESpace2D& CoarseSpace, const TFESpace2D& FineSpace,
    const double* CoarseFunction, size_t n_coarse_dofs,
    double* FineFunction, size_t n_fine_dofs
)
{
  int i,j,k,l;
  const TBaseCell *cell, *parent;
  FE_type CoarseId, FineId;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_FineDOFs, N_CoarseDOFs;
  int FineNumber, CoarseNumber;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
  double s;
  double Val[MaxN_BaseFunctions2D];
  double Val2[MaxN_BaseFunctions2D];
  int Index;
  double *entry;

  // begin code
  auto CoarseColl = CoarseSpace.GetCollection();
  N_CoarseCells = CoarseColl->GetN_Cells();
  N_CoarseDOFs = CoarseSpace.get_n_dof();

  auto FineColl = FineSpace.GetCollection();
  N_FineCells = FineColl->GetN_Cells();
  N_FineDOFs = FineSpace.get_n_dof();

  //Check if the dimensions fit the vector lengths.
  if ((int)n_coarse_dofs != N_CoarseDOFs)
    ErrThrow("Incorrect length of CoarseFunction: ",
             n_coarse_dofs, " != ", N_CoarseDOFs);
  if ((int)n_fine_dofs != N_FineDOFs)
    ErrThrow("Incorrect length of FineFunction: ",
             n_fine_dofs, " != ", N_FineDOFs);

  double* aux = new double[N_FineDOFs];
  memset(aux, 0, sizeof(double)*N_FineDOFs);

  memset(FineFunction, 0, sizeof(double)*N_FineDOFs);

  // set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);
  }

  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }

  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }

  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if (k == -2)
    {
      // cell was already handled
      continue;
    }

    if (k<=-10)
    {
      parent = cell->GetParent();
      N_Children = parent->GetN_Children();
      CoarseNumber = parent->GetClipBoard();
      CoarseId = CoarseSpace.get_fe_type(CoarseNumber);

      auto CoarseElement = CoarseSpace.get_fe(CoarseNumber);
      N_Coarse = CoarseElement.GetN_DOF();

      Ref = parent->GetRefDesc()->GetType();

      auto CoarseDOF = CoarseSpace.GetGlobalDOF(CoarseNumber);

      for(l=0;l<N_Coarse;l++)
        Val[l] = CoarseFunction[CoarseDOF[l]];

      CoarseElement.GetBaseFunct()->ChangeBF(CoarseColl, parent, Val);

      for(j=0;j<N_Children;j++)
      {
        // cout << "child: " << j << endl;
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        FineId = FineSpace.get_fe_type(FineNumber);
        auto FineElement = FineSpace.get_fe(FineNumber);
        N_Fine = FineElement.GetN_DOF();

        // do prolongation
        QQ = FEDatabase::GetProlongationMatrix2D
                (CoarseId, Ref, FineId, j);

        auto FineDOF = FineSpace.GetGlobalDOF(FineNumber);

        for(k=0;k<N_Fine;k++)
        {
          s = 0;
          entry = QQ+k*MaxN_BaseFunctions2D;
          for(l=0;l<N_Coarse;l++)
          {
            // s += QQ[k*MaxN_BaseFunctions2D+l]*Val[l];
            s += entry[l] * Val[l];
          } // endfor l
          Val2[k] = s;
        } // endfor k

        FineElement.GetBaseFunct()->ChangeBF(FineColl, cell, Val2);

        for(k=0;k<N_Fine;k++)
        {
          Index = FineDOF[k];
          FineFunction[Index] += Val2[k];
          aux[Index] += 1;
        }
      } // endfor j
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      FineId = FineSpace.get_fe_type(i);
      auto FineElement = FineSpace.get_fe(i);
      N_Fine = FineElement.GetN_DOF();

      Ref = NoRef;

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace.get_fe_type(CoarseNumber);

      auto CoarseElement = CoarseSpace.get_fe(CoarseNumber);
      N_Coarse = CoarseElement.GetN_DOF();

      // do prolongation
      QQ = FEDatabase::GetProlongationMatrix2D
              (CoarseId, Ref, FineId, 0);

      auto FineDOF = FineSpace.GetGlobalDOF(FineNumber);
      auto CoarseDOF = CoarseSpace.GetGlobalDOF(CoarseNumber);

      for(l=0;l<N_Coarse;l++)
        Val[l] = CoarseFunction[CoarseDOF[l]];

      CoarseElement.GetBaseFunct()->ChangeBF(CoarseColl, cell, Val);

      for(k=0;k<N_Fine;k++)
      {
        s = 0;
        for(l=0;l<N_Coarse;l++)
        {
          s += QQ[k*MaxN_BaseFunctions2D+l]*Val[l];
        } // endfor l
        Val2[k] = s;
      } // endfor k

      FineElement.GetBaseFunct()->ChangeBF(FineColl, cell, Val2);
      for(k=0;k<N_Fine;k++)
      {
        Index = FineDOF[k];
        FineFunction[Index] += Val2[k];
        aux[Index] += 1;
      }
    } // endelse
  } // endfor i

  for(i=0;i<N_FineDOFs;i++)
  {
    FineFunction[i] /= aux[i];
  }

  delete[] aux;
}


/** defect restriction from level+1 to level */
void GridTransfer::DefectRestriction(
    const TFESpace2D& CoarseSpace, const TFESpace2D& FineSpace,
    double* CoarseFunction, size_t n_coarse_dofs,
    const double* FineFunction, size_t n_fine_dofs)
{
  int i,j,k,l;
  const TBaseCell *cell, *parent;
  FE_type CoarseId, FineId;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_FineDOFs, N_CoarseDOFs;
  int FineNumber, CoarseNumber;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
  double s;
  double Val[MaxN_BaseFunctions2D];
  double Val2[MaxN_BaseFunctions2D];
  int Index;

  // begin code
  auto CoarseColl = CoarseSpace.GetCollection();
  N_CoarseCells = CoarseColl->GetN_Cells();
  N_CoarseDOFs = CoarseSpace.get_n_dof();

  auto FineColl = FineSpace.GetCollection();
  N_FineCells = FineColl->GetN_Cells();
  N_FineDOFs = FineSpace.get_n_dof();

  //Check if the dimensions fit the vector lengths.
  if ((int)n_coarse_dofs != N_CoarseDOFs)
    ErrThrow("Incorrect length of CoarseFunction: ",
             n_coarse_dofs, " != ", N_CoarseDOFs);
  if ((int)n_fine_dofs != N_FineDOFs)
    ErrThrow("Incorrect length of FineFunction: ",
             n_fine_dofs, " != ", N_FineDOFs);

  //Make a working copy of FineFunction
  std::vector<double> FineFunction_copy(N_FineDOFs);
  std::copy(FineFunction, FineFunction+N_FineDOFs, FineFunction_copy.begin());

  double* aux = new double[N_FineDOFs];
  memset(aux, 0, sizeof(double)*N_FineDOFs);

  memset(CoarseFunction, 0, sizeof(double)*N_CoarseDOFs);

  // set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);
    auto DOF = FineSpace.GetGlobalDOF(i);
    N_Fine = FineSpace.get_n_local_dof(i);
    for(j=0;j<N_Fine;j++)
      aux[DOF[j]] += 1;
  }

  // modify fine function values, will be repaired at end
  for(i=0;i<N_FineDOFs;i++)
  {
    FineFunction_copy[i] /= aux[i];
  }

  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }

  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }

  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if (k == -2)
    {
      // cell was already handled
      continue;
    }

    if (k<=-10)
    {
      parent = cell->GetParent();
      N_Children = parent->GetN_Children();
      CoarseNumber = parent->GetClipBoard();
      CoarseId = CoarseSpace.get_fe_type(CoarseNumber);

      auto CoarseElement = CoarseSpace.get_fe(CoarseNumber);
      N_Coarse = CoarseElement.GetN_DOF();

      Ref = parent->GetRefDesc()->GetType();

      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        FineId = FineSpace.get_fe_type(FineNumber);
        auto FineElement = FineSpace.get_fe(FineNumber);
        N_Fine = FineElement.GetN_DOF();

        // do restriction
        QQ = FEDatabase::GetProlongationMatrix2D
                (CoarseId, Ref, FineId, j);

        auto FineDOF = FineSpace.GetGlobalDOF(FineNumber);
        auto CoarseDOF = CoarseSpace.GetGlobalDOF(CoarseNumber);

        for(l=0;l<N_Fine;l++)
          Val[l] = FineFunction_copy[FineDOF[l]];

        FineElement.GetBaseFunct()->ChangeBF(FineColl, cell, Val);

        for(k=0;k<N_Coarse;k++)
        {
          s = 0;
          for(l=0;l<N_Fine;l++)
          {
            s += QQ[l*MaxN_BaseFunctions2D+k] * Val[l];
          } // endfor l
          Val2[k] = s;
        } // endfor k

        CoarseElement.GetBaseFunct()->ChangeBF(CoarseColl, parent, Val2);

        for(k=0;k<N_Coarse;k++)
        {
          Index = CoarseDOF[k];
          CoarseFunction[Index] += Val2[k];
        }
      } // endfor j
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      auto FineElement = FineSpace.get_fe(i);
      FineId = FineElement.GetID();
      N_Fine = FineElement.GetN_DOF();

      CoarseNumber = k;
      FineNumber = i;

      auto CoarseElement = CoarseSpace.get_fe(CoarseNumber);
      CoarseId = CoarseElement.GetID();
      N_Coarse = CoarseElement.GetN_DOF();

      Ref = NoRef;

      // do restriction
      QQ = FEDatabase::GetProlongationMatrix2D
              (CoarseId, Ref, FineId, 0);

      auto FineDOF = FineSpace.GetGlobalDOF(FineNumber);
      auto CoarseDOF = CoarseSpace.GetGlobalDOF(CoarseNumber);

      for(l=0;l<N_Fine;l++)
        Val[l] = FineFunction_copy[FineDOF[l]];

      FineElement.GetBaseFunct()->ChangeBF(FineColl, cell, Val);

      for(k=0;k<N_Coarse;k++)
      {
        s = 0;
        for(l=0;l<N_Fine;l++)
        {
          s += QQ[l*MaxN_BaseFunctions2D+k]*Val[l];
        } // endfor l
        Val2[k] = s;
      } // endfor k

      CoarseElement.GetBaseFunct()->ChangeBF(CoarseColl, cell, Val2);

      for(k=0;k<N_Coarse;k++)
        CoarseFunction[CoarseDOF[k]] += Val2[k];
    } // endelse
  } // endfor i

  delete[] aux;
}

/** function restriction from level+1 to level */
void GridTransfer::RestrictFunction(
    const TFESpace2D& CoarseSpace, const TFESpace2D& FineSpace,
    double* CoarseFunction, size_t n_coarse_dofs,
    const double* FineFunction, size_t n_fine_dofs)
{
  int i,j,k,l;
  const TBaseCell *cell, *parent;
  FE_type CoarseId, FineId;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_FineDOFs, N_CoarseDOFs;
  int FineNumber, CoarseNumber;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
  double s;
  double Val[MaxN_BaseFunctions2D];
  double Val2[MaxN_BaseFunctions2D];

  // begin code
  auto CoarseColl = CoarseSpace.GetCollection();
  N_CoarseCells = CoarseColl->GetN_Cells();
  N_CoarseDOFs = CoarseSpace.get_n_dof();

  auto FineColl = FineSpace.GetCollection();
  N_FineCells = FineColl->GetN_Cells();
  N_FineDOFs = FineSpace.get_n_dof();

  //Check if the dimensions fit the vector lengths.
  if ((int)n_coarse_dofs != N_CoarseDOFs)
    ErrThrow("Incorrect length of CoarseFunction: ",
             n_coarse_dofs, " != ", N_CoarseDOFs);
  if ((int)n_fine_dofs != N_FineDOFs)
    ErrThrow("Incorrect length of FineFunction: ",
             n_fine_dofs, " != ", N_FineDOFs);

  double* aux = new double[N_CoarseDOFs];
  memset(aux, 0, sizeof(double)*N_CoarseDOFs);

  memset(CoarseFunction, 0, sizeof(double)*n_coarse_dofs);

  // set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);
  }

  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }

  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }

  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    // cout << "i= " << i << "    ";
    // cout << "k= " << k << endl;
    if (k == -2)
    {
      // cell was already handled
      continue;
    }

    if (k<=-10)
    {
      parent = cell->GetParent();
      N_Children = parent->GetN_Children();
      CoarseNumber = parent->GetClipBoard();

      auto CoarseElement = CoarseSpace.get_fe(CoarseNumber);
      CoarseId = CoarseElement.GetID();
      N_Coarse = CoarseElement.GetN_DOF();

      Ref = parent->GetRefDesc()->GetType();

      memset(Val2, 0, MaxN_BaseFunctions2D*sizeof(double));

      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        auto FineElement = FineSpace.get_fe(FineNumber);
        FineId = FineElement.GetID();
        N_Fine = FineElement.GetN_DOF();

        // do restriction
        QQ = FEDatabase::GetRestrictionMatrix2D
                (CoarseId, Ref, FineId, j);

        auto FineDOF = FineSpace.GetGlobalDOF(FineNumber);

        for(l=0;l<N_Fine;l++)
          Val[l] = FineFunction[FineDOF[l]];

        FineElement.GetBaseFunct()->ChangeBF(FineColl, cell, Val);

        for(k=0;k<N_Coarse;k++)
        {
          s = 0;
          for(l=0;l<N_Fine;l++)
          {
            s += QQ[k*MaxN_BaseFunctions2D+l] * Val[l];
          } // endfor l
          Val2[k] += s;
        } // endfor k
      } // endfor j

      CoarseElement.GetBaseFunct()->ChangeBF(CoarseColl, parent, Val2);

      auto CoarseDOF = CoarseSpace.GetGlobalDOF(CoarseNumber);
      for(k=0;k<N_Coarse;k++)
      {
        l=CoarseDOF[k];
        aux[l] += 1;
        CoarseFunction[l] += Val2[k];
      } // endfor k
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      auto FineElement = FineSpace.get_fe(i);
      FineId = FineElement.GetID();
      N_Fine = FineElement.GetN_DOF();
      
      CoarseNumber = k;
      FineNumber = i;

      auto CoarseElement = CoarseSpace.get_fe(CoarseNumber);
      CoarseId = CoarseElement.GetID();
      N_Coarse = CoarseElement.GetN_DOF();

      Ref = NoRef;

      // do restriction
      QQ = FEDatabase::GetRestrictionMatrix2D
              (CoarseId, Ref, FineId, 0);

      auto FineDOF = FineSpace.GetGlobalDOF(FineNumber);
      auto CoarseDOF = CoarseSpace.GetGlobalDOF(CoarseNumber);

      for(l=0;l<N_Fine;l++)
        Val[l] = FineFunction[FineDOF[l]];

      FineElement.GetBaseFunct()->ChangeBF(FineColl, cell, Val);

      for(k=0;k<N_Coarse;k++)
      {
        s = 0;
        for(l=0;l<N_Fine;l++)
        {
          s += QQ[k*MaxN_BaseFunctions2D+l]*Val[l];
        } // endfor l
        Val2[k] = s;
      } // endfor k

      CoarseElement.GetBaseFunct()->ChangeBF(CoarseColl, cell, Val2);

      for(k=0;k<N_Coarse;k++)
      {
        l=CoarseDOF[k];
        CoarseFunction[l] += Val2[k];
        aux[l] += 1;
      } // endfor k
    } // endelse
  } // endfor i

  for(i=0;i<N_CoarseDOFs;i++)
    CoarseFunction[i] /= aux[i];

  delete[] aux;
} // RestrictFunction

#endif
#ifdef __3D__
/** prolongate */
void GridTransfer::Prolongate(
    const TFESpace3D& CoarseSpace, const TFESpace3D& FineSpace,
    const double* CoarseFunction, size_t n_coarse_dofs,
    double* FineFunction, size_t n_fine_dofs)

{
  int i,j,k,l;
  const TBaseCell *cell, *parent;
  FE_type CoarseId, FineId;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_FineDOFs; // N_CoarseDOFs;
  int FineNumber, CoarseNumber;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
//  double *CurrentCoarseFct, *CurrentFineFct;
  double s;
  double Val[MaxN_BaseFunctions3D];
  double Val2[MaxN_BaseFunctions3D];
  int Index;
  double *entry;
  // begin code
  auto CoarseColl = CoarseSpace.GetCollection();
  N_CoarseCells = CoarseColl->GetN_Cells();

  auto FineColl = FineSpace.GetCollection();
  N_FineCells = FineColl->GetN_Cells();
  N_FineDOFs = FineSpace.get_n_dof();

  //Check if the dimensions fit the vector lengths.
  if ((int)n_coarse_dofs != CoarseSpace.get_n_dof())
    ErrThrow("Incorrect length of CoarseFunction: ",
             n_coarse_dofs, " != ", CoarseSpace.get_n_dof());
  if ((int)n_fine_dofs != N_FineDOFs)
    ErrThrow("Incorrect length of FineFunction: ",
             n_fine_dofs, " != ", N_FineDOFs);

  //cout << "N_FineCells: " << N_FineCells << endl;
  //cout << "N_CoarseCells: " << N_CoarseCells << endl;

  double* aux = new double[N_FineDOFs];
  memset(aux, 0, sizeof(double)*N_FineDOFs);

  memset(FineFunction, 0, sizeof(double)*N_FineDOFs);

#ifdef _OMP
#pragma omp parallel default(shared) private(i,j,k,l,cell,DOF,FineElement,FineBF,N_Fine, \
                                             parent, N_Children, CoarseNumber, CoarseElement, \
                                             CoarseBF, BaseFunctions, Ref, FineNumber,\
                                             QQ, FineDOF, Val, s, Val2, \
                                             Index, N_Coarse, entry)
{
#pragma omp for schedule(guided)
#endif
  // set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);

    auto DOF = FineSpace.GetGlobalDOF(i);
    auto FineElement = FineSpace.get_fe(i);
    N_Fine = FineElement.GetN_DOF();
    for(j=0;j<N_Fine;j++)
#ifdef _OMP
      #pragma omp atomic
#endif
      aux[DOF[j]] += 1;
  }

#ifdef _OMP
#pragma omp for schedule(guided)
#endif
  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }


#ifdef _OMP
#pragma omp for schedule(guided)
#endif
  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }

#ifdef _OMP
#pragma omp for schedule(guided)
#endif
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
#ifdef _MPI
    if(cell->IsHaloCell())   continue;
#endif

    k = cell->GetClipBoard();
    if (k == -2)
    {
      // cell was already handled
      continue;
    }

    if (k<=-10)
    {
      parent = cell->GetParent();
      N_Children = parent->GetN_Children();
      CoarseNumber = parent->GetClipBoard();
      CoarseId = CoarseSpace.get_fe_type(CoarseNumber);

      auto CoarseElement = CoarseSpace.get_fe(CoarseNumber);
      auto BaseFunctions = CoarseElement.GetBaseFunct();
      N_Coarse = BaseFunctions->GetDimension();

      Ref = parent->GetRefDesc()->GetType();

      auto CoarseDOF = CoarseSpace.GetGlobalDOF(CoarseNumber);

      for(l=0;l<N_Coarse;l++)
        Val[l] = CoarseFunction[CoarseDOF[l]];

      BaseFunctions->ChangeBF(CoarseColl, parent, Val);

      for(j=0;j<N_Children;j++)
      {
        // cout << "child: " << j << endl;
        cell = parent->GetChild(j);
// #ifdef _OMP
//  #pragma omp critical
// #endif
        {
         k = cell->GetClipBoard();
         cell->SetClipBoard(-2);
        }

#ifdef _OMP
    if(k==-2) continue;
#endif
  FineNumber = -(k+10);
        FineId = FineSpace.get_fe_type(FineNumber);
        auto FineElement = FineSpace.get_fe(FineNumber);
        N_Fine = FineElement.GetN_DOF();

#ifdef _OMP
       #pragma omp critical
#endif
  {
    QQ = FEDatabase::GetProlongationMatrix3D(CoarseId, Ref, FineId, j);
  }

        auto FineDOF = FineSpace.GetGlobalDOF(FineNumber);

        for(k=0;k<N_Fine;k++)
        {
          s = 0;
          entry = QQ+k*MaxN_BaseFunctions3D;
          for(l=0;l<N_Coarse;l++)
          {
            // s += QQ[k*MaxN_BaseFunctions3D+l]*Val[l];
            s += entry[l] * Val[l];
            // cout << k << " " << l << " " << entry[l] << endl;
          } // endfor l
          Val2[k] = s;
        } // endfor k

        FineElement.GetBaseFunct()->ChangeBF(FineColl, cell, Val2);

        for(k=0;k<N_Fine;k++)
        {
          Index = FineDOF[k];
   {
#ifdef _OMP
       #pragma omp atomic
#endif
          FineFunction[Index] += Val2[k];
   }
        }
      } // endfor j
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      FineId = FineSpace.get_fe_type(i);
      auto FineElement = FineSpace.get_fe(i);
      N_Fine = FineElement.GetN_DOF();

      Ref = NoRef;

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace.get_fe_type(CoarseNumber);

      auto CoarseElement = CoarseSpace.get_fe(CoarseNumber);
      auto BaseFunctions = CoarseElement.GetBaseFunct();
      N_Coarse = BaseFunctions->GetDimension();

#ifdef _OMP
       #pragma omp critical
#endif
  {
    QQ = FEDatabase::GetProlongationMatrix3D(CoarseId, Ref, FineId, 0);
  }

      auto FineDOF = FineSpace.GetGlobalDOF(FineNumber);
      auto CoarseDOF = CoarseSpace.GetGlobalDOF(CoarseNumber);

      for(l=0;l<N_Coarse;l++)
        Val[l] = CoarseFunction[CoarseDOF[l]];

      BaseFunctions->ChangeBF(CoarseColl, cell, Val);

      for(k=0;k<N_Fine;k++)
      {
        s = 0;
        for(l=0;l<N_Coarse;l++)
        {
          s += QQ[k*MaxN_BaseFunctions3D+l]*Val[l];
        } // endfor l
        Val2[k] = s;
      } // endfor k

      FineElement.GetBaseFunct()->ChangeBF(FineColl, cell, Val2);

      for(k=0;k<N_Fine;k++)
      {
        Index = FineDOF[k];
  {
#ifdef _OMP
       #pragma omp atomic
#endif
         FineFunction[Index] += Val2[k];
  }
      }
    } // endelse
  } // endfor i

#ifdef _OMP
#pragma omp for schedule(guided)
#endif
  for(i=0;i<N_FineDOFs;i++)
  {
    FineFunction[i] /= aux[i];
  }
#ifdef _OMP
}
#endif

  delete[] aux;
  
#ifdef _MPI
  /// @todo why is this needed here? 
  FineSpace.get_communicator().CommUpdateReduce(FineFunction);
#endif
}

void GridTransfer::DefectRestriction(
    const TFESpace3D& CoarseSpace, const TFESpace3D& FineSpace,
    double* CoarseFunction, size_t n_coarse_dofs,
    const double* FineFunction, size_t n_fine_dofs)
{
  int i,j,k,l;
  const TBaseCell *cell, *parent;
  FE_type CoarseId, FineId;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_FineDOFs, N_CoarseDOFs;
  int FineNumber, CoarseNumber;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
//  double *CurrentCoarseFct, *CurrentFineFct;
  double s;
  double Val[MaxN_BaseFunctions3D];
  double Val2[MaxN_BaseFunctions3D];
  int Index;

  // begin code
  auto CoarseColl = CoarseSpace.GetCollection();
  N_CoarseCells = CoarseColl->GetN_Cells();
  N_CoarseDOFs = CoarseSpace.get_n_dof();

  auto FineColl = FineSpace.GetCollection();
  N_FineCells = FineColl->GetN_Cells();
  N_FineDOFs = FineSpace.get_n_dof();

  //Check if the dimensions fit the vector lengths.
  if ((int)n_coarse_dofs != N_CoarseDOFs)
    ErrThrow("Incorrect length of CoarseFunction: ",
             n_coarse_dofs, " != ", N_CoarseDOFs);
  if ((int)n_fine_dofs != N_FineDOFs)
    ErrThrow("Incorrect length of FineFunction: ",
             n_fine_dofs, " != ", N_FineDOFs);

  //Make a working copy of FineFunction
  std::vector<double> FineFunction_copy(N_FineDOFs);
  std::copy(FineFunction, FineFunction+N_FineDOFs, FineFunction_copy.begin());

  // cout << "N_FineCells: " << N_FineCells << endl;
  //cout << "N_CoarseCells: " << N_CoarseCells << endl;

  double* aux = new double[N_FineDOFs];
  memset(aux, 0, sizeof(double)*N_FineDOFs);

  memset(CoarseFunction, 0, sizeof(double)*N_CoarseDOFs);


#ifdef _OMP
#pragma omp parallel default(shared) private(i,j,k,l,cell,DOF,FineId,FineElement,N_Fine, \
                                             parent, N_Children, CoarseNumber, CoarseId, CoarseElement, \
                                             BaseFunctions, Ref, FineNumber,\
                                             QQ, FineDOF, CoarseDOF, Val, s, Val2, \
                                             Index, N_Coarse)
{
#pragma omp for schedule(guided)
#endif
//   set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);

    auto DOF = FineSpace.GetGlobalDOF(i);
    FineId = FineSpace.get_fe_type(i);
    auto FineElement = FineSpace.get_fe(i);
    N_Fine = FineElement.GetN_DOF();
    for(j=0;j<N_Fine;j++)
#ifdef _OMP
      #pragma omp atomic
#endif
      aux[DOF[j]] += 1;
  }

#ifdef _OMP
#pragma omp for schedule(guided) nowait
#endif

  for(i=0;i<N_FineDOFs;i++)
  {
    FineFunction_copy[i] /= aux[i];
  }

#ifdef _OMP
#pragma omp for schedule(guided)
#endif
  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }

#ifdef _OMP
#pragma omp for schedule(guided)
#endif
  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }

#ifdef _OMP
#pragma omp for schedule(guided)
#endif
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
#ifdef _MPI
    if(cell->IsHaloCell())  continue;
#endif

    k = cell->GetClipBoard();
    if (k == -2)
    {
      // cell was already handled
      continue;
    }

    if (k<=-10)
    {
      parent = cell->GetParent();
      N_Children = parent->GetN_Children();
      CoarseNumber = parent->GetClipBoard();
      CoarseId = CoarseSpace.get_fe_type(CoarseNumber);
      auto CoarseElement = CoarseSpace.get_fe(CoarseNumber);
      N_Coarse = CoarseSpace.get_n_local_dof(CoarseNumber);

      Ref = parent->GetRefDesc()->GetType();

      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
// #ifdef _OMP
//  #pragma omp critical
// #endif
  {
         k = cell->GetClipBoard();
         cell->SetClipBoard(-2);
  }

#ifdef _OMP
    if(k==-2) continue;
#endif

  FineNumber = -(k+10);
        FineId = FineSpace.get_fe_type(FineNumber);
        auto FineElement = FineSpace.get_fe(FineNumber);
        N_Fine = FineSpace.get_n_local_dof(FineNumber);
#ifdef _OMP
       #pragma omp critical
#endif
  {
    QQ = FEDatabase::GetProlongationMatrix3D(CoarseId, Ref, FineId, j);
  }

        auto FineDOF = FineSpace.GetGlobalDOF(FineNumber);
        auto CoarseDOF = CoarseSpace.GetGlobalDOF(CoarseNumber);

        for(l=0;l<N_Fine;l++)
          Val[l] = FineFunction_copy[FineDOF[l]];

        FineElement.GetBaseFunct()->ChangeBF(FineColl, cell, Val);

        for(k=0;k<N_Coarse;k++)
        {
          s = 0;
          for(l=0;l<N_Fine;l++)
          {
            s += QQ[l*MaxN_BaseFunctions3D+k] * Val[l];
          } // endfor l
          Val2[k] = s;
        } // endfor k

        CoarseElement.GetBaseFunct()->ChangeBF(CoarseColl, parent, Val2);

        for(k=0;k<N_Coarse;k++)
        {
          Index = CoarseDOF[k];
#ifdef _OMP
       #pragma omp atomic
#endif
           CoarseFunction[Index] += Val2[k];
        }
      } // endfor j
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      FineId = FineSpace.get_fe_type(i);
      auto FineElement = FineSpace.get_fe(i);
      N_Fine = FineElement.GetN_DOF();

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace.get_fe_type(CoarseNumber);
      auto CoarseElement = CoarseSpace.get_fe(CoarseNumber);
      N_Coarse = CoarseElement.GetN_DOF();

      Ref = NoRef;
#ifdef _OMP
       #pragma omp critical
#endif
  {
    QQ = FEDatabase::GetProlongationMatrix3D(CoarseId, Ref, FineId, 0);
  }

      auto FineDOF = FineSpace.GetGlobalDOF(FineNumber);
      auto CoarseDOF = CoarseSpace.GetGlobalDOF(CoarseNumber);

      for(l=0;l<N_Fine;l++)
        Val[l] = FineFunction_copy[FineDOF[l]];

      FineElement.GetBaseFunct()->ChangeBF(FineColl, cell, Val);

      for(k=0;k<N_Coarse;k++)
      {
        s = 0;
        for(l=0;l<N_Fine;l++)
        {
          s += QQ[l*MaxN_BaseFunctions3D+k]*Val[l];
        } // endfor l
        Val2[k] = s;
      } // endfor k

      CoarseElement.GetBaseFunct()->ChangeBF(CoarseColl, cell, Val2);

      for(k=0;k<N_Coarse;k++)
      {
#ifdef _OMP
       #pragma omp atomic
#endif
  CoarseFunction[CoarseDOF[k]] += Val2[k];
      }
    } // endelse
  } // endfor i


#ifdef _OMP
}
#endif
  delete[] aux;
}

void GridTransfer::RestrictFunction(
    const TFESpace3D& CoarseSpace, const TFESpace3D& FineSpace,
    double* CoarseFunction, size_t n_coarse_dofs,
    const double* FineFunction, size_t n_fine_dofs)
{
  int i,j,k,l;
  const TBaseCell *cell, *parent;
  FE_type CoarseId, FineId;
  int N_CoarseCells, N_FineCells, N_Children;
  int N_CoarseDOFs; //N_FineDOFs
  int FineNumber, CoarseNumber;
  int N_Fine, N_Coarse;
  Refinements Ref;
  double *QQ;
//  double *CurrentCoarseFct, *CurrentFineFct;
  double s;
  double Val[MaxN_BaseFunctions3D];
  double Val2[MaxN_BaseFunctions3D];
//  int *DOF, Index;
//  double *entry;
  // begin code
  auto CoarseColl = CoarseSpace.GetCollection();
  N_CoarseCells = CoarseColl->GetN_Cells();
  N_CoarseDOFs = CoarseSpace.get_n_dof();

  auto FineColl = FineSpace.GetCollection();
  N_FineCells = FineColl->GetN_Cells();
  
  //Check if the dimensions fit the vector lengths.
  if ((int)n_coarse_dofs != N_CoarseDOFs)
    ErrThrow("Incorrect length of CoarseFunction: ",
             n_coarse_dofs, " != ", N_CoarseDOFs);
  if ((int)n_fine_dofs != FineSpace.get_n_dof())
    ErrThrow("Incorrect length of FineFunction: ",
             n_fine_dofs, " != ", FineSpace.get_n_dof());

  double* aux = new double[N_CoarseDOFs];
  memset(aux, 0, sizeof(double)*N_CoarseDOFs);

  memset(CoarseFunction, 0, sizeof(double)*N_CoarseDOFs);

#ifdef _OMP
#pragma omp parallel default(shared) private(i,cell,k)
{
#pragma omp for schedule(static) nowait
#endif
  // set fine grid clipboard to -1
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    cell->SetClipBoard(-1);
  }
#ifdef _OMP
#pragma omp for schedule(static) nowait
#endif
  // set coarse grid clipboard to implicit number
  for(i=0;i<N_CoarseCells;i++)
  {
    cell = CoarseColl->GetCell(i);
    cell->SetClipBoard(i);
  }
#ifdef _OMP
#pragma omp for schedule(static) nowait
#endif
  // if a cell with clipboard==-1 is found
  // => this cell is only on the fine grid
  // set clipboard to "-number-10"
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
    k = cell->GetClipBoard();
    if(k==-1) cell->SetClipBoard(-i-10);
  }
#ifdef _OMP
}
#endif
  for(i=0;i<N_FineCells;i++)
  {
    cell = FineColl->GetCell(i);
#ifdef _MPI
    // Just skip fine cells which lie in the halo
    // - they will be treated on that process, where they are own cells
    if(cell->IsHaloCell())
      continue;
#endif
    k = cell->GetClipBoard();
    // cout << "i= " << i << "    ";
    // cout << "k= " << k << endl;
    if (k == -2)
    {
      // cell was already handled
      continue;
    }

    if (k<=-10)
    {
      parent = cell->GetParent();
      N_Children = parent->GetN_Children();
      CoarseNumber = parent->GetClipBoard();
      CoarseId = CoarseSpace.get_fe_type(CoarseNumber);

      auto CoarseElement = CoarseSpace.get_fe(CoarseNumber);
      N_Coarse = CoarseElement.GetN_DOF();

      Ref = parent->GetRefDesc()->GetType();

      memset(Val2, 0, MaxN_BaseFunctions3D*sizeof(double));
#ifdef _OMP
// #pragma omp parallel default(shared) private(j,k,s,l,cell,FineNumber,FineId,FineElement,FineBF,N_Fine,QQ,FineDOF,Val2,Index)
// #pragma omp for schedule(guided) nowait
#endif
      auto CoarseDOF = CoarseSpace.GetGlobalDOF(CoarseNumber);
      for(j=0;j<N_Children;j++)
      {
        cell = parent->GetChild(j);
        k = cell->GetClipBoard();
        FineNumber = -(k+10);
        cell->SetClipBoard(-2);
        FineId = FineSpace.get_fe_type(FineNumber);
        auto FineElement = FineSpace.get_fe(FineNumber);
        N_Fine = FineElement.GetN_DOF();

        // do restriction
/*
        cout << "CoarseId: " << CoarseId << endl;
        cout << "Ref: " << Ref << endl;
        cout << "FineId: " << FineId << endl;
        cout << "j: " << j << endl;
*/
        QQ = FEDatabase::GetRestrictionMatrix3D(CoarseId, Ref, FineId, j);

        auto FineDOF = FineSpace.GetGlobalDOF(FineNumber);

        for(l=0;l<N_Fine;l++)
          Val[l] = FineFunction[FineDOF[l]];

        FineElement.GetBaseFunct()->ChangeBF(FineColl, cell, Val);

        for(k=0;k<N_Coarse;k++)
        {
          s = 0;
          for(l=0;l<N_Fine;l++)
          {
            s += QQ[k*MaxN_BaseFunctions3D+l] * Val[l];
          } // endfor l
          Val2[k] += s;
        } // endfor k
      } // endfor j

      CoarseElement.GetBaseFunct()->ChangeBF(CoarseColl, parent, Val2);

#ifdef _OMP
// #pragma omp parallel default(shared) private(k,l)
// #pragma omp for schedule(guided) nowait
#endif
      for(k=0;k<N_Coarse;k++)
      {
        l=CoarseDOF[k];
        aux[l] += 1;
//   #pragma omp critical
        CoarseFunction[l] += Val2[k];
      } // endfor k
    } // endif
    else
    {
      // number in clipboard is number of fine cell in coarse grid
      FineId = FineSpace.get_fe_type(i);
      auto FineElement = FineSpace.get_fe(i);
      N_Fine = FineElement.GetN_DOF();

      CoarseNumber = k;
      FineNumber = i;
      CoarseId = CoarseSpace.get_fe_type(CoarseNumber);

      auto CoarseElement = CoarseSpace.get_fe(CoarseNumber);
      N_Coarse = CoarseElement.GetN_DOF();

      Ref = NoRef;

      // do restriction
/*
      cout << "CoarseId: " << CoarseId << endl;
      cout << "Ref: " << Ref << endl;
      cout << "FineId: " << FineId << endl;
      cout << "j: " << j << endl;
      cout << endl;
*/
      QQ = FEDatabase::GetRestrictionMatrix3D(CoarseId, Ref, FineId, 0);

      auto FineDOF = FineSpace.GetGlobalDOF(FineNumber);
      auto CoarseDOF = CoarseSpace.GetGlobalDOF(CoarseNumber);

#ifdef _OMP
// #pragma omp parallel default(shared) private(k,s,l,Val2,Index)
{
// #pragma omp for schedule(guided) nowait
#endif
      for(l=0;l<N_Fine;l++)
        Val[l] = FineFunction[FineDOF[l]];

      FineElement.GetBaseFunct()->ChangeBF(FineColl, cell, Val);

#ifdef _OMP
// #pragma omp for schedule(guided) nowait
#endif
      for(k=0;k<N_Coarse;k++)
      {
        s = 0;
        for(l=0;l<N_Fine;l++)
        {
          s += QQ[k*MaxN_BaseFunctions3D+l]*Val[l];
        } // endfor l
        Val2[k] = s;
      } // endfor k

      CoarseElement.GetBaseFunct()->ChangeBF(CoarseColl, cell, Val2);

#ifdef _OMP
// #pragma omp for schedule(guided) nowait
#endif
      for(k=0;k<N_Coarse;k++)
      {
        l=CoarseDOF[k];
#ifdef _OMP
  #pragma omp critical
#endif
  {
         CoarseFunction[l] += Val2[k];
         aux[l] += 1;
  }
      } // endfor k
#ifdef _OMP
}
#endif
    } // endelse
  } // endfor i
#ifdef _OMP
#pragma omp parallel default(shared) private(i)
#pragma omp for schedule(static) nowait
#endif

#ifdef _MPI
    // Updates and numbers are stored additive. Make the vectors consistent.
    const TParFECommunicator3D& comm = CoarseSpace.get_communicator();
    comm.update_from_additive_to_consistent_storage(CoarseFunction,0);
    comm.update_from_additive_to_consistent_storage(aux,0);
#endif

  for(i=0;i<N_CoarseDOFs;i++)
    CoarseFunction[i] /= aux[i];

/*
  for(i=0;i<N_CoarseDOFs;i++)
    cout << "CoarseFunction[" << i << "]: " << CoarseFunction[i] << endl;

  for(i=0;i<N_FineDOFs;i++)
    cout << "FineFunction[" << i << "]: " << FineFunction[i] << endl;
*/

  delete[] aux;
} // RestrictFunction

#endif


void GridTransfer::RestrictFunctionRepeatedly(
#ifdef __2D__
  const std::vector<const TFESpace2D*>& space_hierarchy,
#elif __3D__
  const std::vector<const TFESpace3D*>& space_hierarchy,
#endif
  const std::vector<double*>& function_entries,
  const std::vector<size_t>& function_n_dofs)
{
  size_t n_levels =space_hierarchy.size();
  if (n_levels < 2)
    ErrThrow("At least 2 levels needed to perform a grid restriction!");
  if(function_entries.size() != n_levels || function_n_dofs.size() != n_levels)
    ErrThrow("Number of functions does not match number of spaces!");

  for(size_t i =0; i < n_levels - 1; ++i)
  {
#ifdef __2D__
    const TFESpace2D& coarse_space = *space_hierarchy.at(i+1);
    const TFESpace2D& fine_space = *space_hierarchy.at(i);
#elif __3D__
    const TFESpace3D& coarse_space = *space_hierarchy.at(i+1);
    const TFESpace3D& fine_space = *space_hierarchy.at(i);
#endif


    double* coarse_function = function_entries.at(i+1);
    size_t coarse_n_dofs = function_n_dofs.at(i+1);


    double* fine_function = function_entries.at(i);
    size_t fine_n_dofs = function_n_dofs.at(i);

#ifdef _MPI
#ifdef __3D__
    // put the fine function into level 3 (=full) consistency
    TParFECommunicator3D fine_comm = fine_space.get_communicator();
    fine_comm.consistency_update(fine_function, 3);
#endif
#endif

    //do the restriction
    RestrictFunction(coarse_space, fine_space,
                     coarse_function, coarse_n_dofs,
                     fine_function, fine_n_dofs);

#ifdef _MPI
#ifdef __3D__
    //// Restore Level 3 consistency of the coarse solution
    // TODO Doing this here is somewhat against the concept - do it where the result of this is needed!!!
    TParFECommunicator3D coarse_comm = coarse_space.get_communicator();
    coarse_comm.consistency_update(coarse_function, 3);
#endif
#endif
  }
}
