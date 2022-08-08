// =======================================================================
// @(#)Vertex.C        1.1 10/30/98
// 
// Class:       TVertex
// Purpose:     a vertex in a grid 
//
// Author:      Volker Behns  09.07.97
//              Sashikumaar Ganesan 05.11.09 (added parallel methods)
// =======================================================================
#ifdef _MPI
#  include "mpi.h"
#endif

#include <Vertex.h>
#include <string.h>


#ifdef _MPI
int TVertex::NbPeriodicVert = 0;
#endif

// Constructors

TVertex::TVertex(double initX, double initY, double initZ)
{
  X = initX;
  Y = initY;
  Z = initZ;

#ifdef _MPI
  N_Cells = 0;
  SubDomainVert = false;
  N_SubDomains = 0;
  SubDomain_Ranks = nullptr;
  SubDomainGlobalCellNo = nullptr;
  SubDomainLocVertNo = nullptr;
  N_CrossNeibCells = 0;
  Cells = nullptr;
  physical_id = 0;
#endif
  BoundVert=false; 
}
TVertex::TVertex(double initX, double initY)
{
  X = initX;
  Y = initY;
}

// Methods

// set coordinates
void TVertex::SetCoords(double initX, double initY, double initZ)
{
  X = initX;
  Y = initY;
  Z = initZ;
}
void TVertex::SetCoords(double initX, double initY)
{
  X = initX;
  Y = initY;
}

#ifndef __3D__
std::ostream& operator << (std::ostream& s, const TVertex *v)
{
  return s << " X= " << setw(12) << v->X << ", Y= " << setw(12) << v->Y;
}
#else
std::ostream& operator << (std::ostream& s, const TVertex *v)
{
  return s << " X= " << setw(12) << v->X << ", Y= " << setw(12) << v->Y
           << ", Z= " << setw(12) << v->Z;
}
#endif

bool operator < (const TVertex& V, const TVertex& W)
{

  double V_X = V.GetX();
  double W_X = W.GetX();
  double V_Y = V.GetY();
  double W_Y = W.GetY();
#ifdef __3D__
  double V_Z = V.GetZ();
  double W_Z = W.GetZ();
#endif

  double tol = 1e-8;

  if(V_X - W_X < -tol)
    return true;
  else if(V_X-W_X > tol)
    return false;
  else
  {//x regarded as same, compare y
    if(V_Y - W_Y < -tol)
      return true;
    else if(V_Y - W_Y > tol)
      return false;
#ifdef __2D__
    else
      return false;
#endif
#ifdef __3D__
    else
    {
      if(V_Z - W_Z < -tol)
        return true;
      else if(V_Z - W_Z > tol)
        return false;
      else
        return false;
    }
#endif

  }


}

#ifdef _MPI

void TVertex::SetVertexCells(int n_Cells, TBaseCell **cells)
{

 int i;

  N_Cells = n_Cells;
  Cells = new TBaseCell*[N_Cells];

  for(i=0;i<N_Cells;i++)
   Cells[i] = cells[i];
}


/** add only one cell (lowest index cell) from each subdomain as a Hallo cell for 
this vertex */
void TVertex::SetSubDomainInfo(int n_SubDomains, int *subDomain_Ranks, int *subDomainGlobalCellNo, 
                               int *subDomainLocVertNo)
{
 int i;

   if(N_SubDomains)
    delete [] SubDomain_Ranks;

   N_SubDomains = n_SubDomains;
   SubDomain_Ranks = new int[N_SubDomains];
   SubDomainGlobalCellNo = new int[N_SubDomains];
   SubDomainLocVertNo = new int[N_SubDomains];


   for(i=0;i<N_SubDomains;i++)
    {
     SubDomain_Ranks[i] = subDomain_Ranks[i];
     SubDomainGlobalCellNo[i] = subDomainGlobalCellNo[i];
     SubDomainLocVertNo[i] = subDomainLocVertNo[i];
    }
}

void TVertex::AddCrossNeib(int Neib_ID)
{
 int i, tmp;
 
  CrossVert = true;
 
    for(i=N_CrossNeibCells;i<N_SubDomains;i++)
     if(SubDomain_Ranks[i] == Neib_ID) // swap, cross ID will be at the begning
     {
      tmp = SubDomain_Ranks[N_CrossNeibCells];
      SubDomain_Ranks[N_CrossNeibCells] = SubDomain_Ranks[i];
      SubDomain_Ranks[i] = tmp;

      tmp = SubDomainGlobalCellNo[N_CrossNeibCells];
      SubDomainGlobalCellNo[N_CrossNeibCells] = SubDomainGlobalCellNo[i];
      SubDomainGlobalCellNo[i] = tmp;

      tmp = SubDomainLocVertNo[N_CrossNeibCells];
      SubDomainLocVertNo[N_CrossNeibCells] = SubDomainLocVertNo[i];
      SubDomainLocVertNo[i] = tmp;

      N_CrossNeibCells++;
     }
    
}

#endif



TVertex::~TVertex()
{

#ifdef _MPI
//  if(N_Cells)
//    delete [] Cells; // cells itself will be deleted seperately
// 
//  if(N_SubDomains)
//   {
//    delete [] SubDomain_Ranks;
//    delete [] SubDomainGlobalCellNo;
//    delete [] SubDomainLocVertNo;
//   }
#endif

}
