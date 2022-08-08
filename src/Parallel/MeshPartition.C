#ifdef _MPI
// =======================================================================

// @(#)MeshPartition.C
//
// Purpose:     partition the domain into "npart" parts for parallel computing
//
// Author:      Sashikumaar Ganesan
// History:      start of implementation  07/09/09 (Sashikumaar Ganesan)
// =======================================================================
#include "mpi.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <Database.h>
#include <Domain.h>
#include <MeshPartition.h>
#include <MooNMD_Io.h>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <BaseCell.h>
#include <Edge.h>
#include <Joint.h>
#include <MacroCell.h>
#include <PeriodicDomain.h>
#include <Quadrangle.h>
#include <SubDomainHaloJoint.h>
#include <SubDomainJoint.h>
#include <Vertex.h>

#include <MeshPartitionInOut.h> // this includes metis.h

#include <algorithm>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>

// return all vertices ordered by the address in memory
std::vector<const TVertex*> get_sorted_vertices(const TCollection* coll)
{
  // this assumes that all cells have the same geometry:
  int n_cells = coll->GetN_Cells();
  int n_vertices_in_cell = coll->GetCell(0)->GetN_Vertices();
  std::vector<const TVertex*> all_vertices(n_vertices_in_cell * n_cells);
  for(int i = 0, N = 0; i < n_cells; i++)
  {
    auto cell = coll->GetCell(i);
    for(int j = 0; j < n_vertices_in_cell; j++)
    {
      all_vertices.at(N) = cell->GetVertex(j);
      N++;
    }
  }
  // we sort in reverse order here, some other program parts rely on this!
  std::sort(all_vertices.rbegin(), all_vertices.rend());
  return all_vertices;
}

// This is an implementation of finding an Element
// in a sorted array (sorted according to address, i.e., pointer arithmetic)
int get_index_of_vertex(const std::vector<const TVertex*>& v, const TVertex* e)
{
  auto cmp_fct = [](const TVertex* v1, const TVertex* v2) { return v1 > v2; };
  auto it = std::lower_bound(v.begin(), v.end(), e, cmp_fct);
  return std::distance(v.begin(), it);
}


int get_number_vertex(const std::vector<const TVertex*>& all_vertices,
                      const TCollection* coll, std::vector<int>& VertexNumbers)
{
  auto n_all_vertices = all_vertices.size();
  int n_vertices_per_cell = coll->GetCell(0)->GetN_Vertices();
  std::vector<int> number_vertex(n_all_vertices);
  const TVertex* previous_vertex = nullptr;
  int n_root_vertices = -1;
  for(size_t i = 0; i < n_all_vertices; i++)
  {
    const TVertex* current_vertex = all_vertices[i];
    if(current_vertex != previous_vertex)
    {
      n_root_vertices++;
      previous_vertex = current_vertex;
    }
    number_vertex[i] = n_root_vertices;
  }
  n_root_vertices++;

  int n_cells = coll->GetN_Cells();
  VertexNumbers.resize(n_vertices_per_cell * n_cells);
  for(int i = 0, m = 0; i < n_cells; i++)
  {
    const TBaseCell* cell = coll->GetCell(i);
    for(int j = 0; j < n_vertices_per_cell; j++)
    {
      const TVertex* current_vertex = cell->GetVertex(j);
      int N = get_index_of_vertex(all_vertices, current_vertex);
      VertexNumbers[m] = number_vertex[N];
      m++;
    }
  }
  return n_root_vertices;
}

int get_max_n_cells_per_vertex(const std::vector<int>& VertexNumbers,
                                int N_RootVertices,
                                bool is_periodic,
                                const std::vector<std::vector<int>>&
                                                          PeriodicVertexConnect)
{
  int N_AllLocVert = VertexNumbers.size();
  std::vector<int> number_of_cells_sharing_vertex(N_RootVertices, 0);
  for(int i = 0; i < N_AllLocVert; i++)
    number_of_cells_sharing_vertex.at(VertexNumbers[i]) += 1;
  
  // find maximum cells per vertex ( MaxCpV )
  int MaxCpV = *std::max_element(number_of_cells_sharing_vertex.begin(),
                                 number_of_cells_sharing_vertex.end());

  // correction in case of periodic boundaries
  // TODO improve the upper bound (it can currently be to large)
  if (is_periodic)
  {
    int max_sharing_periodic_vertex = 0;
    
    for (int i = 0 ; i < (int)PeriodicVertexConnect.size() ; i++)
    {
      int ncs_periodic_vertex = 0; // number of cells sharing a periodic vertex
      
      for(int j = 0 ; j < (int)PeriodicVertexConnect.at(i).size() ; j++)
      {
        ncs_periodic_vertex += number_of_cells_sharing_vertex.at(
                                 PeriodicVertexConnect.at(i).at(j));
      }
      
      max_sharing_periodic_vertex = std::max(max_sharing_periodic_vertex,
                                             ncs_periodic_vertex);
    }
    
    MaxCpV = std::max(MaxCpV, max_sharing_periodic_vertex);
  }

  return MaxCpV;
}

// the returned vector's first column contains number of neighbor cells
// associated with each vertex, further columns contain the cell numbers
// associated with this vertex
std::vector<int> get_neighbors_of_vertices(int N_RootVertices, int MaxCpV,
                                           int N_Cells, int N_VertInCell,
                                           const std::vector<int>& VertexNumbers,
                                           bool is_periodic,
                                           const std::vector<std::vector<int>>& PeriodicVertexConnect)
{
  std::vector<int> vertex_cell_neighbors_info(N_RootVertices * MaxCpV, 0);
  
  for(int i = 0; i < N_Cells; i++)
  {
    for(int j = 0; j < N_VertInCell; j++)
    {
      int M = VertexNumbers[i * N_VertInCell + j] * MaxCpV;
      vertex_cell_neighbors_info[M]++;
      vertex_cell_neighbors_info[M + vertex_cell_neighbors_info[M]] = i;
    }
  }

  // correction for periodic boundaries
  if(is_periodic)
  {
    for(int i=0 ; i<(int)PeriodicVertexConnect.size() ; i++)
    { 
      int nb_cells = 0;
      std::vector<int> neighbors_cells;

      // collect cells neighboring periodically connected vertices
      for(int j:PeriodicVertexConnect.at(i))
      {
        int idx = j * MaxCpV;
        nb_cells = vertex_cell_neighbors_info.at(idx);

        neighbors_cells.insert(neighbors_cells.end(),
                               &vertex_cell_neighbors_info[idx+1],
                               &vertex_cell_neighbors_info[idx+nb_cells+1]);
      }
      sort_and_select(neighbors_cells);
      nb_cells = (int)neighbors_cells.size();

      // and correct vertex_cell_neighbors_info
      for(int j:PeriodicVertexConnect.at(i))
      {
        int idx = j * MaxCpV;
        vertex_cell_neighbors_info.at(idx) = nb_cells;
        std::copy(neighbors_cells.begin(), neighbors_cells.end(),
                  vertex_cell_neighbors_info.begin()+idx+1);
      }
    }
  }

  return vertex_cell_neighbors_info;
}

void set_all_vertices_to_minus1(const TCollection* coll, int rank,
                                int N_VertInCell)
{
  int N_Cells = coll->GetN_Cells();
  for(int i = 0; i < N_Cells; i++)
  {
    auto cell = coll->GetCell(i);
    if(cell->GetSubDomainNo() == rank)
    {
      for(int j = 0; j < N_VertInCell; j++)
        cell->GetVertex(j)->SetClipBoard(-1);
    }
  }
}

std::vector<int>
get_vertex_neigh_ranks(int N_RootVertices, int MaxCpV, int rank, int size,
                       const TCollection* coll,
                       const std::vector<int>& VertexNumbers,
                       const std::vector<int>& vertex_cell_neighbors_info,
                       int MaxRankPerV)
{
  std::vector<int> VertNeibRanks(N_RootVertices * MaxCpV, 0);

  std::vector<int> HaloCellIndex(size, -1);
  std::vector<int> HaloCellGlobalNo(MaxCpV);
  std::vector<int> HaloCellLocVertNo(MaxCpV);
  int N_Cells = coll->GetN_Cells();
  int N_VertInCell = coll->GetCell(0)->GetN_Vertices();

  for(int i = 0; i < N_Cells; i++)
  {
    auto cell = coll->GetCell(i);
    int ID = cell->GetSubDomainNo();
    /** run only through own cells */
    if(ID == rank)
    {
      cell->SetAsOwnCell();
      /** set SubDomainVert if any vert in this cell is so needed for moving
       * meshes and setting cross vertex */
      for(int j = 0; j < N_VertInCell; j++)
      {
        auto CurrVert = cell->GetVertex(j);
        if(CurrVert->GetClipBoard() != -1)
          continue;

        CurrVert->SetClipBoard(5);
        int N = VertexNumbers[i * N_VertInCell + j] * MaxCpV;
        int N_CellsInThisVert = vertex_cell_neighbors_info.at(N); // nb cells containing this vertex

        // check! any subdomain cell containg this vert has to be added
        // Setting all halo cells
        bool UPDATE = true;
        for(int ii = 1; ii <= N_CellsInThisVert; ii++)
        {
          int VertexCellNo = vertex_cell_neighbors_info.at(N + ii); // loop over cells
          auto Vertexcell = coll->GetCell(VertexCellNo);
          int VertexcellID = Vertexcell->GetSubDomainNo();

          /** Hallo Cell */
          if(VertexcellID != rank)
          {
            // vertex belong to diff subdomain
            if(UPDATE)
            {
              cell->SetAsDependentCell();
              CurrVert->SetAsSubDomainVert();
              HaloCellIndex[VertexcellID] = VertexCellNo;
              UPDATE = false;
            }
            Vertexcell->SetAsHaloCell();

            /** changed on 4 Feb 2012 - by Sashi */
            if(HaloCellIndex[VertexcellID] == -1)
              HaloCellIndex[VertexcellID] = VertexCellNo;
          }
        } //  for(ii=1;ii<

        if(CurrVert->IsSubDomainVert())
        {
          for(int ii = 1; ii <= N_CellsInThisVert; ii++)
          {
            int VertexCellNo = vertex_cell_neighbors_info.at(N + ii);
            auto Vertexcell = coll->GetCell(VertexCellNo);
            int VertexcellID = Vertexcell->GetSubDomainNo();

            // all cells associated with this vert are DepCells
            Vertexcell->SetAsDependentCell();

            if(VertexcellID == rank)
              continue;

            int N_CellsIn_b = VertNeibRanks[N];
            UPDATE = true;
            for(int jj = 1; jj <= N_CellsIn_b; jj++)
            {
              if(VertNeibRanks[N + jj] == VertexcellID)
              {
                UPDATE = false;
                break;
              }
            }

            if(UPDATE)
            {
              int N2 = HaloCellIndex[VertexcellID];
              HaloCellGlobalNo[VertNeibRanks[N]] = N2;

              // find the local index of this vertex in the neib cell
              int jj = 0;
              if(CurrVert->IsPeriodicVert())
              {
                while(CurrVert->GetPeriodicVertIndex()
                  != (coll->GetCell(N2))->GetVertex(jj)->GetPeriodicVertIndex())
                  jj++;
              }
              else
              {
                while(CurrVert != (coll->GetCell(N2))->GetVertex(jj))
                  jj++;
              }
              HaloCellLocVertNo[VertNeibRanks[N]] = jj;

              VertNeibRanks[N]++;
              VertNeibRanks[N + VertNeibRanks[N]] = VertexcellID;

              if(MaxRankPerV < VertNeibRanks[N])
                MaxRankPerV = VertNeibRanks[N];
            } //   if(UPDATE)
          }   //  for(ii=1;ii<=N_CellsInT
        }

        int N_SubDomIn_a = VertNeibRanks[N];
        int* Temp = &VertNeibRanks[N + 1];

        /** set only one cell (lowest index cell) from each subdomain as the
         * Halo cell for this vertex */
        CurrVert->SetSubDomainInfo(N_SubDomIn_a, Temp, HaloCellGlobalNo.data(),
                                   HaloCellLocVertNo.data());
        // reset
        for(int jj = 0; jj < size; jj++)
          HaloCellIndex[jj] = -1;
      }          // for(j=0;j<N_VertInCell
    }            // if(ID==rank
  }              //  for(i=0;i<N_Cells
  MaxRankPerV++; // including own rank
  return VertNeibRanks;
}

void get_n_own_cells(const TCollection* coll, int rank, int& N_LocalCells,
                     int& N_OwnCells)
{
  N_LocalCells = 0;
  N_OwnCells = 0;
  int N_Cells = coll->GetN_Cells();
  for(int i = 0; i < N_Cells; i++)
  {
    auto cell = coll->GetCell(i);
    int ID = cell->GetSubDomainNo();
    if(ID == rank || cell->IsHaloCell())
      N_LocalCells++;
    if(ID == rank)
      N_OwnCells++;
  } //  for (i=0;i<N_Cells;i++)
}

// we need pointers here as domain will take ownership
void get_subdomain_cells_and_global_cell_index(const TCollection* coll,
                                               int*& GlobalCellIndex,
                                               TBaseCell**& SubDomainCells,
                                               int N_LocalCells, int rank,
                                               int N_OwnCells)
{
  GlobalCellIndex = nullptr;
  SubDomainCells = nullptr;
  if(N_LocalCells)
  {
    // collect the own collection
    // we use 'new' here as the domain will take ownership
    SubDomainCells = new TBaseCell*[N_LocalCells];
    GlobalCellIndex = new int[N_LocalCells];
  }

  int N_Cells = coll->GetN_Cells();
  for(int i = 0, N = 0, M = 0; i < N_Cells; i++)
  {
    auto cell = coll->GetCell(i);
    int ID = cell->GetSubDomainNo();
    // in the collection, first own cells then Halo cells
    if(ID == rank)
    {
      SubDomainCells[N] = cell;
      GlobalCellIndex[N] = i;
      N++;
    }
    else if(cell->IsHaloCell())
    {
      SubDomainCells[N_OwnCells + M] = cell;
      GlobalCellIndex[N_OwnCells + M] = i;
      M++;
    }
    // else: there used to be code here which collected cells, edges, and
    // vertices which can be deleted. However these were never deleted, so I
    // (Ulrich) keep it this way
  } //  for (i=0;i<N_Cells;i++)
}

void fill_subdomain_info(TBaseCell** SubDomainCells, int N_OwnCells, int rank)
{
  for(int i = 0; i < N_OwnCells; i++)
  {
    auto cell = SubDomainCells[i];
    if(cell->IsDependentCell())
    {
      int N_Edges = cell->GetN_Edges();
      for(int j = 0; j < N_Edges; j++)
        cell->GetEdge(j)->SetClipBoard(-1);
    }
  }

  for(int i = 0; i < N_OwnCells; i++)
  {
    auto cell = SubDomainCells[i];
    if(cell->IsDependentCell())
    {
      int N_Edges = cell->GetN_Edges();
      for(int j = 0; j < N_Edges; j++)
      {
        TEdge* edge = cell->GetEdge(j);
        if(edge->GetClipBoard() == -1)
        {
          edge->SetClipBoard(5);
          edge->InitSubDomainInfo(rank);
        } // if(edge->GetCli
      }   // for(j=0;j<N_Edges;
    }     // if(cell->IsDependentCell
  }       // for(i=0;i<N_Cells;i++
}

void replace_subdomain_joints(int N_OwnCells, TBaseCell** SubDomainCells,
                              int* GlobalCellIndex, int N_JointsInCell,
                              int rank, const std::vector<int>& VertexNumbers,
                              std::vector<int>& VertNeibRanks, int N_VertInCell,
                              int MaxCpV)
{
  for(int i = 0; i < N_OwnCells; i++)
  {
    auto cell = SubDomainCells[i];
    /** run only through Dependent cells */
    if(cell->IsDependentCell())
    {
      int GblCellNr = GlobalCellIndex[i];
      int ID = cell->GetSubDomainNo();
      /** change the faces, also edges in 3D */
      for(int j = 0; j < N_JointsInCell; j++)
      {
        TJoint* Joint = cell->GetJoint(j);

        if(Joint->GetType() == JointEqN || Joint->GetType() == PeriodicJoint)
        {
          TBaseCell* neib_cell = Joint->GetNeighbour(cell);
          int Neib_ID = neib_cell->GetSubDomainNo();
          if(ID != Neib_ID)
          {
            // this joint belongs to two SubDomains
            cell->SetAsSubDomainInterfaceCell();

            if(Joint->GetType() == JointEqN)
            {
              int GlobCellNo = neib_cell->GetGlobalCellNo();
              // find neib cell local face number
              int m = 0;
              while(Joint != neib_cell->GetJoint(m))
                m++;

              delete Joint;
              TJoint* NewJoint = new TSubDomainJoint(cell, neib_cell,
                                                     Neib_ID, GlobCellNo, m);
              cell->SetJoint(j, NewJoint);
              neib_cell->SetJoint(m, NewJoint);
              NewJoint->SetMapType();
            }

            // set all edges in this face as SubDomainEdges
            const TShapeDesc* ShapeDesc = cell->GetShapeDesc();
            int MaxLen;
            const int *TmpFV, *TmpLen;
            ShapeDesc->GetFaceEdge(TmpFV, TmpLen, MaxLen);
            int N_FaceEdges = TmpLen[j];
            for(int n = 0; n < N_FaceEdges; n++)
            {
              TEdge* edge = cell->GetEdge(TmpFV[j * MaxLen + n]);
              edge->SetAsNotCrossEdgeFor(rank, Neib_ID);
            } // for(n=0;n<N_Edge

            // remove the Neib_ID from VertNeibRanks
            // these face vertices cannot be cross vertices for processors ID
            // and Neib_ID
            ShapeDesc->GetFaceVertex(TmpFV, TmpLen, MaxLen);
            for(int ii = 0; ii < TmpLen[j]; ii++)
            {
              int jj = TmpFV[j * MaxLen + ii];
              int M = VertexNumbers[GblCellNr * N_VertInCell + jj] * MaxCpV;
              int N = VertNeibRanks[M];
              for(int jj = 1; jj <= N; jj++)
              {
                if(VertNeibRanks[M + jj] == Neib_ID)
                {
                  VertNeibRanks[M + jj] = VertNeibRanks[M + N];
                  VertNeibRanks[M]--;
                  break;
                } // for(jj=1;jj<=N;jj++)
              }
            } // for(ii=0;ii<TmpL
          }   // if(ID!=Neib_ID
        }     // if(Joint->GetT
      }       // for (j=0;j< ;
    }         // if(cell->IsDependentCell()
  }           // for(i=0;i<N_OwnCells
}

void change_joint_type(const TCollection* coll, int rank, int N_HaloCells,
                       TBaseCell** SubDomainCells, int N_OwnCells,
                       int N_JointsInCell)
{
  /** set the Halo BD face as BD face: 07 April 2015 - by Sashi  */
  int N_Cells = coll->GetN_Cells();
  for(int i = 0; i < N_Cells; i++)
  {
    auto cell = coll->GetCell(i);
    // cell in this collection
    if(cell->GetSubDomainNo() == rank || cell->IsHaloCell())
    {
      cell->SetClipBoard(i);
    }
    else
    {
      cell->SetClipBoard(-1);
    }
  } //

  for(int i = 0; i < N_HaloCells; i++)
  {
    auto cell = SubDomainCells[N_OwnCells + i];
    for(int j = 0; j < N_JointsInCell; j++)
    {
      auto Joint = cell->GetJoint(j);
      if(Joint->GetType() == JointEqN)
      // if(Joint->GetType() == SubDomainHaloJoint)
      {
        auto neib_cell = Joint->GetNeighbour(cell);
        // j is a Halo BD face
        if(neib_cell->GetClipBoard() == -1)
        {
          Joint->ChangeType(SubDomainHaloJoint); // BoundaryFace);
        }                                        // if(neib_cell->GetCl
      }                                          //   if(Joint->GetTy
    }                                            // for(j=0;j<N_JointsInCel
  }                                              //   for(i=0;i<N_HaloCells;i+
}

void handle_cross_vertices(int N_OwnCells, TBaseCell** SubDomainCells,
                           int* GlobalCellIndex,
                           const std::vector<int>& VertexNumbers,
                           std::vector<int>& VertNeibRanks, int N_VertInCell,
                           int MaxCpV)
{
  /** find cross edges (if any) */
  for(int i = 0; i < N_OwnCells; i++)
  {
    auto cell = SubDomainCells[i];
    if(cell->IsDependentCell())
    {
      int N_Edges = cell->GetN_Edges();
      for(int j = 0; j < N_Edges; j++)
        cell->GetEdge(j)->SetClipBoard(-1);
    } // if(cell->IsDependentCell())
  }   // for(i=0;i<N_OwnCel

  for(int i = 0; i < N_OwnCells; i++)
  {
    auto cell = SubDomainCells[i];
    if(cell->IsDependentCell())
    {
      int M = GlobalCellIndex[i];
      int N_Edges = cell->GetN_Edges();
      const TShapeDesc* ShapeDesc = cell->GetShapeDesc();
      const int* EdgeVertex;
      ShapeDesc->GetEdgeVertex(EdgeVertex);

      for(int j = 0; j < N_Edges; j++)
      {
        TEdge* edge = cell->GetEdge(j);
        if(edge->GetClipBoard() == -1 && edge->IsSubDomainEdge())
        {
          edge->SetClipBoard(5);
          int N_CrossNeibs = edge->GetN_CrossNeibs();
          if(N_CrossNeibs == 0)
            continue;

          cell->SetAsCrossEdgeCell();

          int a = EdgeVertex[2 * j];
          int b = EdgeVertex[2 * j + 1];

          TVertex* Vert_a = cell->GetVertex(a);
          edge->SetCrossNeibInfo(Vert_a);
          int* CrossNeibsRank;
          edge->GetCrossEdgeNeibs(N_CrossNeibs, CrossNeibsRank);

          int M1 = VertexNumbers[M * N_VertInCell + a] * MaxCpV;
          int N1 = VertNeibRanks[M1];
          int M2 = VertexNumbers[M * N_VertInCell + b] * MaxCpV;
          int N2 = VertNeibRanks[M2];

          for(int jj = 0; jj < N_CrossNeibs; jj++)
          {
            int Neib_ID = CrossNeibsRank[jj];
            for(int kk = 1; kk <= N1; kk++)
            {
              if(VertNeibRanks[M1 + kk] == Neib_ID)
              {
                VertNeibRanks[M1 + kk] = VertNeibRanks[M1 + N1];
                VertNeibRanks[M1]--;
                break;
              }
            }

            for(int kk = 1; kk <= N2; kk++)
            {
              if(VertNeibRanks[M2 + kk] == Neib_ID)
              {
                VertNeibRanks[M2 + kk] = VertNeibRanks[M2 + N2];
                VertNeibRanks[M2]--;
                break;
              }
            }
          } // for(jj=0;jj<N_CrossNei
        }   //   if(edge->GetClipBoard()==-1 && edge->IsSubDomainEdge())
      }     // for(j=0;j<N_Edge
    }       // if(cell->IsDependentCell())
  }         // for(i=0;i<N_OwnCel
}

void set_vertex_cells(int N_OwnCells, TBaseCell** SubDomainCells,
                      int N_HaloCells, int* GlobalCellIndex, int N_VertInCell,
                      int MaxCpV, const std::vector<int>& VertexNumbers,
                      const std::vector<int>& vertex_cell_neighbors_info,
                      const TCollection* coll, int rank)
{
  for(int i = 0; i < N_OwnCells; i++)
  {
    auto cell = SubDomainCells[i];
    if(cell->IsDependentCell())
    {
      for(int j = 0; j < N_VertInCell; j++)
        cell->GetVertex(j)->SetClipBoard(-1);
    }
  }

  for(int i = 0; i < N_HaloCells; i++)
  {
    SubDomainCells[N_OwnCells + i]->SetClipBoard(-1);
    // store GlobalCellIndex
    SubDomainCells[N_OwnCells + i]->SetClipBoard_Par(
        GlobalCellIndex[N_OwnCells + i]);
  }

  /** set incident cell list for all vertices */
  std::vector<TBaseCell*> IncCells(MaxCpV);
  std::vector<TBaseCell*> IncHalloCells(MaxCpV);
  for(int i = 0; i < N_OwnCells; i++)
  {
    auto cell = SubDomainCells[i];
    int M = GlobalCellIndex[i];
    if(cell->IsDependentCell())
    {
      for(int j = 0; j < N_VertInCell; j++)
      {
        // set the vertexcells for all vertices in dep. cells
        auto CurrVert = cell->GetVertex(j);
        if(CurrVert->GetClipBoard() != -1)
          continue;

        CurrVert->SetClipBoard(5);

        int N = VertexNumbers[M * N_VertInCell + j] * MaxCpV;
        int N_CellsInThisVert = vertex_cell_neighbors_info.at(N);

        for(int ii = 1; ii <= N_CellsInThisVert; ii++)
          IncCells[ii - 1]
              = coll->GetCell(vertex_cell_neighbors_info.at(N + ii));

        CurrVert->SetVertexCells(N_CellsInThisVert, IncCells.data());


        /** added, halo cells vertCells neib info on 11 Mar 2015 - Sashi */
        for(int ii = 0; ii < N_CellsInThisVert; ii++)
        {
          // fill only for Halo cells, since own Dependent cells' vertices got
          // filled above
          if(IncCells[ii]->IsHaloCell())
          {
            /** this Halo cell is already handled */
            if(IncCells[ii]->GetClipBoard() != -1)
              continue;

            IncCells[ii]->SetClipBoard(5);

            for(int jj = 0; jj < N_VertInCell; jj++)
            {
              auto NeibCurrVert = IncCells[ii]->GetVertex(jj);

              /** even in Halo cells, set neib info only for unset vertices */
              if(NeibCurrVert->GetNNeibs() != 0)
                continue;

              int MM = IncCells[ii]->GetClipBoard_Par();
              int NN = VertexNumbers[MM * N_VertInCell + jj] * MaxCpV;
              int N_CellsInHaloVert = vertex_cell_neighbors_info.at(NN);
              int ll = 0;
              for(int kk = 1; kk <= N_CellsInHaloVert; kk++)
              {
                auto HalloCell
                    = coll->GetCell(vertex_cell_neighbors_info.at(NN + kk));

                /** include only own and Halo in the neib list */
                if(HalloCell->GetSubDomainNo() == rank
                   || HalloCell->IsHaloCell())
                {
                  IncHalloCells[ll] = HalloCell;
                  ll++;
                } // if(HalloCell->GetSubDomainNo(
              }   // for(kk=1;kk<=N_CellsInHaloVert;kk++)
              NeibCurrVert->SetVertexCells(ll, IncHalloCells.data());
            } //   for(jj=0;jj<N_VertInCell;jj++)
          }   // if(IncCells[ii]->IsHalo
        }
      } //  for(jj=0;jj<N_VertInCell;jj++)
    }   //  for (j=0;j<N_VertInCell;
  }     //  for(i=0;i<N_OwnCells
}

void set_cross_cells(int N_OwnCells, TBaseCell** SubDomainCells,
                     int* GlobalCellIndex, int N_VertInCell, int MaxCpV,
                     const std::vector<int>& VertexNumbers,
                     const std::vector<int>& VertNeibRanks)
{
  for(int i = 0; i < N_OwnCells; i++)
  {
    auto cell = SubDomainCells[i];
    if(cell->IsDependentCell())
    {
      for(int j = 0; j < N_VertInCell; j++)
        cell->GetVertex(j)->SetClipBoard(-1);
    }
  }

  for(int i = 0; i < N_OwnCells; i++)
  {
    auto cell = SubDomainCells[i];
    if(cell->IsDependentCell())
    {
      int M = GlobalCellIndex[i];
      for(int j = 0; j < N_VertInCell; j++)
      {
        auto CurrVert = cell->GetVertex(j);
        // continue if the vert is already handled or not a subdomain vert
        if((CurrVert->GetClipBoard() != -1) || !(CurrVert->IsSubDomainVert()))
          continue;

        CurrVert->SetClipBoard(5);

        int M1 = VertexNumbers[M * N_VertInCell + j] * MaxCpV;
        int N1 = VertNeibRanks[M1];
        if(N1 != 0)
          cell->SetAsCrossVertexCell();
        for(int ii = 1; ii <= N1; ii++)
          CurrVert->AddCrossNeib(VertNeibRanks[M1 + ii]);
      } // for(j=0;j<N_VertInCell;j
    }   // if(cell->IsDependentCell())
  }     // for(i=0;i<N_OwnCells
}

void update_domain(TDomain* Domain, int N_LocalCells, int N_OwnCells,
                   TBaseCell** SubDomainCells, int* GlobalCellIndex)
{
  if(N_LocalCells)
  {
    Domain->ReplaceTreeInfo(N_LocalCells, SubDomainCells, GlobalCellIndex,
                            N_OwnCells);
  }
  else
  {
    Domain->SetN_OwnCells(N_OwnCells);
  }

  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
  TDatabase::IteratorDB[It_LE]->SetParam(Domain);
  TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
  TDatabase::IteratorDB[It_Between]->SetParam(Domain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);
}


void prepare_everything(const TCollection* coll, int rank, int size,
                        int N_VertInCell, int N_RootVertices, int MaxCpV,
                        const std::vector<int>& VertexNumbers,
                        const std::vector<int>& vertex_cell_neighbors_info,
                        int MaxRankPerV, int N_JointsInCell, TDomain* Domain,
                        bool flag)
{
  /** *********************************************/
  /** STEP 8.2 : SET ALL VERTICES TO -1 */
  /** *********************************************/
  set_all_vertices_to_minus1(coll, rank, N_VertInCell);

  /** *********************************************/
  /** STEP 9 : Fill the VertNeibRanks info        */
  /** first column contains how many ranks contain this vertex further columns
   * contain the rank ID of the subdomains */
  /** *********************************************/
  std::vector<int> VertNeibRanks = get_vertex_neigh_ranks(
      N_RootVertices, MaxCpV, rank, size, coll, VertexNumbers,
      vertex_cell_neighbors_info, MaxRankPerV);

  /** *************************************************/
  /** STEP 10 : Get own SubDomain collection of cells */
  /** *************************************************/
  int N_LocalCells;
  int N_OwnCells;
  get_n_own_cells(coll, rank, N_LocalCells, N_OwnCells);

  if(N_OwnCells == 0)
    ErrThrow("No mesh cells assigned to this processor.");

  int N_HaloCells = N_LocalCells - N_OwnCells;
  int* GlobalCellIndex = nullptr;
  TBaseCell** SubDomainCells = nullptr;
  get_subdomain_cells_and_global_cell_index(
      coll, GlobalCellIndex, SubDomainCells, N_LocalCells, rank, N_OwnCells);

  // **************************************************************************/
  //** STEP 13 : Edges will be flaged as SubDomainEdges and subdomain infos are
  // filled    */
  //**
  //***************************************************************************/
  fill_subdomain_info(SubDomainCells, N_OwnCells, rank);

  /** change all subdomain joint from JointEqN into TSubDomainJoint */
  replace_subdomain_joints(N_OwnCells, SubDomainCells, GlobalCellIndex,
                           N_JointsInCell, rank, VertexNumbers, VertNeibRanks,
                           N_VertInCell, MaxCpV);

  if(flag)
    change_joint_type(coll, rank, N_HaloCells, SubDomainCells, N_OwnCells,
                      N_JointsInCell);

  /** find cross edges (if any) */
  handle_cross_vertices(N_OwnCells, SubDomainCells, GlobalCellIndex,
                        VertexNumbers, VertNeibRanks, N_VertInCell, MaxCpV);

  /** find cross vertex (if any) */
  // set the clipboard
  set_vertex_cells(N_OwnCells, SubDomainCells, N_HaloCells, GlobalCellIndex,
                   N_VertInCell, MaxCpV, VertexNumbers,
                   vertex_cell_neighbors_info, coll, rank);

  set_cross_cells(N_OwnCells, SubDomainCells, GlobalCellIndex, N_VertInCell,
                  MaxCpV, VertexNumbers, VertNeibRanks);

  update_domain(Domain, N_LocalCells, N_OwnCells, SubDomainCells,
                GlobalCellIndex);
}

#ifdef __2D__


#else // 3D

#define MPI_PARTITION_TAG_N_TOTAL_VERTICES 75
#define MPI_PARTITION_TAG_ALL_VERTICES 80
#define MPI_PARTITION_TAG_MAX_CPV 85
#define MPI_PARTITION_TAG_VERTEX_CELL_NEIGHBORS 90

int Partition_Mesh3D(MPI_Comm comm, TDomain* Domain, int metis_type,
                     int& MaxRankPerV)
{
  MaxRankPerV = -1; // will be set later
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  TCollection* coll = Domain->GetCollection(It_Finest, 0);
  int N_Cells = coll->GetN_Cells();
  // this assumes that all cells have the same geometry:
  TBaseCell* first_cell = coll->GetCell(0);
  int N_VertInCell = first_cell->GetN_Vertices();
  int N_JointsInCell = first_cell->GetN_Joints();
  int N_AllLocVert = N_VertInCell * N_Cells;
  std::vector<int> VertexNumbers(N_AllLocVert);
  std::vector<idx_t> Cell_Rank(N_Cells);

  bool is_periodic = Domain->is_periodic();
  const std::vector<std::vector<int>> PeriodicVertexConnect =
                                            Domain->get_PeriodicVertexConnect();

  // if(rank==0)
  //  printf("Number of  ranks: %d\n",  size);

  int global_flag, flag = 0;

  if (N_Cells < size)
  {
    flag = 1;
  }

  MPI_Allreduce(&flag, &global_flag, 1, MPI_INT, MPI_SUM, comm);
  if (global_flag != 0)
  {
    ErrThrow("More ranks than cells!");
  }

  if (size == 1)
  {
    ErrThrow("Total number of process should be greater than 1 (or 2 if root "
             "is not involved in computation), but you have ", size);
  }

  // check, all cell have same number of vertices
  for (int i = 1; i < N_Cells; i++)
  {
    if (coll->GetCell(i)->GetN_Vertices() != N_VertInCell
       || coll->GetCell(i)->GetN_Joints() != N_JointsInCell)
    {
      ErrThrow("Mesh partition for heterogeneous cells are not yet "
               "implemented ",
               N_JointsInCell, " ", N_VertInCell);
    }
  }

  std::vector<int> vertex_cell_neighbors_info{};
  int N_RootVertices;
  int MaxCpV = 0;

  // partition the mesh in the root and send info to all processors
  if (rank == 0) // Code run by root
  {
    if (N_VertInCell != 4 && N_VertInCell != 8)
    {
      ErrThrow("only Tetra or Hexa mesh can be partitioned !!");
    }

    /** *********************************************/
    /** STEP 1 : STORE VERTEX POINTERS CELL-WISE */
    /** *********************************************/
    std::vector<const TVertex*> Vertices = get_sorted_vertices(coll);

    /** ***************************************************/
    /**STEP 3: STORE THE SORTED POINTER ARRAY AS INDICES */
    /** ***************************************************/
    N_RootVertices = get_number_vertex(Vertices, coll, VertexNumbers);

    /** 1. SEND INFO TO ALL PROCESSORS */
    for (int i = 1; i < size; i++)
    {
      MPI_Send(&N_RootVertices, 1, MPI_INT, i, MPI_PARTITION_TAG_N_TOTAL_VERTICES, comm);
    }

    /** 2. SEND INFO TO ALL PROCESSORS */
    for (int i = 1; i < size; i++)
    {
      MPI_Send(VertexNumbers.data(), N_AllLocVert, MPI_INT, i, MPI_PARTITION_TAG_ALL_VERTICES, comm);
    }

    /** *********************************************/
    /** STEP 5 : MESH PARTITION */
    /** *********************************************/
    bool read_metis = Domain->get_database()["read_metis"];
    bool write_metis = Domain->get_database()["write_metis"];

    if (read_metis)
    {
      // read the partition from a file

      MeshPartitionInOut::read_file(*Domain, size, Cell_Rank);
    }
    else
    {
      // we have to compute our own partition - call METIS if we have it

#ifdef PARMOON_WITH_METIS
      // set up data to pass into METIS

      idx_t ne = (idx_t)N_Cells;
      idx_t nn = (idx_t)N_RootVertices;
      idx_t nparts = (idx_t)size;
      idx_t edgecut = 0;
      idx_t ncommon = 3;

      // each cell's index into the vertex numbers array
      std::vector<idx_t> eptr(N_Cells + 1);
      eptr[0] = 0;

      for (int i = 1; i <= N_Cells; ++i)
      {
        eptr[i] = eptr[i - 1] + N_VertInCell;
      }

      // copy the vertex numbers array to a new vector
      // ('int' might be different from 'idx_t')
      std::vector<idx_t> MetisVertexNumbers(N_AllLocVert);
      std::copy(VertexNumbers.begin(), VertexNumbers.end(),
                MetisVertexNumbers.begin());

      // configure METIS
      idx_t options[METIS_NOPTIONS];
      METIS_SetDefaultOptions(options);
      options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
      options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
      options[METIS_OPTION_NUMBERING] = 0;
      options[METIS_OPTION_CONTIG] = 1;

      // will receive vertex ranks - not used below
      std::vector<idx_t> Vert_Rank(N_AllLocVert);
      double t1 = MPI_Wtime();

      if(metis_type == 0)
      {
        METIS_PartMeshNodal(&ne, &nn, eptr.data(), MetisVertexNumbers.data(),
                            nullptr, nullptr, &nparts, nullptr, options,
                            &edgecut, Cell_Rank.data(), Vert_Rank.data());
      }
      else if(metis_type == 1)
      {
        METIS_PartMeshDual(&ne, &nn, eptr.data(), MetisVertexNumbers.data(),
                           nullptr, nullptr, &ncommon, &nparts, nullptr,
                           options, &edgecut, Cell_Rank.data(),
                           Vert_Rank.data());
      }
      else
      {
        ErrThrow("Partition_Mesh3D implemented for metis_type = 0 or 1 !!");
      }
      t1 = MPI_Wtime() - t1;
      Output::root_info("Domain Decomposition",
                        "Time taken for METIS mesh partitioning ", t1, " sec");
#else
      ErrThrow("ParMooN has been compiled without metis, therefore it can not "
               "be called. You have to either read a partitioning or compile "
               "with metis.");
#endif
    }

    if (write_metis)
    {
      // write the partition to a file

      MeshPartitionInOut::write_file(*Domain, size, Cell_Rank);
    }

    /** *********************************************/
    /** STEP 6 : FIND AND COMMUNICATE MAXIMUM CELLS PER VERTEX */
    /** *********************************************/
    MaxCpV = get_max_n_cells_per_vertex(VertexNumbers, N_RootVertices,
                                        is_periodic, PeriodicVertexConnect);

    // count it up by 1, because the array vertex_cell_neighbors_info needs room
    // for the number of cells sharing each vertex
    MaxCpV++;

    // broadcast MaxCpV to all processes
    for (int i = 1; i < size; i++)
    {
      MPI_Send(&MaxCpV, 1, MPI_INT, i, MPI_PARTITION_TAG_MAX_CPV, comm);
    }

    /** ***********************************************************************/
    /** STEP 7 : CREATE AN ARRAY CONTAINING INFO REGARDING THE NEIGHBOURS OF  */
    /**          VERTICES (CELL INDEX)                                        */
    /** vertex_cell_neighbors_info's first column contains number of neib     */
    /** cells associated with each vertex,                                    */
    /** further columns contain the cell numbers associated with this vertex  */
    /** ***********************************************************************/
    // This vector gets filled with the information, which vertex is shared by
    // which cells.
    vertex_cell_neighbors_info = get_neighbors_of_vertices(N_RootVertices, MaxCpV,
                                                           N_Cells, N_VertInCell,
                                                           VertexNumbers,
                                                           is_periodic,
                                                           PeriodicVertexConnect);
    for (int i = 1; i < size; i++)
    {
      MPI_Send(vertex_cell_neighbors_info.data(), N_RootVertices * MaxCpV,
               MPI_INT, i, MPI_PARTITION_TAG_VERTEX_CELL_NEIGHBORS, comm);
    }
  } // end if root
  else
  {
    // code run by worker processes - receive information form root

    MPI_Status status;
    MPI_Recv(&N_RootVertices, 1, MPI_INT, 0, MPI_PARTITION_TAG_N_TOTAL_VERTICES, comm, &status);
    MPI_Recv(VertexNumbers.data(), N_AllLocVert, MPI_INT, 0, MPI_PARTITION_TAG_ALL_VERTICES, comm, &status);
    MPI_Recv(&MaxCpV, 1, MPI_INT, 0, MPI_PARTITION_TAG_MAX_CPV, comm, &status);

    vertex_cell_neighbors_info.resize(N_RootVertices * MaxCpV, 0);
    MPI_Recv(vertex_cell_neighbors_info.data(), N_RootVertices * MaxCpV,
             MPI_INT, 0, MPI_PARTITION_TAG_VERTEX_CELL_NEIGHBORS, comm, &status);
  } // else if(rank==0)

  if (Domain->get_database().try_get_value("metis_optimize_distributed", false))
  {
    // Try to optimize the distribution of cells across different physical
    // devices to minimize network overhead. Note that this happens *after*
    // (possibly) writing the METIS file because ranks might match to devices
    // differently despite an identical rank count

    std::vector<char> processor_name(MPI_MAX_PROCESSOR_NAME);
    int processor_name_length;

    MPI_Get_processor_name(processor_name.data(), &processor_name_length);

    if (rank == 0)
    {
      // optimization (or not) happens on root.

      // first, gather each rank's processor name

      std::vector<int> processor_name_lengths(size);

      MPI_Gather(&processor_name_length, 1, MPI_INT,
        processor_name_lengths.data(), 1, MPI_INT,
        0, comm);

      std::vector<int> processor_name_offsets(size);

      int offset = 0;
      for (int i = 0; i < size; i++)
      {
        processor_name_offsets[i] = offset;
        offset += processor_name_lengths[i];
      }

      std::vector<char> processor_names(offset);

      MPI_Gatherv(processor_name.data(), processor_name_length, MPI_CHAR,
        processor_names.data(), processor_name_lengths.data(), processor_name_offsets.data(), MPI_CHAR,
        0, comm);

      // count distinct processor ('device') names and track associated ranks

      std::vector<std::string> device_names;
      std::vector<int> device_indices(size);

      for (int i = 0; i < size; i++)
      {
        std::string name_str(processor_names.data() + processor_name_offsets[i], processor_name_lengths[i]);

        int k = -1;
        for (unsigned int j = 0; j < device_names.size(); j++)
        {
          if (device_names[j] == name_str)
          {
            k = j;
            break;
          }
        }

        if (k >= 0)
        {
          device_indices[i] = k;
        }
        else
        {
          device_indices[i] = device_names.size();
          device_names.push_back(name_str);
        }

        // Output::print("Rank ", i, " is running on device ", device_indices[i], ", name '", device_names[device_indices[i]], "'");
      }

      unsigned int num_devices = device_names.size();

      // Output::print("Total: ", num_devices, " distinct devices.");

      if (num_devices < 2)
      {
        Output::print("All ranks are running on the same physical device; "
          "no need to optimize communication.");
      }
      else if (num_devices >= (unsigned)size)
      {
        Output::print("All ranks are running on distinct devices; "
          "no need to optimize communication.");
      }
      else
      {
        // optimize domain decomposition

        // count ranks per device

        std::vector<int> ranks_per_device(num_devices, 0);

        for (int i = 0 ; i < size; i++)
        {
          ranks_per_device[device_indices[i]]++;
        }

        // compute edge weights for each pair of ranks:
        //    halo cells (adjacent cells belonging to the other rank) plus
        //    'reverse halo' cells (own cells adjacent to the other rank)
        // this results in a symmetric weight matrix.

        std::vector<int> edge_weights(size * size, 0);
        std::unordered_set<int> neighboring_cells;

        for (int k = 0; k < N_Cells; k++)
        {
          neighboring_cells.clear();

          int rank_k = Cell_Rank[k];

          for (int i = 0; i < N_VertInCell; i++)
          {
            int M = VertexNumbers[k * N_VertInCell + i] * MaxCpV;
            int num_neighbors = vertex_cell_neighbors_info[M];

            for (int j = 1; j <= num_neighbors; j++)
            {
              int r = vertex_cell_neighbors_info[M + j];

              if (r != k && Cell_Rank[r] != rank_k)
              {
                neighboring_cells.emplace(r);
              }
            }
          }

          for (int r: neighboring_cells)
          {
            int rank_r = Cell_Rank[r];
            edge_weights[rank_k * size + rank_r]++;
            edge_weights[rank_r * size + rank_k]++;
          }
        }

        // for each rank, find adjacent ranks (nonzero edge weight)

        std::vector<std::vector<int>> neighboring_ranks(size);

        for (int i = 0; i < size; i++)
        {
          neighboring_ranks[i] = std::vector<int>();

          /*std::ostringstream os;
          os << "Rank " << i << " neighbors:";

          for (int j = 0; j < size; j++)
          {
            int w = edge_weights[i * size + j];

            if (w > 0)
            {
              neighboring_ranks[i].push_back(j);
              os << "\n\tRank " << j << ", " << w << " connections";
            }
          }

          Output::print(os.str());*/
        }

        // TODO: look for clever reformulations as e.g. an ILP
        //       OR call METIS twice above
        // for now: greedy method

        // each device's list of ranks
        std::vector<std::vector<int>> ranks_assigned(num_devices);

        // ranks not yet assigned to a device
        std::set<int> available_ranks;

        // temp set of ranks assigned to each device
        std::set<int> my_ranks;

        // all ranks are available initially
        for (int i = 0; i < size; i++)
        {
          available_ranks.emplace(i);
        }

        for (unsigned int m = 0; m < num_devices - 1; m++)
        {
          // each device starts with no ranks
          my_ranks.clear();

          int ranks_to_assign = ranks_per_device[m];
          int current_weight = 0;

          while (ranks_to_assign > 0)
          {
            // find the available rank which would result in the best
            // change of weight

            int best_rank = -1;
            int best_weight = 0;

            for (int i: available_ranks)
            {
              int new_weight = current_weight;

              for (int j: neighboring_ranks[i])
              {
                // check all the neighboring interfaces

                if (my_ranks.count(j) == 0)
                {
                  // this would be a new interface
                  new_weight += edge_weights[i * size + j];
                }
                else
                {
                  // this existing interface would become internal
                  new_weight -= edge_weights[i * size + j];
                }
              }

              if (best_rank < 0 || new_weight < best_weight)
              {
                best_weight = new_weight;
                best_rank = i;
              }
            }

            if (best_rank >= 0)
            {
              // assign to this device, update weight, mark unavailable

              --ranks_to_assign;
              current_weight = best_weight;
              my_ranks.emplace(best_rank);
              available_ranks.erase(best_rank);
            }
            else
            {
              ErrThrow("Failed to find assignable rank!");
            }
          }

          // using a set above ensures that this is sorted
          ranks_assigned[m] = std::vector<int>(my_ranks.begin(), my_ranks.end());

          /*std::ostringstream os;

          os << "Ranks to place on device '" << device_names[m] << "':";

          for (int i: ranks_assigned[m])
          {
            os << " " << i;
          }

          os << "\n\tTotal interface weight: " << current_weight;

          Output::print(os.str());*/
        }

        // last device gets all the remaining ranks
        ranks_assigned[num_devices - 1] = std::vector<int>(available_ranks.begin(), available_ranks.end());

        // figure out new rank assignment based on device assignment
        std::vector<int> rank_permutation(size);
        std::vector<int> counter_per_device(num_devices, 0);

        for (int i = 0; i < size; i++)
        {
          int m = device_indices[i];

          rank_permutation[i] = ranks_assigned[m][counter_per_device[m]++];

          // Output::print("Mapping rank ", i, " to rank ", rank_permutation[i]);
        }

        for (int k = 0; k < N_Cells; k++)
        {
          Cell_Rank[k] = rank_permutation[Cell_Rank[k]];
        }
      }
    }
    else
    {
      // non-root ranks only contribute to process name gathering

      MPI_Gather(&processor_name_length, 1, MPI_INT,
        nullptr, 1, MPI_INT,
        0, comm);

      MPI_Gatherv(processor_name.data(), processor_name_length, MPI_CHAR,
        nullptr, nullptr, nullptr, MPI_CHAR,
        0, comm);
    }
  }

  // Send Metis output Cell_Rank to all processes.

  if (std::is_same<idx_t, int32_t>::value)
  {
    MPI_Bcast(Cell_Rank.data(), N_Cells, MPI_INT32_T, 0, comm);
  }
  else if (std::is_same<idx_t, int64_t>::value)
  {
    MPI_Bcast(Cell_Rank.data(), N_Cells, MPI_INT64_T, 0, comm);
  }
  else
  {
    ErrThrow("metis index type not known."); // copy might be possible here
  }

  /** ********************************************************/
  /** STEP 8.1 : SETTING SUBDOMAIN NUMBER AND GLOBAL CELL NO */
  /** ********************************************************/
  for(int i = 0; i < N_Cells; i++)
  {
    auto cell = coll->GetCell(i);
    cell->SetSubDomainNo(Cell_Rank[i]);
    cell->SetGlobalCellNo(i);
  }

  prepare_everything(coll, rank, size, N_VertInCell, N_RootVertices, MaxCpV,
                     VertexNumbers, vertex_cell_neighbors_info, MaxRankPerV,
                     N_JointsInCell, Domain, false);

  delete coll;
  // Barrier is needed, before calling FECommunicator, since neib process must
  // have all info
  MPI_Barrier(comm);
  return 0; // partitioning successful
}
#endif // 3D


void Domain_Crop(MPI_Comm comm, TDomain* Domain)
{
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  TCollection* coll = Domain->GetCollection(It_Finest, 0);
  int N_Cells = coll->GetN_Cells();
  TBaseCell* cell = coll->GetCell(0);
  int N_VertInCell = cell->GetN_Vertices();
  int N_JointsInCell = cell->GetN_Joints();
  int N_AllLocVert = N_VertInCell * N_Cells;

  bool is_periodic = Domain->is_periodic();
  const std::vector<std::vector<int>> PeriodicVertexConnect =
                                            Domain->get_PeriodicVertexConnect();

  /** STEP 1 : STORE VERTEX POINTERS CELL-WISE */
  std::vector<const TVertex*> Vertices = get_sorted_vertices(coll);

  /** STEP 3: STORE THE SORTED POINTER ARRAY AS INDICES */
  std::vector<int> VertexNumbers(N_AllLocVert);
  int N_RootVertices = get_number_vertex(Vertices, coll, VertexNumbers);

  /** STEP 6 : TO FIND MAXIMUM CELLS PER VERTEX */
  int MaxCpV = get_max_n_cells_per_vertex(VertexNumbers, N_RootVertices,
                                          is_periodic, PeriodicVertexConnect);
  int MaxRankPerV = MaxCpV;
  MaxCpV++; // ACCOUNTING FOR AN EXTRA COLUMN GIVING INFO OF NUMBER OF CELLS
            // SHARING THE VERTEX

  /* STEP 7 : CREATE AN ARRAY CONTAINING INFO REGARDING THE NEIGHBOURS OF */
  /*          VERTICES (CELL INDEX)                                       */
  std::vector<int> PointNeighb = get_neighbors_of_vertices(N_RootVertices, MaxCpV,
                                                           N_Cells, N_VertInCell,
                                                           VertexNumbers,
                                                           is_periodic,
                                                           PeriodicVertexConnect);
  /********** VARIABLES FOR MULTIGRID ****************/
  int Nchildren = static_cast<const TBaseCell*>(coll->GetCell(0))
                      ->GetParent()
                      ->GetN_Children();

  /** copy the parent MPI_ID and global cell_no to children */
  for(int i = 0; i < N_Cells; i++)
  {
    auto cell = coll->GetCell(i);
    auto parent_cell = static_cast<const TBaseCell*>(cell)->GetParent();
    cell->SetSubDomainNo(parent_cell->GetSubDomainNo());
    int childn = parent_cell->GetChildNumber(cell);
    int parentglobalno = parent_cell->GetGlobalCellNo();
    cell->SetGlobalCellNo(parentglobalno * Nchildren + childn);
    cell->SetLocalCellNo(i);
  }

  // note: the last argument here was set to true, but with false all the tests
  // pass as well. I don't know what this really does or where it is useful.
  prepare_everything(coll, rank, size, N_VertInCell, N_RootVertices, MaxCpV,
                     VertexNumbers, PointNeighb, MaxRankPerV, N_JointsInCell,
                     Domain, true);

  // Barrier is needed, before calling FECommunicator, since neib process must
  // have all info
  MPI_Barrier(comm);
}

#endif
