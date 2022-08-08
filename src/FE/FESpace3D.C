// =======================================================================
// %W% %G%
// 
// Class:       TFESpace3D
// Purpose:     class for all 3D finite element spaces
//
// Author:      Gunar Matthies (22.11.97)
//
// History:     start of implementation 22.11.97 (Gunar Matthies)
//
// =======================================================================

#include <FESpace3D.h>
#include <Joint.h>
#include <BoundFace.h>
#include "FEDatabase.h"
#include "FEMapperDatabase.h"
#include "FEMapper.h"
#include "FE3DMapper1Reg.h"

#include <Database.h>
#include <IsoInterfaceJoint3D.h>
#include <TetraAffin.h>
#include <TetraIsoparametric.h>
#include <HexaAffin.h>
#include <HexaTrilinear.h>
#include <HexaIsoparametric.h>

#include <HangingNode.h>
#include "HNDesc.h"

#include <MooNMD_Io.h>
#include <numeric> // accumulate
#include <algorithm>

#include <Edge.h>

#ifdef _MPI
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#endif


TFESpace3D::TFESpace3D(const TCollection *coll, const std::string& name,
                       BoundCondFunct3D *BoundaryCondition, int ord)
: TFESpace(coll, name, ord, 3), boundCondition_{BoundaryCondition}
{
  ConstructSpace();
}

TFESpace3D::TFESpace3D(const TCollection *coll, const std::string& name,
                       BoundCondFunct3D *BoundaryCondition, FE_type *fes)
 : TFESpace(coll, name, fes), boundCondition_{BoundaryCondition}
{
  ConstructSpace();
}

TFESpace3D::TFESpace3D(const TCollection *coll, const std::string& name,
                       BoundCondFunct3D *BoundaryCondition, SpaceType type,
                       int ord)
: TFESpace(coll, name, type, ord, 3), boundCondition_{BoundaryCondition}
{
  ConstructSpace();
}

std::array<int, TFESpace::N_DiffBoundNodeTypes>
TFESpace3D::find_boundary_upper_bounds() const
{
  std::array<int, TFESpace::N_DiffBoundNodeTypes> BoundaryUpperBound{};
  int N_Cells = this->Collection->GetN_Cells();
  for(int i=0;i<N_Cells;i++)
  {
    auto cell = Collection->GetCell(i);
    // find number of boundary nodes in this element
    int n_dof_per_joint = get_fe(i).GetFEDesc()->GetN_JointDOF();

    int n_joints = cell->GetN_Joints();
    for(int j=0;j<n_joints;j++) // loop over all Faces of cell
    {
      auto joint=cell->GetJoint(j);
      if(joint->GetType() == BoundaryFace || joint->GetType() == IsoBoundFace)
      {
        // boundary joint
        auto BoundFace = (TBoundFace *)joint;
        if(BoundFace->GetNeighbour(0) == nullptr)
          BoundFace->SetNeighbour(cell);
        BoundFace->set_index_in_neighbour(cell, j);
        BoundCond Cond0 = get_boundary_condition(*BoundFace);
        BoundaryUpperBound[Cond0] += n_dof_per_joint;
#ifdef _MPI
        /** edges on this face are already set, so no need of checking in BD 
         * edges on this face */
        const int *TmpFE, *ETmpLen;
        int EMaxLen;
        cell->GetShapeDesc()->GetFaceEdge(TmpFE, ETmpLen, EMaxLen);
        int N_FaceEdges = ETmpLen[j];
        for(int n=0;n<N_FaceEdges;n++)
        {
          auto edge = cell->GetEdge(TmpFE[j*EMaxLen+n]);
          edge->SetClipBoard(i);
        }

        /** vertices on this face are already set, so no need of checking in BD
         * vert on this face, 30.06.12, sashi */
        const int *TmpFV, *TmpLen;
        int MaxLen;
        cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
        int N_Points = TmpLen[j];
        for(int n=0;n<N_Points;n++)
          cell->GetVertex(TmpFV[j*MaxLen+n])->SetClipBoard(i);
#endif
      } // endif
    } // endfor j

#ifdef _MPI
    int N_Edges=cell->GetN_Edges();
    int N_VertInCell = cell->GetN_Vertices();
    /** identify BD edges but not part of BD faces */
      if(int n_dof_per_edge = get_fe(i).GetFEDesc()->GetN_EdgeDOF())
    {
      const int *EdgeVertex;
      cell->GetShapeDesc()->GetEdgeVertex(EdgeVertex);

      for(int j=0;j<N_Edges;j++)
      {
        auto edge = cell->GetEdge(j);  
        
        if(edge->GetType()==BDEdge3D || edge->GetType()==IsoEdge3D)
        {
          if(edge->GetClipBoard()!=-1)
            continue;
         
          edge->SetClipBoard(1);
          int comp = edge->get_physical_id();
         
          double xp[2], yp[2], zp[2];
          cell->GetVertex(EdgeVertex[2*j])->GetCoords(xp[0], yp[0], zp[0]);
          cell->GetVertex(EdgeVertex[2*j+1])->GetCoords(xp[1], yp[1], zp[1]);
          double X = (xp[0] + xp[1])/2.;
          double Y = (yp[0] + yp[1])/2.;
          double Z = (zp[0] + zp[1])/2.;   

          BoundCond Cond0;
          boundCondition_(comp, X, Y, Z, Cond0);
          BoundaryUpperBound[Cond0] += n_dof_per_edge;
          /** vertices on this edge are already set, so no need of checking in BD vert on this edge, 30.06.12, sashi */             
          cell->GetVertex(EdgeVertex[2*j])->SetClipBoard(i);
          cell->GetVertex(EdgeVertex[2*j+1])->SetClipBoard(i); 
        } //  if(edge->GetType()==BDEdge3D || edge->GetType()==IsoEdge3D)    
      }//  for(j=0;j<N_Edges;j  
    }

    /** identify BD vert but not part of BD faces/edge */
    /** I hope, it is enough to check BD vertices, sice other vertices are anyway inner, so no prob with Dirichlet BD */
    if(int n_dof_per_vertex = get_fe(i).GetFEDesc()->GetN_VertDOF())
    {
      for(int j=0;j<N_VertInCell;j++)
      {
        auto Vert = cell->GetVertex(j);
        if(Vert->IsBoundVert())
        {
          if( Vert->GetClipBoard()!=-1 )
          {
            // dofs connected to this vertex have already
            // been treated elsewhere - Clemens Bartsch
            continue;
          }

          // BD vert is not yet set 
          Vert->SetClipBoard(i);
          int comp = Vert->get_physical_id();

          double X, Y, Z;
          Vert->GetCoords(X, Y, Z);
          BoundCond Cond0;
          boundCondition_(comp, X, Y, Z, Cond0);
          BoundaryUpperBound[Cond0] += n_dof_per_vertex;
        }
      } //   for(j=0;j<N_VertInCell;j++)
    } //     if( FEDesc0_Obj->GetN_VertDOF() > 0) 
#endif   
  } // endfor i
  return BoundaryUpperBound;
}

void TFESpace3D::ConstructSpace()
{
#ifdef _MPI
  int rank;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);
#endif

  std::vector<THangingNode *> VHN{};
  std::vector<int> HNNumbers{};
  FEMapperDatabase mapper_db;

  auto BoundaryUpperBound = find_boundary_upper_bounds();

  // The 'Counters' count how many dofs of a certain kind have been handled
  // during the loop over all cells which follows. The first entries correspond
  // to the boundary dofs in the order of the enumeration `BoundCond` in 
  // Constants.h. The last entry corresponds to all inner dofs.
  std::array<int, TFESpace::N_DiffBoundNodeTypes+1> Counters{};

  const int some_special_number = FIRSTMARK + 1;
  Counters[0] = some_special_number;
  std::transform(BoundaryUpperBound.begin(), BoundaryUpperBound.end(),
                 Counters.begin(), Counters.begin()+1,
                 [](int boundary_upper_bound, int previous_bound_counter)
                 { return previous_bound_counter - boundary_upper_bound; });
  // the following `Counter` is initialized as:
  //   some_special_number - sum of all BoundaryUpperBound
  int & Counter = Counters[TFESpace::N_DiffBoundNodeTypes];

  int N_Cells = this->Collection->GetN_Cells();

  for(int i=0;i<N_Cells;i++)
  {
    auto cell = Collection->GetCell(i);
    int N_Faces=cell->GetN_Joints();

    auto FE0 = this->get_fe(i);
    auto FEDesc0 = FE0.GetFEDesc_ID();
    auto FEDesc0_Obj = FE0.GetFEDesc();
    int I_K0 = BeginIndex[i];

    const int* TmpFV;
    const int* TmpLen;
    int MaxLen;
    cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
    for(int j=0;j<N_Faces;j++)
    {
      auto joint = cell->GetJoint(j);
      int * Indices0 = FEDesc0_Obj->GetJointDOF(j);

      if(joint->GetType() == BoundaryFace ||
         joint->GetType() == IsoBoundFace)
      {
        // boundary joint
        // there is no need to identify dofs from different cells here. Dofs on
        // this boundary joint which are also in another cell must be on another
        // (inner) joint as well, where the identification is handled.
        auto BoundFace=(TBoundFace *)joint;
        BoundCond Cond0 = get_boundary_condition(*BoundFace);
        
        auto mapper = mapper_db.get_fe_mapper(FEDesc0, FEDesc0);

        mapper->map_boundary(&GlobalNumbers[0], I_K0, Indices0,
                             Counters[Cond0], HNNumbers);
#ifdef _MPI
        /** edges on this face are already set, so no need to check in BD edges on this face */
        const int* TmpFE;
        const int* ETmpLen;
        int EMaxLen;
        cell->GetShapeDesc()->GetFaceEdge(TmpFE, ETmpLen, EMaxLen);
        int N_FaceEdges = ETmpLen[j];
        for(int n=0;n<N_FaceEdges;n++)
         {
          auto edge = cell->GetEdge(TmpFE[j*EMaxLen+n]);
          edge->SetClipBoard(i);
         } // fo 
         
         /** vertices on this face are already set, so no need to check in BD vert on this face, 30.06.12, sashi */        
         int N_Points = TmpLen[j];
        for(int n=0;n<N_Points;n++)
          cell->GetVertex(TmpFV[j*MaxLen+n])->SetClipBoard(i);
#endif         
      } // boundary joint
      else
      {
        // no boundary joint
        const auto* neigh = joint->GetNeighbour(cell);
        if (!neigh)
        {
          // there is no neighbour
          // => finer cell in 1 regular grid
          //    will be handle from coarser cell
        } // !neigh
        else
        {
          // there is a neighbour on same level
          // find local edge of neigh on which cell is
          int joint_index_in_neigh = 0;
          while(neigh->GetJoint(joint_index_in_neigh) != joint)
            joint_index_in_neigh++;
          int n = neigh->GetClipBoard();
          if( n == -1 )
          {
            // check for children of neigh on its face joint_index_in_neigh
            auto refdesc=neigh->GetRefDesc();
            if(refdesc->GetType() != NoRef)
            {
              const int *TmpoFnF;
              const int *TmpLen1;
              int MaxLen1;
              const int *TmpFC;
              const int *TmpLen2;
              int MaxLen2;
              const int *TmpoFnlF;
              const int *TmpCTI;
              refdesc->GetOldFaceNewFace(TmpoFnF, TmpLen1, MaxLen1);
              refdesc->GetFaceChild(TmpFC, TmpLen2, MaxLen2);
              refdesc->GetOldFaceNewLocFace(TmpoFnlF);
              refdesc->GetChildTwistIndex(TmpCTI);
              int NFaces=refdesc->GetShapeDesc()->GetN_Joints();
  
              // Output::print("NFaces: " NFaces);
              // Output::print("TmpLen1[joint_index_in_neigh]: ", TmpLen1[joint_index_in_neigh]);
              // Output::print("MaxLen1: ", MaxLen1);
  
              if(TmpLen1[joint_index_in_neigh] == 4)
              {
                int f1=TmpoFnF[joint_index_in_neigh*MaxLen1+0];
                // Output::print("face1: ", f1);
                int chnum1=TmpFC[f1*MaxLen2];
                auto child1=neigh->GetChild(chnum1);
                int c1=child1->GetClipBoard();
                auto FE0 = this->get_fe(c1);
                auto FEDesc1_Obj = FE0.GetFEDesc();
                auto FEDesc1 = FE0.GetFEDesc_ID();
                // Output::print("child: ", 1, " ", c1);
                int I_K1 = BeginIndex[c1];
                int m=TmpoFnlF[chnum1*NFaces+joint_index_in_neigh];
                int * Indices1 = FEDesc1_Obj->GetJointDOF(m);
                int Twist1 = TmpCTI[f1];
                // Output::print("Twist1: ", Twist1);
  
                int f2=TmpoFnF[joint_index_in_neigh*MaxLen1+1];
                // Output::print("face2: ", f2);
                int chnum2=TmpFC[f2*MaxLen2];
                auto child2=neigh->GetChild(chnum2);
                int c2=child2->GetClipBoard();
                // Output::print("child: ", 2, " ", c2);
                int I_K2 = BeginIndex[c2];
                m=TmpoFnlF[chnum2*NFaces+joint_index_in_neigh];
                int * Indices2 = FEDesc1_Obj->GetJointDOF(m);
                int Twist2 = TmpCTI[f2];
                // Output::print("Twist2: ", Twist2);
  
                int f3=TmpoFnF[joint_index_in_neigh*MaxLen1+2];
                // Output::print("face3: ", f3);
                int chnum3=TmpFC[f3*MaxLen2];
                auto child3=neigh->GetChild(chnum3);
                int c3=child3->GetClipBoard();
                // Output::print("child: ", 3, " ", c3);
                int I_K3 = BeginIndex[c3];
                m=TmpoFnlF[chnum3*NFaces+joint_index_in_neigh];
                int * Indices3 = FEDesc1_Obj->GetJointDOF(m);
                int Twist3 = TmpCTI[f3];
                // Output::print("Twist3: ", Twist3);
  
                int f4=TmpoFnF[joint_index_in_neigh*MaxLen1+3];
                // Output::print("face4: ", f4);
                int chnum4=TmpFC[f4*MaxLen2];
                auto child4=neigh->GetChild(chnum4);
                int c4=child4->GetClipBoard();
                // Output::print("child: ", 4, " ", c4);
                int I_K4 = BeginIndex[c4];
                m=TmpoFnlF[chnum4*NFaces+joint_index_in_neigh];
                int * Indices4 = FEDesc1_Obj->GetJointDOF(m);
                int Twist4 = TmpCTI[f4];
                // Output::print("Twist4: ", Twist4);
  
                // Output::print("big: ", i);
                // Output::print("children: ", c1, " ", c2, " ", c3);
                // Output::print(" ", c4);
  
                auto mapper1reg = mapper_db.get_fe_mapper_one_regular3d(FEDesc0,
                                                                        FEDesc1);
                mapper1reg->Map(&GlobalNumbers[0], I_K1, I_K2, I_K3, I_K4,
                                I_K0, Indices1, Indices2, Indices3,
                                Indices4, Indices0,
                                Twist1, Twist2, Twist3, Twist4,
                                joint->GetMapType(),
                                Counters[DIRICHLET],
                                Counter, VHN, HNNumbers, HN_descriptors);
              }
              // There is exactly one children on the neighbour sharing the current face
              else if(TmpLen1[joint_index_in_neigh] == 1)
              {
                // Get Child of Neighbour that contains this face
                int f1 = TmpoFnF[joint_index_in_neigh * MaxLen1];
                int chnum1 = TmpFC[f1 * MaxLen2];
                auto child1 = neigh->GetChild(chnum1);
                int c1 = child1->GetClipBoard();
                
                auto FE1 = this->get_fe(c1);
                auto FEDesc1_Obj = FE1.GetFEDesc();
                auto FEDesc1 = FE1.GetFEDesc_ID();
                
                int I_K1 = BeginIndex[c1];
                int m = TmpoFnlF[chnum1 * NFaces + joint_index_in_neigh];
                if(m == -1)
                {
                  std::cerr << "Error!\n";
                  exit(-1);
                }
                int * Indices1 = FEDesc1_Obj->GetJointDOF(m);
                
                auto mapper = mapper_db.get_fe_mapper(FEDesc0, FEDesc1);
                mapper->map(&GlobalNumbers[0], I_K0, I_K1, Indices0, Indices1,
                            Counter, joint->GetMapType());
              }
              else
              {
                ErrThrow("Only regular refinement are allowed!!");
              }
            }
            else
            {
              // neighbour is not refined

              /*int N_Points = TmpLen[j];
              double X = 0, Y = 0, Z = 0.0;
              double xp[4], yp[4], zp[4];
              for(int k=0;k<N_Points;k++)
              {
                cell->GetVertex(TmpFV[j*MaxLen+k])->GetCoords(xp[k], yp[k], zp[k]);
                X += xp[k];
                Y += yp[k];
                Z += zp[k];
              }
              X /= N_Points;
              Y /= N_Points;
              Z /= N_Points;
  
              // TODO: I don't know why we evaluate the boundary condition for 
              // an inner point. So the correct boundary component is also 
              // unclear, we use 0 here.
              BoundCond Cond0;
              boundCondition_(0, X, Y, Z, Cond0);
              
              auto mapper=mapper_db.get_fe_mapper(FEDesc0, FEDesc0);

              mapper->map_boundary(&GlobalNumbers[0], I_K0, Indices0,
                                   Counters[Cond0], HNNumbers);*/
              
              // HACK: previously we would create spurious boundary DOFs here
              // in the MPI case, causing problems because there are more than
              // the upper bound computed above, which in turn masked other
              // problems created by this approach. Calling map_boundary here
              // with the "inner" DOF counter results in "inner"/natural bc
              // DOFs at the outer edge of a process's halo only when they are
              // not actually on the boundary; the previous approach would
              // in some situations lose the actual boundary conditions on some
              // vertices, resulting in e.g. spurious and mismatched "pinhole"
              // outflows through the wall in TNSE problems.
              // This solution should be examined more closely, though, as it
              // may not interact correctly with other parts of the FEMapper in
              // all situations.
              
              auto mapper = mapper_db.get_fe_mapper(FEDesc0, FEDesc0);
              
              mapper->map_boundary(&GlobalNumbers[0], I_K0, Indices0,
                                   Counter, HNNumbers);
            }
          } // n == -1
          else
          {
            // neighbour is member of this collection
            // => using mappers
            if (n>i)
            {
              // this joint was not handled until now
              auto FE1 = this->get_fe(n);
              auto FEDesc1_Obj = FE1.GetFEDesc();
              auto FEDesc1 = FE1.GetFEDesc_ID();

              int I_K1 = BeginIndex[n];
              int * Indices1 = FEDesc1_Obj->GetJointDOF(joint_index_in_neigh);

              auto mapper = mapper_db.get_fe_mapper(FEDesc0, FEDesc1);
              mapper->map(&GlobalNumbers[0], I_K0, I_K1, Indices0, Indices1,
                          Counter, joint->GetMapType());
            } // n>i
          } // n != -1
        } // end neigh
      } // no boundary joint
    } // endfor j
 
#ifdef _MPI
    if(int N_EdgeDOF = FEDesc0_Obj->GetN_EdgeDOF())
    {
      int N_Edges=cell->GetN_Edges();
      const int * EdgeVertex;
      cell->GetShapeDesc()->GetEdgeVertex(EdgeVertex); 

      for(int j=0;j<N_Edges;j++)
      {
        auto edge = cell->GetEdge(j);
        if(edge->GetType()==BDEdge3D || edge->GetType()==IsoEdge3D)
        {
          if(edge->GetClipBoard()!=-1)
            continue;

          int * Indices0 = FEDesc0_Obj->GetEdgeDOF(j);   
          edge->SetClipBoard(i);
          int comp = edge->get_physical_id();

          double xp[2], yp[2], zp[2];
          cell->GetVertex(EdgeVertex[2*j])->GetCoords(xp[0], yp[0], zp[0]);
          cell->GetVertex(EdgeVertex[2*j+1])->GetCoords(xp[1], yp[1], zp[1]);
          double X = (xp[0] + xp[1])/2.;
          double Y = (yp[0] + yp[1])/2.;
          double Z = (zp[0] + zp[1])/2.;

          BoundCond Cond0;
          boundCondition_(comp, X, Y, Z, Cond0);

          auto mapper=mapper_db.get_fe_mapper(FEDesc0, FEDesc0);
          mapper->map_boundary_edge(N_EdgeDOF, &GlobalNumbers[0], I_K0,
                                    Indices0, Counters[Cond0], HNNumbers);
          /** vertices on this edge are already set, so no need of checking in 
           * BD vert on this edge, 30.06.12, sashi */
          cell->GetVertex(EdgeVertex[2*j])->SetClipBoard(i);
          cell->GetVertex(EdgeVertex[2*j+1])->SetClipBoard(i);
        }    // if(edge->GetType()==BDEdge3D || edge->GetType()==IsoEdge3D)
      }//  for(j=0;j<N_Edges;j 
    } // if(N_EdgeDOF = FEDesc0_Obj->GetN_EdgeDOF()>0)
  
    /** identify BD vert but not part of BD faces/edge */
    /** I hope, it is enough to check BD vertices, since other vertices are 
     * anyway inner, so no prob with Dirichlet BD */
    if(FEDesc0_Obj->GetN_VertDOF() > 0)
    {
      int N_VertInCell = cell->GetN_Vertices();
      for(int j=0;j<N_VertInCell;j++)
      {
        auto Vert = cell->GetVertex(j);
        if(Vert->IsBoundVert() && Vert->GetClipBoard()==-1 )
        {
          // BD vert is not yet set 
          Vert->SetClipBoard(i);
          int comp = Vert->get_physical_id();

          double X, Y, Z;
          Vert->GetCoords(X, Y, Z);
          BoundCond Cond0;
          boundCondition_(comp ,X, Y, Z, Cond0);

          auto mapper=mapper_db.get_fe_mapper(FEDesc0, FEDesc0); 
          mapper->map_boundary_vertex(&GlobalNumbers[0], I_K0,
                                      FEDesc0_Obj->GetVertDOF(j),
                                      Counters[Cond0], HNNumbers);
        }  //   if(Vert->IsBoundVert() && Vert->GetClipBoard()==-1 )
      } //   for(j=0;j<N_VertInCell;j++)
    } //    if(cell->IsBoundCell(
#endif // _MPI

    // handle inner degrees of freedom
    int k = FEDesc0_Obj->GetN_InnerDOF();
    int * Indices0 = FEDesc0_Obj->GetInnerDOF();
    for(int j=0;j<k;j++)
    {
     Counter--;
     GlobalNumbers[BeginIndex[i] + Indices0[j]] = Counter;
    } // endfor j 
  } // endfor i

// ===============================================================================================
// When the domain is partioned, cells in same subdomain may have only an edge(s) or vertex(ies) as common 
// between other cells in the same subcoll, therefore just joint map will not give a unique 
// globaldof on such edges or vertex
//  - sashikumaar Ganesan
// ===============================================================================================


#ifdef _MPI  
  int N_EdgeDOF = -1;
  int N_VertDof = -1;
  for(int i=0;i<N_Cells;i++)
  {
    auto cell = Collection->GetCell(i);
    auto FEDesc0_Obj = get_fe(i).GetFEDesc();

    if(cell->IsDependentCell())
    {
      if(N_EdgeDOF < FEDesc0_Obj->GetN_EdgeDOF())
        N_EdgeDOF = FEDesc0_Obj->GetN_EdgeDOF();

      if(N_VertDof < FEDesc0_Obj->GetN_VertDOF())
        N_VertDof = FEDesc0_Obj->GetN_VertDOF();

      int N_Edges=cell->GetN_Edges();

      for(int j=0;j<N_Edges;j++)
        (cell->GetEdge(j))->SetClipBoard(-1);
    }
  }// for(i=0;i<N_Cells;i++) 

  if (N_EdgeDOF > 0) // only for cont. FESpace
  {
    std::vector<int> w_array (N_EdgeDOF);
    for(int i = 0; i < N_Cells; i++)
    {
      auto cell = Collection->GetCell(i);
      if(cell->IsDependentCell())
      {

        int N_VertInCell = cell->GetN_Vertices();
        for(int j = 0; j < N_VertInCell; j++)
          (cell->GetVertex(j))->SetClipBoard(-1);

        int I_K0 = BeginIndex[i];
        auto FE0 = get_fe(i);
        auto FEDesc0_Obj = FE0.GetFEDesc();
        int N_Edges = cell->GetN_Edges();

        int N_EdgeDOF = FEDesc0_Obj->GetN_EdgeDOF();
        const int * EdgeVertex;
        (cell->GetShapeDesc())->GetEdgeVertex(EdgeVertex);

        for ( int j = 0; j < N_Edges; j++ )
        {
          auto Vert = cell->GetVertex(EdgeVertex[2 * j]); // start vertex of the edge
          auto edge = cell->GetEdge(j);

          if(edge->GetClipBoard() != -1)
            continue;

          edge->SetClipBoard(1);
          int * EdgeDof = FEDesc0_Obj->GetEdgeDOF(j);

          /** fill the already assigned global dof */
          for(int k = 0; k < N_EdgeDOF; k++)
          {
            w_array[k] = I_K0 +  EdgeDof[k];
            int v0;
            while((v0 = GlobalNumbers[w_array[k]]) > -1)
            { w_array[k] = v0; }
          } // for(k=0;k<N_EdgeDOF;k+

          const TBaseCell * const *EdgeNeibs;
          int N_EdgeNeibs;
          edge->GetNeibs(N_EdgeNeibs, EdgeNeibs);

          for(int k = 0; k < N_EdgeNeibs; k++)
          {
            auto neigh = EdgeNeibs[k];

            /** halo cell included, as multigrid fespace will contain Halo cells - Sashi:10-03-15 */
            if((rank != ( neigh->GetSubDomainNo())) && ! (neigh->IsHaloCell()))
              continue;

            /** own cell or neib cell (i.e. Halo cells) is not in this collection */
            if(neigh == cell || neigh->GetClipBoard_Par() == -1)
              continue;

            int l = 0; // find the edge local index in the neib cell
            while(edge != neigh->GetEdge(l)) l++;

            const int * NeibEdgeVertex;
            (neigh->GetShapeDesc())->GetEdgeVertex(NeibEdgeVertex);

            int maptype;
            if(neigh->GetVertex(NeibEdgeVertex[2 * l]) == Vert)
            { maptype =  1; }
            else if(neigh->GetVertex(NeibEdgeVertex[2 * l + 1]) == Vert)
            { maptype =  -1; }
            else
            {
              printf("FESpace3D Error in finding cross edge maptype Edge \n");
              MPI_Abort(MPI_COMM_WORLD, 0);
              ErrThrow("FESpace3D Error in finding cross edge maptype Edge");
            }

            int n = neigh->GetClipBoard();
            auto FE1 = get_fe(n);
            auto FEDesc1_Obj = FE1.GetFEDesc();
            auto NeibEdgeDof = FEDesc1_Obj->GetEdgeDOF(l);

            for(int m = 0; m < N_EdgeDOF; m++)
            {
              int neibdof;
              if(maptype == 1)
              {
                neibdof = BeginIndex[n] + NeibEdgeDof[m];
              }
              else
              {
                neibdof = BeginIndex[n] + NeibEdgeDof[N_EdgeDOF - 1 - m];
              }

              int w1 = neibdof;
              int v1;
              while((v1 = GlobalNumbers[w1]) > -1)
              {
                w1 = v1;
              }

              int w0 = w_array[m];
              if ( GlobalNumbers[w0] != GlobalNumbers[w1] )
              {
                int e = std::max(GlobalNumbers[w0], GlobalNumbers[w1]);
                w_array[m] = std::min(w0, w1);

                GlobalNumbers[w0] = w_array[m];
                GlobalNumbers[w1] = w_array[m];
                GlobalNumbers[w_array[m]] = e;
              }
            } // for(m=0;m<N_EdgeDOF;m++
          } //  for(k=0;k<N_EdgeNeibs;k++)
        } //  for(j=0;j<N_Edges;j+
      } // if(cell->IsDependentCell())
    } // for(i=0;i<N_Cells;i++)
  }//  if(N_EdgeDOF>0)  only for cont. FESpace

  if(N_VertDof > 0) // only for cont. FESpace
  {
    for(int i = 0; i < N_Cells; i++)
    {
      auto cell = Collection->GetCell(i);

      if(cell->IsDependentCell())
      {
        int I_K0 = BeginIndex[i];
        auto FE0 = get_fe ( i );
        auto FEDesc0_Obj = FE0.GetFEDesc();
        int N_VertInCell = cell->GetN_Vertices();

        for(int j = 0; j < N_VertInCell; j++)
        {
          auto Vert = cell->GetVertex(j);

          if(Vert->GetClipBoard() != -1)
            continue;

          Vert->SetClipBoard(5);
          int w0 = I_K0 + FEDesc0_Obj->GetVertDOF(j);

          int v0;
          while((v0 = GlobalNumbers[w0] ) > -1)
          {
            w0 = v0;
          }

          const TBaseCell * const * VertNeibs;
          int N_VertNeibs;
          Vert->GetNeibs (N_VertNeibs, VertNeibs);

          for(int k = 0; k < N_VertNeibs; k++)
          {
            auto neigh = VertNeibs[k];

            /** hallo cells are also need to be considered, so modified 11 Mar 2015 by Sashi */
            if((rank != neigh->GetSubDomainNo()) && ! (neigh->IsHaloCell()))
            {
              continue;
            }

            /** own cell or neib cell (i.e. Hallo cells) is not in this collection */
            if(neigh == cell || neigh->GetClipBoard_Par() == -1)
              continue;

            int n = neigh->GetClipBoard();
            auto FE1 = this->get_fe ( n );
            auto FEDesc1_Obj = FE1.GetFEDesc();

            int l = 0; // find the vert local index in the neib cell
            if(Vert->IsPeriodicVert())
            {
              while(Vert->GetPeriodicVertIndex()
                    != neigh->GetVertex(l)->GetPeriodicVertIndex())
                l++;
            }
            else
            {
              while(Vert != neigh->GetVertex(l))
                l++;
            }
            int w1 = BeginIndex[n] + FEDesc1_Obj->GetVertDOF(l);

            int v1;
            while((v1 = GlobalNumbers[w1] ) > -1)
            {
              w1 = v1;
            }

            if(GlobalNumbers[w0] != GlobalNumbers[w1])
            {
              int e = std::max(GlobalNumbers[w0], GlobalNumbers[w1]);
              int w = std::min(w0, w1);

              GlobalNumbers[w0] = w;
              GlobalNumbers[w1] = w;
              GlobalNumbers[w] = e;
              w0 = w;
            }
          } //for(k=0;k<N_VertNeibs;k++)
        } //  for(j=0;j<N_VertInCell;j++)
      } // if(cell->IsDependentCell())
    } // for(i=0;i<N_Cells;i++)
  }//  if(N_VertDof>0)

#endif // _MPI

  compute_various_numbers(BoundaryUpperBound, VHN, HNNumbers, false);
  
#ifdef _MPI
  // initialize the mapper and communicator for MPI communications
  mapper_.reset(new TParFEMapper3D(1, this));
  comm_ .reset(new TParFECommunicator3D(mapper_.get()));
#endif // _MPI
}

/** return position of all dofs */
void TFESpace3D::GetDOFPosition(int dof, double &x, double &y, double &z) const
{
  const double *xi, *eta, *zeta;
  int N_Points;

  if(dof > N_DegreesOfFreedom)
  {
    Output::warn("TFESpace3D::GetDOFPosition", dof, " dof number is larger "
                 "than total number of degrees of freedom");
    x = -1; y = -1; z = -1;
  }
  
  for(int i = 0, n_cells = Collection->GetN_Cells(); i < n_cells; i++)
  {
    const int *DOF = GetGlobalDOF(i);
    int k = BeginIndex[i+1] - BeginIndex[i];

    int DOFFound = -1;
    for(int j = 0; j < k; j++)
    {
      if(DOF[j] == dof) 
      {
        DOFFound = j;
        break;
      } // endif
    } // endfor


    if(DOFFound>-1) // i.e. dof was found
    {
      //cout << "dof " << dof << " found in cell: " << i << endl;
      auto cell  = Collection->GetCell(i);
      auto fe = get_fe(i);
      ReferenceTransformation_type RefTrans = fe.GetRefTransID();
  
      BFRefElements RefElement = fe.GetBaseFunct()->GetRefElement();
  
      auto nf = fe.GetNodalFunctional();
      nf->GetPointsForAll(N_Points, xi, eta, zeta);
      if(N_Points != fe.GetN_DOF())
      {
        ErrThrow("cannot get the position of a degree of freedom for this "
                 "finite element ", fe.GetID());
      }
  
      bool IsIsoparametric = TDatabase::ParamDB->USE_ISOPARAMETRIC
                             && cell->has_isoparametric_joint();
      if(IsIsoparametric)
      {
        switch(RefElement)
        {
          case BFRefElements::BFUnitHexahedron:
            RefTrans = ReferenceTransformation_type::HexaIsoparametric;
            break;
          case BFRefElements::BFUnitTetrahedron:
            RefTrans = ReferenceTransformation_type::TetraIsoparametric;
            break;
          default:
            ErrThrow("reference element ", RefElement, " not supported in 3D");
            break;
        }
      } // endif IsIsoparametric
  
      TRefTrans3D * rt = FEDatabase::GetRefTrans3D(RefTrans);
      rt->SetCell(cell);
      rt->GetOrigFromRef(xi[DOFFound], eta[DOFFound], zeta[DOFFound], x, y, z);
      break;
    } // endif DOFFound > -1
  } // endfor i
} // end GetDOFPosition

bool TFESpace3D::CheckMesh() const
{
  int N_DOF, found;
  const TBaseCell *Cell;
  int ActiveBound = get_n_active_non_hanging();

  for (int i = 0, n_cells = Collection->GetN_Cells(); i < n_cells; ++i)
  {
    Cell = Collection->GetCell(i);

    auto DOF = GetGlobalDOF(i);
    N_DOF = get_n_local_dof(i);

    found = 0;
    for (int j=0;j<N_DOF;++j)
    {
      if ( DOF[j] < ActiveBound )
      {
        found = 1;
        break;
      }
    }

    if ( found == 0 )
    {
      Output::print("index: " , i, 
                    "\n", Cell->GetVertex(0),
                    "\n", Cell->GetVertex(1),
                    "\n", Cell->GetVertex(2),
                    "\n", Cell->GetVertex(3));
      // Collection->info();      
      return false;
    }
  }

  return true;
}

BoundCond TFESpace3D::get_boundary_condition(const BoundaryJoint& bd_joint)
  const
{
  auto BoundFace = static_cast<const TBoundFace*>(&bd_joint);
  auto cell = BoundFace->GetNeighbour(0);
  int joint_id = BoundFace->get_index_in_neighbour(cell);
  auto BoundComp = BoundFace->GetBoundComp();
  int comp  = BoundComp->get_physical_id();
  const int *TmpFV, *TmpLen;
  int MaxLen;
  const int *TmpFE, *ETmpLen;
  int EMaxLen;
  cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
  cell->GetShapeDesc()->GetFaceEdge(TmpFE, ETmpLen, EMaxLen);
  int N_Points = TmpLen[joint_id];
  double xp[4], yp[4], zp[4];
//        double tp[4], sp[4];
//        double LinComb[4];
//         for(k=0;k<N_Points;k++)
//         {
//           LinComb[k] = 1.0/N_Points;
//           cell->GetVertex(TmpFV[joint_id*MaxLen+k])->GetCoords(xp[k], yp[k], zp[k]);
//         }
//         BoundFace->GetParameters(tp, sp);
// 
//         BoundComp->GetXYZandTS(N_Points, LinComb,
//                                xp, yp, zp, tp, sp,
//                                X, Y, Z, T, S);
  /** modified by sashi, since BD is not regular and T, S are not defined when tetgen is used */
  double X = 0, Y = 0, Z = 0;
  for(int k=0;k<N_Points;k++)
  {
    //LinComb[k] = 1.0/N_Points;
    cell->GetVertex(TmpFV[joint_id*MaxLen+k])->GetCoords(xp[k], yp[k], zp[k]);
    X += xp[k];
    Y += yp[k];
    Z += zp[k];
   }
  X /= N_Points;
  Y /= N_Points;
  Z /= N_Points;
  BoundCond Cond0;
  boundCondition_(comp, X, Y, Z, Cond0);
  return Cond0;
}
