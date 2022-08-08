// =======================================================================
// @(#)FESpace2D.C        1.18 06/27/00
// 
// Class:       TFESpace2D
// Purpose:     class for all 2D finite element spaces
//
// Author:      Gunar Matthies (04.11.97)
//
// History:     start of implementation 04.11.97 (Gunar Matthies)
//
//              split TFESpace into TFESpacexD (15.04.1998) Volker Behns
//              parallel methods  (Sashikumaar Ganesan) 09.09.09
// =======================================================================

#include <FESpace2D.h>
#include <Joint.h>
#include <BoundEdge.h>
#include <BoundComp2D.h>
#include "FEDatabase.h"
#include "BaseCell.h"
#include "HNDesc.h"
#include "FEMapperDatabase.h"
#include "FEMapper.h"
#include "FE2DMapper1Reg.h"

#include <Database.h>
#include <InterfaceJoint.h>
#include <NodalFunctional.h>

#include <MooNMD_Io.h>
#include <numeric> // accumulate
#include <algorithm>


TFESpace2D::TFESpace2D(const TCollection *coll, const std::string& name,
                       BoundCondFunct2D *BoundaryCondition, int ord) :
  TFESpace(coll, name, ord, 2), BoundCondition{BoundaryCondition}
{
  ConstructSpace();
}

TFESpace2D::TFESpace2D(const TCollection *coll,  const std::string& name,
                       BoundCondFunct2D *BoundaryCondition, SpaceType type,
                       int ord) :
  TFESpace(coll, name, type, ord, 2), BoundCondition{BoundaryCondition}
{
  ConstructSpace();
}

TFESpace2D::TFESpace2D(const TCollection *coll, const std::string& name,
                       BoundCondFunct2D *BoundaryCondition, FE_type *fes) :
  TFESpace(coll, name, fes), BoundCondition{BoundaryCondition}
{
  ConstructSpace();
}

std::array<int, TFESpace::N_DiffBoundNodeTypes> 
TFESpace2D::find_boundary_upper_bounds() const
{
  std::array<int, TFESpace::N_DiffBoundNodeTypes> BoundaryUpperBound{};
  int N_Cells = Collection->GetN_Cells();
  for(int i=0;i<N_Cells;i++)
  {
    auto cell = Collection->GetCell(i);
    // find number of boundary dofs for each joint of this cell
    int n_dof_per_joint = get_fe(i).GetFEDesc()->GetN_JointDOF();

    int n_joints = cell->GetN_Joints();
    for(int j = 0; j < n_joints; j++) // loop over all edges of cell
    {
      auto joint = cell->GetJoint(j);
      if(joint->GetType() == BoundaryEdge || joint->GetType() == IsoBoundEdge)
      {
        auto BoundEdge = (TBoundEdge *)joint;
        BoundCond Cond0 = get_boundary_condition(*BoundEdge);
        BoundaryUpperBound[Cond0] += n_dof_per_joint;
      } // endif
    } // endfor j
  } // endfor i
  return BoundaryUpperBound;
}

void TFESpace2D::ConstructSpace()
{
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
  // for i>0 set Counters[i] to Counters[i-1] - BoundaryUpperBound[i-1]
  std::transform(BoundaryUpperBound.begin(), BoundaryUpperBound.end(),
                 Counters.begin(), Counters.begin()+1,
                 [](int boundary_upper_bound, int previous_bound_counter)
                 { return previous_bound_counter - boundary_upper_bound; });
  // the following `Counter` is initialized as:
  //   some_special_number - sum of all BoundaryUpperBound
  int & Counter = Counters[TFESpace::N_DiffBoundNodeTypes];

  int N_Cells = Collection->GetN_Cells();
  // the following loop fills 'GlobalNumbers' with unique numbers which are not
  // yet the final ones. In particular, all numbers are either negative or they
  // are to be understood as the index (within 'GlobalNumbers') of a another dof
  // with which it identifies (i.e., these two dofs will end up being the same 
  // global dof). The initial value -1 means "not yet handled". Boundary dofs 
  // will have negative numbers which start at Counters[i] 
  // (0<i<TFESpace::N_DiffBoundNodeType depends on the type of boundary 
  // condition, starting with Dirichlet) and are smaller by 1 for each such
  // boundary dof. Negative numbers smaller than `Counter` will be inner dofs.
  for(int i=0;i<N_Cells;i++)
  {
    auto cell = Collection->GetCell(i);
    int N_Edges=cell->GetN_Edges();

    auto FE0 = this->get_fe(i);
    auto FEDesc0 = FE0.GetFEDesc_ID();
    auto FEDesc0_Obj = FE0.GetFEDesc();
    // the indices of the local dofs in this cell will be located at 
    // GlobalNumbers[I_K0], GlobalNumbers[I_K0+1], GlobalNumbers[I_K0+2], ...
    int I_K0 = BeginIndex[i];

    // handle cell-boundary degrees of freedom (inner dofs afterwards)
    for(int j=0;j<N_Edges;j++)
    {
      auto joint = cell->GetJoint(j);
      // 'Indices0' is the map from joint-local dofs (on joint j) to cell-local
      // dofs
      int * Indices0 = FEDesc0_Obj->GetJointDOF(j);
      if(joint->GetType() == BoundaryEdge || joint->GetType() == IsoBoundEdge)
      {
        // boundary joint
        // there is no need to identify dofs from different cells here. Dofs on
        // this boundary joint which are also in another cell must be on another
        // (inner) joint as well, where the identification is handled.
        auto BoundEdge = (TBoundEdge *)joint;
        BoundCond Cond0 = get_boundary_condition(*BoundEdge);
        auto mapper=mapper_db.get_fe_mapper(FEDesc0, FEDesc0);
        mapper->map_boundary(&GlobalNumbers[0], I_K0, Indices0,
                             Counters[Cond0], HNNumbers);
      } // boundary joint
      else
      {
        // no boundary joint
        const auto* neigh = joint->GetNeighbour(cell);
        if (!neigh)
        {
          // there is no neighbour
          // => finer cell in 1 regular grid (hanging vertex)
          //    will be handled from coarser cell
        } // !neigh
        else
        {
          // there is a neighbour on same level
          // find local edge index in `neigh` of this joint
          int joint_index_in_neigh=0;
          while(neigh->GetJoint(joint_index_in_neigh) != joint)
            joint_index_in_neigh++;

          int n = neigh->GetClipBoard();
          if(n == -1)
          {
            // neighbor is not member of this collection, this indicates a 
            // hanging vertex
            // check for children of neigh on its joint joint_index_in_neigh
            auto refdesc=neigh->GetRefDesc();
            if(refdesc->GetType() != NoRef)
            {
              const int *TmpoEnE;
              const int *TmpLen1;
              int MaxLen1;
              const int *TmpEC;
              const int *TmpLen2;
              int MaxLen2;
              const int *TmpoEnlE;
              refdesc->GetOldEdgeNewEdge(TmpoEnE, TmpLen1, MaxLen1);
              refdesc->GetEdgeChild(TmpEC, TmpLen2, MaxLen2);
              refdesc->GetOldEdgeNewLocEdge(TmpoEnlE);
              // number of edges in parent cell (i.e., 3 for triangles or 4 for
              // quads), equals neigh->GetN_Edges()
              int NEdges=refdesc->GetShapeDesc()->GetN_Edges();

              // e1 is the edge index (within the local refinement descriptor)
              // of the first child on this edge
              int e1 = TmpoEnE[joint_index_in_neigh*MaxLen1+0];
              // chnum1 is the index (within the local refinement descriptor) of
              // the child cell neighboring `e1`
              int chnum1=TmpEC[e1*MaxLen2];
              // the child cell neighboring `e1`
              auto child1=neigh->GetChild(chnum1);
              // the index of this `child1` within the collection of this space
              int c1=child1->GetClipBoard();
              if(c1 == -1)
              {
                // This should never happen, I think. possibly if the finer
                // cell adjacent to a hanging vertex has been refined with
                // e.g. barycentric refinement.
                ErrThrow("child of neighbor not in the current collection");
              }
              // the finite element on cell `child1`
              auto FE1 = get_fe(c1);
              auto FEDesc1 = FE1.GetFEDesc_ID();
              auto FEDesc1_Obj = FE1.GetFEDesc();
              int I_K1 = BeginIndex[c1];
              // m is the edge index within the cell `chnum1`
              int m=TmpoEnlE[chnum1*NEdges+joint_index_in_neigh];
              // the dofs on this edge
              int * Indices1 = FEDesc1_Obj->GetJointDOF(m);

              // number of inner vertices on joint
              int n_inner_vertices = TmpLen1[joint_index_in_neigh];
              if(n_inner_vertices == 0)
              {
                // there is only ONE child on this edge, no hanging vertex here
                auto mapper = mapper_db.get_fe_mapper(FEDesc0, FEDesc1);
                mapper->map(&GlobalNumbers[0], I_K0, I_K1, Indices0, Indices1,
                            Counter);
              } // n_inner_vertices == 0
              else if(n_inner_vertices == 1)
              {
                // there are exactly two children on this edge, i.e., the is a 
                // hanging vertex
                // e2 is the edge index (within the local refinement descriptor)
                // of the second child on this edge
                int e2=TmpoEnE[joint_index_in_neigh*MaxLen1+1];
                // chnum2 is the index (within the local refinement descriptor)
                // of the child cell neighboring `e2`
                int chnum2=TmpEC[e2*MaxLen2];
                // the child cell neighboring `e2`
                auto child2=neigh->GetChild(chnum2);
                // the index of this `child2` within the collection of this space
                int c2=child2->GetClipBoard();

                if(c2 == -1)
                  ErrThrow("child of neighbor not in the current collection");

                // the finite element on cell `child2`
                auto FE2 = get_fe(c2);
                auto FEDesc2 = FE2.GetFEDesc_ID();
                auto FEDesc2_Obj = FE2.GetFEDesc();
                int I_K2 = BeginIndex[c2];
                // m is the edge index within the cell `chnum2`
                int m2 = TmpoEnlE[chnum2*NEdges+joint_index_in_neigh];
                int *Indices2 = FEDesc2_Obj->GetJointDOF(m2);

                auto mapper1reg = mapper_db.get_fe_mapper_one_regular2d(FEDesc0,
                                                                        FEDesc2);
                mapper1reg->Map_1Reg(i<c1, &GlobalNumbers[0], I_K0, I_K1, I_K2,
                                     Indices0, Indices1, Indices2, Counter,
                                     (c1<c2), VHN, HNNumbers, HN_descriptors);
              } // n_inner_vertices == 1
              else
              {
                ErrThrow("more than two children on one edge are not allowed");
              }
            }
            else
            {
              // neighbour is not refined
              // the dofs on this joints are treated almost like boundary dofs
              // except that we use `Counter` which counts inner dofs.
              // This is a rather rare situation, it could be that the
              // collection only covers a part of the domain
              auto mapper=mapper_db.get_fe_mapper(FEDesc0, FEDesc0);
              mapper->map_boundary(&GlobalNumbers[0], I_K0, Indices0, Counter,
                                   HNNumbers);
            }
          } // n == -1
          else
          {
            // neighbour is member of this collection
            // => using mappers
            if (n>i)
            {
              // this joint was not handled until now
              auto& FE1 = get_fe(n);
              auto FEDesc1 = FE1.GetFEDesc_ID();
              auto FEDesc1_Obj = FE1.GetFEDesc();
              // starting index within `GlobalNumbers` for the dofs of the 
              // neighbor
              int I_K1 = BeginIndex[n];
              // 'Indices1' is the map from joint-local dofs (on joint 
              // `joint_index_in_neigh`) to cell-local dofs in the neighbor
              int * Indices1 = FEDesc1_Obj->GetJointDOF(joint_index_in_neigh);
              auto mapper = mapper_db.get_fe_mapper(FEDesc0, FEDesc1);
              mapper->map(&GlobalNumbers[0], I_K0, I_K1, Indices0, Indices1,
                          Counter);
            } // n>i
          } // n != -1
        } // end neigh
      } // no boundary joint
    } // endfor j

    // handle cell-inner degrees of freedom
    int k = FEDesc0_Obj->GetN_InnerDOF();
    int * Indices0 = FEDesc0_Obj->GetInnerDOF();
    for(int j=0;j<k;j++)
    {
      Counter--;
      GlobalNumbers[BeginIndex[i] + Indices0[j]] = Counter;
    } // endfor j
  } // endfor i

  compute_various_numbers(BoundaryUpperBound, VHN, HNNumbers, true);
}

/** return position of one given dof */
void TFESpace2D::GetDOFPosition(int dof, double &x, double &y) const
{
  const double *xi, *eta;
  int N_Points;

  if(dof > N_DegreesOfFreedom)
  {
    ErrThrow("dof number is larger than total number of degrees of freedom");
  }

  for(int i = 0, n_cells = Collection->GetN_Cells(); i < n_cells; i++)
  {
    const int *DOF = this->TFESpace::GetGlobalDOF(i);
    int k = get_n_local_dof(i);

    int DOFFound = -1;
    for(int j=0;j<k;j++)
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
      nf->GetPointsForAll(N_Points, xi, eta);
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
          case BFRefElements::BFUnitSquare:
           RefTrans = ReferenceTransformation_type::QuadIsoparametric;
            break;
          case BFRefElements::BFUnitTriangle:
            RefTrans = ReferenceTransformation_type::TriaIsoparametric;
            break;
          default:
            ErrThrow("reference element ", RefElement, " not supported in 2D");
            break;
        }
      } // endif IsIsoparametric
  
      auto rt = FEDatabase::GetRefTrans2D(RefTrans);
      rt->SetCell(cell);
      rt->GetOrigFromRef(xi[DOFFound], eta[DOFFound], x, y);
      break;
    } // endif DOFFound > -1
  } // endfor i
} // end GetDOFPosition

BoundCond TFESpace2D::get_boundary_condition(const BoundaryJoint& bd_joint)
  const
{
  constexpr double eps = 1e-6;
  auto BoundEdge = static_cast<const TBoundEdge*>(&bd_joint);
  auto BoundComp = BoundEdge->GetBoundComp();
  double t0, t1;
  BoundEdge->GetParameters(t0, t1);
  int comp=BoundComp->GetID();
  BoundCond Cond0, Cond1;
  if (t0 < t1)
  {
    BoundCondition(comp, t0+eps, Cond0);
    BoundCondition(comp, t1-eps, Cond1);
  }
  else
  {
    BoundCondition(comp, t0-eps, Cond0);
    BoundCondition(comp, t1+eps, Cond1);
  }

  if(Cond0 == Cond1)
  {
    return Cond0;
  }
  else
  {
    ErrThrow("different boundary condition on one edge are not allowed");
  }
}


/** check if FE spaces lhs_space and rhs_space are equal*/
bool operator==(const TFESpace2D &lhs_space, const TFESpace2D &rhs_space)
{
  if(&lhs_space == &rhs_space) // compare pointers
    return true;
  if((lhs_space.N_DegreesOfFreedom == rhs_space.N_DegreesOfFreedom)
     && (lhs_space.UsedElements == rhs_space.UsedElements)
     && (lhs_space.BoundCondition == rhs_space.BoundCondition)
     //&& (lhs_space.OrderOfSpace == rhs_space.OrderOfSpace)
     && (lhs_space.Collection == rhs_space.Collection))
  {
    return true;
  }
  
  return false;
}

/** check if FE spaces lhs_space and rhs_space are not equal*/
bool operator!=(const TFESpace2D &lhs_space, const TFESpace2D &rhs_space)

{
  return !(lhs_space == rhs_space);
}


