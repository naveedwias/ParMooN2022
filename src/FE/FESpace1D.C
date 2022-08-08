#include "FESpace1D.h"
#include "Joint.h"
#include "BaseCell.h"
#include "FEMapper.h"
#include "MooNMD_Io.h"
#include <array> // array for MacOS
#include <cstring> // memset
#include <numeric> // accumulate


TFESpace1D::TFESpace1D(const TCollection *coll, const std::string& name,
                       int ord)
: TFESpace(coll, name, ord, 1)
{
  ConstructSpace();
}

/** constructor for building a space with the given elements */
TFESpace1D::TFESpace1D(const TCollection *coll, const std::string& name,
                         FE_type *fes)
 : TFESpace(coll, name, fes)
{
  ConstructSpace();
}

void TFESpace1D::ConstructSpace()
{
  // find out all used elements
  // boundary structures (no boundary values used)
  N_BoundaryNodes[0]=0;
  N_BoundaryNodes[1]=0;

  int N_Cells = Collection->GetN_Cells();
  // start DOF manager
  int l = FIRSTMARK + 1;
  int Counter = l;

  for(int i=0;i<N_Cells;i++)
  {
    auto cell = Collection->GetCell(i);

    int N_Joints = cell->GetN_Vertices();

    auto FE0 = get_fe(i);
    auto FEDesc0_Obj = FE0.GetFEDesc();
    int I_K0 = BeginIndex[i];
    auto J_K0 = &GlobalNumbers[BeginIndex[i]];

    for(int j=0;j<N_Joints;j++)
    {
      auto joint = cell->GetJoint(j);
      auto Indices0 = FEDesc0_Obj->GetJointDOF(j);
      int N_JointDOF = FEDesc0_Obj->GetN_JointDOF();
      
      // if (N_JointDOF)
      if(joint)
       {
         auto neigh = joint->GetNeighbour(cell);
        if(!neigh)
        {
          // boundary joint
          for(int k=0;k<N_JointDOF;k++)
          {
            if(J_K0[Indices0[k]]==-1)
            {
              // not handled yet
              Counter--;
              J_K0[Indices0[k]] = Counter;
            }
          }
        } // boundary joint
       else
        {
          // no boundary joint
          int n = neigh->GetClipBoard();
          if (n>i && N_JointDOF!=0)
          {
            // this joint was not handled until now

//            I_K1 = BeginIndex[n];
            auto J_K1 = &GlobalNumbers[BeginIndex[n]];
            // find the local edge of neigh on which cell is
            int l=0;
            // while(neigh->GetJoint(l)->GetNeighbour(neigh)!=cell) l++;
            while(neigh->GetJoint(l) != joint) l++;

            auto FE1 = get_fe(n);
            auto FEDesc1_Obj = FE1.GetFEDesc();
            auto Indices1 = FEDesc1_Obj->GetJointDOF(l);

            // cout << "i= " << i << endl;
            // cout << "j= " << j << endl;
            // cout << "n= " << n << endl;
            // cout << "l= " << l << endl;

            int k = J_K0[Indices0[0]];
            if(k<=FIRSTMARK)
            {
              // DOF already handled
              J_K1[Indices1[0]] = I_K0+Indices0[0];
            }
            else
            {
              // DOF not handle yet
              J_K1[Indices1[0]] = I_K0+Indices0[0];
              Counter--;
              J_K0[Indices0[0]] = Counter;
            }
          } // end neigh
        } // no boundary joint
       }
    } // endfor j

    // handle inner degrees of freedom
    int k = FEDesc0_Obj->GetN_InnerDOF();
    auto Indices0 = FEDesc0_Obj->GetInnerDOF();
    for(int j=0;j<k;j++)
    {
      Counter--;
      J_K0[Indices0[j]] = Counter;
    } // endfor j
  } // endfor i

  // find global numbers
  l=0;
  int count = l;

  N_Inner = 0;

  for(int i=0;i<N_Cells;i++)
  {
    auto J_K0 = &GlobalNumbers[BeginIndex[i]];
    int k = get_fe(i).GetN_DOF();
    for(int j=0;j<k;j++)
    {
      l = J_K0[j];
      if (l < -1)
      {
        // node to assign here
        J_K0[j] = count;
        count++;
        N_Inner++;
      } // l < -1
      else
      {
        if (l >= 0)
        {
          J_K0[j] = GlobalNumbers[l];
        }
        else
        {
          cerr << "J_K0[j]==-1" << endl;
        }
      } // l >= -1
    } // endfor j
  } // endfor i

/*
  // print for all elements for global numbers of their local dofs
  for(i=0;i<N_Cells;i++)
  {
    cell = Collection->GetCell(i);
    cout << "cell number: " << i << endl;

    J_K0 = &GlobalNumbers[BeginIndex[i]];
    auto FE0 = get_fe(i, cell);
    k = FE0.GetN_DOF();
    for(j=0;j<k;j++)
    {
      cout << j << ": " << " number: " << J_K0[j] << endl;
    }
    cout << endl;
  } // endfor i
*/

  l = N_Inner;
  N_DegreesOfFreedom = l;
}
