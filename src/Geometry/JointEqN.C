// =======================================================================
// @(#)JointEqN.C        1.8 11/15/99
// 
// Class:       TJointEqN
// Purpose:     connects two cells
//
// Author:      Volker Behns  23.07.97
//
// =======================================================================

#include <BaseCell.h>
#include <JointEqN.h>
#include <RefDesc.h>
#include "Mapper.h"

// Constructors
TJointEqN::TJointEqN(TBaseCell *neighb0) : TJoint()
{
  ID = JointEqN;

  Neighb0 = neighb0;
}

TJointEqN::TJointEqN(TBaseCell *neighb0, TBaseCell *neighb1) : TJoint()
{
  ID = JointEqN;

  Neighb0 = neighb0;
  Neighb1 = neighb1;
}

    // Destructor
TJointEqN::~TJointEqN()
   {
    if(Neighb0)
     { Neighb0 = nullptr;}
  
    if(Neighb1)
     { Neighb1 = nullptr;}   
   }
   
   
// Methods

int TJointEqN::CheckMatchingRef(TBaseCell *Me, int J_i, struct StoreGeom &Tmp)
{
  const TBaseCell *Neighb, *Parent = nullptr;
  int i, J_j, N_, MaxLen1, MaxLen2, auxi, aux;
  const int *TmpValues1, *TmpValues2;
  const TRefDesc *NeighbRefDesc, *ParRefDesc;
  const TJoint *ParJoint;
#ifdef __2D__
  Refinements NeibEdgeRef, MyEdgeRef;
  const int *TmpoEnE, *TmpEC, *TmpECI, *TmpoEnV, *TmpVC, *TmpVCI;
#else
  Refinements NeibFaceRef, MyFaceRef;
  const int *TmpLen1, *TmpLen2, *TmpLen3; // *TmpLen4;
  const int *TmpoFnF, *TmpFC, *TmpFCI, *TmpoFnV, *TmpVC, *TmpVCI;
  const int *MapVerts, *MapFaces, *MapEdges;
  int* MapFacesTmp;
  int MaxLen3;
#endif

  Tmp.Filled = false;

  if (Neighb0 == Me)
    Neighb = Neighb1;
  else
    Neighb = Neighb0;

  if (!Neighb)
  {
    Parent = Me;
    while ((Parent = Parent->GetParent()))
    {
      N_ = Parent->GetN_Children();
      for (i=0;i<N_;i++)
        if (Parent->GetChild(i) == Me) break;

      ParRefDesc = Parent->GetRefDesc();
#ifdef __2D__
      ParRefDesc->GetChildEdge(TmpValues1, MaxLen1);
      ParRefDesc->GetNewEdgeOldEdge(TmpValues2);
#else
      ParRefDesc->GetChildFace(TmpValues1, MaxLen1);
      ParRefDesc->GetNewFaceOldFace(TmpValues2);
#endif

      J_i = TmpValues2[TmpValues1[MaxLen1*i + J_i]];

      // we can go deeper only if nothing happens (NoRef)
#ifdef __2D__
      if (Parent->GetEdgeRef(J_i) != NoRef) return 0;
#else
      if (Parent->GetFaceRef(J_i) != NoRef) return 0;
#endif

      if ((Neighb = Parent->GetJoint(J_i)->GetNeighbour(Parent)))
        break;
    }

    if (!Parent) return 0;
    ParJoint = Parent->GetJoint(J_i);
  }

  if (Neighb->ExistChildren())
  {

#ifdef __2D__
    N_ = Neighb->GetN_Edges();
#else
    N_ = Neighb->GetN_Faces();
#endif
    if (Parent)
    {
      for (J_j=0;J_j<N_;J_j++)
        if (Neighb->GetJoint(J_j) == ParJoint) break;
    }
    else
      for (J_j=0;J_j<N_;J_j++)
        if (Neighb->GetJoint(J_j) == this) break;

#ifdef __2D__
    NeibEdgeRef = Neighb->GetEdgeRef(J_j);
    MyEdgeRef = Me->GetEdgeRef(J_i);

    if (NeibEdgeRef != MyEdgeRef)
       cout << "test JointEqn" << endl;
  
    if (NeibEdgeRef != MyEdgeRef)
      return -2;
   
    Tmp.Filled = true;

    NeighbRefDesc = Neighb->GetRefDesc();

    NeighbRefDesc->GetOldEdgeNewEdge(TmpoEnE, TmpValues1, MaxLen1);
    N_ = TmpValues1[J_j] + 1;

    NeighbRefDesc->GetEdgeChild(TmpEC, TmpValues1, MaxLen2);
    NeighbRefDesc->GetEdgeChildIndex(TmpECI, TmpValues1, MaxLen2);

    for (i=0;i<N_;i++)
    {
      auxi = TmpoEnE[J_j * MaxLen1 + i] * MaxLen2;
      Tmp.Joints[i] = Neighb->GetChild(TmpEC[auxi])->
                              GetJoint(TmpECI[auxi]);
    }

    NeighbRefDesc->GetOldEdgeNewVertex(TmpoEnV, TmpValues1, MaxLen1);
    NeighbRefDesc->GetVertexChild(TmpVC, TmpValues1, MaxLen2);
    NeighbRefDesc->GetVertexChildIndex(TmpVCI, TmpValues1, MaxLen2);

    auxi = J_j * MaxLen1 + 1;
    for (i=1;i<N_;i++)
    {
      aux = TmpoEnV[auxi++] * MaxLen2;
      Tmp.Vertices[i] = Neighb->GetChild(TmpVC[aux])->
                          GetVertex(TmpVCI[aux]);
    }
#else
    NeibFaceRef = Neighb->GetRefDesc()->GetFaceRef(J_j);
    MyFaceRef = Me->GetRefDesc()->GetFaceRef(J_i);

    if((NeibFaceRef == TriReg && MyFaceRef != TriReg) ||
       (NeibFaceRef != TriReg && MyFaceRef == TriReg) )
        return -2;
    if(MyFaceRef >= TriBis0 && MyFaceRef <= TriBis2)
    {
      if(NeibFaceRef != TriBis0 + (2 - (MyFaceRef - TriBis0) + MapType) % 3 )
      {
        std::cerr << "Face Reference Description for Triangle Bisection does not fit\n";
        return -2;
      }
      MapType = MapTriBis00 + 3*(MyFaceRef - TriBis0) + MapType;
    }
    if(MyFaceRef >= TriBis01 && MyFaceRef <= TriBis21)
    {
      int first_bis = (MyFaceRef-TriBis01) / 2;
      int second_bis = (MyFaceRef-TriBis01) % 2;
      if(first_bis <= second_bis) 
        second_bis++;
    
      int neigh_first_bis = ((2-first_bis) + MapType) % 3;
      int neigh_second_bis = ((2-second_bis) + MapType) % 3;
    
      int neigh_type = 2*neigh_first_bis + neigh_second_bis;
      if(neigh_second_bis > neigh_first_bis) 
        neigh_type--;
    
      if(NeibFaceRef != neigh_type + TriBis01)
      {
        std::cerr << "Face Reference Description for Triangle Bisection does not fit\n";
        return -2;
      }
    
      MapType = MapTriBis010 + 3 * (NeibFaceRef - TriBis01) + MapType;
    }


    Tmp.Filled = true;

    NeighbRefDesc = Neighb->GetRefDesc();

    NeighbRefDesc->GetOldFaceNewFace(TmpoFnF, TmpLen1, MaxLen1);
    N_ = TmpLen1[J_j];

    NeighbRefDesc->GetFaceChild(TmpFC, TmpLen2, MaxLen2);
    NeighbRefDesc->GetFaceChildIndex(TmpFCI, TmpLen2, MaxLen2);

    NeighbRefDesc->GetOldFaceNewVertex(TmpoFnV, TmpLen3, MaxLen3);

    if(N_ == 1)
    {
      MapFacesTmp = new int[1];
      MapFacesTmp[0] = 0;
      MapFaces = MapFacesTmp;
      
      Neighb->GetJoint(J_j)->GetMapperOrig(MapVerts, MapEdges);
    }
    else
      Neighb->GetJoint(J_j)->GetMapperRef(MapVerts, MapFaces);

    auxi = J_j * MaxLen1;
    for (i=0;i<N_;i++)
    {
      aux = TmpoFnF[auxi++] * MaxLen2;
      Tmp.Joints[MapFaces[i]] = Neighb->
            GetChild(TmpFC[aux])->GetJoint(TmpFCI[aux]);
    }

    NeighbRefDesc->GetVertexChild(TmpVC, TmpValues1, MaxLen2);
    NeighbRefDesc->GetVertexChildIndex(TmpVCI, TmpValues1, MaxLen2);

    N_ = TmpLen3[J_j];
    auxi = J_j * MaxLen3;
    for (i=0;i<N_;i++)
    {
      aux = TmpoFnV[auxi++] * MaxLen2;
      Tmp.Vertices[MapVerts[i]] = Neighb->GetChild(TmpVC[aux])->
                          GetVertex(TmpVCI[aux]);
    }
#endif // __2D__
  }

  return 0;
}
