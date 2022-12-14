// =======================================================================
// @(#)InterfaceJoint3D.C        1.2 10/18/99
// 
// Class:       TInterfaceJoint3D
// Purpose:     connects two cells on an interface
//
// Author:      Volker Behns  10.03.98
//              Gunar Matthies 05.04.02
//
// =======================================================================

#include <InterfaceJoint3D.h>
#include <BoundComp3D.h>
#include <BaseCell.h>
#include <cmath>

// Constructors
/** initialize the joint with the boundary parameters and one neighbour */
TInterfaceJoint3D::TInterfaceJoint3D(const TBoundComp3D *bdcomp,
                                 double *param1, double *param2,
                                 TBaseCell *neigh0)
   : TJointEqN(neigh0)
{
  int i;
  
  ID = InterfaceJoint3D;

  BoundComp = bdcomp;
  for(i=0;i<4;i++)
  {
    Param1[i] = param1[i];
    Param2[i] = param2[i];
  }

  Neighb0 = neigh0;
}

/** initialize the joint with the boundary parameters and two neighbours */
TInterfaceJoint3D::TInterfaceJoint3D(const TBoundComp3D *bdcomp,
                                 double *param1, double *param2,
                                 TBaseCell *neigh0, TBaseCell *neigh1)
   : TJointEqN(neigh0, neigh1)
{
  int i;
  
  ID = InterfaceJoint3D;

  BoundComp = bdcomp;
  for(i=0;i<4;i++)
  {
    Param1[i] = param1[i];
    Param2[i] = param2[i];
  }

  Neighb0 = neigh0;
  Neighb1 = neigh1;
}

TInterfaceJoint3D::TInterfaceJoint3D(const TBoundComp3D *bdcomp,
                                 TBaseCell *neigh0)
   : TJointEqN(neigh0)
{
  ID = InterfaceJoint3D;

  BoundComp = bdcomp;
  Param1[0] = 0.0;
  Param2[0] = 0.0;
  Param1[1] = 1.0;
  Param2[1] = 0.0;
  Param1[2] = 1.0;
  Param2[2] = 1.0;

  Neighb0 = neigh0;
}

    // Destructor
TInterfaceJoint3D::~TInterfaceJoint3D()
   {
    if(Neighb0)
     { Neighb0 = nullptr;}
  
    if(Neighb1)
     { Neighb1 = nullptr;}   
   }
   
  
  
// Methods
TJoint *TInterfaceJoint3D::NewInst(double, double, TBaseCell *Me)
{
  return new TInterfaceJoint3D(BoundComp, Me);
}

TJoint *TInterfaceJoint3D::NewInst()
{
  return new TInterfaceJoint3D(BoundComp, nullptr);
}

/** return both parameter vectors */
void TInterfaceJoint3D::GetParameters(double *param1, double *param2)
{
  int i;

  for(i=0;i<4;i++)
  {
    param1[i] = Param1[i];
    param2[i] = Param2[i];
  }
}

/** set both parameter vectors */
void TInterfaceJoint3D::SetParameters(double *param1, double *param2)
{
  int i;

  for(i=0;i<4;i++)
  {
    Param1[i] = param1[i];
    Param2[i] = param2[i];
  }
}

int TInterfaceJoint3D::CheckOrientation()
{
  int i,j, N_;
  const int *TmpFV, *TmpLen;
  int MaxLen;
  TVertex *Vert;
  double X, Y, Z, xp, yp, zp;
//  TBaseCell *Neighb2;

  if (Neighb0)
  {
    N_ = Neighb0->GetRefDesc()->GetN_OrigFaces();
    for (i=0;i<N_;i++)
      if (Neighb0->GetJoint(i) == this) break;

    Neighb0->GetRefDesc()->GetShapeDesc()
        ->GetFaceVertex(TmpFV, TmpLen, MaxLen);

    for(j=0;j<TmpLen[i];j++)
    {
      Vert = Neighb0->GetVertex(TmpFV[i*MaxLen+j]);
      Vert->GetCoords(X, Y, Z);
      BoundComp->GetXYZofTS(Param1[j], Param2[j], xp, yp, zp);
      if(std::abs(X-xp)>1e-8 || std::abs(Y-yp)>1e-8 || std::abs(Z-zp)>1e-8)
      {
        BoundComp->GetTSofXYZ(X, Y, Z, Param1[j], Param2[j]);
      }
      // cout << "X, Y, Z: " << X << " " << Y << " " << Z;
      // cout << " P: " << Param1[j] << " " << Param2[j] << endl;
    } // endfor j
  }
  else
  {
    cerr << "Erorr in InterfaceJoint3D: no neighbour given!" << endl;
    return -1;
  }

  return 0;
}
