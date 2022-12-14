// =======================================================================
// @(#)InterfaceJoint.C        1.2 10/18/99
// 
// Class:       TInterfaceJoint
// Purpose:     connects two cells on an interface
//
// Author:      Volker Behns  10.03.98
//
// =======================================================================
#include <cmath>

#include <InterfaceJoint.h>
#include <BoundComp2D.h>
#include <BaseCell.h>

// Constructors
TInterfaceJoint::TInterfaceJoint(const TBoundComp2D *bdcomp, double t_0,
                 double t_1, TBaseCell *neighb0) : TJointEqN(neighb0)
{
  ID = InterfaceJoint;

  BoundComp = bdcomp;

  if (t_0 < t_1)
  {
    T_0 = t_0;
    T_1 = t_1;
  }
  else
  {
    T_0 = t_1;
    T_1 = t_0;
  }
}

TInterfaceJoint::TInterfaceJoint(const TBoundComp2D *bdcomp, double t_0,
                 double t_1, TBaseCell *neighb0, TBaseCell *neighb1) :
                 TJointEqN(neighb0, neighb1)
{
  ID = InterfaceJoint;

  BoundComp = bdcomp;

  if (t_0 < t_1)
  {
    T_0 = t_0;
    T_1 = t_1;
  }
  else
  {
    T_0 = t_1;
    T_1 = t_0;
  }
}

// Methods
TJoint *TInterfaceJoint::NewInst(double newT_0, double newT_1, TBaseCell *Me)
{
  return new TInterfaceJoint(BoundComp, T_0 + newT_0*(T_1 - T_0),
                             T_0 + newT_1*(T_1 - T_0), Me);
}

TJoint *TInterfaceJoint::NewInst()
{
  return new TInterfaceJoint(BoundComp, T_0, T_1, nullptr);
}

int TInterfaceJoint::GetXYofT(double T, double &X, double &Y)
{
  BoundComp->GetXYofT(T, X, Y);

  return 0;
}

#ifdef __2D__
/** update parameters according to the new vertex positions */
void TInterfaceJoint::UpdateParameters(const TVertex *Begin,
                                       const TVertex *End)
{
  double x1, y1, x2, y2;
  double t1, t2;

#ifdef __2D__
  Begin->GetCoords(x1, y1);
  End->GetCoords(x2, y2);
#else
  double z1, z2;
  Begin->GetCoords(x1, y1, z1);
  End->GetCoords(x2, y2, z2);
#endif

  BoundComp->GetTofXY(x1, y1, t1);
  BoundComp->GetTofXY(x2, y2, t2);

  T_0 = t1;
  T_1 = t2;
}
#endif

int TInterfaceJoint::CheckOrientation()
{
  int i, N_;
  const int *TmpEV;
  TVertex *Vert0;
  double X, Y;
  TBaseCell *Neighb2;

  if (Neighb0)
  {
    N_ = Neighb0->GetRefDesc()->GetN_OrigEdges();
    for (i=0;i<N_;i++)
      if (Neighb0->GetJoint(i) == this) break;

    Neighb0->GetRefDesc()->GetShapeDesc()->GetEdgeVertex(TmpEV);

    Vert0 = Neighb0->GetVertex(TmpEV[2*i]);
    BoundComp->GetXYofT(T_0, X, Y);

    if (std::abs(Vert0->GetX() - X) > 1e-6 || std::abs(Vert0->GetY() - Y) > 1e-6)
    {
      Neighb2 = Neighb1;
      Neighb1 = Neighb0;
      Neighb0 = Neighb2;
    }
  }
  else
  {
    cerr << "Erorr in InterfaceJoint: no neighbour given!" << endl;
    return -1;
  }

  return 0;
}
