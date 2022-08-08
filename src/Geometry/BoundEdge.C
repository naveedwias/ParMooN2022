#include "BoundEdge.h"
#include "MooNMD_Io.h"

// Constructors
TBoundEdge::TBoundEdge(const TBoundComp2D *bdcomp, double t_0, double t_1)
{
  ID = BoundaryEdge;

  BoundComp = bdcomp;
  T_0 = t_0;
  T_1 = t_1;
}

// create a new instance of this class
TJoint *TBoundEdge::NewInst(double newT_0, double newT_1, TBaseCell *)
{
  return new TBoundEdge(BoundComp, T_0 + newT_0*(T_1 - T_0),
                        T_0 + newT_1*(T_1 - T_0));
}

TJoint *TBoundEdge::NewInst()
{
  return new TBoundEdge(BoundComp, T_0, T_1);
}

