#include "BoundFace.h"
#include "MooNMD_Io.h"

// Constructors
TBoundFace::TBoundFace(const TBoundComp3D *bdcomp, double *param1, double *param2)
 : BoundaryJoint()
{
  int i;
  ID = BoundaryFace;

  BoundComp = bdcomp;
  for(i=0;i<4;i++)
  {
    Param1[i] = param1[i];
    Param2[i] = param2[i];
  }
}

TBoundFace::TBoundFace(const TBoundComp3D *bdcomp) : BoundaryJoint()
{
  ID = BoundaryFace;

  BoundComp = bdcomp;
  Param1[0] = 0.0;
  Param2[0] = 0.0;
  Param1[1] = 1.0;
  Param2[1] = 0.0;
  Param1[2] = 1.0;
  Param2[2] = 1.0;
  Param1[3] = 0.0;
  Param2[3] = 1.0;
}

// create a new instance of this class
TJoint *TBoundFace::NewInst(double, double, TBaseCell *)
{
  return new TBoundFace(BoundComp);
}

TJoint *TBoundFace::NewInst()
{
  return new TBoundFace(BoundComp);
}

void TBoundFace::SetParameters(double *param1, double *param2)
{
  int i;

  for(i=0;i<4;i++)
  {
    Param1[i] = param1[i];
    Param2[i] = param2[i];
  }
}

void TBoundFace::GetParameters(double *param1, double *param2) const
{
  int i;

  for(i=0;i<4;i++)
  {
    param1[i] = Param1[i];
    param2[i] = Param2[i];
  }
}
