// =======================================================================
// @(#)BdWall.C        1.2 07/16/99
//
// Class:       TBdWall
// Purpose:     a sandwich side wall as a component of a boundary part
//
// Author:      Volker Behns  05.07.99
//
// =======================================================================

#include <BdWall.h>
#include <MooNMD_Io.h>

// Constructor
TBdWall::TBdWall(int id, TBoundComp2D *bdcomp2d) : TBoundComp3D(id)
{
  Type = Wall;

  BdComp2D = bdcomp2d;
}

// Methods
int TBdWall::SetParams(double drx, double dry, double drz, double conic_sc)
{
  DriftX = drx;
  DriftY = dry;
  DriftZ = drz;
  
  ConicScale = conic_sc;

  return 0;
}

int TBdWall::GetXYZofTS(double T, double S,
                        double &X, double &Y, double &Z) const
{
  double x,y; //z;

  BdComp2D->GetXYofT(T, x, y);

  X = x + S * DriftX + x * S * ConicScale;
  Y = y + S * DriftY + y * S * ConicScale;
  Z = S * DriftZ;

  return 0;
}

/** return parameters and coordinates of a given linear
    combination of vertices */
int TBdWall::GetXYZandTS(int N_Points, double *LinComb,
                         double *, double *, double *,
                         double *tp, double *sp,
                         double &X, double &Y, double &Z,
                         double &T, double &S) const
{
  int i;
  double t, s, v;
  double x,y; //z;

  t = s = 0;
  for(i=0;i<N_Points;i++)
  {
    v = LinComb[i];
    t += v*tp[i];
    s += v*sp[i];
  }
  T = t;
  S = s;

  BdComp2D->GetXYofT(T, x, y);

  X = x + S * DriftX + x * S * ConicScale;
  Y = y + S * DriftY + y * S * ConicScale;
  Z = S * DriftZ;

  return 0;
}

int TBdWall::GetTSofXYZ(double X, double Y, double Z,
                        double &T, double &S) const
{
  double x,y,t;

  S = Z / DriftZ;

  x = (X - S * DriftX) / (1+ S * ConicScale);
  y = (Y - S * DriftY) / (1+ S * ConicScale);

  BdComp2D->GetTofXY(x, y, t);

  T = t;
  
  return 0;
}

int TBdWall::ReadIn(std::istream &dat)
{
  return BdComp2D->ReadIn(dat);
}
