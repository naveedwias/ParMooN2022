// =======================================================================
// @(#)BdSphere.C        1.2 07/16/99
//
// Class:       TBdSphere
// Purpose:     a Sphere as a component of a boundary part
//
// Author:      Gunar Matthies 2000/12/04
//
// =======================================================================

#include <BdSphere.h>
#include <cmath>

// Constructor
TBdSphere::TBdSphere(int id) : TBoundComp3D(id)
{
  Type = Sphere;
}

// Methods
void TBdSphere::SetParams (double m_x, double m_y, double m_z,
                         double r)
{
  M_x = m_x;
  M_y = m_y;
  M_z = m_z;

  R = r;
}

int TBdSphere::GetXYZofTS(double T, double S,
                        double &X, double &Y, double &Z) const
{
  X = M_x + R * std::cos(S)*std::cos(T);
  Y = M_y + R * std::cos(S)*std::sin(T);
  Z = M_z + R * std::sin(S);

  return 0;
}

/** return parameters and coordinates of a given linear
    combination of vertices */
int TBdSphere::GetXYZandTS(int N_Points, double *LinComb,
                          double *xp, double *yp, double *zp,
                          double *, double *,
                          double &X, double &Y, double &Z,
                          double &T, double &S) const
{
  double v;
  int i;

  X = 0; Y = 0; Z = 0;
  for(i=0;i<N_Points;i++)
  {
    v = LinComb[i];
    X += v*xp[i];
    Y += v*yp[i];
    Z += v*zp[i];
  }
  X -= M_x;
  Y -= M_y;
  Z -= M_z,

  v = std::sqrt(X*X+Y*Y+Z*Z)/R;
  
  X /= v;
  Y /= v;
  Z /= v;

  S = std::asin(Z/R);
  T = std::atan2(Y, X);

  X += M_x;
  Y += M_y,
  Z += M_z;

  return 0;
}

int TBdSphere::GetTSofXYZ(double X, double Y, double Z,
                        double &T, double &S) const
{
  S = std::asin((Z-M_z)/R);
  T = std::atan2((Y-M_y), (X-M_x));

  return 0;
}

int TBdSphere::ReadIn(std::istream &dat)
{
  char line[100];
  double a, b;

  dat >> M_x >> M_y >> M_z;
  dat.getline (line, 99);
  dat >> R >> a >> b;
  dat.getline (line, 99);

  return 0;
}
