// =======================================================================
// @(#)BdLine.C        1.2 07/16/99
//
// Class:       TBdLine
// Purpose:     a line as a component of a boundary part
//
// Author:      Volker Behns  18.06.97
//
// =======================================================================
#include <cmath>

#include <BdLine.h>
#include <MooNMD_Io.h>

// Constructor
TBdLine::TBdLine(int id) : TBoundComp2D(id)
{
  Type = Line;
}

// Methods
void TBdLine::SetParams (double xstart, double ystart,
                         double delx, double dely)
{
  Xstart = xstart;
  Ystart = ystart;
  delX = delx;
  delY = dely;
}

int TBdLine::GetXYofT(double T, double &X, double &Y) const
{
  if (T < -0.001 || T > 1.001)
  { 
    cerr << "Error in TBdLine: parameter T out of range" << endl;
    cerr << T << endl;
    return -2;
  }

  X = Xstart + T * delX;
  Y = Ystart + T * delY;

  return 0;
}

int TBdLine::GetTofXY(double X, double Y, double &T) const
{
  double T_X = -123, T_Y = -123;
  bool testX = false, testY = false;

  if (std::abs(delX) > 0.0001)
    T = T_X = (X - Xstart) / delX;
  else
    if (std::abs(X - Xstart) < 0.0001) testX = true;
  
  if (std::abs(delY) > 0.0001)
    T = T_Y = (Y - Ystart) / delY;
  else
    if (std::abs(Y - Ystart) < 0.0001) testY = true;

  if (T < -0.001 || T > 1.001 || (!(testX || testY) && std::abs(T_X - T_Y) > 1e-6))
    return -1;
  else
    return 0;
}

int TBdLine::ReadIn(std::istream &dat)
{
  char line[100];

  dat >> Xstart >> Ystart;
  dat.getline (line, 99);
  dat >> delX >> delY;
  dat.getline (line, 99);

  return 0;
}
