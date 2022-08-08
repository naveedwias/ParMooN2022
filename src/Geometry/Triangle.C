// more explicit path due to name conflict with the 'triangle' library on mac
#include "../Geometry/Triangle.h"
#include "Vertex.h"
#include <cmath>

// Constructor
TTriangle::TTriangle()
{
  constexpr int TRIMAXN_EpV = 2;
  static const int DatEdgeVertex[][2] = { {0, 1},  {1, 2},  {2, 0}};
  static const int DatVertexEdge[][TRIMAXN_EpV] = { {2, 0},  {0, 1},  {1, 2}};

  MaxN_EpV = TRIMAXN_EpV;

  EdgeVertex = (const int *) DatEdgeVertex;
  VertexEdge = (const int *) DatVertexEdge;

  Type = Triangle;
  N_Vertices = 3;
  N_Edges = 3;
  N_Joints = 3;
}

// Methods
double TTriangle::GetDiameter(TVertex **Verts) const
{
  double diffX1 = Verts[0]->GetX() - Verts[1]->GetX();
  double diffY1 = Verts[0]->GetY() - Verts[1]->GetY();
  double diffZ1 = Verts[0]->GetZ() - Verts[1]->GetZ();
  double diffX2 = Verts[1]->GetX() - Verts[2]->GetX();
  double diffY2 = Verts[1]->GetY() - Verts[2]->GetY();
  double diffZ2 = Verts[1]->GetZ() - Verts[2]->GetZ();
  double diffX3 = Verts[2]->GetX() - Verts[0]->GetX();
  double diffY3 = Verts[2]->GetY() - Verts[0]->GetY();
  double diffZ3 = Verts[2]->GetZ() - Verts[0]->GetZ();

  return std::sqrt(std::max(diffX1*diffX1 + diffY1*diffY1 + diffZ1*diffZ1,
                   std::max(diffX2*diffX2 + diffY2*diffY2 + diffZ2*diffZ2,
                            diffX3*diffX3 + diffY3*diffY3 + diffZ3*diffZ3)));
}

double TTriangle::GetShortestEdge(TVertex **Verts) const
{
  double diffX1 = Verts[0]->GetX() - Verts[1]->GetX();
  double diffY1 = Verts[0]->GetY() - Verts[1]->GetY();
  double diffZ1 = Verts[0]->GetZ() - Verts[1]->GetZ();
  double diffX2 = Verts[1]->GetX() - Verts[2]->GetX();
  double diffY2 = Verts[1]->GetY() - Verts[2]->GetY();
  double diffZ2 = Verts[1]->GetZ() - Verts[2]->GetZ();
  double diffX3 = Verts[2]->GetX() - Verts[0]->GetX();
  double diffY3 = Verts[2]->GetY() - Verts[0]->GetY();
  double diffZ3 = Verts[2]->GetZ() - Verts[0]->GetZ();
  double len = 1e10;
  
  if (std::sqrt(diffX1*diffX1 + diffY1*diffY1 + diffZ1*diffZ1)< len)
      len = std::sqrt(diffX1*diffX1 + diffY1*diffY1 + diffZ1*diffZ1);
  if (std::sqrt(diffX2*diffX2 + diffY2*diffY2 + diffZ2*diffZ2)< len)
      len = std::sqrt(diffX2*diffX2 + diffY2*diffY2 + diffZ2*diffZ2);
  if (std::sqrt(diffX3*diffX3 + diffY3*diffY3 + diffZ3*diffZ3)< len)
      len = std::sqrt(diffX3*diffX3 + diffY3*diffY3 + diffZ3*diffZ3);

  return(len);
}

double TTriangle::GetLengthWithReferenceMap(TVertex **Verts) const
{
  double x0, x1, x2, y0, y1, y2;
  double xc1, xc2, yc1, yc2;

  x0 = Verts[0]->GetX();
  y0 = Verts[0]->GetY();
  x1 = Verts[1]->GetX();
  y1 = Verts[1]->GetY();
  x2 = Verts[2]->GetX();
  y2 = Verts[2]->GetY();

  //double xc0=x0;
  xc1=x1-x0;
  xc2=x2-x0;

  //double yc0=y0;
  yc1=y1-y0;
  yc2=y2-y0;

  return std::sqrt(std::abs(xc1*yc2-xc2*yc1));
}

double TTriangle::GetMeasure(TVertex **Verts) const
{
  double x1,x2,x3,y1,y2,y3;

  x1 = Verts[0]->GetX();
  y1 = Verts[0]->GetY();
  x2 = Verts[1]->GetX();
  y2 = Verts[1]->GetY();
  x3 = Verts[2]->GetX();
  y3 = Verts[2]->GetY();

  return 0.5*std::abs(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1);
}
