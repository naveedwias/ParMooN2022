#include "Line.h"
#include "Vertex.h"
#include <cmath>

// Constructor
TLine::TLine()
{
  static const int DatEdgeVertex[][2] = { {0, 1}};
  static const int DatVertexEdge[][LINEMAXN_EpV] = { {0},  {0}};

  MaxN_EpV = LINEMAXN_EpV;

  EdgeVertex = (const int *) DatEdgeVertex;
  VertexEdge = (const int *) DatVertexEdge;

  Type = S_Line;
  N_Vertices = 2;
  N_Edges = 2;
  N_Joints = 2;
}

// Methods
double TLine::GetMeasure(TVertex **Verts) const
{
  double x1,x2,y1,y2,z1,z2;

  x1 = Verts[0]->GetX();
  y1 = Verts[0]->GetY();
  z1 = Verts[0]->GetZ();
  x2 = Verts[1]->GetX();
  y2 = Verts[1]->GetY();
  z2 = Verts[1]->GetZ();

  return std::sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
}
