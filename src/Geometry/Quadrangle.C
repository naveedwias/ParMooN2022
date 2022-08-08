#include "Quadrangle.h"
#include "Vertex.h"
#include <cmath>


// Constructor
TQuadrangle::TQuadrangle()
{
  constexpr int QUADMAXN_EpV = 2;
  static const int DatEdgeVertex[][2] = { {0, 1},  {1, 2},  {2, 3},  {3, 0}};
  static const int DatVertexEdge[][QUADMAXN_EpV] =
                 { {3, 0},  {0, 1},  {1, 2},  {2, 3}};

  MaxN_EpV = QUADMAXN_EpV;

  EdgeVertex = (const int *) DatEdgeVertex;
  VertexEdge = (const int *) DatVertexEdge;

  Type = Quadrangle;
  N_Vertices = 4;
  N_Edges = 4;
  N_Joints = 4;
}

// Methods
double TQuadrangle::GetDiameter(TVertex **Verts) const
{
  double diffX1 = Verts[0]->GetX() - Verts[2]->GetX();
  double diffY1 = Verts[0]->GetY() - Verts[2]->GetY();
  double diffZ1 = Verts[0]->GetZ() - Verts[2]->GetZ();
  double diffX2 = Verts[1]->GetX() - Verts[3]->GetX();
  double diffY2 = Verts[1]->GetY() - Verts[3]->GetY();
  double diffZ2 = Verts[1]->GetZ() - Verts[3]->GetZ();

  return std::sqrt(std::max(diffX1*diffX1 + diffY1*diffY1 + diffZ1*diffZ1,
                            diffX2*diffX2 + diffY2*diffY2 + diffZ2*diffZ2));
}

double TQuadrangle::GetShortestEdge(TVertex **Verts) const
{
  double diffX1 = Verts[0]->GetX() - Verts[1]->GetX();
  double diffY1 = Verts[0]->GetY() - Verts[1]->GetY();
  double diffX2 = Verts[1]->GetX() - Verts[2]->GetX();
  double diffY2 = Verts[1]->GetY() - Verts[2]->GetY();
  double diffX3 = Verts[2]->GetX() - Verts[3]->GetX();
  double diffY3 = Verts[2]->GetY() - Verts[3]->GetY();
  double diffX4 = Verts[3]->GetX() - Verts[0]->GetX();
  double diffY4 = Verts[3]->GetY() - Verts[0]->GetY();
  double len;

  if (diffX1*diffX1 + diffY1*diffY1<diffX2*diffX2 + diffY2*diffY2)
      len = diffX1*diffX1 + diffY1*diffY1;
  else
      len = diffX2*diffX2 + diffY2*diffY2;
  if (diffX3*diffX3 + diffY3*diffY3<len)
      len = diffX3*diffX3 + diffY3*diffY3;
  if (diffX4*diffX4 + diffY4*diffY4<len)
      len = diffX4*diffX4 + diffY4*diffY4;

  return std::sqrt(len);
}

// approximation with affine map of vertices 0,1 and 3
double TQuadrangle::GetLengthWithReferenceMap(TVertex **Verts) const
{
    double x0, x1, x3, y0, y1, y3;
    double xc1, xc2, yc1, yc2;
    double detjk, rec_detjk;
    double d11, d12, d21, d22;

    x0 = Verts[0]->GetX();
    y0 = Verts[0]->GetY();
    x1 = Verts[1]->GetX();
    y1 = Verts[1]->GetY();
    x3 = Verts[3]->GetX();
    y3 = Verts[3]->GetY();
    

    //double xc0 = (x1 + x3) * 0.5;
    xc1 = (x1 - x0) * 0.5;
    xc2 = (x3 - x0) * 0.5;
    
    //double yc0 = (y1 + y3) * 0.5;
    yc1 = (y1 - y0) * 0.5;
    yc2 = (y3 - y0) * 0.5;
    
    detjk=xc1*yc2-xc2*yc1;
    rec_detjk=1/detjk;
    
    d11 = (y3-y0) * 0.5 * rec_detjk;  //dxi/dx
    d12 = (x0-x3) * 0.5 * rec_detjk;  //dxi/dy
    d21 = (y0-y1) * 0.5 * rec_detjk;  //deta/dx
    d22 = (x1-x0) * 0.5 * rec_detjk;  //deta/dy

    Output::print(d11, " ", d12, " ", d21, " ", d22);
    ErrThrow("This gives the same result as GetMeasure()");

    // factor 2 = std::sqrt(4.) because of area of unit square
    return 2*std::sqrt(std::abs(detjk));
}

double TQuadrangle::GetMeasure(TVertex **Verts) const
{
  double x1,x2,x3,x4,y1,y2,y3,y4;

  x1 = Verts[0]->GetX();
  y1 = Verts[0]->GetY();
  x2 = Verts[1]->GetX();
  y2 = Verts[1]->GetY();
  x3 = Verts[2]->GetX();
  y3 = Verts[2]->GetY();
  x4 = Verts[3]->GetX();
  y4 = Verts[3]->GetY();

  return 0.5*std::abs(x2*y3-x3*y2-x1*y3+x3*y1+x1*y2-x2*y1)+
    0.5*std::abs(x4*y3-x3*y4-x1*y3+x3*y1+x1*y4-x4*y1);
}

Shapes TQuadrangle::CheckQuad(const TVertex * const * Vertices) const
{
  double test1, test2;

  test1 = Vertices[0]->GetX() - Vertices[1]->GetX() +
          Vertices[2]->GetX() - Vertices[3]->GetX();
  test2 = Vertices[0]->GetY() - Vertices[1]->GetY() +
          Vertices[2]->GetY() - Vertices[3]->GetY();

  if (std::abs(test1) > 1e-8 || std::abs(test2) > 1e-8) return Quadrangle;

// check for rectangle
//  test1 = (Vertices[0]->GetX() - Vertices[1]->GetX()) * 
//          (Vertices[2]->GetX() - Vertices[1]->GetX()) +
//          (Vertices[0]->GetY() - Vertices[1]->GetY()) * 
//          (Vertices[2]->GetY() - Vertices[1]->GetY());
//         
//  test2 = (Vertices[0]->GetX() - Vertices[2]->GetX()) * 
//          (Vertices[3]->GetX() - Vertices[1]->GetX()) +
//          (Vertices[0]->GetY() - Vertices[2]->GetY()) * 
//          (Vertices[3]->GetY() - Vertices[1]->GetY());
//         
//  if (std::abs(test1) > 1e-8 || std::abs(test2) > 1e-8) return Parallelogram;

  test1 = (Vertices[3]->GetX() - Vertices[0]->GetX()) *
          (Vertices[1]->GetX() - Vertices[0]->GetX()) +
          (Vertices[3]->GetY() - Vertices[0]->GetY()) *
          (Vertices[1]->GetY() - Vertices[0]->GetY());

  if (std::abs(test1) > 1e-8) return Parallelogram;
  return Rectangle;
}
