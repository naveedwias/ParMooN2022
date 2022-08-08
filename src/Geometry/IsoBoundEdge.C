// =======================================================================
// @(#)IsoBoundEdge.C        1.2 09/13/99
// 
// Class:       TIsoBoundEdge
// Purpose:     edge on a boundary component with additional vertices
//              for isoparametric reference transformation
//
// Author:      Gunar Matthies  06.08.99
//
// =======================================================================

#include <BoundComp2D.h>
#include <IsoBoundEdge.h>
#include <Vertex.h>

// Constructors
TIsoBoundEdge::TIsoBoundEdge(const TBoundComp2D *bdcomp, double t_0, double t_1)
  : TBoundEdge(bdcomp, t_0, t_1)
{
  ID = IsoBoundEdge;

  N_Vertices = 0;
  Vertices = nullptr;
}

// Methods
int TIsoBoundEdge::CheckMatchingRef(TBaseCell *, int, struct StoreGeom &Tmp)
{
  Tmp.Filled = false;
  return 0;
}

// create a new instance of this class
TJoint *TIsoBoundEdge::NewInst(double newT_0, double newT_1, TBaseCell *)
{
  return new TIsoBoundEdge(BoundComp, T_0 + newT_0*(T_1 - T_0),
                        T_0 + newT_1*(T_1 - T_0));
}

TJoint *TIsoBoundEdge::NewInst()
{
  return new TIsoBoundEdge(BoundComp, T_0, T_1);
}

void TIsoBoundEdge::SetVertices(int n_vertices, TVertex **vertices)
{
  if(Vertices)
    delete Vertices;

  N_Vertices = n_vertices;
  Vertices = vertices;
}

void TIsoBoundEdge::GenerateVertices(int n_vertices)
{
  if(N_Vertices != n_vertices)
  {
    if(Vertices)
      delete Vertices;

    N_Vertices = n_vertices;
    Vertices = new TVertex*[N_Vertices];
  }
#ifdef __2D__
  double t, x, y;
  for(int i=0;i<N_Vertices;i++)
  {
    t = T_0 + ((i+1)*(T_1-T_0))/(N_Vertices+1);
    BoundComp->GetXYofT(t, x, y);
    Vertices[i] = new TVertex(x, y);
  } // endfor i
#endif // __2D__
}

void TIsoBoundEdge::GeneratemidVert(int n_vertices, double* X, double* Y)
{
   if(Vertices)
     delete Vertices;

  N_Vertices = n_vertices;
  Vertices = new TVertex*[N_Vertices];

  for(int i=0;i<N_Vertices;i++)
  {
#ifdef __2D__
    double x = (X[0]+X[1])/2.0;
    double y = (Y[0]+Y[1])/2.0;
    Vertices[i] = new TVertex(x, y);
#else
    (void)X; // avoid compiler warning
    (void)Y; // avoid compiler warning
#endif
    
    //double r = std::sqrt(x*x + y*y);
    
    //cout << r<< " x " << x << " y " << y <<endl;
    //x /=r;
    //y /=r;
    //cout << " x " << x << " y " << y <<endl;
  } // endfor i
}

void TIsoBoundEdge::DeleteVertices()
{
  if(Vertices)
    delete Vertices;
  N_Vertices = 0;
  Vertices = nullptr;
 }
 
 

 
