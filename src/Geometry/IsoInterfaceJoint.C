// =======================================================================
// @(#)IsoInterfaceJoint.C        1.2 09/13/99
// 
// Class:       TIsoInterfaceJoint
// Purpose:     connects two cells on an interface with additional
//              vertices for isoparametric reference transformation
//
// Author:      Gunar Matthies  06.08.99
//
// =======================================================================

#include <IsoInterfaceJoint.h>
#include <BoundComp2D.h>
#include <BaseCell.h>

// Constructors
TIsoInterfaceJoint::TIsoInterfaceJoint(const TBoundComp2D *bdcomp, double t_0,
                 double t_1, TBaseCell *neighb0) 
 : TInterfaceJoint(bdcomp, t_0, t_1, neighb0) 
{
  ID = IsoInterfaceJoint;

  N_Vertices = 0;
  Vertices = nullptr;
}

TIsoInterfaceJoint::TIsoInterfaceJoint(const TBoundComp2D *bdcomp, double t_0,
                 double t_1, TBaseCell *neighb0, TBaseCell *neighb1)
 : TInterfaceJoint(bdcomp, t_0, t_1, neighb0, neighb1) 
{
  ID = IsoInterfaceJoint;

  N_Vertices = 0;
  Vertices = nullptr;
}

// Methods
TJoint *TIsoInterfaceJoint::NewInst(double newT_0, double newT_1, TBaseCell *Me)
{
  return new TIsoInterfaceJoint(BoundComp, T_0 + newT_0*(T_1 - T_0),
                             T_0 + newT_1*(T_1 - T_0), Me);
}

TJoint *TIsoInterfaceJoint::NewInst()
{
  return new TIsoInterfaceJoint(BoundComp, T_0, T_1, nullptr);
}

void TIsoInterfaceJoint::SetVertices(int n_vertices, TVertex **vertices)
{
  delete Vertices;

  N_Vertices = n_vertices;
  Vertices = vertices;
}

void TIsoInterfaceJoint::GenerateVertices(int n_vertices)
{
  if(N_Vertices != n_vertices)
  {
    if(Vertices)
      delete Vertices;

    N_Vertices = n_vertices;
    Vertices = new TVertex*[N_Vertices];
  }

#ifdef __2D__
    int i;
    TVertex *vertex;
    double t, x, y;
  for(i=0;i<N_Vertices;i++)
  {
    t = T_0 + ((i+1)*(T_1-T_0))/(N_Vertices+1);
    
    BoundComp->GetXYofT(t, x, y);

    vertex = new TVertex(x, y);

    Vertices[i] = vertex;
  } // endfor i
#endif // __2D__
}

void TIsoInterfaceJoint::GeneratemidVert(int n_vertices, double* X, double* Y)
{
  if(Vertices)
    delete Vertices;

  N_Vertices = n_vertices;
  Vertices = new TVertex*[N_Vertices];

  int i;
  for(i=0;i<N_Vertices;i++)
  {
#ifdef __2D__
   double x = (X[0]+X[1])/2.0;
   double y = (Y[0]+Y[1])/2.0;
   Vertices[i] = new TVertex(x, y);
//     cout <<N_Vertices << "test ("<< x <<", "<< y << ")"<<endl;
#else
    (void)X; // avoid compiler warning
    (void)Y; // avoid compiler warning
#endif // __2D__
  } // endfor i
}

