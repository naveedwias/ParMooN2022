// =======================================================================
// @(#)InterfaceJoint.C        1.2 01/26/12
// 
// Class:       TInnerInterfaceJoint
// Purpose:     connects two cells of different collection
//              (not treated as a boundary joint)
//
//
// =======================================================================


#include <InnerInterfaceJoint.h>
#include <BoundComp2D.h>
#include <BaseCell.h>
#include <cmath>

// Constructors

TInnerInterfaceJoint::TInnerInterfaceJoint(
    TBaseCell *neighb0, TBaseCell *neighb1) :
  TJointEqN(neighb0, neighb1)
{
  ID = InnerInterfaceJoint;
  NeibGlobalCellNo[0]=neighb0->GetCellIndex();
  NeibGlobalCellNo[1]=neighb1->GetCellIndex();
  SubDomainsID[0]=neighb0->GetReference_ID();
  SubDomainsID[1]=neighb1->GetReference_ID();
  children[0] = nullptr;
  children[1] = nullptr;
  //for 3D
  children[2] = nullptr;
  children[3] = nullptr;
  
  //for 2D starting point and vector pointing to other point of the edge
  Xstart = -4711;
  Ystart = -4711;
  delX = -4711;
  delY = -4711;
  
  //for 3D
  //number of vertices per face
  N_vert = -4711;
  //array of length N_vert containing x-,y- and z-component
  X = nullptr;
  Y = nullptr;
  Z = nullptr;
  //which local face of the cell is this interface
  //which neighboring cell is used depends on the ID given
  //in GetInnerInterface in auxiliaryFunctions.C
  local_face = -4711;

  dim = -4711;
}

TInnerInterfaceJoint::TInnerInterfaceJoint(TBaseCell *neighb0): 
  TJointEqN(neighb0)
{
  ID = InnerInterfaceJoint;
  NeibGlobalCellNo[0]=neighb0->GetCellIndex();
  NeibGlobalCellNo[1]=-1;
  SubDomainsID[0]=neighb0->GetReference_ID();
  SubDomainsID[1]=-1;
  children[0] = nullptr;
  children[1] = nullptr;
  //for 3D
  children[2] = nullptr;
  children[3] = nullptr;
  
  //for 2D starting point and vector pointing to other point of the edge
  Xstart = -4711;
  Ystart = -4711;
  delX = -4711;
  delY = -4711;

  //for 3D
  //number of vertices per face
  N_vert = -4711;
  //array of length N_vert containing x-,y- and z-component
  X = nullptr;
  Y = nullptr;
  Z = nullptr;
  //which local face of the cell is this interface
  //which neighboring cell is used depends on the ID given
  //in GetInnerInterface in auxiliaryFunctions.C
  local_face = -4711;

  dim = -4711;
}


// Methods
TJoint *TInnerInterfaceJoint::NewInst(double, double, TBaseCell *Me)
{
  return new TInnerInterfaceJoint(Me);
}

TJoint *TInnerInterfaceJoint::NewInst()
{
  return new TInnerInterfaceJoint(nullptr);
}

/** Remember which joints are children of this joint */
void TInnerInterfaceJoint::SetChild(TInnerInterfaceJoint *child)
{
  if(children[0] == nullptr)
    children[0] = child;
  else if(children[1] == nullptr)
    children[1] = child;
  else if(children[2] == nullptr && dim==3)
      children[2] = child;
  else if(children[3] == nullptr && dim==3)
      children[3] = child;
  // if two (in 2D) or four (in3D) children are already set, then nothing
  // happens,
  if((children[0] != child && children[1] != child &&
      children[2] != child && children[3] != child && dim==3) ||
     (children[0] != child && children[1] != child && dim==2))
  {
    ErrThrow("TInnerInterfaceJoint::SetChild");
  }
}

/** return one of the two children of this edge */
TInnerInterfaceJoint *TInnerInterfaceJoint::GetChild(int child) const
{
  if(child == 0 || child == 1)
    return children[child];
  else if((child == 2 || child == 3) && dim == 3)
    return children[child];
  else
  {
    ErrThrow("ERROR, TInnerInterfaceJoint::GetChild, no such child available");
  }
}

/** Get either one of the two neighbors */
TBaseCell *TInnerInterfaceJoint::GetNeighbour(int i) const
{
  if (i==1)
    return Neighb1;
  else if(i==0)
    return Neighb0;
  else
  {
    ErrThrow("TInnerInterfaceJoint::GetNeighbour: only two neighbors ", i);
  }
}

/** set the coordinates of the start point and the vector pointing to the     
    second point */
void TInnerInterfaceJoint::SetParams(double xstart, double ystart, 
                                     double delx, double dely)
{
  Xstart = xstart;
  Ystart = ystart;
  delX = delx;
  delY = dely;
  dim = 2;
  N_vert = 2;
}
/** get the coordinates of the start point and the vector pointing to the     
    second point */
void TInnerInterfaceJoint::GetParams(double &xstart, double &ystart, 
                                     double &delx, double &dely) const
{
  xstart = Xstart;
  ystart = Ystart;
  delx = delX;
  dely = delY;
}

void TInnerInterfaceJoint::SetParams(double *x, double *y, double *z, 
                                     int n_vert)
{
  N_vert = n_vert;

  X = new double[n_vert];
  Y = new double[n_vert];
  Z = new double[n_vert];

  for(unsigned int i_vert=0; i_vert < N_vert; i_vert++){
    X[i_vert] = x[i_vert];
    Y[i_vert] = y[i_vert];
    Z[i_vert] = z[i_vert];
  }

  dim = 3;
}

void TInnerInterfaceJoint::GetParams(double *x, double *y, double *z) const
{
  for(unsigned int i_vert=0; i_vert < N_vert; i_vert++){
    x[i_vert] = X[i_vert];
    y[i_vert] = Y[i_vert];
    z[i_vert] = Z[i_vert];
  }
}

double TInnerInterfaceJoint::GetLength() const
{
  return std::sqrt(delX*delX + delY*delY);
}

double TInnerInterfaceJoint::GetArea() const
{
  if(N_vert == 3){
    //using A=0.5*|X x Y|
    double x1 = X[1]-X[0];
    double x2 = Y[1]-Y[0];
    double x3 = Z[1]-Z[0];
    double y1 = X[2]-X[0];
    double y2 = Y[2]-Y[0];
    double y3 = Z[2]-Z[0];
    return 0.5*std::sqrt(std::pow(x2*y3-x3*y2,2) + std::pow(x3*y1-x1*y3,2) + std::pow(x1*y2-x2*y1,2));
  }
  if(N_vert == 4){
  //using A=0.5*std::sqrt(|e|^2*|f|^2-(e,f)^2)
  double norm_e = std::pow((X[0]-X[2]),2)
                 +std::pow((Y[0]-Y[2]),2)
                 +std::pow((Z[0]-Z[2]),2);
  double norm_f = std::pow((X[1]-X[3]),2)
                 +std::pow((Y[1]-Y[3]),2)
                 +std::pow((Z[1]-Z[3]),2);
  double ef = (X[0]-X[2])*(X[1]-X[3])
             +(Y[0]-Y[2])*(Y[1]-Y[3])
             +(Z[0]-Z[2])*(Z[1]-Z[3]);
  return 0.5*std::sqrt(norm_e*norm_f-ef*ef);
  }
  ErrThrow("N_vert = ", N_vert);
}

/** the unit normal of this edge */
void TInnerInterfaceJoint::GetNormal(double &nx, double &ny) const
{
  double length = std::sqrt(delX*delX + delY*delY);
  nx =  delY/length;
  ny = -delX/length;
}

/** the unit normal of this face */
void TInnerInterfaceJoint::GetNormal(double &nx, double &ny, double &nz) const
{
  nx =  ((Y[1]-Y[0])*(Z[N_vert-1]-Z[1]) - (Z[1]-Z[0])*(Y[N_vert-1]-Y[1]));
  ny =  ((Z[1]-Z[0])*(X[N_vert-1]-X[1]) - (X[1]-X[0])*(Z[N_vert-1]-Z[1]));
  nz =  ((X[1]-X[0])*(Y[N_vert-1]-Y[1]) - (Y[1]-Y[0])*(X[N_vert-1]-X[1]));
  double length = std::sqrt(nx*nx + ny*ny + nz*nz);
  nx = nx/length;
  ny = ny/length;
  nz = nz/length;
}

/** the unit tangential of this edge */
void TInnerInterfaceJoint::GetTangent(double &tx, double &ty) const
{
  double length = std::sqrt(delX*delX + delY*delY);
  tx = delX/length;
  ty = delY/length;
}
/** the unit tangentials of this face */
void TInnerInterfaceJoint::GetTangent(double &tx, double &ty, double &tz,
                                      double &tx2, double &ty2, double &tz2) const
{
  // first tangent
  double length = std::sqrt((X[1]-X[0])*(X[1]-X[0])
                     + (Y[1]-Y[0])*(Y[1]-Y[0])
                     + (Z[1]-Z[0])*(Z[1]-Z[0]));
  tx = (X[1]-X[0])/length;
  ty = (Y[1]-Y[0])/length;
  tz = (Z[1]-Z[0])/length;
  // second tangent computed via cross-product with normal
  double nx,ny,nz;
  this->GetNormal(nx,ny,nz);
  length = std::sqrt((ty*nz - tz*ny)*(ty*nz - tz*ny)
              + (tz*nx - tx*nz)*(tz*nx - tx*nz)
              + (tx*ny - ty*nx)*(tx*ny - ty*nx));
  tx2 = (ty*nz - tz*ny)/length;
  ty2 = (tz*nx - tx*nz)/length;
  tz2 = (tx*ny - ty*nx)/length;
}

/** get the index of this joint in given neighbor */
int TInnerInterfaceJoint::GetIndexInNeighbor(const TBaseCell*const neigh) const
{
  if(neigh == Neighb0)
    return IndexInNeighbor[0];
  else if(neigh == Neighb1)
    return IndexInNeighbor[1];
  ErrThrow("TInnerInterfaceJoint::GetIndexInNeighbor !!!!!!!!");
}

/** set the index of this joint in given neighbor */
void TInnerInterfaceJoint::SetIndexInNeighbor(TBaseCell *neigh, int index)
{
  if(neigh == Neighb0)
    IndexInNeighbor[0] = index;
  else if(neigh == Neighb1)
    IndexInNeighbor[1] = index;
  else
  {
    ErrThrow("TInnerInterfaceJoint::SetIndexInNeighbor !!!!!!!!");
  }
}

void TInnerInterfaceJoint::PrintInfo() const
{
  if(dim == 3)
  {
    Output::print("This InnerInterface has ", N_vert, " vertices");
    for(unsigned int iVert = 0; iVert < N_vert; iVert++)
    {
      Output::print("Vertex ", iVert, ": (", X[iVert], ",", Y[iVert], ",",
                    Z[iVert], ")");
    }
    Output::print("Index of this face in the first neighbor: ",
                  IndexInNeighbor[0]);
    Output::print("Index of this face in the second neighbor: ",
                  IndexInNeighbor[1]);
    Output::print("Local face: ", local_face);
  }
  else if(dim == 2)
  {
    Output::print("This InnerInterface has ", N_vert, " vertices");
    Output::print("Vertex 1: (", Xstart, ",", Ystart, ")");
    Output::print("Vertex 2: (", Xstart+delX, ",", Ystart+delY, ")");
    Output::print("Index of this face in the first neighbor: ",
                  IndexInNeighbor[0]);
    Output::print("Index of this face in the second neighbor: ",
                  IndexInNeighbor[1]);
    Output::print("Local face: ", local_face);
  }
  else
    Output::print("dim = ", dim, ". dimension should be 2D or 3D");
}
