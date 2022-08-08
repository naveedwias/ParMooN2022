// =======================================================================
// @(#)Joint.C        1.7 11/18/99
// 
// Class:       TJoint
// Purpose:     superclass for all joints
//
// Author:      Volker Behns  23.09.98
//
// =======================================================================

#include <BaseCell.h>
#include <Database.h>
#include <Joint.h>
#include <stdlib.h>
#include <Line.h>
#include <Triangle.h>
#include <Quadrangle.h>


// Constructors
TJoint::TJoint()
{
  ID = Joint;

  Neighb0 = nullptr;
  Neighb1 = nullptr;

  NeibSubDomainLocalJointNo = -1;
  IndexInNeighbor[0] = -1;
  IndexInNeighbor[1] = -1;
  MapType = -1;
}

// Methods
int TJoint::SetNeighbour(TBaseCell *Neighb)
{
    switch (this->GetType())
    {
        case JointEqN:
        case InterfaceJoint:
        case IsoInterfaceJoint:
#ifdef __3D__
        case InterfaceJoint3D:
        case IsoInterfaceJoint3D:
#endif
        case PeriodicJoint:
        case InnerInterfaceJoint:
            if (Neighb0)
                Neighb1 = Neighb;
            else
                Neighb0 = Neighb;
            
            return 0;
            
        case BoundaryEdge:
        case BoundaryFace:
        case IsoBoundFace:
            Neighb0 = Neighb;
            
            return 0;

       default:
            ErrThrow("Unknown joint type", this->GetType());
    }
}

TBaseCell *TJoint::GetNeighbour(const TBaseCell *Me) const
{
  if (Neighb0 == Me)
    return Neighb1;
  else
    return Neighb0;
}

int TJoint::SetNeighbour(int i, TBaseCell *Neighb)
{
  switch (this->GetType())
  {
    case JointEqN:
    case InterfaceJoint:
    case IsoInterfaceJoint:
#ifdef __3D__
    case InterfaceJoint3D:
    case IsoInterfaceJoint3D:
#endif
    case PeriodicJoint:
      if (i)
        Neighb1 = Neighb;
      else
        Neighb0 = Neighb;

      return 0;
          case BoundaryEdge:
          if (i==0)
              Neighb0 = Neighb;
          else
              Output::print("** Error in TJoint::SetNeighbour: cannot set Neighbour ", i, " for Boundary Edge");
          return 0;    default:
      return -1;
  }
}

TBaseCell *TJoint::GetNeighbour(int i) const
{
  if (i)
    return Neighb1;
  else
    return Neighb0;
}

void TJoint::Delete(TBaseCell *Me)
{
  if (Neighb0 == Me)
    Neighb0 = nullptr;
  else
    Neighb1 = nullptr;
}

//NEW LB 11.10.18
/**
 @brief set the local index of this joint in given neighbor
 @todo write a more general function that does not take index as input
 */
void TJoint::set_index_in_neighbour(const TBaseCell *neigh, int index)
{
  if(neigh == Neighb0)
    IndexInNeighbor[0] = index;
  else if(neigh == Neighb1)
    IndexInNeighbor[1] = index;
  else
  {
    ErrThrow("ERROR in TJoint::set_index_in_neighbour");
  }
}

int TJoint::get_index_in_neighbour(const TBaseCell*const neigh) const
{
  if(neigh == Neighb0)
    return IndexInNeighbor[0];
  else if(neigh == Neighb1)
    return IndexInNeighbor[1];
  else
  {
    ErrThrow("ERROR in TJoint::get_index_in_neighbour");
  }
}

int TJoint::get_joint_nr_in_cell(const TBaseCell* cell) const
{
  int joint_nr = 0;
  // Check if joint belongs to cell
  if(cell != Neighb0 && cell != Neighb1)
  {
    // Joint does not belong to cell
    ErrThrow("This joint does not belong to the given cell.");
  }
  else
  {
    // joint belongs to cell -> Find out the local index
    while( cell->GetJoint(joint_nr) != this )
    {
      joint_nr++;
    }
  }
  return joint_nr;
}

double TJoint::GetDiameter()
{
  // We want to compute the diameter of the joint, i.e. a positive value. We
  // initialize the diameter that is called diam_joint with a negative value to
  // check later if the computation worked.
  double diam_joint = -1;

  // The idea here is to create a TLine, TTriangle or TQuadrangle object
  // depending on the shape of the facet and then use methods of those objects
  // that determine the diameter of the facet.
  // Everything with respect to the Shape is encoded for instance in the shape
  // descriptor of a cell to which the joint belongs to, e.g. the shape
  // descriptor of Neighb0.
  if (Neighb0 == nullptr)
  {
    ErrThrow("Neighb0 not known. Call TJoint::SetNeighbour() before calling ",
        "this method, see documentation of TJoint::GetDiameter().");
  }
  auto shape = Neighb0->GetType();
  // These shape types belong to 2D shapes, see also Shapes.
  if (shape >= 1 && shape <= 4)
  {
    // In 2D, the facet aka line is defined by a starting point and an end
    // point. The information which vertices build the line are encoded in the
    // shape descriptor, see also TShapeDesc.

    // Get the shape descriptor of the cell and the EdgeVertex array
    auto shape_desc = Neighb0->GetShapeDesc();
    const int* edge_vertex;
    shape_desc->GetEdgeVertex(edge_vertex);

    // Get the vertices that are the starting and the end point of the edge. To
    // use TLine::GetMeasure() pointers to the vertices have to be stored in an
    // vertex array.
    std::vector<TVertex*> vertices_vec(2);

    // We further need the index the joint has in the cell Neighb0
    auto joint_index = this->get_joint_nr_in_cell(Neighb0);

    // The first 2 * joint_index vertices belong to the edges with local edge
    // number smaller than joint_index. Therefore, the number of the vertices
    // defining joint_index start at position 2 * joint_index in the EdgeVertex
    // array.
    int starting_index = 2 * joint_index;
    for (int vertex_i = 0; vertex_i < 2; ++vertex_i)
    {
      vertices_vec[vertex_i] = Neighb0->GetVertex( edge_vertex[
          starting_index + vertex_i] );
    }

    // Construct a line and compute the length of the line
    TLine line;
    diam_joint = line.GetMeasure(vertices_vec.data());
  }
  else if (shape >= 5 && shape <= 7)
  {
    // The shape of the cell and hence also the faces is described by the shape
    // descriptor. The idea here is to get the vertices that define the joint
    // and then use TTriangle::GetDiameter and TQuadrange::GetDiameter
    // respectively to compute the diameter of the face.

    // Get the shape descriptor and the FaceVertex array that defines the
    // vertices that belong to the joint. See TShapeDesc for more information
    // about the Shape Descriptor.
    auto shape_desc = Neighb0->GetShapeDesc();
    const int* face_vertex;
    const int* face_vertex_len;
    int max_face_vertex_len;
    shape_desc->GetFaceVertex(face_vertex, face_vertex_len,
        max_face_vertex_len);

    // We further need the index the joint has in the cell Neighb0
    auto joint_index = this->get_joint_nr_in_cell(Neighb0);

    // We later want to use TTriangle::GetDiameter and TQuadrangle::GetDiameter
    // which need the vertices of the facet. Therefore, we compute a vector of
    // pointers to the vertices that define the facet.
    std::vector<TVertex*> vertices_vec( face_vertex_len[joint_index] );

    // Find the position where the vertices of the joint_index-th face start.
    // The ordering in FaceVertex is such that first the vertices of the 0th
    // face are specified, then of the 1st face and so on. Therefore, the
    // starting point for our current face is the sum of all vertices that
    // belong to faces with a smaller local face number.
    int starting_index = 0;
    for (int facet_nr = 0; facet_nr < joint_index; ++facet_nr)
    {
      starting_index += face_vertex_len[facet_nr];
    }
    for (int vertex_i = 0; vertex_i < face_vertex_len[joint_index]; ++vertex_i)
    {
      vertices_vec[vertex_i] = Neighb0->GetVertex( face_vertex[
          starting_index + vertex_i] );
    }

    // The shape of the facet is either a triangle or a quadrilateral depending
    // on the shape of the cell. Now an object of the right shape is initialized
    // and the respective GetDiameter function is called to compute the
    // diameter.
    if (face_vertex_len[joint_index] == 3)
    {
      TTriangle triangle;
      diam_joint = triangle.GetDiameter(vertices_vec.data());
    }
    else if (face_vertex_len[joint_index] == 4)
    {
      TQuadrangle quad;
      diam_joint = quad.GetDiameter(vertices_vec.data());
    }
    else
    {
      ErrThrow("Only tetrahedral or hexahedral elements are implemented yet.");
    }
  }
  else
  {
    ErrThrow("Only TShapes 1-7 are implemented yet.");
  }

  // Check if the computation worked and if so return the diameter
  if (diam_joint < 0)
  {
    ErrThrow("Diameter of facet was not correctly calculated.");
  }
  return diam_joint;
}

#ifdef __3D__
void TJoint::SetMapType()
{
  int N_, LocJoint0, LocJoint1, MaxLen, aux;
  const int *TmpFV, *TmpLen;
  TVertex *Vert;

  if (Neighb0 && Neighb1)
  {
    N_ = Neighb0->GetN_Faces();
    for (LocJoint0=0;LocJoint0<N_;LocJoint0++)
      if (Neighb0->GetJoint(LocJoint0) == this) break;

    N_ = Neighb1->GetN_Faces();
    for (LocJoint1=0;LocJoint1<N_;LocJoint1++)
      if (Neighb1->GetJoint(LocJoint1) == this) break;

    Neighb0->GetRefDesc()->GetShapeDesc()->
             GetFaceVertex(TmpFV, TmpLen, MaxLen);

    Vert = Neighb0->GetVertex(TmpFV[LocJoint0 * MaxLen]);

    Neighb1->GetRefDesc()->GetShapeDesc()->
             GetFaceVertex(TmpFV, TmpLen, MaxLen);

    N_ = TmpLen[LocJoint1];
    aux = LocJoint1 * MaxLen;

    for (MapType=0;MapType<N_;MapType++)
      if (Neighb1->GetVertex(TmpFV[aux + MapType]) == Vert) break;

    if (MapType == N_)
    {
      /*
      int i;
      N_ = Neighb0->GetN_Vertices();
      for (i=0;i<N_;i++)
        cout << " test 0:" << i << Neighb0->GetVertex(i) << "  " << (int)
                Neighb0->GetVertex(i) << endl;

      N_ = Neighb0->GetN_Vertices();
      for (i=0;i<N_;i++)
        cout << " test 1:" << i << Neighb1->GetVertex(i) << "  " << (int)
                Neighb1->GetVertex(i) << endl;
      */

      ErrThrow("Error in SetMapType: could not find vertex");
    }
  }
}

void TJoint::GetMapperRef(const int *&MapVerts, const int *&MapFaces) const
{
  if (MapType != -1)
    switch (Neighb0->GetType())
    {
      case Tetrahedron: TDatabase::MapperDB[MapTriReg0 + MapType]->
                          GetMapperRef(MapVerts, MapFaces);
                        break;

      case Brick:
      case Hexahedron: TDatabase::MapperDB[MapQuadReg0 + MapType]->
                         GetMapperRef(MapVerts, MapFaces);
                       break;
      default:
      break;	
    }
  else
  {
    ErrThrow("Error in GetMapperRef: wrong MapType ", MapType);
  }
}

void TJoint::GetMapperOrig(const int *&MapVerts, const int *&MapEdges) const
{
  if (MapType != -1)
    switch (Neighb0->GetType())
    {
      case Tetrahedron: TDatabase::MapperDB[MapTriReg0 + MapType]->
                          GetMapperOrig(MapVerts, MapEdges);
                        break;

      case Brick:
      case Hexahedron: TDatabase::MapperDB[MapQuadReg0 + MapType]->
                         GetMapperOrig(MapVerts, MapEdges);
                       break;
      default:
      break;
    }
  else
  {
    ErrThrow("Error in GetMapperOrig: wrong MapType ", MapType);
  }
}
int TJoint::GetNeighbourEdgeIndex(const TBaseCell* me, int LocEdge) const
{
  int N_Edges; //MaxLen;
  const int *TmpEV;
  const TVertex *Vert0, *Vert1;
  const TBaseCell* neigh;

  // Set neigh as the neighbour of me
  if(me == Neighb0) neigh=Neighb1;
  else neigh = Neighb0;

  if (me && neigh)
  {
    // Get Vertices of the edge
    me->GetShapeDesc()->GetEdgeVertex(TmpEV);
    Vert0 = me->GetVertex(TmpEV[LocEdge*2]);
    Vert1 = me->GetVertex(TmpEV[LocEdge*2+1]);

    neigh->GetShapeDesc()->GetEdgeVertex(TmpEV);
    N_Edges = neigh->GetShapeDesc()->GetN_Edges();
    // Iterate over all edges of neighbour
    for(int i=0; i<N_Edges; ++i)
    {
        // Check if the vertices of that edge are the ones we searched for
        if((neigh->GetVertex(TmpEV[2*i]) == Vert0 && neigh->GetVertex(TmpEV[2*i+1]) == Vert1) ||
           (neigh->GetVertex(TmpEV[2*i]) == Vert1 && neigh->GetVertex(TmpEV[2*i+1]) == Vert0))
           return i;
    }

    return -1;
  }
  else
  {
    return -1;
  }
}

#endif // __3D__

// Destructor
TJoint::~TJoint()
{
  if(Neighb0)
  { Neighb0 = nullptr;}
  
  if(Neighb1)
  { Neighb1 = nullptr;}
  
}



