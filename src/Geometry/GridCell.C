
#include <Database.h>
#include <GridCell.h>
#include <Vertex.h>
#include <stdlib.h>
#include <MooNMD_Io.h>

#ifdef __2D__
  #include <IsoBoundEdge.h>
  #include <IsoInterfaceJoint.h>
#else
  #include <BoundFace.h>
  #include "FEDatabase.h"
  #include <TetraAffin.h>
  #include "HexaTrilinear.h"
#endif
#include <cmath>

// Constructors
TGridCell::TGridCell(const TRefDesc *refdesc, int reflevel) : TBaseCell(refdesc)
{
  Vertices = new TVertex*[RefDesc->GetN_OrigVertices()];

  Children = nullptr;
  Parent = nullptr;

  RefLevel = reflevel;
}

// Destructor
TGridCell::~TGridCell()
{
  if (Children)
  {
    cerr << "Error: Deleting Cell, but there are children!!!" << endl;
    exit (-1);
  }

  delete [] Vertices;
}

// Methods
int TGridCell::SetVertex(int Vert_i, TVertex *Vert)
{
  Vertices[Vert_i] = Vert;
  return 0;
}

TVertex *TGridCell::GetVertex(int Vert_i)
{
  return Vertices[Vert_i];
}
const TVertex *TGridCell::GetVertex(int Vert_i) const
{
   return Vertices[Vert_i];
}

void TGridCell::ComputeCentroid(double& x, double& y) const
{
    // curent centroid and area
    x = 0;
    y = 0;

    auto n_vertices = this->GetN_Vertices();

    // compute the centroid and the area by dividing the area into trapezoids,
    // compute the centroids of the trapezoids, and then weight them
    // appropriately. This is based on
    // https://en.wikipedia.org/wiki/Centroid#Of_a_polygon
    // and
    // https://en.wikipedia.org/wiki/Shoelace_formula
    for (int pt_i = 0; pt_i < n_vertices; ++pt_i)
    {
      double x_i;
      double y_i;
      double z_i; // dummy value to compile in 3D
      this->GetVertex(pt_i % n_vertices)->GetCoords(x_i, y_i, z_i);
      double x_i_next;
      double y_i_next;
      this->GetVertex((pt_i+1) % n_vertices)->GetCoords(x_i_next, y_i_next, z_i);

      double volume_trapezoid = (x_i * y_i_next - x_i_next * y_i);
      x += (x_i + x_i_next) * volume_trapezoid;
      y += (y_i + y_i_next) * volume_trapezoid;
    }
    // weight appropriately and store into cell_centroids
    double area = this->GetMeasure();
    x /= 6*area;
    y /= 6*area;
}

parmoon::Point TGridCell::getCentroid() const
{
#ifdef __3D__
  int nVertices = this->GetN_Vertices();
  parmoon::Point point(0.0, 0.0, 0.0);
  
  for(int i = 0; i < nVertices; i++)
  {
    const TVertex* pVertex = this->GetVertex(i);
    point.x += pVertex->GetX();
    point.y += pVertex->GetY();
    point.z += pVertex->GetZ();
  }
  
  point.x /= nVertices;
  point.y /= nVertices;
  point.z /= nVertices;
#else
  parmoon::Point point(0.0, 0.0);
  int nVertices = this->GetN_Vertices();
  double doubledArea = 0.0;
  
  for(int i = 0; i < nVertices; i++)
  {
    double x_i;
    double y_i;
    this->GetVertex(i)->GetCoords(x_i, y_i);
    
    double x_i_next;
    double y_i_next;
    int iNext = (i+1) < nVertices ? i+1 : 0;
    this->GetVertex(iNext)->GetCoords(x_i_next, y_i_next);
    
    double commonBracket = (x_i * y_i_next - x_i_next * y_i);
    point.x += (x_i + x_i_next) * commonBracket;
    point.y += (y_i + y_i_next) * commonBracket;
    doubledArea += commonBracket;
  }
  
  point.x /= 3.0*doubledArea;
  point.y /= 3.0*doubledArea;
#endif  
  
  return point;
}

void TGridCell::ComputeReflectedPoint(int edge_j, double x_orig, double y_orig,
        double& x_reflected, double& y_reflected) const
{
  if (TDatabase::ParamDB->USE_ISOPARAMETRIC == 1)
  {
    // The mirroring relies on straight lines, so isoparametric transformations
    // are forbidden
    ErrThrow("This method is only implemented for non-isoparametric ",
        "transformation.");
  }
  // For this purpose we need the point of intersection of the two
  // straights given by
  // (x_orig,y_orig)^T + s * normal_j
  // and
  // Vj0 + t * (Vj1 - Vj0)
  // where s,t are real numbers, normal_j the outer unit normal of
  // edge_j and Vj0 and Vj1 are the Vertices related to edge_j.
  // In other words,
  // (s,t)^t = A_j^-1  b_j,
  // where
  // A_j = (normal_j, Vj0 - Vj1)
  // b_j = Vj0 - (x_orig, y_orig)^T
  //
  // Below we compute only s, since t is not needed in our case.
  //
  // With a slightly abuse of notation, we adapt the notation from
  // https://en.wikipedia.org/wiki/Invertible_matrix#Inversion_of_2_%C3%97_2_matrices
  // That is, A = ([a, b], [c, d]) =
  //            = ([ normal_j[0], Vj0[0] - Vj1[0] ],
  //               [ normal_j[1], Vj0[1] - Vj1[1] ])

  auto n_edges = this->GetN_Joints();
  double Vj0_x, Vj0_y, Vj1_x, Vj1_y, z0;
  this->GetVertex(edge_j)->GetCoords(Vj0_x, Vj0_y, z0);
  this->GetVertex((edge_j+1) % n_edges)->GetCoords(Vj1_x, Vj1_y, z0);
  auto hE = std::sqrt( (Vj1_x - Vj0_x) * (Vj1_x - Vj0_x) +
      (Vj1_y - Vj0_y) * (Vj1_y - Vj0_y) );
  double normal [2] = {(Vj1_y-Vj0_y) / hE, (Vj0_x-Vj1_x) / hE};
  double a = normal[0];
  double b = Vj0_x - Vj1_x;
  double c = normal[1];
  double d = Vj0_y - Vj1_y;

  double detA = a * d - b * c;
  if (detA == 0)
  {
    ErrThrow("The system for finding the parameters to reflect the ",
        "centroid is ill posed!");
  }
  double s;
  s = ( d * (Vj0_x - x_orig) - b * (Vj0_y - y_orig) ) / detA;

  // The reflected point we looked for is given by
  // (x_rig, y_orig)^t + 2 * s * normal_j, since 1 * s projects only to edge_j.
  x_reflected = x_orig + 2 * s * a;
  y_reflected = y_orig + 2 * s * c;
}

int TGridCell::GetN_Children() const
{
  return Children ? RefDesc->GetN_Children() : 0;
}

int TGridCell::GetN_Parents() const
{
  if ( Parent == nullptr)
    return 0;
  else
    return 1;
}

TBaseCell *TGridCell::GetChild(int C_i)
{
  return Children[C_i];
}
const TBaseCell *TGridCell::GetChild(int C_i) const
{
  return Children[C_i];
}

const TBaseCell *TGridCell::GetNeighbor(int joint_nr) const
{
  return this->GetJoint(joint_nr)->GetNeighbour(this);
}

TBaseCell *TGridCell::GetParent()
{
  return Parent;
}
const TBaseCell *TGridCell::GetParent() const
{
  return Parent;
}

int TGridCell::SetParent(TBaseCell *parent)
{
  Parent = parent;
  return 0;
}

int TGridCell::GetChildNumber(TBaseCell *Me) const
{
  int N_ = GetN_Children();
  for(int i = 0; i < N_; i++)
    if (Children[i] == Me)
      return i;
  ErrThrow("Error: could not find child");
}

void TGridCell::PS(std::ofstream& dat, double scale, double StartX,
                   double StartY, int cell_index) const
{
  bool gridcell_with_numbers = cell_index >=0;
  int i;
  double text_x1=0,text_x2=0,text_y1=0,text_y2=0;
  double x, y;
  dat<<"newpath"<<endl;
  dat << (30 + (Vertices[0]->GetX() - StartX) * scale + .5) << " " <<
         (30 + (Vertices[0]->GetY() - StartY) * scale + .5) <<
         " M" << endl;
  text_x1 = Vertices[0]->GetX();
  text_y1 = Vertices[0]->GetY();

  x = Vertices[0]->GetX();
  y = Vertices[0]->GetY();
  //if(!IsHaloCell() && IsDependentCell())

  for (i=1;i<RefDesc->GetN_OrigVertices();i++)
  {
    dat << (30 + (Vertices[i]->GetX() - StartX) * scale + .5) << " " <<
           (30 + (Vertices[i]->GetY() - StartY) * scale + .5) <<
           " L" << endl;
    if(Vertices[i]->GetX() != text_x1)
      text_x2 = Vertices[i]->GetX();
    if(Vertices[i]->GetY() != text_y1)
      text_y2 = Vertices[i]->GetY();

    x += Vertices[i]->GetX();
    y += Vertices[i]->GetY(); 
  }

  text_x1 = (text_x1 + text_x2)/2;
  text_y1 = (text_y1 + text_y2)/2;

  x /= RefDesc->GetN_OrigVertices();
  y /= RefDesc->GetN_OrigVertices();

  dat << (30 + (Vertices[0]->GetX() - StartX) * scale + .5) << " " <<
         (30 + (Vertices[0]->GetY() - StartY) * scale + .5) <<
         " L" << endl;
  dat<<"closepath"<<endl;
  if(gridcell_with_numbers)
  {
    dat << (30 + (x - StartX) * scale + .5) << " " <<
           (30 + (y - StartY) * scale + .5);
    dat << " moveto" << endl;
    dat << "(" << cell_index << ") show" << endl;
  }
  
#ifdef _MPI
if(IsHaloCell())
{
dat<<"gsave"<<endl;
//dat<< r<<" setgray"<<endl;
dat<<"0 1 0 setrgbcolor"<<endl;
dat<<"fill"<<endl;
dat<<"grestore"<<endl;
}
else if(IsDependentCell())
{
dat<<"gsave"<<endl;
dat<<"1 0 0 setrgbcolor"<<endl;
dat<<"fill"<<endl;
dat<<"grestore"<<endl;
}
else
{
dat<<"stroke"<<endl;
}
#else
 dat<<"stroke"<<endl;
#endif

#ifdef _MPI
for (i=0;i<RefDesc->GetN_OrigVertices();i++)
  {
    if(Vertices[i]->IsCrossVert())
    {
     //printf(" I'M HERE\n");
    dat << "newpath" <<endl;
    dat << (30 + (Vertices[i]->GetX() - StartX) * scale + .5) << " " <<
           (30 + (Vertices[i]->GetY() - StartY) * scale + .5) <<" "<<"M"<< endl;
    dat << (30 + (Vertices[i]->GetX() - StartX) * scale + .5) << " " <<
           (30 + (Vertices[i]->GetY() - StartY) * scale + .5) <<" "<<"5 0 360 arc fill"<< endl;
    dat<<"closepath"<<endl;
    }
  }
#endif
dat << (30 + (text_x1 - StartX) * scale + .5) << " " <<
         (30 + (text_y1 - StartY) * scale + .5) <<
         " M" << endl;
#ifdef _MPI
//dat <<"("<<GlobalCellNo<<") show"<<endl;
//dat <<"("<<CellIndex<<") show"<<endl;
dat <<"("<<GetSubDomainNo()<<") show"<<endl;
#endif
}


#ifdef __2D__
int TGridCell::Gen1RegGrid()
{
  int i, N_, ChildNumber, LocJointNum;
  const int *TmpCE, *TmpnEoE;
  int MaxLen;
  const TRefDesc *ParentRefDesc = Parent->GetRefDesc();
  TBaseCell *RefCell;

  ParentRefDesc->GetChildEdge(TmpCE, MaxLen);
  ParentRefDesc->GetNewEdgeOldEdge(TmpnEoE);

  N_ = GetN_Edges();

  for (i=0;i<N_;i++)
  {
    if (Joints[i]->GetType() != BoundaryEdge &&
        Joints[i]->GetType() != IsoBoundEdge)

      if (!Joints[i]->GetNeighbour(this))
      {
        ChildNumber = Parent->GetChildNumber(this);

        LocJointNum = TmpnEoE[TmpCE[ChildNumber * MaxLen + i]];

        if (!Parent->GetJoint(LocJointNum)->GetNeighbour(Parent))
          Parent->Ref1Reg(LocJointNum, RefCell);
      }

    if (Joints[i]->GetType() == InterfaceJoint ||
        Joints[i]->GetType() == IsoInterfaceJoint)
      if (!Joints[i]->GetNeighbour(this))
        Ref1Reg(i, RefCell);
  }

  return 0;
}

#else

int TGridCell::Gen1RegGrid()
{
  if(!Parent)
    return 0;

  int CurrLocEdge;
//  const int *TmpCE, *TmpnEoE;
//  int MaxLen;
//  const TRefDesc *ParentRefDesc = Parent->GetRefDesc(), *GrandParentRefDesc;
//  TBaseCell *RefCell;
  const TBaseCell* Grandfather;
  TBaseCell *CurrCell;
  const TJoint* LastJoint, *Joint;
  const int *TmpEF, *TmpEF2;
  int TmpEFMaxLen, TmpEF2MaxLen;

  if(RefLevel > 1 && (Grandfather = ((const TBaseCell*)Parent)->GetParent()))
  {
//    GrandParentRefDesc = Grandfather->GetRefDesc();
    Grandfather->GetShapeDesc()->GetEdgeFace(TmpEF, TmpEFMaxLen);

    for(int i=0; i<Grandfather->GetN_Edges(); ++i)
    {
      for(int j=0; j<TmpEFMaxLen; ++j)
      {
        LastJoint = Grandfather->GetJoint(TmpEF[2*i+j]);

        if(!(CurrCell = LastJoint->GetNeighbour(Grandfather)))
          continue;

        CurrLocEdge = LastJoint->GetNeighbourEdgeIndex(Grandfather, i);

        while(CurrCell != Grandfather && CurrCell)
        {
          if(CurrLocEdge == -1)
          {
            std::cout << "Error!\n"; exit(-1);
          }
          
          // If neighbour is not refined, refine it regular
          //if(CurrCell->GetRefDesc()->GetType() == NoRef)
          if(CurrCell->GetN_Children() == 0)
          {
            CurrCell->SetRegRefine();
            CurrCell->Refine(RefLevel-1);
            
            if(RefLevel>2)
            {
              for(int k=0; k<CurrCell->GetN_Children(); ++k)
                static_cast<TGridCell*>(CurrCell)->GetChild(k)->Gen1RegGrid();
            }
          }
          
          // Get next joint which contains this edge
          CurrCell->GetShapeDesc()->GetEdgeFace(TmpEF2, TmpEF2MaxLen);
          if(CurrCell->GetJoint(TmpEF2[2*CurrLocEdge]) == LastJoint)
            Joint = CurrCell->GetJoint(TmpEF2[2*CurrLocEdge+1]);
          else
            Joint = CurrCell->GetJoint(TmpEF2[2*CurrLocEdge]);

          // Get new element and the index of our edge in this element
          CurrLocEdge = Joint->GetNeighbourEdgeIndex(CurrCell, CurrLocEdge);
          CurrCell = Joint->GetNeighbour(CurrCell);
          
          LastJoint = Joint;
        }
      } // j=0..N_JointsPerEdge
    }
  }

  return 0;
}
#endif


#ifdef __2D__
int TGridCell::Ref1Reg(int LocJointNum, TBaseCell *&RefCell)
{
  int ChildNumber;
  TBaseCell *LocCell;
  const int *TmpCE, *TmpnEoE;
  int MaxLen;
  const TRefDesc *ParentRefDesc = Parent->GetRefDesc();

  ParentRefDesc->GetChildEdge(TmpCE, MaxLen);
  ParentRefDesc->GetNewEdgeOldEdge(TmpnEoE);

  ChildNumber = Parent->GetChildNumber(this);

  LocCell = Parent->GetJoint(TmpnEoE[TmpCE[ChildNumber * MaxLen +
              LocJointNum]])->GetNeighbour(Parent);

  if (LocCell)
  {
    LocCell->SetRegRefine();
    LocCell->Refine(RefLevel);
  }
  else
  {
    cerr << "Error: generation of 1-regular grid impossible" << endl;
    return -1;
  }

  RefCell = LocCell;
  return 0;
}

#else

int TGridCell::Ref1Reg(int LocJointNum, TBaseCell *&RefCell)
{
  int ChildNumber;
  TBaseCell *LocCell;
  const int *TmpCF, *TmpnFoF;
  int MaxLen;
  const TRefDesc *ParentRefDesc = Parent->GetRefDesc();

  ParentRefDesc->GetChildFace(TmpCF, MaxLen);
  ParentRefDesc->GetNewFaceOldFace(TmpnFoF);

  ChildNumber = Parent->GetChildNumber(this);

  LocCell = Parent->GetJoint(TmpnFoF[TmpCF[ChildNumber * MaxLen +
              LocJointNum]])->GetNeighbour(Parent);

  if (LocCell)
  {
    LocCell->SetRegRefine();
    LocCell->Refine(RefLevel);
  }
  else
  {
    cerr << "Error: generation of 1-regular grid impossible" << endl;
    return -1;
  }

  RefCell = LocCell;
  return 0;
}
#endif


#ifdef __2D__
int TGridCell::MakeConfClosure()
{
  int i;

  int N_ = GetN_Edges();
  // set triangle/quadrangle bit(7/8) (clip<256 means triangle )
  int clip = 1 << (N_ + 4); // mac64 warning bracket added
  // set edge bits(0-3) and number of edges (bits 4-6)
  for (i=0;i<N_;i++)
  {
    // there is a neighbour cell
    if(TBaseCell *Neighb = Joints[i]->GetNeighbour(this))
    {
      // this neighbour cell is refined regularly or it has children 
      if (Neighb->ExistChildren() || Neighb->GetClipBoard() > 512)
      {    
        clip += (1 << i) + 16;
      }
    }
  }
  // set regular refinement bit(9)
  if ((clip > 160 && clip < 256) || clip > 304)
    clip |= 512;

  ClipBoard = clip;

  N_ = clip & 496;
  if (N_ == 160 || N_ == 304)
  {
    if (clip & 128)
      i = (clip & 7) ^ 7;
    else
      i = (clip & 15) ^ 15;

    i--;
    if (i == 3) i--;
    if (i == 7) i = 3;

    if (Joints[i]->InnerJoint())
    {
      TBaseCell *Neighb = Joints[i]->GetNeighbour(this);
      if(Neighb)
        Neighb->MakeConfClosure();
      else
      {
        Ref1Reg(i, Neighb);
        Neighb->Check1Reg();
        Joints[i]->GetNeighbour(this)->MakeConfClosure();
      }
    }
  }

  return 0;
}

#else // 3D
int TGridCell::MakeConfClosure()
{
  int i, j, clip, N_Edges, LocEdge, CurrLocEdge, N_ToRefine=0, CurrClip;
  TBaseCell *CurrCell;
  const int *TmpEF, *TmpEF2;
  int TmpEFMaxLen, TmpEF2MaxLen;
  TJoint *Joint, *LastJoint;
  bool EdgeRefined;

  if(Children)
    {
      Derefine();
      std::cerr << "Element has already been refined\n";
      //exit(-1);
    }

  clip = 0; //GetClipBoard(); 

  GetShapeDesc()->GetEdgeFace(TmpEF, TmpEFMaxLen);
  N_Edges = GetN_Edges();

  for (i=0;i<N_Edges;i++)
    {
      EdgeRefined = false;
      LocEdge = i;
      for(j=0; j<TmpEFMaxLen; ++j)
        {
    LastJoint = GetJoint(TmpEF[2*LocEdge+j]);

    if(!(CurrCell = LastJoint->GetNeighbour(this)))
      continue;

    CurrLocEdge = LastJoint->GetNeighbourEdgeIndex(this, LocEdge);

    while(CurrCell != this && CurrCell)
            {
        if(CurrLocEdge == -1)
                {
      std::cout << "Error!\n"; exit(-1);
                }

        // Check if Edge is refined in CurrCell
        CurrClip = CurrCell->GetClipBoard();
                
        if(CurrClip != -1 && ((CurrClip >> CurrLocEdge) & 1))
                {
      // Check if Cell is already marked for refinement
      if((clip >> LocEdge) & 1)
        std::cout << "Should not happen\n";

      clip |= 1 << LocEdge;
      EdgeRefined = true;
      break;
                }

        // Get next joint which contains this edge
        CurrCell->GetShapeDesc()->GetEdgeFace(TmpEF2, TmpEF2MaxLen);
        if(CurrCell->GetJoint(TmpEF2[2*CurrLocEdge]) == LastJoint)
    Joint = CurrCell->GetJoint(TmpEF2[2*CurrLocEdge+1]);
        else
    Joint = CurrCell->GetJoint(TmpEF2[2*CurrLocEdge]);

        // Get new element and the index of our edge in this element
        CurrLocEdge = Joint->GetNeighbourEdgeIndex(CurrCell, CurrLocEdge);
        CurrCell = Joint->GetNeighbour(CurrCell);

        LastJoint = Joint;
            }

    if(EdgeRefined)
      break;
        } // j=0..N_JointsPerEdge

      //if(EdgeRefined)
      //    std::cout << "Edge " << i << " is refined by an other element\n";
    } // i=0--N_edges

  // Check if a RefDesc exists for this setting
  switch(clip)
    {
      // NoRef
    case 0:
      break;
      std::cout << "Unexpected behaviour!\n";
      exit(-1);
      // Bisections
    case 1:
    case 2:
    case 4:
    case 8:
    case 16:
    case 32:
      // Bis0X
    case 3:
    case 5:
    case 9:
    case 17:
    case 33:
      // Bis 1X
    case 6:
    case 10:
    case 18:
    case 34:
      // Bis2X
    case 12:
    case 20:
    case 36:
      // Bis 3X
    case 24:
    case 40:
      // Bis 4X
    case 48:
      // QuadX
    case 7:
    case 25:
    case 50:
    case 44:
      break;
      // Refinement Rule was found
      // Reg
    case 63:
    default:
      // No Suitable Refinement Rule was found

      N_Edges = GetN_Edges();
      for(i=0; i<N_Edges; ++i)
        {
    // Check if edge is marked for refinement
    if( clip >> i) continue;
//     if( (clip >> i) & 1 == 1) continue;
    
    LocEdge = i;
    for(j=0; j<TmpEFMaxLen; ++j)
            {
        LastJoint = GetJoint(TmpEF[2*LocEdge+j]);

        if(!(CurrCell = LastJoint->GetNeighbour(this)))
    continue;

        CurrLocEdge = LastJoint->GetNeighbourEdgeIndex(this, LocEdge);

        while(CurrCell != this && CurrCell)
                {
      // Add this cell to Gk
      if(CurrCell->GetClipBoard() > 0 && ((CurrCell->GetClipBoard() >> CurrLocEdge) & 1) == 1)
        std::cout << "Edge has been marked for refinement already before\n";

      if(CurrCell->GetClipBoard() != 0)
        {
          N_ToRefine++;         
          CurrCell->SetClipBoard(0);
        }
      
      // Get next joint which contains this edge
      CurrCell->GetShapeDesc()->GetEdgeFace(TmpEF2, TmpEF2MaxLen);
      if(CurrCell->GetJoint(TmpEF2[2*CurrLocEdge]) == LastJoint)
        Joint = CurrCell->GetJoint(TmpEF2[2*CurrLocEdge+1]);
      else
        Joint = CurrCell->GetJoint(TmpEF2[2*CurrLocEdge]);

      // Get new element and the index of our edge in this element
      CurrLocEdge = Joint->GetNeighbourEdgeIndex(CurrCell, CurrLocEdge);
      CurrCell = Joint->GetNeighbour(CurrCell);
      LastJoint = Joint;
                }
            } // j=0..N_JointsPerEdge
        }
      clip = 63;
    }

  if(clip == 0)
    std::cerr << "Fehler!\n";

  SetClipBoard(clip);

  return N_ToRefine;
}
#endif


#ifdef __2D__
int TGridCell::Check1Reg()
{
  int i, N_;
  TBaseCell *Neighb;

  N_ = GetN_Edges();
  for (i=0;i<N_;i++)
    if (Joints[i]->InnerJoint())
      if (!(Neighb = Joints[i]->GetNeighbour(this)))
      {
        Ref1Reg(i, Neighb);
        Neighb->Check1Reg();
      }

  return 0;
}
#else
int TGridCell::Check1Reg()
{
  int i, N_;
  TBaseCell *Neighb;

  N_ = GetN_Faces();
  for (i=0;i<N_;i++)
    if (Joints[i]->InnerJoint())
      if (!(Neighb = Joints[i]->GetNeighbour(this)))
      {
        Ref1Reg(i, Neighb);
        Neighb->Check1Reg();
      }

  return 0;
}
#endif


#ifdef __2D__
int TGridCell::SetRegRefine()
{
  if (ExistChildren()) return -1;

  switch (RefDesc->GetShapeDesc()->GetType())
  {
    case Triangle:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + TriReg];
           break;
    case Parallelogram:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + ParallReg];
           break;
    case Quadrangle:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + QuadReg];
           break;
    case Rectangle:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + RectReg];
           break;
     default:
          ErrThrow("Only 2D cells allowed");
          break;
  }

  return 0;
}

int TGridCell::Set1Refine(int i)
{
  if (ExistChildren()) return -1;

  switch (RefDesc->GetShapeDesc()->GetType())
  {
    case Triangle:
      
        switch (i)
          {
           case 0:  
           RefDesc = TDatabase::RefDescDB[N_SHAPES + TriBis0];
           break;
           
           case 1:  
           RefDesc = TDatabase::RefDescDB[N_SHAPES + TriBis1];
           break;
           
           case 2:  
           RefDesc = TDatabase::RefDescDB[N_SHAPES + TriBis2];
           break;                     

          }   
           
           
           break;
    case Parallelogram:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + ParallReg];
           break;
    case Quadrangle:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + QuadReg];
           break;
    case Rectangle:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + RectReg];
           break;
     default:
          ErrThrow("Only 2D cells allowed");
          break;
  }

  return 0;
}



#endif

#ifdef __3D__
int TGridCell::Set1Refine(int)
{
 
    cerr << "Error: Not yet implemented !!!" << endl;
    exit (-1);
 
}

int TGridCell::SetRegRefine()
{
/*
  double x0, y0, z0;
  double x1, y1, z1;
  double x2, y2, z2;
  double x3, y3, z3;
  double x4, y4, z4;
  double x5, y5, z5;
  double x6, y6, z6;
  double x7, y7, z7;
  double x8, y8, z8;
  double x9, y9, z9;
  double l0, l1, l2;
  double min;
  int index;
*/

  if (ExistChildren()) return -1;

  switch (RefDesc->GetShapeDesc()->GetType())
  {
    case Triangle:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + TriReg];
           break;
    case Parallelogram:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + ParallReg];
           break;
    case Quadrangle:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + QuadReg];
           break;
    case Rectangle:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + RectReg];
           break;
    case Tetrahedron:
         /*
           Vertices[0]->GetCoords(x0, y0, z0);
           Vertices[1]->GetCoords(x1, y1, z1);
           Vertices[2]->GetCoords(x2, y2, z2);
           Vertices[3]->GetCoords(x3, y3, z3);
           x4 = 0.5*(x0+x1);
           y4 = 0.5*(y0+y1);
           z4 = 0.5*(z0+z1);
           x5 = 0.5*(x1+x2);
           y5 = 0.5*(y1+y2);
           z5 = 0.5*(z1+z2);
           x6 = 0.5*(x0+x2);
           y6 = 0.5*(y0+y2);
           z6 = 0.5*(z0+z2);
           x7 = 0.5*(x0+x3);
           y7 = 0.5*(y0+y3);
           z7 = 0.5*(z0+z3);
           x8 = 0.5*(x1+x3);
           y8 = 0.5*(y1+y3);
           z8 = 0.5*(z1+z3);
           x9 = 0.5*(x2+x3);
           y9 = 0.5*(y2+y3);
           z9 = 0.5*(z2+z3);

           l0 = (x6-x8)*(x6-x8) + (y6-y8)*(y6-y8) +(z6-z8)*(z6-z8);
           l1 = (x5-x7)*(x5-x7) + (y5-y7)*(y5-y7) +(z5-z7)*(z5-z7);
           l2 = (x4-x9)*(x4-x9) + (y4-y9)*(y4-y9) +(z4-z9)*(z4-z9);

           min = l0;
           index = 0;
           if(l1 < min)
           {
             min = l1;
             index = 1;
           }
           if(l2 < min)
           {
             min = l2;
             index = 2;
           }

           RefDesc = TDatabase::RefDescDB[N_SHAPES + TetraReg0 + index];
          */
           RefDesc = TDatabase::RefDescDB[N_SHAPES + TetraReg];
           break;
    case Hexahedron:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + HexaReg];
           break;
    case Brick:
           RefDesc = TDatabase::RefDescDB[N_SHAPES + BrickReg];
           break;
    default:
     break;
           
  }

  return 0;
}
#endif

int TGridCell::SetNoRefinement()
{
  switch (GetType())
  {
    case Triangle:
           RefDesc = TDatabase::RefDescDB[Triangle];
           break;
    case Quadrangle:
    case Parallelogram:
    case Rectangle:
           RefDesc = TDatabase::RefDescDB[Quadrangle];
           break;
    case Tetrahedron:
           RefDesc = TDatabase::RefDescDB[Tetrahedron];
           break;
    case Hexahedron:
    case Brick:
           RefDesc = TDatabase::RefDescDB[Hexahedron];
           break;
     default:
          ErrThrow("Only 2D cells allowed. ", GetType());
     break;    
  }

  return 0;
}

int TGridCell::IsToRefine() const
{
  return RefDesc->IsToRefine() ? Children == nullptr ? true : false : false;
}

#ifdef __2D__
// return coordinates of mid point P_j on edge J_i
int TGridCell::LineMidXY(int J_i, int P_j, double &X, double &Y)
{
  const int *TmpVerts, *TmpEV;
  const double *TmpPos;
  int MaxLen;
  double T_0, T;

  TDatabase::RefDescDB[N_SHAPES + RefDesc->GetEdgeRef(J_i)]->
                  GetInnerVerts(TmpVerts, TmpPos, MaxLen);

  RefDesc->GetShapeDesc()->GetEdgeVertex(TmpEV);

  --P_j;

  if (Joints[J_i]->GetType() == BoundaryEdge ||
      Joints[J_i]->GetType() == InterfaceJoint ||
      Joints[J_i]->GetType() == IsoBoundEdge ||
      Joints[J_i]->GetType() == IsoInterfaceJoint)
  {
    T_0 = TmpPos[P_j * MaxLen + 1];

    switch(Joints[J_i]->GetType())
    {
      case IsoBoundEdge:
      T = T_0*((TIsoBoundEdge *) Joints[J_i])->GetEndParameter()
          +(1.-T_0)*((TIsoBoundEdge *) Joints[J_i])->GetStartParameter();

      ((TIsoBoundEdge *) Joints[J_i])->GetXYofT(T, X, Y);
      break;

      case IsoInterfaceJoint:
      T = T_0*((TIsoInterfaceJoint *) Joints[J_i])->GetEndParameter()
        +(1.-T_0)*((TIsoInterfaceJoint *) Joints[J_i])->GetStartParameter();

      ((TIsoInterfaceJoint *) Joints[J_i])->GetXYofT(T, X, Y);
      break;

      case InterfaceJoint:
      T = T_0*((TInterfaceJoint *) Joints[J_i])->GetEndParameter()
          +(1.-T_0)*((TInterfaceJoint *) Joints[J_i])->GetStartParameter();

      ((TInterfaceJoint *) Joints[J_i])->GetXYofT(T, X, Y);
      break;

      case BoundaryEdge:
      T = T_0 * ((TBoundEdge *) Joints[J_i])->GetEndParameter() 
          + (1. - T_0) * ((TBoundEdge *) Joints[J_i])->GetStartParameter();

      ((TBoundEdge *) Joints[J_i])->GetXYofT(T, X, Y);
    break;
     default:
          cerr << "Only edges" << endl;
     exit(-1);
      break;  
    }
  }
  else
  {
    X = Vertices[TmpEV[J_i * 2]]->GetX();
    Y = Vertices[TmpEV[J_i * 2]]->GetY();

    X += TmpPos[P_j * MaxLen + 1] * (Vertices[TmpEV[J_i * 2 + 1]]->GetX() - X);
    Y += TmpPos[P_j * MaxLen + 1] * (Vertices[TmpEV[J_i * 2 + 1]]->GetY() - Y);
  }

  return 0;
}

// return parameters on boundary of subedge SJ_j on edge J_i
int TGridCell::LineMidT(int J_i, int SJ_j, double &T_0, double &T_1)
{
  const int *TmpVerts;
  const double *TmpPos;
  int MaxLen;

  TDatabase::RefDescDB[N_SHAPES + RefDesc->GetEdgeRef(J_i)]->
                  GetInnerVerts(TmpVerts, TmpPos, MaxLen);

  if (SJ_j)
    T_0 = TmpPos[(SJ_j - 1) * MaxLen + 1];
  else
    T_0 = 0.;

  if (SJ_j + 1 != TDatabase::RefDescDB[N_SHAPES + RefDesc->
                    GetEdgeRef(J_i)]->GetN_Edges())
    T_1 = TmpPos[SJ_j * MaxLen + 1];
  else
    T_1 = 1.;

  return 0;
}
#endif

parmoon::Point TGridCell::ComputeMidOfJoint(int joint_i) const
{
  const int *TmpFV;
  int MaxLen;
#ifdef __3D__
  const int *TmpLen;
  this->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
  int n_vert_per_joint = TmpLen[joint_i];
#else // 2D
  MaxLen = 2;
  this->GetShapeDesc()->GetEdgeVertex(TmpFV);
  int n_vert_per_joint = 2;
#endif
  double x = 0;
  double y = 0;
  double z = 0;
  for(int i = 0; i < n_vert_per_joint; ++i)
  {
    auto v = this->GetVertex(TmpFV[joint_i*MaxLen + i]);
    x += v->GetX();
    y += v->GetY();
    z += v->GetY();
  }
  x /= n_vert_per_joint;
  y /= n_vert_per_joint;
  z /= n_vert_per_joint;
#ifdef __3D__
  return parmoon::Point(x, y, z);
#else
  return parmoon::Point(x, y);
#endif
}


bool TGridCell::PointInCell(parmoon::Point p) const
{
#ifdef __2D__
  int N_ = RefDesc->GetN_OrigEdges();
  for(int i = 0; i < N_; i++)
  {
    int j = (i + 1) % N_;

    double DX = Vertices[i]->GetX();
    double DY = Vertices[i]->GetY();
    // normal pointing into cell
    double NX = DY - Vertices[j]->GetY();
    double NY = Vertices[j]->GetX() - DX;

    DX = p.x - DX;
    DY = p.y - DY;

    // check if given point is identical to vertex
    if(DX*DX + DY*DY < 1e-8)
    {
      return true;
    }
    if(NX * DX + NY * DY <= -1e-8)
    {
      return false;
    }
  }
  return true;
#else
  double xi, eta, zeta;
  bool iso = TDatabase::ParamDB->USE_ISOPARAMETRIC
             && this->has_isoparametric_joint();
  if(GetType() == Hexahedron || GetType() == Brick)
  {
    auto type = iso ? ReferenceTransformation_type::HexaIsoparametric
                    : ReferenceTransformation_type::HexaTrilinear;
    auto * rt = FEDatabase::GetRefTrans3D(type);
    rt->SetCell(this); 
    rt->GetRefFromOrig(p.x, p.y, p.z, xi, eta, zeta);
    if(-1.0001 < xi   && xi < 1.0001 &&
       -1.0001 < eta  && eta < 1.0001 &&
       -1.0001 < zeta && zeta < 1.0001)
    {
      return true;
    }
  }
  else
  {
    auto type = iso ? ReferenceTransformation_type::TetraIsoparametric
                    : ReferenceTransformation_type::TetraAffin;
    auto * rt = FEDatabase::GetRefTrans3D(type);
    rt->SetCell(this); 
    rt->GetRefFromOrig(p.x, p.y, p.z, xi, eta, zeta);
    if(-1e-4 < xi && xi < 1.0001 &&
       -1e-4 < eta && eta < 1.0001 &&
       -1e-4 < zeta && zeta < 1.0001 &&
       xi + eta + zeta < 1.0001)
    {
      return true;
    }
  }
  return false;
#endif
}


bool TGridCell::test_constraint(double  A,
                                double  B, 
                                double& range_m,
                                double& range_M) const
{
  bool   ret = false;

  if( std::abs(A) >  0. )
  {
    range_m = (A>0) ? std::max( -B/A,    range_m)
                    : std::max( (1-B)/A, range_m);
    range_M = (A>0) ? std::min( (1-B)/A, range_M)
                    : std::min( -B/A,    range_M);

    ret = (range_m <= range_M);
  }
  else
  {
    ret = ( (B>=0) && (B<=1) );
  }

  return ret;
}


bool TGridCell::IsLineCutingCell(int                  direction,
#ifdef __2D__
                                  std::array<double,2> P,
#else // __3D__
                                  std::array<double,3> P,
#endif
                                 double&              lmin,
                                 double&              lmax) const
{
  double P0_direction = 0.;
  
  double range_m = 0.;
  double range_M = 1.;

#ifdef __2D__
  double detAB;
  double detAC;
  double A[2];
  double B[2];
  double C[2];
  double v0[2];
  double v1[2];
  double v2[2];
  double coeff_a[2];
  double coeff_b[2];
  
  int idx0 = (direction + 1) % 2;
  this->GetVertex(0)->GetCoords(A[0], A[1]);
  this->GetVertex(1)->GetCoords(B[0], B[1]);
  this->GetVertex(2)->GetCoords(C[0], C[1]);

  if( (P[idx0]<std::min(std::min(A[idx0],B[idx0]),C[idx0]))
   || (P[idx0]>std::max(std::max(A[idx0],B[idx0]),C[idx0])) )
  {
    return false;
  }
  else
  {
    // assign v0 v1 v2 (avoiding division by 0) and P0_direction
    detAB = B[idx0] - A[idx0];
    detAC = C[idx0] - A[idx0];
    
    if( detAB != 0. )
    {
      for( int i=0 ; i<2 ; i++ )
      {
        v0[i] = P[i] - A[i];
        v1[i] = B[i] - A[i];
        v2[i] = C[i] - A[i];
      }
    }
    if( detAC != 0. && (std::abs(detAC)>std::abs(detAB)) )
    {
      for( int i=0 ; i<2 ; i++ )
      {
        v0[i] = P[i] - A[i];
        v1[i] = C[i] - A[i];
        v2[i] = B[i] - A[i];
      }
    }
    P0_direction = A[direction];
    
    // solve the convex constrained optimization problem:
    // min/max P[direction] = P0_direction + a1*v1[direction] + a2*v2[direction]
    // with 0<= a1 + a2 <=1,
    // 0<= a2 <=1
    // 0<= a1 <=1
    // P[direction] linear function of a2 => P[direction] extrem if a2 extrem
    // (at the boundary)
    // a1 linear function of a2 => a1 extrem if a2 extrem
    // a1 = v0[idx0]/v1[idx0] + a2 * (-v2[idx0])/v1[idx0];
    // a1 = a2 * coeff_a[1] + coeff_b[1];
    // a1 + a2 = a2 * coeff_a[0] + coeff_b[0]
    // use constraints to restrict range of a2 and solve the problem

    coeff_a[1] = -v2[idx0]/v1[idx0];
    coeff_b[1] = v0[idx0]/v1[idx0];

    coeff_a[0] = coeff_a[1] + 1;
    coeff_b[0] = coeff_b[1];

    for( int i=0 ; i<2 ; i++ )
    {
      if(! TGridCell::test_constraint(coeff_a[i], coeff_b[i], range_m, range_M))
      {
        return false;
      }
    }

    // P[direction] =  P0_direction
    //               + v1[direction] * v0[idx0]/v1[idx0]
    //               + a2 * (v2[direction] - v2[idx0]/v1[idx0])
    // P[direction] = a2 * coeff_a[0] + coeff_b[0];
    coeff_a[0] =  v1[direction]*coeff_a[1]
                + v2[direction];
    coeff_b[0] =  P0_direction
                + v1[direction]*coeff_b[1];

    lmin = std::min(range_m*coeff_a[0]+coeff_b[0],
                    range_M*coeff_a[0]+coeff_b[0]);  
    lmax = std::max(range_m*coeff_a[0]+coeff_b[0],
                    range_M*coeff_a[0]+coeff_b[0]);
    
    return true;
  }  
#else // __3D__
  double detABC;
  double detABD;
  double detBCD;
  double det01;
  double det02;
  double det12;
  double det13;
  double det23;
  double A[3];
  double B[3];
  double C[3];
  double D[3];
  double v0[3];
  double v1[3];
  double v2[3];
  double v3[3];
  double coeff_a[3];
  double coeff_b[3];

  int idx1 = (direction + 1) % 3;
  int idx2 = (direction + 2) % 3;

  this->GetVertex(0)->GetCoords(A[0], A[1], A[2]);
  this->GetVertex(1)->GetCoords(B[0], B[1], B[2]);
  this->GetVertex(2)->GetCoords(C[0], C[1], C[2]);
  this->GetVertex(3)->GetCoords(D[0], D[1], D[2]);

  // choose vectors v1 v2 independent in projection plan O idx1 idx2
  // (to ensure det12 != 0)
  detABC =  (B[idx1]-A[idx1]) * (C[idx2]-B[idx2])
          - (B[idx2]-A[idx2]) * (C[idx1]-B[idx1]);
  detABD =  (B[idx1]-A[idx1]) * (D[idx2]-B[idx2])
          - (B[idx2]-A[idx2]) * (D[idx1]-B[idx1]);
  detBCD =  (C[idx1]-B[idx1]) * (D[idx2]-C[idx2])
          - (C[idx2]-B[idx2]) * (D[idx1]-C[idx1]);

  // assign v0 v1 v2 v3 and P0_direction
  if( detABC != 0.)
  {
    for( int i=0 ; i<3 ; i++ )
    {
      v0[i] = P[i] - A[i];
      v1[i] = B[i] - A[i];
      v2[i] = C[i] - A[i];
      v3[i] = D[i] - A[i];
    }

    P0_direction = A[direction];
  }

  if( (detABD != 0.) && (std::abs(detABD)>std::abs(detABC)) )
  {
    for( int i=0; i<3; i++ )
    {
      v0[i] = P[i] - A[i];
      v1[i] = B[i] - A[i];
      v2[i] = D[i] - A[i];
      v3[i] = C[i] - A[i];
    }

    P0_direction = A[direction];
  }

  if( (detBCD != 0.) && (std::abs(detBCD)>std::max(std::abs(detABC),std::abs(detABD))) )
  {
    for( int i=0 ; i<3 ; i++ )
    {
      v0[i] = P[i] - B[i];
      v1[i] = C[i] - B[i];
      v2[i] = D[i] - B[i];
      v3[i] = A[i] - B[i];
    }

    P0_direction = B[direction];
  }

  // solve the convex constrained optimization problem:
  // min/max P[direction] =  P0_direction + a1*v1[direction] + a2*v2[direction]
  //                       + a3*v3[direction]
  // with 0<= a1 + a2 + a3 <=1,
  // 0<= a3 <=1
  // 0<= a2 <=1
  // 0<= a1 <=1
  // P[direction] linear function of a3 => P[direction] extrem if a3 extrem
  // (at the boundary)
  // a1 and a2 linear functions of a3 => a1 and a2 extrem if a3 extrem
  // a1 = (v0[idx1] + v2[idx1]*det01/det12)/v1[idx1]
  //    + a3 * (det13/det12*v2[idx1] + v3[idx1])/v1[idx1];
  // a1 = a3 * coeff_a[1] + coeff_b[1];
  // a2 = -det01/det12 - a3 * det13/det12;
  // a2 = a3 * coeff_a[2] + coeff_b[2];
  // a1 + a2 + a3 = a3 * coeff_a[0] + coeff_b[0]
  // use constraints to restrict range of a3 and solve the problem

  det01 = (v0[idx1]*v1[idx2] - v0[idx2]*v1[idx1]);
  det02 = (v0[idx1]*v2[idx2] - v0[idx2]*v2[idx1]);
  det12 = (v1[idx1]*v2[idx2] - v1[idx2]*v2[idx1]);
  det13 = (v1[idx1]*v3[idx2] - v1[idx2]*v3[idx1]);  
  det23 = (v2[idx1]*v3[idx2] - v2[idx2]*v3[idx1]);

  coeff_a[2] = -det13/det12;
  coeff_b[2] = -det01/det12;

  coeff_a[1] = det23/det12;
  coeff_b[1] = det02/det12;

  coeff_a[0] = coeff_a[1] + coeff_a[2] + 1;
  coeff_b[0] = coeff_b[1] + coeff_b[2];

  for( int i=0 ; i<3 ; i++ )
  {
    if( ! TGridCell::test_constraint(coeff_a[i], coeff_b[i], range_m, range_M) )
    {
      return false;
    }
  }

  // P[direction] =  P0_direction
  //               + v1[direction]/v1[idx1] * (v0[idx1]+det01/det12*v2[idx1])
  //               - v2[direction] * det01/det12
  //               + a3 * v3[direction]
  //                    * ( 1 - v2[direction]*det13/det12
  //                      + v1[direction]/v1[idx1] * ( det13/det12*v2[idx1]
  //                                                 + v3[idx1]));
  // P[direction] = a3 * coeff_a[0] + coeff_b[0];

  coeff_a[0] =  v1[direction]*coeff_a[1]
              + v2[direction]*coeff_a[2]
              + v3[direction];
  coeff_b[0] =  P0_direction
              + v1[direction]*coeff_b[1]
              + v2[direction]*coeff_b[2];

  lmin = std::min(range_m*coeff_a[0]+coeff_b[0], range_M*coeff_a[0]+coeff_b[0]);  
  lmax = std::max(range_m*coeff_a[0]+coeff_b[0], range_M*coeff_a[0]+coeff_b[0]);

  return true;
#endif
//   double P0_direction;
//   double detABC;
//   double detABD;
//   double detBCD;
//   double det01;
//   double det02;
//   double det12;
//   double det13;
//   double det23;
//   double A[3];
//   double B[3];
//   double C[3];
//   double D[3];
//   double v0[3];
//   double v1[3];
//   double v2[3];
//   double v3[3];
//   double coeff_a[3];
//   double coeff_b[3];
// 
//   int idx1 = (direction + 1) % 3;
//   int idx2 = (direction + 2) % 3;
// 
//   double range_m = 0.;
//   double range_M = 1.;
// 
//   this->GetVertex(0)->GetCoords(A[0], A[1], A[2]);
//   this->GetVertex(1)->GetCoords(B[0], B[1], B[2]);
//   this->GetVertex(2)->GetCoords(C[0], C[1], C[2]);
//   this->GetVertex(3)->GetCoords(D[0], D[1], D[2]);
// 
//   // choose vectors v1 v2 independent in projection plan O idx1 idx2
//   // (to ensure det12 != 0)
//   detABC =  (B[idx1]-A[idx1]) * (C[idx2]-B[idx2])
//           - (B[idx2]-A[idx2]) * (C[idx1]-B[idx1]);
//   detABD =  (B[idx1]-A[idx1]) * (D[idx2]-B[idx2])
//           - (B[idx2]-A[idx2]) * (D[idx1]-B[idx1]);
//   detBCD =  (C[idx1]-B[idx1]) * (D[idx2]-C[idx2])
//           - (C[idx2]-B[idx2]) * (D[idx1]-C[idx1]);
// 
//   // assign v0 v1 v2 v3 and P0_direction
//   if( detABC != 0.)
//   {
//     for( int i=0 ; i<3 ; i++ )
//     {
//       v0[i] = P[i] - A[i];
//       v1[i] = B[i] - A[i];
//       v2[i] = C[i] - A[i];
//       v3[i] = D[i] - A[i];
//     }
// 
//     P0_direction = A[direction];
//   }
// 
//   if( (detABD != 0.) && (std::abs(detABD)>std::abs(detABC)) )
//   {
//     for( int i=0; i<3; i++ )
//     {
//       v0[i] = P[i] - A[i];
//       v1[i] = B[i] - A[i];
//       v2[i] = D[i] - A[i];
//       v3[i] = C[i] - A[i];
//     }
// 
//     P0_direction = A[direction];
//   }
// 
//   if( (detBCD != 0.) && (std::abs(detBCD)>std::max(std::abs(detABC),std::abs(detABD))) )
//   {
//     for( int i=0 ; i<3 ; i++ )
//     {
//       v0[i] = P[i] - B[i];
//       v1[i] = C[i] - B[i];
//       v2[i] = D[i] - B[i];
//       v3[i] = A[i] - B[i];
//     }
// 
//     P0_direction = B[direction];
//   }
// 
//   // solve the convex constrained optimization problem:
//   // min/max P[direction] =  P0_direction + a1*v1[direction] + a2*v2[direction]
//   //                       + a3*v3[direction]
//   // with 0<= a1 + a2 + a3 <=1,
//   // 0<= a3 <=1
//   // 0<= a2 <=1
//   // 0<= a1 <=1
//   // P[direction] linear function of a3 => P[direction] extrem if a3 extrem
//   // (at the boundary)
//   // a1 and a2 linear functions of a3 => a1 and a2 extrem if a3 extrem
//   // a1 = (v0[idx1] + v2[idx1]*det01/det12)/v1[idx1]
//   //    + a3 * (det13/det12*v2[idx1] + v3[idx1])/v1[idx1];
//   // a1 = a3 * coeff_a[1] + coeff_b[1];
//   // a2 = -det01/det12 - a3 * det13/det12;
//   // a2 = a3 * coeff_a[2] + coeff_b[2];
//   // a1 + a2 + a3 = a3 * coeff_a[0] + coeff_b[0]
//   // use constraints to restrict range of a3 and solve the problem
// 
//   det01 = (v0[idx1]*v1[idx2] - v0[idx2]*v1[idx1]);
//   det02 = (v0[idx1]*v2[idx2] - v0[idx2]*v2[idx1]);
//   det12 = (v1[idx1]*v2[idx2] - v1[idx2]*v2[idx1]);
//   det13 = (v1[idx1]*v3[idx2] - v1[idx2]*v3[idx1]);  
//   det23 = (v2[idx1]*v3[idx2] - v2[idx2]*v3[idx1]);
// 
//   coeff_a[2] = -det13/det12;
//   coeff_b[2] = -det01/det12;
// 
//   coeff_a[1] = det23/det12;
//   coeff_b[1] = det02/det12;
// 
//   coeff_a[0] = coeff_a[1] + coeff_a[2] + 1;
//   coeff_b[0] = coeff_b[1] + coeff_b[2];
// 
//   for( int i=0 ; i<3 ; i++ )
//   {
//     if( ! TGridCell::test_constraint(coeff_a[i], coeff_b[i], range_m, range_M) )
//     {
//       return false;
//     }
//   }
// 
//   // P[direction] =  P0_direction
//   //               + v1[direction]/v1[idx1] * (v0[idx1]+det01/det12*v2[idx1])
//   //               - v2[direction] * det01/det12
//   //               + a3 * v3[direction]
//   //                    * ( 1 - v2[direction]*det13/det12
//   //                      + v1[direction]/v1[idx1] * ( det13/det12*v2[idx1]
//   //                                                 + v3[idx1]));
//   // P[direction] = a3 * coeff_a[0] + coeff_b[0];
// 
//   coeff_a[0] =  v1[direction]*coeff_a[1]
//               + v2[direction]*coeff_a[2]
//               + v3[direction];
//   coeff_b[0] =  P0_direction
//               + v1[direction]*coeff_b[1]
//               + v2[direction]*coeff_b[2];
// 
//   lmin = std::min(range_m*coeff_a[0]+coeff_b[0], range_M*coeff_a[0]+coeff_b[0]);  
//   lmax = std::max(range_m*coeff_a[0]+coeff_b[0], range_M*coeff_a[0]+coeff_b[0]);
// 
//   return true;
}


int TGridCell::GetGeoLevel()
{
  int i = 0;
  const TBaseCell *parent = this;

  while ((parent = parent->GetParent()))
    i++;

  return i;
}

int TGridCell::GetSubGridID() const
{
  if (Parent)
    return Parent->GetSubGridID();
  else
    return -1;
}

int TGridCell::get_n_boundary_joints() const
{
  int bdry_joints = 0;
  int n_joints = GetN_Joints();
  for(int i = 0; i < n_joints; i++)
  {
    if (Joints[i]->GetType() == BoundaryEdge ||
        Joints[i]->GetType() == IsoBoundEdge ||
        Joints[i]->GetType() == BoundaryFace ||
        Joints[i]->GetType() == IsoBoundFace)
    {
      bdry_joints++;
    }
  }
  return bdry_joints;
}

void TGridCell::check() const
{
  // check if the joints know 'this' as neighbor
  auto n_joints = this->GetN_Joints();
  auto n_verts = this->GetN_Vertices();
  std::stringstream out;
  out << "cell " << this << "  number of joints " << n_joints << endl;
  for(auto i = 0; i < n_verts; ++i)
  {
    const auto * vert = this->GetVertex(i);
    Output::print(" vertex ", i, " adress ", (void*)vert, "  at ", vert);
  }
  
  for(auto i = 0; i < n_joints; ++i)
  {
    const auto * joint = this->GetJoint(i);
    out << " joint " << i << " address " << joint << " type "
        << joint->GetType() << endl;
    out << "    neighbors " << joint->GetNeighbour(0) << " and "
        << joint->GetNeighbour(1) << endl;
    if((joint->GetNeighbour(0) != this && joint->GetNeighbour(1) != this)
       && (joint->GetType() == JointEqN 
           || joint->GetType() == InnerInterfaceJoint) )
    {
      ErrThrow("One of the joints of this cell does not know the cell as a ",
               "neighbor");
    }
  }
  Output::print(out.str());
}


