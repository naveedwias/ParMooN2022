
#ifdef _MPI
#  include "mpi.h"
#endif

#include <BaseCell.h>
#include <Joint.h>
#include <stdlib.h>

#include <BoundComp3D.h>
#include <BoundEdge.h>
#include <BoundFace.h>
#include <Edge.h>
#include <Point.h>
#include "Database.h"

using namespace parmoon;

// Constructor
TBaseCell::TBaseCell(const TRefDesc *refdesc)
{
 int i, N_;

  RefDesc = refdesc;

  N_ = RefDesc->GetShapeDesc()->GetN_Joints();  
  Joints = new TJoint*[N_];

  for (i=0;i<N_;i++)
   Joints[i] = nullptr;

  ClipBoard = 0;
  Phase_ID = 0;
  Reference_ID = 0;
  CellIndex = -1;
  region = 0;
  LayerCell = 0;  
  
  normalOrientation = nullptr;

#ifdef __3D__   
  N_ = RefDesc->GetN_OrigEdges();
  Edges = new TEdge*[N_];
  for (i=0;i<N_;i++)
   Edges[i] = nullptr;
#endif

#ifdef  _MPI
  ClipBoard_Par  = -1;
  SubDomainNumber = 0;
  GlobalCellNo = -1;
  SubDomainLocalCellNo = -1;

  OwnCell=false;
  HaloCell=false;
  SubDomainInterfaceCell=false;
  DependentCell=false;
  N_NeibProcesses = 0;
  NeibProcessesIds = nullptr;
#endif
}

// Destructor
TBaseCell::~TBaseCell()
{
  int i, N_;
  TJoint *CurrJoint;

#ifdef __2D__
  N_ = RefDesc->GetN_OrigEdges();
#else
  N_ = RefDesc->GetN_OrigFaces();
#endif
  for (i=0;i<N_;i++)
  {
    CurrJoint = Joints[i];
    switch (CurrJoint->GetType())
    {
      case JointEqN:
      case InterfaceJoint:
      case IsoInterfaceJoint:
      case InnerInterfaceJoint:
#ifdef __3D__
      case InterfaceJoint3D:
      case IsoInterfaceJoint3D:
#endif
        if (CurrJoint->GetNeighbour(this))
	{ CurrJoint->Delete(this);}
        else
	{ delete CurrJoint;}
        break;

      default:
        delete CurrJoint;
        break;
    }
  }

  delete[] Joints;
}

double TBaseCell::Get_hK(int cell_measure) const
{
  switch (cell_measure)
  {
    case 0:                                     // diameter
      return this->GetDiameter();
      //case 1: // with reference map
      //Output::print("cell measure ");
      //return this->GetLengthWithReferenceMap();
    case 2:                                     // shortest edge
      return this->GetShortestEdge();
      break;
    case 1:                                     // with reference map
    case 3:                                     // measure
      return std::sqrt(this->GetMeasure());
      break;
    case 4:                                     // mesh size in convection direction, this is just a dummy
      return this->GetDiameter();
      break;
    case 5:                                     // take value from an array
      // this is in general not the diameter but a pw constant value
      // which is needed for some reasons
      return this->GetDiameter();
      break;
    default:                                    // diameter
      Output::print("CELL_MEASURE ", cell_measure, " not available!!!");
      return this->GetDiameter();
      break;
  }
}

bool TBaseCell::has_isoparametric_joint() const
{
  int n_joints = RefDesc->GetShapeDesc()->GetN_Joints();
  for(int i = 0; i < n_joints; i++)
  {
    if(Joints[i]->is_isoparametric())
      return true;
  }
  return false;
}


// added 25.04.2010 for fixing refinement problem
#ifdef __3D__
void TBaseCell::CorrectBoundaryVertices(TVertex **NewVertices, TJoint **NewJoints)
{
  int N_NewFaces = RefDesc->GetN_Faces();
  const int *FaceVertex;
  int MaxN_VpF;
  
  const TBoundComp3D *BoundComp;
  TVertex *Vertex;
  double x, y, z, t, s;

  for (int i=0;i<N_NewFaces;++i)
  {
    if ( NewJoints[i]->GetType() == BoundaryFace )
    {
      BoundComp = ((const TBoundFace*) NewJoints[i])->GetBoundComp();
      
      if ( BoundComp->GetType() == Plane ) continue;
      
      RefDesc->GetFaceVertex(FaceVertex, MaxN_VpF);
      
      for (int j=0;j<MaxN_VpF;++j)
      {
	Vertex = NewVertices[FaceVertex[MaxN_VpF*i+j]];
	
	Vertex->GetCoords(x ,y, z);
	
	BoundComp->GetTSofXYZ(x, y, z, t, s);
	BoundComp->GetXYZofTS(t, s, x, y, z);
	
	Vertex->SetCoords(x, y, z);
      }     
    }
    else continue;
  }
  
}
#endif

// on each joint, decide whether the 
// global normal is outgoing (+1) or ingoing (-1)
void TBaseCell::SetNormalOrientation()
{
  TBaseCell *neighbCell;
  int nEdges = RefDesc->GetN_OrigEdges();
  #ifdef __3D__
  nEdges = RefDesc->GetN_OrigFaces();
  #endif
  if(normalOrientation != nullptr) // nothing more to do
    return;
  normalOrientation = new int[nEdges];
  for (int i=0; i<nEdges;i++)
  {
    normalOrientation[i] = 1;
    const TJoint *joint = Joints[i];
    if(joint->InnerJoint())
    {
      neighbCell = joint->GetNeighbour(this);
      // note: if neighbCell==nullptr, then this joint has a hanging vertex and
      // the normal points from the fine to the coarse side. Also note that the
      // parent of the fine cells (i.e., the neighbor of the coarse side which
      // is not part of the collection) should have been assigned a cell index
      // -1.
      if(neighbCell && neighbCell->GetCellIndex() < GetCellIndex()&& 
         Reference_ID == neighbCell->GetReference_ID()) 
        normalOrientation[i] = -1;
    }
  }
}

//LB ====================================================
bool TBaseCell::IsBoundaryCell( int BoundComp_id ) const
{   int num_inner_joints = 0;
    for(int j = 0;  j < this->GetN_Joints(); j++)
    {   const TJoint *joint = this->GetJoint(j);
        if( joint->InnerJoint())
        {
            num_inner_joints=num_inner_joints+1;
            continue;
        }
        else if(!(joint->InnerJoint()) && BoundComp_id == -4711)
            {
                return true;  
            }
        else if(!(joint->InnerJoint()) && BoundComp_id != -4711)
            {
#ifdef __2D__
               const TBoundEdge *boundedge = (const TBoundEdge *)joint;
                const TBoundComp *BoundComp = boundedge->GetBoundComp();
#elif __3D__
               const TBoundFace *boundface = (const TBoundFace *)joint;
                const TBoundComp *BoundComp = boundface->GetBoundComp();
#endif //__3D__
                if (BoundComp->GetID() == BoundComp_id)
                {
                    return true;
                }
            }
        if(this->GetN_Joints()-1  == num_inner_joints )
            return false;
    }
return false;   
}
//LB ====================================================

#ifdef __3D__
void TBaseCell::computeNormalAndTransformationData(int m,
					std::vector<double>& normal,
					double &transformationDeterminant) const
{
    const int *faceVertexMap, *faceVertexMapLength;
    int maxNVerticesPerFace;
    // For the current cell, get information of faces and local vertices
    // faceVertexMap should be seen as an array of arrays, e.g.
    // faceVertexMap = { {a,b,c},{b,c,d},{a,c,d},{a,b,d}}
    // where faceVertexMap[i] contains the id of vertices defining face i
    // faceVertexMapLength is an array specifying the length of each list
    // note: in the case that faces of an element have differennt number of
    // vertices (e.g. a pyramid), the faceVertexMap lists have all lenght equal to
    // maxNVerticesPerFace, and these are filled with 0 for the faces with less vertices
    this->RefDesc->GetShapeDesc()->GetFaceVertex(faceVertexMap,faceVertexMapLength,maxNVerticesPerFace);
    // simplify: number of vertices on face m (m=joint_id)
    size_t nFaceVertices = faceVertexMapLength[ m ];
    std::vector< Point > faceVertices(nFaceVertices,Point((unsigned int) 3));
    for (size_t l1=0; l1<nFaceVertices; l1++) {
        double _x,_y,_z;
        this->GetVertex(faceVertexMap[ m*maxNVerticesPerFace+l1 ])->GetCoords(_x,_y,_z);
        faceVertices[l1].x = _x;
        faceVertices[l1].y = _y;
        faceVertices[l1].z = _z;
    }
    
    normal.resize(3);
    double xc1, yc1, zc1, xc2, yc2, zc2;
    switch(faceVertices.size()) {
       case 3:
            // compute the 2 vectors that span the plane containing the current face
            xc1 = faceVertices[1].x - faceVertices[0].x;
            xc2 = faceVertices[2].x - faceVertices[0].x;
            
            yc1 = faceVertices[1].y - faceVertices[0].y;
            yc2 = faceVertices[2].y - faceVertices[0].y;
            
            zc1 = faceVertices[1].z - faceVertices[0].z;
            zc2 = faceVertices[2].z - faceVertices[0].z;

            // plane spanned by vectors v1=(xc1, yc1, zc1) and v2=(xc2, yc2, zc2)
            // Area of the triangle: 0.5*||v1 x v2||
            // normed Normal vector = v1 x v2/||v1 x v2||
            // Area of reference triangle (0,0)-(0,1)-(1,0): 1/2*g*h=0.5
            // Determinant of tranform.: A(triangle)/A(ref triangle) = ||v1 x v2||
            normal[0] = yc1*zc2 - zc1*yc2;
            normal[1] = zc1*xc2 - xc1*zc2;
            normal[2] = xc1*yc2 - yc1*xc2;
            // determinant of reference trafo in order to get a normed normal vector
            transformationDeterminant =
            std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
            normal[0] /= transformationDeterminant;
            normal[1] /= transformationDeterminant;
            normal[2] /= transformationDeterminant;
            
            break;
            
        case 4:
            // We consider a quadrilateral (P0,P1,P2,P3) as composed by 2 triangles
            // T1: P0,P1,P2
            // T2: P2,P3,P0
            // and we do the same as above (twice)
            // normed normal: ( (P1-P0) x (P2-P0) ) / || (P1-P0) x (P2-P0) ||
            // area: || (P1-P0) x (P2-P0) || / 2 + || (P3-P2) x (P0-P2) || / 2
            // area reference element [-1,1]x[-1,1]: 4
            // first triangle
            xc1 = faceVertices[1].x - faceVertices[0].x;
            xc2 = faceVertices[2].x - faceVertices[0].x;
            
            yc1 = faceVertices[1].y - faceVertices[0].y;
            yc2 = faceVertices[2].y - faceVertices[0].y;
            
            zc1 = faceVertices[1].z - faceVertices[0].z;
            zc2 = faceVertices[2].z - faceVertices[0].z;
            
            // normal vector (the same (except for length) for T1 and T2)
            normal[0] = yc1*zc2 - zc1*yc2;
            normal[1] = zc1*xc2 - xc1*zc2;
            normal[2] = xc1*yc2 - yc1*xc2;
            
            // determinant of reference transformation in order to get a normed normal vector
            double areaT1 =
            std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
            areaT1 /= 2.0;
            // second triangle
            xc1 = faceVertices[3].x - faceVertices[2].x;
            xc2 = faceVertices[0].x - faceVertices[2].x;
            
            yc1 = faceVertices[3].y - faceVertices[2].y;
            yc2 = faceVertices[0].y - faceVertices[2].y;
            
            zc1 = faceVertices[3].z - faceVertices[2].z;
            zc2 = faceVertices[0].z - faceVertices[2].z;
            
            
            // normal vector (the same (except for length) for T1 and T2)
            normal[0] = yc1*zc2 - zc1*yc2;
            normal[1] = zc1*xc2 - xc1*zc2;
            normal[2] = xc1*yc2 - yc1*xc2;
            
            // determinant of reference trasformation in order to get a normed normal vector
            double areaT2 =
            std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
            normal[0] /= areaT2;
            normal[1] /= areaT2;
            normal[2] /= areaT2;
            
            areaT2 /= 2.0;
            
            // note: the reference element is [-1,1] x [-1,1]
            transformationDeterminant = (areaT1+areaT2)/4.;
            
            break;
            
    } // tria or quads
    
}
#endif

// Methods
#ifdef  _MPI
void TBaseCell::SetNeibProcessesIds(int *Neiblist)
{
 int i;

 if(N_NeibProcesses==0)
  {
   printf(" Set the N_NeibProcesses for this cell first !!!! \n");
   MPI_Finalize();
   exit(0);
  }

 if(NeibProcessesIds) delete [] NeibProcessesIds;

 NeibProcessesIds = new int[N_NeibProcesses];

 for(i=0;i<N_NeibProcesses;i++)
  NeibProcessesIds[i] = Neiblist[i];

}
#endif

double norm(double x,double y, double z)
{
	return std::sqrt(x*x+y*y+z*z);
}
#ifdef __3D__
//Returns normed triple product
double triple_product(const TVertex * vert0, const TVertex * vert1, const TVertex * vert2, const TVertex * vert3)
{
	double x0,x1,x2,x3,y0,y1,y2,y3,z0,z1,z2,z3;
	vert0->GetCoords(x0,y0,z0);
	vert1->GetCoords(x1,y1,z1);
	vert2->GetCoords(x2,y2,z2);
	vert3->GetCoords(x3,y3,z3);

	//calculate the triple product (volume of the hexahedron)
	//TODO: Normiere den Spass, damit sich ablesen laesst, ob Zelle entartet ist.
	double product =((x2-x0)*(y3-y0)-(x3-x0)*(y2-y0))*(z1-z0)
			+((x3-x0)*(y1-y0)-(x1-x0)*(y3-y0))*(z2-z0)
			+((x1-x0)*(y2-y0)-(x2-x0)*(y1-y0))*(z3-z0);
	//double norm = norm(x1-x0,y1-0,z1-z0)*norm(x2-x0,y2-y0,z2-z0)*norm(x3-x0,y3-y0,z3-z0);
	//if (norm==0)
	//	ErrThrow("Two vertices");
	return product;
}

//Check if the edge between vert0 and vert1 and the edge between vert2 and vert3 cross
//Solve LGS
/*{ (x0-x1)  (x2-x3)	{ s		(x3-x1)
 *  (y0-y1)  (y2-y3)	  t}  =	(y3-y1)
 *  (z0-z1)  (z2-z3) }	  		(z3-z1)
 * */
//This method may be numerical instable.
bool find_cross_point(const TVertex * vert0, const TVertex * vert1,
		const TVertex * vert2, const TVertex * vert3)
		//double s,double t
{
	//Input check: Are the vertices pairwise different
	const TVertex* Vertices[4]={vert0 , vert1 , vert2 , vert3};
	for (int iter1 =0; iter1 <4; iter1++)
	{
		for (int iter2=iter1+1; iter2 <4; iter2++)
		{
			//Operator < (actually <=) compares vertices coordinate-wise.
			if ((!(*Vertices[iter1] < *Vertices[iter2])) && (!(*Vertices[iter2] < *Vertices[iter1])))
				{
				Output::print(vert0);
				Output::print(vert1);
				Output::print(vert2);
				Output::print(vert3);
				ErrThrow("Two vertices are the same");
				}
		}
	}

	//Compute cross point (Solve LGS)
	double s,t;

	double x0,x1,x2,x3,y0,y1,y2,y3,z0,z1,z2,z3;
	vert0->GetCoords(x0,y0,z0);
	vert1->GetCoords(x1,y1,z1);
	vert2->GetCoords(x2,y2,z2);
	vert3->GetCoords(x3,y3,z3);

	double a = (x0-x1);
	double b = (x2-x3);
	double x = (x3-x1);
	double c = (y0-y1);
	double d = (y2-y3);
	double y = (y3-y1);
	double e = (z0-z1);
	double f = (z2-z3);
	double z = (z3-z1);

	double det = a*d-c*b;
	if (det!=0)
	{
		s=(d*x-b*y)/det;
		t=(a*y-c*x)/det;

		if (std::abs(e*s+f*t-z)>1e-10)
			return false; //Edges are parallel
	}
	else
	{
		det = c*f-e*d;
			if (det!=0)
			{
				s=(f*y-d*z)/det;
				t=(c*z-e*y)/det;

				if (std::abs(a*s+b*t-x)>1e-10)
					return false; //Edges are parallel
			}
			else
				return false; //Edges are parallel or Vertices are the same...

	}
	if (0<=s && s<=1 && 0<=t && t<=1)
		{
		Output::print("s ",s, " t ", t);
		return true;
		}

	return false;
}



bool TBaseCell::check_orientation() const
{
	#ifdef __3D__
	if (this->GetShapeDesc()->GetType()==Tetrahedron)
	{

		//calculate the triple product (volume of the hexahedron)
		double product =triple_product(
				this->GetVertex(0),this->GetVertex(1),
				this->GetVertex(2),this->GetVertex(3));
		return (product>0);
	}
	else if(this->GetShapeDesc()->GetType()==Hexahedron || this->GetShapeDesc()->GetType()==Brick)
	{
		//calculate the triple product (volume of the hexahedron)
		double product =triple_product(
				this->GetVertex(0),this->GetVertex(1),
				this->GetVertex(2),this->GetVertex(4));
		return (product>0);
	}
	else
	{
		Output::print("cell shape: ", this->GetShapeDesc()->GetType());
		ErrThrow("check_orientation: This cell shape is not implemented!");
		return false;
	}
    #endif

}

bool TBaseCell::check_shape() const
{
	if (this->GetShapeDesc()->GetType()==Tetrahedron)
	{
		//calculate the triple product (volume of the hexahedron)
		//if it is not equal 0 then the vertices are not in a plane
		double product =triple_product(
				this->GetVertex(0),this->GetVertex(1),
				this->GetVertex(2),this->GetVertex(3));
		return (std::abs(product)>1e-10);
	}
	else if(this->GetShapeDesc()->GetType()==Hexahedron|| this->GetShapeDesc()->GetType()==Brick)
	//TODO: Bricks have still more properties to be checked... use CheckHexa()
	{
		const int *vert_per_face;
		const int *n_vert_per_face;
	    int max_n;
		this->GetShapeDesc()->GetFaceVertex(vert_per_face, n_vert_per_face, max_n);

		for(int face_id=0; face_id<6; face_id++)
		{
			int vert0_id =vert_per_face[max_n*face_id +0];
			int vert1_id =vert_per_face[max_n*face_id +1];
			int vert2_id =vert_per_face[max_n*face_id +2];
			int vert3_id =vert_per_face[max_n*face_id +3];

			double product =triple_product(
				this->GetVertex(vert0_id),this->GetVertex(vert1_id),
				this->GetVertex(vert2_id),this->GetVertex(vert3_id));
			if(std::abs(product) >1e-10)
			{
				Output::print("Computing errors may cause that the cell is considered degenerated.");//TODO:
				return false;
			}

			//check if opposite edges cross
			if (find_cross_point(
					this->GetVertex(vert0_id),this->GetVertex(vert1_id),
					this->GetVertex(vert2_id),this->GetVertex(vert3_id)))
				ErrThrow("Opposite edges of a face cross.");

			if (find_cross_point(
					this->GetVertex(vert0_id),this->GetVertex(vert3_id),
					this->GetVertex(vert1_id),this->GetVertex(vert2_id)))
				ErrThrow("Opposite edges of a face cross.");

		} //end face
	}
	else
	{
		ErrThrow("check_shape: For this cell shape has not implemented a check method!");
		return false;
	}

	return true;
}

#endif //__3D__


double TBaseCell::ComputeDiameterOfJoint( const int& joint_index )
{
  auto joint = this->GetJoint(joint_index);
  double diam;
  if (joint->GetNeighb(0) == nullptr)
  {
    joint->SetNeighbour(this);
    diam = joint->GetDiameter();
    joint->Delete(this);
  }
  else
  {
    diam = joint->GetDiameter();
  }
  return diam;
}


std::array<double, 3> TBaseCell::ComputeOuterJointNormal( const int& joint_index
    )
{
  if (TDatabase::ParamDB->USE_ISOPARAMETRIC == 1)
  {
    ErrThrow("This function can only handle non-isoparametric joints.");
  }
  std::array<double, 3> outer_joint_normal;
  // First we have to determine if the cell is a 2D cell or a 3D cell. This is
  // encoded in the shape the call has.
  auto shape = this->GetType();
  // These shape types belong to 2D shapes, see also Shapes.
  if (shape >= 1 && shape <= 4)
  {
    auto n_joint = this->GetN_Joints();
    double x0, y0, x1, y1, z0, z1;
    this->GetVertex(joint_index)->GetCoords(x0, y0, z0);
    this->GetVertex((joint_index+1) % n_joint)->GetCoords(x1, y1, z1);
    if (z0 != 0 || z1 != 0)
    {
      ErrThrow("This method works only in the 2D domain lies in plane with the",
          " z coordinate.");
    }
    auto joint_diam = this->ComputeDiameterOfJoint(joint_index);
    outer_joint_normal[0] = (y1-y0) / joint_diam;
    outer_joint_normal[1] = (x0-x1) / joint_diam;
    outer_joint_normal[2] = 0;
  }
  // These shape types belong to 3D shapes, see also Shapes.
  else if (shape >= 5 && shape <= 7)
  {
    // We compute the unit normal using the cross product rule, i.e. for any two
    // linearly independent vectors in the face the cross product of these two
    // vectors is a (non-unit) normal of the face. Let V_0, V_1, V_2 be
    // three vertices of the facet. Then the vectors for which the cross product
    // is computed are given by vec_0:=(V_1 - V_0) and vec_1:=(V_2 - V_0).

    // Get information about the vertices that build the face with number
    // joint_indext. This information is precisely given by FaceVertex, see also
    // TShapeDesc.
    auto shape_desc = this->GetShapeDesc();
    const int* face_vertex;
    const int* face_vertex_len;
    int max_face_vertex_len;
    shape_desc->GetFaceVertex(face_vertex, face_vertex_len,
        max_face_vertex_len);

    // Find the position in FaceVertex where the vertices of the joint_index-th face
    // start. The ordering in FaceVertex is such that first the vertices of the
    // 0th face are specified, then of the 1st face and so on. Therefore, the
    // starting point for our current face is the sum of all vertices that
    // belong to faces with a smaller local face number.
    int starting_index = 0;
    for (int joint_nr = 0; joint_nr < joint_index; ++joint_nr)
    {
      starting_index += face_vertex_len[joint_nr];
    }

    // Get the coordinates of three vertices of the face. I try to use the same
    // notation as in the beginning of this if block. To remind you, V_i is not
    // the ith vertex in the cell but the i-th vertex in the face.
    std::array<double, 3> V_0;
    std::array<double, 3> V_1;
    std::array<double, 3> V_2;
    this->GetVertex(face_vertex[starting_index + 0])
      ->GetCoords(V_0[0], V_0[1], V_0[2]);
    this->GetVertex(face_vertex[starting_index + 1])
      ->GetCoords(V_1[0], V_1[1], V_1[2]);
    this->GetVertex(face_vertex[starting_index + 2])
      ->GetCoords(V_2[0], V_2[1], V_2[2]);

    // Compute the vectors vec_0:=(V_1 - V_0) and vec_1:=(V_2 - V_0). Those
    // vectors are linearly independent and lie within the face if the face is
    // not degenerated.
    std::array<double, 3> vec_0;
    std::array<double, 3> vec_1;
    for (int component_i = 0; component_i < 3; ++component_i)
    {
      vec_0[component_i] = V_1[component_i] - V_0[component_i];
      vec_1[component_i] = V_2[component_i] - V_0[component_i];
    }

    // Compute the cross product of vec_0 and vec_1 which is a normal vector to
    // the face defined by vec_0 and vec_1, i.e. the face of interest.
    outer_joint_normal[0] = vec_0[1] * vec_1[2] - vec_0[2] * vec_1[1];
    outer_joint_normal[1] = vec_0[2] * vec_1[0] - vec_0[0] * vec_1[2];
    outer_joint_normal[2] = vec_0[0] * vec_1[1] - vec_0[1] * vec_1[0];

    // Normalize this vector to get an unit normal vector, i.e. compute the
    // length of the vector and then dividing the entries by the length.
    auto length = std::sqrt(
        outer_joint_normal[0]*outer_joint_normal[0] +
        outer_joint_normal[1]*outer_joint_normal[1] +
        outer_joint_normal[2]*outer_joint_normal[2]
        );
    for (int dimension_i = 0; dimension_i < 3; ++dimension_i)
    {
      outer_joint_normal[dimension_i] /= length;
    }

    // Up to now it is not clear if the vector points to the outside or
    // the inside of the polygon (even though it is stored in
    // outer_joint_normal). This depends on the ordering of the vertices in the
    // array defining the face.
    // In ParMooN the shape descriptors of 3D shapes shall be defined such that
    // each face is numbered in counterclockwise ordering seen from the INSIDE
    // of the cell. Unfortunately there is no ShapeDesc test that checks this,
    // but this part of the code relies on that assumption. Under this
    // assumption the previously compute unit normal points to the inside of the
    // domain. Therefore, we have to turn it around.
    // Below is a code, that would check if vector has to be turned. I leave it
    // commented out for debugging purposes.
    for (int dimension_i = 0; dimension_i < 3; ++dimension_i)
    {
      outer_joint_normal[dimension_i] *= -1;
    }

    // // To check whether the vector points to the inside
    // // or the outside, another point V* is used, that is a vertex of the cell
    // // but NOT a vertex of the face of interest.  Then we consider the vector
    // // pointing from V_0 to V* and compute the dot product of this vector and
    // // the previously computed normal vector. For two vectors a and b, and the
    // // dot product it holds
    // // a dot b = C cos(theta)
    // // where C is a positive constant and theta is the angle between the
    // // vectors. Let n be the previously compute unit normal vector. Then it
    // // holds
    // // n dot (V* - V_0) > 0   if n points to the inside of the cell and
    // // n dot (V* - V_0) < 0   if n points to the outside of the cell.
    // // The case
    // // n dot (V* - V_0) = 0
    // // would occur if n and (V* - V_0) were orthogonal. But this is impossible
    // // due to construction of n and the choice of V*. This is checked and the
    // // vector is multiplied by -1 if n points to the inside of the cell.
    //
    // // Revert the previous flip
    // for (int dimension_i = 0; dimension_i < 3; ++dimension_i)
    // {
    //   outer_joint_normal[dimension_i] *= -1;
    // }

    // // First, let's find the viewpoint V*, i.e. a vertex of the cell that is not
    // // in the face.
    // std::array<double, 3> V_star;
    // int index_V_star = -1;
    // // Loop over all vertices in the cell and check if the vertex is also in
    // // FaceVertex[joint_index], i.e. a vertex of the face with number joint_index.
    // // We do this until we find a vertex that is not contained in the face.
    // for (int vertex_i = 0; vertex_i < this->GetN_Vertices(); ++vertex_i)
    // {
    //   for (int index_vertex_fc_j = 0; index_vertex_fc_j
    //       < face_vertex_len[joint_index]; ++index_vertex_fc_j)
    //   {
    //     // Check if the vertex is also a vertex of the facet. If not we found
    //     // our vertex.
    //     if (vertex_i == face_vertex[starting_index + index_vertex_fc_j])
    //     {
    //       // Here it is clear it is contained in the face. Therefore we can
    //       // continue with another vertex
    //       break;
    //     }
    //     else if ( index_vertex_fc_j == face_vertex_len[joint_index] - 1 )
    //     {
    //       // Here it is clear that the vertex is not contained in the face since
    //       // we arrived at the last vertex of the face.
    //       index_V_star = vertex_i;
    //     }
    //   }
    //   // If a vertex was found we don't have to check the remaining vertices.
    //   if (index_V_star >=0 )
    //   {
    //     break;
    //   }
    // }
    // // Check if an index was found and if so, get the coordinates of V*
    // if (index_V_star < 0)
    // {
    //   ErrThrow("No index of a vertex lying in the cell but not in the face was",
    //       " found.");
    // }
    // this->GetVertex(index_V_star)->GetCoords(V_star[0], V_star[1], V_star[2]);

    // // Check if the outer_unit_normal vector really points to outward direction.
    // // If this is not the case, then the dot product between this vector and
    // // V*-V_0 is > 0. In this case the unit normal vector has to be flipped
    // // around, i.e. multiplied by -1.
    // double points_inward = 0;
    // for (int dimension_i = 0; dimension_i < 3; ++dimension_i)
    // {
    //   points_inward += outer_jointnormal[dimension_i] * (V_star[dimension_i]
    //       - V_0[dimension_i]);
    // }
    // if (std::abs(points_inward) < 1e-10)
    // {
    //   ErrThrow("V*-V0 is orthogonal to the unit normal vector computed which ",
    //       "should have been forbidden.");
    // }
    // if (points_inward > 0 )
    // {
    //   for (int dimension_i = 0; dimension_i < 3; ++dimension_i)
    //   {
    //     outer_jointnormal[dimension_i] *= -1;
    //   }
    // }
    // else
    // {
    //   ErrThrow("Ohoh, unit normal points outward.")
    // }
  }
  else
  {
    ErrThrow("Only 2D or 3D are implemented.");
  }
  return outer_joint_normal;
}
