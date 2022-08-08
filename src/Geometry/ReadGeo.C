// =======================================================================
// @(#)ReadGeo.C        1.8 11/15/99
//
// Purpose:     read geometry file and
//              generate coarse grid
//
// Author:      Volker Behns  20.01.98
// Version:     1.0
//
// =======================================================================

#include <Database.h>
#include <Domain.h>
#include <MacroCell.h>
#include <JointEqN.h>
#include <PeriodicJoint.h>
#include <Quadrangle.h>
#include <MooNMD_Io.h>
#include <algorithm>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <functional>
#include <limits>
#include <cmath>

#include <InnerInterfaceJoint.h>
#ifdef __2D__
  #include <IsoBoundEdge.h>
  #include <IsoInterfaceJoint.h>
#else
  #include <BoundFace.h>
  #include <IsoBoundFace.h>
  #include <BdWall.h>
  #include <Hexahedron.h>
  #include <BdPlane.h>
  #include <IsoInterfaceJoint3D.h>
  #include <BdNoPRM.h>
  #include <BdCircle.h>
#endif

int TDomain::ReadGeo(std::istream& dat, bool readXgeo)
{
  char line[100];
  int i, j, N_Vertices, NVpF, NVE, NBCT;
  double *DCORVG;
  int *KVERT, *KNPR;
  #ifdef __3D__
    int NBF, *BoundFaces, *FaceParam;
    int *InterfaceParam, N_Interfaces;
  #endif

  // physical references
  int *ELEMSREF;
  
  dat.getline (line, 99);
  dat.getline (line, 99);
  
  // determine dimensions for creating arrays
  dat >> N_RootCells >> N_Vertices >> NVpF >> NVE >> NBCT;
  dat.getline (line, 99);
  dat.getline (line, 99);

  // allocate auxillary fields
  #ifdef __2D__
    DCORVG =  new double[2*N_Vertices];
  #else
    DCORVG =  new double[3*N_Vertices];
  #endif
  KVERT = new int[NVE*N_RootCells];
  KNPR = new int[N_Vertices];
  
  // additional array of physical properties of elements (only if readXgeo)
  ELEMSREF = new int[N_RootCells];

  
  // read fields
  for (i=0;i<N_Vertices;i++)
  {
    #ifdef __2D__
      dat >> DCORVG[2*i] >> DCORVG[2*i + 1];
    #else
      dat >> DCORVG[3*i] >> DCORVG[3*i + 1] >> DCORVG[3*i + 2];
    #endif
    dat.getline (line, 99);
  }

  dat.getline (line, 99);
   
  for (i=0;i<N_RootCells;i++)
  {
    for (j=0;j<NVE;j++)
      dat >> KVERT[NVE*i + j];
    if(readXgeo)
      dat >> ELEMSREF[i];
    else
      ELEMSREF[i] = 0;
    dat.getline (line, 99);
  }

  dat.getline (line, 99);
  
  for (i=0;i<N_Vertices;i++)
    dat >> KNPR[i];

  #ifdef __3D__
    dat.getline (line, 99);
    dat.getline (line, 99);

    NBF = NBCT;
    BoundFaces = new int[NBF*NVpF];
    FaceParam = new int[4*NBF];

    for (i=0;i<NBF;i++)
    {
      for (j=0;j<NVpF;j++)
        dat >> BoundFaces[i*NVpF+j];
      dat.getline (line, 99);
    }

    dat.getline (line, 99);

    for (i=0;i<NBF;i++)
    {
	dat >> FaceParam[i*4] >> FaceParam[i*4+1] >>
	    FaceParam[i*4+2] >>  FaceParam[i*4+3];
	dat.getline (line, 99);
    }

    dat.getline (line, 99);
    N_Interfaces = -12345;
    if (dat.eof())
      N_Interfaces = 0;
    else
      dat >> N_Interfaces;

    if(N_Interfaces == -12345)
      N_Interfaces = 0;

    dat.getline (line, 99);

    if(N_Interfaces)
    {
      InterfaceParam = new int[N_Interfaces*6];

      dat.getline (line, 99);
      for(i=0;i<N_Interfaces;i++)
      {
        dat >> InterfaceParam[i*6]   >> InterfaceParam[i*6+1] >>
               InterfaceParam[i*6+2] >> InterfaceParam[i*6+3] >>
               InterfaceParam[i*6+4] >> InterfaceParam[i*6+5];
        dat.getline (line, 99);
      } // endfor i
    } // endif N_Interfaces
    else
      InterfaceParam = nullptr;

  #endif

  #ifdef __2D__
    MakeGrid(DCORVG, KVERT, KNPR, ELEMSREF, N_Vertices, NVE);
  #else
    MakeGrid(DCORVG, KVERT, KNPR, ELEMSREF, N_Vertices, NVE,
             BoundFaces, FaceParam, NBF, NVpF,
             InterfaceParam, N_Interfaces);
  delete [] BoundFaces;
  delete [] FaceParam;
  #endif

  delete [] DCORVG;
  delete [] KVERT;
  delete [] KNPR;
  delete [] ELEMSREF;

  return 0;
}

#ifdef __2D__
// extended version
// Alfonso, 2011
int TDomain::MakeGrid(double *DCORVG, int *KVERT, int *KNPR, int *ELEMSREF,
                      int N_Vertices, int NVE)
{
  int a, b, j, k, l, comp, Part, NeighborID, N_Edges,maxElpV = 0;
  int aux1, aux2, aux3;
  double T_a, T_b, T, X, Y;
  double Xmin = 1e10, Xmax = -1e10, Ymin = 1e10, Ymax = -1e10;
  int *KVEL;
  TVertex **NewVertices, *LocVerts[4];
  TJoint **KMT, *Joint;
  TBaseCell *JNeib1, *JNeib2;
  Shapes CellType;


  double bd_parameter_a, bd_parameter_b;

  Output::info<4>("Domain::MakeGrid()", "Creating 2D grid.");
 
  // generate vertices, edges and cells
  // search neighbours
  KVEL = new int[N_Vertices];
  
  memset(KVEL, 0, N_Vertices * sizeof(int));
  
  //KVEL(j) = how many times vertices j appears in the list
  // Note: KVERT[i] can be = 0, in case of mixed (tria+quad meshes). In this
  // case, the vertex does not exist.
  for (int i=0;i<NVE*N_RootCells;i++)
  {
    if (KVERT[i])
    {
      KVEL[KVERT[i]-1]++;
    }
  }
  
  // get the maximum
  for (int i=0;i<N_Vertices;i++)
    if (KVEL[i] > maxElpV) maxElpV = KVEL[i];
  
  delete [] KVEL;
  //KVEL = new int[++maxElpV * N_Vertices];
  maxElpV = maxElpV+1;
  KVEL = new int[maxElpV * N_Vertices];
  memset(KVEL, 0, maxElpV * N_Vertices * sizeof(int));
  // first column contains the number of following elements
  for (int i=0;i<NVE*N_RootCells;i++)
  {
    if (KVERT[i])
    {
      j = (KVERT[i] - 1)*maxElpV;
      KVEL[j]++;
      KVEL[j + KVEL[j]] = i / NVE;
    }
  }
  
  // generate vertices
  //cout << " Domain::MakeGrid() generate " << N_Vertices << " vertices " << endl;
  NewVertices = new TVertex*[N_Vertices];
  
  for (int i=0;i<N_Vertices;i++)
  {
    //cout << " i = " << i << endl;
    if (KNPR[i])
    {
      // vertices on a boundary (inner or outer) joint
      // described by bd parametrization of 
      // part BdParts[ KNPR[i] ] [ (int) DCORVG[2*i]] ;

      T = DCORVG[2*i];
      comp = (int) T;
      if (GetLastLocalComp(KNPR[i]-1) == comp - 1)
      {
        comp--;
        T = 1.0;
      }
      else
        T -= comp;

      BdParts[KNPR[i] - 1]->GetXYofT(comp, T, X, Y);

      if (X > Xmax) Xmax = X;
      if (X < Xmin) Xmin = X;
      if (Y > Ymax) Ymax = Y;
      if (Y < Ymin) Ymin = Y;
  
      NewVertices[i] = new TVertex(X, Y);
      //Output::print(" vertex: ", i, " comp: ", comp, " X,Y= ", X, ",", Y);
      //Output::print("bd ", i, " ", X, " ", Y);
    }
    else
    { 
      // inner vertex (described by coordinates)
      NewVertices[i] = new TVertex(DCORVG[2*i], DCORVG[2*i+1]);
    } 
  } // for (i=0;i<N_Vertices;i++) {
  // set bounding box
  StartX = Xmin;
  StartY = Ymin;
  BoundX = Xmax - Xmin;
  BoundY = Ymax - Ymin;
  //StartZ = 0; // this is 2D
  //BoundZ = 0; // this is 2D
  
  // create the CellTree and set references
  CellTree = new TBaseCell*[N_RootCells];
  for (int i=0;i<N_RootCells;i++)
  {
    CellType = Quadrangle;
    if (NVE == 3)
      CellType = Triangle;
    else
      if (!KVERT[NVE*i + 3]) CellType = Triangle;

    if (CellType == Quadrangle)
    {
      LocVerts[0] = NewVertices[KVERT[NVE*i    ]-1];
      LocVerts[1] = NewVertices[KVERT[NVE*i + 1]-1];
      LocVerts[2] = NewVertices[KVERT[NVE*i + 2]-1];
      LocVerts[3] = NewVertices[KVERT[NVE*i + 3]-1];

      CellType = ((TQuadrangle *) TDatabase::RefDescDB[Quadrangle]->
                 GetShapeDesc())->CheckQuad(LocVerts);

      CellTree[i] = new TMacroCell(TDatabase::RefDescDB[CellType],
                                   RefLevel);

      CellTree[i]->SetVertex(0, LocVerts[0]);
      CellTree[i]->SetVertex(1, LocVerts[1]);
      CellTree[i]->SetVertex(2, LocVerts[2]);
      CellTree[i]->SetVertex(3, LocVerts[3]);

    }
    else
    {
      CellTree[i] = new TMacroCell(TDatabase::RefDescDB[
                                   Triangle], RefLevel);
      CellTree[i]->SetVertex(0, NewVertices[KVERT[NVE*i    ]-1]);
      CellTree[i]->SetVertex(1, NewVertices[KVERT[NVE*i + 1]-1]);
      CellTree[i]->SetVertex(2, NewVertices[KVERT[NVE*i + 2]-1]);
    }
    // ReferenceID: used to specify physical/geometrical reference
    CellTree[i]->SetReference_ID(ELEMSREF[i]);
    // set the index in the cell list
    CellTree[i]->SetCellIndex(i);
  }
  
  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);


  // generate edges
  KMT = new TJoint*[N_RootCells*4];
  for (int i=0; i<N_RootCells*4; i++)
    KMT[i] = nullptr;

  for (int i=0;i<N_RootCells;i++)
  {
    N_Edges = CellTree[i]->GetN_Edges();
 
    for (j=0;j<N_Edges;j++)
    {
      a = KVERT[NVE*i + j] - 1;
      b = KVERT[NVE*i + (j+1) % N_Edges] - 1;
      // Part: -1 for inner, 0,1,2,... ofr bdpart 0,1,2,...
      Part = KNPR[a] - 1; // Bdpart
      // find neighbor ID
      NeighborID = -1;
      aux1 = KVEL[a*maxElpV];
      aux2 = KVEL[b*maxElpV];

      for (k=1;k<=aux1;k++)
      {
        aux3 = KVEL[a*maxElpV + k];
        if (aux3 == i) continue;

        for (l=1;l<=aux2;l++)
          if (aux3 == KVEL[b*maxElpV + l])
          {
            NeighborID = aux3;
            break;
          }
        if (NeighborID >= 0) break;
      }
      
      //Output::print(NeighborID);
      if (NeighborID > i)
      {
        if( (KNPR[a]) && (KNPR[b]) )
        {
          if( (KNPR[a] == KNPR[b]) && (Interfaces[KNPR[a]-1] < 0) )
          {
            KMT[i*4 + j] = new TJointEqN(CellTree[i], CellTree[NeighborID]);
            // bd.comp ID of the point is given as integer part of first coordinate
            comp = (int) DCORVG[2*a]; 
            Part = KNPR[a] - 1;
            bd_parameter_a = DCORVG[2*a] - comp;
            // correction for the endpoint of nonclosed interfaces
            if(BdParts[Part]->GetN_BdComps() == comp)
            {
              bd_parameter_a = 1.0;
              comp--;
            }
      
            // check if vertex belongs to boundary component
            if(BdParts[Part]->GetBdComp(comp)->GetTofXY(NewVertices[a]->GetX(), 
                                                        NewVertices[a]->GetY(),
                                                        bd_parameter_a))
            {
              // if not: 
              if (comp)
                comp--;
              else
                comp = GetLastLocalComp(Part);

              int check_bd_node = 
                BdParts[Part]->GetBdComp(comp)->GetTofXY(NewVertices[a]->GetX(),
                                                         NewVertices[a]->GetY(),
                                                         bd_parameter_a);
              if (!check_bd_node)
              {
                cout << " MakeGrid(), line " << __LINE__ << ": Error: vertex " 
                     << a << " does not belong to a bd. component " << endl;
              }
            }

            comp = (int) DCORVG[2*b];
            Part = KNPR[b] - 1;

            if(comp > GetLastLocalComp(Part) ) 
              comp = GetLastLocalComp(Part);

            bd_parameter_b = DCORVG[2*b] - comp;
            if(BdParts[Part]->GetBdComp(comp)->GetTofXY(NewVertices[b]->GetX(), 
                                                        NewVertices[b]->GetY(),
                                                        bd_parameter_b))
            {
              if (comp)
                comp--;
              else
                comp = GetLastLocalComp(Part);
              
              int check_bd_node = BdParts[Part]->GetBdComp(comp)->GetTofXY(
                NewVertices[b]->GetX(), NewVertices[b]->GetY(), bd_parameter_b);
              if (!check_bd_node)
              {
                cout << " MakeGrid(), line " << __LINE__ << ": Error: vertex " 
                     << b << " does not belong to a bd. component " << endl;
              }
            }
      
            T_a = bd_parameter_a;
            T_b = bd_parameter_b;
      
            if (BdParts[Part]->GetBdComp(comp)->GetType() != Line)
            {
              if (std::abs(T_a) < 1e-4 || std::abs(T_a) > 0.9999 ||
                  std::abs(T_b) < 1e-4 || std::abs(T_b) > 0.9999)
              {
                X = (NewVertices[a]->GetX() + NewVertices[b]->GetX()) / 2.;
                Y = (NewVertices[a]->GetY() + NewVertices[b]->GetY()) / 2.;

                BdParts[Part]->GetBdComp(comp)->GetTofXY(X, Y, T);

                if ((T_a - T)*(T - T_b) < 0)
                {
                  cout << __FILE__ << " line: " << __LINE__ 
                       << " changing boundary parameters of vertices " 
                       << a << " and " << b << endl;

                  if (std::abs(T_a) < 1e-4)
                    T_a = 1.0;
                  else
                    if (std::abs(T_a) > 0.9999)
                      T_a = 0.0;

                  if (std::abs(T_b) < 1e-4)
                    T_b = 1.0;
                  else
                    if (std::abs(T_b) > 0.9999)
                      T_b = 0.0;
                }
              }
            }

            if(BdParts[Part]->GetBdComp(comp)->IsFreeBoundary())
              KMT[i*4 + j] = new TIsoInterfaceJoint(BdParts[Part]->
                 GetBdComp(comp), T_a, T_b,CellTree[i], CellTree[NeighborID]);
            else
              KMT[i*4 + j] = new TInterfaceJoint(BdParts[Part]->
                 GetBdComp(comp), T_a, T_b,CellTree[i], CellTree[NeighborID]);

            CellTree[i]->SetJoint(j, KMT[i*4 + j]);
            ((TInterfaceJoint *)KMT[i*4 + j])->CheckOrientation();
          } // two vertices on the same Interfaces
          else
          {
            // two vertices on the boundary
            //KMT[i*4 + j] = new TJointEqN(CellTree[i], CellTree[NeighborID]);
            
            // check if neighbor belongs to different subdomain:
            if (CellTree[i]->GetReference_ID() 
              == CellTree[NeighborID]->GetReference_ID())
            {
              KMT[i*4 + j] = new TJointEqN(CellTree[i], CellTree[NeighborID]);
            }
            else
            {
              //KMT[i*4+j] = new TJointEqN(CellTree[i], CellTree[NeighborID]);
              KMT[i*4+j] = new TInnerInterfaceJoint(CellTree[i], 
                                                    CellTree[NeighborID]);
              ((TInnerInterfaceJoint*)KMT[i*4 + j])->SetParams(
                 NewVertices[a]->GetX(),
                 NewVertices[a]->GetY(),
                 NewVertices[b]->GetX()-NewVertices[a]->GetX(),
                 NewVertices[b]->GetY()-NewVertices[a]->GetY());
              ((TInnerInterfaceJoint*)KMT[i*4 + j])->SetIndexInNeighbor(
                 CellTree[i],j);
            }
          }
        }
        else
        {
          // at least one vertex is inside the domain
          // check if neighbor belongs to different subdomain:
          if(CellTree[i]->GetReference_ID() 
             == CellTree[NeighborID]->GetReference_ID())
          {
            KMT[i*4 + j] = new TJointEqN(CellTree[i], CellTree[NeighborID]);
          }
          else
          {
            //KMT[i*4+j] = new TJointEqN(CellTree[i], CellTree[NeighborID]);
            KMT[i*4+j] = new TInnerInterfaceJoint(CellTree[i], 
                                                  CellTree[NeighborID]);
            ((TInnerInterfaceJoint*)KMT[i*4 + j])->SetParams(
                   NewVertices[a]->GetX(),
                   NewVertices[a]->GetY(),
                   NewVertices[b]->GetX()-NewVertices[a]->GetX(),
                   NewVertices[b]->GetY()-NewVertices[a]->GetY());
            ((TInnerInterfaceJoint*)KMT[i*4 + j])->SetIndexInNeighbor(
                   CellTree[i],j);
          }
        }
        CellTree[i]->SetJoint(j, KMT[i*4 + j]);
      }
      else 
      {
        if (NeighborID == -1)
        {
          if (Interfaces[KNPR[a]-1] < 0)
          {
            comp = (int) DCORVG[2*b];
            T_b = DCORVG[2*b] - comp;
            Part = KNPR[b] - 1;
            
            if (BdParts[Part]->GetBdComp(comp)->GetTofXY(NewVertices[a]->GetX(), 
                                                         NewVertices[a]->GetY(), 
                                                         T_a))
            {
              if (comp)
                comp--;
              else
                comp = GetLastLocalComp(Part);
        
              BdParts[Part]->GetBdComp(comp)->GetTofXY(NewVertices[a]->GetX(), 
                                                       NewVertices[a]->GetY(), 
                                                       T_a);
            }
          }
          else if (Interfaces[KNPR[b]-1] < 0)
          {
            comp = (int) DCORVG[2*a];
            T_a = DCORVG[2*a] - comp;
            
            if (BdParts[Part]->GetBdComp(comp)->GetTofXY(NewVertices[b]->GetX(), 
                                                         NewVertices[b]->GetY(), 
                                                         T_b))
            {
              if (comp)
                comp--;
              else
                comp = GetLastLocalComp(Part);
              
              BdParts[Part]->GetBdComp(comp)->GetTofXY(NewVertices[b]->GetX(), 
                                                       NewVertices[b]->GetY(), 
                                                       T_b);
            }
          }
          else
          {
            //cout << "a=" << a << endl;
            comp = (int) DCORVG[2*a];
            T_a = DCORVG[2*a] - comp;
            T_b = DCORVG[2*b] - comp;
            //cout << "comp: " << comp << endl;
          }
          //cout << "Part = " << Part << endl;
          //cout << " test = " << BdParts[Part]->GetBdComp(comp)->IsFreeBoundary() << endl;
          if (T_b < T_a) 
            T_b = 1.;
          if(BdParts[Part]->GetBdComp(comp)->IsFreeBoundary())
            Joint = new TIsoBoundEdge(BdParts[Part]->GetBdComp(comp), T_a, T_b);
          else 
          {
            Joint = new TBoundEdge(BdParts[Part]->GetBdComp(comp), T_a, T_b);
          }
          CellTree[i]->SetJoint(j, Joint);
        }
        else
        {
          // joint already created (from the neighboring cell)
          TJoint *nJoint;
          for (int jj=0;jj<CellTree[NeighborID]->GetN_Edges();jj++)
          {
            nJoint = CellTree[NeighborID]->GetJoint(jj);
            JNeib1 = nJoint->GetNeighbour(0);
            JNeib2 = nJoint->GetNeighbour(1);
            if((JNeib1 == CellTree[NeighborID] && JNeib2 == CellTree[i]) 
               || (JNeib1 == CellTree[i] && JNeib2 == CellTree[NeighborID])) 
              break;
          }
          CellTree[i]->SetJoint(j, nJoint);
          
          if (CellTree[i]->GetReference_ID()
            != CellTree[NeighborID]->GetReference_ID())
          {
            ((TInnerInterfaceJoint*)nJoint)->SetIndexInNeighbor(CellTree[i],j);
          }
        }
      }
    }
  } //for (i=0;i<N_RootCells;i++) {
  
  // free memory
  delete [] KVEL;
  delete [] NewVertices;
  delete [] KMT;
  
  return 0;
}

int TDomain::MakeGrid(double *DCORVG, int *KVERT, int *KNPR, int N_Vertices, 
                      int NVE)
{
   int a, b, i, j, k, l, comp, Part, Neib, N_E, maxElpV = 0;
  int aux1, aux2, aux3;
  double T_a, T_b, T, X, Y;
  double Xmin = 1e10, Xmax = -1e10, Ymin = 1e10, Ymax = -1e10;
  int *KVEL;
  TVertex **NewVertices, *LocVerts[4];
  TJoint **KMT, *Joint;
  TBaseCell *JNeib1, *JNeib2;
  Shapes CellType;
 
  // generate vertices, edges and cells
  // search neighbours
  KVEL = new int[N_Vertices];

  memset(KVEL, 0, N_Vertices * sizeof(int));
 
  for (i=0;i<NVE*N_RootCells;i++)
    if (KVERT[i]) KVEL[KVERT[i]-1]++;

   for (i=0;i<N_Vertices;i++)
    if (KVEL[i] > maxElpV) maxElpV = KVEL[i];

  delete [] KVEL;
  KVEL = new int[++maxElpV * N_Vertices];

  memset(KVEL, 0, maxElpV * N_Vertices * sizeof(int));
  
  // first column contains the number of following elements
   for (i=0;i<NVE*N_RootCells;i++)
    if (KVERT[i])
    {
      j = (KVERT[i] - 1)*maxElpV;
      KVEL[j]++;
      KVEL[j + KVEL[j]] = i / NVE;
    }

  // generate vertices
    NewVertices = new TVertex*[N_Vertices];
  
  for (i=0;i<N_Vertices;i++)
    if (KNPR[i])
    {
      T = DCORVG[2*i];
      comp = (int) T;
      if (GetLastLocalComp(KNPR[i]-1) == comp - 1)
      {
        comp--;
        T = 1.0;
      }
      else
        T -= comp;
      
      BdParts[KNPR[i] - 1]->GetXYofT(comp, T, X, Y);

      if (X > Xmax) Xmax = X;
      if (X < Xmin) Xmin = X;
      if (Y > Ymax) Ymax = Y;
      if (Y < Ymin) Ymin = Y;

      NewVertices[i] = new TVertex(X, Y);
      //Output::print("bd ", i, " ", X, " ", Y);
    }
    else
      NewVertices[i] = new TVertex(DCORVG[2*i], DCORVG[2*i+1]);

   // set bounding box
  StartX = Xmin;
  StartY = Ymin;
  BoundX = Xmax - Xmin;
  BoundY = Ymax - Ymin;

  // generate cells
  CellTree = new TBaseCell*[N_RootCells];

  for (i=0;i<N_RootCells;i++)
  {
    CellType = Quadrangle;
    if (NVE == 3)
      CellType = Triangle;
    else
      if (!KVERT[NVE*i + 3]) CellType = Triangle;

    if (CellType == Quadrangle)
    {
      LocVerts[0] = NewVertices[KVERT[NVE*i    ]-1];
      LocVerts[1] = NewVertices[KVERT[NVE*i + 1]-1];
      LocVerts[2] = NewVertices[KVERT[NVE*i + 2]-1];
      LocVerts[3] = NewVertices[KVERT[NVE*i + 3]-1];

      CellType = ((TQuadrangle *) TDatabase::RefDescDB[Quadrangle]->
                 GetShapeDesc())->CheckQuad(LocVerts);

      CellTree[i] = new TMacroCell(TDatabase::RefDescDB[CellType],
                                   RefLevel);

      CellTree[i]->SetVertex(0, LocVerts[0]);
      CellTree[i]->SetVertex(1, LocVerts[1]);
      CellTree[i]->SetVertex(2, LocVerts[2]);
      CellTree[i]->SetVertex(3, LocVerts[3]);

    }
    else
    {
      CellTree[i] = new TMacroCell(TDatabase::RefDescDB[
                                   Triangle], RefLevel);
      CellTree[i]->SetVertex(0, NewVertices[KVERT[NVE*i    ]-1]);
      CellTree[i]->SetVertex(1, NewVertices[KVERT[NVE*i + 1]-1]);
      CellTree[i]->SetVertex(2, NewVertices[KVERT[NVE*i + 2]-1]);
    }
  }

  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);

  // generate edges
  KMT = new TJoint*[N_RootCells*4];
  for (i=0;i<N_RootCells*4;i++)
    KMT[i] = nullptr;
  
  for (i=0;i<N_RootCells;i++)
  {
    switch (CellTree[i]->GetType())
    {
      case Triangle: 
           N_E = 3;
           break;
      case Parallelogram: 
      case Quadrangle: 
      case Rectangle: 
           N_E = 4;
           break;
      default:
           cerr << "Error in ReadGeo" << endl;
           return -1;
    }

    //Output::print("edges ", i, " ", N_RootCells);
   for (j=0;j<N_E;j++)
    {
      a = KVERT[NVE*i + j] - 1;
      b = KVERT[NVE*i + (j+1) % N_E] - 1;
      Part = KNPR[a] - 1;
      Neib = -1;

      aux1 = KVEL[a*maxElpV];
      aux2 = KVEL[b*maxElpV];

      for (k=1;k<=aux1;k++)
      {
        aux3 = KVEL[a*maxElpV + k];
        if (aux3 == i) continue;

        for (l=1;l<=aux2;l++)
          if (aux3 == KVEL[b*maxElpV + l])
          {
            Neib = aux3;
            break;
          }
        if (Neib >= 0) break;
      }
      if (Neib > i)
      {
        if( (KNPR[a]) && (KNPR[b]) )
        {
          if( (KNPR[a] == KNPR[b]) && (Interfaces[KNPR[a]-1] < 0) )
          {
            KMT[i*4 + j] = new TJointEqN(CellTree[i], CellTree[Neib]);
            comp = (int) DCORVG[2*a];
            Part = KNPR[a] - 1;
            T_a = DCORVG[2*a] - comp;

            // correction for the endpoint of nonclosed interfaces
            if(BdParts[Part]->GetN_BdComps() == comp)
            {
              T_a = 1.0;
              comp--;
            }

            if (BdParts[Part]->GetBdComp(comp)->GetTofXY(
                  NewVertices[a]->GetX(), NewVertices[a]->GetY(), T_a))
            {
              if (comp)
                comp--;
              else
                comp = GetLastLocalComp(Part);

              BdParts[Part]->GetBdComp(comp)->GetTofXY(
                NewVertices[a]->GetX(), NewVertices[a]->GetY(), T_a);
            }

            comp = (int) DCORVG[2*b];
            Part = KNPR[b] - 1;

            // correction needed for not closed interfaces
            if(comp > GetLastLocalComp(Part) ) comp = GetLastLocalComp(Part);

            T_b = DCORVG[2*b] - comp;
            if (BdParts[Part]->GetBdComp(comp)->GetTofXY(
                  NewVertices[b]->GetX(), NewVertices[b]->GetY(), T_b))
            {
              if (comp)
                comp--;
              else
                comp = GetLastLocalComp(Part);

              BdParts[Part]->GetBdComp(comp)->GetTofXY(
                NewVertices[b]->GetX(), NewVertices[b]->GetY(), T_a);
            }

            if (BdParts[Part]->GetBdComp(comp)->GetType() != Line)
              if (std::abs(T_a) < 1e-4 || std::abs(T_a) > 0.9999 ||
                  std::abs(T_b) < 1e-4 || std::abs(T_b) > 0.9999)
              {
                X = (NewVertices[a]->GetX() + NewVertices[b]->GetX()) / 2;
                Y = (NewVertices[a]->GetY() + NewVertices[b]->GetY()) / 2;

                BdParts[Part]->GetBdComp(comp)->GetTofXY(X, Y, T);

                if ((T_a - T)*(T - T_b) < 0)
                {
                  if (std::abs(T_a) < 1e-4)
                    T_a = 1.0;
                  else
                    if (std::abs(T_a) > 0.9999)
                      T_a = 0.0;

                  if (std::abs(T_b) < 1e-4)
                    T_b = 1.0;
                  else
                    if (std::abs(T_b) > 0.9999)
                      T_b = 0.0;
                }
              }

            if(BdParts[Part]->GetBdComp(comp)->IsFreeBoundary())
              KMT[i*4 + j] = new TIsoInterfaceJoint(BdParts[Part]->
                 GetBdComp(comp), T_a, T_b,CellTree[i], CellTree[Neib]);
            else
              KMT[i*4 + j] = new TInterfaceJoint(BdParts[Part]->
                 GetBdComp(comp), T_a, T_b,CellTree[i], CellTree[Neib]);

            CellTree[i]->SetJoint(j, KMT[i*4 + j]);
            ((TInterfaceJoint *)KMT[i*4 + j])->CheckOrientation();
          }
          else
          {
            KMT[i*4 + j] = new TJointEqN(CellTree[i], CellTree[Neib]);
          }
        }
        else
        {
          KMT[i*4 + j] = new TJointEqN(CellTree[i], CellTree[Neib]);
        }
        CellTree[i]->SetJoint(j, KMT[i*4 + j]);
      }
      else
        if (Neib == -1)
        {
          if (Interfaces[KNPR[a]-1] < 0)
          {
	      comp = (int) DCORVG[2*b];
            T_b = DCORVG[2*b] - comp;
            Part = KNPR[b] - 1;

            if (BdParts[Part]->GetBdComp(comp)->GetTofXY(
                  NewVertices[a]->GetX(), NewVertices[a]->GetY(), T_a))
            {
              if (comp)
                comp--;
              else
                comp = GetLastLocalComp(Part);

              BdParts[Part]->GetBdComp(comp)->GetTofXY(
                NewVertices[a]->GetX(), NewVertices[a]->GetY(), T_a);
            }
          }
          else
            if (Interfaces[KNPR[b]-1] < 0)
            {
		comp = (int) DCORVG[2*a];
              T_a = DCORVG[2*a] - comp;

              if (BdParts[Part]->GetBdComp(comp)->GetTofXY(
                    NewVertices[b]->GetX(), NewVertices[b]->GetY(), T_b))
              {
                if (comp)
                  comp--;
                else
                  comp = GetLastLocalComp(Part);

                BdParts[Part]->GetBdComp(comp)->GetTofXY(
                  NewVertices[b]->GetX(), NewVertices[b]->GetY(), T_b);
              }
            }
            else
            {
              comp = (int) DCORVG[2*a];
	      T_a = DCORVG[2*a] - comp;
              T_b = DCORVG[2*b] - comp;
            }

          if (T_b < T_a) T_b = 1.;
          if(BdParts[Part]->GetBdComp(comp)->IsFreeBoundary())
            Joint = new TIsoBoundEdge(BdParts[Part]->
                 GetBdComp(comp), T_a, T_b);
          else
            Joint = new TBoundEdge(BdParts[Part]->
                 GetBdComp(comp), T_a, T_b);
          CellTree[i]->SetJoint(j, Joint);
        }
        else
        {
          for (k=0;k<4;k++)
          {
            Joint = KMT[Neib*4 + k];
            if (Joint)
            {
              aux1 = Neib*4 + k;
              JNeib1 = KMT[aux1]->GetNeighbour(0);
              JNeib2 = KMT[aux1]->GetNeighbour(1);

              if ((JNeib1 == CellTree[Neib] && JNeib2 == CellTree[i]) ||
                  (JNeib1 == CellTree[i] && JNeib2 == CellTree[Neib]) ){ break;}
            }
          }

          CellTree[i]->SetJoint(j, KMT[Neib*4 + k]);
        }
    }
   //Output::print("edges done ", i, " ", N_RootCells);
  }

   // free memory
  delete [] KVEL;
  delete [] NewVertices;
  delete [] KMT;
  
  return 0;
}

#else


// this makes a 3D coarse grid
int TDomain::MakeGrid(double *DCORVG, int *KVERT, int *KNPR, int *ELEMSREF,
                      int N_Vertices, int NVE, int *,
                      int *FaceParam, int NBF, int,
                      int *InterfaceParam, int N_Interfaces)
{
  int i, j, k, l=0, m, maxElpV = 0, comp, auxi, auxj;
  double S, T, X, Y, Z;
  double Xmin = 1e10, Xmax = -1e10, Ymin = 1e10, Ymax = -1e10;
  double Zmin = 1e10, Zmax = -1e10;
  int *KVEL, *NewJoints_aux1, *NewJoints_aux2, NewJoints_aux3;
  int aux1, aux2, aux3;
  TVertex **NewVertices, *CurrVert, *Vert[4];
  std::vector<TJoint *> NewJoints;
  Shapes CellType;
  TBaseCell *CurrCell, *CurrCell_aux;
  int maxlen, maxlen_aux, N_Faces, N_Verts, N_Faces_aux;
  const int *TmpFV, *TmpFV_aux, *TmpLen, *TmpLen_aux;

  TBoundFace *bdface;
  const TBoundComp3D *bdcomp;
  double Param1[4], Param2[4];
  JointType CurrJointType;
//  TInterfaceJoint3D *IFJoint;
//  TIsoInterfaceJoint3D *IIJoint;

  // generate vertices, faces and cells
  // search neighbours
  KVEL = new int[N_Vertices];

  memset(KVEL, 0, N_Vertices * sizeof(int));

  for (i=0;i<NVE*N_RootCells;i++)
    if (KVERT[i]) KVEL[KVERT[i]-1]++;

  for (i=0;i<N_Vertices;i++)
    if (KVEL[i] > maxElpV) maxElpV = KVEL[i];

  delete[] KVEL;
  KVEL = new int[++maxElpV * N_Vertices];

  memset(KVEL, 0, maxElpV * N_Vertices * sizeof(int));
  
  // first colomn contains the number of following elements
  for (i=0;i<NVE*N_RootCells;i++)
    if (KVERT[i])
    {
      j = (KVERT[i] - 1)*maxElpV;
      KVEL[j]++;
      KVEL[j + KVEL[j]] = i / NVE;
    }

  // generate vertices
  NewVertices = new TVertex*[N_Vertices];

  for (i=0;i<N_Vertices;i++)
    if (KNPR[i])
    {
	// parametrisation of the boundary is given
	if (KNPR[i]>0)
	{
	    // component of boundary
	    comp = (int) DCORVG[3*i];
	    // value of first parameter
	    T = DCORVG[3*i + 1];
	    // value of second parameter
	    S = DCORVG[3*i + 2];
	    // get coordinates from parameters
	    BdParts[KNPR[i] - 1]->GetXYZofTS(comp, T, S, X, Y, Z);
	}
	else
	{
	    X = DCORVG[3*i];
	    Y = DCORVG[3*i+1];
	    Z = DCORVG[3*i+2];
	    //cout << "read " << X << " " << Y << " " << Z << endl;
	}
      if (X > Xmax) Xmax = X;
      if (X < Xmin) Xmin = X;
      if (Y > Ymax) Ymax = Y;
      if (Y < Ymin) Ymin = Y;
      if (Z > Zmax) Zmax = Z;
      if (Z < Zmin) Zmin = Z;

      NewVertices[i] = new TVertex(X, Y, Z);
    }
    else
      NewVertices[i] = new TVertex(DCORVG[3*i], DCORVG[3*i + 1],
                                   DCORVG[3*i + 2]);

  // set bounding box
  StartX = Xmin;
  StartY = Ymin;
  StartZ = Zmin;
  BoundX = Xmax - Xmin;
  BoundY = Ymax - Ymin;
  BoundZ = Zmax - Zmin;

  // generate cells
  CellTree = new TBaseCell*[N_RootCells];

  for (i=0;i<N_RootCells;i++)
  {
    CellType = Hexahedron;
    if (NVE == 4)
      CellType = Tetrahedron;
    else
      if (!KVERT[NVE*i + 7]) CellType = Tetrahedron;

    if (CellType == Hexahedron)
    {
      CellTree[i] = new TMacroCell(TDatabase::RefDescDB[Hexahedron], RefLevel);

      auxi = NVE*i;
      for (j=0;j<8;j++)
        CellTree[i]->SetVertex(j, NewVertices[KVERT[auxi++]-1]);

      // check if Hexahedron is even a brick and if yes, change the refinement
      // descriptor accordingly. This allows us to use an affine reference 
      // mapping instead of the (more expensive) trilinear one.
      CellType = ((THexahedron *) TDatabase::RefDescDB[Hexahedron]->
                 GetShapeDesc())->CheckHexa(((TGridCell *)(CellTree[i]))
                        ->GetVertices());
      if(CellType != Hexahedron)
      {
        int ret = CellTree[i]->SetRefDesc(TDatabase::RefDescDB[CellType]);
        if(ret==-1)
        {
          ErrThrow("setRefDesc(", CellType, ") was not sucessfull");
        }
      }
    }
    else
    {
      CellTree[i] = new TMacroCell(TDatabase::RefDescDB[Tetrahedron],
                                   RefLevel);

      auxi = NVE*i;
      for (j=0;j<4;j++)
        CellTree[i]->SetVertex(j, NewVertices[KVERT[auxi++]-1]);
    }
  }

  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);

  // generate faces
  NewJoints.resize(N_RootCells*6, nullptr);

  //Output::print("NBF a ", NBF);
  // first generate boundary joints
  // NBF -- number of faces on boundaries
  for (i=0;i<NBF;i++)
  {
    auxi = i*4;
    //cout << FaceParam[auxi + 2] <<  endl;
    bdcomp = BdParts[FaceParam[auxi + 2] - 1]->GetBdComp(FaceParam[auxi+3] - 1);
    //cout << "free " << bdcomp->IsFreeBoundary()  << endl;
    if(bdcomp->IsFreeBoundary())
    {
      NewJoints.at(FaceParam[auxi]*6 + FaceParam[auxi + 1] - 7) =
        new TIsoBoundFace(bdcomp);
    }
    else
    {
      NewJoints.at(FaceParam[auxi]*6 + FaceParam[auxi + 1] - 7) =
        new TBoundFace(bdcomp);
    }
  }
  // Output::print("NBF b ", NBF);

  // generate interface joints
  for (i=0;i<N_Interfaces;i++)
  {
    auxi = i*6;

    bdcomp = BdParts[InterfaceParam[auxi + 4] - 1]->
                GetBdComp(InterfaceParam[auxi + 5] - 1);
    if(bdcomp->IsFreeBoundary())
    {
      NewJoints.at(InterfaceParam[auxi    ]*6 + InterfaceParam[auxi + 1] - 7) =
      NewJoints.at(InterfaceParam[auxi + 2]*6 + InterfaceParam[auxi + 3] - 7) =
        new TIsoInterfaceJoint3D(bdcomp, Param1, Param2,
                CellTree[InterfaceParam[auxi]-1],
                CellTree[InterfaceParam[auxi+2]-1]);
    }
    else
    {
      NewJoints.at(InterfaceParam[auxi    ]*6 + InterfaceParam[auxi + 1] - 7) =
      NewJoints.at(InterfaceParam[auxi + 2]*6 + InterfaceParam[auxi + 3] - 7) =
        new TInterfaceJoint3D(bdcomp, Param1, Param2,
                CellTree[InterfaceParam[auxi]-1],
                CellTree[InterfaceParam[auxi+2]-1]);
    }
  }

  // new generate EqN-joints
  NewJoints_aux1 = new int[maxElpV];
  NewJoints_aux2 = new int[maxElpV];

  for (i=0;i<N_RootCells;i++)
  {
    CurrCell = CellTree[i];
    N_Faces = CurrCell->GetN_Faces();
    CurrCell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, maxlen);

    auxi = i*6;
    for (j=0;j<N_Faces;j++)
    {
      if (!NewJoints.at(auxi++))
      {
        N_Verts = TmpLen[j];
        aux2 = (KVERT[NVE*i + TmpFV[maxlen*j]] - 1)*maxElpV;
        aux1 = KVEL[aux2++];
        for (k=0;k<aux1;k++)
          NewJoints_aux1[k] = KVEL[aux2++];

        for (k=1;k<N_Verts;k++)
        {
          aux3 = (KVERT[NVE*i + TmpFV[maxlen*j+k]] - 1)*maxElpV;
          aux2 = KVEL[aux3++];
          for (l=0;l<aux2;l++)
            NewJoints_aux2[l] = KVEL[aux3++];

          for (aux3=l=0;l<aux1;l++)
          {
            NewJoints_aux3 = NewJoints_aux1[l];
            for (m=0;m<aux2;m++)
              if (NewJoints_aux3 == NewJoints_aux2[m])
              {
                NewJoints_aux1[aux3++] = NewJoints_aux3;
                break;
              }
          }

          aux1 = aux3;
        }

        if (aux1 == 2)
        {
          if (i == NewJoints_aux1[0])
          {
            aux1 = i;
            aux2 = NewJoints_aux1[1];
          }
          else
          {
            aux1 = i;
            aux2 = NewJoints_aux1[0];
          }
          
          auxj = j * maxlen;
          for (k=0;k<N_Verts;k++)
            Vert[k] = CurrCell->GetVertex(TmpFV[auxj++]);

          CurrCell_aux = CellTree[aux2];

          N_Faces_aux = CurrCell_aux->GetN_Faces();
          CurrCell_aux->GetShapeDesc()->GetFaceVertex(TmpFV_aux,
                          TmpLen_aux, maxlen_aux);

          for (k=0;k<N_Faces_aux;k++)
            if (N_Verts == TmpLen_aux[k])
            {
              auxj = k * maxlen_aux;
              CurrVert = CurrCell_aux->GetVertex(TmpFV_aux[auxj]);
              for (m=l=0;l<N_Verts;l++)
                if (Vert[l] == CurrVert) break;

              if (l != N_Verts)
                for (m=1;m<N_Verts;m++)
                  if (Vert[(N_Verts + l - m) % N_Verts] != CurrCell_aux->
                        GetVertex(TmpFV_aux[auxj + m]))
                    break;

              if (m == N_Verts) break;
            }

          if (k == N_Faces_aux)
          {
            cerr << "Error in ReadGeo: could not find local face" << endl;
            exit (-1);
          }

          NewJoints.at(aux1*6 + j) = NewJoints.at(aux2*6 + k) = new
            TJointEqN(CellTree[aux2], CellTree[aux1]);
          NewJoints.at(aux1*6 + j)->SetMapType(l);
        }
        else
        {
          cerr << "Error in ReadGeo: no element on face" << endl;
          exit (-1);
        }
      } // endif (!NewJoints[auxi++])
    } // endfor j

    // copy joints to cells
    auxi = i*6;
    for (j=0;j<N_Faces;j++)
    {
      CurrCell->SetJoint(j, NewJoints.at(auxi));
      //if(NewJoints[auxi]->GetType() == InterfaceJoint3D)
      //  cout << auxi << " cell: " << i << " joint: " << j << endl;
      auxi++;
    }
  }

  for(i=0;i<N_RootCells;i++)
  {
    CurrCell = CellTree[i];
    N_Faces = CurrCell->GetN_Faces();
    for(j=0;j<N_Faces;j++)
    {
      if(CurrCell->GetJoint(j)->GetType() == InterfaceJoint3D)
       ((TInterfaceJoint3D *)(CurrCell->GetJoint(j)))->CheckOrientation();
    } // endfor j
  } // endfor i

  for(i=0;i<N_RootCells;i++)
  {
      //cout << "Cell number: " << i << endl;
    CurrCell = CellTree[i];
    N_Faces = CurrCell->GetN_Faces();
    CurrCell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, maxlen);

    for(j=0;j<N_Faces;j++)
    {
      CurrJointType = CurrCell->GetJoint(j)->GetType();
      if(CurrJointType == BoundaryFace ||
         CurrJointType == IsoBoundFace ||
         CurrJointType == InterfaceJoint3D ||
         CurrJointType == IsoInterfaceJoint3D)
      {
        if(CurrJointType == BoundaryFace ||
           CurrJointType == IsoBoundFace)
        {
          bdface = (TBoundFace*)(CurrCell->GetJoint(j));
          bdcomp = bdface->GetBoundComp();
        }

        if(CurrJointType == InterfaceJoint3D ||
           CurrJointType == IsoInterfaceJoint3D )
        {
          bdcomp = ((const TInterfaceJoint3D*)(CurrCell->GetJoint(j)))->GetBoundComp();
          Output::print("Interface ",
                        ((TInterfaceJoint3D*)(CurrCell->GetJoint(j)))
                        ->CheckInside(CurrCell));
        }

        //Output::print("Face number: ", j, " vertices:");
        for(k=0;k<TmpLen[j];k++)
        {
          //cout << TmpFV[maxlen*j+k] << ": ";
          CurrCell->GetVertex(TmpFV[maxlen*j+k])->GetCoords(X, Y, Z);
          //cout << "X: " << X << " Y: " << Y << " Z: " << Z;
          bdcomp->GetTSofXYZ(X, Y, Z, T, S);
          //cout << " T: " << T << " S: " << S << endl;
          Param1[k] = T;
          Param2[k] = S;
        } // endfor k


        if(CurrJointType == BoundaryFace ||
           CurrJointType == IsoBoundFace)
        {
          bdface->SetParameters(Param1, Param2);
        }
        else
        {
          ((TInterfaceJoint3D*)(CurrCell->GetJoint(j)))
                ->SetMapType();
          ((TInterfaceJoint3D*)(CurrCell->GetJoint(j)))
                ->SetParameters(Param1, Param2);
        }
      } // endif BoundaryFace
    } // endfor j
  } // endfor i

  for(int i=0; i<N_RootCells; i++)
  {
    CurrCell = CellTree[i];
    CurrCell->SetReference_ID(ELEMSREF[i]);
  }
  
  // free memory
  delete[] KVEL;
  delete[] NewVertices;
  delete[] NewJoints_aux1;
  delete[] NewJoints_aux2;

  return 0;
}

enum class ElementShape {TETS, QUADS};

/**
 * Check in which order three given vertices are, according
 * to the alphanumeric order given by the operator < of TVertex.
 */
int alphanumeric_order(const TVertex& vert0, const TVertex& vert1, const TVertex& vert2)
{
  if( (vert0 < vert1) && (vert1 < vert2) )
    return 0; // 0-1-2
  if( (vert0 < vert2) && (vert2 < vert1) )
    return 1; // 0-2-1
  if( (vert1 < vert0) && (vert0 < vert2) )
    return 2; // 1-0-2
  if( (vert1 < vert2) && (vert2 < vert0) )
    return 3;  // 1-2-0
  if( (vert2 < vert0) && (vert0 < vert1) )
    return 4; // 2-0-1
  if( (vert2 < vert1) && (vert1 < vert0) )
    return 5; // 2-1-0
  else
  {
    ErrThrow("This should not happen!");
    return -1;
  }
}
/*
 * Given a certain 2D base element by 'current_element_id'
 * and one of its edges by the ids of its two vertices,
 * find that other element, if any, which neighbors the current
 * element, sharing the given vertex with it, and return its ID,
 * or -1 if no such element was found, indicating we are dealing
 * with a boundary edge.
 *
 * This assumes, as does the entire process of creating a
 * sandwich grid, NO HANGING NODES.
 *
 * Performance is of lesser concern still.
 */
size_t find_neighbor_base_element_id(
		const std::vector<std::vector<int>>& all_elements_verts,
		size_t first_element_id,
		size_t vertex_1_id, size_t vertex_2_id)
{
	//declare a lambda
	std::function<bool (const std::vector<int>&, size_t, size_t)> find_two_in_vec = [](
			const std::vector<int>& vec,
			int a, int b)
			{
				bool both_contained = std::find(vec.begin(),vec.end(),a) != vec.end()
								   && std::find(vec.begin(),vec.end(),b) != vec.end();
				return both_contained;
			};

	//check if input make sense
	const std::vector<int>& first_element_verts =
			all_elements_verts.at(first_element_id);
	if(!find_two_in_vec(first_element_verts, vertex_1_id, vertex_2_id))
		ErrThrow("Wrong input in find_neighbor_base_element_id, could"
				" not find vertices in given element.")


	size_t elem_id = 0;
	for(auto e : all_elements_verts)
	{

		if(find_two_in_vec(e, vertex_1_id, vertex_2_id) && // both vertices found...
				(elem_id != first_element_id))				   // ...in another element
			return elem_id;

		elem_id++;
	}
	return -1; //the not found case
}

void sort_vertices(TVertex*& A, TVertex*& B, TVertex*& C)
{
	// sort the three vertices consistently (alphanumeric)
	std::function<bool (const TVertex*, const TVertex*)> less_or_eq
			= [](const TVertex* V, const TVertex*W) { return *V < *W; };
	std::vector<TVertex*> sorted = {A,B,C};
	std::sort(sorted.begin(),sorted.end(), less_or_eq);
	A = sorted[0];
	B = sorted[1];
	C = sorted[2];
}

//For given four vertices and an integer index the method sets the tetrahedron at *tet
//with positive orientation (positive triple product)
void set_tetrahedron(int index, TBaseCell *tet, TVertex* vert0, TVertex* vert1,
                     TVertex* vert2, TVertex* vert3, int &orientation)
{
  tet->SetVertex(0, vert0);
  tet->SetVertex(1, vert1);
  tet->SetVertex(2, vert2);
  tet->SetVertex(3, vert3);
  tet->SetClipBoard(index);
  //Change cell orientation (positive triple product, right hand rule)
  orientation = int(tet->check_orientation());
  if(!orientation)
  {
    tet->SetVertex(0, vert1);
    tet->SetVertex(1, vert0);
  }
}

/* A function recieving three cells forming one prism and an edge (two vertices)
 * and returns which (two) cells contain the edge and which of them is the top 
 * cell which is the bottom cell as well as the correct joint IDs.
 * Joint ID has to match with the face ID (ShapeDesc) for later use (SetMaptype)
 */
void find_low_and_high_face_element(bool orientation_low, bool orientation_mid,
                                    TGridCell* const low_element,
                                    TGridCell* const mid_element,
                                    TGridCell* const high_element,
                                    TVertex* const first_bot,
                                    TVertex* const second_bot,
                                    TGridCell*& low_joint_element,
                                    int& low_joint_id,
                                    TGridCell*& high_joint_element,
                                    int& high_joint_id)
{

  TVertex* A_bot = low_element->GetVertex(0); //get verts 0,1,2, by creation they should be the bot verts
  TVertex* B_bot = low_element->GetVertex(1);
  TVertex* C_bot = low_element->GetVertex(2);

  // just for safety, the verts should still be ordered
  sort_vertices(A_bot, B_bot, C_bot);
  // now A < B < C

  if((first_bot == A_bot && second_bot == B_bot) || (first_bot == B_bot
      && second_bot == A_bot))
  { //this is the AB - edge
    low_joint_element = mid_element;
    low_joint_id = 0;
    high_joint_element = high_element;
    high_joint_id = 0;
    return;
  }
  else if((first_bot == B_bot && second_bot == C_bot)
      || (first_bot == C_bot && second_bot == B_bot))
  { //this is the BC - edge
    low_joint_element = low_element;
    if(orientation_low)
      low_joint_id = 2;
    else
      low_joint_id = 3;
    high_joint_element = mid_element;
    if(orientation_mid)
      high_joint_id = 2;
    else
      high_joint_id = 3;
    return;
  }
  else if((first_bot == C_bot && second_bot == A_bot)
      || (first_bot == A_bot && second_bot == C_bot))
  { //this is the AC - edge
    low_joint_element = low_element;
    if(orientation_low)
      low_joint_id = 3;
    else
      low_joint_id = 2;
    high_joint_element = high_element;
    high_joint_id = 1;
    return;
  }
  else
    ErrThrow("Something is terribly wrong here!");
}

//Treat the special case of a closed boundary component (closed circle)
//uses that vert1 and vert2 are part of the bottom plane !!!
void correct_parameterlist(int n_vert, TBoundComp3D * BdComp, TVertex* vert1,
                           TVertex * vert2, double *parameterlist)
{
  if(BdComp->GetType() == Wall)
  {
    TBdWall * BdWall = dynamic_cast<TBdWall *>(BdComp);
    if(BdWall->GetBdComp2D()->GetType() == Circle)
    {
      auto BdCircle = dynamic_cast<TBdCircle *>(BdWall->GetBdComp2D());

      if(BdCircle->is_full_circle())
      {
        //uses that vert1 and vert2 are part of the bottom plane
        double T1, T2;
        BdCircle->GetTofXY(vert1->GetX(), vert1->GetY(), T1);
        BdCircle->GetTofXY(vert2->GetX(), vert2->GetY(), T2);

        if(std::abs(T1 - T2) > 0.5)
        {
          for(int vert_id = 0; vert_id < n_vert; vert_id++)
          {
            if(1 - parameterlist[vert_id] < 1e-4)
            {
              parameterlist[vert_id] = 0;
            }
            else if(parameterlist[vert_id] < 1e-4)
            {
              parameterlist[vert_id] = 1;
            }
          }
        }
      }
    }
  }
}

//Defines boundary joints given the Cell, the face (int joint_id) where the joint is located
//the correct(!) three vertices of this face and the boundary component matching the joint
//If the component is a BdWall, then vert1 and vert2 are part of the bottom plane and
//vert3 and vert4 are part of the top plane (implicitly used)
void set_bdr_joints(TBaseCell *tet_cell, TVertex *vert1, TVertex *vert2,
                    TVertex *vert3, int joint_id, TBoundComp3D * BdComp)
{
  const int *vert_per_face;
  const int *n_vert_per_face;
  int max_n;
  tet_cell->GetShapeDesc()->GetFaceVertex(vert_per_face, n_vert_per_face,
                                          max_n);

  double T_parameter[4];
  double Z_parameter[4];
  double x, y, z;
  for(int vert_id = 0; vert_id < 3; vert_id++)
  {
    int vert = vert_per_face[max_n * joint_id + vert_id];

    if(tet_cell->GetVertex(vert) == vert1)
    {
      vert1->GetCoords(x, y, z);
      BdComp->GetTSofXYZ(x, y, z, T_parameter[vert_id], Z_parameter[vert_id]);
    }
    else if(tet_cell->GetVertex(vert) == vert2)
    {
      vert2->GetCoords(x, y, z);
      BdComp->GetTSofXYZ(x, y, z, T_parameter[vert_id], Z_parameter[vert_id]);
    }
    else if(tet_cell->GetVertex(vert) == vert3)
    {
      vert3->GetCoords(x, y, z);
      BdComp->GetTSofXYZ(x, y, z, T_parameter[vert_id], Z_parameter[vert_id]);
    }
    else
      ErrThrow(
          "Error in set_bdr_joints (Tetra) in BdWall joint \n  joint does not match face: ",
          joint_id);
  }

  //Treat the special case of a closed boundary component (closed circle)
  //uses that vert1 and vert2 are part of the bottom plane
  correct_parameterlist(3, BdComp, vert1, vert2, T_parameter);

  TJoint* bdry_joint = new TBoundFace(BdComp, T_parameter, Z_parameter);
  tet_cell->SetJoint(joint_id, bdry_joint);
  bdry_joint->SetMapType();
}

//Defines boundary joints given the Cell, the face (int joint_id) where the joint is located
//the correct(!) four vertices of this face and the boundary component matching the joint
//If the component is a BdWall, then vert1 and vert2 are part of the bottom plane and
//vert3 and vert4 are part of the top plane (implicitly used)
void set_bdr_joints(TBaseCell *hexa_cell, TVertex *vert1, TVertex *vert2,
                    TVertex *vert3, TVertex *vert4, int joint_id,
                    TBoundComp3D * BdComp)
{
  const int *vert_per_face;
  const int *n_vert_per_face;
  int max_n;
  hexa_cell->GetShapeDesc()->GetFaceVertex(vert_per_face, n_vert_per_face,
                                           max_n);

  double T_parameter[4];
  double Z_parameter[4];
  double x, y, z;
  for(int vert_id = 0; vert_id < 4; vert_id++)
  {
    int vert = vert_per_face[joint_id * max_n + vert_id];
    if(hexa_cell->GetVertex(vert) == vert1)
    {
      vert1->GetCoords(x, y, z);
      BdComp->GetTSofXYZ(x, y, z, T_parameter[vert_id], Z_parameter[vert_id]);
    }
    else if(hexa_cell->GetVertex(vert) == vert2)
    {
      vert2->GetCoords(x, y, z);
      BdComp->GetTSofXYZ(x, y, z, T_parameter[vert_id], Z_parameter[vert_id]);
    }
    else if(hexa_cell->GetVertex(vert) == vert3)
    {
      vert3->GetCoords(x, y, z);
      BdComp->GetTSofXYZ(x, y, z, T_parameter[vert_id], Z_parameter[vert_id]);
    }
    else if(hexa_cell->GetVertex(vert) == vert4)
    {
      vert4->GetCoords(x, y, z);
      BdComp->GetTSofXYZ(x, y, z, T_parameter[vert_id], Z_parameter[vert_id]);
    }
    else
      ErrThrow(
          "Error in set_bdr_joints (Hexa) in BdComp joint \n  joint does not match face: ",
          joint_id)
  }
  //Treat the special case of a closed boundary component (closed circle)

  correct_parameterlist(4, BdComp, vert1, vert2, T_parameter);

  TJoint* outer_joint = new TBoundFace(BdComp, T_parameter, Z_parameter);
  hexa_cell->SetJoint(joint_id, outer_joint);
  outer_joint->SetMapType();
}

//Given a GridCell and an edge (first_bot,second_bot) of bottom face
//this method returns the wall face which contains the edge (both vertices)
int find_joint_face(TGridCell* const hexa_cell, TVertex* const first_bot,
                    TVertex* const second_bot)
{
  const int* vert_per_face;
  const int* n_vert_per_face; //={4,4,4,4,4,4}
  int max_n_vert_per_face;    //=4
  hexa_cell->GetShapeDesc()->GetFaceVertex(vert_per_face, n_vert_per_face,
                                           max_n_vert_per_face);

  //face 0 is the bottom face and face 5 is the top face
  //thus, both faces are not considered
  for(int face_id = 1; face_id < 5; face_id++)
  {
    for(int vert_1_id = 0; vert_1_id < 4; vert_1_id++)
    {
      if(hexa_cell->GetVertex(
          vert_per_face[max_n_vert_per_face * face_id + vert_1_id])
         == first_bot)
      {
        for(int vert_2_id = 0; vert_2_id < 4; vert_2_id++)
        {
          if(hexa_cell->GetVertex(
              vert_per_face[max_n_vert_per_face * face_id + vert_2_id])
             == second_bot)
          {
            return face_id;
          }
        }
      }
    }
  }
  ErrThrow("In find_joint_face the given cell doesn't contain the edge!");
  return -1;
}

/**@brief constructs a 3D grid from a given 2D grid
 * \ref page_sandwich_grid , provides further information on this method
 */
void TDomain::MakeSandwichGrid(double *DCORVG, int *KVERT, int *KNPR,
                                  int N_Vertices, int NVE, double drift_x,
                                  double drift_y, double drift_z,
                                  double conic_scale,
                                  const std::vector<double>& lambda)
{
  size_t n_vert_layers = lambda.size();
  size_t n_elem_layers = lambda.size() - 1;
  size_t n_base_vertices = N_Vertices;
  size_t n_base_elements = this->N_RootCells;

  // TODO NVE is superfluous, given KVERT were a 2D vector
  ElementShape elements_shape;
  if(NVE == 3)
    elements_shape = ElementShape::TETS;
  else if(NVE == 4)
    elements_shape = ElementShape::QUADS;
  else
    ErrThrow(
        "What kind of element should that be? Enter 3 for TETS or 4 for quads.");

  // put input into more suitable forms
  // This contains coordinates of the base (i.e., 2D) vertices.
  // If it is an inner vertex, it contains (x,y) coordinates.
  // If it is a boundary vertex, it contains (t,0) coordinates,
  // with t the parameter on the respective boundary. The containing
  // boundary component is given by the next vector.
  std::vector<std::pair<double, double>> base_vert_coords(n_base_vertices);
  for(size_t v = 0; v < n_base_vertices; ++v)
    base_vert_coords.at(v) =
    { DCORVG[2*v], DCORVG[2*v + 1]};

    // Contains IDs of the boundary (component!) on which each vertex lies.
    // For inner vertices, it is 0.
  std::vector<int> base_vert_bdry(n_base_vertices);
  for(size_t v = 0; v < n_base_vertices; ++v)
    base_vert_bdry.at(v) = KNPR[v] - 1; //add 1, input KNPR has "-1" for inner vertex

  // This data structure gathers which base vertices form which base element.
  std::vector<std::vector<int>> base_verts_per_base_elem(n_base_elements,
                                                         std::vector<int>(NVE));

  for(size_t e = 0; e < n_base_elements; ++e)
  {
    std::vector<int> verts(NVE);
    for(int v = 0; v < NVE; ++v)
      verts.at(v) = KVERT[e * NVE + v] - 1;

    base_verts_per_base_elem.at(e) = verts;
  }
  
  double BdBox_xmin = std::numeric_limits<double>::infinity();
  double BdBox_ymin = std::numeric_limits<double>::infinity();
  double BdBox_xmax = -std::numeric_limits<double>::infinity();
  double BdBox_ymax = -std::numeric_limits<double>::infinity();

  std::vector<TVertex*> sandw_vertices(n_base_vertices * n_vert_layers,
                                       nullptr);
  for(size_t bv = 0; bv < n_base_vertices; ++bv)
  {
    double base_x = 0;
    double base_y = 0;

    if(base_vert_bdry.at(bv) == -1)
    { //this is an inner base vertex - coordinates are cartesian (x,y)
      base_x = base_vert_coords.at(bv).first;
      base_y = base_vert_coords.at(bv).second;
    }
    else
    { //this is a boundary base vertex - coordinates are (tau, 0),
      //tau from a parametrization of the boundary given by base_vert_bdry
      double comp_plus_tau = base_vert_coords.at(bv).first; //second is irrelevant
      int part = base_vert_bdry.at(bv);
      int comp = std::floor(comp_plus_tau); //floor tau to the nearest integer
      double tau = comp_plus_tau - comp;

      //The boundary on which the vertex lies is encoded in the input as follows:
      // BOUNDARY PART		: base_vert_bdry.at(v)
      // BOUNDARY COMPONENT : std::floor(comp_plus_tau)
      // TAU (ON PART, COMP): comp_plus_tau - comp

      if(GetLastLocalComp(part) == comp - 1)
      {      //the very last vertex - correct "comp" not reached by std::floor
        comp--;
        tau = 1.0;
      }

      TBdWall* BdWall = dynamic_cast<TBdWall*>(BdParts[part]->GetBdComp(comp));
      BdWall->GetBdComp2D()->GetXYofT(tau, base_x, base_y);
    }

    //compute parameters for the bounding box
    if(base_x < BdBox_xmin)
      BdBox_xmin = base_x;
    if(base_y < BdBox_ymin)
      BdBox_ymin = base_y;
    if(base_x > BdBox_xmax)
      BdBox_xmax = base_x;
    if(base_y > BdBox_ymax)
      BdBox_ymax = base_y;
    
    // ...in both cases, create all corresponding 3D vertices
    for(size_t l = 0; l < n_vert_layers; ++l)
    {
      double sandw_x = base_x + drift_x * lambda.at(l) + base_x * conic_scale * lambda.at(l);
      double sandw_y = base_y + drift_y * lambda.at(l) + base_y * conic_scale * lambda.at(l);
      double sandw_z = drift_z * lambda.at(l);
      // create a new 3D vertex and put it into the list
      int new_vertex_id = l * n_base_vertices + bv;
      sandw_vertices.at(new_vertex_id) = new TVertex(sandw_x, sandw_y, sandw_z);
      if(base_vert_bdry.at(bv) != -1 || l == 0 || l == n_vert_layers-1)
        sandw_vertices.at(new_vertex_id)->SetAsBoundVert();
    }
  }      //end 2D vertices

         //Set the bounding box.
  StartX = BdBox_xmin;
  StartY = BdBox_ymin;
  StartZ = 0;
  BoundX = (BdBox_xmax - BdBox_xmin)*std::max(conic_scale,1.) + drift_x; //coarse bounds
  BoundY = (BdBox_ymax - BdBox_ymin)*std::max(conic_scale,1.) + drift_y; //to be improved
  BoundZ = drift_z;

  //Set BdWalls
  int N_BdWalls = 0;
  for(int part = 0; part < N_BoundParts; part++)
  {
    int N_BoundComp = BdParts[part]->GetN_BdComps();
    N_BdWalls += N_BoundComp;

    for(int comp = 0; comp < N_BoundComp; comp++)
    {
      TBdWall* BdWall = dynamic_cast<TBdWall*>(BdParts[part]->GetBdComp(comp));
      BdWall->SetParams(drift_x, drift_y, drift_z,conic_scale);
    }
  }

  // The BdComp.IDs (N_BdWalls and N_BdWalls+1) may not be important,
  // probably they aren't used later anyway...
  TBdPlane* Bottom = new TBdPlane(N_BdWalls);
  Bottom->SetParams(0, 0, 0, 1, 0, 0, 0, 0, -1);
  TBdPlane* Top = new TBdPlane(N_BdWalls + 1);
  Top->SetParams(0, 0, drift_z, 1, 0, 0, 0, 0, 1);

  if(elements_shape == ElementShape::QUADS)
  {
    size_t n_sandw_elements = n_base_elements * n_elem_layers;
    CellTree = new TBaseCell*[n_sandw_elements];

    //Declare Cells (Hexahedrons)
    //for (size_t c=0; c<n_sandw_elements; c++)
    //  CellTree[c] = new TMacroCell(TDatabase::RefDescDB[Hexahedron], RefLevel);

    for(size_t e = 0; e < n_base_elements; e++)
    {
      int base_vert_1_id = base_verts_per_base_elem.at(e)[0];
      int base_vert_2_id = base_verts_per_base_elem.at(e)[1];
      int base_vert_3_id = base_verts_per_base_elem.at(e)[2];
      int base_vert_4_id = base_verts_per_base_elem.at(e)[3];

      //Set hexahedrons and top-bottom-joints
      for(size_t l = 0; l < n_elem_layers; l++)
      {
        TVertex* A_bot = sandw_vertices.at(
            l * n_base_vertices + base_vert_1_id);
        TVertex* B_bot = sandw_vertices.at(
            l * n_base_vertices + base_vert_2_id);
        TVertex* C_bot = sandw_vertices.at(
            l * n_base_vertices + base_vert_3_id);
        TVertex* D_bot = sandw_vertices.at(
            l * n_base_vertices + base_vert_4_id);

        TVertex* A_top = sandw_vertices.at(
            (l + 1) * n_base_vertices + base_vert_1_id);
        TVertex* B_top = sandw_vertices.at(
            (l + 1) * n_base_vertices + base_vert_2_id);
        TVertex* C_top = sandw_vertices.at(
            (l + 1) * n_base_vertices + base_vert_3_id);
        TVertex* D_top = sandw_vertices.at(
            (l + 1) * n_base_vertices + base_vert_4_id);

        TVertex * locVertices[8] = {A_bot, B_bot, C_bot, D_bot, A_top, B_top,
                                    C_top, D_top};

        auto CellType = ((THexahedron *) TDatabase::RefDescDB[Hexahedron]
            ->GetShapeDesc())->CheckHexa(locVertices);

        //Set Cell
        CellTree[l * n_base_elements + e] = new TMacroCell(
            TDatabase::RefDescDB[CellType], RefLevel);
        TBaseCell* hexa = CellTree[l * n_base_elements + e];
        hexa->SetVertex(0, A_bot);
        hexa->SetVertex(1, B_bot);
        hexa->SetVertex(2, C_bot);
        hexa->SetVertex(3, D_bot);
        hexa->SetVertex(4, A_top);
        hexa->SetVertex(5, B_top);
        hexa->SetVertex(6, C_top);
        hexa->SetVertex(7, D_top);

        //Correct orientation of the grid cell
        //Should never enter this case...
        if(!hexa->check_orientation())
        {
          Output::print(
              "Something is wrong with the order of vertices in the given quadrilateral");
          hexa->SetVertex(0, B_bot);
          hexa->SetVertex(1, A_bot);
          hexa->SetVertex(2, D_bot);
          hexa->SetVertex(3, C_bot);
          hexa->SetVertex(4, B_top);
          hexa->SetVertex(5, A_top);
          hexa->SetVertex(6, D_top);
          hexa->SetVertex(7, C_top);
        }

        if(l == 0)
          set_bdr_joints(hexa, A_bot, B_bot, C_bot, D_bot, 0, Bottom);

        if(l > 0)
        {
          TBaseCell *hexa_below = CellTree[(l - 1) * n_base_elements + e];
          TJointEqN* top_bot_outer_joint = new TJointEqN(hexa, hexa_below);
          hexa_below->SetJoint(5, top_bot_outer_joint);
          hexa->SetJoint(0, top_bot_outer_joint);
          top_bot_outer_joint->SetMapType();
        }

        if(l == n_elem_layers - 1)
          set_bdr_joints(hexa, A_top, B_top, C_top, D_top, 5, Top);

      } //end layers

    } //end base_elements

    for(size_t e = 0; e < n_base_elements; e++)
    {
      int base_vert_0_id = base_verts_per_base_elem.at(e)[0];
      int base_vert_1_id = base_verts_per_base_elem.at(e)[1];
      int base_vert_2_id = base_verts_per_base_elem.at(e)[2];
      int base_vert_3_id = base_verts_per_base_elem.at(e)[3];

      //treat all three edges, inserting all side joints in all sandwich elements
      //first edge
      std::vector<std::pair<int, int>> edges = {
          {base_vert_0_id, base_vert_1_id}, {base_vert_1_id, base_vert_2_id}, {
              base_vert_2_id, base_vert_3_id},
          {base_vert_3_id, base_vert_0_id}};

      //keep track of which elements were treated already
      std::vector<std::vector<int>> treated_edges(n_base_vertices,
                                                  std::vector<int>());

      for(int edge = 0; edge < 4; edge++)
      {
        int vert_1_id = edges.at(edge).first;
        int vert_2_id = edges.at(edge).second;

        //Check on the element containing the edge
        if(std::find(treated_edges.at(vert_1_id).begin(),
                     treated_edges.at(vert_1_id).end(), vert_2_id)
           != treated_edges.at(vert_1_id).end())
          continue;

        treated_edges.at(vert_1_id).push_back(vert_2_id);
        treated_edges.at(vert_2_id).push_back(vert_1_id);

        int e_neighbor = find_neighbor_base_element_id(base_verts_per_base_elem,
                                                       e, vert_1_id, vert_2_id);

        //treat boundary edges
        if(e_neighbor == -1)
        {
          // check conditions
          if(Interfaces[KNPR[vert_1_id] - 1] < 0 || Interfaces[KNPR[vert_2_id]
              - 1]
                                                    < 0)
          {
            ErrThrow("Error in ReadGeo");
          }

          //treatment of boundary edges
          int edge_bdr_part = base_vert_bdry.at(vert_1_id);
          int edge_bdr_comp = static_cast<int>(base_vert_coords.at(vert_1_id)
              .first);

          if(BdParts[edge_bdr_part]->GetBdComp(edge_bdr_comp)->IsFreeBoundary())
          {
            ErrThrow("Error in ReadGeo");
          }

          for(size_t l = 0; l < n_elem_layers; l++)
          {
            TGridCell* my_hexa = dynamic_cast<TGridCell*>(CellTree[l
                * n_base_elements
                                                                   + e]);

            TVertex* vert1(sandw_vertices.at(l * n_base_vertices + vert_1_id));
            TVertex* vert2(sandw_vertices.at(l * n_base_vertices + vert_2_id));
            TVertex* vert3(
                sandw_vertices.at((l + 1) * n_base_vertices + vert_1_id));
            TVertex* vert4(
                sandw_vertices.at((l + 1) * n_base_vertices + vert_2_id));
            double T_vert_1 = base_vert_coords.at(vert_1_id).first
                - edge_bdr_comp;
            double T_vert_2 = base_vert_coords.at(vert_2_id).first
                - edge_bdr_comp;

            //prevents confusion in vertices of two different boundary-parts
            //guarantees 0<=T<=1
            if(T_vert_2 < T_vert_1)
              T_vert_2 = 1.;

            int my_joint_id = find_joint_face(my_hexa, vert1, vert2);

            set_bdr_joints(my_hexa, vert1, vert2, vert3, vert4, my_joint_id,
                           BdParts[edge_bdr_part]->GetBdComp(edge_bdr_comp));
          }

        }	        		  //end boundary treatment
        else
        {
          for(size_t l = 0; l < n_elem_layers; l++)
          {
            TGridCell* my_hexa = dynamic_cast<TGridCell*>(CellTree[l
                * n_base_elements
                                                                   + e]);
            TGridCell* yo_hexa =
                dynamic_cast<TGridCell*>(CellTree[l * n_base_elements
                    + e_neighbor]);
            TVertex* vert1(sandw_vertices.at(l * n_base_vertices + vert_1_id));
            TVertex* vert2(sandw_vertices.at(l * n_base_vertices + vert_2_id));

            int my_joint_id = find_joint_face(my_hexa, vert1, vert2);
            int yo_joint_id = find_joint_face(yo_hexa, vert1, vert2);

            TJointEqN* inner_joint = new TJointEqN(my_hexa, yo_hexa);
            my_hexa->SetJoint(my_joint_id, inner_joint);
            yo_hexa->SetJoint(yo_joint_id, inner_joint);
            inner_joint->SetMapType();
          }
        }

      }
      // N_RootCells auf die Zahl 3D (nicht: 2D) Zellen setzen
      this->N_RootCells = n_sandw_elements;

    }

    //ErrThrow("NewMakeSandwichGrid only deals with triangular/"
    //"tetrehedral meshes so far. ");
  }

  else if(elements_shape == ElementShape::TETS)
  {
    // Allocate the CellTree, i.e., the tree structure of 3D cells.
    // Each 2D base element will lead to 3 times n_elem_layers cells
    // in the CellTree (each prism with triangular base area will be split into three tets)
    size_t n_sandw_elements = 3 * n_base_elements * n_elem_layers;
    CellTree = new TBaseCell*[n_sandw_elements];
    for(size_t c = 0; c < n_sandw_elements; ++c)
      CellTree[c] = new TMacroCell(TDatabase::RefDescDB[Tetrahedron], RefLevel);

    //For further use the orientation (sign of the triple product of the vertices)
    //has to be positive. If not it has to be changed and the numbering of faces/joints has to be adjusted.
    std::vector<int> Cell_orientations(n_sandw_elements);

    for(size_t e = 0; e < n_base_elements; e++)
    {

      int base_vert_1_id = base_verts_per_base_elem.at(e)[0];
      int base_vert_2_id = base_verts_per_base_elem.at(e)[1];
      int base_vert_3_id = base_verts_per_base_elem.at(e)[2];

      for(size_t l = 0; l < n_elem_layers; l++)
      {
        TVertex* A_bot = sandw_vertices.at(
            l * n_base_vertices + base_vert_1_id);
        TVertex* B_bot = sandw_vertices.at(
            l * n_base_vertices + base_vert_2_id);
        TVertex* C_bot = sandw_vertices.at(
            l * n_base_vertices + base_vert_3_id);
        sort_vertices(A_bot, B_bot, C_bot);
        TVertex* A_top = sandw_vertices.at(
            (l + 1) * n_base_vertices + base_vert_1_id);
        TVertex* B_top = sandw_vertices.at(
            (l + 1) * n_base_vertices + base_vert_2_id);
        TVertex* C_top = sandw_vertices.at(
            (l + 1) * n_base_vertices + base_vert_3_id);
        sort_vertices(A_top, B_top, C_top);

        size_t my_prism_index = 3 * (l * n_base_elements + e);

        // Set up "low" Tetrahedron
        TBaseCell* low_tet = CellTree[my_prism_index];
        set_tetrahedron(my_prism_index, low_tet, A_bot, B_bot, C_bot, C_top,
                        Cell_orientations.at(my_prism_index));

        // Set up "mid" Tetrahedron
        TBaseCell* mid_tet = CellTree[my_prism_index + 1];
        set_tetrahedron(my_prism_index + 1, mid_tet, A_bot, B_bot, B_top, C_top,
                        Cell_orientations.at(my_prism_index + 1));

        // Set up "high" Tetrahedron
        TBaseCell* high_tet = CellTree[my_prism_index + 2];
        set_tetrahedron(my_prism_index + 2, high_tet, A_bot, A_top, B_top,
                        C_top, Cell_orientations.at(my_prism_index + 2));

        // put both inner joints into place
        // Joint ID has to match with the face ID (ShapeDesc) for later use (SetMaptype)
        TJointEqN* low_mid_inner_joint = new TJointEqN(low_tet, mid_tet);
        low_tet->SetJoint(1, low_mid_inner_joint);
        mid_tet->SetJoint(1, low_mid_inner_joint);

        low_mid_inner_joint->SetMapType();

        TJointEqN* mid_high_inner_joint = new TJointEqN(mid_tet, high_tet);
        //Joint face depends on orientation
        if(Cell_orientations[my_prism_index + 1])
          mid_tet->SetJoint(3, mid_high_inner_joint);
        else
          mid_tet->SetJoint(2, mid_high_inner_joint);

        //Joint face depends on orientation
        if(Cell_orientations[my_prism_index + 2])
          high_tet->SetJoint(3, mid_high_inner_joint);
        else
          high_tet->SetJoint(2, mid_high_inner_joint);

        mid_high_inner_joint->SetMapType();

        //set bottom plane joint
        if(l == 0)
          set_bdr_joints(low_tet, A_bot, B_bot, C_bot, 0, Bottom);

        // put outer joint between layers into place
        if(l > 0)
        {
          TBaseCell *top_tet = CellTree[3 * ((l - 1) * n_base_elements + e) + 2];
          TJointEqN* top_bot_outer_joint = new TJointEqN(low_tet, top_tet);
          if(Cell_orientations[3 * ((l - 1) * n_base_elements + e) + 2])
            top_tet->SetJoint(2, top_bot_outer_joint);
          else
            top_tet->SetJoint(3, top_bot_outer_joint);
          low_tet->SetJoint(0, top_bot_outer_joint);
          top_bot_outer_joint->SetMapType();
        }

        //set top plane joint
        if(l == n_elem_layers - 1)
        {
          if(Cell_orientations[my_prism_index + 2])
            set_bdr_joints(high_tet, A_top, B_top, C_top, 2, Top);
          else
            set_bdr_joints(high_tet, A_top, B_top, C_top, 3, Top);
        }
      }      //end layers
    }      //end prisms

           //to keep track of which elements were treated already
    std::vector<std::vector<int>> treated_edges(n_base_vertices,
                                                std::vector<int>());

    // set boundary joints
    for(size_t e = 0; e < n_base_elements; e++)
    {
      int base_vert_0_id = base_verts_per_base_elem.at(e)[0];
      int base_vert_1_id = base_verts_per_base_elem.at(e)[1];
      int base_vert_2_id = base_verts_per_base_elem.at(e)[2];
      //treat all three edges, inserting all side joints in all sandwich elements
      //first edge
      std::vector<std::pair<int, int>> edges = {
          {base_vert_0_id, base_vert_1_id}, {base_vert_1_id, base_vert_2_id}, {
              base_vert_2_id, base_vert_0_id}};

      for(int edge = 0; edge < 3; ++edge)
      {

        int vert_1_id = edges.at(edge).first;
        int vert_2_id = edges.at(edge).second;

        //Check on the element containing the edge
        if(std::find(treated_edges.at(vert_1_id).begin(),
                     treated_edges.at(vert_1_id).end(), vert_2_id)
           != treated_edges.at(vert_1_id).end())
          continue;

        treated_edges.at(vert_1_id).push_back(vert_2_id);
        treated_edges.at(vert_2_id).push_back(vert_1_id);

        int e_neighbor = find_neighbor_base_element_id(base_verts_per_base_elem,
                                                       e, vert_1_id, vert_2_id);

        if(e_neighbor == -1)
        {    //this is a boundary edge, treat it accordingly
             // check conditions
          if(Interfaces[KNPR[vert_1_id] - 1] < 0
              || Interfaces[KNPR[vert_2_id]- 1] < 0)
          {
            ErrThrow("Error in ReadGeo" );
          }

          //treatment of boundary edges
          int edge_bdr_part = base_vert_bdry.at(vert_1_id);
          int edge_bdr_comp = static_cast<int>(base_vert_coords.at(vert_1_id)
              .first);

          if(BdParts[edge_bdr_part]->GetBdComp(edge_bdr_comp)->IsFreeBoundary())
          {
            ErrThrow("Error in ReadGeo " );
          }

          double T_vert_1 = base_vert_coords.at(vert_1_id).first
              - edge_bdr_comp;
          double T_vert_2 = base_vert_coords.at(vert_2_id).first
              - edge_bdr_comp;

          //prevents confusion in vertices of two different boundary-parts
          //guarantees 0<=T<=1
          //runs without, too, but causes error in BdLine
          if(T_vert_2 < T_vert_1)
            T_vert_2 = 1.;

          for(size_t l = 0; l < n_elem_layers; ++l)
          {
            size_t my_prism_index = 3 * (l * n_base_elements + e);

            TGridCell* my_low_tet =
                dynamic_cast<TGridCell*>(CellTree[my_prism_index]);
            TGridCell* my_mid_tet =
                dynamic_cast<TGridCell*>(CellTree[my_prism_index + 1]);
            TGridCell* my_high_tet =
                dynamic_cast<TGridCell*>(CellTree[my_prism_index + 2]);

            bool orientation_low = Cell_orientations.at(my_prism_index);
            bool orientation_mid = Cell_orientations.at(my_prism_index + 1);

            TVertex* vert1(sandw_vertices.at(l * n_base_vertices + vert_1_id));
            TVertex* vert2(sandw_vertices.at(l * n_base_vertices + vert_2_id));
            TVertex* vert1_high(
                sandw_vertices.at((l + 1) * n_base_vertices + vert_1_id));
            TVertex* vert2_high(
                sandw_vertices.at((l + 1) * n_base_vertices + vert_2_id));

            int my_lower_joint_id, my_upper_joint_id;
            TGridCell* my_lower_joint_element;
            TGridCell* my_upper_joint_element;

            find_low_and_high_face_element(orientation_low, orientation_mid,
                                           my_low_tet, my_mid_tet, my_high_tet,
                                           vert1, vert2, my_lower_joint_element,
                                           my_lower_joint_id,
                                           my_upper_joint_element,
                                           my_upper_joint_id);

            //there exists an edge between (smaller vertex, low) and (larger vertex, up)
            if(*vert2 < *vert1)
            {
              //Set joint for lower tet
              set_bdr_joints(my_lower_joint_element, vert1, vert2, vert1_high,
                             my_lower_joint_id,
                             BdParts[edge_bdr_part]->GetBdComp(edge_bdr_comp));
              //Set joint for upper tet
              set_bdr_joints(my_upper_joint_element, vert1_high, vert2_high,
                             vert2, my_upper_joint_id,
                             BdParts[edge_bdr_part]->GetBdComp(edge_bdr_comp));
            }
            else
            {
              //Set joint for lower tet
              set_bdr_joints(my_lower_joint_element, vert1, vert2, vert2_high,
                             my_lower_joint_id,
                             BdParts[edge_bdr_part]->GetBdComp(edge_bdr_comp));
              //Set joint for upper tet
              set_bdr_joints(my_upper_joint_element, vert1_high, vert2_high,
                             vert1, my_upper_joint_id,
                             BdParts[edge_bdr_part]->GetBdComp(edge_bdr_comp));
            }
          }    		  //end layers
        }    		  //end set boundary joints

        //set inner joints between prisms
        else
        {	  // there are two tets sharing a face to both sides of all sandwich
          // repetitions of the current edge - create their joints!
          for(size_t l = 0; l < n_elem_layers; ++l)
          {
            size_t my_prism_index = 3 * (l * n_base_elements + e);
            size_t yo_prism_index = 3 * (l * n_base_elements + e_neighbor);

            TGridCell* my_low_tet =
                dynamic_cast<TGridCell*>(CellTree[my_prism_index]);
            TGridCell* my_mid_tet =
                dynamic_cast<TGridCell*>(CellTree[my_prism_index + 1]);
            TGridCell* my_high_tet =
                dynamic_cast<TGridCell*>(CellTree[my_prism_index + 2]);

            bool my_orientation_low = Cell_orientations.at(my_prism_index);
            bool my_orientation_mid = Cell_orientations.at(my_prism_index + 1);

            TGridCell* yo_low_tet =
                dynamic_cast<TGridCell*>(CellTree[yo_prism_index]);
            TGridCell* yo_mid_tet =
                dynamic_cast<TGridCell*>(CellTree[yo_prism_index + 1]);
            TGridCell* yo_high_tet =
                dynamic_cast<TGridCell*>(CellTree[yo_prism_index + 2]);

            bool yo_orientation_low = Cell_orientations.at(yo_prism_index);
            bool yo_orientation_mid = Cell_orientations.at(yo_prism_index + 1);

            TVertex* vert1(sandw_vertices.at(l * n_base_vertices + vert_1_id)); //global numbering
            TVertex* vert2(sandw_vertices.at(l * n_base_vertices + vert_2_id));

            { //set the lower of the two joints
              int my_lower_joint_id, yo_lower_joint_id;
              TGridCell* my_lower_joint_element;
              TGridCell* yo_lower_joint_element;

              int my_upper_joint_id, yo_upper_joint_id;
              TGridCell* my_upper_joint_element;
              TGridCell* yo_upper_joint_element;

              find_low_and_high_face_element(my_orientation_low,
                                             my_orientation_mid, my_low_tet,
                                             my_mid_tet, my_high_tet, vert1,
                                             vert2, my_lower_joint_element,
                                             my_lower_joint_id,
                                             my_upper_joint_element,
                                             my_upper_joint_id);

              find_low_and_high_face_element(yo_orientation_low,
                                             yo_orientation_mid, yo_low_tet,
                                             yo_mid_tet, yo_high_tet, vert1,
                                             vert2, yo_lower_joint_element,
                                             yo_lower_joint_id,
                                             yo_upper_joint_element,
                                             yo_upper_joint_id);

              TJointEqN* low_joint = new TJointEqN(my_lower_joint_element,
                                                   yo_lower_joint_element);
              my_lower_joint_element->SetJoint(my_lower_joint_id, low_joint);
              yo_lower_joint_element->SetJoint(yo_lower_joint_id, low_joint);
              low_joint->SetMapType();

              TJointEqN* high_joint = new TJointEqN(my_upper_joint_element,
                                                    yo_upper_joint_element);
              my_upper_joint_element->SetJoint(my_upper_joint_id, high_joint);
              yo_upper_joint_element->SetJoint(yo_upper_joint_id, high_joint);
              high_joint->SetMapType();
            }
          }
        } //end Set joints between prisms
      }

    }
    // N_RootCells auf die Zahl 3D (nicht: 2D) Zellen setzen
    this->N_RootCells = n_sandw_elements;
  } //END construct tetrahedral mesh

  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);

}

#endif // __2D__

