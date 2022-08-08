
// @(#)Tests.C        1.3 07/20/99

#include <Database.h>
#include <Domain.h>
#include <JointEqN.h>
#include <PeriodicJoint.h>
#include <MacroCell.h>
#include <BdLine.h>
#ifndef __3D__
  #include <BoundEdge.h>
#else
  #include <BoundFace.h>
#endif


#ifndef __3D__

void TDomain::SetBoundBox(double boundx, double boundy)
{
  BoundX = boundx;
  BoundY = boundy;
}

void TDomain::SetBoundBoxstart(double startx , double starty)
{
  StartX = startx;
  StartY = starty;
}

void TDomain::TestGrid1()
{
  TVertex *v[5];
  TBoundEdge *b[5];
  TJointEqN *j;

  v[0] = new TVertex(0.0, 0.0);
  v[1] = new TVertex(1.0, 0.0);
  v[2] = new TVertex(1.5, 0.5);
  v[3] = new TVertex(1.0, 1.0);
  v[4] = new TVertex(0.0, 1.0);

  b[0] = new TBoundEdge(BdParts[0]->GetBdComp(0), 0.0, 1.0);
  b[1] = new TBoundEdge(BdParts[0]->GetBdComp(1), 0.0, 1.0);
  b[2] = new TBoundEdge(BdParts[0]->GetBdComp(2), 0.0, 1.0);
  b[3] = new TBoundEdge(BdParts[0]->GetBdComp(3), 0.0, 1.0);
  b[4] = new TBoundEdge(BdParts[0]->GetBdComp(4), 0.0, 1.0);

  CellTree = new TBaseCell*[2];
  N_InitVCells = 2;
  N_RootCells = 2;

  SetBoundBox(1.5, 1.);

  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);

  CellTree[0] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadToTri0], RefLevel);
  CellTree[1] = new TGridCell(TDatabase::RefDescDB[Triangle],
                      RefLevel);

  CellTree[0]->SetClipBoard(Refinement);
  CellTree[1]->SetClipBoard(Refinement);

  j = new TJointEqN(CellTree[0], CellTree[1]);

  CellTree[0]->SetVertex(0, v[0]);
  CellTree[0]->SetVertex(1, v[1]);
  CellTree[0]->SetVertex(2, v[3]);
  CellTree[0]->SetVertex(3, v[4]);
  CellTree[1]->SetVertex(0, v[1]);
  CellTree[1]->SetVertex(1, v[2]);
  CellTree[1]->SetVertex(2, v[3]);

  CellTree[0]->SetJoint(0, b[0]);
  CellTree[0]->SetJoint(1, j);
  CellTree[0]->SetJoint(2, b[3]);
  CellTree[0]->SetJoint(3, b[4]);
  CellTree[1]->SetJoint(0, b[1]);
  CellTree[1]->SetJoint(1, b[2]);
  CellTree[1]->SetJoint(2, j);
}

/******************************************************************************/
/*                                                                            */
/* creates an inital triangulation of the unit square consisting of the two   */
/* triangles with the vertices (0,0), (1,1), (0,1) and (0,0), (1,0), (1,1)    */
/*                                                                            */
/******************************************************************************/

void TDomain::TwoTriangles(bool bottomleft_to_topright)
{
  TVertex *v[4];
  TBoundEdge *b[4];
  TJointEqN *j;

  // vertices of the unit square
  v[0] = new TVertex(0.0, 0.0);
  v[1] = new TVertex(1.0, 0.0);
  v[2] = new TVertex(1.0, 1.0);
  v[3] = new TVertex(0.0, 1.0);

  // edges of the triangles which are on the boundary
  b[0] = new TBoundEdge(BdParts[0]->GetBdComp(0), 0.0, 1.0);
  b[1] = new TBoundEdge(BdParts[0]->GetBdComp(1), 0.0, 1.0);
  b[2] = new TBoundEdge(BdParts[0]->GetBdComp(2), 0.0, 1.0);
  b[3] = new TBoundEdge(BdParts[0]->GetBdComp(3), 0.0, 1.0);

  // allocate cell tree with the two triangles
  CellTree = new TBaseCell*[2];
  N_InitVCells = 2;
  N_RootCells = 2;
  
  // set bounding box for output
  SetBoundBoxstart(0., 0.);
  SetBoundBox(1, 1);

  // set the iterators, the pointer this stands for the 
  // TDomain which calls this routine
  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);

  // include two triangles into the cell tree
  CellTree[0] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[1] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);

  CellTree[0]->SetClipBoard(Refinement);
  CellTree[1]->SetClipBoard(Refinement);

  // allocate joint (edge) between the triangles
  j = new TJointEqN(CellTree[0], CellTree[1]);

  // set the vertices of triangle 0 at (0,0), (1,0), (1,1)
  CellTree[0]->SetVertex(0, v[0]);
  CellTree[0]->SetVertex(1, v[1]);
  CellTree[0]->SetVertex(2, v[bottomleft_to_topright ? 2 : 3]);

  // set the vertices of triangle 1 at (0,0), (1,1), (0,1)
  CellTree[1]->SetVertex(0, v[bottomleft_to_topright ? 0 : 1]);
  CellTree[1]->SetVertex(1, v[2]);
  CellTree[1]->SetVertex(2, v[3]);

  // set the edges of triangle 0 (bottom boundary, right boundary, inner)
  if(bottomleft_to_topright)
  {
    CellTree[0]->SetJoint(0, b[0]);
    CellTree[0]->SetJoint(1, b[1]);
    CellTree[0]->SetJoint(2, j);
  }
  else
  {
    CellTree[0]->SetJoint(0, b[0]);
    CellTree[0]->SetJoint(1, j);
    CellTree[0]->SetJoint(2, b[3]);
  }

  // set the edges of triangle 0 (inner, upper boundary, left boundary)
  if(bottomleft_to_topright)
  {
    CellTree[1]->SetJoint(0, j);
    CellTree[1]->SetJoint(1, b[2]);
    CellTree[1]->SetJoint(2, b[3]);
  }
  else
  {
    CellTree[1]->SetJoint(0, b[1]);
    CellTree[1]->SetJoint(1, b[2]);
    CellTree[1]->SetJoint(2, j);
  }
}

/******************************************************************************/
/*                                                                            */
/* creates an inital triangulation of the unit square consisting of the eight */
/* triangles                                                                  */
/* initially a triangulation with the vertices (0,0), (1,1), (0,1)            */
/* and (0,0), (1,0), (1,1) is generated which is then refined once            */
/*                                                                            */
/******************************************************************************/

void TDomain::TwoTrianglesRef()
{
  TVertex *v[4];
  TBoundEdge *b[4];
  TJointEqN *j;

  // vertices of the unit square
  v[0] = new TVertex(0.0, 0.0);
  v[1] = new TVertex(1.0, 0.0);
  v[2] = new TVertex(1.0, 1.0);
  v[3] = new TVertex(0.0, 1.0);

  // edges of the triangles which are on the boundary
  b[0] = new TBoundEdge(BdParts[0]->GetBdComp(0), 0.0, 1.0);
  b[1] = new TBoundEdge(BdParts[0]->GetBdComp(1), 0.0, 1.0);
  b[2] = new TBoundEdge(BdParts[0]->GetBdComp(2), 0.0, 1.0);
  b[3] = new TBoundEdge(BdParts[0]->GetBdComp(3), 0.0, 1.0);

  // allocate cell tree with the two triangles
  CellTree = new TBaseCell*[2];
  N_InitVCells = 2;
  N_RootCells = 2;

  // set bounding box for output
  SetBoundBoxstart(0., 0.);
  SetBoundBox(1, 1);

  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);

  // include two triangles into the cell tree
  CellTree[0] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[1] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);

  CellTree[0]->SetClipBoard(Refinement);
  CellTree[1]->SetClipBoard(Refinement);

  // allocate joint (edge) between the triangles
  j = new TJointEqN(CellTree[0], CellTree[1]);

  // set the vertices of triangle 0 at (0,0), (1,0), (1,1)
  CellTree[0]->SetVertex(0, v[0]);
  CellTree[0]->SetVertex(1, v[1]);
  CellTree[0]->SetVertex(2, v[2]);

  // set the vertices of triangle 1 at (0,0), (1,1), (0,1)
  CellTree[1]->SetVertex(0, v[0]);
  CellTree[1]->SetVertex(1, v[2]);
  CellTree[1]->SetVertex(2, v[3]);

  // set the edges of triangle 0 (bottom boundary, right boundary, inner)
  CellTree[0]->SetJoint(0, b[0]);
  CellTree[0]->SetJoint(1, b[1]);
  CellTree[0]->SetJoint(2, j);

   // set the edges of triangle 0 (inner, upper boundary, left boundary)
 CellTree[1]->SetJoint(0, j);
  CellTree[1]->SetJoint(1, b[2]);
  CellTree[1]->SetJoint(2, b[3]);

  // refine once
  CellTree[0]->Refine(RefLevel);
}

void TDomain::SquareInSquare()
{
  TVertex *v[8];
  TBoundEdge *b[8];
  TJointEqN *j[4];

  v[0] = new TVertex(0.0, 0.0);
  v[1] = new TVertex(0.5, 0.0);
  v[2] = new TVertex(1.0, 0.0);
  v[3] = new TVertex(0.0, 0.5);
  v[4] = new TVertex(1.0, 0.5);
  v[5] = new TVertex(0.0, 1.0);
  v[6] = new TVertex(0.5, 1.0);
  v[7] = new TVertex(1.0, 1.0);

  b[0] = new TBoundEdge(BdParts[0]->GetBdComp(0), 0.0, 0.5);
  b[1] = new TBoundEdge(BdParts[0]->GetBdComp(0), 0.5, 1.0);
  b[2] = new TBoundEdge(BdParts[0]->GetBdComp(1), 0.0, 0.5);
  b[3] = new TBoundEdge(BdParts[0]->GetBdComp(1), 0.5, 1.0);
  b[4] = new TBoundEdge(BdParts[0]->GetBdComp(2), 0.0, 0.5);
  b[5] = new TBoundEdge(BdParts[0]->GetBdComp(2), 0.5, 1.0);
  b[6] = new TBoundEdge(BdParts[0]->GetBdComp(3), 0.0, 0.5);
  b[7] = new TBoundEdge(BdParts[0]->GetBdComp(3), 0.5, 1.0);

  CellTree = new TBaseCell*[5];
  N_InitVCells = 5;
  N_RootCells = 5;

  // set bounding box for output
  SetBoundBoxstart(0., 0.);
  SetBoundBox(1, 1);

  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);

  CellTree[0] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[1] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[2] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[3] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[4] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);

  CellTree[0]->SetClipBoard(Refinement);
  CellTree[1]->SetClipBoard(Refinement);
  CellTree[2]->SetClipBoard(Refinement);
  CellTree[3]->SetClipBoard(Refinement);
  CellTree[4]->SetClipBoard(Refinement);

  j[0] = new TJointEqN(CellTree[0], CellTree[4]);
  j[1] = new TJointEqN(CellTree[1], CellTree[4]);
  j[2] = new TJointEqN(CellTree[2], CellTree[4]);
  j[3] = new TJointEqN(CellTree[3], CellTree[4]);

  CellTree[0]->SetVertex(0, v[0]);
  CellTree[0]->SetVertex(1, v[1]);
  CellTree[0]->SetVertex(2, v[3]);

  CellTree[1]->SetVertex(0, v[1]);
  CellTree[1]->SetVertex(1, v[2]);
  CellTree[1]->SetVertex(2, v[4]);

  CellTree[2]->SetVertex(0, v[4]);
  CellTree[2]->SetVertex(1, v[7]);
  CellTree[2]->SetVertex(2, v[6]);

  CellTree[3]->SetVertex(0, v[3]);
  CellTree[3]->SetVertex(1, v[6]);
  CellTree[3]->SetVertex(2, v[5]);

  CellTree[4]->SetVertex(0, v[1]);
  CellTree[4]->SetVertex(1, v[4]);
  CellTree[4]->SetVertex(2, v[6]);
  CellTree[4]->SetVertex(3, v[3]);


  CellTree[0]->SetJoint(0, b[0]);
  CellTree[0]->SetJoint(1, j[0]);
  CellTree[0]->SetJoint(2, b[7]);

  CellTree[1]->SetJoint(0, b[1]);
  CellTree[1]->SetJoint(1, b[2]);
  CellTree[1]->SetJoint(2, j[1]);

  CellTree[2]->SetJoint(0, b[3]);
  CellTree[2]->SetJoint(1, b[4]);
  CellTree[2]->SetJoint(2, j[2]);

  CellTree[3]->SetJoint(0, j[3]);
  CellTree[3]->SetJoint(1, b[5]);
  CellTree[3]->SetJoint(2, b[6]);

  CellTree[4]->SetJoint(0, j[1]);
  CellTree[4]->SetJoint(1, j[2]);
  CellTree[4]->SetJoint(2, j[3]);
  CellTree[4]->SetJoint(3, j[0]);

}

void TDomain::SquareInSquareRef()
{
  TVertex *v[8];
  TBoundEdge *b[8];
  TJointEqN *j[4];

  v[0] = new TVertex(0.0, 0.0);
  v[1] = new TVertex(0.5, 0.0);
  v[2] = new TVertex(1.0, 0.0);
  v[3] = new TVertex(0.0, 0.5);
  v[4] = new TVertex(1.0, 0.5);
  v[5] = new TVertex(0.0, 1.0);
  v[6] = new TVertex(0.5, 1.0);
  v[7] = new TVertex(1.0, 1.0);

  b[0] = new TBoundEdge(BdParts[0]->GetBdComp(0), 0.0, 0.5);
  b[1] = new TBoundEdge(BdParts[0]->GetBdComp(0), 0.5, 1.0);
  b[2] = new TBoundEdge(BdParts[0]->GetBdComp(1), 0.0, 0.5);
  b[3] = new TBoundEdge(BdParts[0]->GetBdComp(1), 0.5, 1.0);
  b[4] = new TBoundEdge(BdParts[0]->GetBdComp(2), 0.0, 0.5);
  b[5] = new TBoundEdge(BdParts[0]->GetBdComp(2), 0.5, 1.0);
  b[6] = new TBoundEdge(BdParts[0]->GetBdComp(3), 0.0, 0.5);
  b[7] = new TBoundEdge(BdParts[0]->GetBdComp(3), 0.5, 1.0);

  CellTree = new TBaseCell*[5];
  N_InitVCells = 5;
  N_RootCells = 5;

  // set bounding box for output
  SetBoundBoxstart(0., 0.);
  SetBoundBox(1, 1);

  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);

  CellTree[0] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[1] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[2] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[3] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[4] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);

  CellTree[0]->SetClipBoard(Refinement);
  CellTree[1]->SetClipBoard(Refinement);
  CellTree[2]->SetClipBoard(Refinement);
  CellTree[3]->SetClipBoard(Refinement);
  CellTree[4]->SetClipBoard(Refinement);

  j[0] = new TJointEqN(CellTree[0], CellTree[4]);
  j[1] = new TJointEqN(CellTree[1], CellTree[4]);
  j[2] = new TJointEqN(CellTree[2], CellTree[4]);
  j[3] = new TJointEqN(CellTree[3], CellTree[4]);

  CellTree[0]->SetVertex(0, v[0]);
  CellTree[0]->SetVertex(1, v[1]);
  CellTree[0]->SetVertex(2, v[3]);

  CellTree[1]->SetVertex(0, v[1]);
  CellTree[1]->SetVertex(1, v[2]);
  CellTree[1]->SetVertex(2, v[4]);

  CellTree[2]->SetVertex(0, v[4]);
  CellTree[2]->SetVertex(1, v[7]);
  CellTree[2]->SetVertex(2, v[6]);

  CellTree[3]->SetVertex(0, v[3]);
  CellTree[3]->SetVertex(1, v[6]);
  CellTree[3]->SetVertex(2, v[5]);

  CellTree[4]->SetVertex(0, v[1]);
  CellTree[4]->SetVertex(1, v[4]);
  CellTree[4]->SetVertex(2, v[6]);
  CellTree[4]->SetVertex(3, v[3]);


  CellTree[0]->SetJoint(0, b[0]);
  CellTree[0]->SetJoint(1, j[0]);
  CellTree[0]->SetJoint(2, b[7]);

  CellTree[1]->SetJoint(0, b[1]);
  CellTree[1]->SetJoint(1, b[2]);
  CellTree[1]->SetJoint(2, j[1]);

  CellTree[2]->SetJoint(0, b[3]);
  CellTree[2]->SetJoint(1, b[4]);
  CellTree[2]->SetJoint(2, j[2]);

  CellTree[3]->SetJoint(0, j[3]);
  CellTree[3]->SetJoint(1, b[5]);
  CellTree[3]->SetJoint(2, b[6]);

  CellTree[4]->SetJoint(0, j[1]);
  CellTree[4]->SetJoint(1, j[2]);
  CellTree[4]->SetJoint(2, j[3]);
  CellTree[4]->SetJoint(3, j[0]);

  CellTree[0]->Refine(RefLevel);
  CellTree[1]->Refine(RefLevel);
  CellTree[2]->Refine(RefLevel);
  CellTree[3]->Refine(RefLevel);
}

void TDomain::UnitSquare()
{
  TVertex *v[4];
  TBoundEdge *b[4];

  v[0] = new TVertex(0.0, 0.0);
  v[1] = new TVertex(1.0, 0.0);
  v[2] = new TVertex(1.0, 1.0);
  v[3] = new TVertex(0.0, 1.0);

  b[0] = new TBoundEdge(BdParts[0]->GetBdComp(0), 0.0, 1.0);
  b[1] = new TBoundEdge(BdParts[0]->GetBdComp(1), 0.0, 1.0);
  b[2] = new TBoundEdge(BdParts[0]->GetBdComp(2), 0.0, 1.0);
  b[3] = new TBoundEdge(BdParts[0]->GetBdComp(3), 0.0, 1.0);

  CellTree = new TBaseCell*[1];
  N_InitVCells = 1;
  N_RootCells = 1;

  // set bounding box for output
  SetBoundBoxstart(0., 0.);
  SetBoundBox(1, 1);

  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);

  CellTree[0] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);

  CellTree[0]->SetClipBoard(Refinement);

  CellTree[0]->SetVertex(0, v[0]);
  CellTree[0]->SetVertex(1, v[1]);
  CellTree[0]->SetVertex(2, v[2]);
  CellTree[0]->SetVertex(3, v[3]);

  CellTree[0]->SetJoint(0, b[0]);
  CellTree[0]->SetJoint(1, b[1]);
  CellTree[0]->SetJoint(2, b[2]);
  CellTree[0]->SetJoint(3, b[3]);

}

void TDomain::DrivenCavitySquareQuads()
{
  TVertex *v[4];
  TBoundEdge *b[4];

  v[0] = new TVertex(-1.0, -1.0);
  v[1] = new TVertex(1.0, -1.0);
  v[2] = new TVertex(1.0, 1.0);
  v[3] = new TVertex(-1.0, 1.0);

  b[0] = new TBoundEdge(BdParts[0]->GetBdComp(0), 0.0, 1.0);
  b[1] = new TBoundEdge(BdParts[0]->GetBdComp(1), 0.0, 1.0);
  b[2] = new TBoundEdge(BdParts[0]->GetBdComp(2), 0.0, 1.0);
  b[3] = new TBoundEdge(BdParts[0]->GetBdComp(3), 0.0, 1.0);

  CellTree = new TBaseCell*[1];
  N_InitVCells = 1;
  N_RootCells = 1;

  // set bounding box for output
  SetBoundBoxstart(-1., -1.);
  SetBoundBox(1, 1);

  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);

  CellTree[0] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);

  CellTree[0]->SetClipBoard(Refinement);

  CellTree[0]->SetVertex(0, v[0]);
  CellTree[0]->SetVertex(1, v[1]);
  CellTree[0]->SetVertex(2, v[2]);
  CellTree[0]->SetVertex(3, v[3]);

  CellTree[0]->SetJoint(0, b[0]);
  CellTree[0]->SetJoint(1, b[1]);
  CellTree[0]->SetJoint(2, b[2]);
  CellTree[0]->SetJoint(3, b[3]);

}

void TDomain::UnitSquareRef()
{
  TVertex *v[4];
  TBoundEdge *b[4];

  v[0] = new TVertex(0.0, 0.0);
  v[1] = new TVertex(1.0, 0.0);
  v[2] = new TVertex(1.0, 1.0);
  v[3] = new TVertex(0.0, 1.0);

  b[0] = new TBoundEdge(BdParts[0]->GetBdComp(0), 0.0, 1.0);
  b[1] = new TBoundEdge(BdParts[0]->GetBdComp(1), 0.0, 1.0);
  b[2] = new TBoundEdge(BdParts[0]->GetBdComp(2), 0.0, 1.0);
  b[3] = new TBoundEdge(BdParts[0]->GetBdComp(3), 0.0, 1.0);

  CellTree = new TBaseCell*[1];
  N_InitVCells = 1;
  N_RootCells = 1;

  // set bounding box for output
  SetBoundBoxstart(0., 0.);
  SetBoundBox(1, 1);

  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);

  CellTree[0] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);

  CellTree[0]->SetClipBoard(Refinement);

  CellTree[0]->SetVertex(0, v[0]);
  CellTree[0]->SetVertex(1, v[1]);
  CellTree[0]->SetVertex(2, v[2]);
  CellTree[0]->SetVertex(3, v[3]);

  CellTree[0]->SetJoint(0, b[0]);
  CellTree[0]->SetJoint(1, b[1]);
  CellTree[0]->SetJoint(2, b[2]);
  CellTree[0]->SetJoint(3, b[3]);

  CellTree[0]->Refine(RefLevel);

  CellTree[0]->GetChild(1)->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + QuadReg]);

  CellTree[0]->GetChild(1)->Refine(RefLevel);

}


void TDomain::UnitSquare_US22()
{
  TVertex *v[9];
  TBoundEdge *b[8];
  TJointEqN *j[4];
  // vertices of the unit square
  v[0] = new TVertex(0.0, 0.0);
  v[1] = new TVertex(0.5, 0.0);
  v[2] = new TVertex(1.0, 0.0);
  v[3] = new TVertex(0.0, 0.5);
  v[4] = new TVertex(0.6, 0.6);
  v[5] = new TVertex(1.0, 0.5);
  v[6] = new TVertex(0.0, 1.0);
  v[7] = new TVertex(0.5, 1.0);
  v[8] = new TVertex(1.0, 1.0);

   // edges of the square which are on the boundary
  b[0] = new TBoundEdge(BdParts[0]->GetBdComp(0), 0.0, 0.5);
  b[1] = new TBoundEdge(BdParts[0]->GetBdComp(0), 0.5, 1.0);
  
  b[2] = new TBoundEdge(BdParts[0]->GetBdComp(1), 0.0, 0.5);
  b[3] = new TBoundEdge(BdParts[0]->GetBdComp(1), 0.5, 1.0);  
  
  b[4] = new TBoundEdge(BdParts[0]->GetBdComp(2), 0.0, 0.5);
  b[5] = new TBoundEdge(BdParts[0]->GetBdComp(2), 0.5, 1.0);  
  
  b[6] = new TBoundEdge(BdParts[0]->GetBdComp(3), 0.0, 0.5);
  b[7] = new TBoundEdge(BdParts[0]->GetBdComp(3), 0.5, 1.0);

  CellTree = new TBaseCell*[4];
  N_InitVCells = 4;
  N_RootCells = 4;

  // set bounding box for output
  SetBoundBoxstart(0., 0.);
  SetBoundBox(1, 1);

  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);
  
  CellTree[0] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);
  CellTree[1] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);
  CellTree[2] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);
  CellTree[3] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);
  
  // include four squares into the cell tree
  CellTree[0]->SetClipBoard(Refinement);
  CellTree[1]->SetClipBoard(Refinement);
  CellTree[2]->SetClipBoard(Refinement);
  CellTree[3]->SetClipBoard(Refinement);
  

  // allocate joint (edge) between the triangles  
  j[0] = new TJointEqN(CellTree[0], CellTree[1]);
  j[1] = new TJointEqN(CellTree[1], CellTree[2]);
  j[2] = new TJointEqN(CellTree[2], CellTree[3]);
  j[3] = new TJointEqN(CellTree[3], CellTree[0]);

  // set the vertices of square 0 at (0,0), (0.5,0)
  CellTree[0]->SetVertex(0, v[0]);
  CellTree[0]->SetVertex(1, v[1]);
  CellTree[0]->SetVertex(2, v[4]);
  CellTree[0]->SetVertex(3, v[3]);
  // set the vertices of square 1 
  CellTree[1]->SetVertex(0, v[2]);
  CellTree[1]->SetVertex(1, v[5]);
  CellTree[1]->SetVertex(2, v[4]);
  CellTree[1]->SetVertex(3, v[1]);
  // set the vertices of square 2 
  CellTree[2]->SetVertex(0, v[8]);
  CellTree[2]->SetVertex(1, v[7]);
  CellTree[2]->SetVertex(2, v[4]);
  CellTree[2]->SetVertex(3, v[5]);
  // set the vertices of square 3
  CellTree[3]->SetVertex(0, v[6]);
  CellTree[3]->SetVertex(1, v[3]);
  CellTree[3]->SetVertex(2, v[4]);
  CellTree[3]->SetVertex(3, v[7]);

  // s0
  CellTree[0]->SetJoint(0, b[0]);
  CellTree[0]->SetJoint(1, j[0]);
  CellTree[0]->SetJoint(2, j[3]);
  CellTree[0]->SetJoint(3, b[7]);
  // s1 
  CellTree[1]->SetJoint(0, b[2]);
  CellTree[1]->SetJoint(1, j[1]);
  CellTree[1]->SetJoint(2, j[0]);
  CellTree[1]->SetJoint(3, b[1]);
  // s2 
  CellTree[2]->SetJoint(0, b[4]);
  CellTree[2]->SetJoint(1, j[2]);
  CellTree[2]->SetJoint(2, j[1]);
  CellTree[2]->SetJoint(3, b[3]);
  // s3
  CellTree[3]->SetJoint(0, b[6]);
  CellTree[3]->SetJoint(1, j[3]);
  CellTree[3]->SetJoint(2, j[2]);
  CellTree[3]->SetJoint(3, b[5]);
}

void TDomain::TestGrid2()
{
  TBaseCell *LocCell;

  CellTree[0]->GetChild(1)->SetRefDesc(TDatabase::
                  RefDescDB[N_SHAPES + QuadReg]);
  CellTree[1]->GetChild(0)->SetRefDesc(TDatabase::
                  RefDescDB[N_SHAPES + TriReg]);

  LocCell = CellTree[1]->GetChild(3);
  LocCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TriReg]);
  LocCell->Refine(RefLevel);

  LocCell->GetChild(2)->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TriReg]);
}

void TDomain::TestGrid3()
{
  CellTree[1]->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TriReg]);
}

void TDomain::TestShishkin()
{
  CellTree[0]->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + QuadToTri0]);
  CellTree[1]->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + QuadToTri0]);
  CellTree[2]->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + QuadToTri0]);
  CellTree[3]->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + QuadToTri1]);

  Refine();

  CellTree[0]->GetChild(0)->SetClipBoard(0);
  CellTree[0]->GetChild(1)->SetClipBoard(1);
  CellTree[1]->GetChild(0)->SetClipBoard(2);
  CellTree[1]->GetChild(1)->SetClipBoard(3);
  CellTree[2]->GetChild(0)->SetClipBoard(4);
  CellTree[2]->GetChild(1)->SetClipBoard(5);
  CellTree[3]->GetChild(0)->SetClipBoard(6);
  CellTree[3]->GetChild(1)->SetClipBoard(7);

  int i, j;

  for (i=0;i<4;i++)
  for (j=0;j<2;j++)
  cout << " 1 " << CellTree[i]->GetChild(j)->GetJoint(0)->GetNeighbour(CellTree[i]->GetChild(j))->GetClipBoard() << CellTree[i]->GetChild(j)->GetJoint(1)->GetNeighbour(CellTree[i]->GetChild(j))->GetClipBoard() << CellTree[i]->GetChild(j)->GetJoint(2)->GetNeighbour(CellTree[i]->GetChild(j))->GetClipBoard() << endl;
}

/*
double indicator(double x, double y, double a)
{
  if (x < (a-0.1)*3.3) return -1;
  if (x > (a-0.1)*3.3) return 1;

  return 0;
}
*/
///*
double indicator(double x, double y, double a)
{
  double xa=x-0.5;
  double ya=(1-y)-(0.5-0.875*a);

  return (xa*xa+ya*ya)*(xa*xa+ya*ya-2*a*ya)-a*a*xa*xa;
}
//*/

void TDomain::RefCardioide(double A)
{
  TBaseCell *CurrCell;
  TVertex *CurrVertex;
  int i, info, N_, inner, outer;
  double indi;

  TDatabase::IteratorDB[It_Between]->Init(0);

  // loop over all cells
  while ((CurrCell = TDatabase::IteratorDB[It_Between]->Next(info)))
  {
    N_ = CurrCell->GetN_Vertices();
    for (inner=outer=i=0;i<N_;i++)
    {
      CurrVertex = CurrCell->GetVertex(i);
      indi = indicator(CurrVertex->GetX(), CurrVertex->GetY(), A);
      if (indi >= 0 ) inner++;
      if (indi <= 0 ) outer++;
    }

    if (inner && outer)
    {
      if (info != 1)
        CurrCell->SetClipBoard(Refinement);
      else
        CurrCell->SetClipBoard(NoRefinement);
    }
    else
      if (info != -1)
        CurrCell->SetClipBoard(DeRefinement);
      else
        CurrCell->SetClipBoard(NoRefinement);
  }

  TDatabase::IteratorDB[It_Finest]->Init(0);
  for (i=TDatabase::IteratorDB[It_Finest]->GetMaxLevel();i>0;i--)
  {
    TDatabase::IteratorDB[It_EQ]->Init(i);
    while ((CurrCell = TDatabase::IteratorDB[It_EQ]->Next(info)))
      CurrCell->Gen1RegMarks();
  }

  TDatabase::IteratorDB[It_Finest]->Init(0);
  while ((CurrCell = TDatabase::IteratorDB[It_Finest]->Next(info)))
  {
    switch (CurrCell->GetClipBoard())
    {
      case DeRefinement:
        CurrCell = CurrCell->GetParent();
        N_ = CurrCell->GetN_Children();
        for (info=i=0;i<N_;i++)
          if (CurrCell->GetChild(i)->GetClipBoard() != DeRefinement)
          {
            info = 1;
            break;
          }

        if (info)
        {
          for (i=0;i<N_;i++)
            if (CurrCell->GetChild(i)->GetClipBoard() == DeRefinement)
              CurrCell->GetChild(i)->SetClipBoard(NoRefinement);
        }
        else
        {
          CurrCell->SetClipBoard(NoRefinement);
          CurrCell->Derefine();
        }
        break;

      case Refinement:
        CurrCell->SetRegRefine();
        CurrCell->Refine(0);
        break;
    }
  }
}

void TDomain::PeriodicSquares()
{
  TVertex *v[9];
  TBoundEdge *b[4];
  TJoint *j[6];

  TDatabase::ParamDB->INTERNAL_PERIODIC_IDENTITY = 1;

  v[0] = new TVertex(0.0, 0.0);
  v[1] = new TVertex(0.5, 0.0);
  v[2] = new TVertex(1.0, 0.0);
  v[3] = new TVertex(0.0, 0.5);
  v[4] = new TVertex(0.5, 0.5);
  v[5] = new TVertex(1.0, 0.5);
  v[6] = new TVertex(0.0, 1.0);
  v[7] = new TVertex(0.5, 1.0);
  v[8] = new TVertex(1.0, 1.0);

  b[0] = new TBoundEdge(BdParts[0]->GetBdComp(0), 0.0, 0.5);
  b[1] = new TBoundEdge(BdParts[0]->GetBdComp(0), 0.5, 1.0);
  b[2] = new TBoundEdge(BdParts[0]->GetBdComp(2), 0.0, 0.5);
  b[3] = new TBoundEdge(BdParts[0]->GetBdComp(2), 0.5, 1.0);

  CellTree = new TBaseCell*[4];
  N_InitVCells = 4;
  N_RootCells = 4;

  SetBoundBox(1, 1);

  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);

  CellTree[0] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);
  CellTree[1] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);
  CellTree[2] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);
  CellTree[3] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);

  CellTree[0]->SetClipBoard(Refinement);
  CellTree[1]->SetClipBoard(Refinement);
  CellTree[2]->SetClipBoard(Refinement);
  CellTree[3]->SetClipBoard(Refinement);

  j[0] = new TJointEqN(CellTree[0], CellTree[1]);
  j[1] = new TJointEqN(CellTree[1], CellTree[2]);
  j[2] = new TJointEqN(CellTree[2], CellTree[3]);
  j[3] = new TJointEqN(CellTree[3], CellTree[0]);
  j[4] = new TPeriodicJoint(CellTree[1], CellTree[0]);
  j[5] = new TPeriodicJoint(CellTree[2], CellTree[3]);

  CellTree[0]->SetVertex(0, v[0]);
  CellTree[0]->SetVertex(1, v[1]);
  CellTree[0]->SetVertex(2, v[4]);
  CellTree[0]->SetVertex(3, v[3]);

  CellTree[1]->SetVertex(0, v[2]);
  CellTree[1]->SetVertex(1, v[5]);
  CellTree[1]->SetVertex(2, v[4]);
  CellTree[1]->SetVertex(3, v[1]);

  CellTree[2]->SetVertex(0, v[8]);
  CellTree[2]->SetVertex(1, v[7]);
  CellTree[2]->SetVertex(2, v[4]);
  CellTree[2]->SetVertex(3, v[5]);

  CellTree[3]->SetVertex(0, v[6]);
  CellTree[3]->SetVertex(1, v[3]);
  CellTree[3]->SetVertex(2, v[4]);
  CellTree[3]->SetVertex(3, v[7]);

  CellTree[0]->SetJoint(0, b[0]);
  CellTree[0]->SetJoint(1, j[0]);
  CellTree[0]->SetJoint(2, j[3]);
  CellTree[0]->SetJoint(3, j[4]);

  CellTree[1]->SetJoint(0, j[4]);
  CellTree[1]->SetJoint(1, j[1]);
  CellTree[1]->SetJoint(2, j[0]);
  CellTree[1]->SetJoint(3, b[1]);

  CellTree[2]->SetJoint(0, b[2]);
  CellTree[2]->SetJoint(1, j[2]);
  CellTree[2]->SetJoint(2, j[1]);
  CellTree[2]->SetJoint(3, j[5]);

  CellTree[3]->SetJoint(0, j[5]);
  CellTree[3]->SetJoint(1, j[3]);
  CellTree[3]->SetJoint(2, j[2]);
  CellTree[3]->SetJoint(3, b[3]);

}

void TDomain::PeriodicTriangles()
{
  TVertex *v[9];
  TBoundEdge *b[4];
  TJoint *j[10];

  TDatabase::ParamDB->INTERNAL_PERIODIC_IDENTITY = 1;

  v[0] = new TVertex(0.0, 0.0);
  v[1] = new TVertex(0.5, 0.0);
  v[2] = new TVertex(1.0, 0.0);
  v[3] = new TVertex(0.0, 0.5);
  v[4] = new TVertex(0.5, 0.5);
  v[5] = new TVertex(1.0, 0.5);
  v[6] = new TVertex(0.0, 1.0);
  v[7] = new TVertex(0.5, 1.0);
  v[8] = new TVertex(1.0, 1.0);

  b[0] = new TBoundEdge(BdParts[0]->GetBdComp(0), 0.0, 0.5);
  b[1] = new TBoundEdge(BdParts[0]->GetBdComp(0), 0.5, 1.0);
  b[2] = new TBoundEdge(BdParts[0]->GetBdComp(2), 0.0, 0.5);
  b[3] = new TBoundEdge(BdParts[0]->GetBdComp(2), 0.5, 1.0);

  CellTree = new TBaseCell*[8];
  N_InitVCells = 8;
  N_RootCells = 8;

  SetBoundBox(2, 2);

  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);

  CellTree[0] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[1] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[2] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[3] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[4] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[5] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[6] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[7] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);

  CellTree[0]->SetClipBoard(Refinement);
  CellTree[1]->SetClipBoard(Refinement);
  CellTree[2]->SetClipBoard(Refinement);
  CellTree[3]->SetClipBoard(Refinement);
  CellTree[4]->SetClipBoard(Refinement);
  CellTree[5]->SetClipBoard(Refinement);
  CellTree[6]->SetClipBoard(Refinement);
  CellTree[7]->SetClipBoard(Refinement);

  j[0] = new TJointEqN(CellTree[0], CellTree[2]);
  j[1] = new TJointEqN(CellTree[0], CellTree[1]);
  j[2] = new TJointEqN(CellTree[1], CellTree[3]);
  j[3] = new TJointEqN(CellTree[2], CellTree[4]);
  j[4] = new TJointEqN(CellTree[3], CellTree[5]);
  j[5] = new TJointEqN(CellTree[4], CellTree[5]);
  j[6] = new TJointEqN(CellTree[4], CellTree[6]);
  j[7] = new TJointEqN(CellTree[5], CellTree[7]);
  j[8] = new TPeriodicJoint(CellTree[7], CellTree[1]);
  j[9] = new TPeriodicJoint(CellTree[6], CellTree[0]);

  // assign vertices to the mesh cells
  CellTree[0]->SetVertex(0, v[0]);
  CellTree[0]->SetVertex(1, v[4]);
  CellTree[0]->SetVertex(2, v[3]);
 
  CellTree[1]->SetVertex(0, v[3]);
  CellTree[1]->SetVertex(1, v[4]);
  CellTree[1]->SetVertex(2, v[6]);

  CellTree[2]->SetVertex(0, v[0]);
  CellTree[2]->SetVertex(1, v[1]);
  CellTree[2]->SetVertex(2, v[4]);

  CellTree[3]->SetVertex(0, v[4]);
  CellTree[3]->SetVertex(1, v[7]);
  CellTree[3]->SetVertex(2, v[6]);

  CellTree[4]->SetVertex(0, v[1]);
  CellTree[4]->SetVertex(1, v[5]);
  CellTree[4]->SetVertex(2, v[4]);
 
  CellTree[5]->SetVertex(0, v[4]);
  CellTree[5]->SetVertex(1, v[5]);
  CellTree[5]->SetVertex(2, v[7]);

  CellTree[6]->SetVertex(0, v[1]);
  CellTree[6]->SetVertex(1, v[2]);
  CellTree[6]->SetVertex(2, v[5]);

  CellTree[7]->SetVertex(0, v[5]);
  CellTree[7]->SetVertex(1, v[8]);
  CellTree[7]->SetVertex(2, v[7]);
 
  // assign links to the edges
  CellTree[0]->SetJoint(0, j[0]);
  CellTree[0]->SetJoint(1, j[1]);
  CellTree[0]->SetJoint(2, j[9]);
 
  CellTree[1]->SetJoint(0, j[1]);
  CellTree[1]->SetJoint(1, j[2]);
  CellTree[1]->SetJoint(2, j[8]);

  CellTree[2]->SetJoint(0, b[0]);
  CellTree[2]->SetJoint(1, j[3]);
  CellTree[2]->SetJoint(2, j[0]);

  CellTree[3]->SetJoint(0, j[4]);
  CellTree[3]->SetJoint(1, b[3]);
  CellTree[3]->SetJoint(2, j[2]);
 
  CellTree[4]->SetJoint(0, j[6]);
  CellTree[4]->SetJoint(1, j[5]);
  CellTree[4]->SetJoint(2, j[3]);
 
  CellTree[5]->SetJoint(0, j[5]);
  CellTree[5]->SetJoint(1, j[7]);
  CellTree[5]->SetJoint(2, j[4]);

  CellTree[6]->SetJoint(0, b[1]);
  CellTree[6]->SetJoint(1, j[9]);
  CellTree[6]->SetJoint(2, j[6]);

  CellTree[7]->SetJoint(0, j[8]);
  CellTree[7]->SetJoint(1, b[2]);
  CellTree[7]->SetJoint(2, j[7]);
}

void TDomain::PeriodicSquaresLarge()
{
  TVertex *v[9];
  TBoundEdge *b[4];
  TJoint *j[6];

  TDatabase::ParamDB->INTERNAL_PERIODIC_IDENTITY = 2;

  v[0] = new TVertex(-1.0, -1.0);
  v[1] = new TVertex(0, -1.0);
  v[2] = new TVertex(1.0, -1.0);
  v[3] = new TVertex(-1.0, 0.0);
  v[4] = new TVertex(0, 0.0);
  v[5] = new TVertex(1.0, 0.0);
  v[6] = new TVertex(-1.0, 1.0);
  v[7] = new TVertex(0, 1.0);
  v[8] = new TVertex(1.0, 1.0);

  b[0] = new TBoundEdge(BdParts[0]->GetBdComp(0), 0.0, 0.5);
  b[1] = new TBoundEdge(BdParts[0]->GetBdComp(0), 0.5, 1.0);
  b[2] = new TBoundEdge(BdParts[0]->GetBdComp(2), 0.0, 0.5);
  b[3] = new TBoundEdge(BdParts[0]->GetBdComp(2), 0.5, 1.0);

  CellTree = new TBaseCell*[4];
  N_InitVCells = 4;
  N_RootCells = 4;

  SetBoundBox(2, 2);

  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);

  CellTree[0] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);
  CellTree[1] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);
  CellTree[2] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);
  CellTree[3] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);

  CellTree[0]->SetClipBoard(Refinement);
  CellTree[1]->SetClipBoard(Refinement);
  CellTree[2]->SetClipBoard(Refinement);
  CellTree[3]->SetClipBoard(Refinement);

  j[0] = new TJointEqN(CellTree[0], CellTree[1]);
  j[1] = new TJointEqN(CellTree[1], CellTree[2]);
  j[2] = new TJointEqN(CellTree[2], CellTree[3]);
  j[3] = new TJointEqN(CellTree[3], CellTree[0]);
  j[4] = new TPeriodicJoint(CellTree[1], CellTree[0]);
  j[5] = new TPeriodicJoint(CellTree[2], CellTree[3]);

  CellTree[0]->SetVertex(0, v[0]);
  CellTree[0]->SetVertex(1, v[1]);
  CellTree[0]->SetVertex(2, v[4]);
  CellTree[0]->SetVertex(3, v[3]);

  CellTree[1]->SetVertex(0, v[2]);
  CellTree[1]->SetVertex(1, v[5]);
  CellTree[1]->SetVertex(2, v[4]);
  CellTree[1]->SetVertex(3, v[1]);

  CellTree[2]->SetVertex(0, v[8]);
  CellTree[2]->SetVertex(1, v[7]);
  CellTree[2]->SetVertex(2, v[4]);
  CellTree[2]->SetVertex(3, v[5]);

  CellTree[3]->SetVertex(0, v[6]);
  CellTree[3]->SetVertex(1, v[3]);
  CellTree[3]->SetVertex(2, v[4]);
  CellTree[3]->SetVertex(3, v[7]);

  CellTree[0]->SetJoint(0, b[0]);
  CellTree[0]->SetJoint(1, j[0]);
  CellTree[0]->SetJoint(2, j[3]);
  CellTree[0]->SetJoint(3, j[4]);

  CellTree[1]->SetJoint(0, j[4]);
  CellTree[1]->SetJoint(1, j[1]);
  CellTree[1]->SetJoint(2, j[0]);
  CellTree[1]->SetJoint(3, b[1]);

  CellTree[2]->SetJoint(0, b[2]);
  CellTree[2]->SetJoint(1, j[2]);
  CellTree[2]->SetJoint(2, j[1]);
  CellTree[2]->SetJoint(3, j[5]);

  CellTree[3]->SetJoint(0, j[5]);
  CellTree[3]->SetJoint(1, j[3]);
  CellTree[3]->SetJoint(2, j[2]);
  CellTree[3]->SetJoint(3, b[3]);

}

void TDomain::PeriodicTrianglesLarge()
{
  TVertex *v[9];
  TBoundEdge *b[4];
  TJoint *j[10];

  TDatabase::ParamDB->INTERNAL_PERIODIC_IDENTITY = 2;

  v[0] = new TVertex(-1.0, -1.0);
  v[1] = new TVertex(0, -1.0);
  v[2] = new TVertex(1.0, -1.0);
  v[3] = new TVertex(-1.0, 0.0);
  v[4] = new TVertex(0, 0.0);
  v[5] = new TVertex(1.0, 0.0);
  v[6] = new TVertex(-1.0, 1.0);
  v[7] = new TVertex(0, 1.0);
  v[8] = new TVertex(1.0, 1.0);

  b[0] = new TBoundEdge(BdParts[0]->GetBdComp(0), 0.0, 0.5);
  b[1] = new TBoundEdge(BdParts[0]->GetBdComp(0), 0.5, 1.0);
  b[2] = new TBoundEdge(BdParts[0]->GetBdComp(2), 0.0, 0.5);
  b[3] = new TBoundEdge(BdParts[0]->GetBdComp(2), 0.5, 1.0);

  CellTree = new TBaseCell*[8];
  N_InitVCells = 8;
  N_RootCells = 8;

  SetBoundBox(2, 2);

  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);

  CellTree[0] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[1] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[2] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[3] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[4] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[5] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[6] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);
  CellTree[7] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      TriReg], RefLevel);

  CellTree[0]->SetClipBoard(Refinement);
  CellTree[1]->SetClipBoard(Refinement);
  CellTree[2]->SetClipBoard(Refinement);
  CellTree[3]->SetClipBoard(Refinement);
  CellTree[4]->SetClipBoard(Refinement);
  CellTree[5]->SetClipBoard(Refinement);
  CellTree[6]->SetClipBoard(Refinement);
  CellTree[7]->SetClipBoard(Refinement);

  j[0] = new TJointEqN(CellTree[0], CellTree[2]);
  j[1] = new TJointEqN(CellTree[0], CellTree[1]);
  j[2] = new TJointEqN(CellTree[1], CellTree[3]);
  j[3] = new TJointEqN(CellTree[2], CellTree[4]);
  j[4] = new TJointEqN(CellTree[3], CellTree[5]);
  j[5] = new TJointEqN(CellTree[4], CellTree[5]);
  j[6] = new TJointEqN(CellTree[4], CellTree[6]);
  j[7] = new TJointEqN(CellTree[5], CellTree[7]);
  j[8] = new TPeriodicJoint(CellTree[7], CellTree[1]);
  j[9] = new TPeriodicJoint(CellTree[6], CellTree[0]);

  // assign vertices to the mesh cells
  CellTree[0]->SetVertex(0, v[0]);
  CellTree[0]->SetVertex(1, v[4]);
  CellTree[0]->SetVertex(2, v[3]);
 
  CellTree[1]->SetVertex(0, v[3]);
  CellTree[1]->SetVertex(1, v[4]);
  CellTree[1]->SetVertex(2, v[6]);

  CellTree[2]->SetVertex(0, v[0]);
  CellTree[2]->SetVertex(1, v[1]);
  CellTree[2]->SetVertex(2, v[4]);

  CellTree[3]->SetVertex(0, v[4]);
  CellTree[3]->SetVertex(1, v[7]);
  CellTree[3]->SetVertex(2, v[6]);

  CellTree[4]->SetVertex(0, v[1]);
  CellTree[4]->SetVertex(1, v[5]);
  CellTree[4]->SetVertex(2, v[4]);
 
  CellTree[5]->SetVertex(0, v[4]);
  CellTree[5]->SetVertex(1, v[5]);
  CellTree[5]->SetVertex(2, v[7]);

  CellTree[6]->SetVertex(0, v[1]);
  CellTree[6]->SetVertex(1, v[2]);
  CellTree[6]->SetVertex(2, v[5]);

  CellTree[7]->SetVertex(0, v[5]);
  CellTree[7]->SetVertex(1, v[8]);
  CellTree[7]->SetVertex(2, v[7]);
 
  // assign links to the edges
  CellTree[0]->SetJoint(0, j[0]);
  CellTree[0]->SetJoint(1, j[1]);
  CellTree[0]->SetJoint(2, j[9]);
 
  CellTree[1]->SetJoint(0, j[1]);
  CellTree[1]->SetJoint(1, j[2]);
  CellTree[1]->SetJoint(2, j[8]);

  CellTree[2]->SetJoint(0, b[0]);
  CellTree[2]->SetJoint(1, j[3]);
  CellTree[2]->SetJoint(2, j[0]);

  CellTree[3]->SetJoint(0, j[4]);
  CellTree[3]->SetJoint(1, b[3]);
  CellTree[3]->SetJoint(2, j[2]);
 
  CellTree[4]->SetJoint(0, j[6]);
  CellTree[4]->SetJoint(1, j[5]);
  CellTree[4]->SetJoint(2, j[3]);
 
  CellTree[5]->SetJoint(0, j[5]);
  CellTree[5]->SetJoint(1, j[7]);
  CellTree[5]->SetJoint(2, j[4]);

  CellTree[6]->SetJoint(0, b[1]);
  CellTree[6]->SetJoint(1, j[9]);
  CellTree[6]->SetJoint(2, j[6]);

  CellTree[7]->SetJoint(0, j[8]);
  CellTree[7]->SetJoint(1, b[2]);
  CellTree[7]->SetJoint(2, j[7]);

}

void TDomain::PeriodicRectangle_2_4()
{
  TVertex *v[15];
  TBoundEdge *b[4];
  TJoint *j[14];

  TDatabase::ParamDB->INTERNAL_PERIODIC_IDENTITY = 3;

  // set vertices of initial grid
  v[0] = new TVertex(-1.0, -2.0);
  v[1] = new TVertex(0, -2.0);
  v[2] = new TVertex(1.0, -2.0);
  v[3] = new TVertex(-1.0, -1.0);
  v[4] = new TVertex(0, -1.0);
  v[5] = new TVertex(1.0, -1.0);
  v[6] = new TVertex(-1.0, 0.0);
  v[7] = new TVertex(0, 0.0);
  v[8] = new TVertex(1.0, 0.0);
  v[9] = new TVertex(-1.0, 1.0);
  v[10] = new TVertex(0, 1.0);
  v[11] = new TVertex(1.0, 1.0);
  v[12] = new TVertex(-1.0, 2.0);
  v[13] = new TVertex(0, 2.0);
  v[14] = new TVertex(1.0, 2.0);

  // create boundary edges at y=-1 and y=2
  b[0] = new TBoundEdge(BdParts[0]->GetBdComp(0), 0.0, 0.5);
  b[1] = new TBoundEdge(BdParts[0]->GetBdComp(0), 0.5, 1.0);
  b[2] = new TBoundEdge(BdParts[0]->GetBdComp(2), 0.0, 0.5);
  b[3] = new TBoundEdge(BdParts[0]->GetBdComp(2), 0.5, 1.0);

  CellTree = new TBaseCell*[8];
  N_InitVCells = 8;
  N_RootCells = 8;

  SetBoundBox(2, 4);

  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);

  CellTree[0] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);
  CellTree[1] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);
  CellTree[2] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);
  CellTree[3] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);
  CellTree[4] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);
  CellTree[5] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);
  CellTree[6] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);
  CellTree[7] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);

  CellTree[0]->SetClipBoard(Refinement);
  CellTree[1]->SetClipBoard(Refinement);
  CellTree[2]->SetClipBoard(Refinement);
  CellTree[3]->SetClipBoard(Refinement);
  CellTree[4]->SetClipBoard(Refinement);
  CellTree[5]->SetClipBoard(Refinement);
  CellTree[6]->SetClipBoard(Refinement);
  CellTree[7]->SetClipBoard(Refinement);

  // make joints between the meshcells
  j[0] = new TJointEqN(CellTree[0], CellTree[1]);
  j[1] = new TJointEqN(CellTree[1], CellTree[2]);
  j[2] = new TJointEqN(CellTree[2], CellTree[3]);
  j[3] = new TJointEqN(CellTree[3], CellTree[0]);
  j[4] = new TJointEqN(CellTree[3], CellTree[4]);
  j[5] = new TJointEqN(CellTree[2], CellTree[5]);
  j[6] = new TJointEqN(CellTree[4], CellTree[5]);
  j[7] = new TJointEqN(CellTree[4], CellTree[7]);
  j[8] = new TJointEqN(CellTree[5], CellTree[6]);
  j[9] = new TJointEqN(CellTree[6], CellTree[7]);
  j[10] = new TPeriodicJoint(CellTree[1], CellTree[0]);
  j[11] = new TPeriodicJoint(CellTree[2], CellTree[3]);
  j[12] = new TPeriodicJoint(CellTree[4], CellTree[5]);
  j[13] = new TPeriodicJoint(CellTree[6], CellTree[7]);

  // assign vertices to  mesh cells
  CellTree[0]->SetVertex(0, v[0]);
  CellTree[0]->SetVertex(1, v[1]);
  CellTree[0]->SetVertex(2, v[4]);
  CellTree[0]->SetVertex(3, v[3]);

  CellTree[1]->SetVertex(0, v[2]);
  CellTree[1]->SetVertex(1, v[5]);
  CellTree[1]->SetVertex(2, v[4]);
  CellTree[1]->SetVertex(3, v[1]);

  CellTree[2]->SetVertex(0, v[8]);
  CellTree[2]->SetVertex(1, v[7]);
  CellTree[2]->SetVertex(2, v[4]);
  CellTree[2]->SetVertex(3, v[5]);

  CellTree[3]->SetVertex(0, v[6]);
  CellTree[3]->SetVertex(1, v[3]);
  CellTree[3]->SetVertex(2, v[4]);
  CellTree[3]->SetVertex(3, v[7]);

  CellTree[4]->SetVertex(0, v[6]);
  CellTree[4]->SetVertex(1, v[7]);
  CellTree[4]->SetVertex(2, v[10]);
  CellTree[4]->SetVertex(3, v[9]);

  CellTree[5]->SetVertex(0, v[7]);
  CellTree[5]->SetVertex(1, v[8]);
  CellTree[5]->SetVertex(2, v[11]);
  CellTree[5]->SetVertex(3, v[10]);

  CellTree[6]->SetVertex(0, v[10]);
  CellTree[6]->SetVertex(1, v[11]);
  CellTree[6]->SetVertex(2, v[14]);
  CellTree[6]->SetVertex(3, v[13]);

  CellTree[7]->SetVertex(0, v[9]);
  CellTree[7]->SetVertex(1, v[10]);
  CellTree[7]->SetVertex(2, v[13]);
  CellTree[7]->SetVertex(3, v[12]);

  // assign joints to mesh cells
  // 0-th joint is joint with vertices no. 0 and 1
  // 1-st joint is joint with vertices no. 1 and 2 ...
  CellTree[0]->SetJoint(0, b[0]);
  CellTree[0]->SetJoint(1, j[0]);
  CellTree[0]->SetJoint(2, j[3]);
  CellTree[0]->SetJoint(3, j[10]);

  CellTree[1]->SetJoint(0, j[10]);
  CellTree[1]->SetJoint(1, j[1]);
  CellTree[1]->SetJoint(2, j[0]);
  CellTree[1]->SetJoint(3, b[1]);

  CellTree[2]->SetJoint(0, j[5]);
  CellTree[2]->SetJoint(1, j[2]);
  CellTree[2]->SetJoint(2, j[1]);
  CellTree[2]->SetJoint(3, j[11]);

  CellTree[3]->SetJoint(0, j[11]);
  CellTree[3]->SetJoint(1, j[3]);
  CellTree[3]->SetJoint(2, j[2]);
  CellTree[3]->SetJoint(3, j[4]);

  CellTree[4]->SetJoint(0, j[4]);
  CellTree[4]->SetJoint(1, j[6]);
  CellTree[4]->SetJoint(2, j[7]);
  CellTree[4]->SetJoint(3, j[12]);

  CellTree[5]->SetJoint(0, j[5]);
  CellTree[5]->SetJoint(1, j[12]);
  CellTree[5]->SetJoint(2, j[8]);
  CellTree[5]->SetJoint(3, j[6]);

  CellTree[6]->SetJoint(0, j[8]);
  CellTree[6]->SetJoint(1, j[13]);
  CellTree[6]->SetJoint(2, b[2]);
  CellTree[6]->SetJoint(3, j[9]);

  CellTree[7]->SetJoint(0, j[7]);
  CellTree[7]->SetJoint(1, j[9]);
  CellTree[7]->SetJoint(2, b[3]);
  CellTree[7]->SetJoint(3, j[13]);
}

void TDomain::QuadShishkin(double tau1, double tau2)
{
  TVertex *v[9];
  TBoundEdge *b[8];
  TJointEqN *j[4];

  v[0] = new TVertex( 0.0,  0.0);
  v[1] = new TVertex(tau1,  0.0);
  v[2] = new TVertex( 1.0,  0.0);
  v[3] = new TVertex( 0.0, tau2);
  v[4] = new TVertex(tau1, tau2);
  v[5] = new TVertex( 1.0, tau2);
  v[6] = new TVertex( 0.0,  1.0);
  v[7] = new TVertex(tau1,  1.0);
  v[8] = new TVertex( 1.0,  1.0);

  b[0] = new TBoundEdge(BdParts[0]->GetBdComp(0), 0.0, tau1);
  b[1] = new TBoundEdge(BdParts[0]->GetBdComp(0), tau1, 1.0);
  b[2] = new TBoundEdge(BdParts[0]->GetBdComp(1), 0.0, tau2);
  b[3] = new TBoundEdge(BdParts[0]->GetBdComp(1), tau2, 1.0);
  b[4] = new TBoundEdge(BdParts[0]->GetBdComp(2), 0.0, 1-tau1);
  b[5] = new TBoundEdge(BdParts[0]->GetBdComp(2), 1-tau1, 1.0);
  b[6] = new TBoundEdge(BdParts[0]->GetBdComp(3), 0.0, 1-tau2);
  b[7] = new TBoundEdge(BdParts[0]->GetBdComp(3), 1-tau2, 1.0);

  CellTree = new TBaseCell*[4];
  N_InitVCells = 4;
  N_RootCells = 4;

  SetBoundBox(1, 1);

  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);

  CellTree[0] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);
  CellTree[1] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);
  CellTree[2] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);
  CellTree[3] = new TGridCell(TDatabase::RefDescDB[N_SHAPES +
                      QuadReg], RefLevel);

  CellTree[0]->SetClipBoard(Refinement);
  CellTree[1]->SetClipBoard(Refinement);
  CellTree[2]->SetClipBoard(Refinement);
  CellTree[3]->SetClipBoard(Refinement);

  j[0] = new TJointEqN(CellTree[0], CellTree[1]);
  j[1] = new TJointEqN(CellTree[1], CellTree[2]);
  j[2] = new TJointEqN(CellTree[2], CellTree[3]);
  j[3] = new TJointEqN(CellTree[3], CellTree[0]);

  CellTree[0]->SetVertex(0, v[0]);
  CellTree[0]->SetVertex(1, v[1]);
  CellTree[0]->SetVertex(2, v[4]);
  CellTree[0]->SetVertex(3, v[3]);

  CellTree[1]->SetVertex(0, v[2]);
  CellTree[1]->SetVertex(1, v[5]);
  CellTree[1]->SetVertex(2, v[4]);
  CellTree[1]->SetVertex(3, v[1]);

  CellTree[2]->SetVertex(0, v[8]);
  CellTree[2]->SetVertex(1, v[7]);
  CellTree[2]->SetVertex(2, v[4]);
  CellTree[2]->SetVertex(3, v[5]);

  CellTree[3]->SetVertex(0, v[6]);
  CellTree[3]->SetVertex(1, v[3]);
  CellTree[3]->SetVertex(2, v[4]);
  CellTree[3]->SetVertex(3, v[7]);


  CellTree[0]->SetJoint(0, b[0]);
  CellTree[0]->SetJoint(1, j[0]);
  CellTree[0]->SetJoint(2, j[3]);
  CellTree[0]->SetJoint(3, b[7]);

  CellTree[1]->SetJoint(0, b[2]);
  CellTree[1]->SetJoint(1, j[1]);
  CellTree[1]->SetJoint(2, j[0]);
  CellTree[1]->SetJoint(3, b[1]);

  CellTree[2]->SetJoint(0, b[4]);
  CellTree[2]->SetJoint(1, j[2]);
  CellTree[2]->SetJoint(2, j[1]);
  CellTree[2]->SetJoint(3, b[3]);

  CellTree[3]->SetJoint(0, b[6]);
  CellTree[3]->SetJoint(1, j[3]);
  CellTree[3]->SetJoint(2, j[2]);
  CellTree[3]->SetJoint(3, b[5]);

}

void TDomain::Rectangular(int n, int m)
{
  int i,j;
  TVertex **Vertices;
  TJoint **Horizontal;
  TJoint **Vertical;
  TMacroCell *MCell;
  double x, y;

  Vertices = new TVertex*[(n+1)*(m+1)];
  Horizontal = new TJoint*[n*(m+1)];
  Vertical = new TJoint*[(n+1)*m];

  CellTree = new TBaseCell*[n*m];
  N_InitVCells = n*m;
  N_RootCells = n*m;

  for(i=0;i<n*m;i++)
    CellTree[i] = new TMacroCell(TDatabase::RefDescDB[N_SHAPES +
                          QuadReg], RefLevel);

  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);

  StartX = 0;
  StartY = -4;
  SetBoundBox(1,8);

  for(i=0;i<=n;i++)
    for(j=0;j<=m;j++)
    {
      x = ((double)i)/n;
      y = 8*((double)j)/m-4;
      Vertices[i+j*(n+1)] = new TVertex(x, y);
    }

  for(i=0;i<n;i++)
    for(j=1;j<m;j++)
    {
      Horizontal[i+j*n] = new TJointEqN(CellTree[i+(j-1)*n],
                                        CellTree[i+ j   *n]);
    }

  for(i=0;i<n;i++)
  {
    Horizontal[i] = new TBoundEdge(BdParts[0]->GetBdComp(0), 
                           ((double)i)/n, ((double)(i+1))/n );
    Horizontal[i+m*n] = new TBoundEdge(BdParts[0]->GetBdComp(2), 
                           1-((double)(i+1))/n, 1-((double)i)/n );
  }

  for(i=1;i<n;i++)
    for(j=0;j<m;j++)
    {
      Vertical[i+j*(n+1)] = new TJointEqN(CellTree[i-1+j*n],
                                      CellTree[i  +j*n]);
    }

  for(j=0;j<m;j++)
  {
    Vertical[n+j*(n+1)] = new TBoundEdge(BdParts[0]->GetBdComp(1), 
                           ((double)j)/m, ((double)(j+1))/m );
    Vertical[j*(n+1)] = new TBoundEdge(BdParts[0]->GetBdComp(3), 
                           1-((double)(j+1))/m, 1-((double)j)/m );
  }

  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
    {
      MCell = (TMacroCell *)CellTree[i+j*n];
      MCell->SetVertex(0, Vertices[i  + j   *(n+1)]);
      MCell->SetVertex(1, Vertices[i+1+ j   *(n+1)]);
      MCell->SetVertex(2, Vertices[i+1+(j+1)*(n+1)]);
      MCell->SetVertex(3, Vertices[i  +(j+1)*(n+1)]);

      MCell->SetJoint(0, Horizontal[i+j*n]);
      MCell->SetJoint(1, Vertical[i+1+j*(n+1)]);
      MCell->SetJoint(2, Horizontal[i+(j+1)*n]);
      MCell->SetJoint(3, Vertical[i+j*(n+1)]);

      if(j<m/2)
        MCell->SetSubGridID(0);
      else
        MCell->SetSubGridID(1);
    }
}

void TDomain::initializeDefaultUnitSquareBdry()
{
  Output::root_info("Domain","Loading default domain description: Default_UnitSquare.");

  //create an i/o string stream
  std::stringstream bdryStream;

  //...and fill it with input in .PRM-format
  bdryStream << "NBCT \n";
  bdryStream << "  1 \n";
  bdryStream << "IBCT \n";
  bdryStream << "1 \n";
  bdryStream << "NCOMP \n";
  bdryStream << "4 \n";
  bdryStream << "ITYP NSPLINE NPAR \n";
  bdryStream << "1     1     2 \n";
  bdryStream << "1     1     2 \n";
  bdryStream << "1     1     2 \n";
  bdryStream << "1     1     2 \n";
  bdryStream << "PARAMETERS \n";
  bdryStream << "0.000000E+0000    0.000000E+0000 \n";
  bdryStream << "1.000000E+0000    0.000000E+0000 \n";
  bdryStream << "1.000000E+0000    0.000000E+0000 \n";
  bdryStream << " 0.000000E+0000    1.000000E+0000 \n";
  bdryStream << "1.000000E+0000    1.000000E+0000 \n";
  bdryStream << "-1.000000E+0000    0.000000E+0000 \n";
  bdryStream << "0.000000E+0000    1.000000E+0000 \n";
  bdryStream << "0.000000E+0000   -1.000000E+0000 \n";

  //call the readin method
  ReadBdParam( bdryStream );

}

void TDomain::initializeDrivenCavitySquareBdry()
{
  Output::root_info("Domain","Loading default domain description: DrivenCavitySquare.");

  //create an i/o string stream
  std::stringstream bdryStream;

  //...and fill it with input in .PRM-format
  bdryStream << "NBCT \n";
  bdryStream << "  1 \n";
  bdryStream << "IBCT \n";
  bdryStream << "1 \n";
  bdryStream << "NCOMP \n";
  bdryStream << "4 \n";
  bdryStream << "ITYP NSPLINE NPAR \n";
  bdryStream << "1     1     2 \n";
  bdryStream << "1     1     2 \n";
  bdryStream << "1     1     2 \n";
  bdryStream << "1     1     2 \n";
  bdryStream << "PARAMETERS \n";
  bdryStream << "-1.000000E+0000    -1.000000E+0000 \n";
  bdryStream << "2.000000E+0000    0.000000E+0000 \n";
  bdryStream << "1.000000E+0000    -1.000000E+0000 \n";
  bdryStream << " 0.000000E+0000    2.000000E+0000 \n";
  bdryStream << "1.000000E+0000    1.000000E+0000 \n";
  bdryStream << "-2.000000E+0000    0.000000E+0000 \n";
  bdryStream << "-1.000000E+0000    1.000000E+0000 \n";
  bdryStream << "0.000000E+0000   -2.000000E+0000 \n";

  //call the readin method
  ReadBdParam( bdryStream );

}

#else

void TDomain::SetBoundBox(double boundx, double boundy, double boundz)
{
  BoundX = boundx;
  BoundY = boundy;
  BoundZ = boundz;
}

void TDomain::TestGrid3D()
{
//  TVertex *v[9];
//  TBoundFace *b[9];
//  TJointEqN *j1;

//  v[0] = new TVertex(0.0, 0.0, 0.0);
//  v[1] = new TVertex(1.0, 0.0, 0.0);
//  v[2] = new TVertex(1.0, 1.0, 0.0);
//  v[3] = new TVertex(0.0, 1.0, 0.0);
//  v[4] = new TVertex(0.0, 0.0, 1.0);
//  v[5] = new TVertex(1.0, 0.0, 1.0);
//  v[6] = new TVertex(1.0, 1.0, 1.0);
//  v[7] = new TVertex(0.0, 1.0, 1.0);
//  v[8] = new TVertex(1.5, 0.5, 0.5);
         
}

int TDomain::initializeDefaultCubeBdry()
{
  Output::root_info("Domain","Loading default domain description: Default_UnitCube.");

  //create an i/o string stream
  std::stringstream bdryStream;

  //...and fill it with input in .PRM-format
  bdryStream << "NBCT \n";
  bdryStream << "1 \n";
  bdryStream << "IBCT \n";
  bdryStream << "1 \n";
  bdryStream << "NCOMP \n";
  bdryStream << "6 \n";
  bdryStream << "ITYP NSPLINE NPAR \n";
  bdryStream << "10     0     3 \n";
  bdryStream << "10     0     3 \n";
  bdryStream << "10     0     3 \n";
  bdryStream << "10     0     3 \n";
  bdryStream << "10     0     3 \n";
  bdryStream << "10     0     3 \n";
  bdryStream << "PARAMETERS \n";
  bdryStream << "0.000000E+0000    0.000000E+0000    0.000000E+0000 \n";
  bdryStream << "1.000000E+0000    0.000000E+0000    0.000000E+0000 \n";
  bdryStream << "0.000000E+0000    0.000000E+0000   -1.000000E+0000 \n";
  bdryStream << "0.000000E+0000    0.000000E+0000    0.000000E+0000 \n";
  bdryStream << "0.000000E+0000    1.000000E+0000    0.000000E+0000 \n";
  bdryStream << "-1.000000E+0000    0.000000E+0000    0.000000E+0000 \n";
  bdryStream << "0.000000E+0000    0.000000E+0000    0.000000E+0000 \n";
  bdryStream << "1.000000E+0000    0.000000E+0000    0.000000E+0000 \n";
  bdryStream << "0.000000E+0000   -1.000000E+0000    0.000000E+0000 \n";
  bdryStream << "0.000000E+0000    0.000000E+0000    1.000000E+0000 \n";
  bdryStream << "1.000000E+0000    0.000000E+0000    0.000000E+0000 \n";
  bdryStream << "0.000000E+0000    0.000000E+0000    1.000000E+0000 \n";
  bdryStream << "1.000000E+0000    0.000000E+0000    0.000000E+0000 \n";
  bdryStream << "0.000000E+0000    1.000000E+0000    0.000000E+0000 \n";
  bdryStream << "1.000000E+0000    0.000000E+0000    0.000000E+0000 \n";
  bdryStream << "1.000000E+0000    1.000000E+0000    1.000000E+0000 \n";
  bdryStream << "-1.000000E+0000    0.000000E+0000    0.000000E+0000 \n";
  bdryStream << "0.000000E+0000    1.000000E+0000    0.000000E+0000 \n";

  //call the readin method
  bool isSandwich;
  ReadBdParam( bdryStream, isSandwich );

  return isSandwich;
}

int TDomain::initialize_block_1x2x3_bdry()
{
  Output::root_info("Domain","Loading domain description: Default_Block_1x2x3.");

  std::stringstream bdryStream;

  // input in .PRM-format
  bdryStream << "NBCT \n";
  bdryStream << "1 \n";
  bdryStream << "IBCT \n";
  bdryStream << "1 \n";
  bdryStream << "NCOMP \n";
  bdryStream << "6 \n";
  bdryStream << "ITYP NSPLINE NPAR \n";
  bdryStream << "10     0     3 \n";
  bdryStream << "10     0     3 \n";
  bdryStream << "10     0     3 \n";
  bdryStream << "10     0     3 \n";
  bdryStream << "10     0     3 \n";
  bdryStream << "10     0     3 \n";
  bdryStream << "PARAMETERS \n";
  bdryStream << "0.000000E+0000    0.000000E+0000    0.000000E+0000 \n";
  bdryStream << "1.000000E+0000    0.000000E+0000    0.000000E+0000 \n";
  bdryStream << "0.000000E+0000    0.000000E+0000   -1.000000E+0000 \n";
  bdryStream << "0.000000E+0000    0.000000E+0000    0.000000E+0000 \n";
  bdryStream << "0.000000E+0000    2.000000E+0000    0.000000E+0000 \n";
  bdryStream << "-1.000000E+0000    0.000000E+0000    0.000000E+0000 \n";
  bdryStream << "0.000000E+0000    0.000000E+0000    0.000000E+0000 \n";
  bdryStream << "1.000000E+0000    0.000000E+0000    0.000000E+0000 \n";
  bdryStream << "0.000000E+0000   -1.000000E+0000    0.000000E+0000 \n";
  bdryStream << "0.000000E+0000    0.000000E+0000    3.000000E+0000 \n";
  bdryStream << "1.000000E+0000    0.000000E+0000    0.000000E+0000 \n";
  bdryStream << "0.000000E+0000    0.000000E+0000    1.000000E+0000 \n";
  bdryStream << "1.000000E+0000    0.000000E+0000    0.000000E+0000 \n";
  bdryStream << "0.000000E+0000    2.000000E+0000    0.000000E+0000 \n";
  bdryStream << "1.000000E+0000    0.000000E+0000    0.000000E+0000 \n";
  bdryStream << "1.000000E+0000    2.000000E+0000    0.000000E+0000 \n";
  bdryStream << "-1.000000E+0000    0.000000E+0000    0.000000E+0000 \n";
  bdryStream << "0.000000E+0000    1.000000E+0000    0.000000E+0000 \n";

  bool isSandwich;
  ReadBdParam( bdryStream, isSandwich );

  return isSandwich;
}

void TDomain::initialize_cube_hexa_mesh()
{
  Output::root_info("Domain","Loading default initial mesh: Default_UnitCube_Hexa.");

  //create an i/o string stream
  std::stringstream geoStream;

  //...and fill it with input in .PRM-format

  geoStream << "Grobgitter 3D \n";
  geoStream << "Parametrisierung \n";
  geoStream << "1    8    4    8    6         NEL NVT NVpF NVpEL NBF \n";
  geoStream << "DCORVG (comp, t, s  or  x, y, z) \n";
  geoStream << "0.000000E+0000     0.000000E+0000     0.000000E+0000 \n";
  geoStream << "0.000000E+0000     1.000000E+0000     0.000000E+0000 \n";
  geoStream << "0.000000E+0000     1.000000E+0000     1.000000E+0000 \n";
  geoStream << "0.000000E+0000     0.000000E+0000     1.000000E+0000 \n";
  geoStream << "3.000000E+0000     0.000000E+0000     0.000000E+0000 \n";
  geoStream << "3.000000E+0000     1.000000E+0000     0.000000E+0000 \n";
  geoStream << "3.000000E+0000     1.000000E+0000    -1.000000E+0000 \n";
  geoStream << "3.000000E+0000     0.000000E+0000    -1.000000E+0000 \n";
  geoStream << "KVERT \n";
  geoStream << "1   2   3   4   5   6   7   8 \n";
  geoStream << "KNPR \n";
  geoStream << "1   1   1   1   1   1   1   1 \n";
  geoStream << "BoundFaces \n";
  geoStream << "1   2   3   4 \n";
  geoStream << "1   2   6   5 \n";
  geoStream << "2   3   7   6 \n";
  geoStream << "3   4   8   7 \n";
  geoStream << "4   1   5   8 \n";
  geoStream << "5   6   7   8 \n";
  geoStream << "FaceParam (cell, locface, boundary part/component) \n";
  geoStream << "1   1   1   1 \n";
  geoStream << "1   2   1   3 \n";
  geoStream << "1   3   1   5 \n";
  geoStream << "1   4   1   6 \n";
  geoStream << "1   5   1   2 \n";
  geoStream << "1   6   1   4 \n";

  //call the readin method
  ReadGeo( geoStream, false );
}

void TDomain::initialize_cube_tetra_mesh()
{
  Output::root_info("Domain","Loading default initial mesh: Default_UnitCube_Tetra.");

  //create an i/o string stream
  std::stringstream geoStream;

  //...and fill it with input in .PRM-format
  geoStream << "Grobgitter 3D \n";
  geoStream << "Parametrisierung \n";
  geoStream << "6    8    3    4    12         NEL NVT NVpF NVpEL NBF \n";
  geoStream << "DCORVG (comp, t, s  or  x, y, z) \n";
  geoStream << "0.000000E+0000     0.000000E+0000     0.000000E+0000 \n";
  geoStream << "0.000000E+0000     1.000000E+0000     0.000000E+0000 \n";
  geoStream << "0.000000E+0000     1.000000E+0000     1.000000E+0000 \n";
  geoStream << "0.000000E+0000     0.000000E+0000     1.000000E+0000 \n";
  geoStream << "3.000000E+0000     0.000000E+0000     0.000000E+0000 \n";
  geoStream << "3.000000E+0000     1.000000E+0000     0.000000E+0000 \n";
  geoStream << "3.000000E+0000     1.000000E+0000    -1.000000E+0000 \n";
  geoStream << "3.000000E+0000     0.000000E+0000    -1.000000E+0000 \n";
  geoStream << "KVERT \n";
  geoStream << "1  2  3  7 \n";
  geoStream << "1  6  2  7 \n";
  geoStream << "1  3  4  7 \n";
  geoStream << "1  4  8  7 \n";
  geoStream << "1  5  6  7 \n";
  geoStream << "1  8  5  7 \n";
  geoStream << "KNPR  \n";
  geoStream << "1   1   1   1   1   1   1   1   \n";
  geoStream << "BoundFaces \n";
  geoStream << "1   2   3 \n";
  geoStream << "3   2   7 \n";
  geoStream << "1   6   2 \n";
  geoStream << "2   6   7 \n";
  geoStream << "1   3   4 \n";
  geoStream << "4   3   7 \n";
  geoStream << "1   4   8 \n";
  geoStream << "8   4   7 \n";
  geoStream << "1   5   6 \n";
  geoStream << "6   5   7 \n";
  geoStream << "1   8   5 \n";
  geoStream << "5   8   7 \n";
  geoStream << "FaceParam (cell, locface, boundary part/component) \n";
  geoStream << "1   1   1   1 \n";
  geoStream << "1   3   1   5 \n";
  geoStream << "2   1   1   3 \n";
  geoStream << "2   3   1   5 \n";
  geoStream << "3   1   1   1 \n";
  geoStream << "3   3   1   6 \n";
  geoStream << "4   1   1   2 \n";
  geoStream << "4   3   1   6 \n";
  geoStream << "5   1   1   3 \n";
  geoStream << "5   3   1   4 \n";
  geoStream << "6   1   1   2 \n";
  geoStream << "6   3   1   4 \n";

  //call the readin method
  ReadGeo( geoStream, false );
}

void TDomain::initialize_block_1x2x3_tetra_mesh()
{
  Output::root_info("Domain","Loading initial mesh: Default_Block_1x2x3_Tetra.");

  std::stringstream geoStream;
  // Brief description by Ondrej Partl
  // NEL = number of elements
  // NVT = number of nodes
  // NVpF = number of nodes per face of each element (I hope)
  // NVpEL = number of nodes per element
  // NBF = number of element faces on boundary
  // DCORVG: for KNPR == 0, list of coordinates (x, y, z) of nodes
  // KVERT: list of element nodes; first node is 1
  // KNPR: see DCORVG
  // BoundFaces: list of nodes of element faces lying on boundary
  // FaceParam: information about element faces lying on boundary
  
  geoStream << "Coarse grid 3D \n";
  geoStream << "Parametrization \n";
  geoStream << "6    8    3    4    12         NEL NVT NVpF NVpEL NBF \n";
  geoStream << "DCORVG (comp, t, s  or  x, y, z) \n";
  geoStream << "0.000000E+0000     0.000000E+0000     0.000000E+0000 \n";
  geoStream << "0.000000E+0000     2.000000E+0000     0.000000E+0000 \n";
  geoStream << "1.000000E+0000     0.000000E+0000     0.000000E+0000 \n";
  geoStream << "1.000000E+0000     2.000000E+0000     0.000000E+0000 \n";
  geoStream << "1.000000E+0000     0.000000E+0000     3.000000E+0000 \n";
  geoStream << "0.000000E+0000     0.000000E+0000     3.000000E+0000 \n";
  geoStream << "1.000000E+0000     2.000000E+0000     3.000000E+0000 \n";
  geoStream << "0.000000E+0000     2.000000E+0000     3.000000E+0000 \n";
  geoStream << "KVERT \n";
  geoStream << "1  3  4  7 \n";
  geoStream << "1  5  3  7 \n";
  geoStream << "1  4  2  7 \n";
  geoStream << "1  2  8  7 \n";
  geoStream << "1  6  5  7 \n";
  geoStream << "1  8  6  7 \n";
  geoStream << "KNPR  \n";
  geoStream << "0   0   0   0   0   0   0   0   \n";
  geoStream << "BoundFaces \n";
  geoStream << "1   3   4 \n";
  geoStream << "3   4   7 \n";
  geoStream << "1   5   3 \n";
  geoStream << "3   7   5 \n";
  geoStream << "1   4   2 \n";
  geoStream << "4   2   7 \n";
  geoStream << "1   2   8 \n";
  geoStream << "2   8   7 \n";
  geoStream << "1   6   5 \n";
  geoStream << "6   5   7 \n";
  geoStream << "1   8   6 \n";
  geoStream << "8   6   7 \n";
  geoStream << "FaceParam (element, local face, boundary part, boundary face corresponding to .PRM file) \n";
  geoStream << "1   1   1   1 \n";
  geoStream << "1   3   1   5 \n";
  geoStream << "2   1   1   3 \n";
  geoStream << "2   3   1   5 \n";
  geoStream << "3   1   1   1 \n";
  geoStream << "3   3   1   6 \n";
  geoStream << "4   1   1   2 \n";
  geoStream << "4   3   1   6 \n";
  geoStream << "5   1   1   3 \n";
  geoStream << "5   3   1   4 \n";
  geoStream << "6   1   1   2 \n";
  geoStream << "6   3   1   4 \n";

  ReadGeo( geoStream, false );
}

void TDomain::initialize_block_1x2x3_hexa_mesh()
{
  Output::root_info("Domain","Loading initial mesh: Default_Block_1x2x3_Hexa.");

  std::stringstream geoStream;
  // Brief description by Ondrej Partl
  // NEL = number of elements
  // NVT = number of nodes
  // NVpF = number of nodes per face of each element (I hope)
  // NVpEL = number of nodes per element
  // NBF = number of element faces on boundary
  // DCORVG: for KNPR == 0, list of coordinates (x, y, z) of nodes
  // KVERT: list of element nodes; first node is 1
  // KNPR: see DCORVG
  // BoundFaces: list of nodes of element faces lying on boundary
  // FaceParam: information about element faces lying on boundary
  
  geoStream << "Coarse grid 3D \n";
  geoStream << "Parametrization \n";
  geoStream << "1    8    4    8    6         NEL NVT NVpF NVpEL NBF \n";
  geoStream << "DCORVG (comp, t, s  or  x, y, z) \n";
  geoStream << "0.000000E+0000     0.000000E+0000     0.000000E+0000 \n";
  geoStream << "1.000000E+0000     0.000000E+0000     0.000000E+0000 \n";
  geoStream << "1.000000E+0000     2.000000E+0000     0.000000E+0000 \n";
  geoStream << "0.000000E+0000     2.000000E+0000     0.000000E+0000 \n";
  geoStream << "0.000000E+0000     0.000000E+0000     3.000000E+0000 \n";
  geoStream << "1.000000E+0000     0.000000E+0000     3.000000E+0000 \n";
  geoStream << "1.000000E+0000     2.000000E+0000     3.000000E+0000 \n";
  geoStream << "0.000000E+0000     2.000000E+0000     3.000000E+0000 \n";
  geoStream << "KVERT \n";
  geoStream << "1   2   3   4   5   6   7   8 \n";
  geoStream << "KNPR \n";
  geoStream << "0   0   0   0   0   0   0   0   \n";
  geoStream << "BoundFaces \n";
  geoStream << "1   2   3   4 \n";
  geoStream << "1   2   6   5 \n";
  geoStream << "2   3   7   6 \n";
  geoStream << "3   4   8   7 \n";
  geoStream << "4   1   5   8 \n";
  geoStream << "5   6   7   8 \n";
  geoStream << "FaceParam (element, local face, boundary part, boundary face corresponding to .PRM file) \n";
  geoStream << "1   1   1   1 \n";
  geoStream << "1   2   1   3 \n";
  geoStream << "1   3   1   5 \n";
  geoStream << "1   4   1   6 \n";
  geoStream << "1   5   1   2 \n";
  geoStream << "1   6   1   4 \n";

  ReadGeo( geoStream, false );
}

#endif
