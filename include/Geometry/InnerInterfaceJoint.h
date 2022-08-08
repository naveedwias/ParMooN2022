// =======================================================================
// @(#)InnerInterfaceJoint.h        1.3 01/26/12
// 
// Class:       TInnerInterfaceJoint
// Purpose:     connects two cells of different collection
//              (not treated as a boundary joint)
//
//
// =======================================================================

#ifndef __INNERINTERFACEJOINT__
#define __INNERINTERFACEJOINT__

#include <JointEqN.h>
#include <Vertex.h>

/** connects two cells on an interface */
class TInnerInterfaceJoint : public TJointEqN
{
  protected:
    /** during refinement the two (in 2D) or four (in 3D) children are marked in the following list */
    TInnerInterfaceJoint *children[4];
    
    /** x coordinate of the begin of line */
    double Xstart;
    /** y coordinate of the begin of line */
    double Ystart;
    /** x progress of line */
    double delX;
    /** y progress of line */
    double delY;

    /** dimension of the interface */
    int dim;
    
    //only in 3D
    //x-, y- and z- component of the n_vert vertices
    unsigned int N_vert;
    double *X, *Y, *Z;
    //local face number of the face in neighboring cell
    //since the local face number can be different depending on the neighboring cell
    //always use the same Cell_Id, e.g. for Stokes-Darcy, always the Stokes or the Darcy cell
    //which cell is used can be found in GetInnerInterfaceJoints in auxiliaryFunctions.C
    int local_face;


 public:
  /** global cell number of the neibs' cell, which contains this joint */
  int NeibGlobalCellNo[2];
  int SubDomainsID[2];
  
  // Constructors
  /** constructor with one initial neighbour */
  explicit TInnerInterfaceJoint(TBaseCell *neighb0);
  TInnerInterfaceJoint(TBaseCell *neighb0, TBaseCell *neighb1);
  
  // Methods
    /** create a new instance of the same class */
  virtual TJoint *NewInst(double T_0, double T_1, TBaseCell *Me);
  virtual TJoint *NewInst();
  
  
  virtual bool InnerJoint() const
  { return true; }
  
  /** Remember which joints are children of this joint */
  void SetChild(TInnerInterfaceJoint *child);
  /** return one of the two children of this edge */
  TInnerInterfaceJoint *GetChild(int child) const;
  
  /** Get either one of the two neighbors */
  TBaseCell *GetNeighbour(int i) const;
  
  /** set/get the coordinates of the start point and the vector pointing to the     
    second point */
  void SetParams (double xstart, double ystart, double delx, double dely);
  void GetParams (double &xstart, double &ystart, double &delx, double &dely) const;
  
  /**
   * WARNING points must be ordered such that point i and point i-1 and i+1
   * are an edge of the face
   */
  void SetParams (double *x, double *y, double *z, int n_vert);

  int GetN_Vertices() const {return N_vert;}

  void GetParams (double *x, double *y, double *z) const;

  void SetLocalface(int l_face) {local_face = l_face;}

  int GetLocalface() const  {return local_face;}
  
  /** Compute the length of this edge */
  double GetLength() const;
  
  /** Compute the area of this face
   *  Warning: order of points must be clockwise or counter clockwise*/
  double GetArea() const;
  
  /** the unit normal of this edge */
  void GetNormal(double &nx, double &ny) const;
  
  /** the unit normal of this face */
  void GetNormal(double &nx, double &ny, double &nz) const;

  /** the unit tangential of this edge */
  void GetTangent(double &tx, double &ty) const;
  
  /** the unit tangentials of this face */
  void GetTangent(double &tx, double &ty, double &tz,
                  double &tx2, double &ty2, double &tz2) const;
  
  /** set/get the index of this joint in given neighbor */
  int GetIndexInNeighbor(TBaseCell const * const neigh) const;
  void SetIndexInNeighbor(TBaseCell *neigh, int index);
  
  /**
   * Prints all infos for debugging/testing
   */
  void PrintInfo() const;
  
    // Destructor
  virtual ~TInnerInterfaceJoint(){};
  
};

#endif
