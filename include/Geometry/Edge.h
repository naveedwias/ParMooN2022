// =======================================================================
// @(#)Edge.h  
// 
// Class:       TEdge
// Purpose:     superclass for edges in 3D
//
// Author:      Sashikumaar Ganesan  03.09.2010
//
// History:
//
// ======================================================================= 

#ifndef __EDGE__
#define __EDGE__

class TEdge;

#include <Vertex.h>
#include <Constants.h>
#include <BaseCell.h>

/** an edge in a 3D grid */
class TEdge
{
  protected:
    JointType EdgeID;

   /** Number of 3D cells containing this edge **/
    int N_Neibs;

    /** neighbours */
    TBaseCell **Neighbs;

    /** an integer for storing clipboard information, Bad style! */
    mutable int ClipBoard;

#ifdef  _MPI
  /** a bool to check this edge belongs to two or more SubDomains **/
    bool SubDomainEdge;

  /** a bool to check this SubDomainEdge is a cross edge **/
  /** i.e, this edge in their cell is the only connection between two subdomains */
    bool SubDomainCrossEdge;

   /** Number of 3D SubDomain cells containing this edge **/
    int N_SubDomainNeibs;

    /** SubDomain neighbours' rank */
    int *SubDomainNeibsRank;

    /** global cell number of the neibs' cell, which contains this edge */
    int *CrossNeibsGlobalNo;

   /** Number of 3D SubDomain cells containing this edge as cross egde**/
    int N_CrossNeibs;

    /** SubDomain neighbours' rank */
    int *CrossNeibsRank;

    /** local edge number of this joint in the neib cell */
    int *CrossNeibsLocalEdgeNo;

    /** local edge number of this joint in the neib cell */
    int *CrossEdgeMaptype;

    /**  identifier of the physical part */
    int physical_id;
#endif

  public:
    // Constructors
    TEdge(int n_Neibs, TBaseCell **neighbs);

 // methods
    /** return type */
    JointType GetType() const
    { return EdgeID; }

    /** set value in ClipBoard */
    void SetClipBoard(int value) const
    { ClipBoard=value; }
    /** get value from ClipBoard */
    int GetClipBoard() const
    { return ClipBoard; }

#ifdef  _MPI

  void InitSubDomainInfo(int rank);

  void SetAsNotCrossEdgeFor(int rank, int Neib_ID);

  void SetCrossNeibInfo(TVertex *Vert_a);

  bool IsSubDomainEdge()
   {return SubDomainEdge;}

//   void SetAsSubDomainCrossEdge()
//   {SubDomainCrossEdge = true;}

  bool IsSubDomainCrossEdge()
   {return SubDomainCrossEdge;}

  void GetNeibs(int &n_Neibs, const TBaseCell * const *&neighbs) const
   {
    n_Neibs = N_Neibs;
    neighbs = Neighbs;
   }

  int GetN_CrossNeibs()
   {return N_CrossNeibs; }

  void GetNeibSubDomainRanks(int &n_SubDomainNeibs, int *&subDomainNeibsRank)
   {
     n_SubDomainNeibs = N_SubDomainNeibs;
     subDomainNeibsRank = SubDomainNeibsRank;
   }

  void GetCrossEdgeNeibs(int &n_CrossNeibs, int *&crossNeibsRank)
   {
    n_CrossNeibs = N_CrossNeibs;
    crossNeibsRank = CrossNeibsRank;
   }

  void GetCrossEdgeNeibs(int &n_CrossNeibs, int *&crossNeibsRank, int *&crossNeibsGlobalNo,
                         int *&crossNeibsLocalEdgeNo, int *&crossEdgeMaptype);

  /** return physical id */
  void set_physical_id(int _ref)
  { physical_id = _ref; }

  /** return physical id */
  int get_physical_id() const
  { return physical_id; }
#endif

  // Destructor
  ~TEdge();

};

#endif
