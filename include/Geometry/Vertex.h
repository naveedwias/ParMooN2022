// =======================================================================
// @(#)Vertex.h        1.1 10/30/98
// 
// Class:       TVertex
// Purpose:     a vertex in a grid
//
// Author:      Volker Behns  09.07.97
//              Sashikumaar Ganesan 05.11.09 (added parallel methods)
//              Sashikumaar Ganesan 08.09.2010 (added 3D parallel methods)
// =======================================================================

#ifndef __VERTEX__
#define __VERTEX__

#include <MooNMD_Io.h>
class TBaseCell;

/** a vertex in a grid */
class TVertex
{
  protected:
    /** first coordinate */
    double X;
    /** second coordinate */
    double Y;
    /** third coordinate (3D) */
    double Z = 0;

    /** an integer for storing clipboard information, Bad style*/
    mutable int ClipBoard;

#ifdef _MPI

    /** Number of 3D cells containing this cells **/
    /** Note !this info only set for dependent cells !!!!!!! */
    int N_Cells;

    /** cells */
    /** Note ! this info only set for dependent cells !!!!!!!!!!*/
    TBaseCell **Cells;

    /** marking this vertex as subdomain vertex */
    bool SubDomainVert;

    /** marking this vertex as cross vertex */
    bool CrossVert;

    /** an integer which stores the number of ranks (SubDomains) contain this vertex */
    int N_SubDomains;

    /** an integer which stores the rank of SubDomains, which contain this vertex */
    int *SubDomain_Ranks;

    /** list of neib cell Globalnumbers, which incident this vertex */
    int *SubDomainGlobalCellNo;

    /** list of neib cell local vert no */
    int *SubDomainLocVertNo;

    /** an integer which stores the number of Cross neib cells, which incident this vertex */
    int N_CrossNeibCells;

    /**  identifier of the physical part */
    int physical_id;
#endif


    /** marking this vertex as Bound vertex */
    bool BoundVert;

    /** if the value is >-1: marking this vertex as Periodic vertex
     *  the same value is shared by the vertices to be identified,
     *  a same value is shared by 2^i vertices (in case of many periodic
     *  boundaries, some vertices belong to i different periodic boundaries */
    int PeriodicVertIndex = -1;

  public:
    // Constructors

      /** 3D vertex */
      TVertex(double initX, double initY, double initZ);
      /** 2D vertex */ 
      TVertex(double initX, double initY);

    // Destructor
    ~TVertex();
    TVertex& operator=(const TVertex&) = delete;
    TVertex(const TVertex&) = delete;

    // Methods

    // set coordinates
    /** set the coordinates in 3D */
    void SetCoords(double initX, double initY, double initZ);
    /** set the coordinate in 2D */
    void SetCoords(double initX, double initY);

    /** return the x coordinate */
    double GetX() const
    { return X; }
    /** return the y coordinate */
    double GetY() const
    { return Y; }

    /** return the z coordinate (3D) */
    double GetZ() const
    { return Z; }
    /** return all three coordinates */
    void GetCoords(double& x, double& y, double& z) const
    {
      x = X;
      y = Y;
      z = Z; // 0 in 2D
    }
    /** return all two coordinates */
    void GetCoords(double& x, double& y) const
    {
      x = X;
      y = Y;
    }

    /** write some information of the vertex in stream s */
    friend std::ostream& operator << (std::ostream& s, const TVertex *v);
    /**
     * This operator introduces an alphanumeric order on the vertices. It will
     * compare first the x component, then the y and finally the z component, if the
     * other ones are equal.
     * The vertices are regarded as equal with a tolerance of 1e-8.
     */
    friend bool operator < (const TVertex& V, const TVertex& W);

    /** set value in ClipBoard */
    void SetClipBoard(int value) const
    { ClipBoard=value; }
    /** get value from ClipBoard */
    int GetClipBoard() const
    { return ClipBoard; }

    void SetAsBoundVert()
    { BoundVert = true; }

    bool IsBoundVert() const
    { return BoundVert; }

    bool IsPeriodicVert() const
    { return PeriodicVertIndex >= 0; }

    void SetPeriodicVertIndex(int idx)
    { PeriodicVertIndex = idx; }

    int GetPeriodicVertIndex() const
    { return PeriodicVertIndex; }


#ifdef _MPI
    static int NbPeriodicVert;

    /** Note ! this info only set for dependent cells !!!!!! */
    void SetVertexCells(int n_Cells, TBaseCell **cells);

    void SetSubDomainInfo(int n_SubDomains, int *subDomain_Ranks, int *subDomainGlobalCellNo, 
                          int *subDomainLocVertNo);

    void AddCrossNeib(int Neib_ID);

    void SetAsSubDomainVert()
    { SubDomainVert = true; }

    bool IsSubDomainVert()
    { return SubDomainVert; }

    void SetAsCrossVert()
    { CrossVert = true; }

    bool IsCrossVert()
    { return CrossVert; }

    void GetCrossNeibs(int &n_VertCrossNeibs, int *&vertCrossNeibs) 
    {
      n_VertCrossNeibs = N_CrossNeibCells;
      vertCrossNeibs = SubDomain_Ranks;
    }

    void GetCrossNeibsInfo(int &N_NeibCells, int *&NeibCellRank, 
                          int *&GlobalNo, int *&LocVertNo)
    {
      N_NeibCells = N_CrossNeibCells;
      NeibCellRank = SubDomain_Ranks;
      GlobalNo = SubDomainGlobalCellNo;
      LocVertNo = SubDomainLocVertNo;
    }

    void GetNeibs(int &n_Neibs, const TBaseCell * const *&neighbs) const
    {
      n_Neibs = N_Cells;
      neighbs = Cells;
    }

    int GetNNeibs() const
    { return N_Cells; }
    /** return physical id */
    void set_physical_id(int _ref)
    { physical_id = _ref; }

    /** return physical id */
    int get_physical_id() const
    { return physical_id; }
#endif
};


#endif
