// =======================================================================
// @(#)Joint.h        1.7 11/15/99
// 
// Class:       TJoint
// Purpose:     superclass for edges and faces
//
// Author:      Volker Behns  23.07.97
//
// History:     Add method for getting a certain neighbour
//              (Gunar Matthies 17.10.97)
//
// =======================================================================

#ifndef __JOINT__
#define __JOINT__

#include <Constants.h>

#ifndef __3D__
  #define MAXN_nVpoJ  3
  #define MAXN_nJpoJ  2
#else
  #define MAXN_nVpoJ  9
  #define MAXN_nJpoJ  4
#endif

class TJoint;

struct StoreGeom
{
  const TVertex *Vertices[MAXN_nVpoJ];
  const TJoint *Joints[MAXN_nJpoJ];
  bool Filled;
};

/** supercall for edges and faces */
class TJoint
{
  protected:
    JointType ID;

    /** first neighbour */
    TBaseCell *Neighb0;
    /** second neighbour */
    TBaseCell *Neighb1;

    /** value in ClipBoard (dG method)*/
    int ClipBoard;
    
    /** */
    int NeibSubDomainLocalJointNo;
    
    /** The index of this joint in the two neighbors */
    ///@todo set the size of the array as a function of the Joint type (e.g. 1 for boundaries)
    int IndexInNeighbor[2];
    
    int MapType = 0;

  public:
    // Constructors
    TJoint();

    // Methods
    /** return type */
    JointType GetType() const
    { return ID; }
    
    void ChangeType(JointType New_ID)
     {ID = New_ID;}

#ifdef __3D__
    /** set mapper type automatically */
    void SetMapType();

    /** set mapper type */
    void SetMapType(int maptype)
    { MapType = maptype; }

    /** Function is used to get local edge index on neighboured element */
    int GetNeighbourEdgeIndex(const TBaseCell*, int) const;
#endif

    /** return mapper type */
    int GetMapType() const
    { return MapType; }

    /** check the refinement pattern on both sides for matching,
        return already existing object on the joint in Tmp */
    virtual int CheckMatchingRef(TBaseCell *Me, int J_i,
                  struct StoreGeom &Tmp) = 0;

    /** create a new instance of the same class */
    virtual TJoint *NewInst(double T_0, double T_1, TBaseCell *Me) = 0;
    virtual TJoint *NewInst() = 0;

    /** set the neighbour to Neighb */
    int SetNeighbour(TBaseCell *Neighb);
    /** return the neighbour of this joint which is not equal to Me */
    TBaseCell *GetNeighbour(const TBaseCell *Me) const;

    /** set neighbour i to Neighb */
    int SetNeighbour(int i, TBaseCell *Neighb);
    /** return neighbour with number i */
    TBaseCell *GetNeighbour(int i) const;

    /** remove a neighbour */
    void Delete(TBaseCell *Neighb);

    /** function for debug purpose only */
    TBaseCell *GetNeighb(int i) const
    { return(i ? Neighb1 : Neighb0); }

    /** return whether this is an interior joint */
    virtual bool InnerJoint() const = 0;
    
    /** return whether this is a joint on an isoparametric */
    virtual bool is_isoparametric() const = 0;
    
    /** @brief set the local index of this joint in given neighbor
     * @todo write a more general function that does not take index as input
     */
    void set_index_in_neighbour(const TBaseCell *neigh, int index);
    /** get the index of this joint in given neighbor */
    int get_index_in_neighbour(const TBaseCell*const neigh) const;

    /** @brief get the number that the joint has in a given cell
     *
     * @details If the joint is a facet of the given cell the local index that
     * the joint has in the cell is returned. Else an error is thrown.
     * @param[in] cell a cell for which the joint is a facet
     * @returns the local index the joint has in the cell
     */
    int get_joint_nr_in_cell(const TBaseCell* cell) const;


    /**
     * @brief computes the diameter of the joint.
     * @details This function computes the diameter of a (d-1) dimensional
     * joint. The facet can either be a line or an area. In the former case this
     * function returns the length of the line, and in the latter case it
     * returns the diameter of the area.
     * @note For this function to work Neighb0 has to be set. You might want
     * to call the method like this:
     * \code{.cpp}
     *    TBaseCell cell_for_joint;  // cell where this joint is a part of
     *    auto joint = cell_for_joint(joint_nr) // joint_nr is the local number of the joint of interest
     *    double diam;
     *    if (joint->GetNeighb(0) == nullptr)
     *    {
     *      joint->SetNeighbour(cell_for_joint);  // set Neighb0
     *      diam = joint->GetDiameter();  // determine the diameter
     *      joint->Delete(cell_for_joint); // delete Neighb0 to not infer with later calls of SetNeighbor
     *    }
     *    else
     *    {
     *      diam = joint->GetDiameter();
     *    }
     * \endcode
     * @return the diameter of the (d-1) dimensional joint.
     */
    double GetDiameter();

    #ifdef __3D__
      /** return mapper of refined vertices and faces */
      void GetMapperRef(const int *&MapVerts, const int *&MapFaces) const;

      /** return mapper of original vertices and edges */
      void GetMapperOrig(const int *&MapVerts, const int *&MapEdges) const;
    #endif
            
    /** set value in ClipBoard */
    void SetClipBoard(int value)
    { ClipBoard=value; }
    /** get value from ClipBoard */
    int GetClipBoard()
    { return ClipBoard; }
    
    void SetNeibSubDomainLocalJointNo(int value)
    { NeibSubDomainLocalJointNo = value; }
    
    int GetNeibSubDomainLocalJointNo()
    { return NeibSubDomainLocalJointNo; }  
    
    // Destructor
    virtual ~TJoint();

};

#endif
