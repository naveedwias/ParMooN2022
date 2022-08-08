// =======================================================================
// @(#)JointEqN.h        1.3 11/15/99
// 
// Class:       TJointEqN
// Purpose:     connects two cells
//
// Author:      Volker Behns  23.07.97
//
// History:     Add methods for manipulating a certain neighbour 
//              (Gunar Matthies 17.10.97)
//
// =======================================================================

#ifndef __JOINTEQN__
#define __JOINTEQN__

#include <Joint.h>

/** connects two cells */
class TJointEqN : public TJoint
{
  public:
    // Constructors
    /** constructor with one initial neighbour */
    explicit TJointEqN(TBaseCell *neighb0);
    /** constructor with two initial neighbours */
    TJointEqN(TBaseCell *neighb0, TBaseCell *neighb1);

    // Methods
    /** checking for matching joints */
    virtual int CheckMatchingRef(TBaseCell *Me, int J_i,
                     struct StoreGeom &Tmp);

    /** create a new instance of this class */
    virtual TJoint *NewInst(double, double, TBaseCell *Me)
    {
      return new TJointEqN(Me);
    }
    virtual TJoint *NewInst()
    { return new TJointEqN(nullptr); }

    /** return whether this is an interior joint */
    virtual bool InnerJoint() const
    { return true; }
    
    virtual bool is_isoparametric() const
    { return false; }
    
    //Destructor
    virtual ~TJointEqN();   
};

#endif
