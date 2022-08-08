// =======================================================================
// @(#)BDEdge3D.h  
// 
// Class:       TBDEdge3D
// Purpose:     class for boundary edges in 3D
//
// Author:      Sashikumaar Ganesan  03.09.2010
//
// History:
//
// ======================================================================= 

#ifndef __BDEDGE3D__
#define __BDEDGE3D__

#include <Edge.h>


/** an edge in a 3D grid */
class TBDEdge3D : public TEdge
{
  protected:

  public:

    /** constructor with neighbours */
    TBDEdge3D(int n_Neibs, TBaseCell **neighbs);

};

#endif
