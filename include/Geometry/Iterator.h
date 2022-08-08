#ifndef __ITERATOR__
#define __ITERATOR__

class TIterator;

constexpr int MAX_ItLevel = 50;

#include "Enumerations_geometry.h"
#include <Domain.h>

/** iterator to produce a series of cells with some
    special properties */
class TIterator
{
  protected:
    /** current level */
    int Level;

    /** domain on which the iterator works */
    TDomain *Domain;
    /** a copy of tree of cells */
    TBaseCell **CellTree;
    /** number of cells on trees root */
    int N_RootCells;

  public:
    // Constructors
    virtual ~TIterator(){};

    // Methods
    /** set all parameters to the given values */
    int SetParam(TDomain *domain);

    /** return the next cell */
    virtual TBaseCell *Next(int &info) = 0;
    /** return the previous cell */
    virtual TBaseCell *Prev() = 0;

    /** Initialize at level "level" */
    virtual int Init(int level) = 0;

    /** return the maximum level */
    virtual int GetMaxLevel() = 0;
};

#endif
