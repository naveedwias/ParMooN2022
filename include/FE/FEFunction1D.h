#ifndef __FEFUNCTION1D__
#define __FEFUNCTION1D__

#include "FESpace1D.h"
#include "memory"

/** a function from a finite element space */
class TFEFunction1D
{
  protected:
    /** name of the function */
    std::string Name;

    /** space to which this function belongs to */
    std::shared_ptr<TFESpace1D> FESpace1D;

    /** double vector according to FE isomorphism */
    double *Values;

    /** length of vector */
    int Length;

  public:
    /** constructor with vector initialization */
    TFEFunction1D(std::shared_ptr<TFESpace1D> fespace1D,
                  const std::string& name, double *values, int length);

    /** destructor */
    ~TFEFunction1D();

    /** return name */
    std::string GetName()
    { return Name; }

    /** return fe space */
    std::shared_ptr<TFESpace1D> GetFESpace1D()
    { return FESpace1D; }

    /** return length */
    int GetLength()
    { return Length; }

    /** return vector of data */
    double *GetValues()
    { return Values; }
};

#endif
