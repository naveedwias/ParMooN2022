#ifndef __SQUAREMATRIX1D__
#define __SQUAREMATRIX1D__

#include <SquareMatrix.h>

/** @brief deprecated, use the FEMatrix class instead */
class TSquareMatrix1D : public TSquareMatrix
{
  public:
    /** generate the matrix */
    explicit TSquareMatrix1D(std::shared_ptr<const TFESpace1D> space);

    /** destructor: free Entries array */
    ~TSquareMatrix1D() = default;
};

#endif

