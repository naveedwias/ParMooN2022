#ifndef __SQUAREMATRIX2D__
#define __SQUAREMATRIX2D__

#include <SquareMatrix.h>

/** @brief deprecated, use the FEMatrix class instead */
class TSquareMatrix2D : public TSquareMatrix
{
  public:
    /** generate the matrix */
    explicit TSquareMatrix2D(std::shared_ptr<const TFESpace2D> space);
    
    TSquareMatrix2D(const TSquareMatrix2D & m) = default;

    /** destructor: free Entries array */
    ~TSquareMatrix2D() = default;
};

#endif
