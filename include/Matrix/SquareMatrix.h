#ifndef __SQUAREMATRIX__
#define __SQUAREMATRIX__

#include <FEMatrix.h>

/** @brief deprecated, use the FEMatrix class instead */
class TSquareMatrix : public FEMatrix
{
  protected:
    /** generate the matrix, called from derived classes */
    explicit TSquareMatrix(std::shared_ptr<const TFESpace1D> space);
    explicit TSquareMatrix(std::shared_ptr<const TFESpace2D> space);
    #ifdef __3D__
    explicit TSquareMatrix(std::shared_ptr<const TFESpace3D> space);
    #endif // 3D
  public:

    /** destructor: free Entries array */
    ~TSquareMatrix() = default;
    
    TSquareMatrix(const TSquareMatrix & m) = default;
};

#endif // __SQUAREMATRIX__
