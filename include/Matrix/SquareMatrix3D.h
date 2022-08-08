#ifndef __SQUAREMATRIX3D__
#define __SQUAREMATRIX3D__

#include <SquareMatrix.h>

/** @brief deprecated, use the FEMatrix class instead */
class TSquareMatrix3D : public TSquareMatrix
{
  public:
    /** generate the matrix */
    explicit TSquareMatrix3D(std::shared_ptr<const TFESpace3D> space);
    
    TSquareMatrix3D(const TSquareMatrix3D & m) = default;
    
    /** destructor: free Entries array */
    ~TSquareMatrix3D() = default;
};

#endif // __SQUAREMATRIX3D__
