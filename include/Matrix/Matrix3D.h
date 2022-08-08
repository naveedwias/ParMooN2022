#ifndef __MATRIX3D__
#define __MATRIX3D__

#include <FEMatrix.h>

/** @brief deprecated, use the FEMatrix class instead */
class TMatrix3D : public FEMatrix
{
  public:
    /** generate the matrix */
    TMatrix3D(std::shared_ptr<const TFESpace3D> testspace,
              std::shared_ptr<const TFESpace3D> ansatzspace);
    
    TMatrix3D(const TMatrix3D & m) = default;
    
    /** destructor: free Entries array */
    ~TMatrix3D() = default;
};

#endif // __MATRIX3D__
