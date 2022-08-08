#ifndef __MATRIX2D__
#define __MATRIX2D__

#include <FEMatrix.h>

/** @brief deprecated, use the FEMatrix class instead */
class TMatrix2D : public FEMatrix
{
  public:
    /** generate the matrix */
    TMatrix2D(std::shared_ptr<const TFESpace2D> testspace,
              std::shared_ptr<const TFESpace2D> ansatzspace);
    
    TMatrix2D(const TMatrix2D & m) = default;
    
    /** destructor: free Entries array */
    ~TMatrix2D() = default;
};


#endif // __MATRIX2D__
