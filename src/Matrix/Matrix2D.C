#include <Matrix2D.h>

TMatrix2D::TMatrix2D(std::shared_ptr<const TFESpace2D> testspace,
                     std::shared_ptr<const TFESpace2D> ansatzspace)
 : FEMatrix(testspace, ansatzspace)
{
  
}
