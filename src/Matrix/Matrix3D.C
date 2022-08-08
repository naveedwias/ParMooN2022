#include <Matrix3D.h>

TMatrix3D::TMatrix3D(std::shared_ptr<const TFESpace3D> testspace,
                     std::shared_ptr<const TFESpace3D> ansatzspace)
 : FEMatrix(testspace, ansatzspace)
{
  
}
