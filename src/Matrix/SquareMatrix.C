#include <SquareMatrix.h>

TSquareMatrix::TSquareMatrix(std::shared_ptr<const TFESpace1D> space)
 : FEMatrix(space)
{
  
}
TSquareMatrix::TSquareMatrix(std::shared_ptr<const TFESpace2D> space)
 : FEMatrix(space)
{
  
}
#ifdef __3D__
TSquareMatrix::TSquareMatrix(std::shared_ptr<const TFESpace3D> space)
 : FEMatrix(space)
{
  
}
#endif // 3D

