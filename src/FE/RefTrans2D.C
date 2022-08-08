#include <RefTrans2D.h>
#include <MooNMD_Io.h>

void TRefTrans2D::GetTransformationDerivatives(double, double, double*) const
{
  ErrThrow("Reference transformation derivatives not implemented for cell type!");
}

void TRefTrans2D::GetInverseTransformationDerivatives(double xi, double eta, double* matrix) const
{
  GetTransformationDerivatives(xi, eta, matrix);

  double a11 = matrix[0];
  double a21 = matrix[1];
  double a12 = matrix[2];
  double a22 = matrix[3];

  double det = a11 * a22 - a12 * a21;

  if (det == 0.0)
  {
    ErrThrow("Degenerate cell found!");
  }

  matrix[0] = a22 / det;
  matrix[1] = -a21 / det;
  matrix[2] = -a12 / det;
  matrix[3] = a11 / det;
}