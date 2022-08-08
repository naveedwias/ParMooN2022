#include <RefTrans3D.h>
#include <MooNMD_Io.h>

void TRefTrans3D::GetTransformationDerivatives(double, double, double,
	double*) const
{
	ErrThrow("Reference transformation derivatives not implemented for cell type!");
}

void TRefTrans3D::GetInverseTransformationDerivatives(
  double xi, double eta, double zeta, double* matrix) const
{
  double J[9];

  GetTransformationDerivatives(xi, eta, zeta, J);

  // J is column-major:
  //
  // 0 3 6
  // 1 4 7
  // 2 5 8

  // first row: cross-product of second and third column
  matrix[0] = J[4] * J[8] - J[5] * J[7];
  matrix[3] = J[5] * J[6] - J[3] * J[8];
  matrix[6] = J[3] * J[7] - J[4] * J[6];

  // second row: cross-product of third and first column
  matrix[1] = J[7] * J[2] - J[8] * J[1];
  matrix[4] = J[8] * J[0] - J[6] * J[2];
  matrix[7] = J[6] * J[1] - J[7] * J[0];

  // third row: cross-product of first and second column
  matrix[2] = J[1] * J[5] - J[2] * J[4];
  matrix[5] = J[2] * J[3] - J[0] * J[5];
  matrix[8] = J[0] * J[4] - J[1] * J[3];

  double det = J[0] * matrix[0] + J[1] * matrix[3] + J[2] * matrix[6];

  if (det == 0.0)
  {
    ErrThrow("Degenerate cell found!");
  }

  for (int j = 0; j < 9; j++)
  {
    matrix[j] /= det;
  }
}