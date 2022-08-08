/**
 * Little unit test of wrapper class for LAPACK LU factorization and solve.
 * Puts up a system found here
 *   https://www.nag.co.uk/numeric/fl/manual/pdf/F07/f07adf.pdf,
 * computes the LU factorization and solves. Solution is checked against
 * hard-coded known (MATLAB confirmed...) solution.
 *
 * @date 2016/06/21
 * @author Clemens Bartsch
 */
#include <DenseMatrix.h>
#include <MooNMD_Io.h>

#include <vector>
#include <cmath>

int main(int, char**)
{

  DenseMatrix mat(4,4); //construct a zero-filled 4-by-4 matrix

  //set entries
  mat.setEntry(0, 0, 1.8); //1st row
  mat.setEntry(0, 1, 2.88);
  mat.setEntry(0, 2, 2.05);
  mat.setEntry(0, 3, -0.89);

  mat.setEntry(1, 0, 5.25); //2nd row
  mat.setEntry(1, 1, -2.95);
  mat.setEntry(1, 2, -0.95);
  mat.setEntry(1, 3, -3.80);

  mat.setEntry(2, 0, 1.58); //3rd row
  mat.setEntry(2, 1, -2.69);
  mat.setEntry(2, 2, -2.90);
  mat.setEntry(2, 3, -1.04);

  mat.setEntry(3, 0, -1.11); //4th row
  mat.setEntry(3, 1, -0.66);
  mat.setEntry(3, 2, -0.59);
  mat.setEntry(3, 3, 0.80);

  mat.decomposeLU();

  double vec[4] = {1, 2, 3, 4};

  mat.solve(vec);

  double sol_known[4] = {22.4384, -6.62493, 5.01208 , 34.3641};
  for(int i =0 ; i <4; ++i)
  {
    if(std::abs(vec[i] - sol_known[i]) > 1e-4) //tolerance high due to coarse output...
      ErrThrow("No-no!");
  }
}
