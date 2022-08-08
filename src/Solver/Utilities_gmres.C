#include <Utilities_gmres.h>
#include <BlockVector.h>
#include <MooNMD_Io.h>
#include <cmath>

HessenbergMatrix::HessenbergMatrix(const unsigned int size) : entries(0)
{
  resize(size);
}

void HessenbergMatrix::resize(const unsigned int size)
{
  if (size != entries.size())
  {
    entries.resize(size);

    if (size > 1)
    {
      // the first two rows have `size` many columns
      entries[0] = std::vector<double>(size, 0.0);
      entries[1] = std::vector<double>(size, 0.0);

      // the following rows have one column less the their previous row, the
      // last row has two entries
      for (unsigned int row = 2; row < size; ++row)
      {
        size_t n = size - row + 1; // number of entries in this row
        entries[row] = std::vector<double>(n, 0.0);
      }
    }
  }
}

double& HessenbergMatrix::operator()(int i, int j)
{
  if (i > j + 1)
  {
    ErrThrow("HessenbergMatrix: index (", i, ",", j, ") is not upper "
             "triangular");
  }

  // notice that index 0 in  entries[i] corresponds to column i - 1
  j = 1 + j - i;
  return entries.at(i).at(j);
}

const double& HessenbergMatrix::operator()(int i, int j) const
{
  if (i > j + 1)
  {
    ErrThrow("HessenbergMatrix: index (", i, ",", j, ") is not upper "
             "triangular");
  }

  // notice that index 0 in  entries[i] corresponds to column i-1
  j = 1 + j - i;
  return entries.at(i).at(j);
}

void GenerateGivensRotation(const double &dx, const double &dy,
    double &cs, double &sn)
{
  if (dy == 0.0)
  {
    cs = 1.0;
    sn = 0.0;
  }
  else if (std::abs(dy) > std::abs(dx))
  {
    double temp = dx / dy;
    sn = 1.0 / std::sqrt( 1.0 + temp * temp);
    cs = temp * sn;
  }
  else
  {
    double temp = dy / dx;
    cs = 1.0 / std::sqrt( 1.0 + temp * temp);
    sn = temp * cs;
  }
}

void ApplyGivensRotation(double &dx, double &dy,
  const double &cs, const double &sn)
{
  double temp = cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx = temp;
}

template <class Matrix, class Vector>
void Update(Vector &x, const int k, const Matrix &h, std::vector<double> y,
            std::vector<Vector>& v)
{
  // Hessenberg backsolve:
  for (int i = k; i >= 0; i--)
  {
    y.at(i) /= h(i, i);
    for (int j = i - 1; j >= 0; j--)
    {
      y[j] -= h(j, i) * y[i];
    }
  }

  for (int j = 0; j <= k; j++)
  {
    x.add_scaled(v[j], y[j]);
  }
}

template void Update<HessenbergMatrix, BlockVector>(BlockVector&, const int,
  const HessenbergMatrix&, std::vector<double>, std::vector<BlockVector>&);