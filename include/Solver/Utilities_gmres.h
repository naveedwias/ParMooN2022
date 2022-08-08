#ifndef __UTILITIES_GMRES__
#define __UTILITIES_GMRES__

#include <vector>

/* ************************************************************************** */
/// @brief a simple class to store an upper triangular matrix and one off 
/// diagonal
///
/// Example:
///     ( * * * * * )
///     ( * * * * * )
///     (   * * * * )
///     (     * * * )
///     (       * * )
/// This is needed for gmres.
class HessenbergMatrix
{
  public:

    explicit HessenbergMatrix(const unsigned int size);

    double& operator()(int i, int j);
    const double& operator()(int i, int j) const;
    void resize(const unsigned int size);

  private:
    // one vector for each row
    std::vector<std::vector<double>> entries;
};

void GenerateGivensRotation(const double &dx, const double &dy,
    double &cs, double &sn);

void ApplyGivensRotation(double &dx, double &dy,
  const double &cs, const double &sn);

template <class Matrix, class Vector>
void Update(Vector &x, const int k, const Matrix &h, std::vector<double> y,
            std::vector<Vector>& v);

#endif