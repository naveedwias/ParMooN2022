#include <limits> // std::numeric_limits<double>::quiet_NaN
#include <cmath>  // std::isnan
#include "Point.h"
#include "MooNMD_Io.h"

using namespace parmoon;
/* ************************************************************************** */
Point::Point(double x)
 : x(x), y(std::numeric_limits<double>::quiet_NaN()), 
   z(std::numeric_limits<double>::quiet_NaN())
{
  
}

/* ************************************************************************** */
Point::Point(double x, double y)
 : x(x), y(y), z(std::numeric_limits<double>::quiet_NaN())
{
  
}

/* ************************************************************************** */
Point::Point(double x, double y, double z)
 : x(x), y(y), z(z)
{
  
}

/* ************************************************************************** */
Point::Point(const std::vector<double>& p)
{
  const unsigned int p_size = p.size();
  if(p_size == 2)
  {
    this->x = p[0];
    this->y = p[1];
    this->z = std::numeric_limits<double>::quiet_NaN();
  }
  else if(p_size == 3)
  {
    this->x = p[0];
    this->y = p[1];
    this->z = p[2];
  }
  else if(p_size == 1)
  {
    this->x = p[0];
    this->y = std::numeric_limits<double>::quiet_NaN();
    this->z = std::numeric_limits<double>::quiet_NaN();
  }
  else
  {
    ErrThrow("unable to create a Point with a vector of size ", p_size);
  }
}

/* ************************************************************************** */
Point::Point(unsigned int dimension)
 : x(0.0), y(0.0), z(0.0)
{
  if(dimension == 0 || dimension > 3)
  {
    ErrThrow("cannot create a Point of dimension ", dimension);
  }
  if(dimension < 3)
  {
    this->z = std::numeric_limits<double>::quiet_NaN();
    if(dimension < 2)
      this->y = std::numeric_limits<double>::quiet_NaN();
  }
}

/* ************************************************************************** */
unsigned int Point::dimension() const
{
  if(!std::isnan(this->z))
    return 3;
  else if(!std::isnan(this->y))
    return 2;
  else
    return 1;
}

/* ************************************************************************** */
bool Point::is_equal(const Point& B, double tol) const
{
  unsigned int dim = this->dimension();
  if(dim != B.dimension())
  {
    return false;
  }
  else
  {
    if(std::abs(this->x - B.x) > tol)
    {
      return false;
    }
    if(dim>1 && (std::abs(this->y - B.y) > tol))
    {
      return false;
    }
    if(dim>2 && (std::abs(this->z - B.z) > tol))
    {
      return false;
    }
    return true;
  }
}

/* ************************************************************************** */
Point Point::operator- (const Point& B)
{
  unsigned int dim = this->dimension();
  if(dim != B.dimension())
  {
    ErrThrow("unable to create a vector from Points of different dimensions");
  }
  Point C(dim);
  C.x = this->x - B.x;
  if(dim>1)
  {
    C.y = this->y - B.y;
    if(dim>2)
    {
      C.z = this->z - B.z;
    }
  }
  return C;
}

/* ************************************************************************** */
Point Point::operator+ (const Point& B)
{
  unsigned int dim = this->dimension();
  if(dim != B.dimension())
  {
    ErrThrow("unable to translate Point, dimensions mismatch");
  }
  Point C(dim);
  C.x = this->x + B.x;
  if(dim>1)
  {
    C.y = this->y + B.y;
    if(dim>2)
    {
      C.z = this->z + B.z;
    }
  }
  return C;
}

/* ************************************************************************** */
Point operator*(double lambda, const Point &p)
{
  unsigned int dim = p.dimension();

  if (dim == 1)
  {
    return Point(lambda * p.x);
  }
  else if (dim == 2)
  {
    return Point(lambda * p.x, lambda * p.y);
  }
  else
  {
    return Point(lambda * p.x, lambda * p.y, lambda * p.z);
  }
}

/* ************************************************************************** */
Point operator*(const Point &p, double lambda)
{
  unsigned int dim = p.dimension();

  if (dim == 1)
  {
    return Point(lambda * p.x);
  }
  else if (dim == 2)
  {
    return Point(lambda * p.x, lambda * p.y);
  }
  else
  {
    return Point(lambda * p.x, lambda * p.y, lambda * p.z);
  }
}

/* ************************************************************************** */
Point::operator std::vector<double>()
{
  const unsigned int p_size = this->dimension();
  if(p_size == 1)
  {
    std::vector<double> vector_cast{this->x};
    return vector_cast;
  }
  else if(p_size == 2)
  {
    std::vector<double> vector_cast{this->x, this->y};
    return vector_cast;
  }
  else if(p_size == 3)
  {
    std::vector<double> vector_cast{this->x, this->y, this->z};
    return vector_cast;
  }
  else
  {
    std::vector<double> vector_cast;
    return vector_cast;
  }
}
