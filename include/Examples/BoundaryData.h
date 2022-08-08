#ifndef __BOUNDARYDATA__
#define __BOUNDARYDATA__

#include <functional> // std::function
#include <utility> // std::pair
#include <Joint.h>
#include <Point.h>
#include <BoundaryCondition.h>
#include <AnalyticalFunction.h>

// forward declaration
enum class Problem_type;

using namespace parmoon;

/** ************************************************************************
*
* @class      BoundaryData
* @brief      handle where on the boundary are which boundary data
*
* There some different methods to find the boundary data at a specific point on
* the boundary. This is always related to the BoundaryComponent of a Mesh.
*
* This class represents scalar data only, however you can use a vector of
* objects of this class to represent components (x/y or normal/tangential).
*
* As a simplification there is a possibility to construct a BoundaryData object
* which is constant, e.g. all-homogeneous or all-one.
*
* Either the from_parametrization*d (*=2 or 3), the from_coordinates, or the
* from_component is used to compute the boundary data. Depending on the
* constructor only one of them is set to a callable function.
*
* In general boundary data can depend on the boundary component, the location
* in space (coordinates or parametrization). and the point in time. This class 
* enables simpler constructors if the resulting boundary data is not supposed to
* depend on some of these.
*
* @author     Ulrich Wilbrandt
* @date       03.11.2014
*
****************************************************************************/
class BoundaryData
{
 public:
  /** @brief constructor for BoundaryData which only depends on the boundary
   * component.
   * 
   * The function 'f' returns boundary data for a given boundary component.
   */
  explicit BoundaryData(std::function<double(unsigned int)> f);
  
  /** @brief constructor in 2D using the parametrization on the boundary
   *
   * The function 'f' returns boundary data for a given boundary component
   * and parameter (often called t). The last argument is the point in time.
   *
   * For stationary BoundaryData simply define `f` such that it does not depend
   * on its last argument.
   */
  explicit BoundaryData(std::function<double(unsigned int, double, double)> f);

  /** @brief constructor in 3D using the parametrization on the boundary
   *
   * The function 'f' returns boundary data for a given boundary component
   * and parameters (often called t and s). The last argument is the point in
   * time.
   *
   * For stationary BoundaryData simply define `f` such that it does not depend
   * on its last argument.
   */
  explicit BoundaryData(
    std::function<double(unsigned int, double, double, double)> f);

  /** @brief constructor using the coordinates on the boundary (2D and 3D)
   *
   * The function 'f' returns boundary data for given coordinates(x, y). The
   * last argument is the point in time.
   *
   * For stationary BoundaryData simply define `f` such that it does not depend
   * on its last argument.
   */
  explicit BoundaryData(std::function<double(const Point&, double)> f);

  /** @brief constructor for constant boundary data (2D and 3D)
   *
   * Use this if your boundary data is constant on the entire boundary, for
   * example homogeneous.
   */
  explicit BoundaryData(double);

  /** @brief compute the boundary data given boundary conditions and an exact
   * solution.
   *
   * Whenever an exact solution is known, this is the easiest way to define
   * matching boundary data.
   */
  BoundaryData(const BoundaryCondition& bc, const AnalyticalFunction& exact,
               Problem_type type);

  /** @brief get the boundary data on a given joint and parameter s (2D)
   *
   * The parameter s is in \f$[-1,1]\f$ and determines where on the joint the
   * boundary data is supposed to be computed. s=-1 means on first vertex, s=1
   * means on second vertex, s=0 means in the center.
   *
   * For time dependent boundary data the third argument is used, otherwise it
   * is ignored.
   */
  double get_data(const TJoint& joint, double s, double time = 0.) const;

  /** @brief get the boundary data on a given joint and parameters t and s (3D)
   *
   * The parameter range for t and s depend on the Mesh and are a priori unknown
   * here. This is why using this method is not recommended. Use
   * get_data(const Point&, double)
   */
  double
  get_data(const TJoint& joint, std::pair<double, double> ts, double time = 0.)
  const;

  /** @brief get the boundary data at given coordinates (2D)
   *
   * Note that is is not (always) checked if the specified point is really on
   * the boundary.
   */
  double get_data(const TJoint& joint, const Point&, double time = 0.) const;
  
  /// @brief indicate if the BoundaryData depends on time. 
  ///
  /// This can avoid reassembling of a right hand side vector
  bool depends_on_time = true;

 private:
  /** @brief function to compute the boundary data on a given boundary component
   * 
   * These boundary data therefore do not depend on time and space.
   */
  std::function<double(unsigned int)> from_component;
   
  /** @brief function to compute boundary data at a given boundary part,
   *         component and parameter (2D)
   *
   * The unsigned int is the boundary component index, the first double value
   * is the parameter (often called t). The last argument is the point in time.
   */
  std::function<double(unsigned int, double, double)> from_parametrization2d;

  /** @brief function to compute boundary data at a given boundary part,
   *         component and parameters (3D)
   *
   * The unsigned int is the component index, the double values are the
   * parameters (often called t, s). The last argument is the point in time.
   */
  std::function<double(unsigned int, double, double, double)>
      from_parametrization3d;

  /** @brief function to compute boundary data at given coordinates (2D and 3D)
   */
  std::function<double(const Point&, double)> from_coordinates;
};

#endif // __BOUNDARYDATA__
