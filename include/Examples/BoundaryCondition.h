#ifndef __BOUNDARYCONDITION__
#define __BOUNDARYCONDITION__

#include <functional> // std::function
#include <Constants.h>
#include <Joint.h>
#include <Point.h>

using namespace parmoon;

class BoundaryData;

/** ************************************************************************
*
* @class      BoundaryCondition
* @brief      handle where on the boundary are which boundary conditions
*
* There some different methods to find the boundary condition at a specific
* point on the boundary. The is always related to the boundary component of a
* Mesh. Preferably on each boundary component the example poses a different
* boundary condition. However it is also possible to have different boundary
* conditions on one boundary component. However note that on each joint only
* one boundary condition can be prescribed which is evaluated at the center of
* the joint. In that case it is possible to use the parametrization or
* coordinates on the boundary.
*
* This class represents scalar boundary conditions only, however you can use a
* vector of objects of this class to represent components (x/y or
* normal/tangential). Furthermore it is not possible to define boundary
* conditions which depend on time, however boundary data can.
*
* As a simplification there is a possibility to construct a BoundaryCondition
* which is constant, e.g. all-Dirichlet or all-Neumann.
*
*
* @author     Ulrich Wilbrandt
* @date       03.11.2014
*
****************************************************************************/
class BoundaryCondition
{
 public:
  /** @brief recommended constructor in 2D and 3D
   *
   * This sets the 'from_component' member variable. Whenever the boundary
   * condition only depends on the boundary component this is the preferred
   * constructor. If on the other hand if the boundary conditions differ within
   * one boundary component, this constructor is not suitable.
   *
   * @param f function returning the boundary condition at a given component
   */
  BoundaryCondition(std::function<BoundCond(unsigned int)> f);

  /** @brief constructor using the coordinates on the boundary in 2D and 3D
   *
   * The function 'f' returns a BoundCond for given coordinates.
   * Use this constructor only if the boundary condition changes within one
   * boundary component. Otherwise use
   * BoundaryCondition(std::function<BoundCond(unsigned int)> f).
   */
  BoundaryCondition(std::function<BoundCond(const Point&)> f);

  /** @brief constructor in 2D using the parametrization on the boundary
   *
   * The function 'f' returns a BoundCond for a given component and parameter
   * (often called t). Use this constructor only if the boundary condition
   * changes within one boundary component. Otherwise use
   * BoundaryCondition(std::function<BoundCond(unsigned int)> f).
   */
  explicit BoundaryCondition(std::function<BoundCond(unsigned int, double)> f);

  /** @brief constructor in 3D using the parametrization on the boundary
   *
   * The function 'f' returns a BoundCond for a given component, part
   * and parameters (often called t and s). Use this constructor only if the
   * boundary condition changes within one boundary component. Otherwise use
   * BoundaryCondition(std::function<BoundCond(unsigned int)> f).
   */
  explicit BoundaryCondition(
    std::function<BoundCond(unsigned int, double, double)> f);

  /** @brief Use this for constant boundary condition on all boundary
   * components
   *
   * Note that this is essentially equivalent to calling
   * BoundaryCondition([](unsigned int component){ return bc;}). So this is
   * suitable for both 2D and 3D.
   */
  explicit BoundaryCondition(BoundCond bc);

  /** @brief copy constructor */
  BoundaryCondition(const BoundaryCondition&) = default;

  /** @brief move constructor */
  BoundaryCondition(BoundaryCondition&&) = default;

  /** @brief set the dimension of this example
   *
   * This is meaningful only if you used one of the constructors which are valid
   * for both 2D and 3D. If you used a constructor for 2D and try to change the
   * dimension to 3D (or vice versa), an exception is thrown as this is probably
   * not intended.
   *
   * @param d new dimension, must be 2 or 3, otherwise an exception is thrown
   */
  void set_dim(unsigned int d);

  /** @brief get the boundary condition on a given joint
   *
   * If 'this' does not depend on the parametrization or coordinates, only
   * the component of the BoundaryComponent in the joint is used. Otherwise
   * either coordinates or parameters are computed first and then used to get
   * the boundary condition.
   *
   * If this is not a boundary joint, an exception is thrown.
   *
   * @todo this method should take a BoundaryJoint. This class does not yet
   * exist, but it should be derived from TJoint, furthermore TBoundEdge and
   * TBoundFace should both derive from it.
   */
  BoundCond get_condition(const TJoint& joint) const;

 private:
  /** @brief dimension of this example (either 0, 2 or 3)
   *
   * dim = 0 means that this boundary condition does not depend on the
   * dimension, i.e. the dimension is not specified. This is usually the case if
   * the boundary condition only depends on the boundary component, or is even
   * constant on the entire boundary.
   */
  unsigned int dim;

  /** @brief function returning the boundary condition in 2D and 3D
   *
   * The unsigned parameter is the boundary component index. So if the boundary
   * conditions do not change within one boundary component this is probably
   * the best way to go.
   */
  std::function<BoundCond(unsigned int)> from_component;

  /** @brief function returning the boundary condition depending on the
   *         coordinates
   *
   * This function can be used if within one boundary component there are
   * different boundary conditions, depending on the coordinates.
   *
   * This is an alternative to 'from_parametrization2d' or
   * 'from_parametrization3d' respectively, which essentially does the same
   * using the parametrization of the boundary component.
   *
   * It is recommended to use 'from_component'. You can for example change your
   * Mesh such that there are two boundary components where there could be only
   * one.
   */
  std::function<BoundCond(const Point&)> from_coordinates;

  /** @brief function returning the boundary condition in 2D
   *
   * The unsigned parameter is the boundary component index. The double value
   * is the parameter (often called t).
   * This function can be used if within one boundary component there are
   * different boundary conditions, depending on the parameter. However note
   * that this is a bit dangerous, because on each edge of the Mesh only one
   * condition can be set, so the transition parameter also needs to represent
   * a vertex of the Mesh.
   *
   * It is recommended to use 'from_component'. You can for example change your
   * Mesh such that there are two boundary components where there could be only
   * one.
   */
  std::function<BoundCond(unsigned int, double)> from_parametrization2d;

  /** @brief function returning the boundary condition in 3D
   *
   * The unsigned parameter is the component component index. The double values
   * are the parameters (often called t and s). This function can be used if
   * within one boundary component there are different boundary conditions,
   * depending on the parameter. However note that this is a bit dangerous,
   * because on each face of the Mesh only one condition can be set, so the
   * transition parameters also needs to represent edges of the Mesh.
   *
   * It is recommended to use 'from_component'. You can for example change your
   * Mesh such that there are two boundary components where there could be only
   * one.
   */
  std::function<BoundCond(unsigned int, double, double)> from_parametrization3d;
  
  friend BoundaryData;
};

#endif // __BOUNDARYCONDITION__
