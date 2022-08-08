#include <limits> // std::numeric_limits<double>::quiet_NaN
#include <BoundaryCondition.h>
#include <MooNMD_Io.h>
#ifdef __2D__
#include <BoundEdge.h>
#else // __3D__
#include <BoundFace.h>
#endif // 3D

using namespace parmoon;

/** ************************************************************************ */
BoundaryCondition::BoundaryCondition(std::function<BoundCond(unsigned int)> f)
  : dim(0), from_component(f), from_coordinates(), from_parametrization2d(),
    from_parametrization3d()
{
}

/** ************************************************************************ */
BoundaryCondition::BoundaryCondition(
    std::function<BoundCond(unsigned int, double)> f)
  : dim(2), from_component(), from_coordinates(), from_parametrization2d(f),
    from_parametrization3d()
{
}

/** ************************************************************************ */
BoundaryCondition::BoundaryCondition(std::function<BoundCond(const Point&)> f)
  : dim(0), from_component(), from_coordinates(f), from_parametrization2d(),
    from_parametrization3d()
{
}

/** ************************************************************************ */
BoundaryCondition::BoundaryCondition(
    std::function<BoundCond(unsigned int, double, double)> f)
  : dim(3), from_component(), from_coordinates(), from_parametrization2d(),
    from_parametrization3d(f)
{
}

/** ************************************************************************ */
BoundaryCondition::BoundaryCondition(BoundCond bc)
  : BoundaryCondition([bc](unsigned int)
                      {
                        return bc;
                      })
{
}

/** ************************************************************************ */
void BoundaryCondition::set_dim(unsigned int d)
{
  if(this->dim == 0)
  {
    if(d == 2 || d == 3)
      this->dim = d;
    else
    {
      // error, only dimensions 2 and 3 are allowed
      ErrThrow("cannot change dimension in BoundaryCondition to ", d);
    }
  }
  else if(this->dim == d)
    return; // nothing needs to be done here
  else
  {
    // error, cannot change the dimension unless it was unspecified before (0)
    ErrThrow("cannot change the dimension in BoundaryCondition from ",
             this->dim, " to ", d);
  }
}

/** ************************************************************************ */
BoundCond BoundaryCondition::get_condition(const TJoint& joint) const
{
  if(joint.InnerJoint())
  {
    ErrThrow("cannot find a boundary condition for an inner joint");
  }

#ifdef __2D__
  TBoundEdge* bd_joint = (TBoundEdge*)&joint;
  bool two_d = true; // indicate space dimension
  if(this->dim == 3)
    ErrThrow("3D BoundaryCondition called on a TBoundEdge (which is 2D)");
#else
  TBoundFace* bd_joint = (TBoundFace*)&joint;
  bool two_d = false; // indicate space dimension
  if(this->dim == 2)
    ErrThrow("2D BoundaryCondition called on a TBoundFace (which is 3D)");
#endif
  unsigned int component = bd_joint->GetBoundComp()->GetID();

  if(this->from_component) // is callable
  {
    return this->from_component(component);
  }

  // the parameter (in the parametetrization of the boundary) for the center of
  // the joint
  double t_center = 0.0;
#ifdef __2D__
  {
    double t0, t1;
    bd_joint->GetParameters(t0, t1);
    t_center = 0.5 * (t0 + t1);
  }
#else // 3D
  ErrThrow("BoundaryCondition::get_condition not implemented in 3D in the "
           "case where the boundary condition does not only depend on the "
           "boundary component");
#endif

  if(this->from_coordinates) // is callable
  {
    // create a all-zero Point where third coordinate is set to nan in 2D
    Point p((unsigned int)(two_d ? 2 : 3));
    /// @todo the joint should be able to compute its center point, even for
    /// isoparametric boundaries, we do it by hand here.

    // fill Point p with correct values as center of Joint.
    if(bd_joint->GetType() == IsoBoundEdge)
    {
      // for curved boundaries, the center of the joint is in general not on
      // the boundary. So in that case we get the center parameter and compute
      // the corresponding coordinates
      ErrThrow("iso parametric boundaries are not supported currently");
    }
#ifdef __2D__
    bd_joint->GetXYofT(t_center, p.x, p.y);
#else
// ??
#endif
    return this->from_coordinates(p);
  }
  else if(this->from_parametrization2d) // is callable
  {
    // evaluate at center of joint
    return this->from_parametrization2d(component, t_center);
  }
  else if(this->from_parametrization3d) // is callable
  {
    ErrThrow("boundary conditions depending on the boundary parametrization "
             " in 3D are not yet supported");
  }
  else
  {
    ErrThrow("no function specified in BoundaryCondition");
  }
}


/** ************************************************************************ */
