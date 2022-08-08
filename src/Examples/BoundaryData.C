#include <BoundaryData.h>
#include <MooNMD_Io.h>
#include <PDECoefficients.h> // for enum class Problem_type

#include <BoundEdge.h>
#include <BoundFace.h>

using namespace parmoon;

/* ************************************************************************* */
BoundaryData::BoundaryData(std::function<double(unsigned int)> f)
 : depends_on_time(false), from_component(f), from_parametrization2d(),
   from_parametrization3d(), from_coordinates()
{
}

/* ************************************************************************* */
#ifdef __2D__
BoundaryData::BoundaryData(
    std::function<double(unsigned int, double, double)> f)
  : from_component(), from_parametrization2d(f), from_parametrization3d(),
    from_coordinates()
{
}
#endif // 2D

/* ************************************************************************* */
#ifdef __3D__
BoundaryData::BoundaryData(
    std::function<double(unsigned int, double, double, double)> f)
  : from_component(), from_parametrization2d(), from_parametrization3d(f),
    from_coordinates()
{
}
#endif // 3D

/* ************************************************************************* */
BoundaryData::BoundaryData(std::function<double(const Point&, double)> f)
  : from_component(), from_parametrization2d(), from_parametrization3d(),
    from_coordinates(f)
{
}

/* ************************************************************************* */
BoundaryData::BoundaryData(double v)
  : BoundaryData([v](unsigned int){return v;})
{
}

/* ************************************************************************* */
BoundaryData::BoundaryData(const BoundaryCondition& bc,
                           const AnalyticalFunction& exact, Problem_type type)
  : from_component(), from_parametrization2d(), from_parametrization3d(),
    from_coordinates()
{
  /// @todo we need a description of the boundary to compute normals and 
  /// possibly convert parametrizations to coordinates and vice versa.
  if(!exact.exists())
  {
    ErrThrow("can not create BoundaryData from an AnalyticalFunction which has "
             "been default constructed. That means it is not known.");
  }
  if(bc.from_component && exact.is_constant()) // callable
  {
    switch(type)
    {
      case Problem_type::ConvDiffReac:
        this->from_component = [bc, exact](unsigned int c)
          {
            auto bd_condition = bc.from_component(c);
            unsigned int space_dim = 2;
            FunctionEvaluation f_eval(space_dim, 1);
            Point p(space_dim); // exact is constant so p can be anything
            exact.get(p, f_eval);
            switch(bd_condition)
            {
              case DIRICHLET:
                return f_eval.get(0, 0);
                break;
              case NEUMANN:
                ErrThrow("BoundaryCondition other DIRICHLET not yet supported");
                break;
              default:
                ErrThrow("BoundaryCondition other than NEUMANN or DIRICHLET "
                         "not yet supported");
            }
            return 0.;
          };
        break;
      default:
        ErrThrow("unable to generate BoundaryData from BoundaryCondition and "
                 "AnalyticalFunction for type ", type);
    }
  }
  else
    ErrThrow("automatic computation of boundary data from boundary conditions "
             "and exact solution is not yet implemented");
}


/* ************************************************************************* */
double BoundaryData::get_data(const TJoint& joint, double s, double time) const
{
#ifdef __3D__
  ErrThrow("BoundaryData::get_data(const Joint&, double, double) is only for "
           "2D");
  TBoundFace* bd_joint = (TBoundFace*)&joint; // to avoid compiler errors
#endif // 3D
#ifdef __2D__
  TBoundEdge* bd_joint = (TBoundEdge*)&joint;
#endif
  unsigned int component = bd_joint->GetBoundComp()->GetID();
  
  if(this->from_component) // callable
  {
    return this->from_component(component);
  }
  
  // the parameter (in the parametetrization of the boundary) for the desired
  // point on the joint
  double t = 0.0;
  {
    double t0, t1;
#ifdef __2D__
    bd_joint->GetParameters(t0, t1);
#endif
    t = t0 + (s+1.)*0.5 * (t1 - t0);
  }
  
  if(this->from_parametrization2d) // callable
  {
    return this->from_parametrization2d(component, t, time);
  }
  else if(this->from_coordinates) // callable
  {
    // compute x,y, from s and joint
    Point p(0.0, 0.0);
#ifdef __2D__
    bd_joint->GetXYofT(t, p.x, p.y);
#endif
    return this->from_coordinates(p, time);
  }
  else if(this->from_parametrization3d) // callable
  {
    ErrThrow("BoundaryData::get_data(const Joint&, double, double) is only for "
             "2D");
  }
  else
  {
    ErrThrow("no 2D BoundaryData function available");
  }
}

/* ************************************************************************* */
double BoundaryData::get_data(const TJoint& joint, std::pair<double, double> ts,
                              double time) const
{
  double t = ts.first;
  double s = ts.second;
#ifdef __2D__
  ErrThrow("BoundaryData::get_data(const Joint&, std::pair<double,double>, "
           "double) is only for 3D");
  TBoundEdge* bd_joint = (TBoundEdge*)&joint; // to avoid compiler errors
#endif // 2D
#ifdef __3D__
  TBoundFace* bd_joint = (TBoundFace*)&joint;
#endif
  unsigned int component = bd_joint->GetBoundComp()->GetID();
  
  if(this->from_component) // callable
  {
    return this->from_component(component);
  }

  if(this->from_parametrization3d) // callable
  {
    return this->from_parametrization3d(component,t, s, time);
  }
  else if(this->from_coordinates) // callable
  {
    // calculate point from parametrization t,s
    // point_on_boundary
    Point point((unsigned int)3); // space dimension
#ifdef __3D__
    bd_joint->GetXYZofTS(t, s, point.x, point.y, point.z);
#endif
    return this->from_coordinates(point, time);
  }
  else if(this->from_parametrization2d)
  {
    ErrThrow("BoundaryData::get_data(const Joint&, double, double, double) is "
             "only for 3D");
  }
  else
  {
    ErrThrow("no 3D BoundaryData function available");
  }
}

/* ************************************************************************* */
double BoundaryData::get_data(const TJoint& joint, const Point& point, 
                              double time) const
{
#ifdef __2D__
  TBoundEdge* bd_joint = (TBoundEdge*)&joint;
#else
  TBoundFace* bd_joint = (TBoundFace*)&joint;
#endif
  unsigned int component = bd_joint->GetBoundComp()->GetID();
  
  if(this->from_component) // callable
  {
    return this->from_component(component);
  }
  else if(this->from_coordinates) // callable
  {
    return this->from_coordinates(point, time);
  }
  else if(this->from_parametrization2d) // callable
  {
    // compute t from point
    double t = 0;
#ifdef __2D__
    bd_joint->GetTofXY(point.x, point.y, t);
#endif
    return this->from_parametrization2d(component, t, time);
  }
  else if(this->from_parametrization3d)
  {
    // compute t,s from Point 
    double t = 0.0, s = 0.0;
#ifdef __3D__
    bd_joint->GetTSofXYZ(point.x, point.y, point.z, t, s);
#endif
    return this->from_parametrization3d(component, t, s, time);
  }
  else
  {
    ErrThrow("no BoundaryData function using coordinates available");
  }
}

/* ************************************************************************* */
