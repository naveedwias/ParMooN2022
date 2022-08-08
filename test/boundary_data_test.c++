#include <MooNMD_Io.h>
#include <cmath>
#include <list>
#include <BoundaryData.h>
#include <PDECoefficients.h> // for enum class Problem_type

#include <BoundEdge.h>
#include <BdCircle.h>
#include <BdLine.h>
#include <BoundFace.h>
#include <BdPlane.h>
#include <BdWall.h>

using namespace parmoon;

void test_BdData_fromBdCondition_and_AnalyticalFunction()
{
  double value = 1.0;
  BoundaryCondition bc_component(DIRICHLET);
  AnalyticalFunction af(value);
  
  // This is currently the only case where this seems to work
  BoundaryData bd(bc_component, af, Problem_type::ConvDiffReac);
  
  
}

int main(int, char **)
{
  // at first we create a few joints in 2D and 3D to have object on which we 
  // can later call BoundaryCondition::get_condition.
  
  // create a TBdCircle object
  TBdCircle circle(0);
  // center at (1,1), Radius 0.5 (in both directions making it a circle rather 
  // than an ellipse), starting at an angle Phi=pi, ending at Phi=2pi. This 
  // describes a lower half circle.
  circle.SetParams(1.0, 1.0, 0.5, 0.5, M_PI, 2*M_PI);
  
  // create a TBdLine object
  TBdLine line(1);
  // horizontal line starting at (1.5,1) ending at (0.5,1). These are the end 
  // points of the half circle defined above
  line.SetParams(1.5, 1, -1., 0.);
  
  std::list<TBoundEdge> bd_edges;
  // some edges with end point on the circle
  bd_edges.push_back(TBoundEdge(&circle, 0.4, 0.6));
  bd_edges.push_back(TBoundEdge(&circle, 0.6, 0.8));
  // some edges with end points on the line
  bd_edges.push_back(TBoundEdge(&line, 0.0, 0.4));
  bd_edges.push_back(TBoundEdge(&line, 0.4, 0.9));
  
#ifdef __3D__
  // create a TBdCylinder
  TBdPlane plane(0);
  // plane through point (1,1,1) and (2,1,1) and normal (0, 0.5,2.)
  plane.SetParams(1.,1.,1., 1,0.,0., 0., 0.5, 2.);
  
  // create a TBdWall using the circle
  TBdWall wall_circle(1, &circle);
  wall_circle.SetParams(2., 1., 4., 0.);
  
  // create a TBdWall using the line
  TBdWall wall_line(1, &line);
  wall_line.SetParams(2., 1., 4., 0.);
  
  std::list<TBoundFace> bd_faces;
  // some faces on the TBoundComp3D objects plane, wall_circle, wall_line
  bd_faces.push_back(TBoundFace(&plane));
  bd_faces.push_back(TBoundFace(&wall_circle));
  bd_faces.push_back(TBoundFace(&wall_line));
#endif // 3D
  
  
  /////////////////////////////////////////////////////////////////////////////
  // done creating boundary joints, now the tests start
  
  
  // constant boundary data
  Output::print("testing constant boundary data");
  double value = 1.0;
  BoundaryData bd_constant(value);
  Point p(1., 2.);
#ifdef __2D__
  for(auto e : bd_edges)
  {
    if(bd_constant.get_data(e, p) != value)
      ErrThrow("wrong boundary data ", bd_constant.get_data(e, p), "  ", value);
    if(bd_constant.get_data(e, 0.5) != value)
      ErrThrow("wrong boundary data ", bd_constant.get_data(e, p, 0.5), "  ",
               value);
  }
#elif defined __3D__
  for(auto e : bd_edges)
  {
    if(bd_constant.get_data(e, p) != value)
      ErrThrow("wrong boundary data ", bd_constant.get_data(e, p), "  ", value);
    try
    {
      if(bd_constant.get_data(e, 0.5) != value)
        ErrThrow("wrong boundary data ", bd_constant.get_data(e, p, 0.5), "  ",
                 value);
      Output::warn("BoundaryData::get_data(TJoint, double)", 
                   "calling get_data with a 2D parametrization in 3D should "
                   "throw an exception");
      return 1;
    }
    catch(...)
    {
      //this was expected
    }
    auto ts = std::make_pair(0.5, 0.1);
    if(bd_constant.get_data(e, ts, 0.) != value)
      ErrThrow("wrong boundary data ", bd_constant.get_data(e, ts, 0.),
               "  ", value);
  }
  p = Point(1., 2., 3.);
  for(auto e : bd_faces)
  {
    if(bd_constant.get_data(e, p) != value)
      ErrThrow("wrong boundary data ", bd_constant.get_data(e, p), "  ", value);
  }
#endif // 3D
  
  
  Output::print("testing boundary data which only depends on the component");
  auto from_component = [](unsigned int c){ return c == 0 ? 1. : -1.; };
  BoundaryData bd_component(from_component);
  for(auto e : bd_edges)
  {
    unsigned int c = e.GetBoundComp()->GetID(); // boundary component index
    if(bd_component.get_data(e, p) != from_component(c))
      ErrThrow("wrong boundary data ", bd_component.get_data(e, p), "  ",
               from_component(c), " component ", c);
#ifdef __2D__
    if(bd_component.get_data(e, 0.5) != from_component(c))
      ErrThrow("wrong boundary data ", bd_component.get_data(e, p, 0.5), "  ",
               from_component(c));
#endif // 2D
  }
  
  
  Output::print("testing boundary data which depends on the parametrization");
#ifdef __2D__
  auto from_parametrization2d = [](unsigned int c, double t, double /*time*/)
                                 { return c == 0 ? t : 1-t; };
  BoundaryData bd_param2d(from_parametrization2d);
  p = Point(1.5, 1.);
  for(auto e : bd_edges)
  {
    unsigned int c = e.GetBoundComp()->GetID(); // boundary component index
    double param = (c == 0 ? 1. : 0.);
    if(bd_param2d.get_data(e, p) != from_parametrization2d(c, param, 0.0))
      ErrThrow("wrong boundary data ", bd_param2d.get_data(e, p), "  ",
               from_parametrization2d(c, param, 0.0), " component ", c);
  }
#endif // 2D
#ifdef __3D__
  auto from_parametrization3d = [](unsigned int c, double t, double s, double)
                                 { return c == 0 ? t+s : 1-t-s; };
  BoundaryData bd_param2d(from_parametrization3d);
  p = Point(1.5, 1., 0.5);
  // still to do
#endif
  
  
  Output::print("testing boundary data which depends on the coordinates");
  auto from_coordinates = [](const Point & p, double)
                          { return p.x > 0.3 ? 1. : -1.; };
  BoundaryData bd_coord(from_coordinates);
  p = Point(1.5, 1., 0.5);
  for(auto e : bd_edges)
  {
    unsigned int c = e.GetBoundComp()->GetID(); // boundary component index
    if(bd_coord.get_data(e, p) != from_coordinates(p, 0.0))
      ErrThrow("wrong boundary data ", bd_coord.get_data(e, p), "  ",
               from_coordinates(p, 0.0), " component ", c);
  }
  
  
  Output::print("testing boundary data originating from a BoundaryCondition "
                "and an AnalyticalFunction");
  { // local scope
    double value = 1.0;
    BoundaryCondition bc_component(DIRICHLET);
    AnalyticalFunction af(value);
    
    // This is currently the only case where this seems to work
    BoundaryData bd(bc_component, af, Problem_type::ConvDiffReac);
    for(auto e : bd_edges)
    {
      if(bd.get_data(e, p) != value)
        ErrThrow("wrong boundary data ", bd.get_data(e, p), "  ", value);
#ifdef __2D__
      if(bd.get_data(e, 0.5) != value)
        ErrThrow("wrong boundary data ", bd.get_data(e, p, 0.5), "  ", value);
#elif defined __3D__
      auto ts = std::make_pair(0.5, 0.1);
      if(bd.get_data(e, ts, 0.) != value)
        ErrThrow("wrong boundary data ", bd.get_data(e, ts, 0.), "  ", value);
#endif
    }
  }
  test_BdData_fromBdCondition_and_AnalyticalFunction();
}
