#include <MooNMD_Io.h>
#include <cmath>
#include <list>
#include <BoundaryCondition.h>

#include <BoundEdge.h>
#include <BdCircle.h>
#include <BdLine.h>
#include <BoundFace.h>
#include <BdPlane.h>
#include <BdWall.h>

using namespace parmoon;

int main(int, char**)
{
  double tol = 1.e-13;
  // at first we create a few joints in 2D and 3D to have object on which we 
  // can later call BoundaryCondition::get_condition.
  
  // create a TBdCircle object
  TBdCircle circle(0);
  // center at (1,1), Radius 0.5 (in both directions making it a circle rather 
  // than an ellipse), starting at an angle Phi=pi, ending at Phi=2pi. This 
  // describes a lower half circle.
  circle.SetParams(1.0, 1.0, 0.5, 0.5, M_PI, 2*M_PI);
  // testing the circle a little bit
  {
    double t = 0.5;
    double x, y;
    circle.GetXYofT(t, x, y);
    if(std::abs(x-1.) > tol || std::abs(y-0.5) > tol)
      ErrThrow("TBdCircle::GetXYofT is not working correctly ", x, "  ", y);
    circle.GetTofXY(x, y, t);
    if(std::abs(t-0.5) > tol)
      ErrThrow("TBdCircle::GetTofXY is not working correctly ", t);
  }
  
  // create a TBdLine object
  TBdLine line(1);
  // horizontal line starting at (1.5,1) ending at (0.5,1). These are the end 
  // points of the half circle defined above
  line.SetParams(1.5, 1, -1., 0.);
  // testing the line a little bit
  {
    double t = 0.5;
    double x, y;
    line.GetXYofT(t, x, y);
    if(std::abs(x-1.) > tol || std::abs(y-1.) > tol)
      ErrThrow("TBdLine::GetXYofT is not working correctly ", x, "  ", y);
    line.GetTofXY(x, y, t);
    if(std::abs(t-0.5) > tol)
      ErrThrow("TBdLine::GetTofXY is not working correctly ", t);
  }
  
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
  // testing the plane a little bit
  {
    double t = 0.5, s = -1.;
    double x,y,z;
    plane.GetXYZofTS(t, s, x, y, z);
    if(std::abs(x-1.5) > tol || std::abs(y-3.) > tol || std::abs(z-0.5) > tol)
      ErrThrow("TBdPlane::GetXYZofTS does not work correctly ", x, " ", y, " ",
               z);
    plane.GetTSofXYZ(x, y, z, t, s);
    if(std::abs(t-0.5) > tol || std::abs(s+1.) > tol)
      ErrThrow("TBdPlane::GetTSofXYZ not correctly working ", t, "  ", s);
  }
  
  // create a TBdWall using the circle
  TBdWall wall_circle(1, &circle);
  wall_circle.SetParams(2., 1., 4., 0.);
  //testing the wall_circle a little bit
  {
    double t = 0.5, s = 0.6;
    double x,y,z;
    wall_circle.GetXYZofTS(t, s, x, y, z);
    if(std::abs(x-2.2) > tol || std::abs(y-1.1) > tol || std::abs(z-2.4) > tol)
      ErrThrow("TBdWall::GetXYZofTS does not work correctly ", x, " ", y, " ", 
               z);
    wall_circle.GetTSofXYZ(x, y, z, t, s);
    if(std::abs(t-0.5) > tol || std::abs(s-0.6) > tol)
      ErrThrow("TBdWall::GetTSofXYZ not correctly working ", t, "  ", s);
  }
  
  // create a TBdWall using the line
  TBdWall wall_line(1, &line);
  wall_line.SetParams(2., 1., 4., 0.);
  //testing the wall_line a little bit
  {
    double t = 0.5, s = 0.6;
    double x,y,z;
    wall_line.GetXYZofTS(t, s, x, y, z);
    if(std::abs(x-2.2) > tol || std::abs(y-1.6) > tol || std::abs(z-2.4) > tol)
      ErrThrow("TBdWall::GetXYZofTS does not work correctly ", x, " ", y, " ", 
               z);
    wall_line.GetTSofXYZ(x, y, z, t, s);
    if(std::abs(t-0.5) > tol || std::abs(s-0.6) > tol)
      ErrThrow("TBdWall::GetTSofXYZ not correctly working ", t, "  ", s);
  }
  
  std::list<TBoundFace> bd_faces;
  // some faces on the TBoundComp3D objects plane, wall_circle, wall_line
  bd_faces.push_back(TBoundFace(&plane));
  bd_faces.push_back(TBoundFace(&wall_circle));
  bd_faces.push_back(TBoundFace(&wall_line));
#endif // 3D
  
  
  /////////////////////////////////////////////////////////////////////////////
  // done creating boundary joints, now the tests start
  
  Output::print("Testing constant BoundaryCondition");
  BoundaryCondition bc_constant(NEUMANN);
  for(auto e : bd_edges)
  {
    if(bc_constant.get_condition(e) != NEUMANN)
      ErrThrow("wrong boundary condition ", e.GetType());
  }
  bc_constant.set_dim(2);
  // try to reset the dimension
  try
  {
    bc_constant.set_dim(3);
    Output::warn("BoundaryCondition::set_dim", 
                 "switching the dimension should not be allowed");
    return 1;
  }
  catch(std::runtime_error err)
  {
    // this is what we expected
  }
  catch(...)
  {
    Output::warn("BoundaryCondition::set_dim", "caught unexpected exception");
  }
  
  
  Output::print("Testing BoundaryCondition depending on the boundary "
                "component");
  auto from_component = [](unsigned int c){return c==0 ? NEUMANN : DIRICHLET;};
  BoundaryCondition bc_component(from_component);
  for(auto e : bd_edges)
  {
    unsigned int component_id = e.GetBoundComp()->GetID();
    if(bc_component.get_condition(e) != from_component(component_id))
      ErrThrow("wrong boundary condition ", e.GetType());
  }
#ifdef __3D__
  for(auto e : bd_faces)
  {
    unsigned int component_id = e.GetBoundComp()->GetID();
    if(bc_component.get_condition(e) != from_component(component_id))
      ErrThrow("wrong boundary condition ", e.GetType());
  }
#endif // 3D
  
  
#ifdef __3D__
  /// @todo test BoundaryCondition depending on coordinates and parametrization
  /// in 3D
  return 0;
#endif // 3D
  
  Output::print("Testing BoundaryCondition depending on the coordinates");
  auto from_point = [](const Point& p)
                    {return p.x > 0.3 ? NEUMANN : DIRICHLET;};
  BoundaryCondition bc_point(from_point);
  for(auto e : bd_edges)
  {
    if(bc_point.get_condition(e) != NEUMANN)
      ErrThrow("wrong boundary condition ", e.GetType());
  }
  
  Output::print("Testing BoundaryCondition depending on the parametrization");
  auto from_parametrization = [](unsigned int component, double t)
                              { Output::print("HALLO ", component, " ", t);
                                return component == 0 ? NEUMANN : DIRICHLET; };
  BoundaryCondition bc_param(from_parametrization);
  for(auto e : bd_edges)
  {
    unsigned int component_id = e.GetBoundComp()->GetID();
    auto bc = component_id== 0 ? NEUMANN : DIRICHLET;
    if(bc_param.get_condition(e) != bc)
      ErrThrow("wrong boundary condition ", e.GetType());
  }
}
