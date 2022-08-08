#include <Example.h>
#include <MooNMD_Io.h>

/* ************************************************************************* */
Example::Example(std::vector<BoundaryCondition>&& bc,
                 std::vector<BoundaryData>&& bd, PDECoefficients&& coeffs,
                 std::vector<AnalyticalFunction>&& exact)
  : bc(std::move(bc)), bd(std::move(bd)), coeffs(std::move(coeffs)),
    exact(std::move(exact))
{
}


/* ************************************************************************* */
/** a method to be used for the unit square, used in multiple examples */
bool is_unit_square(const TDomain&)
{
  ErrThrow("to be implemented");
  //   if(bdesc.dimension() != 2)
  //   {
  //     Output::warn<>("Example", "The given domain should be 2D. ");
  //     return false;
  //   }
  //   if(bdesc.n_parts() != 1)
  //   {
  //     Output::warn<>("Example", "There should be only one boundary part, not
  //     ",
  //                    bdesc.n_parts());
  //     return false;
  //   }
  //   if(bdesc.part(0)->n_bd_comps() != 4)
  //   {
  //     Output::warn("Example", "There should be four boundary components, not
  //     ",
  //              bdesc.part(0)->n_bd_comps());
  //     return false;
  //   }
  //   std::shared_ptr<BoundaryComponent> bc;
  //   {
  //     bc = bdesc.part(0)->get_bd_comp(0);
  //     if(bc->GetType() != Boundary_type::Line)
  //     {
  //       Output::warn("Example", "The 0-th boundary component should be a
  //       line,"
  //                      " not a ", bc->GetType());
  //       return false;
  //     }
  //     Point p((unsigned int) 2); // 2D
  //     bc->get_point_of_t(0.0, p);
  //     if(std::abs(p.x()) + std::abs(p.y()) > 1e-14)
  //     {
  //       Output::warn<>("Example", "The first point of the first boundary "
  //                      "component should be zero");
  //       return false;
  //     }
  //     bc->get_point_of_t(1.0, p);
  //     if(std::abs(p.x()-1.0) + std::abs(p.y()) > 1e-14)
  //     {
  //       Output::warn<>("Example", "The second point of the first boundary "
  //                      "component should be (1.0, 0.0)");
  //       return false;
  //     }
  //   }
  //   {
  //     bc = bdesc.part(0)->get_bd_comp(1);
  //     if(bc->GetType() != Boundary_type::Line)
  //     {
  //       Output::warn("Example", "The 1-st boundary component should be a
  //       line,"
  //                      " not a ", bc->GetType());
  //       return false;
  //     }
  //     Point p((unsigned int) 2); // 2D
  //     bc->get_point_of_t(0.0, p);
  //     if(std::abs(p.x()-1.0) + std::abs(p.y()) > 1e-14)
  //     {
  //       Output::warn<>("Example", "The first point of the second boundary "
  //                      "component should be (1.0, 0.0)");
  //       return false;
  //     }
  //     bc->get_point_of_t(1.0, p);
  //     if(std::abs(p.x()-1.0) + std::abs(p.y()-1.0) > 1e-14)
  //     {
  //       Output::warn<>("Example", "The second point of the second boundary "
  //                      "component should be (1.0, 1.0)");
  //       return false;
  //     }
  //   }
  //   {
  //     bc = bdesc.part(0)->get_bd_comp(2);
  //     if(bc->GetType() != Boundary_type::Line)
  //     {
  //       Output::warn("Example", "The 2-nd boundary component should be a
  //       line,"
  //                      " not a ", bc->GetType());
  //       return false;
  //     }
  //     Point p((unsigned int) 2); // 2D
  //     bc->get_point_of_t(0.0, p);
  //     if(std::abs(p.x()-1.0) + std::abs(p.y()-1.0) > 1e-14)
  //     {
  //       Output::warn<>("Example", "The first point of the third boundary "
  //                      "component should be (1.0, 1.0)");
  //       return false;
  //     }
  //     bc->get_point_of_t(1.0, p);
  //     if(std::abs(p.x()) + std::abs(p.y()-1.0) > 1e-14)
  //     {
  //       Output::warn<>("Example", "The second point of the third boundary "
  //                      "component should be (0.0, 1.0)");
  //       return false;
  //     }
  //   }
  //   {
  //     bc = bdesc.part(0)->get_bd_comp(3);
  //     if(bc->GetType() != Boundary_type::Line)
  //     {
  //       Output::warn("Example", "The 3-rd boundary component should be a
  //       line,"
  //                      " not a ", bc->GetType());
  //       return false;
  //     }
  //     Point p((unsigned int) 2); // 2D
  //     bc->get_point_of_t(0.0, p);
  //     if(std::abs(p.x()) + std::abs(p.y()-1.0) > 1e-14)
  //     {
  //       Output::warn<>("Example", "The first point of the fourth boundary "
  //                      "component should be (0.0, 1.0)");
  //       return false;
  //     }
  //     bc->get_point_of_t(1.0, p);
  //     if(std::abs(p.x()) + std::abs(p.y()) > 1e-14)
  //     {
  //       Output::warn<>("Example", "The second point of the fourth boundary "
  //                      "component should be (0.0, 0.0)");
  //       return false;
  //     }
  //   }
  return true;
}

/* ************************************************************************* */
bool Example::check_suitable_domain(const TDomain&) const
{
  Output::warn<>("Example::check_suitable_domain", "not yet implemented");
  
  return true;
}

/* ************************************************************************* */
int Example::get_example_from_database(const ParameterDatabase& param_db)
{
  try
  {
    return param_db["example"];
  }
  catch(...)
  {
    ErrThrow("could not find a parameter name 'example' in the "
             "ParameterDatabase. I am unable to create a suitable Example "
             "object without it.");
  }
}
