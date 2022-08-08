#include "Domain.h"
#include "MooNMD_Io.h"
#include "NavierStokes.h"
#include "Database.h"
#include "ParMooN.h"

template <int d>
void check(ParameterDatabase db, TDomain &domain, int velocity_order)
{
  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  // automatically choose inf-sup stable pressure space
  TDatabase::ParamDB->PRESSURE_SPACE = -4711;
  
  NavierStokes<d> nse(domain, db);
  Output::print("created nse2d object");
  nse.assemble_linear_terms();
  nse.solve();
  nse.output();
  
  auto div_error = nse.get_errors()[1];
  if(div_error > 1.e-12)
  {
    ErrThrow("Scott-Vogelius elements should lead to zero divergence, but ",
             div_error, " has been computed");
  }
}

int main(int , char** )
{
  ParameterDatabase db = parmoon::parmoon_initialize();
  db.merge(TDomain::default_domain_parameters(), true);
#ifdef __2D__
  db.merge(Example2D::default_example_database(), true);
#else // 3D
  db.merge(Example3D::default_example_database(), true);
#endif // 2D
  db["refinement_n_initial_steps"] = 2;
  db["refinement_final_step_barycentric"] = true;
  db["example"] = 1;
  Output::setVerbosity(3);
  
#ifdef __2D__
  db["boundary_file"] = "Default_UnitSquare";
  db["geo_file"] = "TwoTriangles";
#else // 3D
  db["boundary_file"] = "Default_UnitCube";
  db["geo_file"] = "Default_UnitCube_Tetra";
  db["refinement_n_initial_steps"] = 1;
#endif // 2D
  //db.add("output_write_vtk", true, "");
  
  TDomain domain(db);
  // refinement
  domain.refine_and_get_hierarchy_of_collections(db);
  
#ifdef __2D__
  check<2>(db, domain, 12);
  check<2>(db, domain, 13);
  check<2>(db, domain, 14);
#else // 2D -> 3D
  check<3>(db, domain, 13);
  // check<3>(db, domain, 14); // P4 on tetrahedra not implemented!
#endif // 2D
  parmoon::parmoon_finalize();
}

