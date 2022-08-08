#include <Database.h>
#include <Domain.h>
#include <MooNMD_Io.h>
#include <FEFunction3D.h>
#include <BlockVector.h>
#include <Database.h>
#include <MainUtilities.h>
#include "ParMooN.h"

#ifdef _MPI
#include <mpi.h>
#endif
// some functions dummy given as an analytical representation
void analytic_function(double x, double y, double z, double * values)
{
  values[0] = x*y*z;
}


int main(int argc, char** argv)
{
  ParameterDatabase db = parmoon::parmoon_initialize(argc, argv);
  Output::setVerbosity(2);

  // set a few parameters in the database to construct a unit cube with one
  // cell. then refine uniformly twice
  db.add("boundary_file", "Default_UnitCube", "");
   db.add("geo_file", "Default_UnitCube_Hexa", "",
          {"Default_UnitCube_Hexa","Default_UnitCube_Tetra"});
   db.add("refinement_n_initial_steps", (size_t) 1,"", (size_t) 0, (size_t) 2);
  // construct a domain object
  TDomain domain(db);
  //Refinement of grid in 3D for MPI ist different than for the sequential case in 2D of read_write_fe.
  // refine grid
  std::list<TCollection*> grid_collections
     = domain.refine_and_get_hierarchy_of_collections(db);
  // create finite element space, for that we need a
  // collection of grid cells ( = grid)
  // type of finite element on each cell (here: Q2)
  size_t ansatz_order = 2;
  std::shared_ptr<const TFESpace3D> fe_space(
    new TFESpace3D(grid_collections.back(), "feSpace_",
                   BoundConditionNoBoundCondition, ansatz_order));

  //fe_space->info();

  // create some vector which could be seen as a finite element function
  // representation
  BlockVector v1(fe_space->get_n_dof());
  {
    // fill the vector with some values, according to some given function
    TFEFunction3D f(fe_space, "dummy", v1.get_entries());
    f.Interpolate(analytic_function);
  }

  // write the vector
  std::stringstream out_in_stream;
  v1.write_to_stream(out_in_stream);

// create second vector, read into it
  BlockVector v2(fe_space->get_n_dof());
  v2.read_from_stream(out_in_stream);

  // substract first vector (difference should be zero now)
  v2.add_scaled(v1, -1.);

  double err = v2.norm();
  if(err != 0.0)
  {
    ErrThrow("writing and then reading a BlockVector failed with error ", err);
  }
  Output::print("test successful");
  parmoon::parmoon_finalize();
}
