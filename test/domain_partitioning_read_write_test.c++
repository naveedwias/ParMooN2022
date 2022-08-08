#include <Database.h>
#include <Domain.h>
#include <MooNMD_Io.h>
#include <ParameterDatabase.h>
#include "BaseCell.h"
#include "ParMooN.h"

#include <mpi.h>
#include <stdio.h>

int main(int /*argc*/, char **/*argv*/)
{
  ParameterDatabase db = parmoon::parmoon_initialize();

  int my_rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  Output::setVerbosity(5);

  // add some parameters to the database.
  db.add("read_metis", false, "");
  db.add("read_metis_file", "metis_domain_decomp.out", "");
  db.add("write_metis", true, "");
  db.add("write_metis_file", "metis_domain_decomp.out", "");

  db.add("refinement_n_initial_steps", (size_t) 2, "");
  db.add("boundary_file", "Default_UnitCube", "");
  db.add("geo_file", "Default_UnitCube_Hexa", "");

  // Construct the "original" domain.
  TDomain domain_original(db);

  // Intial refinement.
  std::list<TCollection* > grids
    = domain_original.refine_and_get_hierarchy_of_collections(db);

  domain_original.print_info("original domain");

  // Construct the "copy" domain - first change the database accordingly.
  db["read_metis"] = true;
  db["write_metis"] = false;

  TDomain domain_copycat(db);
  std::list<TCollection* > grids_copy
    = domain_copycat.refine_and_get_hierarchy_of_collections(db);

  domain_copycat.print_info("copied domain");

  // compare the two domain objects:
  if(grids.size() != grids_copy.size())
  {
    ErrThrow("wrong number of collections");
  }
  auto n_cells = grids.back()->GetN_Cells();
  if(n_cells != grids_copy.back()->GetN_Cells())
  {
    ErrThrow("wrong number of cells after partitioning ",
             grids.back()->GetN_Cells(), " != ",
             grids_copy.back()->GetN_Cells());
  }
  for(int i = 0; i < n_cells; ++i)
  {
    const TBaseCell* c_orig = grids.back()->GetCell(i);
    const TBaseCell* c_copy = grids_copy.back()->GetCell(i);
    // in 3D there are only either tetrahedra or hexahedra
    // so we can assume that the number of vertices in the two cell is equal
    auto n_vert = c_orig->GetN_Vertices();
    for(int v = 0; v < n_vert; ++v)
    {
      const TVertex* v_orig = c_orig->GetVertex(v);
      const TVertex* v_copy = c_copy->GetVertex(v);
      if(  v_orig->GetX() != v_copy->GetX()
        || v_orig->GetY() != v_copy->GetY()
        || v_orig->GetZ() != v_copy->GetZ() )
      {
        Output::print("Found two different vertices, Cell ", i, " vertex ", v,
                      v_orig, "  != ", v_copy);
        ErrThrow("wrong mesh partitioning");
      }
    }
  }

  // tidy up in the build directory
  if(my_rank == 0)
	  remove("metis_domain_decomp.out");

  parmoon::parmoon_finalize();
}

