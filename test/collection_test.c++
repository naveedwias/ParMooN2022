#include "Domain.h"
#include "Database.h"
#include <ParMooN_repository_info.h>
#include "MooNMD_Io.h"
#include "ParMooN.h"

const std::string path_to_repo = parmoon::source_directory;
const std::string path_to_meshes = path_to_repo + "/data/mesh/";


struct TestObject
{
    std::string Geo_file;
    std::string Boundary_file;
    bool has_simplex_cell;
};

TestObject prepare_test0()
{
#ifdef __2D__
  std::string geo_file = path_to_meshes + "unit_square/unit_square_quad6.mesh";
  std::string bd_file = path_to_meshes + "unit_square/unit_square.PRM";
#else
  std::string geo_file = "Default_UnitCube_Hexa";
  std::string bd_file = "Default_UnitCube";
#endif
  return TestObject{geo_file, bd_file, false};
}
TestObject prepare_test1()
{
#ifdef __2D__
  std::string geo_file = path_to_meshes + "unit_square/unit_square_tria6.mesh";
  std::string bd_file = path_to_meshes + "unit_square/unit_square.PRM";
#else
  std::string geo_file = "Default_UnitCube_Tetra";
  std::string bd_file = "Default_UnitCube";
#endif
  return TestObject{geo_file, bd_file, true};
}
TestObject prepare_test2()
{
#ifdef __2D__
  std::string geo_file = path_to_meshes + "Hemker_tria.mesh";
  std::string bd_file = path_to_meshes + "Hemker.PRM";
#else
  std::string geo_file = path_to_meshes + "cylinder.3d.3K.mesh";
  std::string bd_file = "";
#endif
  return TestObject{geo_file, bd_file, true};
}
TestObject prepare_test3()
{
#ifdef __2D__
  std::string geo_file = path_to_meshes + "Hemker_quad.mesh";
  std::string bd_file = path_to_meshes + "Hemker.PRM";
  auto has_simplex_cell = false;
#else
  std::string geo_file = path_to_meshes + "channel.3d.mesh";
  std::string bd_file = "";
  auto has_simplex_cell = true;
#endif
  return TestObject{geo_file, bd_file, has_simplex_cell};
}


void check_collection(TCollection& coll, bool has_simplex_cell = true)
{
  auto n_cells = coll.GetN_Cells();
  Output::print("check collection with ", n_cells, " cells");
  for(int i = 0; i < n_cells; ++i)
  {
    auto cell = coll.GetCell(i);
    if(i != coll.get_cell_index(cell))
    {
      ErrThrow("wrong cell index, ", i, " ", coll.get_cell_index(cell));
    }
  }
  if (has_simplex_cell)
  {
    if (!coll.has_tria_cell() or !coll.has_tetra_cell())
    {
      ErrThrow("has_tria_cell() or has_tetra_cell() failed");
    }
  }
  else
  {
    if (!coll.has_quad_cell() or !coll.has_quad_cell())
    {
      ErrThrow("has_tria_cell() or has_tetra_cell() failed");
    }
  }
}

int main()
{
  parmoon::parmoon_initialize();
  std::vector<TestObject> TestObjects = {prepare_test0(), prepare_test1(),
                                         prepare_test2(), prepare_test3()};
  for(auto &Test : TestObjects)
  {
    Output::print("\ntesting geo_file, ", Test.Geo_file);

    // Construct the ParMooN Databases.
    ParameterDatabase parmoon_db =
        ParameterDatabase::parmoon_default_database();
    parmoon_db.merge(TDomain::default_domain_parameters());
    parmoon_db["boundary_file"].set(Test.Boundary_file, false);
    parmoon_db["geo_file"].set(Test.Geo_file, false);
    
    // Construct domain, thereby read in controls from the input file.
    TDomain domain(parmoon_db);
    
    domain.RegRefineAll();

    auto coll = domain.GetCollection(It_Finest, 0);
    check_collection(*coll, Test.has_simplex_cell);
    delete coll;
  }
  parmoon::parmoon_finalize();
}
