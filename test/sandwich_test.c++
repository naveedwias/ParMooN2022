#include <Domain.h>
#include <Database.h>
#include <BaseCell.h>
#include <ParMooN_repository_info.h>
#include <FESpace3D.h>
#include "ParMooN.h"

const std::string path = parmoon::source_directory;
const std::string path_to_repo = path + "/data/mesh/";

struct ReferenceData
{
  public:
    int N_Cells;
    int N_DOFs;
    int N_actDOFs;
};

struct TestObject
{
  public:
    std::string Name;
    std::string Geo_file;
    std::string Boundary_file;
    ReferenceData Unrefined_Mesh;
    ReferenceData Refined_Mesh;
};

TestObject prepare_test0()
{
  std::string name = "UnitSquare_tria";

  std::string geo_file = path_to_repo + "UnitSquare_tria.GEO";
  std::string bd_file = path_to_repo + "UnitSquare.PRM";

  auto unrefined = ReferenceData{12, 45, 3};
  auto refined = ReferenceData{96, 225, 63};

  return TestObject{name, geo_file, bd_file, unrefined, refined};
}
TestObject prepare_test1()
{
  std::string name = "UnitSquareCrissCross";

  std::string geo_file = path_to_repo + "UnitSquareCrissCross.GEO";
  std::string bd_file = path_to_repo + "UnitSquare.PRM";

  auto unrefined = ReferenceData{24, 65, 15};
  auto refined = ReferenceData{192, 369, 175};

  return TestObject{name, geo_file, bd_file, unrefined, refined};
}
TestObject prepare_test2()
{
  std::string name = "UnitSquare_quads";

  std::string geo_file = path_to_repo + "UnitSquare_quads.mesh";
  std::string bd_file = path_to_repo + "UnitSquare.PRM";

  auto unrefined = ReferenceData{8, 125, 27};
  auto refined = ReferenceData{64, 729, 343};

  return TestObject{name, geo_file, bd_file, unrefined, refined};
}
TestObject prepare_test3()
{
  std::string name = "UnitSquareIrregular";

  std::string geo_file = path_to_repo + "UnitSquareIrregular.GEO";
  std::string bd_file = path_to_repo + "UnitSquare.PRM";

  auto unrefined = ReferenceData{54, 120, 42};
  auto refined = ReferenceData{432, 747, 441};

  return TestObject{name, geo_file, bd_file, unrefined, refined};
}
TestObject prepare_test4()
{
  std::string name = "Hemker_tria";

  std::string geo_file = path_to_repo + "Hemker_tria.GEO";
  std::string bd_file = path_to_repo + "Hemker.PRM";

  auto unrefined = ReferenceData{624, 1220, 516};
  auto refined = ReferenceData{4992, 8136, 5320};

  return TestObject{name, geo_file, bd_file, unrefined, refined};
}
TestObject prepare_test5()
{
  std::string name = "Hemker_quad";

  std::string geo_file = path_to_repo + "Hemker_quad.mesh";
  std::string bd_file = path_to_repo + "Hemker.PRM";

  auto unrefined = ReferenceData{104, 1220, 516};
  auto refined = ReferenceData{832, 8136, 5320};

  return TestObject{name, geo_file, bd_file, unrefined, refined};
}

void all_dirichlet_boundary_condition(int, double, double, double,
                                      BoundCond & bc)
{
  bc = DIRICHLET;
}

void get_input_data(const std::string geo_file, const std::string boundary_file,
                    ParameterDatabase& db)
{
  db.merge(TDomain::default_domain_parameters(), true, true);
  auto nested_db = TDomain::default_sandwich_grid_parameters();
  nested_db["n_layers"].set(2);
  db.add_nested_database(nested_db);

  db["boundary_file"].set(boundary_file, false);
  db["geo_file"].set(geo_file, false);
  db["refinement_n_initial_steps"].set(1);
  db["sandwich_grid"] = true;

}

bool space_info(const ReferenceData ref_data, TCollection* gridCollections)
{
  auto space = TFESpace3D(gridCollections, "refined_mesh",
                          all_dirichlet_boundary_condition, 2);
  int n_cells = space.GetCollection()->GetN_Cells();
  Output::print("N_Cells    : ", setw(13), n_cells, ", while expected ",
                ref_data.N_Cells);
  int n_DOFS = space.get_n_dof();
  Output::print("dofs all   : ", setw(13), n_DOFS, ", while expected ",
                ref_data.N_DOFs);
  int n_actDOFS = space.get_n_active();
  Output::print("dof active : ", setw(13), n_actDOFS, ", while expected ",
                ref_data.N_actDOFs);

  if(n_cells != ref_data.N_Cells)
  {
    return false;
  }
  if(n_DOFS != ref_data.N_DOFs)
  {
    return false;
  }
  if(n_actDOFS != ref_data.N_actDOFs)
  {
    return false;
  }

  return true;
}

int main()
{
  parmoon::parmoon_initialize();
  std::vector<TestObject> TestObjects = {prepare_test0(), prepare_test1(),
                                         prepare_test2(), prepare_test3(),
                                         prepare_test4(), prepare_test5()};

  for(auto &Test : TestObjects)
  {
    Output::print(Test.Name);

    // Construct the ParMooN Databases.
    ParameterDatabase parmoon_db =
        ParameterDatabase::parmoon_default_database();
    get_input_data(Test.Geo_file, Test.Boundary_file, parmoon_db);

    // Construct domain, thereby read in controls from the input file.
    TDomain domain(parmoon_db);

    if(!domain.check())
      ErrThrow("The Domain constructed by MakeSandwichGrid is not consistent.");

    std::list<TCollection*> gridCollections;
    gridCollections.push_front(domain.GetCollection(It_Finest, 0));

    domain.RegRefineAll();
    gridCollections.push_front(domain.GetCollection(It_Finest, 0));

    if(!space_info(Test.Unrefined_Mesh, gridCollections.back()))
      ErrThrow(
          "The Domain constructed by MakeSandwichGrid() leads to a wrong FESpace");

    if(!space_info(Test.Refined_Mesh, gridCollections.front()))
      ErrThrow("The Domain after Refine() leads to a wrong FESpace");
  }
  parmoon::parmoon_finalize();
}
