#include <BlockFEMatrix.h>
#include <FEMatrix.h>


#include "Domain.h"
#include "Database.h"
#include <ParMooN_repository_info.h>
#include <MainUtilities.h>
#include <Structure.h>
#include <array>
#include <map>
#include <string>
#include <templateNames.h>
#include "MooNMD_Io.h"
#include "ParMooN.h"
#ifdef __3D__
constexpr int d = 3;
#else
constexpr int d = 2;
#endif
using FESpace = typename Template_names<d>::FESpace;
using BoundCondFunct = typename Template_names<d>::BoundaryConditionFunction;


/* ########################################################################## */
/* Declare functions and constant expressions, see below for definition */
/* ########################################################################## */
void test_all(bool testall, ParameterDatabase parmoon_db);
void test_all_values(TStructure structure, int additional_face_integrals,
    int problem_number, int order, int other_order, bool has_hanging_nodes);
void test_constructors(bool testall, ParameterDatabase parmoon_db);
  std::array<std::string, 6> list_of_checked_values = {
  "Number hanging rows", "Number rows","Number columns",
  "Number entries", "Number entries in row", "Number hanging entries"};
std::array<unsigned int, 4> get_reference_values(int
    has_additional_face_integrals, int problem_number, int order,
    int other_order, bool has_hanging_nodes);
void test_default_constructor();
void test_remaining_constructors();
void test_is_square();
void test_equalequal();
void test_fortran_shift();
void test_transposed();
void test_multiplication();
void print_error(std::string error_string, int has_additional_face_integrals,
    int problem_number, int order, int other_order, bool has_hanging_nodes);
void indicator(double x, double y, double* values);
void indicator(double x, double y, double z, double* values);
const std::string path_to_repo = parmoon::source_directory;
const std::string path_to_meshes = path_to_repo + "/data/mesh/";


/* ########################################################################## */
/* Main program */
/* ########################################################################## */
int main(int, char* argv[])
{
  ParameterDatabase parmoon_db = parmoon::parmoon_initialize();
  bool testall = false;
  if (argv[1])
  {
    testall = (std::string(argv[1]).compare("testall") == 0);
  }
  parmoon_db.merge(TDomain::default_domain_parameters());

  Output::setVerbosity(3);    // set verbosity [ 1,5 ] higher for more output
  test_all(testall, parmoon_db);   // do the actual test
  Output::print(d, "D test passed.");
  parmoon::parmoon_finalize();
}


/* ########################################################################## */
/* Declaration of functions */
/* ########################################################################## */
void test_all(bool testall, ParameterDatabase parmoon_db)
{ // executes all the tests

  /* Check all constructors */
  test_constructors(testall, parmoon_db);

  /* Test is square */
  test_is_square();

  /* Test == */
  test_equalequal();

  /* Test Fortran shift */
  test_fortran_shift();

  /* Check transposed structure and multiplication */
  test_transposed();

  /* Check structure after several multiplications */
  test_multiplication();
} // end test_all


void test_constructors(bool testall, ParameterDatabase parmoon_db)
{ // tests all the constructors
  std::vector<short int> has_additional_face_integrals = {0, 1};
  std::vector<std::string> geometries;   // list of tested geometries
  std::vector<std::string> boundary_files; // corresponding boundary
  if (d == 2)
  {
    geometries ={"UnitSquare_quad.mesh", "UnitSquareIrregular.GEO","disk.mesh",
      "geothermal2d.mesh", "backward_facing_step/backward_facing_step_quad1."
        "mesh", "backward_facing_step/backward_facing_step_tria1.mesh",
      "doubleRiverbed.mesh", "flow_around_cylinder2D/flow_around_cylinder_"
        "quad1.mesh", "flow_around_cylinder2D/flow_around_cylinder_tria6.mesh",
      "Hemker_quad.mesh", "Hemker_tria.mesh", "ring.mesh", "Tshape_meshCodina_"
        "coarse.mesh", "UnitSquareCrissCross.GEO", "UnitSquareWithSpherical"
        "InscribedRegion.mesh"};
    boundary_files = {"UnitSquare.PRM", "UnitSquare.PRM", "disk.PRM",
      "geothermal2d.PRM", "backward_facing_step/backward_facing_step.PRM",
      "backward_facing_step/backward_facing_step.PRM", "doubleRiverbed.PRM",
      "flow_around_cylinder2D/flow_around_cylinder.PRM",
      "flow_around_cylinder2D/flow_around_cylinder.PRM", "Hemker.PRM",
      "Hemker.PRM","ring.PRM","Tshape.PRM", "UnitSquare.PRM","UnitSquare.PRM"};}
  else
  {
    geometries = {"Default_UnitCube_Hexa", "Default_UnitCube_Tetra",
      "channel.3d.mesh", "cylinder.3d.3K.mesh"};
    boundary_files = {"Default_UnitCube", "Default_UnitCube", "", ""};
  }
  if (!testall)
  {
    geometries.resize(4);
    boundary_files.resize(4);
  }

  // check everything twice, once with additional face integrals, once without
  for (int additional_face_integrals : has_additional_face_integrals)
  {
    /* Construct domain with various boundary conditions and hanging nodes*/
    for (unsigned int problem_number = 0; problem_number < geometries.size();
        ++problem_number)
    { // try all problems
      std::string geo_file;
      std::string boundary_file;
      if (d == 3)
      { // default geometries are directly coded in ParMooN
        if (problem_number < 2)
        {
          geo_file = geometries[problem_number];
        }
        else
        {
          geo_file = path_to_meshes + geometries[problem_number];
        }
        boundary_file = boundary_files[problem_number];
      }
      else
      {
        geo_file = path_to_meshes + geometries[problem_number];
        boundary_file = path_to_meshes+boundary_files[problem_number];
      }
      parmoon_db["geo_file"].set(geo_file, false);
      parmoon_db["boundary_file"].set(boundary_file, false);
      // construct domain, FE spaces and structure from database
      TDomain domain(parmoon_db);
      // Boundary conditions
      BoundCondFunct* BoundaryCondition;
      std::string boundcond;  // used for printing information
      if (problem_number % 2 == 0)
      { // half problems with Neumann boundary
        BoundaryCondition = BoundConditionNoBoundCondition;
        boundcond = "Neumann";
      }
      else
      { // half with Dirichlet boundary
#ifdef __2D__
        BoundaryCondition = BoundConditionNSE;
#else
        BoundaryCondition = BoundaryConditionNewton;
#endif
        boundcond = "Dirichlet";
      }

      auto coll = domain.GetCollection(It_Finest, 0); // collect particular mesh

      Output::print("\n\n---------- NEW TEST ----------\n\nTest with following",
          " properties:",
          "\nDimension:                 ", d,
          "\nProblem:                   ", geometries[problem_number]," (",
          problem_number,")",
          "\nBoundary Condition:        ", boundcond,
          "\nAdditional face integrals: ", additional_face_integrals,'\n');

      // Introduce Hanging nodes for first four geometries and try again. Due to
      // time only first four geometries are tested.
      parmoon_db["conforming_closure"] = false; // introduce hanging nodes
      TDomain domain_hn(parmoon_db);
      TCollection* coll_hn = nullptr;
      bool test_hanging_nodes = (problem_number < 4 && d != 3);
      if (test_hanging_nodes)
      { // refine the domain
        unsigned int max_refinements = 3; // maximal number of refinements
        for (unsigned int refinement_level = 0; refinement_level <
            max_refinements; ++refinement_level)
        {
          domain_hn.RefineByIndicator(indicator);
        }
        coll_hn = domain_hn.GetCollection(It_Finest,0);
        domain_hn.PS("hn_domain.ps", coll_hn);
      }

      std::vector<short int> list_order_spaces = {1, 2, -1, -12}; // see src/FE/FESpace2D.C
      if (!testall)
      {
        list_order_spaces = {1, -1};
      }
      for (auto order : list_order_spaces)
      {   // try all orders
        std::shared_ptr<FESpace> fespace (new FESpace(coll, "Test Space",
              BoundaryCondition, order)); // FE-space
        // structure that will be tested
        TStructure structure(fespace, false, additional_face_integrals);

        /* Test all values */
        test_all_values(structure, additional_face_integrals, problem_number,
            order, order, false);

        if(test_hanging_nodes && order != -12 && order != -1)
        { // order -12 (2D/3D) and -1(3D) are not implemented in the descriptor
          std::shared_ptr<FESpace> fespace_hn(new FESpace(
                coll_hn, "Test Space", BoundaryCondition, order)); // FE-space
          if(!additional_face_integrals)
          {
            TStructure structure_hn(fespace_hn, false,
                                    additional_face_integrals);
            test_all_values(structure_hn, additional_face_integrals,
                            problem_number, order, order, true);
          }
        }

        /* Test constructor for different FESpaces */
        if (problem_number < 4)
        { // check only the first 4 geometries due to time
          for (auto other_order : list_order_spaces)
          { // try all combinations
            if (other_order != order)
            {
              std::shared_ptr<FESpace> other_fespace (new FESpace(coll,
                    "Test Space", BoundaryCondition, other_order)); // FE-space
              TStructure structure2 = TStructure(fespace, other_fespace, false,
                                                 additional_face_integrals);

              /* Test all values */
              test_all_values(structure2, additional_face_integrals,
                  problem_number, order, other_order, false);

            } // endif other_order!=order
          } // end loop over other_order
        } // endif additional_face_integrals == 0
      } // end loop over orders

      delete coll;
      if (problem_number < 4)
      {
        delete coll_hn;
      }
    } // end loop over problems
  } // end loop over has_additional_face_integrals

  /* Check default constructor */
  test_default_constructor();

  /* Check (nRows, nCols)-constructor and (nRows, nCols, nActive, N_entries,
   *col_ptr, *row_ptr)-constructor */
  test_remaining_constructors();

} // end test_constructors

void test_all_values(TStructure structure, int additional_face_integrals,
    int problem_number, int order, int other_order, bool has_hanging_nodes)
{ // test all values for constructor with two different FE spaces
  long unsigned int test_values[] = {
    structure.get_n_rows(), structure.get_n_columns(),
    structure.get_n_entries(),structure.get_n_entries_in_row(3)};
  auto reference_values = get_reference_values(additional_face_integrals,
      problem_number, order, other_order, has_hanging_nodes);
  std::cout <<
    "\n\nDimension:                 "<< d<<
    "\nAdditional face integrals: "<< additional_face_integrals<<
    "\nHanging Nodes:             "<< has_hanging_nodes<<
    "\nProblem:                   "<< problem_number<<
    "\nOrder:                     "<< order<<
    "\nOther Order:               "<< other_order<<
    '\n';
  Output::print("\nrows,", " columns, ", "entries, ", "entries in row, ");
  for (unsigned int j = 0; j < reference_values.size()-1; j++)
  {
    std::cout << test_values[j] << ", ";
  }
  std::cout << test_values[reference_values.size()-1]<<'\n' << std::endl;
  for (unsigned int j = 0; j < reference_values.size(); j++)
  {
    if (test_values[j] != reference_values[j])
    {
      Output::print(list_of_checked_values[j], ":\nNow: ", test_values[j],
          ", expected: ", reference_values[j]);
      std::string error_string = list_of_checked_values[j] + " has changed.";
      print_error(error_string, additional_face_integrals, problem_number,
          order, other_order, has_hanging_nodes);
    }
  }
} // end test_all_values

void test_is_square()
{ // check is_square()
  TStructure structure_square = TStructure(3);
  if (!structure_square.is_square())
  {
    ErrThrow("is_square() of a square structure did not return true.");
  }
  TStructure structure_not_square = TStructure(3,4);
  if (structure_not_square.is_square())
  {
    ErrThrow("is_square() of a non-square structure returned true.");
  }
} // end check is_square

void test_equalequal()
{ //check == operator
  TStructure structure = TStructure(3, 4);
  TStructure structure2 = TStructure(4, 4);
  if (structure != structure)
  {
    ErrThrow("structure == structure does not return true.");
  }
  if (structure == structure2)
  {
    ErrThrow("structure == structure2 does not return false.");
  }
} // end check == operator

void test_fortran_shift()
{ // check fortran_shift()
  TStructure structure = TStructure(3, 4);
  if (structure.is_fortran_shifted())
  {   // check Fortran shift
    ErrThrow("Structure is Fortran shifted but should NOT be.");
  }
  else
  {
    structure.fortran_shift();
    if (!structure.is_fortran_shifted())
    {
      ErrThrow("Structure is NOT Fortran shifted but should be.");
    }
    structure.fortran_shift();
    if (structure.is_fortran_shifted())
    {
      ErrThrow("Structure is Fortran shifted but should NOT be.");
    }
  }
} // end test_fortran_shift

void test_default_constructor()
{
  TStructure structure0 = TStructure();
  std::array<unsigned int, 4> test_values_0 = { structure0.get_n_rows(),
    structure0.get_n_columns(), structure0.get_n_entries(), 0};
  for (unsigned int j = 0; j < test_values_0.size(); j++)
  {
    if (test_values_0[j] != 0)
    {
      Output::print("Now: ", test_values_0[j], ", should be: ", 0);
      ErrThrow("Default constructor TStructure() did not set ",
          list_of_checked_values[j], " to 0.");
    }
  }
} // end test_default_constructor

void test_remaining_constructors()
{ // test (nRows, nCols)-constructor and (nRows, nCols, nActive, N_entries,
  // *col_ptr, *row_ptr)-constructor
  TStructure structure_5x12 = TStructure(5, 12);
  std::array<unsigned int, 4> test_values_5x12 = {
    structure_5x12.get_n_rows(), structure_5x12.get_n_columns(),
    structure_5x12.get_n_entries(), 0};
  unsigned int reference_values_5x12[] = {5, 12, 0, 0};
  for (unsigned int j = 0; j < test_values_5x12.size(); j++)
  {
    if (test_values_5x12[j] != reference_values_5x12[j])
    {
      Output::print("Now: ", test_values_5x12[j], ", should be: ",
          reference_values_5x12[j]);
      ErrThrow("Constructor TStructure(5, 12) did not set ",
          list_of_checked_values[j], " to ", reference_values_5x12[j] ,".");
    }
  }

  std::array<int, 10> col_vector = {1, 2, 3, 0, 1, 2, 0, 1, 2, 3};
  std::array<int, 5> row_vector = {0, 2, 3, 6, 10};
  TStructure structure = TStructure(4, 6, 10, &col_vector[0],&row_vector[0]);
  if (structure.get_index_of_entry(0,0) != -1)
  {
    ErrThrow("get_index_of_entry did not return -1 even though the tested entry"
        " is not in the sparsity pattern.");
  }
  if (structure.get_index_of_entry(0,1) != 0)
  {
    ErrThrow("get_index_of_entry did not return 0 even though the tested entry"
        " is the first entry in the column vector.");
  }
  auto test_row_ptr = structure.get_row_ptr();
  auto test_rows = structure.get_row_array();
  for (int j = 0; j < 5; j++)
  {
    if (test_row_ptr[j] != row_vector[j])
    {
      ErrThrow("TStructure(int nRows, int nCols, int nActive, ",
          "int N_entries, int *col_ptr, int *row_ptr) failed to construct the "
          "structure correctly (or get_row_ptr() returned wrong values).\n"
          "test_row_ptr[", j,"] = ", test_row_ptr[j],", row_vector[",j,"] = ",
          row_vector[j],".");
    }
    if (test_rows[j] != row_vector[j])
    {
      ErrThrow("TStructure(int nRows, int nCols, int nActive, ",
          "int N_entries, int *col_ptr, int *row_ptr) failed to construct the "
          "structure correctly (or get_row_array() returned wrong values)."
          "\ntest_rows[", j,"] = ",test_rows[j], ", row_vector[",j,"] = ",
          col_vector[j],".");
    }
  }
  auto test_vector_columns = structure.get_vector_columns();
  auto test_columns = structure.get_columns();
  for (int j = 0; j < 10; j++)
  {
    if (test_vector_columns[j] != col_vector[j])
    {
      ErrThrow("TStructure(int nRows, int nCols, int nActive, ",
          "int N_entries, int *col_ptr, int *row_ptr) failed to construct the "
          "structure correctly (or get_vector_columns() returned wrong values)."
          "\ntest_vector_columns[", j,"] = ", test_vector_columns[j],
          ", col_vector[",j,"] = ", col_vector[j],".");
    }
    if (test_columns[j] != col_vector[j])
    {
      ErrThrow("TStructure(int nRows, int nCols, int nActive, ",
          "int N_entries, int *col_ptr, int *row_ptr) failed to construct the "
          "structure correctly (or get_columns() returned wrong values)."
          "\ntest_columns[", j,"] = ",test_columns[j], ", col_vector[",j,"] = ",
          col_vector[j],".");
    }
  }
} // end test_remaining_constructors

void test_transposed()
{
  std::array<int, 9> col_vector_3x4 = {1, 2, 3, 0, 2, 0, 1, 2, 3};
  std::array<int, 4> row_vector_3x4 = {0, 3, 5, 9};
  TStructure structure_3x4 = TStructure(3, 4, 9, &col_vector_3x4[0],
      &row_vector_3x4[0]);
  // transposed matrix
  std::shared_ptr<TStructure> structure_3x4T = structure_3x4.get_transposed();
  std::array<unsigned int, 4> test_values_3x4T = {
    structure_3x4T->get_n_rows(), structure_3x4T->get_n_columns(),
    structure_3x4T->get_n_entries(), 0};    // transposed values
  std::array<unsigned int, 4> reference_values_3x4T = {4, 3, 9, 0};
  for (int j = 0; j < 4; j++)
  {
    if (test_values_3x4T[j] != reference_values_3x4T[j])
    {
      ErrThrow("get_transposed() returned wrong structure. The value for \"",
          list_of_checked_values[j],"\" is incorrect.\ntest_values_3x4T[", j,
          "] = ", test_values_3x4T[j], ", reference_values_3x4T[", j,"] = ",
          reference_values_3x4T[j]);
    }
  }
  auto test_columns_3x4T = structure_3x4T->get_columns();
  std::array<int, 9> reference_cols_3x4T = {1, 2, 0, 2, 0, 1, 2, 0, 2};
  for (unsigned int j = 0; j < test_columns_3x4T.size(); j++)
  {
    if (test_columns_3x4T[j] != reference_cols_3x4T[j])
    {
      Output::print(test_columns_3x4T[j], ", ", reference_cols_3x4T[j]);
      ErrThrow("get_transposed() returned wrong structure. The column vector ",
          "is incorrect at position ",j,".");
    }
  }
  auto test_rows_3x4T = structure_3x4T->get_row_array();
  std::array<int, 5> reference_rows_3x4T = {0, 2, 4, 7, 9};
  for (unsigned int j = 0; j < test_rows_3x4T.size(); j++)
  {
    if (test_rows_3x4T[j] != reference_rows_3x4T[j])
    {
      Output::print(test_rows_3x4T[j], ", ", reference_rows_3x4T[j]);
      ErrThrow("get_transposed() returned wrong structure. The row vector ",
          "is incorrect at position ",j,".");
    }
  }
} // end test_transposed

void test_multiplication()
{ // Test the multiplication methods: get_product_structure and
  // get_structure_of_product_with_transpose_from_right
  std::array<int, 9> col_vector_3x4 = {1, 2, 3, 0, 2, 0, 1, 2, 3};
  std::array<int, 4> row_vector_3x4 = {0, 3, 5, 9};
  TStructure structure_3x4 = TStructure(3, 4, 9, &col_vector_3x4[0],
      &row_vector_3x4[0]);
  std::array<int, 3> col_vector_4x2 = {0, 1, 0};
  std::array<int, 5> row_vector_4x2 = {0, 1, 2, 2, 3};
  TStructure structure_4x2 = TStructure(4, 2, 3, &col_vector_4x2[0],
      &row_vector_4x2[0]);
  // Test multiplication
  auto structure_product = get_product_structure(structure_3x4, structure_4x2);
  std::array<unsigned int, 4> test_values_product = {
    structure_product->get_n_rows(), structure_product->get_n_columns(),
    structure_product->get_n_entries(), 0};    // transposed values
  std::array<unsigned int, 4> reference_values_product = {3, 2, 5, 0};
  for (int j = 0; j < 4; j++)
  {
    if (test_values_product[j] != reference_values_product[j])
    {
      ErrThrow("get_product_structure() returned wrong structure. The value for"
          "\"", list_of_checked_values[j],"\" is incorrect.\n",
          "test_values_product[", j, "] = ", test_values_product[j],
          ", reference_values_product[", j,"] = ", reference_values_product[j]);
    }
  }
  auto test_columns_product = structure_product->get_columns();
  std::array<int, 5> reference_cols_product = {0, 1, 0, 0, 1};
  for (unsigned int j = 0; j < test_columns_product.size(); j++)
  {
    if (test_columns_product[j] != reference_cols_product[j])
    {
      Output::print(test_columns_product[j], ", ", reference_cols_product[j]);
      ErrThrow("get_product_structure() returned wrong structure. The column ",
          "vector is incorrect at position ",j,".");
    }
  }
  auto test_rows_product = structure_product->get_row_array();
  std::array<int, 4> reference_rows_product = {0, 2, 3, 5};
  for (unsigned int j = 0; j < test_rows_product.size(); j++)
  {
    if (test_rows_product[j] != reference_rows_product[j])
    {
      Output::print(test_rows_product[j], ", ", reference_rows_product[j]);
      ErrThrow("get_product_structure() returned wrong structure. The row ",
          "vector is incorrect at position ",j,".");
    }
  }

  std::shared_ptr<TStructure> structure_3x4T = structure_3x4.get_transposed();
  auto product_structure =get_product_structure(structure_3x4, *structure_3x4T);
  auto product_with_transposed_structure =
    structure_3x4.get_structure_of_product_with_transpose_from_right();
  if (*product_structure != *product_with_transposed_structure)
  {
    ErrThrow("structure3x4.get_structure_of_product_with_transpose_from_right()"
        "is NOT equal to get_product_structure(structure_3x4, structure_3x4T)")
  }

  std::array<int, 6> col_vector_4x4 = {0, 3, 2, 0, 1, 2};
  std::array<int, 5> row_vector_4x4 = {0, 2, 3, 4, 6};
  TStructure structure_4x4 = TStructure(4, 4, 6, &col_vector_4x4[0],
      &row_vector_4x4[0]);
  auto double_product_structure = get_product_structure(
      *get_product_structure(structure_3x4,structure_4x4), *structure_3x4T);
  if (
      *double_product_structure !=
      *structure_3x4.get_structure_of_product_with_transpose_from_right(
        structure_4x4)
     )
  {
    ErrThrow("structure3x4.get_structure_of_product_with_transpose_from_right(",
        "structure_4x4) is NOT equal to get_product_structure(get_product_",
        "structure(*get_product_structure(structure_3x4,structure_4x4), ",
        "*structure_3x4T), structure_3x4T)")
  }
} // end test_multiplication

void print_error(std::string error_string, int additional_face_integrals,
    int problem_number, int order, int other_order, bool has_hanging_nodes)
{ // Print some error information
  ErrThrow(error_string,
      "\nDimension:                 ", d,
      "\nAdditional face integrals: ", additional_face_integrals,
      "\nHanging Nodes              ", has_hanging_nodes,
      "\nProblem:                   ", problem_number,
      "\nOrder:                     ", order,
      "\nOther order:               ", other_order
      );
} // end print_error

void indicator(double x, double y, double* values)
{ // Indicator functions for adaptive refinement
  indicator(x, y, 0, values);
}
void indicator(double x, double y, double z, double* values)
{
  values[0] = -0.5*x*x + y + 2*z + 0.5;
}


/* ########################################################################## */
/* Reference values */
/* ########################################################################## */
std::array<unsigned int, 4> get_reference_values(int
    has_additional_face_integrals, int problem_number, int order,
    int other_order, bool has_hanging_nodes)
{ // Reference values to check weather something has changed. This may
  // does not check corner cases.
  std::array<unsigned int, 4> reference_values;    // return value
  // Order of return values:
  // Number active rows, number hanging rows, number active entries,number rows,
  // number columns, number entries, number entries in row 3, number hanging
  // entries
  if (has_hanging_nodes == false)
  {
    if (d == 2)
    { // 2D test
      switch (has_additional_face_integrals)
      {
        case 0: // no additional face integrals
          switch (problem_number)
          { // different geometries
            case 0: // UnitSquare_quad.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {4, 4, 16, 4};
                      break;
                    case 2:
                      reference_values = {4, 9, 36, 9};
                      break;
                    case -1:
                      reference_values = {4, 4, 16, 4};
                      break;
                    case -12:
                      reference_values = {4, 9, 36, 9};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {9, 4, 36, 4};
                      break;
                    case 2:
                      reference_values = {9, 9, 81, 9};
                      break;
                    case -1:
                      reference_values = {9, 4, 36, 4};
                      break;
                    case -12:
                      reference_values = {9, 9, 81, 9};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {4, 4, 16, 4};
                      break;
                    case 2:
                      reference_values = {4, 9, 36, 9};
                      break;
                    case -1:
                      reference_values = {4, 4, 16, 4};
                      break;
                    case -12:
                      reference_values = {4, 9, 36, 9};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {9, 4, 36, 4};
                      break;
                    case 2:
                      reference_values = {9, 9, 81, 9};
                      break;
                    case -1:
                      reference_values = {9, 4, 36, 4};
                      break;
                    case -12:
                      reference_values = {9, 9, 81, 9};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 1: // UnitSquareIrregular.GEO
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {8, 8, 40, 5};
                      break;
                    case 2:
                      reference_values = {8, 24, 99, 12};
                      break;
                    case -1:
                      reference_values = {8, 16, 59, 7};
                      break;
                    case -12:
                      reference_values = {8, 54, 162, 18};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {24, 8, 99, 4};
                      break;
                    case 2:
                      reference_values = {24, 24, 228, 9};
                      break;
                    case -1:
                      reference_values = {24, 16, 129, 5};
                      break;
                    case -12:
                      reference_values = {24, 54, 324, 12};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {16, 8, 59, 4};
                      break;
                    case 2:
                      reference_values = {16, 24, 129, 9};
                      break;
                    case -1:
                      reference_values = {16, 16, 70, 5};
                      break;
                    case -12:
                      reference_values = {16, 54, 162, 12};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {54, 8, 162, 3};
                      break;
                    case 2:
                      reference_values = {54, 24, 324, 6};
                      break;
                    case -1:
                      reference_values = {54, 16, 162, 3};
                      break;
                    case -12:
                      reference_values = {54, 54, 324, 6};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 2: // disk.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {102, 102, 652, 7};
                      break;
                    case 2:
                      reference_values = {102, 377, 1724, 19};
                      break;
                    case -1:
                      reference_values = {102, 275, 1072, 12};
                      break;
                    case -12:
                      reference_values = {102, 1044, 3132, 36};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {377, 102, 1724, 4};
                      break;
                    case 2:
                      reference_values = {377, 377, 4115, 9};
                      break;
                    case -1:
                      reference_values = {377, 275, 2391, 5};
                      break;
                    case -12:
                      reference_values = {377, 1044, 6264, 12};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {275, 102, 1072, 4};
                      break;
                    case 2:
                      reference_values = {275, 377, 2391, 9};
                      break;
                    case -1:
                      reference_values = {275, 275, 1319, 5};
                      break;
                    case -12:
                      reference_values = {275, 1044, 3132, 12};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {1044, 102, 3132, 3};
                      break;
                    case 2:
                      reference_values = {1044, 377, 6264, 6};
                      break;
                    case -1:
                      reference_values = {1044, 275, 3132, 3};
                      break;
                    case -12:
                      reference_values = {1044, 1044, 6264, 6};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 3: // geothermal2d.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {17114,17114,119472,5};
                      break;
                    case 2:
                      reference_values = {17114,68293,324028,13};
                      break;
                    case -1:
                      reference_values = {17114,51179,204556,8};
                      break;
                    case -12:
                      reference_values = {17114, 204396, 613188, 24};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {68293,17114,324028,4};
                      break;
                    case 2:
                      reference_values = {68293,68293,784159,9};
                      break;
                    case -1:
                      reference_values = {68293,51179,460131,5};
                      break;
                    case -12:
                      reference_values = {68293, 204396, 1226376, 12};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {51179,17114,204556,4};
                      break;
                    case 2:
                      reference_values = {51179,68293,460131,9};
                      break;
                    case -1:
                      reference_values = {51179,51179,255575,5};
                      break;
                    case -12:
                      reference_values = {51179,204396,613188,12};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {204396, 17114, 613188, 3};
                      break;
                    case 2:
                      reference_values = {204396, 68293, 1226376, 6};
                      break;
                    case -1:
                      reference_values = {204396, 51179, 613188, 3};
                      break;
                    case -12:
                      reference_values = {204396, 204396, 1226376, 6};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 4: // backward_facing_step_quad1.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {128, 128, 882, 8};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 2:
                      reference_values = {423, 423, 5685, 12};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case -1:
                      reference_values = {213, 213, 1221, 6};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case -12:
                      reference_values = {762, 762, 6786, 6};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 5: // backward_facing_step_tria1.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {126, 126, 708, 7};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 2:
                      reference_values = {417, 417, 4155, 9};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case -1:
                      reference_values = {291, 291, 1287, 5};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case -12:
                      reference_values = {996, 996, 5976, 6};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 6: // dobuleRiverbed.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {59, 59, 375, 7};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 2:
                      reference_values = {217, 217, 2365, 9};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case -1:
                      reference_values = {158, 158, 758, 5};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case -12:
                      reference_values = {600, 600, 3600, 6};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 7: // flow_around_cylinder_quad1.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {60, 60, 444, 9};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 2:
                      reference_values = {208, 208, 2944, 15};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case -1:
                      reference_values = {104, 104, 632, 7};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case -12:
                      reference_values = {396, 396, 3564, 9};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 8: // flow_around_cylinder_tria6.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {899, 899, 5955, 7};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 2:
                      reference_values = {3427, 3427,38143,9};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case -1:
                      reference_values = {2528, 2528,12302,5};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case -12:
                      reference_values = {9774, 9774,58644,6};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 9: // Hemker_quad.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {70, 70, 522, 9};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 2:
                      reference_values = {244, 244, 3472, 25};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case -1:
                      reference_values = {122, 122, 746, 7};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case -12:
                      reference_values = {468, 468, 4212, 9};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 10:  // Hemker_tria.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {70, 70, 418, 8};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 2:
                      reference_values = {244, 244, 2536, 9};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case -1:
                      reference_values = {174, 174, 798, 5};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case -12:
                      reference_values = {624, 624, 3744, 6};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 11:  // ring.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {248, 248, 1664, 8};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 2:
                      reference_values = {956, 956, 10724, 9};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case -1:
                      reference_values = {708, 708, 3468, 5};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case -12:
                      reference_values = {2760, 2760,16560,6};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 12:  // Tshape_meshCodina_coarse.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {407, 407, 2659, 7};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 2:
                      reference_values = {1533, 1533,16929,9};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case -1:
                      reference_values = {1126, 1126, 5446, 5};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case -12:
                      reference_values = {4320, 4320,25920,6};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 13:  // UnitSquareCrissCross.GEO
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {5, 5, 21, 4};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 2:
                      reference_values = {13, 13, 109, 9};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case -1:
                      reference_values = {8, 8, 32, 5};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case -12:
                      reference_values = {24, 24, 144, 6};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 14:  // UnitSquareWithSphericalInscribedRegion.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {190, 190, 1244, 8};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 2:
                      reference_values = {717, 717, 7935, 9};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case -1:
                      reference_values = {527, 527, 2555, 5};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case -12:
                      reference_values = {2028, 2028,12168,6};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
          } // end cases for geometries
          break;
        case 1: // with additional face integrals
          switch (problem_number)
          { // different geometries
            case 0: // UnitSquare_quad.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {4, 4, 16, 4};
                      break;
                    case 2:
                      reference_values = {4, 9, 36, 9};
                      break;
                    case -1:
                      reference_values = {4, 4, 16, 4};
                      break;
                    case -12:
                      reference_values = {4, 9, 36, 9};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {9, 4, 36, 4};
                      break;
                    case 2:
                      reference_values = {9, 9, 81, 9};
                      break;
                    case -1:
                      reference_values = {9, 4, 36, 4};
                      break;
                    case -12:
                      reference_values = {9, 9, 81, 9};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {4, 4, 16, 4};
                      break;
                    case 2:
                      reference_values = {4, 9, 36, 9};
                      break;
                    case -1:
                      reference_values = {4, 4, 16, 4};
                      break;
                    case -12:
                      reference_values = {4, 9, 36, 9};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {9, 4, 36, 4};
                      break;
                    case 2:
                      reference_values = {9, 9, 81, 9};
                      break;
                    case -1:
                      reference_values = {9, 4, 36, 4};
                      break;
                    case -12:
                      reference_values = {9, 9, 81, 9};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 1: // UnitSquareIrregular.GEO
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {8, 8, 58, 7};
                      break;
                    case 2:
                      reference_values = {8, 24, 157, 19};
                      break;
                    case -1:
                      reference_values = {8, 16, 99, 12};
                      break;
                    case -12:
                      reference_values = {8, 54, 294, 36};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {24, 8, 157, 7};
                      break;
                    case 2:
                      reference_values = {24, 24, 410, 18};
                      break;
                    case -1:
                      reference_values = {24, 16, 253, 11};
                      break;
                    case -12:
                      reference_values = {24, 54, 720, 30};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {16, 8, 99, 7};
                      break;
                    case 2:
                      reference_values = {16, 24, 253, 19};
                      break;
                    case -1:
                      reference_values = {16, 16, 154, 12};
                      break;
                    case -12:
                      reference_values = {16, 54, 426, 36};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {54, 8, 294, 5};
                      break;
                    case 2:
                      reference_values = {54, 24, 720, 12};
                      break;
                    case -1:
                      reference_values = {54, 16, 426, 7};
                      break;
                    case -12:
                      reference_values = {54, 54, 1116, 18};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 2: // disk.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {102, 102, 1146, 13};
                      break;
                    case 2:
                      reference_values = {102, 377, 3206, 37};
                      break;
                    case -1:
                      reference_values = {102, 275, 2060, 24};
                      break;
                    case -12:
                      reference_values = {102, 1044, 6096, 72};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {377, 102, 3206, 8};
                      break;
                    case 2:
                      reference_values = {377, 377, 8561, 21};
                      break;
                    case -1:
                      reference_values = {377, 275, 5355, 13};
                      break;
                    case -12:
                      reference_values = {377, 1044, 15156,36};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {275, 102, 2060, 8};
                      break;
                    case 2:
                      reference_values = {275, 377, 5355, 21};
                      break;
                    case -1:
                      reference_values = {275, 275, 3295, 13};
                      break;
                    case -12:
                      reference_values = {275, 1044, 9060, 36};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {1044, 102, 6096, 6};
                      break;
                    case 2:
                      reference_values = {1044, 377,15156,15};
                      break;
                    case -1:
                      reference_values = {1044, 275, 9060, 9};
                      break;
                    case -12:
                      reference_values = {1044,1044,24048,24};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 3: // geothermal2d.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {17114,17114,220778,9};
                      break;
                    case 2:
                      reference_values = {17114,68293,628678,25};
                      break;
                    case -1:
                      reference_values = {17114,51179,407900,16};
                      break;
                    case -12:
                      reference_values = {17114, 204396, 1225416, 48};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {68293,17114,628678,7};
                      break;
                    case 2:
                      reference_values = {68293, 68293, 1699573, 19};
                      break;
                    case -1:
                      reference_values = {68293, 51179, 1070895, 12};
                      break;
                    case -12:
                      reference_values = {68293, 204396, 3063060, 36};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {51179,17114,407900,7};
                      break;
                    case 2:
                      reference_values = {51179, 68293, 1070895, 19};
                      break;
                    case -1:
                      reference_values = {51179, 51179, 662995, 12};
                      break;
                    case -12:
                      reference_values = {51179, 204396, 1837644, 36};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {204396, 17114, 1225416, 6};
                      break;
                    case 2:
                      reference_values = {204396, 68293, 3063060, 15};
                      break;
                    case -1:
                      reference_values = {204396, 51179, 1837644, 9};
                      break;
                    case -12:
                      reference_values = {204396, 204396, 4899744, 24};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 4: // backward_facing_step_quad1.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {128, 128, 1652, 13};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 2:
                      reference_values = {423, 423, 13529, 22};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case -1:
                      reference_values = {213, 213, 3203, 13};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case -12:
                      reference_values = {762, 762, 26928, 24};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 5: // backward_facing_step_tria1.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {126, 126, 1114, 11};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 2:
                      reference_values = {417, 417, 7849, 18};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case -1:
                      reference_values = {291, 291, 2935, 11};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case -12:
                      reference_values = {996, 996, 20880, 18};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 6: // dobuleRiverbed.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {59, 59, 655, 12};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 2:
                      reference_values = {217, 217, 4905, 21};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case -1:
                      reference_values = {158, 158, 1890, 13};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case -12:
                      reference_values = {600, 600, 13824, 24};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 7: // flow_around_cylinder_quad1.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {60, 60, 860, 15};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 2:
                      reference_values = {208, 208, 7320, 31};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case -1:
                      reference_values = {104, 104, 1726, 18};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case -12:
                      reference_values = {396, 396, 15228, 36};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 8: // flow_around_cylinder_tria6.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {899, 899, 10633, 13};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 2:
                      reference_values = {3427,3427,80445,21};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case -1:
                      reference_values = {2528,2528,31134,13};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case -12:
                      reference_values = {9774,9774,228492,24};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 9: // Hemker_quad.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {70, 70, 1090, 19};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 2:
                      reference_values = {244, 244, 9184, 45};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case -1:
                      reference_values = {122, 122, 2174, 20};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case -12:
                      reference_values = {468, 468, 18144, 27};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 10:  // Hemker_tria.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {70, 70, 694, 14};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 2:
                      reference_values = {244, 244, 5020, 18};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case -1:
                      reference_values = {174, 174, 1902, 11};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case -12:
                      reference_values = {624, 624, 13680, 18};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 11:  // ring.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {248, 248, 2976, 13};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 2:
                      reference_values = {956, 956, 22692, 18};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case -1:
                      reference_values = {708, 708, 8812, 13};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case -12:
                      reference_values = {2760,2760,64944,24};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 12:  // Tshape_meshCodina_coarse.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {407, 407, 4727, 13};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 2:
                      reference_values = {1533,1533,35541,21};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case -1:
                      reference_values = {1126,1126,13718,13};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case -12:
                      reference_values = {4320,4320,100368,18};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 13:  // UnitSquareCrissCross.GEO
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {5, 5, 25, 5};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 2:
                      reference_values = {13, 13, 165, 13};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case -1:
                      reference_values = {8, 8, 60, 8};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case -12:
                      reference_values = {24, 24, 432, 18};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 14:  // UnitSquareWithSphericalInscribedRegion.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {190, 190, 2202, 15};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 2:
                      reference_values = {717, 717, 16637, 21};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case -1:
                      reference_values = {527, 527, 6435, 13};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case -12:
                      reference_values = {2028,2028,47232,24};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
          } // end cases for geometries
          break;
      } // end cases for face_integrals
    } // end 2D case
    else
    { // 3D test
      switch (has_additional_face_integrals)
      {
        case 0: // no additional face integrals
          switch (problem_number)
          { // different geometries
            case 0: // Default_UnitCube_Hexa
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {8, 8, 64, 8};
                      break;
                    case 2:
                      reference_values = {8, 27, 216, 27};
                      break;
                    case -1:
                      reference_values = {8, 6, 48, 6};
                      break;
                    case -12:
                      reference_values = {8, 27, 216, 27};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {27, 8, 216, 8};
                      break;
                    case 2:
                      reference_values = {27, 27, 729, 27};
                      break;
                    case -1:
                      reference_values = {27, 6, 162, 6};
                      break;
                    case -12:
                      reference_values = {27, 27, 729, 27};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {6, 8, 48, 8};
                      break;
                    case 2:
                      reference_values = {6, 27, 162, 27};
                      break;
                    case -1:
                      reference_values = {6, 6, 36, 6};
                      break;
                    case -12:
                      reference_values = {6, 27, 162, 27};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {27, 8, 216, 8};
                      break;
                    case 2:
                      reference_values = {27, 27, 729, 27};
                      break;
                    case -1:
                      reference_values = {27, 6, 162, 6};
                      break;
                    case -12:
                      reference_values = {27, 27, 729, 27};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 1: // Default_UnitCube_Tetra
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {8, 8, 46, 8};
                      break;
                    case 2:
                      reference_values = {8, 27, 138, 27};
                      break;
                    case -1:
                      reference_values = {8, 18, 78, 18};
                      break;
                    case -12:
                      reference_values = {8, 60, 240, 60};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {27, 8, 138, 5};
                      break;
                    case 2:
                      reference_values = {27, 27, 393, 14};
                      break;
                    case -1:
                      reference_values = {27, 18, 204, 7};
                      break;
                    case -12:
                      reference_values = {27, 60, 600, 20};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {18, 8, 78, 5};
                      break;
                    case 2:
                      reference_values = {18, 27, 204, 14};
                      break;
                    case -1:
                      reference_values = {18, 18, 90, 7};
                      break;
                    case -12:
                      reference_values = {18, 60, 240, 20};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {60, 8, 240, 4};
                      break;
                    case 2:
                      reference_values = {60, 27, 600, 10};
                      break;
                    case -1:
                      reference_values = {60, 18, 240, 4};
                      break;
                    case -12:
                      reference_values = {60, 60, 600, 10};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 2: // channel.3d.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {96, 96, 882, 6};
                      break;
                    case 2:
                      reference_values = {96, 489, 3183, 18};
                      break;
                    case -1:
                      reference_values = {96, 505, 2343, 10};
                      break;
                    case -12:
                      reference_values = {96, 2070, 8280, 30};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {489, 96, 3183, 6};
                      break;
                    case 2:
                      reference_values = {489, 489, 10149, 19};
                      break;
                    case -1:
                      reference_values = {489, 505, 6342, 12};
                      break;
                    case -12:
                      reference_values = {489, 2070, 20700,40};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {505, 96, 2343, 5};
                      break;
                    case 2:
                      reference_values = {505, 489, 6342, 14};
                      break;
                    case -1:
                      reference_values = {505, 505, 2989, 7};
                      break;
                    case -12:
                      reference_values = {505, 2070, 8280, 20};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {2070, 96, 8280, 4};
                      break;
                    case 2:
                      reference_values = {2070, 489,20700,10};
                      break;
                    case -1:
                      reference_values = {2070, 505, 8280, 4};
                      break;
                    case -12:
                      reference_values = {2070,2070,20700,10};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 3: // cylinder.3d.3K.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {3328,3328,45326,15};
                      break;
                    case 2:
                      reference_values = {3328,24327,190290,65};
                      break;
                    case -1:
                      reference_values = {3328,34322,169566,60};
                      break;
                    case -12:
                      reference_values = {3328,166500,666000,240};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {24327,3328,190290,16};
                      break;
                    case 2:
                      reference_values = {24327,24327,662085,70};
                      break;
                    case -1:
                      reference_values = {24327,34322,472332,65};
                      break;
                    case -12:
                      reference_values = {24327, 166500, 1665000, 260};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {34322,3328,169566,5};
                      break;
                    case 2:
                      reference_values = {34322,24327,472332,14};
                      break;
                    case -1:
                      reference_values = {34322,34322,234122,7};
                      break;
                    case -12:
                      reference_values = {34322,166500,666000,20};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {166500,3328,666000,4};
                      break;
                    case 2:
                      reference_values = {166500, 24327, 1665000, 10};
                      break;
                    case -1:
                      reference_values = {166500, 34322, 666000, 4};
                      break;
                    case -12:
                      reference_values = {166500, 166500, 1665000, 10};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
          } // end cases for geometries
          break;
        case 1: // with additional face integrals
          switch (problem_number)
          { // different geometries
            case 0: // Default_UnitCube_Hexa
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {8, 8, 64, 8};
                      break;
                    case 2:
                      reference_values = {8, 27, 216, 27};
                      break;
                    case -1:
                      reference_values = {8, 6, 48, 6};
                      break;
                    case -12:
                      reference_values = {8, 27, 216, 27};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {27, 8, 216, 8};
                      break;
                    case 2:
                      reference_values = {27, 27, 729, 27};
                      break;
                    case -1:
                      reference_values = {27, 6, 162, 6};
                      break;
                    case -12:
                      reference_values = {27, 27, 729, 27};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {6, 8, 48, 8};
                      break;
                    case 2:
                      reference_values = {6, 27, 162, 27};
                      break;
                    case -1:
                      reference_values = {6, 6, 36, 6};
                      break;
                    case -12:
                      reference_values = {6, 27, 162, 27};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {27, 8, 216, 8};
                      break;
                    case 2:
                      reference_values = {27, 27, 729, 27};
                      break;
                    case -1:
                      reference_values = {27, 6, 162, 6};
                      break;
                    case -12:
                      reference_values = {27, 27, 729, 27};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 1: // Default_UnitCube_Tetra
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {8, 8, 58, 8};
                      break;
                    case 2:
                      reference_values = {8, 27, 186, 27};
                      break;
                    case -1:
                      reference_values = {8, 18, 114, 18};
                      break;
                    case -12:
                      reference_values = {8, 60, 360, 60};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {27, 8, 186, 7};
                      break;
                    case 2:
                      reference_values = {27, 27, 585, 22};
                      break;
                    case -1:
                      reference_values = {27, 18, 348, 13};
                      break;
                    case -12:
                      reference_values = {27, 60, 1080, 40};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {18, 8, 114, 7};
                      break;
                    case 2:
                      reference_values = {18, 27, 348, 22};
                      break;
                    case -1:
                      reference_values = {18, 18, 198, 13};
                      break;
                    case -12:
                      reference_values = {18, 60, 600, 40};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {60, 8, 360, 6};
                      break;
                    case 2:
                      reference_values = {60, 27, 1080, 18};
                      break;
                    case -1:
                      reference_values = {60, 18, 600, 10};
                      break;
                    case -12:
                      reference_values = {60, 60, 1800, 30};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 2: // channel.3d.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {96, 96, 1280, 9};
                      break;
                    case 2:
                      reference_values = {96, 489, 4986, 30};
                      break;
                    case -1:
                      reference_values = {96, 505, 3957, 19};
                      break;
                    case -12:
                      reference_values = {96, 2070, 14350, 60};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {489, 96, 4986, 13};
                      break;
                    case 2:
                      reference_values = {489, 489, 18029, 47};
                      break;
                    case -1:
                      reference_values = {489, 505, 13083, 33};
                      break;
                    case -12:
                      reference_values = {489, 2070,45370,110};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {505, 96, 3957, 8};
                      break;
                    case 2:
                      reference_values = {505, 489, 13083, 26};
                      break;
                    case -1:
                      reference_values = {505, 505, 8401, 16};
                      break;
                    case -12:
                      reference_values = {505, 2070, 27270,50};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {2070, 96, 14350, 6};
                      break;
                    case 2:
                      reference_values = {2070, 489,45370,18};
                      break;
                    case -1:
                      reference_values = {2070, 505,27270,10};
                      break;
                    case -12:
                      reference_values = {2070,2070,85300,30};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
            case 3: // cylinder.3d.3K.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {3328,3328,86466,30};
                      break;
                    case 2:
                      reference_values = {3328,24327,378474,134};
                      break;
                    case -1:
                      reference_values = {3328,34322,340026,123};
                      break;
                    case -12:
                      reference_values = {3328, 166500, 1311560, 480};
                      break;
                  } // end cases for other order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {24327,3328,378474,30};
                      break;
                    case 2:
                      reference_values = {24327, 24327, 1485277, 136};
                      break;
                    case -1:
                      reference_values = {24327, 34322, 1177380, 127};
                      break;
                    case -12:
                      reference_values = {24327, 166500, 4247240, 500};
                      break;
                  } // end cases for other order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {34322,3328,340026,10};
                      break;
                    case 2:
                      reference_values = {34322, 24327, 1177380, 35};
                      break;
                    case -1:
                      reference_values = {34322, 34322, 791918, 24};
                      break;
                    case -12:
                      reference_values = {34322, 166500, 2602680, 80};
                      break;
                  } // end cases for other order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {166500, 3328, 1311560, 7};
                      break;
                    case 2:
                      reference_values = {166500, 24327, 4247240, 22};
                      break;
                    case -1:
                      reference_values = {166500, 34322, 2602680, 13};
                      break;
                    case -12:
                      reference_values = {166500, 166500, 8120600, 40};
                      break;
                  } // end cases for other order
                  break;
              } // end cases for order
              break;
          } // end cases for geometries
          break;
      } // end cases for face_integrals
    } // end 2D case
  } // end has_hanging_nodes == false
  else
  { // has_hanging_nodes == true
    if (d == 2)
    { // 2D case
      switch (has_additional_face_integrals)
      {
        case 0:
          switch (problem_number)
          {
            case 0: // UnitSquare_quad.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {19, 19, 131, 4};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values = {57, 57, 761, 37};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {32, 32, 152, 4};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
              } // endswitch order
              break;
            case 1: // UnitSquareIrregular.GEO
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {29, 29, 211, 9};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values = {93, 93, 1017, 19};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {73, 73, 289, 5};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
              } // endswitch order
              break;
            case 2: // disk.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {538, 538, 4462,7};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values ={1918,1918,22579,9};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {1570,1570,6628,5};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
              } // endswitch order
              break;
            case 3: // geothermal2d.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {18523, 18523, 131943, 5};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values = {73295, 73295, 844348, 9};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {55394, 55394, 272894, 5};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
              } // endswitch order
              break;
          } // endswitch problem_number
          break;
      } // endswitch has_additional_face_integrals
    } // endif d==2
    else
    { // 3D case
      switch (has_additional_face_integrals)
      {
        case 0:
          switch (problem_number)
          {
            case 0: // Default_UnitCube_Hexa
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {137, 137, 2503,27};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values = {771,771,36717,75};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
              } // endswitch order
              break;
            case 1: // Default_UnitCube_Tetra
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {62, 62, 614, 7};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values = {298, 298, 6144, 72};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
              } // endswitch order
              break;
            case 2: // channel.3d.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {2793, 2793, 40641, 38};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values = {17027, 17027, 435951, 23};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
              } // endswitch order
              break;
            case 3: // cylinder.3d.3K.mesh
              switch (order)
              {
                case 1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {26426, 26426, 351660, 13};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case 2:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values = {168191, 168191, 4344612, 27};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case -1:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
                case -12:
                  switch (other_order)
                  {
                    case 1:
                      reference_values = {};
                      break;
                    case 2:
                      reference_values = {};
                      break;
                    case -1:
                      reference_values = {};
                      break;
                    case -12:
                      reference_values = {};
                      break;
                  } // endswitch other_order
                  break;
              } // endswitch order
              break;
          } // endswitch problem_number
          break;

      } // endswitch has_additional_face_integrals
    } // endif d == 3
  } // end has_hanging_nodes == true
  return reference_values;
} // end check_all_values for constructors
