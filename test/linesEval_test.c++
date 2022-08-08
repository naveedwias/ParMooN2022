#include <Domain.h>
#include <Database.h>
#ifdef __2D__
#include <FESpace2D.h>
#include <FEFunction2D.h>
#else // __3D__
#include <FESpace3D.h>
#include <FEFunction3D.h>
#endif
#include <LinesEval.h>
#include "MooNMD_Io.h"
#include "ParMooN.h"

#ifdef _MPI
#include <mpi.h>
#endif

/* ************************************************************************** */
void analytic_function(
#ifdef __2D__  
                       double x, double y, double *values)
{ values[0] = x + y; }
#else // __3D__
                       double x, double y, double z, double *values)
{ values[0] = x + y + z; }
#endif


/* ************************************************************************** */
void all_dirichlet_boundary_condition(int, double,
#ifdef __2D__
                                                   BoundCond & bc)
#else // __3D__
                                                   double, double,
                                      BoundCond & bc)
#endif
{ bc = DIRICHLET; }


/* ************************************************************************** */
template<int d>
bool additional_tests(TDomain & domain)
{
  auto lines_eval_db = LinesEval<d>::default_lineseval_parameters();
  lines_eval_db["line_direction"]  = 1;
  lines_eval_db["line_position"]   = std::vector<double>({0.,0.,0., 0.,0.5,0.});
  lines_eval_db["line_refinement"] = 1;

  try
  {
    lines_eval_db["line_position"] = std::vector<double>({0.,0.,0., 0.5,0.});
    LinesEval<d> lines(domain, lines_eval_db);
    Output::print("there should be an error with invalid 'line_position'");
    return false;
  }
  catch(...)
  {
    // ok, this should have been an error
    lines_eval_db["line_position"] = std::vector<double>({0.,0.,0., 0.,0.5,0.});
  }

  try
  {
    lines_eval_db["line_direction"] = 5;
    //LinesEval lines(domain, lines_eval_db);
    Output::print("there should be an error with invalid 'line_direction'");
    return false;
  }
  catch(...)
  {
    // ok, this should have been an error
    lines_eval_db["line_direction"] = 1;
  }

  return true;
}


/* ************************************************************************** */
int main(int argc, char* argv[])
{
  auto db = parmoon::parmoon_initialize(argc, argv);
  db.merge(TDomain::default_domain_parameters());

  bool root = true;
#ifdef _MPI
  db["refinement_n_initial_steps"] = 4;
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  root = (my_rank == 0);
#else
  db["refinement_n_initial_steps"] = 0;
#endif
  Output::setVerbosity(3);
  
#ifdef __2D__
  const int d = 2;
  db["boundary_file"] = "Default_UnitSquare";
  db["geo_file"]      = "TwoTriangles";
#else // __3D__
  const int d = 3;
  db["boundary_file"] = "Default_UnitCube";
  db["geo_file"]      = "Default_UnitCube_Tetra";
#endif

  using FESpace    = typename Template_names<d>::FESpace;
  using FEFunction = typename Template_names<d>::FEFunction;

  TDomain domain(db);
  
  // refinement
  domain.refine_and_get_hierarchy_of_collections(db);
  
  if(!additional_tests<d>(domain))
  {
    Output::print("additional_tests failed");
    return 1;
  }
  
  auto coll = domain.get_grid_collections().front();
  
  std::shared_ptr<const FESpace> fespace(
    new FESpace(coll, "dummy", all_dirichlet_boundary_condition, 2));
  int length = fespace->get_n_dof();
  std::vector<double> entries(length, 0);
  FEFunction fefunction(fespace, "testfunction", &entries[0]);
  fefunction.Interpolate(analytic_function);

  auto lines_eval_db = LinesEval<d>::default_lineseval_parameters();
  std::vector<double> posi = {0.0,0.0,0.0, 0.5,0.5,0.5};
//  lines_eval_db["position_file"];
  lines_eval_db["line_direction"]  = 1;
  lines_eval_db["line_position"]   = posi;
  lines_eval_db["line_refinement"] = 1;
  db.add_nested_database(lines_eval_db);

  LinesEval<d> lines(domain, lines_eval_db);
  LinesEval<d> lines_test(domain, db);

#ifdef _MPI
  double position[66];
  double space_avg[2] = {0.5, 1.5};
  int n_lines = 2;
  int n_points[2] = {33, 33};
  for(int i=0;i<33;i++)
  {
    position[i]    = i*0.03125;
    position[i+33] = i*0.03125;
  }
#else
  // expected values
#ifdef __2D__
  double position[11]  = {0., 0.5, 1., 0., 0.5, 1., 0., 0.25, 0.5, 0.75, 1.};
  double space_avg[3] = {0.5, 0.5, 1.};
  int n_lines = 3;
  int n_points[3] = {3, 3, 5};
#else // __3D__
  double position[8]  = {0., 0.5, 1., 0., 0.25, 0.5, 0.75, 1.};
  double space_avg[2] = {0.5, 1.5};
  int n_lines = 2;
  int n_points[2] = {3, 5};
#endif
#endif // _MPI

  if( n_lines != lines.GetLength() )
  {
    ErrThrow("wrong number of lines computed, ", n_lines, " != ",
             lines.GetLength());
    return 1;
  }

  int idx = 0;

  for( int i=0 ; i<n_lines ; i++ )
  {
    auto line_i = lines.GetLine(i);

    if( n_points[i] != line_i.GetNbPoints() )
    {
      ErrThrow("wrong number of points computed, ", n_points[i], " != ",
              line_i.GetNbPoints());
      return 1;
    }
  
    for( int j=0 ; j<n_points[i] ; j++ )
    {
      if( std::abs(position[idx] - line_i.GetPosition(j) ) > 1.e-10 )
      {
        ErrThrow("wrong position computed, ", position[idx], " != ",
                 line_i.GetPosition(j));
        return 1;
      }

      idx++;
    }

    if( std::abs(space_avg[i]-line_i.space_average_value(fefunction)) > 1.e-10 )
    {
      ErrThrow("wrong mean value computed, ", space_avg[i], " != ",
               line_i.space_average_value(fefunction));
      return 1;
    }
  }

  // delete temporary file used for this test
  std::string dir_name = lines_eval_db["directory_name"];
  if(root && std::remove(dir_name.c_str()) != 0)
  {
    ErrThrow("linesEval_test: Error deleting temporary directory ", dir_name);
  }

  parmoon::parmoon_finalize();
}
