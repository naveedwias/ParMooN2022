#include <Domain.h>
#include <Database.h>
#include "ConvectionDiffusion.h"
#include <Example_CD3D.h>
#include <Chrono.h>
#include "ParMooN.h"

// main program
// =======================================================================
int main(int argc, char* argv[])
{
  auto parmoon_db = parmoon::parmoon_initialize(argc, argv);
#ifdef _MPI
  TDatabase::ParamDB->Comm = MPI_COMM_WORLD;
  int my_rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if(my_rank==0)
  {
    Output::print("<<<<< Running ParMooN: CD3D Main Program >>>>>");
    Output::info("CD3D", "MPI, using ", size, " processes");
  }
#else
  int my_rank = 0;
  Output::print("<<<<< Running ParMooN: CD3D Main Program >>>>>");
#endif

  Chrono timer;

  TDomain domain(parmoon_db);

  // Intial refinement
  domain.refine_and_get_hierarchy_of_collections(parmoon_db);
  
  //print information on the mesh partition on the finest grid
  domain.print_info("cd3d domain");

  timer.restart_and_print("setup(domain, database)");
  // Construct the cd3d problem object.
  ConvectionDiffusion<3> cd3d(domain, parmoon_db);
  timer.restart_and_print("constructing CD3D object");
  
  //=========================================================================
  //Start the actual computations.
  //=========================================================================

  cd3d.assemble(); // assemble matrix and rhs
  timer.restart_and_print("Assembling");
  
  cd3d.solve();    // solve the system
  timer.restart_and_print("Solving");
  
  cd3d.output();   // produce nice output
  timer.restart_and_print("output");
  
  //=========================================================================

  if(my_rank==0)
    Output::print("<<<<< ParMooN Finished: CD3D Main Program >>>>>");

  timer.print_total_time("CD3D_ParMooN program");
  parmoon::parmoon_finalize();
} // end main


