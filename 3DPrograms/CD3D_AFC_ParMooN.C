/**
 * @brief Main program for solving a 3D stationary scalar equation using AFC 
 *        in ParMooN.
 * @author Sashikumaar Ganesan, Clemens Bartsch, Abhinav Jha
 *
 * Implementation started on 2015/01/23. Rework since 2015/10/19.
 *
 */
#include <Domain.h>
#include <Database.h>
#include "ConvectionDiffusion_AFC.h"
#include <Example_CD3D.h>
#include <Chrono.h>
#include "ParMooN.h"

#include <sys/stat.h>

// main program
// =======================================================================
int main(int argc, char* argv[])
{
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
  Output::info("CD3D", "SEQUENTIAL (or OMP...)");
#endif

  Chrono timer;

  TDomain domain(parmoon_db);

  // Intial refinement
  domain.refine_and_get_hierarchy_of_collections(parmoon_db);
  
  //print information on the mesh partition on the finest grid
  domain.print_info("cd3d domain");

  // Choose and construct example.
  Example_CD3D example(parmoon_db);

  timer.restart_and_print("setup(domain, example, database)");
  // Construct the cd3d problem object.
  ConvectionDiffusion_AFC<3> cd3d(domain, parmoon_db);
  timer.restart_and_print("constructing CD3D object");
  
  //=========================================================================
  //Start the actual computations.
  //=========================================================================
   // assemble and solve poisson equation with right-hand side
  cd3d.assemble_poisson();
  cd3d.solve_poisson();
  cd3d.assemble(0); // assemble matrix and rhs
  timer.restart_and_print("Assembling");
  
  cd3d.solve(0);    // solve the system
  timer.restart_and_print("Solving");
  if( cd3d.get_db()["algebraic_flux_correction"].is("afc") )
  {//nonlinear loop necessary
    size_t Max_It = cd3d.get_db()["afc_nonlinloop_maxit"];
    for(unsigned int k = 1;; k++)
    {
      bool converged;
      converged = cd3d.solve(k);

      if ((converged)||(k>= Max_It))
        break;
    }
  }
  
  cd3d.output();   // produce nice output
  timer.restart_and_print("output");
  
  //=========================================================================

  if(my_rank==0)
    Output::print("<<<<< ParMooN Finished: CD3D Main Program >>>>>");

  timer.print_total_time("CD3D_ParMooN program");
  parmoon::parmoon_finalize();
  }
} // end main


