// =======================================================================
//
// Purpose:  main program for solving a stationary scalar equation using ParMooN
//
// Author:   Sashikumaar Ganesan
//
// History:  Implementation started on 22.08.2014
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include "ConvectionDiffusion_AFC.h"
#include "ParMooN.h"

#include <sys/stat.h>
#include <sys/types.h>

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  auto parmoon_db = parmoon::parmoon_initialize(argc, argv);
  
  TDomain domain(parmoon_db);
  // refine grid
  domain.refine_and_get_hierarchy_of_collections(parmoon_db);
  
  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    domain.PS("Domain.ps", It_Finest, 0);
   
  //=========================================================================  
  //Create an object of CD2D AFC type
  ConvectionDiffusion_AFC<2> cd2d(domain, parmoon_db);
  //If we have MCL for Balance laws then we need to assemble for Poisson problem with
  //homogenous Neumann condition. This "fairly" unused parameter (only used in AFC Adaptive)
  //tells whether to use Neumann or Dirichlet BCs
  if(cd2d.get_db()["afc_limiter"].is("monolithic_steady"))
    TDatabase::ParamDB->INTERNAL_COERCIVITY = 27;

  //Assemble an object with Neumann BCs
  ConvectionDiffusion_AFC<2> cd2d_poisson(domain, parmoon_db);
  // assemble and solve poisson equation with right-hand side
  cd2d_poisson.assemble_poisson();
  cd2d_poisson.solve_poisson();
  //Reset the unused parameter to -1
  TDatabase::ParamDB->INTERNAL_COERCIVITY = -1;

  //Assemble the Diffusion matrix (this differs from previous assembly as this is full
  //Neumann whereas in the last case for uniqueness we have one point as Dirichlet)
  cd2d.diffusion_matrix = cd2d_poisson.diffusion_matrix;
  //Copy the Poisson solution
  cd2d.poisson_sol = cd2d_poisson.poisson_sol;
  //Here the computation for AFC Scheme start
  //Assemble the CDR equation
  cd2d.assemble(0);
  cd2d.solve(0);
  
  if( cd2d.get_db()["algebraic_flux_correction"].is("afc") )
  {//nonlinear loop necessary
    size_t Max_It = cd2d.get_db()["afc_nonlinloop_maxit"];
    for(unsigned int k = 1;; k++)
    {
      bool converged;
      converged = cd2d.solve(k);

      if ((converged)||(k>= Max_It))
        break;
    }
  }


  cd2d.output();
  //=========================================================================
  parmoon::parmoon_finalize();
} // end main
