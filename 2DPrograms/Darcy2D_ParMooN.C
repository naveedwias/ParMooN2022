#include <Domain.h>
#include <Database.h>
#include <Darcy.h>
#include "MainUtilities.h"
#include "ParMooN.h"

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
  Darcy<2> darcy2d(domain, parmoon_db);
  darcy2d.assemble();
  darcy2d.solve();
  darcy2d.output();
  //=========================================================================
  
  parmoon::parmoon_finalize();
}
