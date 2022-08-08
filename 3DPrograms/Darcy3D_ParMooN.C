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
  
  //=========================================================================
  Darcy<3> darcy3d(domain, parmoon_db);
  darcy3d.assemble();
  darcy3d.solve();
  darcy3d.output();
  //=========================================================================
  
  parmoon::parmoon_finalize();
}
