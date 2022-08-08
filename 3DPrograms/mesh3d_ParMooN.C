// =======================================================================
//
// Purpose:  main program for testing the Mesh class (3D) and the 
//           initialization of Domain class using different format
//           (mainly: .mesh, .smesh through TetGen)
//
// Author:   Alfonso Caiazzo
//
// History:  Implemetation started on 12.08.2016
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <Mesh.h>
#include <ParameterDatabase.h>
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

  //Output::print(" ... writing 3d mesh ... ");
  //TCollection *coll = domain.GetCollection(It_Finest, 0);
  //coll->writeMesh("test.mesh");
  

  
  //=========================================================================
  parmoon::parmoon_finalize();
} // end main
