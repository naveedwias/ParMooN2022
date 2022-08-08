// =======================================================================
//
// Purpose:  main program for converting .(x)GEO files into .mesh
//
// Author:   Alfonso Caiazzo
//
// History:  Implementation started on 10.05.2016
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <Mesh.h>
#include <ParameterDatabase.h>
#include "MooNMD_Io.h"
#include "ParMooN.h"

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  auto parmoon_db = parmoon::parmoon_initialize(argc, argv);
  TDomain domain(parmoon_db);
  
  if(!parmoon_db.contains("mesh_file"))
  {
    parmoon_db.add("mesh_file", "converted.mesh", "");
  }
  
  // write .mesh file
  TCollection *coll = domain.GetCollection(It_Finest, 0);
  coll->writeMesh(parmoon_db["mesh_file"]);
  
  // create mesh object and display info for checking
  Mesh m(parmoon_db["mesh_file"]);
  m.info();
  // set boundary description
  m.setBoundary(parmoon_db["boundary_file"]);
  // display boundary info
  m.boundary.info();
  
  //=========================================================================
  parmoon::parmoon_finalize();
} // end main
