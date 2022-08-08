// =======================================================================
//
// Purpose:  main program for testing the Mesh class and the 
//           mesh-conversion routines
//
// Author:   Alfonso Caiazzo
//
// History:  Implementation started on 22.03.2016
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

  unsigned int testType = 4;
  if (testType==0) {
    Output::print(" Test: ");
    Output::print("   * read .PRM and .GEO file and write the corresponding .mesh");

    // write .mesh file
    TCollection *coll = domain.GetCollection(It_Finest, 0);
    coll->writeMesh("test.mesh");
    
  } else if (testType==1) {
    //=========================================================================
    Output::print(" Test: ");
    Output::print("   * read .PRM and .GEO file and write the corresponding .mesh");
    Output::print("   * create a Mesh and display info");
    Output::print("   * rewirte a new .GEO file");

    // initialize domain
    domain.Init(parmoon_db["boundary_file"], parmoon_db["geo_file"]);

    // write .mesh file
    TCollection *coll = domain.GetCollection(It_Finest, 0);
    coll->writeMesh("output/meshFromPRMandGEO.mesh");

    // create mesh object and display info
    Mesh m("output/meshFromPRMandGEO.mesh");
    m.info();
    
    // set boundary description
    m.setBoundary(parmoon_db["boundary_file"]);
    // display boundary info
    m.boundary.info();
    // write new .GEO
    m.writeToGEO("output/GEOFromMesh.xGEO");

    Output::print(" ... done. Compare the new and the old .GEO files. ");
    
    } else if (testType==2) {

    Output::print(" Test: ");
    Output::print("   * read .PRM and the new .GEO file and write the corresponding .mesh");
    domain.Init(parmoon_db["boundary_file"], "output/GEOFromMesh.xGEO");
    TCollection *coll = domain.GetCollection(It_Finest, 0);
    coll->writeMesh("output/meshFromMesh.mesh");

    Output::print(" ... done. Compare the mesh file with the one generated with testType = 1. ");

     } else if (testType==3) {

    Output::print(" Test: ");
    Output::print("   * read .PRM and the .mesh to initialize the domain");
    Output::print("   * write out the mesh file corresponding to the grid");
    domain.InitFromMesh(parmoon_db["boundary_file"],
		"output/meshFromMesh.mesh"); 
    TCollection *coll = domain.GetCollection(It_Finest, 0);
    coll->writeMesh("output/meshFromPRMandMesh.mesh");

    Output::print(" ... done. Compare now meshFromPRMandGEO.mesh and meshFromPRMandMesh.mesh");

  } else if (testType==4) {
    
    Output::print(" Test: ");
    Output::print("   * read .PRM and .mes file and write the corresponding .GEO");

    // create mesh object and display info
    if(!parmoon_db.contains("mesh_file"))
    {
      ErrThrow("please provide a mesh file via a parameter 'mesh_file'.");
    }
    Mesh m(parmoon_db["mesh_file"]);
    m.info();
    m.setBoundary(parmoon_db["boundary_file"]);
    // write new .GEO
    m.writeToGEO(parmoon_db["geo_file"]);
  } 
  
  
  //=========================================================================
  parmoon::parmoon_finalize();
} // end main
