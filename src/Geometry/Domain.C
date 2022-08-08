#ifdef _MPI
#include "mpi.h"
#include <MeshPartition.h>
#endif

#include <Database.h>
#include <Domain.h>
#include <JointEqN.h>
#include <MacroCell.h>
#include <MooNMD_Io.h>
#include <BdSphere.h>

#include <Quadrangle.h>

#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <stdio.h>
#include <limits>

#include <BdCircle.h>
#include <BdLine.h>
#include <BdSpline.h>
#include <BdPolygon.h>
#include <BdNonUniformSpline.h>

#ifdef __2D__
  #include <IsoBoundEdge.h>
  #include <IsoInterfaceJoint.h>
#endif

#ifdef __3D__
  #include <BoundFace.h>
  #include <BdNoPRM.h>
  #include <BdPlane.h>
  #include <BdWall.h>
  #include <BdSphere.h>
  #include <BdCylinder.h>
  #include <IsoInterfaceJoint3D.h>
  #include <InterfaceJoint3D.h>
  #include <IsoBoundFace.h>

  #include <BDEdge3D.h>
  #include <InnerEdge.h>
  #include <IsoEdge3D.h>
  #include <ChannelTauRoutines.h>
#endif

#include <PeriodicDomain.h>
#include "RefinementStrategy.h"

#include <string.h>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <cmath>
#include "Multigrid.h"

ParameterDatabase TDomain::default_domain_parameters()
{
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();

  db.add("boundary_file", "Default_UnitSquare",
         "This is a file describing the boundary of the computational domain. "
         "You probably want to adjust this to be the path to some file which "
         "typically has the extension 'PRM'. See the documentation for GEO and "
         "PRM files.",
         {"Default_UnitSquare", "Default_UnitCube",
          "Default_Block_1x2x3", ""});

  db.add("geo_file", "UnitSquare",
         "This files describes the computational mesh. You probably want to "
         "adjust this to be the path to some file which typically has the "
         "extension 'GEO', 'xGEO', or '.mesh'. See the documentation for GEO "
         "and PRM files.",
         {"UnitSquare", "TwoTriangles", "TwoTriangles2",
          "Default_UnitCube_Hexa", "Default_UnitCube_Tetra",
          "Default_Block_1x2x3_Tetra", "Default_Block_1x2x3_Hexa"});
#ifdef __3D__
  // change the default values in 3D
  db["boundary_file"] = "Default_UnitCube";
  db["geo_file"] = "Default_UnitCube_Hexa";
#endif

  db.add("refinement_n_initial_steps", 0u,
         "This is the number of refinement steps before any computation "
         "starts. Usually the mesh is uniformly refined. In a multigrid "
         "program, this determines the number of uniform refinements until the "
         "finest mesh.", 0u, 20u);

  db.add("refinement_n_steps_before_partitioning", 100u,
         "The number of refinement steps before the domain is partitioned. "
         "This parameter is only useful in mpi mode. It has to be smaller or "
         "equal to 'refinement_n_initial_steps', larger values mean equal to "
         "'refinement_n_initial_steps'. The mesh needs to be fine enough to be "
         "partitioned. Therefore the default is to partition on the finest "
         "mesh. Note that if you use multigrid this value is not considered "
         "and the partitioning is done on the coarsest multigrid level",
         0u, 100u);

  db.add("refinement_max_n_adaptive_steps", 0u,
         "A maximum number of adaptive refinement steps "
         "which may be applied to this domain.",
         0u, 1000u);

  db.add("refinement_max_n_uniform_steps", 0u,
         "A maximum number of uniform refinement steps "
         "which may be applied to this domain.",
         0u, 1000u);

  db.add("refinement_final_step_barycentric", false,
         "If this is set to true, the domain will be refined as usual except "
         "the last refinement step is not uniform but barycentric. Then the "
         "Scott-Vogelius finite element pair is inf-sup stable. If "
         "'refinement_n_initial_steps' is zero, this parameter has no effect. "
         " NOTE: This is not yet correctly implemented, you have to do this by "
         "hand.");

  db.add("conforming_closure", true, "Chose whether conforming grid or not");

  db.add("fix_boundary_elements", false,
    "Try to fix cells with more than one boundary face");

  db.add("fix_boundary_elements_write_mesh", false,
    "Write the mesh to a new file after fixing boundary elements.");

  db.add("fix_boundary_elements_write_mesh_file",
    std::string("fixed_mesh.mesh"),
    "Write the mesh to a new file after fixing boundary elements. "
    "(Path to write to)");

  db.add("write_adaptive_mesh", false, "If this is set true then the mesh is written"
          "after refinment in .ps format. Useful in adaptive refinments");

  db.add("read_metis", false , "This Boolean will state if you read a file "
         "which contains the partition of the cells on the processors.");

  db.add("read_metis_file", std::string("mesh_partitioning_file.txt"),
        "The Mesh-file will be read here.");

  db.add("write_metis", false , "This Boolean will state if you write out "
         "which cell belongs to which processor into a file (see parameter "
         "'write_metis_file'.");

  db.add("write_metis_file", std::string("mesh_partitioning_file.txt"), 
         "The partitioning of the mesh will be written here.");

  db.add("metis_optimize_distributed", false,
         "If true, will try to optimize network throughput when running on "
         "multiple physical devices by reassigning ranks.");

  db.add("sandwich_grid", false, "If 'true', then this domain should"
         "be constructed as a sandwich grid from a 2d geometry and mesh.");

  db.add("periodic_translations", std::vector<double>(1, 0.),
         "List of vectors defining the translation between periodic "
         "boundaries.");

  db.merge(Multigrid::default_multigrid_database());

  return db;
}

ParameterDatabase TDomain::default_sandwich_grid_parameters()
{
  ParameterDatabase db("Sandwich Grid Database");

  db.add("drift_x", 0.0, "Grid stretch in x direction.", 0.0, 1.e10);

  db.add("drift_y", 0.0, "Grid stretch in y direction.", 0.0, 1.e10);

  db.add("drift_z", 1.0, "Grid stretch in z direction.", 0.0, 1.e10);

  db.add("conic_scale", 0.0, "Conic scale in x and y direction.", 0.0, 1.e10);

  db.add("n_layers", 1, "Number of cell layers to be stacked "
      "in the sandwich geometry. If this is specified and 'lambda' "
      "is at its default value {0,1}, a total number of 'n_layer'"
      "cell layers are equally distributed across the height of the"
      "sandwich grid. Note that 'lambda' is prioritised over 'n_layers'.",
      1, 1000000000 );

  db.add("lambda", {0.0,1.0}, "The vector of seperators of the cell "
      "layers. Minimum length is 2, the first value must be 0 and "
      "the last value must be 1. The values in between give the "
      "relative positions of the layer seperators. If this is "
      "specified other than {0,1} it is given priority over n_layers.");

  return db;
}

// A little helper function, copied from
// https://stackoverflow.com/questions/874134/find-if-string-ends-with-another-string-in-c,
// which makes it easier to decide, whether a string has a specific ending.
bool ends_with (std::string const &fullString, std::string const &ending)
{
  if (fullString.length() >= ending.length())
  {
    return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
  }
  else
  {
    return false;
  }
}

TDomain::TDomain(const ParameterDatabase& param_db)
  : Interfaces(nullptr), RefLevel(0), db(default_domain_parameters())
{
  // get the relevant parameters from the new database
  db.merge(param_db, false);

#ifdef __3D__
  // Check if this should be a sandwich grid
  if (db["sandwich_grid"])
  {
    // Try to find a nested database with the correct name

    ParameterDatabase sandwich_db = default_sandwich_grid_parameters();
    ParameterDatabase sandwich_in_db("");

    try
    {
      sandwich_in_db = param_db.get_nested_database("Sandwich Grid Database");
    }
    catch(...)
    {
      Output::warn("DOMAIN","If the domain should be built as sandwich grid, you"
          " must give a nested database named '[Sandwich Grid Database]'.");
      Output::warn("CONT'D","Did not find such a database, proceeding with defaults.");

    }

    // Merge the sandwich input database into the default one. Here we allow
    // for the creation of new parameters, for easy include of example specific
    // parameters, which should not appear in the default database.
    sandwich_db.merge(sandwich_in_db, true);

    //The sandwich database is then nested into the domain database.
    db.add_nested_database(sandwich_db);
  }
#endif

  Output::info<4>("Domain" "Starting domain initialization.");

  std::string geoname = db["geo_file"];
  std::string boundname = db["boundary_file"];

  // Find out what kind of geometry input files are given.
  if( ends_with(geoname, "GEO") && ends_with(boundname, ".PRM") )
  {//GEO and PRM file given (old-school MooNMD)
    Output::root_info("Domain", "Initializing Domain using GEO file ", geoname,
          " and .PRM file ", boundname);
    Init(boundname, geoname);
  }
  else if( ends_with(geoname, ".mesh") && ends_with(boundname, ".PRM") )
  {//.mesh file (medit format) and PRM file given
    Output::root_info("Domain", "Initializing Domain using .mesh file ", geoname,
             " and .PRM file ", boundname);
    InitFromMesh(boundname, geoname);
  }
#ifdef __3D__
  else if( ends_with(geoname, ".mesh") )
  {// in 3D, we allow domain initialization with .mesh file and without PRM
    Output::root_info("Domain", "Initializing Domain using .mesh file ", geoname,
             " and no .PRM file.");
    InitFromMesh(boundname, geoname);
  }
#endif
  else 
  {// If nothing else fits, the program supposes a default geometry.
    Output::root_info("Domain", "Trying to find a default geometry and mesh "
        " for ", geoname, " and ", boundname);
    // default cases for the tests
    Init(boundname, geoname);
  }
}

TDomain::~TDomain()
{
  delete [] StartBdCompID;
  if (Interfaces)
    delete [] Interfaces;
  for(int i = 0; i < N_BoundParts; ++i)
    delete BdParts[i];
  delete [] BdParts;
  for(auto coll: gridCollections)
    delete coll;
  gridCollections.clear();
  /// @todo delete cells, joints, vertices here
  delete [] CellTree;
// #ifdef __3D__
//   for(auto p : meshVertices)
//     delete p;
//   for(auto p : meshJoints)
//     delete p;
// #endif
}

// Methods
int TDomain::GetBdPartID(int BdCompID)
{
  int i;

  for (i=N_BoundParts-1;i>0;i--)
    if (BdCompID >= StartBdCompID[i]) break;

  return i;
}

int TDomain::GetLocalBdCompID(int BdCompID)
{
  int i;

  for (i=N_BoundParts-1;i>0;i--)
    if (BdCompID >= StartBdCompID[i]) break;

  return BdCompID - StartBdCompID[i];
}

#ifdef __2D__

void TDomain::Init(const std::string& PRM, const std::string& GEO)
{

  // Interpret the PRM string
  if (PRM == "Default_UnitSquare")
  {// default case - boundary of the [0,1] unit square
    initializeDefaultUnitSquareBdry();
  }
  else if (PRM == "DrivenCavitySquare")
  {// default case - boundary of the [-1,1] unit square
    initializeDrivenCavitySquareBdry();
  }
  else
  {
    // non-default: read in from .PRM file
    std::ifstream bdryStream(PRM);
    if (!bdryStream)
    {
      ErrThrow("Cannot open .PRM file ", PRM);
    }

    //do the actual read in
    ReadBdParam(bdryStream);
  }

  // Interpret the GEO string - start with default meshes
  if (GEO == "TwoTriangles") // diagonal from bottom left to top right
  {
    TwoTriangles();
  }
  else if (GEO == "TwoTriangles2") // diagonal from top left to bottom right
  {
    TwoTriangles(false);
  }
  else if (GEO == "TwoTrianglesRef")
  {
    TwoTrianglesRef();
  }
  else if (GEO ==  "UnitSquare")
  {
    UnitSquare();
  }
  else if (GEO == "DrivenCavitySquareQuads")
  {
    DrivenCavitySquareQuads();
  }
  else if (GEO ==  "UnitSquareRef")
  {
    UnitSquareRef();
  }
  else if (GEO == "SquareInSquare")
  {
    SquareInSquare();
  }
  else if (GEO == "UnitSquare_US22")
  {
    UnitSquare_US22();
  }
  else if (GEO == "SquareInSquareRef")
  {
    SquareInSquareRef();
  }
  else if (GEO == "PeriodicSquares")
  {
    PeriodicSquares();
  }
  else if (GEO == "PeriodicTriangles")
  {
    PeriodicTriangles();
  }
  else if (GEO ==  "PeriodicSquaresLarge")
  {
    PeriodicSquaresLarge();
  }
  else if (GEO == "PeriodicTrianglesLarge")
  {
    PeriodicTrianglesLarge();
  }
  else if (GEO == "PeriodicRectangle_2_4")
  {
    PeriodicRectangle_2_4();
  }
  else
  {
    // check if the GEO file actually is an .xGEO-file
    bool isxGEO = ends_with(GEO, ".xGEO");

    // make an input file string from the file "GEO"
    std::ifstream geo_stream(GEO);
    if (!geo_stream)
    {
      ErrThrow("Cannot open GEO file ", GEO);
    }

    // and call the read-in method
    ReadGeo(geo_stream,isxGEO);
  }
}


#else // 3D
void TDomain::Init(const std::string& PRM, const std::string& GEO)
{

  // start with read in of boundary description
  bool prm_is_sandwich = false;
  if (PRM == "Default_UnitCube")
  {
    // one implemented default case: unit cube
    initializeDefaultCubeBdry();
  }
  else if (PRM == "Default_Block_1x2x3")
  {
    initialize_block_1x2x3_bdry();
  }
  else
  {
    // "PRM" interpreted as file path to PRM-file.
    std::ifstream bdryStream(PRM);
    if (!bdryStream)
      ErrThrow("Cannot open PRM file ", PRM);

    ReadBdParam(bdryStream, prm_is_sandwich);
  }

  if(db["sandwich_grid"])
  {
    //this ought to be a sandwich grid
    if(!prm_is_sandwich)
      ErrThrow("The specified .PRM file ", PRM ," is not suitable for a "
          "sandwich geometry. It seems to lack a TBdWall boundary.");

    // make an input file string from the file "GEO"
    ReadSandwichGeo(GEO);

    return;
  }

  if(prm_is_sandwich)
    Output::root_warn("DOMAIN", "Your .PRM file is adapted to sandwich grid. This might "
        "cause problems, since 'sandwich_grid' is false. (Untested!)");

  // Check for default geometries.
  if (GEO == "TestGrid3D")
  {
    TestGrid3D();
  }
  else if (GEO == "Default_UnitCube_Hexa")
  {
    initialize_cube_hexa_mesh();
  }
  else if (GEO == "Default_UnitCube_Tetra")
  {
    initialize_cube_tetra_mesh();
  }
  else if (GEO == "Default_Block_1x2x3_Tetra")
  {
    initialize_block_1x2x3_tetra_mesh();
  }
  else if (GEO == "Default_Block_1x2x3_Hexa")
  {
    initialize_block_1x2x3_hexa_mesh();
  }
  else
  { // non-default case, try to read an initial geometry from file
    // make an input file string from the file "GEO"
    std::ifstream geo_stream(GEO);
    if (!geo_stream)
    {
      ErrThrow("Cannot open GEO file ", GEO);
    }
    bool isxGEO = ends_with(GEO, ".xGEO"); // check if the GEO file actually is an .xGEO-file
    ReadGeo( geo_stream, isxGEO );
  }
}

void TDomain::ReadSandwichGeo(const std::string& file_name,
                              const std::string& prm_file_name)
{

  // The constructor took care of this nested database being present.
  const ParameterDatabase& sw_db = db.get_nested_database("Sandwich Grid Database");
  double drift_x = sw_db["drift_x"];
  double drift_y = sw_db["drift_y"];
  double drift_z = sw_db["drift_z"];
  double conic_scale = sw_db["conic_scale"];
  int n_layers = sw_db["n_layers"]; //layers of cells (!)

  //Three C-style arrays needed for the call to MakeSandwichGrid.
  double *DCORVG = nullptr;
  int* KVERT = nullptr;
  int* KNPR = nullptr;
  int N_Vertices;
  int NVE;

  // First case: the file_name is a .GEO file.
  if(ends_with(file_name, ".GEO"))
  {
    std::ifstream dat(file_name);
    if (!dat)
      ErrThrow("Cannot open .GEO file ", file_name);
    char line[100];
    int NVpF, NBCT;

    dat.getline (line, 99);
    dat.getline (line, 99);

    // determine dimensions for creating arrays
    dat >> N_RootCells >> N_Vertices >> NVpF >> NVE >> NBCT;
    dat.getline (line, 99);
    dat.getline (line, 99);

    // allocate auxillary fields
    DCORVG =  new double[2*N_Vertices];
    KVERT = new int[NVE*N_RootCells];
    KNPR = new int[N_Vertices];
    // read fields
    for (int i=0;i<N_Vertices;i++)
    {
      dat >> DCORVG[2*i] >> DCORVG[2*i + 1];
      dat.getline (line, 99);
    }

    dat.getline (line, 99);

    for (int i=0;i<N_RootCells;i++)
    {
      for (int j=0;j<NVE;j++)
        dat >> KVERT[NVE*i + j];
      dat.getline (line, 99);
    }

    dat.getline (line, 99);

    for (int i=0;i<N_Vertices;i++)
      dat >> KNPR[i];
  }
  // Second case: the file_name is a .mesh file.
  else if (ends_with(file_name, ".mesh"))
  {
    // read the mesh file and prepare input for MakeSandwichGrid
     Mesh m(file_name);
     m.setBoundary(prm_file_name);

     size_t numberOfElements = m.triangle.size() + m.quad.size();
     size_t maxNVertexPerElem = 3;
     // make the ParMooN-grid
     if (m.quad.size()) {
       maxNVertexPerElem = 4;
     }
     if (numberOfElements==0) {
       ErrThrow(" ** Error(Domain::Init) the mesh has no elements");
     }

     // vertices data
     size_t n_vertices = m.vertex.size();
     KNPR = new int[n_vertices];
     DCORVG =  new double[2*n_vertices];

     // fill the DCORVG array with GEO-like coordinates of two-dimensional points
     for (size_t i=0; i<n_vertices; i++) {
       // check if the vertex is on the boundary
       double localParam;
       int partID = m.boundary.isOnComponent(m.vertex[i].x,m.vertex[i].y,localParam);
       if (partID>=0) {
         DCORVG[2*i] = localParam;
         DCORVG[2*i+1] = 0.;
         KNPR[i] = partID+1;
       } else {
         DCORVG[2*i] = m.vertex[i].x;
         DCORVG[2*i+1] = m.vertex[i].y;
         KNPR[i] = 0;
       }
     }
     KVERT = new int[maxNVertexPerElem * numberOfElements];

     // store triangles (+ 0 when using a mixed mesh)
     for (size_t i=0;i<m.triangle.size();i++)
     {
       for (size_t  j=0; j<3; j++)
         KVERT[maxNVertexPerElem*i + j] = m.triangle[i].nodes[j];
       if (maxNVertexPerElem==4)
         KVERT[maxNVertexPerElem*i + 3] = 0;
     }

     // store quadrilaterals
     for (size_t i=0;i<m.quad.size();i++)
     {
       for (size_t  j=0; j<4; j++)
         KVERT[maxNVertexPerElem* (m.triangle.size() + i) + j] = m.quad[i].nodes[j];
     }
     // N_RootCells is an internal Domain variable used in other functions
     N_RootCells = numberOfElements;
     N_Vertices = m.vertex.size();
     NVE = maxNVertexPerElem;
  }

  // figure out the relative placing of sandwich grid layers.
  std::vector<double> lambda = sw_db["lambda"];

  bool control_by_lambda = true;
  if(!std::is_sorted(lambda.begin(), lambda.end()))
  {//check if lambda is valid
    Output::root_warn("DOMAIN", "Invalid lambda for sandwich grid given. "
                 "The entries must be in ascending order. "
                 "Opting for control parameter n_layers.");
    control_by_lambda = false;
  }
  if(*std::min_element(lambda.begin(), lambda.end()) != 0 ||
     *std::max_element(lambda.begin(), lambda.end()) != 1)
  {//check if lambda is valid
    Output::root_warn("DOMAIN", "Invalid lambda for sandwich grid given. "
                 "The vector must start with 0 and end with 1. "
                 "Opting for control parameter n_layers.");
    control_by_lambda = false;
  }
  if(lambda == std::vector<double>({0,1}))
  {//default lambda is given - take n_layers instead
    Output::root_info("DOMAIN", "Default lambda for sandwich grid given. "
                 "Opting for control parameter n_layers.");
    control_by_lambda = false;
  }
  if(!control_by_lambda)
  {//n_layers is used as control parameter

    Output::root_info("DOMAIN", "Taking sandwich grid database "
        "parameter 'n_layers' = ", sw_db["n_layers"]);
    int n_node_layers = sw_db["n_layers"];
    ++n_node_layers; //count up to get the number of node (and not cell) layers

    //fill a lambda with equally spaced numbers
    lambda.resize(n_node_layers);
    for(int i=0 ; i<n_node_layers ; i++)
    {
      lambda[i] = i * (1.0/n_layers);
    }
  }
  else
  {//control by lambda
    Output::root_info("DOMAIN", "Taking sandwich grid database "
      "parameter 'lambda'.");
  }

  MakeSandwichGrid(DCORVG, KVERT, KNPR, N_Vertices, NVE,
                   drift_x, drift_y, drift_z, conic_scale, lambda);

  delete[] DCORVG;
  delete[] KVERT;
  delete[] KNPR;

}

#endif

#ifdef __2D__
void TDomain::ReadBdParam(std::istream& dat)
#else
void TDomain::ReadBdParam(std::istream& dat, bool& sandwich_flag)
#endif
{
#ifdef _MPI
  int rank; // out_rank=int(TDatabase::ParamDB->Par_P0);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  char line[100];
  int i, j, CurrentBdPart, N_BdComp, CompType, CompID = 0, N_Spls;
#ifdef __2D__
  TBoundComp2D *BdComp;
#else
  TBoundComp3D *BdComp=nullptr;
  TBoundComp2D *BdComp2D=nullptr;
  sandwich_flag = false;
#endif



  // determine dimensions for creating arrays

  // get the number of boundaries (inner and outer) of the domain
  dat.getline (line, 99);
  dat >> N_BoundParts;
  dat.getline (line, 99);

  BdParts = new TBoundPart*[N_BoundParts];
  Interfaces = new int[N_BoundParts];
  N_BoundComps = 0;
  // StartBdCompID: index in the bdpart list where the new BdPart starts
  StartBdCompID = new int[N_BoundParts + 1];
  StartBdCompID[0] = 0; // set the first to 0

  for (i=0;i<N_BoundParts;i++)
  {
    dat.getline (line, 99); // IBCT
    dat >> CurrentBdPart;
    dat.getline (line, 99);
    Interfaces[i] = CurrentBdPart;
    CurrentBdPart = std::abs(CurrentBdPart); // it can be negative (for orientation)
    if (i+1 != CurrentBdPart)
    {
#ifdef _MPI
      if(rank==0)
#endif
        ErrThrow("different number of boundary part\nCurrentBdPart ", i, "  ",
                 CurrentBdPart);
    }

    // get number of components f the bdpart i
    dat.getline (line, 99);
    dat >> N_BdComp;
    dat.getline (line, 99);

    BdParts[i] = new TBoundPart(N_BdComp);
    N_BoundComps += N_BdComp;
    SetStartBdCompID(N_BoundComps, i+1); // set the start ID of the next bdcomp

    dat.getline (line, 99);
    for (j=0;j<N_BdComp;j++)
    {
      dat >> CompType >> N_Spls;  //ITYP NSPLINE NPAR
      dat.getline (line, 99);

#ifdef __2D__
      // 2D types: Line (1), Circle (2), Spline (3), Poygon (4), NonUnif Spline (5)
      switch (abs(CompType))
      {
        case 1: BdComp = new TBdLine(CompID++);
                break;
        case 2: BdComp = new TBdCircle(CompID++);
                break;
        case 3: BdComp = new TBdSpline(CompID++, N_Spls);
                break;
        case 4: BdComp = new TBdPolygon(CompID++, N_Spls);
                break;
        case 5: BdComp = new TBdNonUniformSpline(CompID++, N_Spls);
                break;
        default:
#ifdef _MPI
                if(rank==0)
#endif
                  ErrThrow("ReadParam.C: Boundary type not implemented");
      }
#else
      // 3D types: Line (1), Circle (2), Spline (3), Poygon (4), NonUnif Spline (5),
      //           Plane (10), Sphere (11)
      switch (abs(CompType))
      {
        case 1: BdComp2D = new TBdLine(CompID++);
                break;
        case 2: BdComp2D = new TBdCircle(CompID++);
                break;
        case 3: BdComp2D = new TBdSpline(CompID++, N_Spls);
                break;
        case 4: BdComp2D = new TBdPolygon(CompID++, N_Spls);
                break;
        case 5: BdComp2D = new TBdNonUniformSpline(CompID++, N_Spls);
                break;
        case 10: BdComp = new TBdPlane(CompID++);
                 break;
        case 11: BdComp = new TBdSphere(CompID++);
                 break;
        case 4711: BdComp = new TBdNoPRM(CompID++); // create grid without PRM file
                   break;
        default:
#ifdef _MPI
                   if(rank==0)
#endif
                     ErrThrow("ReadParam.C: Boundary type (3D) not implemented");
      }

      if(abs(CompType)<10)
      {
        BdComp = new TBdWall(CompID-1, BdComp2D);
        sandwich_flag = true;
      }
#endif // 3D
      BdParts[i]->SetBdComp(j, BdComp);

      if(CompType<0)
      {
        BdComp->SetFreeBoundaryStatus(true);
        Output::root_info("Domain", i, " ReadBdParam : ", j);
      }
    }
  }

  dat.getline (line, 99);
  for (i=0;i<N_BoundParts;i++)
  {
    N_BdComp = BdParts[i]->GetN_BdComps();
    for (j=0;j<N_BdComp;j++)
    {
      BdParts[i]->GetBdComp(j)->ReadIn(dat);
    }
  }

  // read HOLES (if any)
  dat.getline (line, 99);
  N_Holes = -12345;
  if (dat.eof())
    N_Holes = 0;
  else
    dat >> N_Holes;

  if(N_Holes == -12345)
    N_Holes = 0;

  dat.getline (line, 99);

  if (N_Holes)
  {
    // coordinates of a point in a hole
    PointInHole = new double[2*N_Holes];

    dat.getline (line, 99);
    for (i=0;i<N_Holes;i++)
    {
      dat >> PointInHole[2*i] >> PointInHole[2*i+1];
      dat.getline (line, 99);
    }
  }
  else
    PointInHole = nullptr;

  dat.getline (line, 99);
  N_Regions = -12345;
  if (dat.eof())
    N_Regions = 0;
  else
    dat >> N_Regions;

  if(N_Regions == -12345)
    N_Regions = 0;

  dat.getline (line, 99);

  // read REGIONS (if any)
  if (N_Regions)
  {
    PointInRegion = new double[4*N_Regions];

    dat.getline (line, 99);
    for (i=0;i<N_Regions;i++)
    {
      dat >> PointInRegion[4*i] >> PointInRegion[4*i+1];
      PointInRegion[4*i+2] = i;
      PointInRegion[4*i+3] = 10000;
      dat.getline (line, 99);
    }
  }
  else
    PointInRegion = nullptr;
}

// initialize a domain from a mesh and a boundary(PRM) file
void TDomain::InitFromMesh(const std::string& PRM, const std::string& MESHFILE)
{
#ifdef __2D__
  Output::root_info<3>("Domain", "TDomain:: InitFromMesh using ", PRM, " and ", MESHFILE);
  //make an input file string from the file "PRM"
  std::ifstream bdryStream(PRM);
  if (!bdryStream)
  {
    ErrThrow(" ** Error(TDomain::Init) cannot open PRM file ", PRM);
  }
  //do the actual read in
  ReadBdParam(bdryStream);

  // read mesh
  Mesh m(MESHFILE);
  m.setBoundary(PRM);
  unsigned int numberOfElements = m.triangle.size() + m.quad.size();
  unsigned int maxNVertexPerElem = 3;
  // make the ParMooN-grid
  if (m.quad.size()) {
    maxNVertexPerElem = 4;
  }
  if (numberOfElements==0) {
    ErrThrow(" ** Error(Domain::Init) the mesh has no elements");
  }

  // vertices data
  double *DCORVG;
  int *KNPR;
  unsigned int n_vertices = m.vertex.size();
  KNPR = new int[n_vertices];
  DCORVG =  new double[2*n_vertices];

  // fill the DCORVG array with GEO-like coordinates of two-dimensional points
  for (unsigned int i=0; i<n_vertices; i++) {

    // check if the vertex is on the boundary
    double localParam;
    int partID = m.boundary.isOnComponent(m.vertex[i].x,m.vertex[i].y,localParam);
    if (partID>=0) {
      DCORVG[2*i] = localParam;
      DCORVG[2*i+1] = 0.;
      KNPR[i] = partID+1;
    } else {
      DCORVG[2*i] = m.vertex[i].x;
      DCORVG[2*i+1] = m.vertex[i].y;
      KNPR[i] = 0;
    }  
  }
  int *KVERT,*ELEMSREF;
  KVERT = new int[maxNVertexPerElem * numberOfElements];
  ELEMSREF = new int[numberOfElements];

  // store triangles (+ 0 when using a mixed mesh)
  for (unsigned int i=0;i<m.triangle.size();i++)
  {
    for (unsigned int  j=0; j<3; j++) 
      KVERT[maxNVertexPerElem*i + j] = m.triangle[i].nodes[j];
    
    if (maxNVertexPerElem==4) KVERT[maxNVertexPerElem*i + 3] = 0;
    
    ELEMSREF[i] = m.triangle[i].reference;
    
  }

  // store quadrilaterals
  for (unsigned int i=0;i<m.quad.size();i++)
  {
    for (unsigned int  j=0; j<4; j++) 
      KVERT[maxNVertexPerElem* (m.triangle.size() + i) + j] = m.quad[i].nodes[j];

    ELEMSREF[m.triangle.size()+i] = m.quad[i].reference;
    
  }
  // N_RootCells is an internal Domain variable used in other functions
  N_RootCells = numberOfElements;
  
  MakeGrid(DCORVG,KVERT,KNPR,ELEMSREF,m.vertex.size(),maxNVertexPerElem);
  delete [] DCORVG;
  delete [] KVERT;
  delete [] KNPR;
  delete [] ELEMSREF;
#else

  bool prm_is_sandwich = false;

  std::ifstream bdry_stream(PRM);
  if (bdry_stream)
  {
    // read .PRM if available (for example for sandwich grid)
    ReadBdParam(bdry_stream, prm_is_sandwich);

    if(db["sandwich_grid"])
    {
      //this ought to be a sandwich grid
      if(!prm_is_sandwich)
        ErrThrow("The specified .PRM file ", PRM ," is not suitable for a "
                 "sandwich geometry. It seems to lack a TBdWall boundary.");

      ReadSandwichGeo(MESHFILE, PRM);

    }
    else
    {
      ErrThrow("The combination of .mesh file and sandwich grid "
          "requires specification of a .PRM file.")
    }
  } else {
    // if not sandwich, one can proceed without a PRM file
    Mesh m(MESHFILE);
    GenerateFromMesh(m);
  }

#endif
}

void TDomain::PS(const char *name, Iterators iterator, int arg)
{
  TCollection * coll = this->GetCollection(iterator, arg);
  this->PS(name, coll);
  delete coll;
}

void TDomain::PS(const char *name, TCollection *Coll)
{
  std::ofstream dat(name);
  if (!dat)
  {
    Output::root_warn("TDomain::PS", "unable to open ", name, " for output");
    return;
  }
  Output::root_info("PS","Generating postscript file ", name);
  
  // scale = 5350 / BoundX;
  // if (7820 / BoundY < scale) scale = 7820 / BoundY;
  double scale = 535 / BoundX;
  if (782 / BoundY < scale) scale = 782 / BoundY;

  int BX = (int) (BoundX * scale + .5);
  int BY = (int) (BoundY * scale + .5);

  dat << "%!PS-Adobe-2.0" << endl;
  dat << "%%Creator: MooN_MD (Volker Behns)" << endl;
  dat << "%%DocumentFonts: Helvetica" << endl;
  // dat << "%%BoundingBox: 300 300 " << 300+BX << " " << 300+BY << endl;
  dat << "%%BoundingBox: 25 25 " << 35+BX << " " << 35+BY << endl;
  dat << "%%Pages: 1" << endl;
  dat << "%%EndComments" << endl;
  dat << "%%EndProlog" << endl;
  dat << "%%Page: 1 1" << endl;
  dat << "/Helvetica findfont 14 scalefont setfont" << endl;
  dat << "/Cshow { dup stringwidth pop 2 div neg 0 rmoveto show } def" << endl;
  dat << "/M { moveto } def" << endl;
  dat << "/L { lineto } def" << endl;
  // dat << "0.10 0.10 scale" << endl;
  // dat << "10.0 setlinewidth" << endl;
  dat << "0.5 setlinewidth" << endl;

  // loop over all cells
  int n_cells = Coll->GetN_Cells();
  for(int i = 0; i < n_cells; i++)
  {
    Coll->GetCell(i)->PS(dat, scale, StartX, StartY); // no cell indices
    //Coll->GetCell(i)->PS(dat, scale, StartX, StartY, i); // with cell indices
  }

  dat << "stroke" << endl;
  dat << "showpage" << endl;
  dat << "%%Trailer" << endl;
  dat << "%%Pages: 1" << endl;
}


int TDomain::Refine()
{
#ifdef __3D__
  if (db["conforming_closure"] && CellTree[0]->GetN_Vertices() == 6 )
  {
    ErrThrow("Adaptive refinement with conforming closure on Hexahedrons does "
             "not work.");
  }
#endif
  TBaseCell *CurrCell;
  int info;

  RefLevel++;

  TDatabase::IteratorDB[It_Finest]->Init(0);

  // loop over all cells
  while( (CurrCell = TDatabase::IteratorDB[It_Finest]->Next(info)))
    CurrCell->Refine(RefLevel);

  return 0;
}

int TDomain::RegRefineAll()
{
  TBaseCell *CurrCell;
  int info;

  RefLevel++;

  TDatabase::IteratorDB[It_Finest]->Init(0);
  // loop over all cells
  while ((CurrCell = TDatabase::IteratorDB[It_Finest]->Next(info)))
  {
    CurrCell->SetRegRefine();
    CurrCell->Refine(RefLevel);
  }

  gridCollections.push_front(this->GetCollection(It_Finest, 0));
  return 0;
}


void TDomain::RefineByIndicator(
#ifdef __2D__ 
  DoubleFunct2D *Indicator
#else
  DoubleFunct3D *Indicator
#endif
              )
{
#ifdef __2D__
  auto db = RefinementStrategy<2>::default_refinement_strategy_database();
  auto f = [&](std::array<double, 3> point)
           { 
             double val;
             Indicator(point[0], point[1], &val);
             return val;
           };
  RefinementStrategy<2> rs(db);
  TCollection* coll = this->GetCollection(It_Finest, 0);
  rs.apply_indicator(coll, f);
  delete coll;
  this->RefineByRefinementStrategy(rs);
#else
    auto db = RefinementStrategy<3>::default_refinement_strategy_database();
    auto f = [&](std::array<double, 3> point)
           { 
             double val;
             Indicator(point[0], point[1], point[2], &val);
             return val;
           };
  RefinementStrategy<3> rs(db);
  TCollection* coll = this->GetCollection(It_Finest, 0);
  rs.apply_indicator(coll, f);
  delete coll;
  this->RefineByRefinementStrategy(rs);
#endif
//   TBaseCell *CurrCell;
//   TVertex *vert;
//   int j, k, info;
//   int Inner, Outer;
//   double x,y,val;
// 
//   TDatabase::IteratorDB[It_Finest]->Init(0);
// 
//   // loop over all cells
//   while ((CurrCell = TDatabase::IteratorDB[It_Finest]->Next(info)))
//   {
//     Inner=0;
//     Outer=0;
//     k=CurrCell->GetN_Vertices();
//     for(j=0;j<k;j++)
//     {
//       vert=CurrCell->GetVertex(j);
//       x=vert->GetX();
//       y=vert->GetY();
//       Indicator(x,y,&val);
//       if(val<=1e-10) Inner++;
//       if(val>=-1e-10) Outer++;
//     }
//     if((Inner>0) && (Outer>0)) 
//     {
//       // there vertices on both sides
//       CurrCell->SetRegRefine();
//     }
//   }
// 
//   Refine();
//   if(conf_closure)
//   { 
//     MakeConfClosure();
//   }
//   else
//   {
//     Gen1RegGrid();
//   }
//   gridCollections.push_front(this->GetCollection(It_Finest, 0));
}

template <int d>
void TDomain::RefineByRefinementStrategy(RefinementStrategy<d>& strategy)
{
  TCollection* collection = this->GetCollection(It_Finest, 0);
  for(int k = 0; k < collection->GetN_Cells(); k++)
  {
    if(strategy.should_refine_cell(k))
    {
      collection->GetCell(k)->SetRegRefine();
    }
    else
    {
      collection->GetCell(k)->SetNoRefinement();
    }
    //Output::print("cell ", k, " : ", collection->GetCell(k)->IsToRefine());
  }
  if (db["conforming_closure"] &&
      (collection->GetCell(0)->GetShapeDesc()->GetType() == Hexahedron ||
       collection->GetCell(0)->GetShapeDesc()->GetType() == Brick))
  {
    ErrThrow("Conforming closure with Hexahedrons is not possible with this"
        " code.");
  }
  delete collection;
  
  if(db["conforming_closure"])
  { 
    MakeConfClosure();
  }
  else
  {
    Refine();
    //Output::print("after TDomain::Refine");
    Gen1RegGrid();
  }
  gridCollections.push_front(this->GetCollection(It_Finest, 0));
}
#ifdef __3D__
template void TDomain::RefineByRefinementStrategy(
  RefinementStrategy<3>& strategy);
#else
template void TDomain::RefineByRefinementStrategy(
  RefinementStrategy<2>& strategy);
#endif



#ifdef __2D__
int TDomain::MakeConfClosure()
{
  int info;
  // delete existing closures
  TDatabase::IteratorDB[It_Finest]->Init(0);
  while(TBaseCell *CurrCell = TDatabase::IteratorDB[It_Finest]->Next(info))
  {
    if(TBaseCell *parent = CurrCell->GetParent())
    {
      Refinements type = parent->GetRefDesc()->GetType();
      // check if the parent has been refined non-regularly, which (almost only)
      // happens during conforming closures
      if(type >= TriBis0 && type <= Quad2Conf3)
      {
        // check whether a child is marked for refinement
        // TriReg means regular refinement for triangles and quadrangles
        type = CurrCell->GetRefDesc()->GetType();
        CurrCell = TDatabase::IteratorDB[It_Finest]->Next(info);
        if (type < CurrCell->GetRefDesc()->GetType())
          type = TriReg;

        if (parent->GetN_Children() == 3)
        {
          CurrCell = TDatabase::IteratorDB[It_Finest]->Next(info);
          if (type < CurrCell->GetRefDesc()->GetType())
            type = TriReg;
        }

        parent->Derefine();

        /// @todo can type be QuadReg and the other children NoRef? If so, the 
        /// following will not correctly refine the parent, will it?
        if(type == TriReg)
          parent->SetRegRefine();
      }
    }
    // there is no parent cell (coarsest level)
  }
    

  // generate a new 1-regular grid
  Refine();
  Gen1RegGrid();

  // initialize clipboards
  TDatabase::IteratorDB[It_Finest]->Init(0);
  while(TBaseCell *CurrCell = TDatabase::IteratorDB[It_Finest]->Next(info))
    CurrCell->SetClipBoard(0);

  // look for elements which have to be refined
  TDatabase::IteratorDB[It_Finest]->Init(0);
  int MaxLevel = TDatabase::IteratorDB[It_Finest]->GetMaxLevel();
  //this->PS("before_closure.ps", It_Finest, 0);
  //Output::print("MaxLevel " << MaxLevel);
  
  for(int i = MaxLevel; i >= 0; i--)
  {
    // get iterator on level i
    TDatabase::IteratorDB[It_EQ]->Init(i);
    // loop over the mesh cells on level i
    while(TBaseCell *CurrCell = TDatabase::IteratorDB[It_EQ]->Next(info))
    {
      // if current cell does not possess children, compute conforming closure
      if(!CurrCell->ExistChildren())
      {
        CurrCell->MakeConfClosure();
      }
    }
    // get iterator on finest level
    TDatabase::IteratorDB[It_Finest]->Init(0);
    // loop over mesh cells on finest level
    while(TBaseCell *CurrCell = TDatabase::IteratorDB[It_Finest]->Next(info))
    {
      // get clip board
      int clip = CurrCell->GetClipBoard();
      //Output::print(i, " clip ", clip);
      // compute refinement rule
      switch (clip)
      {
        case 145:
          CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TriBis0]);
          CurrCell->Refine(RefLevel);
          break;
        case 146: 
          CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TriBis1]);
          CurrCell->Refine(RefLevel);
          break;
        case 148:
          CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TriBis2]);
          CurrCell->Refine(RefLevel);
          break;
        case 273:
          CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + Quad1Conf2]);
          CurrCell->Refine(RefLevel);
          break;
        case 274:
          CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + Quad1Conf3]);
          CurrCell->Refine(RefLevel);
          break;
        case 276:
          CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + Quad1Conf0]);
          CurrCell->Refine(RefLevel);
          break;
        case 280:
          CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + Quad1Conf1]);
          CurrCell->Refine(RefLevel);
          break;
        case 291:
          CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + Quad2Conf1]);
          CurrCell->Refine(RefLevel);
          break;
        case 293:
          CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + QuadBis0]);
          CurrCell->Refine(RefLevel);
          break;
        case 294:
          CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + Quad2Conf2]);
          CurrCell->Refine(RefLevel);
          break;
        case 297:
          CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + Quad2Conf0]);
          CurrCell->Refine(RefLevel);
          break;
        case 298:
          CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + QuadBis1]);
          CurrCell->Refine(RefLevel);
          break;
        case 300:
          CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + Quad2Conf3]);
          CurrCell->Refine(RefLevel);
          break;
        default:
          if (clip > 512)
          {
            CurrCell->SetRegRefine();
            CurrCell->Refine(RefLevel);
          }
          break;
      }
    }
  }
  return 0;
}

#else

int TDomain::CloseGrid(int level)
{
  TBaseCell *CurrCell, *Cell, *Child, *Neigh;
  int i, j, k, clip;
  int N_Edges, N_Children;
  int info, LocEdge, CurrLocEdge, LocFace=0, NeighLocFace;
  int MapType;
  int first_bis, second_bis, neigh_first_bis, neigh_second_bis, neigh_type;
  TJoint *Joint, *LastJoint;
  const int *TmpEF, *TmpEF2;
  int TmpEFMaxLen, TmpEF2MaxLen;
  bool RefineRegular;
  Refinements NeighFaceRef, MyFaceRef;
  int N_ToRefine = 0;

  // Reset ClipBoards
  TDatabase::IteratorDB[It_EQ]->Init(level);
  k=0;
  while( (CurrCell = TDatabase::IteratorDB[It_EQ]->Next(info)))
  {
    N_Children = CurrCell->GetRefDesc()->GetN_Children();

    // Check if CurrCell is refined irregularly but has a son which has to be refined
    if(2 <= N_Children && N_Children < 8) // CurrCell is refined irregularly
    {
      RefineRegular = false;

      for(i=0; i<N_Children; ++i)
      {
        Child = CurrCell->GetChild(i);
        // Check if Child is refined directly
        
        //if(Child->GetRefDesc()->GetType() != NoRef) // TODO
        if(Child->GetClipBoard() > 0)
          RefineRegular = true;

        // Check if Child contains an edge that is refined by a neighbour
        if(!RefineRegular)
        {         
          // TODO
        }
        
        if(RefineRegular)
        {
          CurrCell->SetRegRefine();
          break;
        }
        else
        CurrCell->SetNoRefinement();
      }
    }
    
    // Initialize Clipboard
    if(CurrCell->GetRefDesc()->GetType() >= TetraReg 
       && CurrCell->GetRefDesc()->GetType() <= TetraReg2)
      CurrCell->SetClipBoard(63);
    else
      CurrCell->SetClipBoard(-1);
  }

  /*
   * Set Clipboard to 0 if element is unrefined but contains an edge that is refined by an other tetrahedron
   */
  k = -1;
  TDatabase::IteratorDB[It_EQ]->Init(level);
  while ( (Cell = TDatabase::IteratorDB[It_EQ]->Next(info)))
    {
      Cell->GetShapeDesc()->GetEdgeFace(TmpEF, TmpEFMaxLen);
      N_Edges = Cell->GetN_Edges();

      k++;
      if(Cell->GetClipBoard() == 63)
  {
    for (i=0;i<N_Edges;i++)
      {
        if(Cell->GetClipBoard() & (std::pow(2, i) == 0))
    continue;

        LocEdge = i;
        for(j=0; j<TmpEFMaxLen; ++j)
    {
      LastJoint = Cell->GetJoint(TmpEF[2*LocEdge+j]);

      if(!(CurrCell = LastJoint->GetNeighbour(Cell)))
        continue;

      CurrLocEdge = LastJoint->GetNeighbourEdgeIndex(Cell, LocEdge);

      while(CurrCell != Cell && CurrCell)
        {
          // Check if Edge is refined in CurrCell
          if(CurrCell->GetClipBoard() == -1)
      {
        CurrCell->SetClipBoard(0);
        N_ToRefine++;
      }

          // Get next joint which contains this edge
          CurrCell->GetShapeDesc()->GetEdgeFace(TmpEF2, TmpEF2MaxLen);
          if(CurrCell->GetJoint(TmpEF2[2*CurrLocEdge]) == LastJoint)
      Joint = CurrCell->GetJoint(TmpEF2[2*CurrLocEdge+1]);
          else
      Joint = CurrCell->GetJoint(TmpEF2[2*CurrLocEdge]);

          // Get new element and the index of our edge in this element
          CurrLocEdge = Joint->GetNeighbourEdgeIndex(CurrCell, CurrLocEdge);
          CurrCell = Joint->GetNeighbour(CurrCell);
          LastJoint = Joint;
        }
    } // j=0..N_JointsPerEdge
      } // i=0..N_edges
  }
    }

  while(N_ToRefine > 0)
    {
      TDatabase::IteratorDB[It_EQ]->Init(level);
      k=0;
      while ( (CurrCell = TDatabase::IteratorDB[It_EQ]->Next(info)) )
  {
    if(CurrCell->GetClipBoard() == 0)
      {
        //Output::print("Scanning Cell ", k);
        N_ToRefine += CurrCell->MakeConfClosure() - 1;
      }
    k++;
  }
    }
      
  k=-1;
  TDatabase::IteratorDB[It_EQ]->Init(level);
  while ((CurrCell = TDatabase::IteratorDB[It_EQ]->Next(info)))
    {
      k++;
      clip = CurrCell->GetClipBoard();

      switch(clip)
  {
    /*
     * Unique Refinements
     */
    // NoRef
  case 0:
    std::cerr << "This should not happen\n";
    break;
  case -1:
    CurrCell->SetNoRefinement();
    break;
    // Bisections
  case 1:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis0]);
    break;
  case 2:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis1]);
    break;
  case 4:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis2]);
    break;
  case 8:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis3]);
    break;
  case 16:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis4]);
    break;
  case 32:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis5]);
    break;
    // Diagonal Bisections
  case 33:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis05]);
    break;
  case 10:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis13]);
    break;
  case 20:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis24]);
    break;
    // QuadX
  case 7:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraQuad0]);
    break;
  case 25:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraQuad1]);
    break;
  case 50:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraQuad2]);
    break;
  case 44:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraQuad3]);
    break;
    // Reg
  case 63:
    CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraReg]);
    break;
  case 3:
  case 5:
  case 9:
  case 17:
  case 6:
  case 18:
  case 34:
  case 12:
  case 36:
  case 24:
  case 40:
  case 48:
    // Determine Local Face Index of face which is refined by double bisection
    switch(clip)
      {
      case 3: LocFace = 0; break;
      case 5: LocFace = 0; break;
      case 9: LocFace = 1; break;
      case 17: LocFace = 1; break;
      case 6: LocFace = 0; break;
      case 18: LocFace = 2; break;
      case 34: LocFace = 2; break;
      case 12: LocFace = 3; break;
      case 36: LocFace = 3; break;
      case 24: LocFace = 1; break;
      case 40: LocFace = 3; break;
      case 48: LocFace = 2; break;
      }

    // Check if Neighbour is already refined
    Joint = CurrCell->GetJoint(LocFace);
    Neigh = Joint->GetNeighbour(CurrCell);

    if(Neigh && Neigh->ExistChildren())
      {
        // Find Local Face Index at Neigh
        for(NeighLocFace=0; NeighLocFace<Neigh->GetN_Faces(); ++NeighLocFace)
    if(Joint == Neigh->GetJoint(NeighLocFace))
      break;

        if(NeighLocFace == Neigh->GetN_Faces())
    {
      std::cerr << "Face was not found at Neighbour\n";

    }

        NeighFaceRef = Neigh->GetRefDesc()->GetFaceRef(NeighLocFace);

        MapType = Joint->GetMapType();

        first_bis = (NeighFaceRef-TriBis01) / 2;
        second_bis = (NeighFaceRef-TriBis01) % 2;
        if(first_bis <= second_bis) second_bis++;

        neigh_first_bis = ((2-first_bis) + MapType) % 3;
        neigh_second_bis = ((2-second_bis) + MapType) % 3;

        neigh_type = 2*neigh_first_bis + neigh_second_bis;
        if(neigh_second_bis > neigh_first_bis) neigh_type--;

        MyFaceRef = Refinements(TriBis01 + neigh_type);

        switch(clip)
    {
      /*
       * Non-unique refinements
       */
      // Bis0X
    case 3:
      if(MyFaceRef == TriBis01)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis01]);
      else if(MyFaceRef == TriBis10)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis10]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
    case 5:
      if(MyFaceRef == TriBis02)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis02]);
      else if(MyFaceRef == TriBis20)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis20]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
    case 9:
      if(MyFaceRef == TriBis20)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis03]);
      else if(MyFaceRef == TriBis02)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis30]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
    case 17:
      if(MyFaceRef == TriBis21)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis04]);
      else if(MyFaceRef == TriBis12)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis40]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
      // Bis 1X
    case 6:
      if(MyFaceRef == TriBis12)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis12]);
      else if(MyFaceRef == TriBis21)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis21]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
    case 18:
      if(MyFaceRef == TriBis01)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis14]);
      else if(MyFaceRef == TriBis10)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis41]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
    case 34:
      if(MyFaceRef == TriBis02)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis15]);
      else if(MyFaceRef == TriBis20)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis51]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
      // Bis2X
    case 12:
      if(MyFaceRef == TriBis02)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis23]);
      else if(MyFaceRef == TriBis20)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis32]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
    case 36:
      if(MyFaceRef == TriBis01)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis25]);
      else if(MyFaceRef == TriBis10)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis52]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
      // Bis 3X
    case 24:
      if(MyFaceRef == TriBis01)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis34]);
      else if(MyFaceRef == TriBis10)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis43]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
    case 40:
      if(MyFaceRef == TriBis21)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis35]);
      else if(MyFaceRef == TriBis12)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis53]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
      // Bis 4X
    case 48:
      if(MyFaceRef == TriBis12)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis45]);
      else if(MyFaceRef == TriBis21)
        CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis54]);
      else
        std::cerr << "RefDesc does not fit to FaceRefDesc!\n";
      break;
    }
      }
    else
      {
        switch(clip)
    {
      /*
       * Non-unique refinements
       */
      // Bis0X
    case 3:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis01]);
      break;
    case 5:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis02]);
      break;
    case 9:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis03]);
      break;
    case 17:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis04]);
      break;
      // Bis 1X
    case 6:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis12]);
      break;
    case 18:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis14]);
      break;
    case 34:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis15]);
      break;
      // Bis2X
    case 12:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis23]);
      break;
    case 36:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis25]);
      break;
      // Bis 3X
    case 24:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis34]);
      break;
    case 40:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis35]);
      break;
      // Bis 4X
    case 48:
      CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBis45]);
      break;
    }
      }


    break;
  default:
    Output::print("Undefined Reference Description: ", CurrCell->GetClipBoard());
    break;
  }

      CurrCell->Refine(level+1);
    }
  return 0;
}

/**
 * This function is used to set the right refinement descriptions to each cell such that the grid is conforming afterwards.
 * After refinement of one ore more selected Cells this method has to be called.
 *
 * Max Winkler (23.02.2012)
 */
int TDomain::MakeConfClosure()
{
  int m;
  int MaxLevel;
  
  // Generate 1-regular grid
  Gen1RegGrid();

  TDatabase::IteratorDB[It_Finest]->Init(0);
  MaxLevel = TDatabase::IteratorDB[It_Finest]->GetMaxLevel();

  for(m=MaxLevel; m>=0; --m)
    {
      CloseGrid(m);
    }
  return 0;
}
#endif

#ifdef __2D__
int TDomain::Gen1RegGrid()
{
  int MaxLevel, CurrLevel, info;
  TBaseCell *CurrCell;

  TDatabase::IteratorDB[It_Finest]->Init(0);
  MaxLevel = TDatabase::IteratorDB[It_Finest]->GetMaxLevel();

  for (CurrLevel=MaxLevel;CurrLevel>0;CurrLevel--)
  {
    TDatabase::IteratorDB[It_EQ]->Init(CurrLevel);

    while ((CurrCell = TDatabase::IteratorDB[It_EQ]->Next(info)))
      if (!CurrCell->ExistChildren())
        CurrCell->Gen1RegGrid();
  }
  return 0;
}
#else
int TDomain::Gen1RegGrid()
{
  TBaseCell* CurrCell;
  int info;

  TDatabase::IteratorDB[It_Finest]->Init(0);
  while((CurrCell = TDatabase::IteratorDB[It_Finest]->Next(info)) )   
    CurrCell->Gen1RegGrid();

  return 0;
}
#endif

int TDomain::ConvertQuadToTri(int type)
{
  TBaseCell *CurrCell;
  Shapes CellType;
  int info;

  if (type)
  {
    RefLevel++;

    TDatabase::IteratorDB[It_Finest]->Init(0);

    // loop over all cells
    while ((CurrCell = TDatabase::IteratorDB[It_Finest]->Next(info)))
    {
      CellType = CurrCell->GetType();
      if (CellType == Parallelogram || CellType == Quadrangle ||
          CellType == Rectangle)
      {
        if (type == 1)
        {
          CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES +
                               QuadToTri0]);
        }
        else
        {
          CurrCell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES +
                               QuadToTri1]);
        }
        CurrCell->Refine(RefLevel);
      }
    }
  }

  return 0;
}

void TDomain::barycentric_refinement()
{
  RefLevel++;
  TDatabase::IteratorDB[It_Finest]->Init(0);
  int info;
  // loop over all cells
  while(auto current_cell = TDatabase::IteratorDB[It_Finest]->Next(info))
  {
    auto cell_type = current_cell->GetType();
    if(cell_type == Triangle)
    {
      current_cell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TriBary]);
    }
    else if(cell_type == Tetrahedron)
    {
      current_cell->SetRefDesc(TDatabase::RefDescDB[N_SHAPES + TetraBary]);
    }
    else
    {
      ErrThrow("unable to do barycentric refinement on a cell of type ",
               cell_type);
    }
    current_cell->Refine(RefLevel);
  }
  gridCollections.push_front(this->GetCollection(It_Finest, 0));
}


/*
 Extract a subcollection from a collection object:
*/
TCollection *TDomain::GetCollection(TCollection *coll, int reference)
{
  Output::root_info<2>("GetCollection", " Domain::GetCollection with reference: ", reference);
  int n_cells;
  TBaseCell **cells, *CurrCell;
  TCollection *subcoll;

  n_cells = 0;
  for (int i=0;i<coll->GetN_Cells();i++)
  {
    CurrCell = coll->GetCell(i);
    if (CurrCell->GetReference_ID()==reference)
    {
      n_cells++;
    }
  }

  // fill array of cells
  cells = new TBaseCell*[n_cells];
  int j = 0;
  for (int i=0;i<coll->GetN_Cells();i++)
  {
    CurrCell = coll->GetCell(i);
    if (CurrCell->GetReference_ID()==reference)
    {
      cells[j] = CurrCell;
      j++;
    }
  }
  Output::root_info<2>("Domain", "TDomain::GetCollection() creating collection, n_cells = ", n_cells);
  // create collection from an array of cells
  subcoll = new TCollection(n_cells, cells);
  
  return subcoll;
}


/*
 Extract a collection from a Domain object:
 For a given Iterator
   - count how many cells it contain
   - add cells to collection
*/
TCollection *TDomain::GetCollection(Iterators it, int level) const
{
  TCollection *coll;
  int i, n_cells, info;
  TBaseCell **cells, *CurrCell;

  // initialize the iterator
  TDatabase::IteratorDB[it]->Init(level);
  n_cells=0;
  // if the pointer to the next item is not empty increase number of cells
  // info: set to the current level
  while (TDatabase::IteratorDB[it]->Next(info)) n_cells++;

  // fill array of cells
  cells=new TBaseCell*[n_cells];
  TDatabase::IteratorDB[it]->Init(level);
  i=0;
#ifdef _MPI //the number of own cells in this collection must be counted
  int own_cell_counter = 0;
#endif
  while (((CurrCell = TDatabase::IteratorDB[it]->Next(info))))
  {
    cells[i]=CurrCell;
    cells[i]->SetCellIndex(i);
    i++;
#ifdef _MPI
    if (!CurrCell->IsHaloCell())
      ++own_cell_counter;
#endif
  }

  // create collection from an array of cells
  coll=new TCollection(n_cells, cells);

  #ifdef  _MPI
  coll->SetN_OwnCells(own_cell_counter);
  #endif

  return coll;
}

TCollection *TDomain::GetCollection(Iterators it, int level, int ID) const
{
  if(ID == -4711) 
    return this->GetCollection(it, level);
  
  int info;
  // initialize the iterator
  TDatabase::IteratorDB[it]->Init(level);
  
  //Output::print(" count cells ");
  int n_cells = 0;
  // loop over all cells, check their reference id
  // info: set to the current level
  while(TBaseCell *CurrCell = TDatabase::IteratorDB[it]->Next(info))
  {
    //Output::print(" id: ", CurrCell->GetPhase_ID());
    if (CurrCell->GetReference_ID() == ID)
      n_cells++;
  }
  
  //Output::print(" fill cells ");
  // fill array of cells
  TBaseCell **cells = new TBaseCell*[n_cells];
  int i=0;
  TDatabase::IteratorDB[it]->Init(level);
  while(TBaseCell *CurrCell = TDatabase::IteratorDB[it]->Next(info))
  {
    if (CurrCell->GetReference_ID() == ID)
    {
      cells[i] = CurrCell;
      cells[i]->SetCellIndex(i);
      i++;
    }
  }
  
  // create collection from an array of cells
  TCollection *coll = new TCollection(n_cells, cells);
  
  #ifdef  _MPI 
  coll->SetN_OwnCells(N_OwnCells);
  #endif
  
  return coll;
}

const std::list<TCollection*> TDomain::get_grid_collections() const
{
  if(gridCollections.empty())
    ErrThrow("There are no collections stored in the TDomain object, you "
             "need to call 'refine_and_get_hierarchy_of_collections'");
  return gridCollections;
}

void TDomain::DeRefine()
{
  CellTree[0]->Derefine();
}

#ifdef __3D__
static void nonRecursiveQuicksort (TVertex **array, int length)
{
  int stackSize = 0, leftEnd = 0, rightEnd = length-1;
  int i, j, m;
  TVertex *mid, *temp;

  const int initialStackLength = (int)( std::log((double)length)/std::log(2.0) + 1.0 );
  int stackLength = initialStackLength;
  std::vector<int> stack(stackLength);

  do
  {
    do
    {
      i = leftEnd;
      j = rightEnd;

      m = ( leftEnd + rightEnd ) / 2;
      mid = array[m];

      do
      {
        while(array[i] > mid)
          i++;

        while(array[j] < mid)
          j--;

        if (i<=j)
        {
          temp = array[i];
          array[i] = array[j];
          array[j] = temp;
          i++; j--;
        }
      } while (i<=j);

      if (leftEnd < j)
      {
        if(stackSize >= stackLength - 1)
        {
          stackLength += initialStackLength;
          stack.resize(stackLength);
        }
        
        stackSize++;
        stack[stackSize] = rightEnd;
        rightEnd = j;
      }
    } while (leftEnd < j);

    if(stackSize>0)
    {
      rightEnd = stack[stackSize];
      stackSize--;
    }

    leftEnd = i;

  } while(i < rightEnd);
}

static int GetIndex(TVertex **Array, int Length, TVertex *Element)
{
  int l=0, r=Length, m=(r+l)/2;
  TVertex *Mid;

  Mid=Array[m];
  while(Mid!=Element)
  {
    if(Mid>Element)
    {
      l=m;
    }
    else
    {
      r=m;
    }
    m=(r+l)/2;
    Mid=Array[m];
  }
  return m;
}

int TDomain::Grape(const char *name, TCollection *coll)
{
  int Format, Version, N_Elements, N_Vertices, N_LocVertices;
  int Header[5];
  int i,j,k,l,m, N_;
  int *Type;
  int *VertexNumbers, *NumberVertex;
  double *Coords;
  TBaseCell *cell;
  TVertex **Vertices, *Last, *Current;

  std::ofstream dat(name);
  if(!dat)
  {
    cerr << "cannot open file for output" << endl;
    return -1;
  }

  Format = 1;
  Version = 0;
  N_Elements = coll->GetN_Cells();

  Type = new int[N_Elements];

  N_LocVertices = 0;
  for(i=0;i<N_Elements;i++)
  {
    cell = coll->GetCell(i);
    N_LocVertices += cell->GetN_Vertices();
    switch(cell->GetType())
    {
      case Tetrahedron:
        Type[i] = 0;
      break;

      case Hexahedron:
      case Brick:
        Type[i] = 1;
      break;
      default:
        ErrThrow("TDomain::Grape only in 3D");
        break;
    }
  }

  Vertices = new TVertex*[N_LocVertices];
  N_ = 0;
  for(i=0;i<N_Elements;i++)
  {
    cell = coll->GetCell(i);
    k = cell->GetN_Vertices();
    for(j=0;j<k;j++)
      Vertices[N_++] = cell->GetVertex(j);
  }

  // sort the Vertices array
    nonRecursiveQuicksort (Vertices, N_);

  Last = nullptr;
  N_Vertices = 0;
  for(i=0;i<N_LocVertices;i++)
    if ((Current = Vertices[i]) != Last)
    {
      N_Vertices++;
      Last = Current;
    }

  Header[0] = Format;
  Header[1] = Version;
  Header[2] = N_Elements;
  Header[3] = N_Vertices;
  Header[4] = N_LocVertices;

  dat.write((char *)Header,sizeof(int)*5);
  dat.write((char *)Type,sizeof(int)*N_Elements);
  delete Type;

  Coords = new double[3*N_Vertices];
  VertexNumbers = new int[N_LocVertices];
  NumberVertex = new int[N_LocVertices];

  Last=nullptr;
  N_=0; k=-1;
  for(i=0;i<N_LocVertices;i++)
  {
    if((Current=Vertices[i])!=Last)
    {
      Vertices[i]->GetCoords(Coords[N_],Coords[N_+1],Coords[N_+2]);
      // Output::print(N_/2, "  ", Coords[N_], "  ", Coords[N_+1], "  ", Coords[N_+2]);
      k++;
      N_ += 3;
      Last=Current;
    }
    NumberVertex[i]=k;
  }

  m=0;
  for(i=0;i<N_Elements;i++)
  {
    cell = coll->GetCell(i);
    k = cell->GetN_Vertices();
    for(j=0;j<k;j++)
    {
      Current=cell->GetVertex(j);
      // Output::print((int)(Current));
      l=GetIndex(Vertices, N_LocVertices, Current);
      VertexNumbers[m]=NumberVertex[l];
      // Output::print(VertexNumbers[m]);
      m++;
    } // endfor j
  } //endfor i

  dat.write((char *)VertexNumbers,sizeof(int)*N_LocVertices);
  dat.write((char *)Coords,sizeof(double)*3*N_Vertices);
  
  delete Coords;
  delete NumberVertex;
  delete VertexNumbers;
  delete Vertices;

  return 0;
}

/** make boundary parameter consistent */
void TDomain::MakeBdParamsConsistent(TCollection *coll)
{
  int i,j,k;
  int N_Cells, N_Joints;
  double Param1[4], Param2[4];
  TBaseCell *cell;
  TJoint *joint;
  TBoundFace *boundface;
  const TBoundComp3D *bdcomp;
  const int *TmpFV, *TmpLen;
  int MaxLen;
  double x,y,z;

  N_Cells = coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);
    N_Joints = cell->GetN_Joints();
    cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
    for(j=0;j<N_Joints;j++)
    {
      joint = cell->GetJoint(j);
      if(joint->GetType() == BoundaryFace ||
                joint->GetType() == IsoBoundFace)
      {
        boundface = (TBoundFace *)joint;
        bdcomp = boundface->GetBoundComp();

        if(bdcomp->GetType() == Plane)
        {
          for(k=0;k<TmpLen[j];k++)
          {
            cell->GetVertex(TmpFV[j*MaxLen+k])->GetCoords(x,y,z);
            bdcomp->GetTSofXYZ(x,y,z,Param1[k], Param2[k]);
          } // endfor k
          boundface->SetParameters(Param1, Param2);
        } // endif Plane

        if(bdcomp->GetType() == Wall) 
        {
          if( ((TBdWall*)bdcomp)->GetBdComp2D()->GetType() == Line  )
          {
            for(k=0;k<TmpLen[j];k++)
            {
              cell->GetVertex(TmpFV[j*MaxLen+k])->GetCoords(x,y,z);
              bdcomp->GetTSofXYZ(x,y,z,Param1[k], Param2[k]);
            } // endfor k
            boundface->SetParameters(Param1, Param2);
          } // endif Line
        } // endif Wall

      } // endif BoundaryFace
    } // endfor j
  } // endfor i
} // MakeBdParamsConsistent

#endif



#ifdef __3D__
/** @brief added by sashi */
int TDomain::GenerateEdgeInfo()
{
  int i, ii, j, k, m, I, n_cells, info, N_RootVertices, *PointNeighb, MaxCpV = 0, EMaxLen;
  int level=0, N, NVert_All, *NumberVertex, *VertexNumbers, N_Edges, N_FaceEdges, N_FaceVert;
  int *RootVer_Loc, *NeibRootVer_Loc, a, b, len1, len2, cell_a, N_Neibs, *Neibs, *Neib_EdgeNo;
  int NeibN_Edges, *VertBeginIndex, Neiba, Neibb, BoundEdgeMarker[MAXN_EDGES3D], N_Faces, MaxLen;
  const int *EdgeVertex, *TmpFV, *TmpLen,  *TmpFE, *ETmpLen;;
 
  TBaseCell **cells, *CurrCell, *NeibCell, **NeibCells;
  Iterators it=It_Finest;
  TVertex **Vertices_All, *Last, *Current;
  const TShapeDesc *ShapeDesc, *NeibShapeDesc;
  TEdge *edge, *Neibedge;
  TJoint *joint;

#ifdef _MPI
  int rank;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);
  int EdgePhysicalID[MAXN_EDGES3D];
#endif

  NVert_All=0;
  n_cells=0;
  TDatabase::IteratorDB[it]->Init(level);
  while((CurrCell=TDatabase::IteratorDB[it]->Next(info)))
  {
    NVert_All +=CurrCell->GetN_Vertices();
    n_cells++;
  }

  cells=new TBaseCell*[n_cells];
  VertBeginIndex = new int[n_cells+1];
  NumberVertex = new int[NVert_All];
  VertexNumbers = new int[NVert_All];
  Vertices_All = new TVertex*[NVert_All];

  NVert_All=0;
  n_cells=0;
  VertBeginIndex[0] = 0;
  TDatabase::IteratorDB[it]->Init(level);
  while((CurrCell = TDatabase::IteratorDB[it]->Next(info)))
  {
    cells[n_cells]=CurrCell;
    N=CurrCell->GetN_Vertices();

    for(i=0;i<N;i++)
    {
      Vertices_All[NVert_All] = CurrCell->GetVertex(i);
      NVert_All++;
    }

    n_cells++;
    VertBeginIndex[n_cells] = NVert_All;

  } // while
// NVert_All--;
  // sort the Vertices array
    nonRecursiveQuicksort (Vertices_All, NVert_All);
// NVert_All++;
  Last=nullptr;
  N_RootVertices=-1;
  for(i=0;i<NVert_All;i++)
   {
    if((Current=Vertices_All[i])!=Last)
    {
      N_RootVertices++;
      Last=Current;
    }
    NumberVertex[i]=N_RootVertices;
   }
  N_RootVertices++;

  // Output::print("N_RootVertices ", N_RootVertices, " NVert_All ", NVert_All);

  m=0;
  for(ii=0;ii<n_cells;ii++)
   {
    CurrCell = cells[ii];
    N=CurrCell->GetN_Vertices();
    for(i=0;i<N;i++)
    {
      Current= CurrCell->GetVertex(i);
      I=GetIndex(Vertices_All, NVert_All, Current);
      VertexNumbers[m]=NumberVertex[I];
      // Output::print(VertexNumbers[m]);
      m++;
    } // endfor j
   } //while

  /** find max No. cells met any vertex in the collection **/
//    VertBound = new int[N_RootVertices];
   PointNeighb = new int[N_RootVertices];
   memset(PointNeighb, 0, N_RootVertices*sizeof(int));
//    memset(VertBound, 0, N_RootVertices*sizeof(int));

  for(i=0;i<NVert_All;i++)
    PointNeighb[VertexNumbers[i]]++;

  for(i=0;i<N_RootVertices;i++)
    if (PointNeighb[i] > MaxCpV) MaxCpV = PointNeighb[i];

   //Output::print("Max number of cells per vertex ", MaxCpV);
   /** No. of edge neibs will not exceed MaxCpV */ 
   Neibs = new int[MaxCpV];
   Neib_EdgeNo = new int[MaxCpV];
   NeibCells = new TBaseCell*[MaxCpV];

  /** PointNeighb's first column contains number of neib cells met with each vertex **/
  /** further columns contain the cell numbers met with this vertex **/
   MaxCpV++;

   delete [] PointNeighb;
   PointNeighb = new int[N_RootVertices*MaxCpV];
   memset(PointNeighb, 0, (N_RootVertices*MaxCpV)*sizeof(int));

  NVert_All=0;
  for(ii=0;ii<n_cells;ii++)
   {
    CurrCell = cells[ii];
    N=CurrCell->GetN_Vertices();
//     CurrCell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);

    for(i=0;i<N;i++)
     {
      I = VertexNumbers[NVert_All] *MaxCpV;
      PointNeighb[I]++;
      PointNeighb[I + PointNeighb[I] ] = ii;
      NVert_All ++;
     }

   } // for(ii=0;ii<n_cells;ii++)

// exit(0);

//   for(i=0;i<N_RootVertices;i++)
//     if ( VertBound[i]==0) printf("Vert %d VertBound %d \n", i, VertBound[i]);
//     //printf("Number of cells in vertex %d \n",  PointNeighb[i*MaxCpV]);
// 
// exit(0);
//   NVert_All=0;
  for(ii=0;ii<n_cells;ii++)
   {
    CurrCell = cells[ii];
    N=CurrCell->GetN_Vertices();
    N_Edges=CurrCell->GetN_Edges();
    N_Faces=CurrCell->GetN_Joints();
    
    ShapeDesc= CurrCell->GetShapeDesc();
    ShapeDesc->GetEdgeVertex(EdgeVertex);
    ShapeDesc->GetFaceEdge(TmpFE, ETmpLen, EMaxLen);
    ShapeDesc->GetFaceVertex(TmpFV, TmpLen, MaxLen);
    
    RootVer_Loc = VertexNumbers+VertBeginIndex[ii];
    
    
    /** set the marker for Edges */
    for(i=0;i<MAXN_EDGES3D;i++)
    {
      BoundEdgeMarker[i] = 0;
#ifdef _MPI
      EdgePhysicalID[i]  = 0;
#endif
    }
    
    /** set bound and isobound Edges */
    for(i=0;i<N_Faces;i++)
    {
      joint=CurrCell->GetJoint(i);

      if(joint->GetType() == BoundaryFace)
      {
#ifdef _MPI
          auto boundface              = (TBoundFace*)joint;
          const TBoundComp* BoundComp = boundface->GetBoundComp();
          int comp                    = BoundComp->get_physical_id();
#endif
        N_FaceEdges = ETmpLen[i];
        for(j=0;j<N_FaceEdges;j++) 
         {
          /** not modified if it is isobound edge */
          if(BoundEdgeMarker[TmpFE[i*EMaxLen+j] ] == 0)
          {
            BoundEdgeMarker[TmpFE[i*EMaxLen+j] ] = 1;
#ifdef _MPI
            EdgePhysicalID[TmpFE[i*EMaxLen+j] ] = comp;
#endif
          }
         }
         
         /** set all vertices on the bound face as bound vert */
         N_FaceVert = TmpLen[i];        
        for(j=0;j<N_FaceVert;j++)
        {
          CurrCell->GetVertex(TmpFV[i*MaxLen + j])->SetAsBoundVert();
#ifdef _MPI
          CurrCell->GetVertex(TmpFV[i*MaxLen + j])->set_physical_id(comp);
#endif
        }
          
      } // if(joint->GetType() == Boundar
      else if(joint->GetType() == IsoBoundFace)
      {
#ifdef _MPI
        auto boundface              = (TBoundFace*)joint;
        const TBoundComp* BoundComp = boundface->GetBoundComp();
        int comp                    = BoundComp->get_physical_id();
#endif
        N_FaceEdges = ETmpLen[i];
        for(j=0;j<N_FaceEdges;j++)
        {
          BoundEdgeMarker[TmpFE[i*EMaxLen+j] ] = 2;
#ifdef _MPI
          EdgePhysicalID[TmpFE[i*EMaxLen+j] ] = comp;
#endif
        }

        /** set all vertices on the bound face as bound vert */
        N_FaceVert = TmpLen[i];        
        for(j=0;j<N_FaceVert;j++)
        {
#ifdef _MPI
          CurrCell->GetVertex(TmpFV[i*MaxLen + j])->SetAsBoundVert();
          CurrCell->GetVertex(TmpFV[i*MaxLen + j])->set_physical_id(comp);
#endif
        }
      } // 
     } // for(i=0;i<N_    
    
    for(i=0;i<N_Edges;i++)
     {
      N_Neibs = 0;
      edge = CurrCell->GetEdge(i);
      
      if(edge==nullptr) // edge not yet created
       {
        a = RootVer_Loc[EdgeVertex[2*i]];
        b = RootVer_Loc[EdgeVertex[2*i+1]];

        len1 = PointNeighb[a*MaxCpV]; // No. cells met with vert a
        len2 = PointNeighb[b*MaxCpV]; // No. cells met with vert b

        /** find No. cells having this edge */
        for(j=1;j<=len1;j++)
         {
          cell_a = PointNeighb[a*MaxCpV + j];
// 	  if(rank==0 && cell_a==16)
//            Output::print(" cell_a ", cell_a);

          for(k=1;k<=len2;k++)
           {

            if(cell_a == PointNeighb[b*MaxCpV +k] )
             {  
              Neibs[N_Neibs] = cell_a;
              N_Neibs ++;
              break;
             }
           } // for(k=1;k<=
         } // for(j=1;j<

        /** find the local index of the edge in all cells */
        for(j=0;j<N_Neibs;j++)
         {

          if(Neibs[j]==ii)
           {
             Neib_EdgeNo[j] = i;
           }
          else 
           {
            NeibCell = cells[Neibs[j]];
            NeibShapeDesc= NeibCell->GetShapeDesc();
            NeibN_Edges=NeibCell->GetN_Edges();
            NeibShapeDesc->GetEdgeVertex(EdgeVertex);
            NeibRootVer_Loc = VertexNumbers+VertBeginIndex[ Neibs[j]];
            //Output::print("  Neib ", Neibs[j]);

            for(k=0;k<NeibN_Edges;k++)
             {
              Neibedge = NeibCell->GetEdge(k);

              if(Neibedge==nullptr) // edge not yet created
               {
                 Neiba = NeibRootVer_Loc[EdgeVertex[2*k]];
                 Neibb = NeibRootVer_Loc[EdgeVertex[2*k+1]];

                 if( (Neiba==a &&  Neibb==b) || (Neiba==b &&  Neibb==a) )
                  {
                   //Output::print("  Neib edge found ");
                   Neib_EdgeNo[j] = k;
                   break;
                  }
               }
             } //for(k=0;k<NeibN_Edges;k

            if(k==NeibN_Edges)
             {
              ErrThrow("Error could not find edge index in the Neib cell !");
             }
            }
         } // for(j=0;j<N_Neib

       /**generate the edge and set in all cells */
       for(j=0;j<N_Neibs;j++)
        {
         NeibCells[j] = cells[Neibs[j]];
         NeibCells[j]->SetEdge(Neib_EdgeNo[j], edge);
        }

//   	  if(rank==0 && ii==35 && i==1)      
// 	  {
//            for(j=0;j<N_Neibs;j++)
// 	     printf("rank %d j %d NeibCells %d BoundEdgeMarker %d \n", rank, j, Neibs[j], BoundEdgeMarker[i]); 
// 	     
// 	  }
//         
        
        
       if(BoundEdgeMarker[i]==2 )
        {
         edge = new TIsoEdge3D(N_Neibs, NeibCells);
        }
       else if(BoundEdgeMarker[i]==1)
        {
         edge = new TBDEdge3D(N_Neibs, NeibCells);
        }
       else
        {
         edge = new TInnerEdge(N_Neibs, NeibCells);
        }
#ifdef _MPI
       edge->set_physical_id(EdgePhysicalID[i]);
#endif
       } //if(edge==nullptr)
//      else
//       {
//          ErrThrow(ii, " Edge already set ", i);
//       }

       for(j=0;j<N_Neibs;j++)
         NeibCells[j]->SetEdge(Neib_EdgeNo[j], edge);

     } // for(i=0;i<N_Edges;i++)
   } //  for(ii=0;ii<n_cel

   delete [] Neibs;
   delete [] Neib_EdgeNo;
   delete [] NeibCells;
   delete [] cells;
   delete [] VertBeginIndex;
   delete [] NumberVertex;
   delete [] VertexNumbers;
   delete [] Vertices_All;
//    delete [] VertBound;
   delete [] PointNeighb;

#ifdef _MPI
   if(rank==0)
#endif
    Output::info<5>("Domain.C","3D Mesh Edges Generated ");

  return 0;
}


#endif

size_t TDomain::get_n_initial_refinement_steps() const
{
  return db["refinement_n_initial_steps"];
}

size_t TDomain::get_max_n_adaptive_steps() const
{
  return db["refinement_max_n_adaptive_steps"];
}

void TDomain::print_info(const std::string& name) const
{
#ifdef _MPI
  int my_rank, size;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &my_rank);
  MPI_Comm_size(TDatabase::ParamDB->Comm, &size);

  int sbuf_own = N_OwnCells;
  int sbuf_halo = N_RootCells - N_OwnCells;

  std::vector<int> ns_own_cells(size,0);
  std::vector<int> ns_halo_cells(size,0);

  {
    MPI_Gather(
      &sbuf_own, 1, MPI_INT,            //send
      &ns_own_cells.at(0), 1, MPI_INT,  //receive
      0, MPI_COMM_WORLD);               //control
    MPI_Gather(
      &sbuf_halo, 1, MPI_INT,            //send
      &ns_halo_cells.at(0), 1, MPI_INT,  //receive
      0, MPI_COMM_WORLD);               //control
  }
  if(my_rank == 0)
  {
    Output::stat("Domain", name);
    size_t sum_cells_total = 0;
    for(int i =0; i < size ;++i)
    {
      Output::dash("Process", i, "\t n_own_cells: ", ns_own_cells.at(i),
                    "\t n_halo_cells: ", ns_halo_cells.at(i));
      sum_cells_total += ns_own_cells.at(i);
    }
    Output::dash("Total number of cells: ", sum_cells_total);

  }

#else
  Output::stat("Domain ", name);
  Output::dash("No domain statistics printout in non-MPI case so far.");
#endif
}


void determine_n_refinement_steps_multigrid(
  const std::string& multigrid_type,
  int n_multigrid_levels,
  int n_initial_refinement_steps,
  int& n_ref_before, int& n_ref_after)
{
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#else
  int my_rank=0;
#endif
  //Check if input multigrid type makes sense.
  if(multigrid_type!=std::string("standard") && multigrid_type!=std::string("mdml"))
    ErrThrow("Unknown multigrid type ", multigrid_type,". Use 'standard' or 'mdml'.");
  //Check if input numbers make sense
  if(n_multigrid_levels > n_initial_refinement_steps + 1)
    ErrThrow("You requested ", n_multigrid_levels, " geometric levels for "
        "multigrid, but only ", n_initial_refinement_steps + 1, " grid levels "
        "will be present after refining ", n_initial_refinement_steps, " times. "
        " Choose a smaller multigrid_n_levels.");

  //split the number of refinement steps
  n_ref_after =  n_multigrid_levels - 1;
  n_ref_before =  n_initial_refinement_steps - n_ref_after;

  //print information
  if(my_rank == 0)
  {
    Output::info("REFINEMENT FOR MULTIGRID","You are using ",multigrid_type,
                 " multigrid with ",n_multigrid_levels," geometric levels. ");
#ifdef _MPI
    Output::dash("Number of refinement steps before domain partitioning: ", n_ref_before);
    Output::dash("Number of refinement steps after domain partitioning: ", n_ref_after);
    Output::dash("Total number of refinement steps: ", n_initial_refinement_steps);
#else
    Output::dash("Number of refinement steps before picking grids: ", n_ref_before);
    Output::dash("Number of refinement steps after picking grids: ", n_ref_after);
    Output::dash("Total number of refinement steps: ", n_initial_refinement_steps);
#endif
  }
}


std::list<TCollection*> TDomain::refine_and_get_hierarchy_of_collections(
                                      const ParameterDatabase& parmoon_db,
#ifdef __2D__
                                      BoundCondFunct2D* BoundaryCondition)
#else
                                      BoundCondFunct3D* BoundaryCondition)
#endif
{
#ifdef _MPI
  int my_rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  gridCollections.clear();
  gridCollections.push_front(this->GetCollection(It_Finest, 0));

  // Split the number of refinement steps - see doc of the method
  // TODO Removing this strange construction is a TODO
  int n_ref_before = this->get_n_initial_refinement_steps();
  int n_ref_after = 0;
  int n_steps_before_partition = db["refinement_n_steps_before_partitioning"];
  if(n_steps_before_partition < n_ref_before)
  {
    n_ref_after = n_ref_before - n_steps_before_partition;
    n_ref_before = n_steps_before_partition;
  }
  if(parmoon_db.contains("preconditioner")
     && parmoon_db["preconditioner"].is(std::string("multigrid"))
     && parmoon_db.contains("solver_type")
     && parmoon_db["solver_type"].is(std::string("iterative")))
  {
    determine_n_refinement_steps_multigrid(
      parmoon_db["multigrid_type"],
      parmoon_db["multigrid_n_levels"],
      this->get_n_initial_refinement_steps(),
      n_ref_before, n_ref_after);
  }
  bool final_barycentric = this->db["refinement_final_step_barycentric"];
  for(int i = 0; i < n_ref_before; i++)
  {
    if(final_barycentric && i+1 == n_ref_before && n_ref_after == 0)
      this->barycentric_refinement();
    else
      this->RegRefineAll();
  }


  // set periodic boundaries, currently only implemented with 3D geometry
  if(BoundaryCondition)
  {
#ifdef _MPI
    if(parmoon_db.contains("preconditioner")
    && parmoon_db["preconditioner"].is(std::string("multigrid")))
    {
      ErrThrow("Periodic boundaries are not tested for MPI with multigrid yet");
    }
#endif
    set_PeriodicJoint(BoundaryCondition, this);
  }


#ifdef _MPI
  // Partition the by now finest grid using Metis and distribute among processes.

  // 1st step: Analyse interfaces and create edge objects,.
  this->GenerateEdgeInfo();

  // 2nd step: Call the mesh partitioning.

  int maxCellsPerVertex;

  //do the actual partitioning, and examine the return value
  if(Partition_Mesh3D(MPI_COMM_WORLD, this, TDatabase::ParamDB->Par_P2,
                      maxCellsPerVertex) == 1)
  {
    ErrThrow("Partitioning did not succeed.");
  }

  // 3rd step: Generate edge info anew
  //(since distributing changed the domain).
  this->GenerateEdgeInfo();

  // with _MPI the coarser grids are no longer used
  gridCollections.clear();
  gridCollections.push_front(this->GetCollection(It_Finest, 0));

#endif

  
  

  //this is only relevant for multigrid
  for(int level=0; level < n_ref_after; ++level)
  {
    if(final_barycentric && level+1 == n_ref_after)
      this->barycentric_refinement();
    else
      this->RegRefineAll();

#ifdef _MPI
    // remove new collection from 'gridCollections', because we need to call two 
    // methods before, then insert it again
    gridCollections.pop_front();
    this->GenerateEdgeInfo();  // has to be called anew after every refinement step
    Domain_Crop(MPI_COMM_WORLD, this); // remove unwanted cells in the halo after refinement
    gridCollections.push_front(this->GetCollection(It_Finest, 0));
#endif
    // Grab collection.
  }
  return gridCollections;
}


#ifdef __3D__


void TDomain::GenerateFromMesh(Mesh& m)
{
  // create inner faces (if needed)
  m.createInnerFaces();

  /** 
      compute number boundary faces:
      This step at the moment is necessary to initialize some useful arrays
      (e.g. boundary faces)
  */
  m.computeNumberOfBoundaryFaces();

  if (db["fix_boundary_elements"])
  {
    m.fixBoundaryElements();

#ifdef _MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
#endif
      if (db["fix_boundary_elements_write_mesh"])
      {
        m.writeToMesh(db["fix_boundary_elements_write_mesh_file"], false);
      }
#ifdef _MPI
    }
#endif
  }

  // build the boundary parts
  this->buildBoundary(m);

  // build the mesh
  this->buildParMooNMesh(m);

  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(this);
  TDatabase::IteratorDB[It_LE]->SetParam(this);
  TDatabase::IteratorDB[It_Finest]->SetParam(this);
  TDatabase::IteratorDB[It_Between]->SetParam(this);
  TDatabase::IteratorDB[It_OCAF]->SetParam(this);  
}

// note this function needs:
// list of faces (triangles or quads)
// number of boundary faces
// adjacency (1,-1 for boundary faces)
// list of points
///@todo write part of the function inside the Mesh class, here set only the BdParts
void TDomain::buildBoundary(Mesh& m)
{
  // we consider only 1 boundary part
  ///@todo do we need to take into account the general case?
  this->N_BoundParts = 1;

  // vector of boundary parts
  this->BdParts = new TBoundPart* [this->N_BoundParts];

  // the number of boundary components is computed from the adjacency
  this->N_BoundComps = m.n_boundary_faces;
  meshBoundComps.resize(this->N_BoundComps);

  Output::root_info("Domain", "TDomain::buildBoundary() - N_BoundComps: ", this->N_BoundComps);
  this->BdParts[0] = new TBoundPart(this->N_BoundComps);

  // StartBdCompID[i] gives the index of the element in BdParts
  // where the BdComp i starts. The last element is equal to to N_BoundParts
  this->StartBdCompID = new int[this->N_BoundParts + 1];

  // we now consider only a boundary components,
  // i.e. StartBdCompID[0] = 0, StartBdCompID[1] = N_BoundComps
  ///@todo check how we can extend this to multiple boundary markers
  this->StartBdCompID[0] = 0;
  this->StartBdCompID[1] = this->N_BoundComps;

  const std::vector<std::vector<int>>* ptr_faceToCell;
  int nb_face;

  if (m.hexa.size() == 0)  // mesh with only tetrahedra
  {
    ptr_faceToCell = &m.faceToTetra;
    nb_face = m.triangle.size();
  }
  else if (m.tetra.size() == 0)  // mesh with only hexahedra
  {
    ptr_faceToCell = &m.faceToHexa;
    nb_face = m.quad.size();
  }
  else
  {
    ErrThrow(" *** Error(TDomain::buildBoundary)");
  }

  // in order to describe the boundary, each boundary face is considered
  // to be a boundary component of the (single) boundary part.
  // Hence,
  // (*) first, we check if the face is on the boundary
  // (*) second, given the vertices, we define a TBdPlane describing the face

  int counter = 0; // count the number of boundary faces found

  if (m.boundaryFacesMarker.size() == 0)
  {
    m.boundaryFacesMarker.resize(nb_face);
  }

  for (int i = 0; i < nb_face; ++i)
  {
    if (ptr_faceToCell->at(i)[0] < 0 || ptr_faceToCell->at(i)[1] < 0)
    {
      double p[3], a[3], b[3], n[3];

      if (m.is_tetramesh)
      {
        // the face is on the boundary
        meshBoundComps.at(counter) = new TBdPlane(counter,
                                                  m.triangle[i].reference);

        /// @attention the lists of nodes starts from 1 (not from 0)
        p[0] = m.vertex[ m.triangle[i].nodes[0] - 1].x;
        a[0] = m.vertex[ m.triangle[i].nodes[1] - 1].x - p[0];
        b[0] = m.vertex[ m.triangle[i].nodes[2] - 1].x - p[0];

        p[1] = m.vertex[ m.triangle[i].nodes[0] - 1].y;
        a[1] = m.vertex[ m.triangle[i].nodes[1] - 1].y - p[1];
        b[1] = m.vertex[ m.triangle[i].nodes[2] - 1].y - p[1];

        p[2] = m.vertex[ m.triangle[i].nodes[0] - 1].z;
        a[2] = m.vertex[ m.triangle[i].nodes[1] - 1].z - p[2];
        b[2] = m.vertex[ m.triangle[i].nodes[2] - 1].z - p[2];
      }
      else //if(m.is_hexamesh)
      {
        // the face is on the boundary
        meshBoundComps.at(counter) = new TBdPlane(counter,
                                                  m.quad[i].reference);

        /// @attention the lists of nodes starts from 1 (not from 0)
        p[0] = m.vertex[ m.quad[i].nodes[0] - 1].x;
        a[0] = m.vertex[ m.quad[i].nodes[1] - 1].x - p[0];
        b[0] = m.vertex[ m.quad[i].nodes[3] - 1].x - p[0];

        p[1] = m.vertex[ m.quad[i].nodes[0] - 1].y;
        a[1] = m.vertex[ m.quad[i].nodes[1] - 1].y - p[1];
        b[1] = m.vertex[ m.quad[i].nodes[3] - 1].y - p[1];

        p[2] = m.vertex[ m.quad[i].nodes[0] - 1].z;
        a[2] = m.vertex[ m.quad[i].nodes[1] - 1].z - p[2];
        b[2] = m.vertex[ m.quad[i].nodes[3] - 1].z - p[2];
      }

      // normalize vector a
      double fac = std::sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
      a[0] /= fac;
      a[1] /= fac;
      a[2] /= fac;

      // normalize vector b
      fac = std::sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
      b[0] /= fac;
      b[1] /= fac;
      b[2] /= fac;

      // cross product of a and b
      n[0] = a[1]*b[2] - a[2]*b[1];
      n[1] = a[2]*b[0] - a[0]*b[2];
      n[2] = a[0]*b[1] - a[1]*b[0];

      // normalize n
      fac = std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
      n[0] /= fac;
      n[1] /= fac;
      n[2] /= fac;

      ((TBdPlane*) meshBoundComps[counter])->SetParams(p[0], p[1], p[2],
                                                       a[0], a[1], a[2],
                                                       n[0], n[1], n[2]);

      ++counter;

      // the id of the boundary face is taken as the attribute of the triangle
      // note: later, trifacemarkerlist=0 is used to identify inner faces (see below)
      m.boundaryFacesMarker[i] = counter;

      if(m.is_tetramesh)
      {
        Output::print<5>("TDomain::buildBoundary() Triangle ", i,
                         " boundaryMarker = ", m.boundaryFacesMarker[i],
                         " attribute = ", m.triangle[i].reference);
      }
      else if(m.is_hexamesh)
      {
        Output::print<5>("TDomain::buildBoundary() Quad ", i,
                         " boundaryMarker = ", m.boundaryFacesMarker[i],
                         " attribute = ", m.quad[i].reference);
      }
    }
    else // not a boundary face
    {
      m.boundaryFacesMarker[i] = 0;
    }
  }

  // check if the number of boundary components is consistent
  assert(counter==N_BoundComps);
  // initialize ParMooN BdParts and the BdComp
  for(int i=0; i<this->N_BoundComps; i++)
    this->BdParts[0]->SetBdComp(i, meshBoundComps[i]);
}


void TDomain::buildParMooNMesh(Mesh& m)
{
  // create a vector<TVertex*> from the list of pointa
  this->setVertices(m);

  // allocate memory for CellTree and set the vertices
  this->allocRootCells(m);

  // create joints and set joints type
  this->distributeJoints(m);
}

/**
 * @brief set the coordinates of the vertices
 *
 * This functions reads the vertices coordinates 
 * and creates TVertex*, stored in the vector meshVertices
 * Moreover, bounds for the domain are computed
 */
void TDomain::setVertices(Mesh& m)
{
 
  double xmin=0, xmax=0, ymin=0, ymax=0, zmin=0, zmax=0;
  meshVertices.resize(m.vertex.size()); 
  Output::root_info("Domain", "TDomain::setVertices() total number of vertices: ",  meshVertices.size());

  for(unsigned int i=0;i<meshVertices.size();++i)
  {
    double x = m.vertex[i].x;
    double y = m.vertex[i].y;
    double z = m.vertex[i].z;

    if(i > 0)
    {
      if(xmin > x) xmin = x;
      if(xmax < x) xmax = x;

      if(ymin > y) ymin = y;
      if(ymax < y) ymax = y;

      if(zmin > z) zmin = z;
      if(zmax < z) zmax = z;
    }
    else
    {
      xmin = xmax = x;
      ymin = ymax = y;
      zmin = zmax = z;
    }

    meshVertices.at(i)=new TVertex (x, y, z);
  }

  StartX = xmin;
  BoundX = xmax - xmin;

  StartY = ymin;
  BoundY = ymax - ymin;

  StartZ = zmin;
  BoundZ = zmax - zmin;

}

void TDomain::allocRootCells(Mesh& m)
{
  TVertex *Vertex;
  TMacroCell *Cell;

  if (m.hexa.size()!=0 && m.tetra.size()!=0)
  {
    ErrThrow(" *** Error(TDomain::allocRootCells)");
  }

  if (m.is_tetramesh)
  {
    // set the number of total cells and allocate the cell array
    this->N_RootCells = m.tetra.size();

    Output::root_info<3>("Domain", "TDomain::allocRootCells() number of tetrahedra: ",
                     this->N_RootCells);

    this->CellTree = new TBaseCell*[this->N_RootCells];

    // set the vertices of each cell, reading from meshVertices
    for (int i = 0; i < this->N_RootCells; ++i)
    {
      Cell = new TMacroCell(TDatabase::RefDescDB[Tetrahedron], 0);
      this->CellTree[i] = Cell;

      for (int j = 0; j < 4; ++j)
      {
        Vertex = meshVertices.at(m.tetra[i].nodes[j] - 1);
        Cell->SetVertex(j, Vertex);
      }

      Cell->SetPhase_ID((int) m.tetra[i].reference);
    }
  }
  else if(m.is_hexamesh)
  {
    // set the number of total cells and allocate the cell array
    this->N_RootCells = m.hexa.size();
    Output::print<3>("TDomain::allocRootCells() number of hexahedra: ",
                     this->N_RootCells);

    this->CellTree = new TBaseCell* [this->N_RootCells];
    // set the vertices of each cell, reading from meshVertices
    for(int i=0;i<this->N_RootCells;++i)
    {
      Cell = new TMacroCell (TDatabase::RefDescDB[Hexahedron], 0);
      this->CellTree[i] = Cell;

      for(int j=0;j<8;++j)
      {
        Vertex = meshVertices.at(m.hexa[i].nodes[j]-1);
        Cell->SetVertex(j, Vertex);
      }
      Cell->SetPhase_ID((int) m.hexa[i].reference);
    }
  }

}

///@attention the functions hashTriFaces() or hashQuadFaces(),
///           and CreateAdjacency() must have been called before
void TDomain::distributeJoints(Mesh& m)
{
  if(m.hexa.size()!=0 && m.tetra.size()!=0)
  {
    ErrThrow(" *** Error(TDomain::allocRootCells)");
  }

  if(m.is_tetramesh)
  {
  // allocate the joints array
  meshJoints.resize(m.triangle.size());
  Output::root_info("Domain", "number of joints: ", meshJoints.size());

  // set the faces for each joint
  for(unsigned int i=0; i<meshJoints.size(); ++i)
  {
    // find element that contain this face
    if(m.boundaryFacesMarker[i] == 0)// inner joints
    {
      int left = m.faceToTetra[i][0];
      int right = m.faceToTetra[i][1];
      Output::print<5>("TDomain::distributeJoints() Joint ", i ," Tetra: ", left, " and ", right);
      meshJoints.at(i) = new TJointEqN (CellTree[left], CellTree[right]);
    }
    else // boundary joints
    {
      Output::print<5>("Joint ", i, " is a boundary joint");
      int bdcomp = m.boundaryFacesMarker[i] - 1;
      TBoundComp3D* BoundComp = meshBoundComps[bdcomp];
      meshJoints.at(i)=new TBoundFace (BoundComp);
    }
  }

  
  for(unsigned int i=0;i<m.tetra.size();++i)
  {
    const int *TmpFV, *TmpLen;
    int MaxLen;
    
    const TShapeDesc *ShapeDesc = CellTree[i]->GetShapeDesc();
    ShapeDesc->GetFaceVertex(TmpFV, TmpLen, MaxLen);

    for(int j=0;j<4;++j)
    {
      ///@attention make sure that the vector meshTrifaceHash has been created 
      int triface = m.findTriFace(m.tetra[i].nodes[TmpFV[j*MaxLen  ]],
				  m.tetra[i].nodes[TmpFV[j*MaxLen+1]],
				  m.tetra[i].nodes[TmpFV[j*MaxLen+2]]);

      // check that a face has been found
      assert (triface != -1);

      CellTree[i]->SetJoint(j, meshJoints.at(triface));

      // correct parameters for boundary faces
      if(meshJoints.at(triface)->GetType() == BoundaryFace)
      {

        TBoundFace *BoundFace = (TBoundFace*) meshJoints.at(triface);
        const TBoundComp3D *BoundComp = BoundFace->GetBoundComp();
        double param1[4], param2[4];

        for(int k=0;k<TmpLen[j];++k)
        {
	  double x, y, z, t, s;
	  CellTree[i]->GetVertex(TmpFV[MaxLen*j+k])->GetCoords(x,y,z);
          BoundComp->GetTSofXYZ(x,y,z, t, s);

          param1[k] = t;
          param2[k] = s;
        }
        BoundFace->SetParameters(param1, param2);
      }
    }
  }
  
  // set map type
  for(unsigned int i=0;i<meshJoints.size();++i)
  {
    meshJoints[i]->SetMapType();
  }
  }
  else if(m.is_hexamesh)
  {
    // allocate the joints array
    meshJoints.resize(m.quad.size());
    Output::print("number of joints: ", meshJoints.size());

    // set the faces for each joint
    for(unsigned int i=0; i<meshJoints.size(); ++i)
    {
      // find element that contain this face
      if(m.boundaryFacesMarker[i] == 0)// inner joints
      {
        int left = m.faceToHexa[i][0];
        int right = m.faceToHexa[i][1];
        Output::print<5>("TDomain::distributeJoints() Joint ", i ,
                         " Hexa: ", left, " and ", right);
        meshJoints.at(i) = new TJointEqN (CellTree[left], CellTree[right]);
      }
      else // boundary joints
      {
        Output::print<5>("Joint ", i, " is a boundary joint");
        int bdcomp = m.boundaryFacesMarker[i] - 1;
        TBoundComp3D* BoundComp = meshBoundComps[bdcomp];
        meshJoints.at(i)=new TBoundFace (BoundComp);
      }
    }

    for(unsigned int i=0;i<m.hexa.size();++i)
    {
      const int *TmpFV, *TmpLen;
      int MaxLen;

      const TShapeDesc *ShapeDesc = CellTree[i]->GetShapeDesc();
      ShapeDesc->GetFaceVertex(TmpFV, TmpLen, MaxLen);

      for(int j=0;j<6;++j)
      {
        ///@attention make sure that the vector meshQuadfaceHash has been created 
        int quadface = m.findQuadFace(m.hexa[i].nodes[TmpFV[j*MaxLen  ]],
                                      m.hexa[i].nodes[TmpFV[j*MaxLen+1]],
                                      m.hexa[i].nodes[TmpFV[j*MaxLen+2]],
                                      m.hexa[i].nodes[TmpFV[j*MaxLen+3]]);

        // check that a face has been found
        assert (quadface != -1);

        CellTree[i]->SetJoint(j, meshJoints.at(quadface));

        // correct parameters for boundary faces
        if(meshJoints.at(quadface)->GetType() == BoundaryFace)
        {

          TBoundFace *BoundFace = (TBoundFace*) meshJoints.at(quadface);
          const TBoundComp3D *BoundComp = BoundFace->GetBoundComp();
          double param1[4], param2[4];

          for(int k=0;k<TmpLen[j];++k)
          {
            double x, y, z, t, s;
            CellTree[i]->GetVertex(TmpFV[MaxLen*j+k])->GetCoords(x,y,z);
            BoundComp->GetTSofXYZ(x,y,z, t, s);

            param1[k] = t;
            param2[k] = s;
          }
          BoundFace->SetParameters(param1, param2);
        }
      }
    }

    // set map type
    for(unsigned int i=0;i<meshJoints.size();++i)
    {
      meshJoints[i]->SetMapType();
    }
  }
}
#endif // 3D


#ifdef  __3D__
bool TDomain::check() const
{
	   for(int cell_id=0;cell_id<N_RootCells;cell_id++)
	   {
		   TBaseCell *CurrCell = CellTree[cell_id];

		   //check cells
		   if(!CurrCell->check_orientation())
			   ErrThrow("Cell ", cell_id," does not meet the right hand rule.");
		   if(!CurrCell->check_shape())
			   ErrThrow("Cell ", cell_id," is not of the given shape.");

		   int n_joints = CurrCell->GetN_Joints();

	     for(int joint_id=0;joint_id<n_joints;joint_id++)
	     {
	    	 TJoint* CurrJoint = CurrCell->GetJoint(joint_id);

    		 const int *TmpFV, *TmpLen;
    		 int MaxLen;
    		 CurrCell->GetRefDesc()->GetShapeDesc()
    	         ->GetFaceVertex(TmpFV, TmpLen, MaxLen);

	    	 if ( !CurrJoint )
	    		 Output::print("Some Joint is not set: Cell ",cell_id, ", Joint ",joint_id);

	    	 //Check outer joints
	    	 else if(!CurrJoint->InnerJoint())
	    	 {
	    		 TBoundFace * CurrFace = dynamic_cast<TBoundFace *>(CurrCell->GetJoint(joint_id));
	    		 double Param1[4];
	    		 double Param2[4];

	    	     for(int vert_id=0;vert_id<TmpLen[joint_id];vert_id++)
	    	     {
	    	       double X , Y , Z , xp , yp , zp;
	    	       CurrFace->GetParameters(Param1, Param2);
	    	       TVertex *Vert = CurrCell->GetVertex(TmpFV[joint_id*MaxLen+vert_id]);
	    	       Vert->GetCoords(X, Y, Z);
	    	       CurrFace->GetXYZofTS(Param1[vert_id], Param2[vert_id], xp, yp, zp);
	    	       //Output::print("Coordinates :", X, " ",Y," ",Z);
	    	       //Output::print("Parametrization Coordinates :", xp, " ",yp," ",zp);
	    	       if(-1e-10>X-xp || 1e-10<X-xp)
	    	    	   ErrThrow("Error in parametrization in Cell ",cell_id, ", Joint ",joint_id);
	    	       if(-1e-10>Y-yp || 1e-10<Y-yp)
	    	    	   ErrThrow("Error in parametrization in Cell ",cell_id, ", Joint ",joint_id);
	    	       if(-1e-10>Z-zp || 1e-10<Z-zp)
	    	    	   ErrThrow("Error in parametrization in Cell ",cell_id, ", Joint ",joint_id);
	    	     }
	    	 }//end outer joints

	    	 //check inner joints
	    	 else
	    	 {
	    		 TBaseCell * NeighCell = CurrJoint->GetNeighbour(CurrCell);

	    		 int neighjoint_id=0;
	    		 for(int joint=0; joint<n_joints;joint++)
	    		 {
	    			 if( NeighCell->GetJoint(joint)==CurrJoint )
	    			 {
	    				 neighjoint_id=joint;
	    				 break;
	    			 }
	    		 }

	    		 for(int vert_id=0;vert_id<TmpLen[joint_id];vert_id++)
	    		 {
		    		 bool vert_match= false;

	    			 for(int neighvert_id=0; neighvert_id<TmpLen[joint_id]; neighvert_id++)
	    				 if(CurrCell->GetVertex(TmpFV[joint_id*MaxLen+vert_id])
	    						 ==NeighCell->GetVertex(TmpFV[neighjoint_id*MaxLen+neighvert_id]))
	    					 vert_match=true;

	    			 if(!vert_match)
	    				 ErrThrow("Some Joint does not match: Cell ",cell_id, ", Joint ",joint_id);
	    		 }

	    	 }//end inner joints

	     }
	   }

	return true;
}

#endif
//END DEBUG

