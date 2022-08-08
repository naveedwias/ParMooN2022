/** ************************************************************************ 
*
* @class TDomain  
* @date  09.07.97
* @brief  contains the boundary description, the virtual cell tree and macro grid
* @author Volker Behns
* @History: collection methods (Gunar Matthies 14.10.97), 
            methods for Refine/Derefine algorithm  (Gunar Matthies 17.10.97)
            mesh class, edge generation, parallel methods and documentation (Sashikumaar Ganesan 08.08.2014)
  
  @todo implement a method to report a summary of all angles in a mesh.
  @todo write a method to avoid too large angles during refinement
************************************************************************  */

#ifndef __DOMAIN__
#define __DOMAIN__

class TDomain;

#include <BoundPart.h>
#include <Collection.h>
#include "Constants.h"
#include <Iterator.h>
#include <Mesh.h>
#include <ParameterDatabase.h>

class TJoint;

template <int d> class RefinementStrategy;

/** contains the boundary description, the virtual cell
    tree and macro grid */
class TDomain
{
  protected:
    
    /** @brief number of boundary parts */
    int N_BoundParts;    
    
    /** @brief boundary parts of domain */
    TBoundPart **BdParts;

    /** @brief number of all boundary components */
    int N_BoundComps;
    
    /** @brief start id of boundary component on each boundary part */
    int *StartBdCompID;

    /** @brief boundary part id's of all Interfaces */
    int *Interfaces;

    /** @brief number of holes */
    int N_Holes;
    
    /** @brief point in each hole */
    double *PointInHole;

    /** @brief number of regions */
    int N_Regions;
    
    /** @brief point in each region */
    double *PointInRegion;

    /** @brief array of all root cells of cell tree */
    TBaseCell **CellTree;
    
    /** @brief number of all root cells of cell tree */
    int N_RootCells;

    /** @brief number of virtuell cells on initial level */
    int N_InitVCells;

    /** @brief x coordinate of the start point (2D) */
    double StartX;
    /** @brief y coordinate of the start point (2D) */
    double StartY;
    /** @brief x length of bounding box */
    double BoundX;
    /** @brief y length of bounding box */
    double BoundY;

#ifdef __3D__
      /** @brief third coordinate of start point (3D) */
      double StartZ;
      /** @brief return number of cell in Y direction */
      int N_InitVCellsY;
      /** @brief z length of the bounding box */
      double BoundZ;
#endif

    /** @brief current refinment level */
    int RefLevel;
    
#ifdef  _MPI
      /** @brief array contains the global cell number of local cells (including Halo cells) */
      int *GlobalCellIndex;

      /** @brief Number of own cells (excluding Halo cells) */
      int N_OwnCells;

      /** @brief true if the domain contains periodic boundaries */
      bool PeriodicDomain = false;

      /** @brief In case of periodic joints, each vector contains the sorted
       *         global indexes of the vertices which are periodically connected
       *         It is used to complete the halo cells in MeshPartition.C. */
      std::vector<std::vector<int>> PeriodicVertexConnect;
#endif

      /** @brief see documentation of refine_and_get_hierarchy_of_collections */
      std::list<TCollection*> gridCollections;
      /**
       * A Database object which holds parameters which are of a certain
       * interest to this domain object.
       */
      ParameterDatabase db;

  public:
    /**
     * @brief Constructor.
     * @param db A database to be merged into the domain's own.
     */
    explicit TDomain(const ParameterDatabase& db);
     
    /** @brief copy constructor, deleted as a precaution */
    TDomain(const TDomain&) = delete;
     
    /** @brief Default copy assignment operator, deleted as a precaution */
    TDomain& operator=(const TDomain&) = delete;
    
    /** @brief destructor */
    ~TDomain();
    
    /**
     * Creates a database filled with default parameters. This database will
     * contain all necessary parameters for the control of a domain.
     */
    static ParameterDatabase default_domain_parameters();
    
    /**
     * Creates a database filled with default parameters. This database will
     * contain all necessary parameters for the control of a domain that should
     * by constructed as a "sandwich grid" from a 2D initial mesh.
     */
    static ParameterDatabase default_sandwich_grid_parameters();

    // Methods
    /** @brief Read in initial mesh from ".GEO" file.
     *
     * @param[in] dat Input stream which contains the initial mesh information
     * in MooNMD-native ".GEO"-format (or .xGEO).
     *
     * @param[in] isXGeoFile Set true if the input is in extended GEO (.xGEO) format.
     *
     * @return Integer O on success, other numbers otherwise.
     */
    int ReadGeo(std::istream& dat, bool isXGeoFile);


#ifdef __3D__
    /**
     *  @brief Read in initial mesh from either ".GEO" oder ".mesh" file as sandwich geometry.
     * @param[in] file_name Input stream which contains the initial mesh information
     * in either MooNMD-native ".GEO"-format or medit style ".mesh"-format.
     */
    void ReadSandwichGeo(const std::string& file_name,
                         const std::string& prm_file_name = "");

    /** @brief make boundary parameter consistent */
    void MakeBdParamsConsistent(TCollection *coll);
    
    int CloseGrid(int level);


    /** @brief Checks if the given Domain consists of  valid cells*/
    bool check() const;

#endif


    /** @brief Reads in boundary parameterization in ".PRM" file format
     *
     * @param[in] dat Input stream containing the boundary information in
     * ".PRM" style (the default MooNMD file format for domain description).
     *
     * @param[out] Flag Used in 3D only, set to 1 if a boundary component
     * of type TBdWall has to be constructed - typical for sandwich grids.
     * Otherwise set to 0.
     */
#ifdef __2D__
    void ReadBdParam(std::istream& dat);
#else
    void ReadBdParam(std::istream& dat, bool& sandwich_flag);
#endif

    /** @brief get boundary part of BdCompID */
    int GetBdPartID(int BdCompID);
    /** @brief get local number of boundary component */
    int GetLocalBdCompID(int BdCompID);
    /** @brief get local number of last boundary component on part BdPartID*/
    int GetLastLocalComp(int BdPartID)
    { return StartBdCompID[BdPartID+1] - StartBdCompID[BdPartID] - 1; }
    /** @brief set start BdCompID on boundary part i */
    void SetStartBdCompID(int BdCompID, int i)
    { StartBdCompID[i] = BdCompID; }

    /** @brief get i-th boundary part */
    TBoundPart *GetBdPart(int i)
    { return BdParts[i]; }

    /** @brief get tree of cells */
    void GetTreeInfo(TBaseCell **&celltree, int &N_rootcells)
    { 
      celltree = CellTree;
      N_rootcells = N_RootCells;
    }

    /** @brief set tree of cells */
    void SetTreeInfo(TBaseCell **celltree, int N_rootcells)
    {
      CellTree = celltree;
      N_RootCells = N_rootcells;

      #ifdef  _MPI
      N_OwnCells = 0;
      #endif 
    }

    /**
      * @brief Chooses in what way to construct the domain's geometry
      * and calls the corresponding method.
      *
      * The strings (or rather: char arrays) handed over to this function
      * determine how the domain description and the initial mesh are constructed.
      * The method itself does non of the initializing work but invokes those
      * functions which do.
      *
      * @param[in] PRM The description of the domain boundaries. Possibilities are
      *
      * "Default_UnitSquare" - Default domain boundary description (2D only).
      * "Default_UnitCube" - Default domain boundary description (3D only).
      *
      * If not this and not nullptr, PRM is interpreted as a filepath to a .PRM file
      * and gets handed over to TDomain::ReadBdParam() to be read in.
      *
      *	@param[in] GEO The description of the initial mesh. Possibilities are:
      *
      *	"InitGrid" - Grid Generator. (2D only) Call TDomain::GenInitGrid() (creates initial grid with TRIANGLE)
      * "TwoTriangles" - Default mesh. Call TDomain::TwoTriangles. (2D only)
      * "TwoTrianglesRef" - Default mesh. Call TDomain::TwoTrianglesRef. (2D only)
      * "UnitSquare" - Default mesh. ... (2D only)
      * "UnitSquareRef" - Default mesh. ... (2D only)
      * "SquareInSquare" - Default mesh. ... (2D only)
      * "SquareInSquareRef" - Default mesh. ... (2D only)
      * "PeriodicSquares" - Default mesh. ... (2D only)
      * "PeriodicSquaresLarge" - Default mesh. ... (2D only)
      * "PeriodicTrianglesLarge" - Default mesh. ... (2D only)
      * "PeriodicRectangle_2_4" - Default mesh. ... (2D only)
      *
      * "TestGrid3D" - Default mesh. ... (3D only)
      * "Default_UnitCube_Hexa" - Default regular hexahedron mesh of the unit cube. (3D only)
      * "Default_UnitCube_Tetra" - Default regular triangular mesh of the unit cube. (3D only)
      *
      * If none of these, the string is considered to be the path to a .GEO
      * file and thus is handed over to TDomain::ReadGeo().
      *
      * @note: The combination of the default domain and mesh descriptions
      * is neither checked nor tested. So use them carefully and be prepared for the worst!
      *
      */
      void Init(const std::string& PRM, const std::string& GEO);

      /**
       * @brief Initialize the domain starting from a boundary file and a mesh
       *
       * @param[in] PRM filepath to the boundary description
       * @param[in] m filepath to mesh file
       * @attention this function uses only strings (new convention, 05.2016)
       */
      void InitFromMesh(const std::string& PRM, const std::string& m);
      
#ifdef __3D__

      void GenerateFromMesh(Mesh& m);
      void buildBoundary(Mesh& m);
      void buildParMooNMesh(Mesh& m);
      void setVertices(Mesh& m);
      void allocRootCells(Mesh& m);
      void distributeJoints(Mesh& m);      

      // auxiliary vector of boundary components 
      std::vector<TBoundComp3D*> meshBoundComps;
      // auxiliary vector of vertices to store the mesh vertices
      std::vector<TVertex*> meshVertices;
      // auxiliary vector of Joints
      std::vector<TJoint*> meshJoints;
#endif
    /** @brief write mesh into a postscript file */
    void PS(const char *name, Iterators iterator, int arg);
    /** @brief write collection into a postscript file */
    void PS(const char *name, TCollection *Coll);
    
    /** @brief refine the grid according the cell refinement descriptors */
    int Refine();
    /** @brief refine all cells regular */
    int RegRefineAll();
    /** @brief generate a 1-regular grid */
    int Gen1RegGrid();
    /** @brief refine the finest grid according the given indicator function */
    void RefineByIndicator(
#ifdef __2D__ 
  DoubleFunct2D *Indicator
#else
  DoubleFunct3D *Indicator
#endif
  );
    /** @brief refine the finest grid according to the marked cells in 
     *  refinement strategy.
     * @param[in] ConfClosure Do or don't do conforming closures 
     */
    template <int d>
    void RefineByRefinementStrategy(RefinementStrategy<d>& strategy);
    /** @brief refine the finest grid if necessary in order to get a grid with
     *  conforming closures */
    int MakeConfClosure();

    
    /** @brief refine/derefine algorithm for a 1-regular grid, geolevel of all
        cells on the finest grid is between MinLevel and MaxLevel */
    void Refine1Reg(int MinLevel, int MaxLevel);

    /** @brief derefinemnt */
    void DeRefine();

    /** @brief convert all finest quadrangles into two triangles */
    int ConvertQuadToTri(int type);
    
    /// @brief do a barycentric refinement for triangles/tetrahedra
    /// @note this throws an exception whenever a quadrangle/hexahedron is found
    void barycentric_refinement();
    
    int get_ref_level() const
    { return RefLevel; }

    /** @brief produce a collection with all cells returned by iterator it */
    TCollection *GetCollection(Iterators it, int level) const;
    
    /** @brief produce a collection with all cells of a given Collection which
     *         have the given reference as ReferenceID 
     * 
     * This will give you a subcollection.
     */
    TCollection *GetCollection(TCollection *coll, int reference);
    
    /**
     * @brief get collection of cells on a certain level and with certain 
     *        reference id
     */
    TCollection *GetCollection(Iterators it, int level, int ID) const;

    /** @brief get bounding box parameters */
    void GetBoundBox(double &startx, double &starty,
                     double &boundx, double &boundy)
    {
      startx = StartX;
      starty = StartY;
      boundx = BoundX;
      boundy = BoundY;
    }
    
#ifdef __3D__
    /** @brief get bounding box parameters */
    void GetBoundBox(double &startx, double &starty, double &startz,
                     double &boundx, double &boundy, double &boundz)
    {
      startx = StartX;
      starty = StartY;
      startz = StartZ;
      boundx = BoundX;
      boundy = BoundY;
      boundz = BoundZ;
    }
    
    void SetBoundBox(double startx, double starty, double startz,
                     double boundx, double boundy, double boundz)
    {
      StartX = startx;
      StartY = starty;
      StartZ = startz;
      BoundX = boundx;
      BoundY = boundy;
      BoundZ = boundz;
    }
#endif
    
    // test
    #ifndef __3D__
      void TestGrid1();
      void TestGrid2();
      void TestGrid3();
      void TestShishkin();
      void TriangleShishkin();
      void UnitSquare();
      void DrivenCavitySquareQuads();
      void UnitSquareRef();
      void TwoTriangles(bool bottomleft_to_topright = true);
      void TwoTrianglesRef();
      void SquareInSquare();
      void SquareInSquareRef();
      void SetBoundBox(double boundx, double boundy);
      void SetBoundBoxstart(double startx , double starty);
      void RefCardioide(double A);
      void PeriodicSquares();
      void PeriodicTriangles();
      void PeriodicSquaresLarge();
      void PeriodicRectangle_2_4();
      void PeriodicTrianglesLarge();
      void QuadShishkin(double tau1, double tau2);
      void Rectangular(int dimx, int dimy);
      void TestTriaConf();
      void TestTriaConf2();
      void UnitSquare_US22();
      void CheckCells();

      /**
       * @brief Initialize the domain boundary as if data/UnitSquare.PRM
       * would have been read in as .PRM file.
       *
       * This is useful for testing purposes, when the program must be
       * independent of the path from working directory to "data" directory.
       */
      void initializeDefaultUnitSquareBdry();

      /** Boundary for a square stretching between (-1,-1) and (1,1),
       *  as needed for the driven cavity example.*/
      void initializeDrivenCavitySquareBdry();

    #else
      void TestGrid3D();
      /**
       * @brief Initialize the domain boundary as if data/Wuerfel.PRM
       * would have been read in as .PRM file.
       *
       * This method is useful for tests, which are supposed to be independent
       * of whether an extern .PRM-file is available. To let the domain call this
       * method when initializing, pass "Default_UnitCube" as first argument
       * to the init method.
       */
      int initializeDefaultCubeBdry();
      
      /**
      * @brief Initialize the boundary of the domain
      * in the form of the block with the edegs of the lengths 1, 2, 3.
      *
      * To let the domain call this method when initializing,
      * pass "Default_Block_1x2x3" as the parameter boundary_file.
      */
      int initialize_block_1x2x3_bdry();

      /**
       * @brief Initialize the initial mesh as if data/Wuerfel.GEO
       * would have been read in as .GEO file.
       *
       * This method is useful for tests, which are supposed to be independent
       * of whether an extern .PRM-file is available. To let the domain call this
       * method when initializing, pass "Default_UnitCube_Hexa" as second argument
       * to the init method.
       */
      void initialize_cube_hexa_mesh();

      /**
       * @brief Initialize the initial mesh as if data/SixTetras.GEO
       * would have been read in as .GEO file.
       *
       * This method is useful for tests, which are supposed to be independent
       * of whether an extern .PRM-file is available. To let the domain call this
       * method when initializing, pass "Default_UnitCube_Tetra" as second argument
       * to the init method.
       */
      void initialize_cube_tetra_mesh();
      
      /**
      * @brief Initialize the initial mesh of tetrahedra for the domain
      * in the form of the block with the edges of the lengths 1, 2, 3.
      *
      * To let the domain call this method when initializing, pass
      * "Default_Block_1x2x3_Tetra" as the parameter geo_file.
      */
      void initialize_block_1x2x3_tetra_mesh();
      
      /**
      * @brief Initialize the initial mesh of hexahedra for the domain
      * in the form of the block with the edges of the lengths 1, 2, 3.
      *
      * To let the domain call this method when initializing, pass
      * "Default_Block_1x2x3_Hexa" as the parameter geo_file.
      */
      void initialize_block_1x2x3_hexa_mesh();

      void SetBoundBox(double boundx, double boundy, double boundz);
    #endif

    #ifdef __3D__
      int Grape(const char *name, TCollection *coll);

    #endif

    void TetrameshGen();


#ifdef  _MPI
  void ReplaceTreeInfo(int n_cells, TBaseCell **cells, int *GLOB_cellIndex, int n_OwnCells)
    {
    if(CellTree) delete[] CellTree;
    N_RootCells = n_cells;
    CellTree = cells;
    GlobalCellIndex = GLOB_cellIndex;
    N_OwnCells = n_OwnCells;
    }

  void SetN_OwnCells(int n_OwnCells)
    { N_OwnCells = n_OwnCells; }

  int GetN_OwnCells()
    { return N_OwnCells; }

  int GetN_HaloCells()
    { return (N_RootCells - N_OwnCells); }

    void SetPeriodicDomain()
    { PeriodicDomain = true; }

    bool is_periodic()
    { return PeriodicDomain; }

    std::vector<std::vector<int>>* get_PeriodicVertexConnect_ptr()
    { return &PeriodicVertexConnect; }

    const std::vector<std::vector<int>> get_PeriodicVertexConnect()
    { return PeriodicVertexConnect; }
#endif

#ifdef __3D__
  public: 

     /** @brief generate edge info in 3D mesh **/
     // TODO CB Add documentation. This method does mark certain Vertices as boundary vertices
     // (TVertex::SetAsBoundVertex(...) ) and creates edge objects, pointers to which are then
     // stored as members of certain mesh cells (TBaseCell::SetEdge(...)).
     // Does work only on currently finest level.
     int GenerateEdgeInfo();

#endif

  /** @brief adaptive refine  */
  int AdaptRefineAll();   

  /**
   * @return The value of parameter "refinement_n_initial_steps"
   * stored in this object's database.
   */
  size_t get_n_initial_refinement_steps() const;

  /**
   * @return The value of parameter "refinement_max_n_adaptive_steps"
   * stored in this object's database.
   */
  size_t get_max_n_adaptive_steps() const;


  /**
   * Prints info on this domain,
   * to console and outfile.
   * @param name A name for the domain.
   */
  void print_info(const std::string& name) const;

  /** This is a method which wraps together two things which are awful about
   * MPI ParMooN, especially in connection with multigrid.
   * It must be used in main programs whenever one watns to set up a problem which
   * is already adapted to use MPI multigrid (CD3D, NSE3D,...) and therefore gets
   * a list of TCollections and not a Domain as parameter.
   * Read documentation of
   *   determine_n_refinement_steps_multigrid
   * for a description of the problem.
   * @param[in] parmoon_db The input database.
   * @return A hierarchy of geometric grids which can be used for multigrid,
   * finest grid first.
   *
   */
  std::list<TCollection*> refine_and_get_hierarchy_of_collections(
                           const ParameterDatabase& parmoon_db,
#ifdef __2D__
                           BoundCondFunct2D* BoundaryCondition = nullptr);
#else
                           BoundCondFunct3D* BoundaryCondition = nullptr);
#endif


  const std::list<TCollection*> get_grid_collections() const;

  /// Get a const reference to the database of this domain object.
  const ParameterDatabase& get_database() const {return this->db;};

  //TODO reorganize methods, put everything private here!
  private:
    #ifdef __2D__
      /** @brief make initial 2D grid */
      int MakeGrid(double *DCORVG, int *KVERT, int *KNPR, int N_Vertices,
                   int NVE);
      /** @brief make initial 2D grid, extended version which sets ReferenceID
       *         in cells */
      int MakeGrid(double *DCORVG, int *KVERT, int *KNPR, int *ELEMSREF,
                   int N_Vertices, int NVE);
    #else
      /** @brief make initial 3D grid */
      int MakeGrid(double *DCORVG, int *KVERT, int *KNPR, int *ELEMSREF,
                   int N_Vertices, int NVE, int *BoundFaces, int *FaceParam,
                   int NBF, int NVpF,
                   int *Interfaceparam, int N_Interfaces);

      /** @brief make initial sandwich grid */
      void MakeSandwichGrid(
    		  double *DCORVG, int *KVERT, int *KNPR,
    		  int N_Vertices, int NVE,
    		  double DriftX, double DriftY, double DriftZ,
                  double ConicScale,
    		  const std::vector<double>& Lambda);

     #endif

};

/**
 * @todo Dear team geometry: please read this description and maybe find a way to
 * reimplement the domain decompositioning in a way which does not destroy the
 * cell tree hierarchy (and, btw. write tests for the domain decomp).
 *
 * This is a helper methods which must be called in all main programs which
 * intend to make use of a parallelized multigrid solver, i.e. the (T)CD3D and
 * (T)NSE3D mains programs and their test programs.
 *
 * The point is, that due to the modifications which apppear to the domain in
 * Partition_Mesh3D and Domain_Crop, it is a bit tricky to gain a grid hierarchy
 * for multigrid. Actually one can pick a grid (a "TCollection") only, when it
 * currently is the finest one of the domain by using
 *   domain.GetCollection(It_Finest, 0)
 * Other "cell tree iterators" cannot be used succesfully, mainly because
 * Partition_Mesh3D and/or Domain_Crop destroy the cell tree and create a new one,
 * which only consists of one level.
 * If now, as we usually do in ParMoooN, one has an initiali grid given and
 * wants to perform some initial refinement steps to gain the fine computational
 * grid, one has to split these initial steps, into some performed BEFORE partitioning
 * the domain and some AFTER partitioning the Domain. The number of steps AFTER
 * domain partitioning must be such, that one can pick as many grids as one needs
 * for the requested multigrid hierarchy one by one using
 *   domain.GetCollection(It_Finest, 0).
 *
 * @param[in] multigrid_type May have the values "standard" and "mdml". Otherwise
 * the program quits.
 * @param[in] n_multigrid_levels The number of GEOMETRIC multigrid levels.
 * The algebraic hierarchy of mdml will have one extra level, but n_multigrid_levels
 * will be the same as for mdml. Interpreting the parameter in such a way, means,
 * that mdml and standard mg applied to the same problem with the same parameter
 * n_multigrid_levels wil lead to both of them solving the same coarse grid problem,
 * which we find preferable.
 * @param[in] n_initial_refinement_steps The number of refinement steps to be
 * performed in total
 * @param[out] n_ref_before Number of refinement steps to be performed before
 * partitioning the domain and performing
 */
void determine_n_refinement_steps_multigrid(
  const std::string& multigrid_type,
  int n_multigrid_levels,
  int n_initial_refinement_steps,
  int& n_ref_before, int& n_ref_after);

#endif


