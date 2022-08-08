#ifndef __SYSTEM_NAVIERSTOKES_DG__
#define __SYSTEM_NAVIERSTOKES_DG__

#include <BlockFEMatrix.h>
#ifdef __2D__
#include <Example_NSE2D.h>
#include <FEFunction2D.h>
#else
#include <Example_NSE3D.h>
#include <FEFunction3D.h>
#endif
#include <BlockVector.h>
#include <ParameterDatabase.h>
#include <Solver.h>
#include <DataWriter.h>
#include "CheckpointIO.h"

#include <deque>
#include <array>

template <int d>
class NavierStokes_DG
{
  public:
    using FEFunction = typename Template_names<d>::FEFunction;
    using FESpace = typename Template_names<d>::FESpace;
    using Example_NSE = typename Template_names<d>::Example_NSE;

  protected:

    /** @brief store a complete system on a particular grid
     *
     * This combines a matrix, rhs, solution, spaces and functions needed to
     * describe one Darcy problem in 3D.
     */
    struct System_per_grid
    {
      /** @brief Finite Element space for the velocity */
      std::shared_ptr<const FESpace> velocity_space;
      /** @brief Finite Element space for the pressure */
      std::shared_ptr<const FESpace> pressure_space;
      /** @brief the system matrix (here one block)
       *  [ A  BT ]
       *  [ B  C  ]
       */
      BlockFEMatrix matrix;
      /** @brief the right hand side vector */
      BlockVector rhs;
      /** @brief solution vector with two components. */
      BlockVector solution;
      /** @brief Finite Element function for velocity */
      FEFunction u;
      /** @brief Finite Element function for pressure */
      FEFunction p;

      /** @brief constructor */
      System_per_grid(const Example_NSE& example, const TCollection& coll,
                      std::pair<int,int> velocity_pressure_orders);

      //! Delete copy constructor.
      System_per_grid(const System_per_grid&) = delete;

      //! Delete move constructor.
      System_per_grid(System_per_grid&&) = delete;

      //! Delete copy assignment operator.
      System_per_grid& operator=(const System_per_grid&) = delete;

      //! Delete move assignment operator.
      System_per_grid& operator=(System_per_grid&&) = delete;

      //! Default destructor.
      ~System_per_grid() = default;
    };

    /** @brief a complete system on each grid
     *
     * Note that the size of this deque is at least one and larger only in case
     * of multigrid.
     */
    std::deque<System_per_grid> systems;

    /** @brief Definition of the used example */
    const Example_NSE example;

    /** @brief a local parameter database which constrols this class
     *
     * The database given to the constructor will be merged into this one. Only
     * parameters which are of interest to this class are stored (and the
     * defualt ParMooN parameters). Note that this usually does not include
     * other parameters such as solver parameters. Those are only in the
     * Darcy3D::solver object.
     */
    ParameterDatabase db;

    /** @brief class for output handling */
    DataWriter<d> outputWriter;

    /// @brief to read/write the solution as a checkpoint. Of course, reading an
    /// initial solution makes little sense when using a direct solver for a
    /// linear problem.
    CheckpointIO checkpoint_io;

    /** @brief a solver object which will solve the linear system
     *
     * Storing it means that for a direct solver we also store the factorization
     * which is usually not necessary.
     */
    Solver<BlockFEMatrix, BlockVector> solver;

    /** @brief store the errors to access them from outside this class
     *
     * This array is filled during a call to Darcy3D::output if
     * 'output_compute_errors' is set to true. The exact solution is
     * taken from Darcy3D::example. If that example does not provide an exact
     * solution, typically it is set to be zero, so that this array contains
     * the norms of the solution instead of the error.
     *
     * The errors are stored in the following order:
     *
     *  - L2 error of velocity
     *  - L2 error of divergence of velocity
     *  - H1-semi error of velocity
     *  - L2 error of pressure
     *  - H1-semi error of pressure
     */
    std::array<double, 5> errors;

    /** @brief set parameters in database
     *
     * This functions checks if the parameters in the database are meaningful
     * and resets them otherwise. The hope is that after calling this function
     * this class is fully functional.
     *
     * If some parameters are set to unsupported values, an error occurs and
     * throws an exception.
     */
    void set_parameters();

  public:

    /** @brief constructor
     *
     * This constructor calls the other constructor creating an Example_NSE
     * object for you. See there for more documentation.
     */
    NavierStokes_DG(const TDomain& domain, const ParameterDatabase & db);

    /** @brief constructor
     *
     * The domain must have been refined a couple of times already if you want
     * to use multigrid. On the finest level the finite element spaces and
     * functions as well as matrices, solution and right hand side vectors are
     * initialized. This is determined by the values in the parameter database.
     *
     * The reference_id can be used if only the cells with the give reference_id
     * should be used. The default implies all cells.
     */
    NavierStokes_DG(const TDomain& domain, const ParameterDatabase & db,
          const Example_NSE ex);

    //! Delete copy constructor.
    NavierStokes_DG(const NavierStokes_DG&) = delete;

    //! Delete move constructor.
    NavierStokes_DG(NavierStokes_DG&&) = delete;

    //! Delete copy assignment operator.
    NavierStokes_DG& operator=(const NavierStokes_DG&) = delete;

    //! Delete move assignment operator.
    NavierStokes_DG& operator=(NavierStokes_DG&&) = delete;

    /** @brief standard destructor */
    ~NavierStokes_DG() = default;

    /** @brief return a database with all parameters necessary for Darcy.
     *
     * If complete is false, this will only return those parameters which are
     * used by this class directly. If the parameters which are used by of all
     * its members are also desired, set complete to true.
     */
    static ParameterDatabase default_NavierStokes_DG_database(bool complete = false);

    /** @brief assemble matrix,
     *
     * depending on 'this->db["space_discretization_type]' different (local)
     * assembling routines are used. Also in case of multigrid the matrices
     * on all grids are assembled.
     */
    void assemble();
    
    //added 21.02.22
    void assemble_nonlinear_term();
    bool stop_it(unsigned int iteration_counter);

    /** @brief solve the system */
    void solve();

    /**
     * @brief measure errors and write pictures
     *
     * The current errors will be printed out if
     * TDatabase::ParamDB->MEASURE_ERRORS is true. The errors will be stored
     * and are then accessible via the get*Error methods. If desired, further
     * output, e.g., vtk files are created.
     *
     * @param i suffix for output file name, -1 means no suffix
     */
    void output(int i = -1);

    /// @name return computed errors
    ///
    /// You have to call Darcy3D::output for any of these to return a
    /// meaningful value.
    //@{
    /// @brief return the computed L2 error of the velocity
    double getL2VelocityError() const;
    /// @brief return the computed L2 error of the divergence of the velocity
    double getL2DivergenceError() const;
    /// @brief return the computed H1-semi error of the velocity
    double getH1SemiVelocityError() const;
    /// @brief return the computed L2 error of the pressure
    double getL2PressureError() const;
    /// @brief return the computed L2 error of the pressure
    double getH1SemiPressureError() const;
    //@}

    // getters and setters
    const BlockFEMatrix & get_matrix() const
    { return this->systems.front().matrix; }
    BlockFEMatrix & get_matrix()
    { return this->systems.front().matrix; }
    const BlockVector & get_rhs() const
    { return this->systems.front().rhs; }
    BlockVector & get_rhs()
    { return this->systems.front().rhs; }
    const FEFunction & get_velocity() const
    { return this->systems.front().u; }
    const FEFunction & get_pressure() const
    { return this->systems.front().p; }
    std::shared_ptr<const FESpace> get_velocity_space() const
    { return this->systems.front().velocity_space; }
    std::shared_ptr<const FESpace> get_pressure_space() const
    { return this->systems.front().pressure_space; }
    const BlockVector & get_solution() const;
//     { return this->systems.front().solution; }
    BlockVector & get_solution();
//     { return this->systems.front().solution; }
    unsigned int get_size() const
    { return this->systems.front().solution.length(); }
    const Example_NSE& get_example() const
    { return example; }
    const ParameterDatabase & get_db() const
    { return db; }
};

#endif // __SYSTEM_NAVIERSTOKES_DG__
