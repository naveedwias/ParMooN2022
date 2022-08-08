#ifndef __SYSTEM_TIMECONVECTIONDIFFUSION_H__
#define __SYSTEM_TIMECONVECTIONDIFFUSION_H__

#include <BlockFEMatrix.h>
#include <BlockVector.h>
#ifdef __2D__
#include "Example_TimeCD2D.h"
#include "FEFunction2D.h"
#else
#include "Example_TimeCD3D.h"
#include "FEFunction3D.h"
#endif
#include "Solver.h"
#include "DataWriter.h"
#include "LocalAssembling.h"
#include "TimeDiscretizations.h"
#include "templateNames.h"
#include "CheckpointIO.h"

#include <vector>
#include <array>
#include <deque>

#include "AlgebraicFluxCorrectionTCD.h"

class TDomain;

template<int d>
class TimeConvectionDiffusion
{
  public:
    using FEFunction = typename Template_names<d>::FEFunction;
    using MultiIndex_vector = typename Template_names<d>::MultiIndex_vector;
    using FESpace = typename Template_names<d>::FESpace;
    using Example_TimeCD = typename Template_names<d>::Example_TimeCD;
  
    /** @brief constructor
     * This constructor calls the other constructor creating an Example_TimeCD
     * object. 
     */
    TimeConvectionDiffusion(const TDomain& domain,
                            const ParameterDatabase& param_db);
    
    /** @brief The standard constructor, can be used for multigrid and 
     * non-multigrid.
     *
     * @param[in] domain The computational domain providing the grids
     * @param[in] param_db A parameter database with parameters concerning this
     *                     class or any of its members (fe space, solver,
     *                     assemble,...)
     * @param[in] example The example which is to be calculated.
     */
    TimeConvectionDiffusion(const TDomain& domain,
                            const ParameterDatabase &param_db,
                            const Example_TimeCD& example);
    
    
    ~TimeConvectionDiffusion()
    {
      if(afcTcdPtr != nullptr)
        destroy_afcTcdObject();
    }
    
    
    /** @brief return a database with all parameters necessary for 
     * time-dependent convection-diffusion (tcd) probems
     * 
     * If complete is false, this will only return those parameters which are
     * used by this class directly. If the parameters which are used by of all
     * its members are also desired, set complete to true.
     */
    static ParameterDatabase default_tcd_database(bool complete = false);
    
    /** @brief Assemble all the matrices before the time iterations
     * 
     * This includes the assembling of: Stiff_matrix, Mass_Matrix, 
     * (additional matrixK in case of SUPG stabilization), rhs
     */
    void assemble_initial_time(FEFunction* velo1 = nullptr,
                               FEFunction* velo2 = nullptr,
                               FEFunction* velo3 = nullptr);
    
    /** @brief assemble the matrices
     * this function will assemble the stiffness matrix and rhs
     * In addition the system matrix and the rhs which passes to the solver 
     * are also prepared within the function
     */
    void assemble(FEFunction* velo1 = nullptr,
                  FEFunction* velo2 = nullptr,
                  FEFunction* velo3 = nullptr,
                  FEFunction* sources_and_sinks = nullptr);
    
    /** @brief solve the system
     */
    void solve();
    
    void fem_fct_nonlinear_loop();
    
    /** @brief Calculates the new solution using algebraic flux correction
     * via the Runge-Kutta-Heun method. */
     void calculateNewSolutionExplicitly();
    
    /** @brief measure errors and write solution
     */
    void output();
    
    /** @brief Prints solution values at given point into file
     * (using "append" option!).
     * If no file name is provided, it prints them to standard output. */
    void printValuesAtPoint (double x, double y, double z,
                              const std::string& fileName = "");
    
     // getters and setters
    const Example_TimeCD& get_example() const
    { return example; }
    const FEFunction & get_function() const
    { return this->systems.front().fe_function; }
    FEFunction & get_function()
    { return this->systems.front().fe_function; }
    const BlockFEMatrix & get_stiff_matrix() const
    { return this->systems.front().stiffness_matrix; }
    const BlockVector & get_rhs() const
    { return this->systems.front().rhs; }
    BlockVector & get_rhs()
    { return this->systems.front().rhs; }
    const BlockVector & get_solution() const
    { return this->systems.front().solution; }
    const BlockVector & get_solution_m1() const
    { return this->systems.front().solution_m1; }
    const BlockVector & get_solution_m2() const
    { return this->systems.front().solution_m2; }
    std::shared_ptr<const FESpace> get_space() const
    { return this->systems.front().fe_space; }
    const ParameterDatabase & get_db() const
    { return db; }
    const Solver<BlockFEMatrix, BlockVector> & get_solver() const
    { return solver; }
    
    TimeDiscretization& get_time_stepping_scheme()
    { return time_stepping_scheme; }
    const TimeDiscretization& get_time_stepping_scheme() const
    { return time_stepping_scheme; }
    
    /**
     * @brief return the computed errors at each discret time point
     */
    std::array<double, 3> get_errors() const;
    
   protected:
    /** @brief store a complete system on a particular grid.
     * 
     * This combines a matrix, rhs, solution, spaces and functions 
     * needed to describe a Time CD problem
     */
    struct SystemPerGrid
    {
      /** @brief Finite element space */
      std::shared_ptr<const FESpace> fe_space;
      /** @brief Stiffness Matrix */
      BlockFEMatrix stiffness_matrix;
      /** @brief Mass matrix */
      BlockFEMatrix mass_matrix;
      /** @brief right hand side vector */
      BlockVector rhs;
      /** @brief solution vector */
      BlockVector solution;

      /** @brief Finite element function */
      FEFunction fe_function;
      
      /** @brief
       * old solution for the computation of the residual
       * that passes as an FEFunction to the local assembling
       * routines
       * 
       * @todo why don't we store a vector of size 
       * time_stepping_scheme.n_old_solutions() here?
       */
      BlockVector solution_m1;
      FEFunction u_m1;
      BlockVector solution_m2;
      FEFunction u_m2;

      /** @brief constructor in mpi case
       * @param[in] example The current example.
       * @param[in] coll The collection of mesh cells of this grid.
       * @param[in] ansatz_order polynomial order for finite element space
       */
      SystemPerGrid(const Example_TimeCD& example, const TCollection& coll,
                    int ansatz_order);
      /**
       * Special member functions mostly deleted, for struct takes ownership of
       * the bad classes TFEFunction2D/TFEFunction3D and TFESpace2D/TFESpace3D.
       */
      //! Delete copy constructor.
      SystemPerGrid(const SystemPerGrid&) = delete;

      //! Delete move constructor.
      SystemPerGrid(SystemPerGrid&&) = delete;

      //! Delete copy assignment operator.
      SystemPerGrid& operator=(const SystemPerGrid&) = delete;

      //! Delete move assignment operator.
      SystemPerGrid& operator=(SystemPerGrid&&) = delete;

      //! Default destructor. Most likely causes memory leaks.
      ~SystemPerGrid() = default;

    };
    
    /** @brief a local parameter database which controls this class
     *
     * The database given to the constructor will be merged into this one. Only
     * parameters which are of interest to this class are stored (and the
     * default ParMooN parameters). Note that this usually does not include
     * other parameters such as solver parameters. Those are only in the
     * Solver object.
     */
    ParameterDatabase db;
        
    /** @brief a solver object which will solve the linear system
     * 
     * Storing it means that for a direct solver we also store the factorization
     * which is usually not necessary.
     */
    Solver<BlockFEMatrix, BlockVector> solver;
    
    /** @brief a complete system on each grid
     *
     * Note that the size of this deque is at least one and larger than that
     * only in case of multigrid (when it holds as many systems as there are
     * multigrid levels).
     */
    std::deque<SystemPerGrid> systems;
    
    
    /** @brief Definition of the used example */
    const Example_TimeCD example;
    
    /** @brief old right hand side vectior 
     * this will be used to save the right hand side from the 
     * previous time step that will be used for different 
     * time stepping schemes
     */
    BlockVector old_rhs;

    // extended to 8 for supg error. Therefore the last two values are only used
    // in case of SUPG space discretization
    static constexpr int n_errors = 8;

    /** @brief store the errors to compute accumulated error norms
     *
     * The errors are stored in the following order:
     *  - L2(0,T;L2)
     *  - previous L2 error squared
     *  - L2(0,T;H1)
     *  - previous H1 error squared
     *  - L_inf(0,T;L_inf)
     *  - L_inf(0,T;L2)
     *  - L2(0,T;SUPG)
     *  - previous SUPG error squared
     */
    std::array<double, n_errors> errors;
    
    // time at which L_inf(0,T;L_inf) and L_inf(0,T;L2) are reached
    double t_Linf = 0.;
    double t_L2 = 0.;

    /** @brief output object */
    DataWriter<d> outputWriter;
    
    /// @brief time stepping scheme object to access everything
    TimeDiscretization time_stepping_scheme;
    /// @brief system right hand side which passes to the solver
    BlockVector rhs_from_time_disc;
    
    /// @brief to read/write the current solution as a checkpoint
    CheckpointIO checkpoint_io;
    
    /**
     * @brief Check whether the program will be working with the
     * current input parameters.
     *
     * ParMooN is work in progress, and so is this class. This method checks the
     * parameters stored in the database and stops execution of the program, if
     * some of these do not match. The method is a little makeshift and the
     * cases caught here are various, but basically it is intended to stop
     * execution of cases with parameter combinations which are not implemented
     * for Time_ConvectionDiffusion or are currently known to be problematic.
     *
     * This is not yet a guarantee for a functioning program, but is
     * intended to be, someday. Eventually this method and the like
     * will be moved to TDatabase.
     */
    void check_and_set_parameters();
    
    /** @brief write some information (number of cells, dofs, ...) */
    void output_problem_size_info() const;

    /**
     * Scramble together the parameters which Assemble2D/Assemble3D needs and
     * call it. Is only put here to keep the code of assemble() slender.
     * assemble() should take care of the right choice of the LocalAssembling
     * object and whether it fits the system matrix' block structure
     * @param s The sytem where rhs, mass and stiffness matrices are to be assembled.
     * @param la_stiff The local assembling object of choice.
     * @param assemble_both If true, both stiffness (+rhs) and mass matrix are
     * assembled, if false only stiffness matrix and rhs.
     */
    void call_assembling_routine(TimeConvectionDiffusion::SystemPerGrid& system,
                                 LocalAssembling<d>& la_stiff,
                                 bool assemble_both);
    void modify_and_call_assembling_routine(
        SystemPerGrid& s,
        LocalAssembling<d>& la,
        bool assemble_both,
        FEFunction* velo1 = nullptr,
        FEFunction* velo2 = nullptr,
        FEFunction* velo3 = nullptr,
        FEFunction* sources_and_sinks = nullptr);
    
    
     /** @brief Object that carries out algebraic flux correction for TCD. */
     AFC_TCD<d>* afcTcdPtr = nullptr;
     
     /** @brief Supplementary object for algebraic flux correction for TCD.
      Carries lists of active and nonactive master rows of MPI processes.*/
     RowsColsRecord<d>* rowsColsRecordPtr = nullptr;
     
     /** @brief Supplementary object for algebraic flux correction for TCD.
      Sets boundary values for BlockVectors.*/
     TCD_BC_setter<d>* bcSetterPtr = nullptr;
     
     void create_afcTcdObject(AFC_TCD_params<d>& passedParameters);
     
     void destroy_afcTcdObject();
    
     void setAfcTcdObjectParameters(AFC_TCD_params<d>& params);

    /**
     * @brief Compute the defects A.x0 - b and A.xf - b and print their norm
     *
     * A is the current matrix, b is the current right hand side,
     * x0 is the last solution (initial value for the current step) and xf is
     * the current solution.
     * Call this function after solving the system, but before resetting the
     * matrices or updating the solution_m1.
     */
    void compute_residuals();
};

#endif // __SYSTEM_TIMECONVECTIONDIFFUSION_H__
