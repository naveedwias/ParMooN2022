#ifndef __SYSTEM_TIMENAVIERSTOKES_H__
#define __SYSTEM_TIMENAVIERSTOKES_H__

#include "BlockFEMatrix.h"
#ifdef __2D__
#include "Example_TimeNSE2D.h"
#include "FEFunction2D.h"
#include "FEVectFunct2D.h"
#else
#include "Example_TimeNSE3D.h"
#include "FEFunction3D.h"
#include "FEVectFunct3D.h"
#endif
#include "BlockVector.h"
#include "LinesEval.h"
#include "ParameterDatabase.h"
#include "Solver.h"
#include "templateNames.h"
#include "LocalAssembling.h"
#include "DataWriter.h"
#include "Residuals.h"
#include "MainUtilities.h" // FixedSizeQueue
#include "TimeDiscretizations.h"
#include "CheckpointIO.h"
#include "WindkesselBoundary.h"

class PointwiseAssemblyData;

#include <deque>
#include <array>

template <int d>
class TimeNavierStokes
{
  public:
    using FEFunction = typename Template_names<d>::FEFunction;
    using FEVectFunct = typename Template_names<d>::FEVectFunct;
    using FESpace = typename Template_names<d>::FESpace;
    using Example_TimeNSE = typename Template_names<d>::Example_TimeNSE;

    /** @brief see the other constructor */
    TimeNavierStokes(const TDomain& domain, const ParameterDatabase& param_db);

    /** @brief constructor
     *
     * The domain must have been refined a couple of times already. On the
     * finest level the finite element spaces and functions as well as
     * matrices, solution and right hand side vectors are initialized.
     *
     * @param domain the computational domain to get the grid(s)
     * @param param_db parameters controlling this class
     * @param example The example to use
     */
    TimeNavierStokes(const TDomain& domain, const ParameterDatabase& param_db,
                     const Example_TimeNSE& ex);

    virtual ~TimeNavierStokes() {}

    /** @brief return a database with all parameters necessary for 
     * time-dependent Navier-Stokes (tnse) probems
     *
     * If complete is false, this will only return those parameters which are
     * used by this class directly. If the parameters which are used by of all
     * its members are also desired, set complete to true.
     */
    static ParameterDatabase default_tnse_database(bool complete = false);

    /** @brief set the solution to its initial values.
     *
     * This reads a solution from file if the parameter 'read_initial_solution'
     * is true. Otherwise the initial solution as specified in the example is
     * interpolated. If no solution is known, this is zero.
     *
     * @return true if a solution was read from file, false if an interpolation
     * is used.
     *
     * @note This method is called during construction. In case one wants to 
     * solve multiple time-dependent Navier-Stokes problems, this can be used
     * to initialize another time loop with the same object.
     */
    bool initialize_solution();

    /** @brief Assemble all the matrices and rhs before the time iterations
     *
     * This includes the assembling of: Stiff_matrix, Mass_Matrix,
     * (additional matrixK in case of SUPG stabilization), rhs
     */
    void assemble_initial_time();

    /** @brief This will assemble the right-hand side, will prepare the 
     * right-hand side using the class TimeDiscretization object. 
     * The pressure blocks are scaled properly and the nonlinear matrices
     * are assembled again. The system matrix is assemble using the 
     * object of class TimeDiscretization and modify the matrices if needed.
     */
    void assemble_matrices_rhs(unsigned int it_counter);

    /** @brief solve the system */
    void solve();

    /**
     * @brief Compute the defect Ax-b and store its norm.
     *
     * A is the current matrix, x is the current solution and b is the
     * right hand side. Call this function after assembling the nonlinear
     * matrix with the current solution.
     */
    void compute_residuals();

    Residuals compute_residuals_unstored();

    /** @brief check if one of the stopping criteria is fulfilled
     *
     * either converged, maximum number of iterations reached, or slow
     * convergence
     * 
     * This version does not reset on convergence failure.
     *
     * @param it_counter current iterate
     */
    bool stop_it(unsigned int it_counter);

    /** @brief check if one of the stopping criteria is fulfilled
     *
     * either converged, maximum number of iterations reached, or slow
     * convergence
     *
     * This version resets on convergence failure.
     *
     * @param it_counter current iterate
     * @param iteration_failed receives failure flag
     */
    bool stop_it(unsigned int it_counter, bool& iteration_failed);

    /** @brief check if one of the stopping criteria is fulfilled
     *
     * either converged, maximum number of iterations reached, or slow
     * convergence
     *
     * @param it_counter current iterate
     * @param reset_on_failure Reset on convergence failure?
     * @param iteration_failed receives failure flag
     */
    bool stop_it(unsigned int it_counter, bool reset_on_failure, bool& iteration_failed);

    /** @brief
     * compute errors and write solution
     */
    void output();

    /** @brief check if the nonlinear term is replaced by a linear one.
     * 
     * This is the case whenever the imex or the fully explicit scheme is used,
     * set by the parameter 'time_discretization_nonlinear_term'.
     */
    bool linearized();

    void get_time_derivative(BlockVector& dt);

    // getters and setters
    const BlockFEMatrix & get_matrix() const
    { return this->systems.front().matrix; }
    BlockFEMatrix & get_matrix()
    { return this->systems.front().matrix; }

    const BlockVector & get_rhs() const
    { return this->systems.front().rhs; }
    BlockVector & get_rhs()
    { return this->systems.front().rhs; }

    const FEVectFunct & get_velocity() const
    { return this->systems.front().u; }

    const FEVectFunct & get_velocity_old() const
    { return this->systems.front().u_m1;}
    // try not to use this as it is not const
    std::unique_ptr<FEFunction> get_velocity_component(int i);

    const FEFunction & get_pressure() const
    { return this->systems.front().p; }

    FEFunction & get_pressure()
    { return this->systems.front().p; }

    std::shared_ptr<const FESpace> get_velocity_space() const
    { return this->systems.front().velocity_space; }
    std::shared_ptr<const FESpace> get_pressure_space() const
    { return this->systems.front().pressure_space; }
    const BlockVector & get_solution() const
    { return this->systems.front().solution; }
    BlockVector & get_solution()
    { return this->systems.front().solution; }
    const BlockVector & get_solution_m1() const
    { return this->systems.front().solution_m1; }
    const BlockVector & get_solution_m2() const
    { return this->systems.front().solution_m2; }
    unsigned int get_size() const
    { return this->systems.front().solution.length(); }
    const Example_TimeNSE & get_example() const
    { return example; }
    const ParameterDatabase & get_db() const
    { return db; }
    /// @brief return the computed errors at each discre time point
    std::array<double, int(10)> get_errors() const;

    /// @brief get the current residuals 
    /// @details updated in TimeNavierStokes<d>::compute_residuals() which in 
    /// turn is called from TimeNavierStokes<d>::stop_it().
    const Residuals& get_residuals() const;
    /// @brief get the current impuls residual
    /// @details updated in TimeNavierStokes<d>::compute_residuals() which in 
    /// turn is called from TimeNavierStokes<d>::stop_it().
    double get_impuls_residual() const;
    /// @brief get the current mass residual
    /// @details updated in TimeNavierStokes<d>::compute_residuals() which in 
    /// turn is called from TimeNavierStokes<d>::stop_it().
    double get_mass_residual() const;
    /// @brief get the current residual
    /// @details updated in TimeNavierStokes<d>::compute_residuals() which in 
    /// turn is called from TimeNavierStokes<d>::stop_it().
    double get_full_residual() const;
    /// @brief reset the residuals.
    /// use this if you want to use this object again for a second nonlinear 
    /// iteration.
    void reset_residuals();

    TimeDiscretization& get_time_stepping_scheme()
    { return time_stepping_scheme; }
    const TimeDiscretization& get_time_stepping_scheme() const
    { return time_stepping_scheme; }

    void add_to_output(const FEVectFunct* fe_vector_fct)
    { outputWriter.add_fe_vector_function(fe_vector_fct); }
    void add_to_output(const FEFunction* fe_fct)
    { outputWriter.add_fe_function(fe_fct); }
    void add_to_output(const std::string& name, std::function<double()> func)
    { outputWriter.add_global_output_variable(name, func); }

  protected:

    /// @brief default copy constructor (useful in derived classes)
    TimeNavierStokes(const TimeNavierStokes<d> &) = default;

    virtual void adjust_assembled_rhs();

    /** @brief store a complete system on a paticular grid.
     *
     * This combines a matrix, rhs, solution, spaces and functions
     * needed to describe a Time Navier-Stokes problem in 2D
     */
    struct System_per_grid
    {
      /** @brief Finite element space for velocity*/
      std::shared_ptr<FESpace> velocity_space;

      /** @brief Finite element space for pressure*/
      std::shared_ptr<FESpace> pressure_space;

      /** @brief system matrix                     |    [ A11  A12  A13  B1T ]
       *                       [ A11  A12  B1 ]    |    [ A21  A22  A23  B2T ]
       *                       [ A21  A22  B2 ]    |    [ A31  A32  A33  B3T ]
       *                       [ B3   B4   C  ]    |    [ B1   B2   B3   C   ]
      */
      BlockFEMatrix matrix;

      /** @brief mass matrix: this will be the standard mass matrix
       * for the standard Galerkin scheme. However, for the SUPG/RBVMS
       * schems, this includes all the terms which are related to the
       * time derivatives in the fully discrete scheme Eq: (45)
       * Ahmed, Rebollo, John and Rubino (2015)
       *  [ M11  M12  0 ]
       *  [ M21  M22  0 ]
       *  [ 0     0   0 ]
       * This additionaly created due to the struture of the
       * residual based Variational Multiscale Method:
       */
      BlockFEMatrix mass_matrix;

      /** @brief right hand side vector*/
      BlockVector rhs;

      /** @brief solution vector*/
      BlockVector solution;

      /** @brief Finite element function for velocity*/
      FEVectFunct u;

      /** @brief Finite element function for pressure*/
      FEFunction p;

      std::shared_ptr<PointwiseAssemblyData> persistent_data;

      /** @brief
       * old solution for the computation of the residual
       * that passes as an FEFunction to the local assembling
       * routines
       * 
       * @todo why don't we store a vector of size 
       * time_stepping_scheme.n_old_solutions() here?
       */
      BlockVector solution_m1;
      FEVectFunct u_m1;
      FEFunction p_m1;

      BlockVector solution_m2;
      FEVectFunct u_m2;
      FEFunction p_m2;

      double tau_m1;
      double tau_m2;

      BlockVector time_avg_sol;
      FEVectFunct u_time_avg;
      FEFunction  p_time_avg;

      BlockVector combined_old_sols;
      FEVectFunct comb_old_u;

      BlockVector extrapolate_sol;
      FEVectFunct extrapolate_u;
      FEFunction extrapolate_p;

      /** @brief constructor*/
      System_per_grid(const Example_TimeNSE& example, const TCollection& coll,
                      std::pair<int,int> order);

      /** @brief copy constructor*/
      System_per_grid(const System_per_grid&);

      //! Delete move constructor. No moves allowed.
      System_per_grid( System_per_grid&& ) = delete;

      //! Delete copy assignment operator. No copies allowed.
      System_per_grid& operator=( const System_per_grid& ) = delete;

      //! Default move assignment operator. No moves allowed.
      System_per_grid& operator=( System_per_grid&& ) = delete;

      //! Default destructor. Does most likely cause memory leaks.
      ~System_per_grid() = default;
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

    /** @brief a complete system on each grid
     *
     * Note that the size of this deque is at least one and larger only in case
     * of multigrid.
     */
    std::deque<System_per_grid> systems;

    /** @brief class for output handling */
    DataWriter<d> outputWriter;


    /** @brief Definition of the used example */
    Example_TimeNSE example;

    /** @brief a solver object which will solve the linear system
     *
     * Storing it means that for a direct solver we also store the factorization
     * which is usually not necessary.
     */
    Solver<BlockFEMatrix, BlockVector> solver;

    /** @brief an array to store defect, so that we don't have to reallocate
     *         so often
     */
    BlockVector defect;

    ///@brief The norms of residuals from up to 10 previous iterations
    FixedSizeQueue<10, Residuals> old_residuals;

    /** @brief store the initial residual so that the nonlinear iteration can
    *         be stopped as soon as a desired reduction is achieved
    */
    double initial_residual;

    int solve_count;

    /** @brief store errors  */
    std::array<double, 10> errors;

    BlockVector base_rhs;

    /** @brief right hand side vector from previous time step (on finest mesh) */
    BlockVector base_rhs_m1;

    /** @brief right hand side vector from two time steps ago (on finest mesh) */
    BlockVector base_rhs_m2;

    /// @brief time stepping scheme object to access everything
    TimeDiscretization time_stepping_scheme;
    /// @brief is the mass matrix and right hand side is solution
    /// dependent: for example the SUPG method
    bool is_rhs_and_mass_matrix_nonlinear;
    /// @brief system right hand side which passes to the solver
    BlockVector rhs_from_time_disc;

    /// this is for setting the discretization type globally, since the DISCTYPE 
    /// is removed fromt the global database but it's not simple/possible to use
    /// different space discretization in the LocalAssembling2D
    int space_disc_global;

    /// @brief to read/write the current solution as a checkpoint
    CheckpointIO checkpoint_io;

    /// @brief to write the current averaged solution as a checkpoint
    CheckpointIO time_average_io;

    /** @brief set parameters in database
    *
    * This functions checks if the parameters in the database are meaningful
    * and resets them otherwise. The hope is that after calling this function
    * this class is fully functional.
    *
    * If some parameters are set to unsupported values, an error occurs and
    * throws an exception.
    */
    void check_and_set_parameters();

    /** @brief get velocity and pressure space*/
    void get_velocity_pressure_orders(
      std::pair<int,int> &velocity_pressure_orders);

    /** @brief write some information (number of cells, dofs, ...) */
    void output_problem_size_info() const;

    /// @brief this routines wraps up the call to Assemble2D
    void call_assembling_routine(TimeNavierStokes<d>::System_per_grid& s,
                                 LocalAssembling_type type,
                                 bool first_nonlinear_step = false,
                                 bool do_upwinding = false,
                                 bool assemble_dirichlet_rows = false);

    using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
    using MatrixD = typename Template_names<d>::MatrixD; 
    /// @brief set the matrices and right hand side depending on the
    /// assemling routines, nstypes and the methods
    void set_matrices_rhs(TimeNavierStokes<d>::System_per_grid& s,
                          LocalAssembling_type type,
                          std::vector<SquareMatrixD*> &sqMat,
                          std::vector<MatrixD*> &reMat,
                          std::vector<double*> &rhs);

    /// @brief set the spaces depending on disc types
    void set_arrays(TimeNavierStokes<d>::System_per_grid& s,
                    std::vector<const FESpace*> &spaces,
                    std::vector<const FESpace*>& spaces_rhs,
                    std::vector<const FEFunction*>& functions,
                    bool first_nonlinear_step);  
    /// @brief restrict the function to every (coarse) grid
    /// nonlinear assembling requires an approximate velocity on every grid
    void restrict_function();

    /// update matrices for local projection stabilization
    void update_matrices_lps(TimeNavierStokes<d>::System_per_grid& s);

    /// @brief Special case for the SUPG and RBVMS method. The right-hand side is 
    /// nonlinear and need a complete re-assembling before passing to the solver 
    /// for solving.
    void assemble_rhs_nonlinear();

    /** @brief Assemble terms on boundary edges (2D) or faces (3D).
     */
    void assemble_boundary_terms(bool first_nonlinear_step);

    /** @brief modify matrices according to the Slip type boundary
     * conditions
     * If the mass matrix (also BT's are independent of solution )
     * is independent of solution then it only needs to modify only
     * once:
     * NOTE: mass matrix and BT's are solution dependent for
     * residual-VMS, and SUPG case.
     * Nonlinear matrices needs to be modify within each
     * time step and also within non-linear iteration
     */
    void modify_slip_bc(bool BT_Mass = false, bool slip_A_nl = false);

    void prepare_postprocessing(TCollection *coll);

    /** @brief LineEval object to store the lines where evalation has to be done
     */
    LinesEval<d> Lines;

    /** @brief calculate the time averaging of the solution in time_avg_sol */
    void time_averaging();

    /** @brief modify computed pressure due to the used nonlinear form
     * 
     *  Using the rotational or emac form of the nonlinear term leads to a 
     * pressure which includes an additional term 0.5*|u|^2. This modification 
     * is substracted, so that the resulting pressure is comparable to the ones
     * from other nonlinear forms.
     */
    void adjust_pressure();

    /** @brief projection space used for VMS method*/
    std::shared_ptr<FESpace> projection_space_;

    /** @brief finite element function for vms projection*/
    // can we rename it to large scales?? also check BlockVector!! currently just vector
    std::vector<double> vms_small_resolved_scales; 
    std::shared_ptr<FEVectFunct> vms_small_resolved_scales_fefct;

    /** matrices for turbulence model*/
    std::vector<std::shared_ptr<FEMatrix>> matrices_for_turb_mod;

    // piecewise constant space containing the labels of the local projection space
    std::shared_ptr<FESpace> label_for_local_projection_space_;

    // vector used for indicating local projection space 
    std::vector<double> label_for_local_projection;

    // finite element function for local projection space
    std::shared_ptr<FEFunction> label_for_local_projection_fefct;

    /** @brief method that will prepare the needed spaces and matrices
     * used for the Variational method
     */
    void prepare_vms(const TDomain& domain);

    /** @brief some common declaration for postprocessing
      * these are only on the finest grid
     */
    std::shared_ptr<FESpace> vorticity_space;
    unsigned int n_vort_dofs;
    std::vector<double> vorticity;
    std::shared_ptr<FEFunction> vorticity_funct;
    std::shared_ptr<FEFunction> divergence;
    double zero_vorticity;

    std::shared_ptr<FESpace> stream_function_space;
    unsigned int n_psi;
    std::vector<double> psi;
    std::shared_ptr<FEFunction> stream_function;

    std::shared_ptr<FESpace> wall_shear_stress_space;
    unsigned int n_wall_shear_stress_dofs;
    std::vector<double> wall_shear_stress;
    std::shared_ptr<FEVectFunct> wall_shear_stress_funct;

    BlockVector time_derivative;
    std::shared_ptr<FEFunction> p_time_derivative_funct;
    std::shared_ptr<FEVectFunct> u_time_derivative_funct;

    std::shared_ptr<FESpace> tke_space;
    unsigned int n_tke_dofs;
    std::vector<double> tke;
    std::shared_ptr<FEFunction> tke_funct;

    std::shared_ptr<FESpace> nu_eff_space;
    unsigned int n_nu_eff_dofs;
    std::vector<double> nu_eff;
    std::shared_ptr<FEFunction> nu_eff_funct;

    std::vector<WindkesselBoundary> wk_boundary;
    std::vector<double> wk_p_out;

    double last_kinetic_energy;
    double last_tke;
    double last_max_velocity;

 public:
  const FEFunction & get_vorticity_funct() const
    {return *vorticity_funct.get();}
};


#endif // __SYSTEM_TIMENAVIERSTOKES_H__
