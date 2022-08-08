/** ****************************************************************************
*
* @name      TNSE_ROM
* @brief     solve ROM for TNSE
*
* @note      only implemented for time-independent boundary conditions
*
* Velocity POD basis functions are supposed to be (discretely) divergence-free,
* thus the pressure term in the momentum equation drops out and the computation
* of the velocity is decoupled of the pressure. The pressure is recovered with a
* stabilization-based equation (SM-ROM approach).
* First the reduced velocity is computed from the velocity POD basis functions,
* then the reduced pressure is computed from the pressure POD basis functions.
*
* For details, see:
* "Numerical and Analytical Aspects of POD-Based Reduced-Order Modeling in
* Computational Fluid Dynamics", S. Giere, 2016.
*
*******************************************************************************/

#ifndef TNSE_ROM_H
#define TNSE_ROM_H

#include "ROM.h"

#include "BlockFEMatrix.h"
#include "BlockVector.h"
#include "DataWriter.h"
#include "deque"
#include "LocalAssembling.h"
#include "Solver.h"
#include "templateNames.h"
#include "TimeDiscretizations.h"
#include "TimeNavierStokes.h"

#ifdef __2D__
#include "Example_TimeNSE2D.h"
#include "FEFunction2D.h"
#else
#include "Example_TimeNSE3D.h"
#include "FEFunction3D.h"
#endif


template<int d>
class TNSE_ROM : public TimeNavierStokes<d>
{
  public:
    using FEFunction = typename Template_names<d>::FEFunction;
    using FEVectFunct = typename Template_names<d>::FEVectFunct;
    using FESpace = typename Template_names<d>::FESpace;
    using Example_TimeNSE = typename Template_names<d>::Example_TimeNSE;
    using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
    using MatrixD = typename Template_names<d>::MatrixD;
    using BdValFunct = typename Template_names<d>::BoundaryValuesFunction;
    using BdCondFunct = typename Template_names<d>::BoundaryConditionFunction;

    /** @brief constructor
     *
     * This constructor calls the other constructor creating an
     * Example_TimeNSE object
     */
    TNSE_ROM(const TDomain& domain, const ParameterDatabase& param_db);

    /** @brief The standard constructor
     *
     * @param[in] Domain   The computational domain providing the grid
     * @param[in] param_db A parameter database with parameters concerning this
     *                     class or any of its members (fe space, solver,
     *                     assemble,...)
     * @param[in] example  The example which is to be calculated.
     */
    TNSE_ROM(const TDomain&           domain,
             const ParameterDatabase& param_db,
             const Example_TimeNSE&   ex);

    /** @brief return a database with all parameters necessary for
     *         time-dependent Navier Stokes problems with ROM
     *
     * @param[in] param_db A parameter database with parameters concerning this
     *                     class or any of its members (fe space, solver,
     *                     assemble,...)
     * @param[in] is_full  flag is set to true when used with constructor
     *                     TimeNavierStokes() in order to take account of
     *                     all parameters
     */
    static ParameterDatabase set_databaseROM(const ParameterDatabase& param_db,
                                             bool is_full=false);

    /** @brief assemble all the matrices and rhs before the time iterations
     *         and call the compute_initial_condition() routine
     */
    void assemble_initial_time();

    /** @brief assemble the matrices and rhs
     *
     * This function will assemble the rhs at each first nonlinear step.
     * In addition it calls the function prepare_system_velocity() where the
     * convective nonlinear term is updated.
     *
     * @param it_counter nonlinear iteration number
     */
    void assemble_matrices_rhs(unsigned int it_counter);

    /** @brief solve the system
     *
     * @param[in] is_velocity flag to select if velocity (true)
     *                        or pressure (false) has to be computed
     */
    void solve(bool is_velocity);

    /** @brief check if one of the stopping criteria is fulfilled
     *
     * either converged or maximun number of iterations reached
     *
     * @param it_counter current iterate
     */
    bool stop_it(unsigned int it_counter);

    /** @brief get the current reduced velocity residual
     */
    double get_velocity_residual() const
    { return residual_u_r; }

    /** @brief print out the residual of the corresponding fem system computed
     *         from the parent class TimeNavierStokes<d>
     */
    void print_fem_residuals(int loop_index);


  private:

    /** @brief ROM related to the velocity components */
    ROM ROM_u;
    /** @brief ROM related to the pressure */
    ROM ROM_p;

    /** @brief a local parameter database which controls this class
     *
     * The database given to the constructor will be merged into this one. Only
     * parameters which are of interest to this class are stored (and the
     * default ParMooN parameters). Note that this usually does not include
     * other parameters such as solver parameters. Those are only in the
     * Solver object.
     */
    ParameterDatabase db;

    /** @brief rank of the velocity POD basis*/
    int rank_u;

    /** @brief Reduced system matrix */
    DenseMatrix sys_mat_u_r;

    /** @brief Linear part of the reduced system matrix */
    DenseMatrix lhs_mat_u_r;

    /** @brief Reduced system rhs */
    std::vector<double> sys_rhs_u_r;

    /** @brief Linear parts of the reduced system rhs */
    DenseMatrix rhs_mat_u_r;

    /** @brief Reduced mass matrix */
    DenseMatrix mass_mat_u_r;

    /** @brief Reduced diffusion matrix*/
    DenseMatrix diff_mat_u_r;
    /** @brief Contribution due to the mean fluctuation correction*/
    std::vector<double> diff_mean_u_r;

    /** @brief Reduced convection matrix (linear part)*/
    DenseMatrix conv_lin_mat_u_r;
    /** @brief Contribution due to the mean fluctuation correction*/
    std::vector<double> conv_lin_mean_u_r;

    /** @brief Reduced convection matrix (nonlinear part)*/
    std::vector<DenseMatrix> conv_non_lin_mat_u_r;
    /** @brief Contribution due to the mean fluctuation correction*/
    std::vector<std::vector<double>> conv_non_lin_mean_u_r;

    /** @brief Reduced source term*/
    std::vector<double> rhs_src_u_r;

    /** @brief Reduced solution, i.e. ROM solution */
    std::vector<double> sol_u_r;

    /** @brief Reduced solution at the last time step */
    std::vector<double> sol_m1_u_r;

    /** @brief Reduced source term at the last time step for Crank Nicolson*/
    std::vector<double> rhs_src_m1_u_r;

    /** @brief Reduced solution at the penultimate time step for BDF2*/
    std::vector<double> sol_m2_u_r;

    /** @brief rank of the pressure POD basis*/
    int rank_p;

    /** @brief Reduced system matrix */
    DenseMatrix sys_mat_p_r;

    /** @brief Reduced system rhs */
    std::vector<double> sys_rhs_p_r;

    /** @brief Reduced diffusion matrix*/
    DenseMatrix diff_mat_p_r;

    /** @brief Reduced matrix (u, grad q)*/
    DenseMatrix div_mat_p_r;
    /** @brief Contribution due to the mean fluctuation correction*/
    std::vector<double> div_mean_p_r;

    /** @brief Reduced convection matrix (linear part)*/
    DenseMatrix conv_lin_mat_p_r;
    /** @brief Contribution due to the mean fluctuation correction*/
    std::vector<double> conv_lin_mean_p_r;

    /** @brief Reduced convection matrix (nonlinear part)*/
    std::vector<DenseMatrix> conv_non_lin_mat_p_r;
    /** @brief Contribution due to the mean fluctuation correction*/
    std::vector<std::vector<double>> conv_non_lin_mean_p_r;

    /** @brief Reduced source term*/
    std::vector<double> rhs_src_p_r;

    /** @brief Reduced solution, i.e. ROM solution */
    std::vector<double> sol_p_r;

    /** @brief reduced velocity residual */
    double residual_u_r;

    /** @brief coeffcients of the BDF2 scheme */
    std::vector<double> bdf2_coeffs;

    /** @brief For the time discretization in case of BDF2 scheme at the first
     *         step which need the use of an other time scheme (BDF1 or CN)
     */
    bool pre_stage_bdf = false;

    /** @brief length of the last time step, used in prepare_system_velocity()
     *         to check if the linear parts have to be actualized
     */
    double tau_old = 0.;

    /** @brief flag indicating if the corresponding FEM residuals have
     *         to be computed and printed out.
     */
    bool fem_residuals = false;

    /** @brief check and set parameters in database
     *
     * This functions checks if the parameters in the database are meaningful
     * and resets them otherwise. The hope is that after calling this function
     * this class is fully functional.
     *
     * If some parameters are set to unsupported values, an error occurs and
     * throws an exception.
     */
    void check_and_set_parameters();

    /** @brief Compute ROM initial condition
     *
     * The ROM initial condition is computed in the standard way: euclidean, L2
     * or H1 projection into the POD space. In case of L2 or H1 projection,
     * the gramian matrix (inner product matrix) is temporary stored in the
     * BlockFEMatrix mass_matrix (matrix respectively) from the parent class
     * TimeNavierStokes for the velocity
     * or a temporary BlockFEMatrix for the pressure
     */
    void compute_initial_condition();

    /** @brief Assemble gramian matrice (inner products) of velocity or pressure
     *
     * For velocity:
     *   L2: LocalAssembling_type::TimeNavierStokesMass
     *   H1: LocalAssembling_type::NavierStokesLinear
     *
     * For pressure:
     *   L2: LocalAssembling_type::TCDMassOnly
     *   H1: LocalAssembling_type::TCDDiffOnly
     * 
     * @param[in] BlockFEMatrix* when this pointer is not nullptr, assemble the
     *                           gramian for pressure in the pointed matrix
     *                           otherwise assemble the gramian for velocity.
     * @return a shared pointer to the TMatrix containing the gramian.
     */
    std::shared_ptr<TMatrix> assemble_gramian(BlockFEMatrix* mat_ptr=nullptr);

    /** @brief set the parameters for calling the function assemble_and_reduce()
     *         according to the LocalAssembling type and the unknown (velocity
     *         or pressure)
     *
     * For each contribution (i.e. diffusion, convection, rhs, ...) different
     * LocalAssembling_type are used in combination with the flag is_velocity:
     * For the velocity:
     *   Assembling mass: LocalAssembling_type::TimeNavierStokesMass
     *   Assembling diff: LocalAssembling_type::NavierStokesLinear
     *   Assembling rhs (source): LocalAssembling_type::TimeNavierStokesRhs
     *   Assembling conv_lin          \
     *   Assembling conv_lin_mean (BC) \
     *                                 >LocalAssembling_type::TimeNavierStokesNL
     *   Assembling conv_non_lin       /
     *   Assembling conv_non_lin_mean /
     * For the pressure:
     *   Assembling diff: LocalAssembling_type::Custom
     *   Assembling div: LocalAssembling_type::NavierStokesLinear
     *   Assembling rhs (source): LocalAssembling_type::TimeNavierStokesRhs
     *   Assembling conv_lin          \
     *   Assembling conv_lin_mean (BC) \
     *                                 >LocalAssembling_type::TimeNavierStokesNL
     *   Assembling conv_non_lin       /
     *   Assembling conv_non_lin_mean /
     * The different LocalAssembling_type are used in combination with the flag
     * is_velocity as parameters selector.
     *
     * @param type       Flag indicating which kind of assembling has to be done
     *                   TimeNavierStokesAll or TimeNavierStokesRhs
     * @param is_velocity Flag indicating if the velocity or pressure system has
     *                    to be assembled
     */
    void assemble_system(LocalAssembling_type type, bool is_velocity);

    /** @brief This wraps up the tedious call to assemble and reduce FE terms
     *
     * The LocalAssembling_type is used in combination with the flag is_velocity
     * as parameters selector for the functions set_arrays(), set_matrices_rhs()
     * and reduce(). It is modified with the function adapt_type() when passed
     * to the LocalAssembling constructor in order to 'approach' the expected
     * LocalAssembling object, in this way calling the LocalAssembling custom
     * constuctor is avoided, but the resulting LocalAssembling object has to be
     * adjusted with the function adapt_local_assembling() afterwards.
     * This proceeding is meant to be more intelligible and reliable as coping a
     * large part of LocalAssembling::set_parameters_for_tcd() and
     * LocalAssembling::set_parameters_for_tnse().
     *
     * @param is_reduced Flag selecting assemble only or assemble and reduce
     * @param type       Flag indicating which kind of assembling has to be done
     * @param is_velocity which unknown (velocity or pressure)
     * @param pressure_matrix pointer to the matrix to be assembled
     */
    void assemble_and_reduce(bool                 is_reduced,
                             LocalAssembling_type type,
                             bool                 is_velocity,
                             BlockFEMatrix*       pressure_matrix=nullptr);

    /** @brief set the spaces and fe functions depending on the type (and the
     *         unknown, but enough information is still contained in the type)
     */
    void set_arrays(std::vector<const FESpace*>&    spaces,
                    std::vector<const FESpace*>&    spaces_rhs,
                    std::vector<const FEFunction*>& functions,
                    std::vector<BdCondFunct*>&      bdCond,
                    std::vector<BdValFunct*>&       bdVal,
                    LocalAssembling_type            type,
                    int                             new_term_idx);

    /** @brief set the matrices and right hand side depending on the
     *         assemling type and the unknown (velocity or pressure)
     */
    void set_matrices_rhs(std::vector<SquareMatrixD*>& sqMat,
                          std::vector<MatrixD*>&       reMat,
                          std::vector<double*>&        rhs_array,
                          LocalAssembling_type         type,
                          bool                         is_velocity,
                          BlockFEMatrix*               pressure_matrix);

    /** @brief adapt the LocalAssembling_type for the LocalAssembling
     *         constructor.
     *
     * Until now LocalAssembling_type was only used as parameters selector in
     * combination with the flag is_velocity.
     * The adapted LocalAssembling_type have been chosen for the proximity of
     * the resulting LocalAssembling objects (obtained with
     * LocalAssembling::set_parameters_for_tcd()
     * or LocalAssembling::set_parameters_for_tnse()) and the expected
     * one, adjustments are done with adapt_local_assembling() afterwards.
     */
    LocalAssembling_type adapt_type(bool is_velocity,LocalAssembling_type type);

    /** @brief Adjust the default LocalAssembling object, according to the
     *         velocity or pressure ROM specificities.
     */
    void adapt_local_assembling(bool                 is_velocity,
                                LocalAssembling_type type,
                                LocalAssembling<d>&  la);

    /** @brief reduce the assembled matrices and rhs.
     */
    bool reduce(LocalAssembling_type type,
                bool                 is_velocity,
                int&                 new_term_idx,
                BlockFEMatrix*       pressure_matrix);

    /** @brief Compute the complete velocity system matrix and rhs according to
     *         the time discretisation
     *
     * @note Functions prepare_system_matrix() and prepare_rhs_from_time_disc()
     *       of TimeDiscretizations.C cannot be used, their parameter are
     *       BlockFEMatrix instead of DenseMatrix.
     *
     * This function calls prepare_offline_velocity() if necessary (once at the
     * beginning, and for each first nonlinear step if the step size changed)
     *
     * @param first_nonlinear_step Flag indicating if the iteration is currently
     *                             the first nonlinear step.
     */
    void prepare_system_velocity(bool first_nonlinear_step);

    /** @brief Compute the velocity linear parts lhs and rhs
     */
    void prepare_offline_velocity(void);

    /** @brief Compute the pressure system matrix and rhs according to the
     *         velocity–pressure ROM based on a stabilization of the coupled
     *         problem (SM-ROM) as described in:
     *         "A numerical investigation of velocity–pressure reduced order
     *         models for incompressible flows",
     *         A. Caiazzo, T. Iliescu, V. John, S. Schyschlowa,
     *         Journal of Computational Physics 259 (2014) 598–616
     */
    void prepare_system_pressure(void);

    /** @brief Compute the defect Ax-b and store its norm.
     */
    void compute_residuals(void);

};

#endif // TNSE_ROM_H
