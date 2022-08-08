/** ****************************************************************************
*
* @name      TCD_ROM
* @brief     solve ROM for TCD (only time-independent BC)
*
*******************************************************************************/

#ifndef TCD_ROM_H
#define TCD_ROM_H

#include "ROM.h"

#include "BlockFEMatrix.h"
#include "BlockVector.h"
#include "DataWriter.h"
#include "deque"
#include "LocalAssembling.h"
#include "Solver.h"
#include "templateNames.h"
#include "TimeConvectionDiffusion.h"
#include "TimeDiscretizations.h"

#ifdef __2D__
#include "Example_TimeCD2D.h"
#include "FEFunction2D.h"
#else
#include "Example_TimeCD3D.h"
#include "FEFunction3D.h"
#endif


template<int d>
class TCD_ROM : public TimeConvectionDiffusion<d>, public ROM
{
  public:
    using FEFunction = typename Template_names<d>::FEFunction;
    using FESpace = typename Template_names<d>::FESpace;
    using Example_TimeCD = typename Template_names<d>::Example_TimeCD;

    /** @brief constructor
     *
     * This constructor calls the other constructor creating an
     * Example_TimeCD object
     */
    TCD_ROM(const TDomain& domain, const ParameterDatabase& param_db);

    /** @brief The standard constructor
     *
     * @param[in] Domain   The computational domain providing the grid
     * @param[in] param_db A parameter database with parameters concerning this
     *                     class or any of its members (fe space, solver,
     *                     assemble,...)
     * @param[in] example  The example which is to be calculated.
     */
    TCD_ROM(const TDomain&           domain,
            const ParameterDatabase& param_db,
            const Example_TimeCD&    example);

    /** @brief return a database with all parameters necessary for
     *        time-dependent convection-diffusion (tcd) problems with ROM
     */
    static ParameterDatabase set_databaseROM(const ParameterDatabase& param_db);

    /** @brief assemble all the matrices before the time iterations and call
     *         the compute_initial_condition() routine
     *
     * This function includes the assembling of: Stiff_matrix, Mass_Matrix,
     * (additional matrixK in case of SUPG stabilization) and rhs
     */
    void assemble_initial_time();

    /** @brief assemble the matrices
     *
     * This function will assemble the stiffness matrix (additional matrixK in
     * case of SUPG stabilization) and rhs.
     * In addition the system matrix and the rhs which passes to the solver
     * are also prepared by calling the function prepare_system_matrix_rhs()
     */
    void assemble();

    /** @brief solve the system
    */
    void solve();

    /** @brief measure errors and write solution
      * 
      */
    void output();

    const std::vector<double> get_sol_r() const
    { return sol_r; }

  protected:
    /** @brief a local parameter database which controls this class
    *
    * The database given to the constructor will be merged into this one. Only
    * parameters which are of interest to this class are stored (and the
    * default ParMooN parameters). Note that this usually does not include
    * other parameters such as solver parameters. Those are only in the
    * Solver object.
    */
    ParameterDatabase db;

    /** @brief For the time discretization in case of BDF2 scheme at the first
    *          step which need the use of an other time scheme (BDF1 or CN)
    */
    bool pre_stage_bdf = false;

    /** @brief coeffcients of the BDF2 scheme */
    std::vector<double> bdf2_coeffs;

    /** @brief flag indicating if stiff and mass matrices are time dependent*/
    bool mat_time_dependent;

    /** @brief Reduced system matrix */
    DenseMatrix sys_mat_r;

    /** @brief Reduced system rhs */
    std::vector<double> sys_rhs_r;

    /** @brief Reduced solution, i.e. ROM solution */
    std::vector<double> sol_r;

    /** @brief Reduced solution at the last time step */
    std::vector<double> sol_m1_r;

    /** @brief Reduced solution at the penultimate time step */
    std::vector<double> sol_m2_r;

    /** @brief Reduced mass matrix */
    DenseMatrix mass_mat_r;

    /** @brief Reduced convection-diffusion-reaction matrix (stiffness matrix)*/
    DenseMatrix cdr_mat_r;

    /** @brief Reduced convection-diffusion-reaction matrix times snaps mean
    *         only needed if the POD modes are computed from the snapshots'
    *         fluctuations (usefull for steady Dirichlet BC)*/
    std::vector<double> cdr_mat_mean_r;

    /** @brief Reduced source term (of its variational formulation)*/
    std::vector<double> rhs_r;

    /** @brief Reduced source term at the last time step */
    std::vector<double> rhs_m1_r;

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
    * With parameter 'rom_init_regularized' set to false, the ROM initial
    * condition is computed in the standard way: euclidean or L2 projection into
    * the POD space. Otherwise, the regularized ROM initial condition (see PhD
    * Thesis of S.Giere, Sec. 3.2.2.) is computed. The filter width in the
    * Helmholtz equation is controlled by the database parameter
    * 'differential_filter_width'.
    */
    void compute_initial_condition();

    /** @brief This wraps up the tedious call to Assemble and reduce FE terms
     *
     * @param type Flag indicating which kind of assembling has to be done.
     *             TCDMassOnly is used to compute the gramian matrix during
     *             assemble_initial_time().
     *             Custom is used to compute Mass and Stiffness matrices and Rhs
     *             during compute_initial_condition() if 'rom_init_regularized'
     *             is true to avoid adding SUPG terms to the Helmholtz equation.
     */
    void assemble_and_reduce(LocalAssembling_type type);

    /** @brief Compute system matrix and system rhs
     *
     * @note Note that assembling in every time step in very inefficient in the
     * ROM context. Try to avoid it, e.g., by storing all FE source terms at all
     * times and reduce them offline and online just use the reduced vectors.
     * Alternatively, for source terms given in the separated time-space form,
     * one could pre-assemble the space part and reduce it offline and within
     * the time loop just multiply it with the corresponding temporal
     * coefficient (e.g., see PhD Thesis of S.Giere, Sec. 4.2.).
     *
     * During the computation of system rhs the system is updated with the
     * function assemble_and_reduce(). It must be done before computing the
     * system matrix.
     * System rhs sys_rhs_r and system matrix sys_mat_r are computed.
     *
     * The functions prepare_system_matrix() and prepare_rhs_from_time_disc()
     * of TimeDiscretizations.C cannot be used, their parameter are
     * BlockFEMatrix instead of DenseMatrix.
     *
     * @param type Flag indicating which kind of assembling has to be done.
     *             Custom is used to only compute the system matrix at the 
     *             initial time in case of mat_time_dependent==false.
     */
    void prepare_system_matrix_rhs(LocalAssembling_type type);

};

#endif // TCD_ROM_H
