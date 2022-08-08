/** ****************************************************************************
*
* @name       TNSE_POD
* @brief      Computation of POD basis - specific routines for TNSE problems
*
* Velocity POD basis functions and pressure POD basis functions are computed
* separately according to the decoupled approach.
* Data corresponding to the velocity or the pressure are respectively stored in
* files with 'u' or 'p' extension.
*
*******************************************************************************/

#ifndef TNSE_POD_H
#define TNSE_POD_H

#include "POD.h"

#include "BlockFEMatrix.h"
#include "BlockVector.h"
#ifdef __2D__
#include "Example_TimeNSE2D.h"
#include "FEFunction2D.h"
#include "FEVectFunct2D.h"
#else
#include "Example_TimeNSE3D.h"
#include "FEFunction3D.h"
#include "FEVectFunct3D.h"
#endif

#include "DataWriter.h"
#include "templateNames.h"
#include "TimeDiscretizations.h"
#include "AlgebraicFluxCorrection.h"
#include "LocalAssembling.h"


template<int d>
class TNSE_POD
{
  public:
    using FEFunction  = typename Template_names<d>::FEFunction;
    using FEVectFunct = typename Template_names<d>::FEVectFunct;
    using FESpace     = typename Template_names<d>::FESpace;
    using Example_TimeNSE = typename Template_names<d>::Example_TimeNSE;

    /** @brief constructor
    * this constructor calls the other constructor creating an
    * Example_TimeNSE object
    */
    TNSE_POD(TCollection&             coll,
             const ParameterDatabase& param_db);

    /** @brief The standard constructor, can be used for multigrid and
    *          non-multigrid.
    *
    * @param[in] collection collection of cells
    * @param[in] param_db   A parameter database with parameters concerning this
    *                       class or any of its members (fe space, solver,
    *                       assemble,...)
    * @param[in] example    The example which is to be calculated.
    */
    TNSE_POD(TCollection&             coll,
             const ParameterDatabase& param_db, 
             const Example_TimeNSE&   ex);

    /**
    * @brief return a database with all parameters necessary for
    * time-dependent Navier-Stokes (tnse) problems
    */
    static ParameterDatabase set_pod_basis_database(
                                             const ParameterDatabase& param_db);

    /**
    * @brief Compute POD basis
    *
    * Compute POD basis from snapshots with the compute_basis() routine
    * from the class POD. Within the function, appropriate gramian matrix
    * according to the 'pod_inner_product' parameter is assembled and
    * incorporated into computation of the POD basis.
    */
    void compute_pod_basis(const ParameterDatabase& param_db);

    /// @brief write POD basis to file and write vtk files
    void output(const std::vector<std::vector<double>>& snaps_avr,
                const std::vector<DenseMatrix>&         basis);

  protected:

    /** @brief suffix to differenciate velocity data from pressure data*/
    std::string data_suffix[2] = {"u", "p"};

    /** @brief Finite element space for velocity */
    std::shared_ptr<FESpace> velocity_space;

    /** @brief Finite element space for pressure */
    std::shared_ptr<FESpace> pressure_space;

    /** @brief used in output() to pass function to outputWriter */
    BlockVector pod_mode;

    /** @brief Finite element function for velocity */
    FEVectFunct velocity_function;

    /** @brief Finite element function for presure */
    FEFunction pressure_function;

    /** @brief a local parameter database which controls this class
    *
    * The database given to the constructor will be merged into this one. Only
    * parameters which are of interest to this class are stored (and the
    * default ParMooN parameters). Note that this usually does not include
    * other parameters such as solver parameters. Those are only in the
    * Solver object.
    */
    ParameterDatabase db;

    /** @brief Definition of the used example */
    const Example_TimeNSE example;

    /** @brief class for output handling */
    DataWriter<d> outputWriter;

    /** @brief check parameters in database
    *
    * This functions checks if the parameters in the database are meaningful
    * and resets them otherwise. The hope is that after calling this function
    * this class is fully functional.
    *
    * If some parameters are set to unsupported values, an error occurs and
    * throws an exception.
    */
    void set_parameters();

    /** @brief Assemble gramian matrix
    *
    * Assemble gramian matrix, which describes the inner product, with respect
    * to which the POD basis will be computed.
    */
    void assemble_gramian(BlockFEMatrix* mat_ptr, bool is_velocity, std::string inner_product);

    /** @brief print the problem infomation
    */
    void output_problem_size_info() const;


};



#ifdef __2D__
    /**
     * Named constructor for a matrix taking the block structure
     *
     * ( A  0  0 )
     * ( 0  A  0 )
     * ( 0  0  C )
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE3D POD-ROM problems
     */
BlockFEMatrix NSE2D_Diag(std::shared_ptr<const TFESpace2D> velocity,
                         std::shared_ptr<const TFESpace2D> pressure);

#else
    /**
     * Named constructor for a matrix taking the block structure
     *
     * ( A  0  0  0 )
     * ( 0  A  0  0 )
     * ( 0  0  A  0 )
     * ( 0  0  0  C )
     *
     * @param velocity The velocity finite element space.
     * @param pressure The pressure finite element space.
     * @return A newly constructed BlockFEMatrix for NSE3D POD-ROM problems
     */
BlockFEMatrix NSE3D_Diag(std::shared_ptr<const TFESpace3D> velocity,
                         std::shared_ptr<const TFESpace3D> pressure);
#endif
#endif // TNSE_POD_H
