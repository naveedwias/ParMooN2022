/** ****************************************************************************
*
* @name       TCD_POD
* @brief      Computation of POD basis - specific routines for TCD problems
*
*******************************************************************************/

#ifndef TCD_POD_H
#define TCD_POD_H

#include "POD.h"

#include "BlockFEMatrix.h"
#include "BlockVector.h"
#ifdef __2D__
#include "Example_TimeCD2D.h"
#include "FEFunction2D.h"
#else
#include "Example_TimeCD3D.h"
#include "FEFunction3D.h"
#endif

#include "DataWriter.h"
#include "templateNames.h"
#include "TimeDiscretizations.h"
#include "AlgebraicFluxCorrection.h"
#include "LocalAssembling.h"


template<int d>
class TCD_POD
{
  public:
    using FEFunction = typename Template_names<d>::FEFunction;
    using FESpace = typename Template_names<d>::FESpace;
    using Example_TimeCD = typename Template_names<d>::Example_TimeCD;

    /** @brief constructor
    * this constructor calls the other constructor creating an
    * Example_TimeCD object
    */
    TCD_POD(TCollection&             coll,
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
    TCD_POD(TCollection&             coll,
            const ParameterDatabase& param_db, 
            const Example_TimeCD&    ex);

    /** TODO: adapt comment (if more than tcd)
    * @brief return a database with all parameters necessary for
    * time-dependent convection-diffusion (tcd) problems
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
    void compute_pod_basis();

    /// @brief write POD basis to file and write vtk files
    void output();

    POD* get_pod_c()
    { return &pod_c; }

  protected:
    /** @brief Finite Element space */
    std::shared_ptr<FESpace> fe_space;

    POD pod_c;

    /** @brief gramian matrix (needed for POD computation)
     *
     * For euclidean-norm the gramian is the identity matrix and does not have
     * to be computed. For L2-norm it is equal to the mass matrix.
     */
    BlockFEMatrix gramian_matrix;

    /** @brief used in output() to pass function to outputWriter */
    BlockVector pod_mode;

    /** @brief Finite Element function */
    FEFunction fe_function;

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
    const Example_TimeCD example;

    /** @brief class for output handling */
    DataWriter<d> outputWriter;

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

    /** @brief Assemble gramian matrix
    *
    * Assemble gramian matrix, which describes the inner product, with respect
    * to which the POD basis will be computed.
    */
    void assemble_gramian();

    /** @brief Assemble stiffness matrix
    */
    void assemble_stiff(BlockFEMatrix& stiff_matrix);

    /** @brief Write spectral norm of the stiffness matrix of the POD basis
    */
    void write_S2();

    /** @brief print the problem infomation
    */
    void output_problem_size_info() const;
};

#endif // TCD_POD_H
