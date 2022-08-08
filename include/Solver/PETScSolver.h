/** ************************************************************************ 
*
* @class PETScSolver
* @brief solving a linear system using PETSc
*
*
* \page page_using_petsc Using PETSc
* \brief Some examples how to use PETSc
*
* PETSc solves your problems with preconditioned Krylov subspace methods.
* You can have an overview of all Krylov methods and preconditioners at
* http://www.mcs.anl.gov/petsc/documentation/linearsolvertable.html
*
* The default Krylov method is GMRES with ILU as preconditioner.
*
* PETSc is configured by setting two things in your input
* *.dat file:
*   - set 'solver_type: petsc'
*   - set 'petsc_arguments:' to the usual PETSc command line
*     arguments as found in the PETSc documentation
*     https://www.mcs.anl.gov/petsc/
*
*
* \section solving_scalar Solving scalar problems
* set 'petsc_arguments:' as
* \code{.sh}
* -ksp_type richardson -ksp_monitor -pc_type lu -pc_factor_mat_solver_package umfpack
* \endcode
*
*
* \section solving_saddle_point Solving saddle point problems
* set 'petsc_arguments:' as
* \code{.sh}
* -ksp_type fgmres -pc_type fieldsplit -pc_fieldsplit_type schur -fieldsplit_0_ksp_atol 1.0e-13 -fieldsplit_0_ksp_rtol 0. -fieldsplit_1_ksp_atol 1.0e-13 -fieldsplit_1_ksp_rtol 0.
* \endcode
*
*
* PETSc is able to handle systems of equations in block form, similar to what
* a BlockMatrix in ParMooN does. It uses a special matrix representation for 
* such a case. This is the only way to used specialized saddle point solvers in
* PETSc. Furthermore for other systems of multiple equations, one has a very 
* fine control of which solvers/preconditioners are used for each equation 
* (subsystem). However there seem to be some limitations. Saddle point problems
* require a two-by-two system which is not what we do in ParMooN for 
* Navier-Stokes equations. And in case of a multiple equations you can use a 
* direct solver for the entire system. If you want to use direct solvers for 
* system of multiple equations, you should use the DirectSolver class in 
* ParMooN.
*
****************************************************************************/

#ifndef __PETSCSOLVER__
#define __PETSCSOLVER__

#include "all_defines_external_libraries.h"
#ifdef PARMOON_WITH_PETSC
#include <petscksp.h> // defines type 'Mat' below and much more
#include <petscpc.h> // defines type 'Mat' below and much more
#endif // PARMOON_WITH_PETSC
#include <ParameterDatabase.h>
#include <vector>
#include <memory>

// forward declaration
class FEMatrix;
class BlockMatrix;
class BlockFEMatrix;
class BlockVector;
template <class V>
class CompositeOperator;

#ifdef _MPI
class TParFECommunicator3D;
#endif // _MPI

class PETScSolver
{
  public:
    PETScSolver(const BlockFEMatrix& matrix, const ParameterDatabase& db);

    /** @brief constructor: an internal copy is created and stored. Ready to 
     * solve.
     */
    PETScSolver(const BlockMatrix& matrix, const ParameterDatabase& db);

    PETScSolver(const CompositeOperator<BlockVector>& matrix,
      const ParameterDatabase& db);

    /** @brief This class is not copy constructible */
    PETScSolver(const PETScSolver&) = delete;

    /** @brief move constructor */
    PETScSolver(PETScSolver&&);

    /** @brief This class is not copy assignable */
    PETScSolver& operator=(const PETScSolver&) = delete;

    /** @brief move assignment operator */
    PETScSolver& operator=(PETScSolver&&);
    
    /** @brief release all memory */
    ~PETScSolver();

    /**
     * @brief Solves the equation A*(solution)=rhs for solution.
     * 
     * @param rhs the right-hand side of the problem Ax=b
     * @param solution vector to store the solution into
     */
    void solve(const BlockVector& rhs, BlockVector& solution);
    
    void set_prefix(std::string s)
    { prefix = s; }
    
  private:
    /** @brief some identifying string to prefix output. */
    std::string prefix;
#ifdef PARMOON_WITH_PETSC
    /** @brief the matrix of the linear equation A*x=b, PETSc format */
    Mat petsc_mat;
    
    /** @brief in case of a saddle point problem these are sub blocks 
     * 
     * This is similar to the BlockMatrix in ParMooN. Some specialized PETSc 
     * solvers require this format rather than one big matrix.
     */
    std::vector<Mat> sub_petsc_mats;
    
    /// @brief PETSc linear solver object (Krylov subspace method)
    KSP ksp;
#endif // PARMOON_WITH_PETSC
#ifdef _MPI
    std::vector<const TParFECommunicator3D*> comms_;
#endif // _MPI
};
#endif // __PETSCSOLVER__

