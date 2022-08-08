#ifndef ALGEBRAIC_FLUX_CORRECTION_H_
#define ALGEBRAIC_FLUX_CORRECTION_H_

#include <ParameterDatabase.h>

#include <vector>
#include <BlockVector.h>
#include <BlockFEMatrix.h>

#ifdef __2D__
#include "FEFunction2D.h"
#else
#include "FEFunction3D.h"
#endif


class FEMatrix;

namespace AlgebraicFluxCorrection {
  
/**
 * control parameter which limiter to use (in steady-state problem so far)
 * KUZMIN: limiter from Zalesak (1979) and Kuzmin (2007)
 * BJK17: limiter from Barrenechea, John, Knobloch (2017)
 * MONOLITHIC: limiter from Kuzmin (2020)
 */
enum class Limiter
{
    KUZMIN, BJK17, MONOLITHIC, MONOLITHIC_STEADY, MUAS, MUAS_Kno21, MUAS_MAX, MUAS_MAX_ABS 
};

/**
 * control parameter for the iteration scheme for the afc problem 
 * (in steady-state problem so far)
 * FIXEDPOINT_RHS: fixed point iterations with changes of the right-hand side
 * FIXEDPOINT_MATRIX: fixed point iterations with changes of the matrix
 * NEWTON: Newton's method
 * NEWTON_REGU: Newton Regualrized method; only implemented for KUZMIN limiter
 */
enum class Iteration_Scheme
{
    FIXEDPOINT_RHS, FIXEDPOINT_MATRIX, NEWTON, NEWTON_REGU
};
/**
 * Sets up and returns a default algebraic flux correction database,
 * which contains all control parameters necessary for an algebraic
 * flux correction. All of them are initialized with default values.
 */
ParameterDatabase default_afc_database();

/**
 * FEM-TVD ('Total Variation Diminishing') algorithm for steady-state
 * convection-diffusion-reaction problems, following
 *
 * Kuzmin, D.(2007): Algebraic Flux Correction for Finite Element
 * Discretizations of Coupled Systems.
 *
 * Note that for the algorithm to operate properly, the system matrix
 * must have been constructed as if it only had active degrees of freedom.
 * This is not checked, because when reaching here the matrix is a purely
 * algebraic object without knowledge of FE spaces.
 * To regain the correct Dirichlet rows, one has to call
 * correct_dirichlet_rows upon the system matrix afterwards.
 *
 * @param[in,out] system_matrix The matrix to undergo algebraic flux correction.
 *                         Must be square.
 * @param[in] diffusion_matrix Diffusion part of the matrix
 * @param[in] sol The current solution vector.
 * @param[out] rhs The current right hand side vector. Gets modified.
 * @param[in] neum_to_diri array which contains the indices of actual
 *  Dirichlet dof which were but treated as Neumann dof. CB: That paramerter
 *  is unused at the moment!
 * @param[in] coeffs Bilinear coeffecients
 * @param[in] afc_matrix_D_entries entries of the matrix D
 * @param[in] afc_matrix_D_B matrix D-B used for Modified Kuzmin
 * @param[in] gamma weights needed in the linearity preserving limiter from 
 * [BJK17]
 * @param[in] compute_D_and_gamma boolean that gives the information  of the 
 *  matrix D and the vector gamma need to be computed 
 * @param[in] limiter Specifies the used limiter
 * @param[in] it_scheme Specifies the used iteration scheme
 * @param[in] is_not_afc_fixed_point_rhs Specifies if the scheme is 
 *            fixed_point_rhs or not; required so as to avoid recomputation of
 *  A+D
 */
void steady_state_algorithm(FEMatrix& system_matrix, 
                            FEMatrix& diffusion_matrix, 
                            const std::vector< double >& sol, 
                            std::vector< double >& rhs, 
                            BlockVector& poisson_sol,
                            const std::vector< int >& neum_to_diri, 
                            double* coeffs,
                            FEMatrix& D, FEMatrix& D_B, 
                            std::vector< double >& gamma, 
                            std::vector< double >& alphas, 
                            bool compute_D_and_gamma, 
                            const ParameterDatabase& db, 
                            Limiter limiter = AlgebraicFluxCorrection::Limiter::KUZMIN, 
                            Iteration_Scheme it_scheme = AlgebraicFluxCorrection::Iteration_Scheme::FIXEDPOINT_MATRIX, 
                            const int is_not_afc_fixed_point_rhs = 1);


/**
 * Sets all Dirichlet rows of the input Matrix to 0 and 1 on diagonal.
 * In the assembling of a matrix which is supposed to be treated with
 * AFC all dofs are treated as active. So after the AFC completed, the
 * actual dirichlet rows have to be reset. This is achieved here.
 *
 * @param[in,out] MatrixA The matrix to be corrected. Must be square.
 */
void correct_dirichlet_hanging_rows(FEMatrix& MatrixA);

void correct_dirichlet_hanging_rhs(FEMatrix& A, BlockVector& RHS);




/** 
 * Computes the new iterate for the AFC schemes
 * @param[in] old_solution The solution from the previous iteration
 * @param[in,out] new_solution Input: solution of the linear system, would be
 * the new iterate if there is no damping
 * 
 * Projection to Admissible values, in case the method doesn't satisfy the DMP.
 * NOTE: Only applicable if the maximum and minimum of the solution is known.
 * 
 */
 void AFC_Compute_New_Iterate(const BlockVector& old_solution, 
                              BlockVector& new_solution, 
                              const ParameterDatabase& db,
                              FEMatrix& matrix);
}

 //Stiffness entries for A
 static std::vector<double> original_stiffness_Entries;
 ///raw fluxes f_{ij} in [BJK16.SINUM]
 static std::vector<double> raw_fluxes;
 ///Convective part of the matrix
 static std::vector<double> convective_matrix;
 ///Artificial diffusion matrix entries
 static std::vector<double> artificial_diffusion_entries;

#endif /* ALGEBRAIC_FLUX_CORRECTION_H_ */
