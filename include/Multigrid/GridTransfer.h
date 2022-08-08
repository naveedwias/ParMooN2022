/**
 * @file GridTransfer.h
 *
 * Gathers grid transfer operations needed for multgrid, formerly these were
 * declared in LinAlg.h and defined in MultigridComponents2D and
 * MultigridComponents3D.
 *
 * @date 2016/05/10
 * @author Clemens Bartsch
 */

#ifndef INCLUDE_MULTIGRID_GRIDTRANSFER_H_
#define INCLUDE_MULTIGRID_GRIDTRANSFER_H_

#include <vector>
#include <cstddef>

// Experimental Macro to avoid double code.
#ifdef __2D__
#define TFESpaceXD TFESpace2D
#elif __3D__
#define TFESpaceXD TFESpace3D
#endif

//forward declarations
class TFESpace2D;
class TFESpace3D;

namespace GridTransfer
{
/**
 * Prolongate a function from finer to coarser level.
 *
 * In MPI case you must make sure that the input CoarseFunction enters in level
 * 3 consistency. Due to its const-ness the update cannot be performed within
 * this function.
 *
 * @note No checks whether the spaces actually form the correct hierarchy are
 * performed - we rely on "garbage in, garbage out" here.
 */
void Prolongate(
    const TFESpaceXD& CoarseSpace, const TFESpaceXD& FineSpace,
    const double* CoarseFunction, size_t n_coarse_dofs,
    double* FineFunction, size_t n_fine_dofs);

/**
 * Restrict function from fine function to coarse function -
 * maps into dual space. Used for restricting the defect from fine grid to
 * coarse grid (where it will form the right hand side).
 * For an in-depth explanation see lecture notes of Volker John.
 *
 * In MPI case you must make sure that the input FineFunction enters in level
 * 3 consistency. Due to its const-ness the update cannot be performed within
 * this function.
 *
 * @note No checks whether the spaces actually form the correct hierarchy are
 * performed - we rely on "garbage in, garbage out" here.
 */
void DefectRestriction(
    const TFESpaceXD& CoarseSpace, const TFESpaceXD& FineSpace,
    double* CoarseFunction, size_t n_coarse_dofs,
    const double* FineFunction, size_t n_fine_dofs);

/**
 * Restricts the function on the finest grid to all other grids
 * successively.
 *
 * @note No checks whether the spaces form an actual hierarchy are
 * performed - we rely on "garbage in, garbage out" here.
 *
 * @param[in] space_hierarchy An ordered hierarchy of FE Spaces, finest first.
 * @param[in,out] function_entries The functions to be filled.
 * @param[in] The lengths of the functions. Must match the FESpaces number of dofs.
 *
 * TODO First function should be const.
 */
void RestrictFunctionRepeatedly(
  const std::vector<const TFESpaceXD*>& space_hierarchy,
  const std::vector<double*>& function_entries,
  const std::vector<size_t>& function_n_dofs);


/** Restrict a function from coarse to fine level.
 * @note No checks whether the spaces actually form the correct hierarchy are
 * performed - we rely on "garbage in, garbage out" here.*/
void RestrictFunction(
    const TFESpaceXD& CoarseSpace, const TFESpaceXD& FineSpace,
    double* CoarseFunction, size_t n_coarse_dofs,
    const double* FineFunction, size_t n_fine_dofs);
}
#undef TFESpaceXD
#endif /* INCLUDE_MULTIGRID_GRIDTRANSFER_H_ */
