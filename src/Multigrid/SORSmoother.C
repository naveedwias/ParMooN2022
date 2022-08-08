/**
 * @file SORSmoother.C
 * Implementation of class SORSmoother.
 *
 * @date 2016/09/09
 * @author Clemens Bartsch
 */

#include <SORSmoother.h>
#include <BlockVector.h>
#include <BlockFEMatrix.h>
#include <Iteration_sor.h>
#ifdef _MPI
#include <ParFECommunicator3D.h>
#endif

#ifdef _MPI
SORSmoother::SORSmoother(double omega, size_t sor_strat, const std::string& sor_par_strat)
: sor_(nullptr), omega_(omega), sor_strat_(sor_strat), parallel_strategy_(sor_par_strat)
{
}
#else
SORSmoother::SORSmoother(double omega, size_t sor_strat)
: sor_(nullptr), omega_(omega), sor_strat_(sor_strat)
{
}
#endif

void SORSmoother::smooth(const BlockVector& rhs, BlockVector& solution)
{
  //MPI: solution and rhs enter as level 3 consistent
  sor_->apply_smoother(rhs, solution); //one sor sweep - like in preconditioning
}

void SORSmoother::update(const BlockFEMatrix& matrix)
{
  //Reset the sor object, w.
  sor_.reset(new Iteration_sor<BlockFEMatrix, BlockVector>(matrix, sor_strat_, omega_));

#ifdef _MPI
  sor_->set_parallel_strategy(parallel_strategy_);
#endif
}


