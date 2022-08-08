/**
 * @file Implementation of class MultigridLevel.
 *
 * @date 2016/05/10
 * @author Clemens Bartsch
 */

#include <BlockFEMatrix.h>
#include <DirectSmoother.h>
#include <JacobiSmoother.h>
#include <MooNMD_Io.h>
#include <MultigridLevel.h>
#include <ParameterDatabase.h>
#include <VankaSmoother.h>
#include <SORSmoother.h>

#ifdef _MPI
#include <ParFECommunicator3D.h>
#endif
#include <cmath>

MultigridLevel::MultigridLevel(BlockFEMatrix* matrix,
                               SmootherCode sm,
                               const ParameterDatabase& db)
:  matrix_(matrix),
   defect_(*matrix, true), residual_(1e10),
   rhs_(*matrix,true), solution_(*matrix, false),
   smoother_(nullptr)
{
  Output::info<5>("MultigridLevel", "Constructed a MultigridLevel object. matrix "
               "dimensions: (", matrix->get_n_total_rows(), ",",
               matrix->get_n_total_columns(), "), n_cells ",
               matrix->get_ansatz_space(0,0)->GetCollection()->GetN_Cells());

  //Determine which smoother to use and construct the object.
  switch(sm)
  {
    case SmootherCode::DIRECT_SOLVE:
      smoother_ = std::make_shared<DirectSmoother>();
      break;
    case SmootherCode::JACOBI:
    {
      double damp = db["jacobi_damp"];
      smoother_ = std::make_shared<JacobiSmoother>(damp);
      break;
    }
    case SmootherCode::SOR:
    { // sor smoother with backward sweep ("0"), overrelaxation parameter omega
      // determined by input db
      double omega = db["sor_omega"];
#ifdef _MPI
      smoother_ = std::make_shared<SORSmoother>(omega, 0, db["sor_parallel_type"]);
#else
      smoother_ = std::make_shared<SORSmoother>(omega, 0);
#endif
      break;
    }
    case SmootherCode::SSOR:
    { // sor smoother with backward sweep ("0"), overrelaxation parameter omega
      // determined by input db
      double omega = db["sor_omega"];
#ifdef _MPI
      smoother_ = std::make_shared<SORSmoother>(omega, 2, db["sor_parallel_type"]);
#else
      smoother_ = std::make_shared<SORSmoother>(omega, 2);
#endif
      break;
    }
    case SmootherCode::NODAL_VANKA:
    {
      double damp = db["multigrid_vanka_damp_factor"];
      smoother_ = std::make_shared<VankaSmoother>(VankaType::NODAL, damp);
      break;
    }
    case SmootherCode::CELL_VANKA:
    {
      double damp = db["multigrid_vanka_damp_factor"];
      smoother_ = std::make_shared<VankaSmoother>(VankaType::CELL, damp);
      break;
    }
    case SmootherCode::PATCH_VANKA:
    {
      double damp = db["multigrid_vanka_damp_factor"];
      smoother_ = std::make_shared<VankaSmoother>(VankaType::PATCH, damp);
      break;
    }
    case SmootherCode::CELL_VANKA_JACOBI:
    {
      double damp = db["multigrid_vanka_damp_factor"];
      smoother_ = std::make_shared<VankaSmoother>(VankaType::CELL_JACOBI, damp);
      break;
    }
    //Vanka smoothers with storage of local systems
    case SmootherCode::NODAL_VANKA_STORE:
    {
      double damp = db["multigrid_vanka_damp_factor"];
      smoother_ = std::make_shared<VankaSmoother>(VankaType::NODAL, damp, true);
      break;
    }
    case SmootherCode::CELL_VANKA_STORE:
    {
      double damp = db["multigrid_vanka_damp_factor"];
      smoother_ = std::make_shared<VankaSmoother>(VankaType::CELL, damp, true);
      break;
    }
    case SmootherCode::PATCH_VANKA_STORE:
    {
      double damp = db["multigrid_vanka_damp_factor"];
      smoother_ = std::make_shared<VankaSmoother>(VankaType::PATCH, damp, true);
      break;
    }
    default:
      ErrThrow("Unknown SmootherCode!");
  }
}


void MultigridLevel::apply_smoother()
{
#ifdef _MPI
  // establish level 3 consistency of rhs and solution
  std::vector<const TParFECommunicator3D*> comms = matrix_->get_communicators();
  for(size_t bl =0; bl < comms.size(); ++bl)
  {
    comms[bl]->consistency_update(rhs_.block(bl), 3);
    comms[bl]->consistency_update(solution_.block(bl), 3);
  }
#endif
  smoother_->smooth(rhs_, solution_);
}

void MultigridLevel::calculate_defect()
{
  defect_ = rhs_;

#ifdef _MPI
  std::vector<const TParFECommunicator3D*> comms = matrix_->get_communicators();
  for(size_t bl =0; bl < comms.size(); ++bl)
  {
    comms[bl]->consistency_update(solution_.block(bl), 3); //restore level 3 consistency of solution_
  }
#endif

  matrix_->apply_scaled_add(solution_, defect_, -1.0);

#ifdef _MPI
  residual_ = defect_.norm(comms);
#else
  residual_ = std::sqrt(dot(defect_,defect_));
#endif
}

void MultigridLevel::update_smoother()
{
  smoother_->update(*matrix_);
}
