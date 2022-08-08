/*
 * VankaSmoother.C
 *
 *  Created on: May 16, 2016
 *      Author: bartsch
 */

#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <DenseMatrix.h>
#include <DirectSolver.h>
#include <MooNMD_Io.h>
#include <VankaSmoother.h>

#include <memory>

#ifdef _MPI
#include <mpi.h>
#include <ParFECommunicator3D.h>
#include "BaseCell.h"
#endif

// Experimental Macro to avoid double code.
#ifdef __2D__
#define TFESpaceXD TFESpace2D
#elif __3D__
#define TFESpaceXD TFESpace3D
#endif

//! Default constructor.
VankaSmoother::VankaSmoother(VankaType type, double damp_factor, bool store)
: type_(type), dimension_(0), damp_factor_(damp_factor),
  matrix_global_(nullptr), press_dofs_local_(0), velo_dofs_local_(0),
  local_systems_(0), store_systems_(store), pressure_space_(nullptr),
  velocity_space_(nullptr)
{
#ifdef _MPI
  if(type != VankaType::CELL_JACOBI && type != VankaType::CELL)
    ErrThrow("Among the Vanka smoothers only cell_vanka and cell_vanka_jacobi"
        " are enabled for MPI so far.");
#endif
}

void VankaSmoother::update(const BlockFEMatrix& matrix)
{
  //Check if matrix looks like saddle point problem (i.e.: all but last block
  // row space are the same)
  size_t n_blocks = matrix.get_n_cell_rows();

  //First space will be interpreted as velocity space, last space as pressure space.
  const TFESpaceXD* first_space = matrix.get_test_space(0,0).get();
  const TFESpaceXD* last_space = matrix.get_test_space(n_blocks - 1 , 0).get();

  for(size_t i = 0; i < n_blocks - 1 ; ++i)
  {
    if(matrix.get_test_space(i,0).get() != first_space)
    {
      ErrThrow("So far VankaSmoother will only operate for matrices, whose"
          " first n-1 blocks have the same fe space."
          " This matrix does not fulfil that requirement! ", i);
    }
  }
  //Find out if spaces changed, and if so: reassort dof patches
  //pressure (when pressure space has changed)
  // @todo TODO the following two "if" dont check correctly if the spaces changed
  // in comparison to the previous one. This needs to be corrected to make a safer
  // comparison. For the moment, I dont comment them out, because then the tests are
  // much longer, especially in tnse3d case: reassort when needed (with "if") ~ 20s,
  // reassort systematically (without "if") ~ 315s (Najib)!!!
  if(last_space != pressure_space_)
  {
    //TODO Check that cell vanka and continuous space are not used together!
    // Trouble is: Space does not know whether it holds discontinuous elements.
    pressure_space_ = last_space;
    set_up_pressure_patches(*pressure_space_);
  }
  //velocity ( when either space changed)
  if(first_space != velocity_space_ || last_space != pressure_space_)
  {
    velocity_space_ = first_space;
    bool is_transposed{true}; //out-variable for get_block
    std::shared_ptr<const TMatrix> coupling_block
      = matrix.get_block(n_blocks - 1, 0, is_transposed); //the block B1
    if(is_transposed)
      ErrThrow("Oy vey! That coupling block is transposed, "
          "that's not what I expected.")
    set_up_velocity_patches(*coupling_block.get(), *velocity_space_);
  }

  //Reset the stored global matrix on which all the work is done
  dimension_ = n_blocks - 1;
  matrix_global_=matrix.get_combined_matrix();

  //Throw away all local systems.
  for(auto sys : local_systems_)
  {
    delete sys; sys = nullptr;
  }
  local_systems_ = std::vector<DenseMatrix*>(press_dofs_local_.size(),nullptr);

#ifdef _MPI
    // store the communicator structure belonging to the matrices' FESpaces
    comms_ = matrix.get_communicators();
#endif
}

/// @todo TODO By now it's too many different strategies in here, it is not anymore
/// intelligible - divide the method into parts!
//The implementation of the smoothing step is very procedural in nature...
void VankaSmoother::smooth(const BlockVector& rhs, BlockVector& solution )
{

  // Check input
  if(rhs.n_blocks() != dimension_ + 1)
    ErrThrow("VankaSmoother: rhs dimension does not fit!");

  if(solution.n_blocks() != dimension_ + 1)
    ErrThrow("VankaSmoother: solution dimension does not fit!");

#ifdef _MPI
  // Update solution to level 3 consistency.
  for (size_t bl = 0; bl < comms_.size();++bl)
  {
    comms_.at(bl)->consistency_update(solution.block(bl),3);
  }
#endif

  // These vectors will be used in the CELL_JACOBI strategy to collect updates.
  BlockVector updates_jacobi(solution);
  updates_jacobi.reset(); //set to 0
  BlockVector update_counts(updates_jacobi); //actually stores only ints
  update_counts.reset(); //set to 0
  // For Jacobi, global residual has to be calculated only once per sweep.
  std::vector<double> residual_glo(rhs.get_entries_vector());
  matrix_global_->multiply(solution.get_entries(), &residual_glo.at(0), - 1.0);

  // Loop over all local systems and do all the work
  for(size_t i = 0; i < press_dofs_local_.size() ; ++i)
  {
    const DofPatch& velo_dofs = velo_dofs_local_.at(i);
    const DofPatch& press_dofs = press_dofs_local_.at(i);
    size_t n_velo_dofs_local = velo_dofs.getSize();
    size_t n_press_dofs_local = press_dofs.getSize();

    size_t n_velo_dofs_global = this->velocity_space_->get_n_dof();

    size_t size_local = dimension_* n_velo_dofs_local + n_press_dofs_local;

    //read: dof_map[i] = global dof corresponding to local dof i (added over all block rows)
    std::vector<int> dof_map(size_local, 0);

    //local solution and rhs
    std::vector<double> rhs_local(size_local, 0.0);
    std::vector<double> solution_local(size_local, 0.0);

    /* ******** Fill dof map. *********** */
    for(size_t k = 0; k < size_local; ++k)
    {
      if(k < dimension_* n_velo_dofs_local) //is a velo dof
      {
        int k_modulo_n_velo_dofs_local = k%n_velo_dofs_local;
        int block = k/n_velo_dofs_local; //integer division
        dof_map.at(k) = block * n_velo_dofs_global
            + velo_dofs.getDof(k_modulo_n_velo_dofs_local);
      }
      else //is a pressure dof
      {
        int block = dimension_;
        int k_minus_dim_times_velo = k - dimension_*n_velo_dofs_local;
        dof_map.at(k) = block * n_velo_dofs_global
            + press_dofs.getDof(k_minus_dim_times_velo);
      }
    }

    if (local_systems_.at(i) == nullptr )
    {
      /* ******** Set up local matrix *********** */
      DenseMatrix* matrix_local = new DenseMatrix(size_local,size_local);
      for(size_t r_loc = 0; r_loc < size_local; ++r_loc) //loop over all local rows
      {
        size_t r_glo = dof_map.at(r_loc); //the corresponding global dof

        int begin_r_glo = matrix_global_->get_row_ptr()[r_glo];
        int end_r_glo = matrix_global_->get_row_ptr()[r_glo + 1];

        size_t c_loc = 0;// start with the 0th local column - exploit that
        // the global KCol Array and the dof_map array are sorted

        for( int j = begin_r_glo ; j < end_r_glo ; ++j )
        {
          double entry = matrix_global_->GetEntries()[j];
          if(entry != 0) //don't copy zeroes.
          {
            int c_glo = matrix_global_->get_vector_columns()[j];
            //find out what the local column is
            while(dof_map.at(c_loc) < c_glo)
            {
              ++c_loc;
              if(c_loc == size_local)
                break;
            }
            if(c_loc == size_local)
              break; //break loop, we're behind the end!
            if(dof_map.at(c_loc) == c_glo) //the dof c_glo is of interest for the local system
            {
              matrix_local->setEntry(r_loc, c_loc, entry);
            }
            //else just go on
          }
        }

      }
      //decompose and store
      matrix_local->decomposeLU();
      local_systems_.at(i) = matrix_local;
    }

    /* ******** Set up local right hand side. *********** */
    /* (Local right hand side is global defect in local rows) */
    for(size_t r_loc = 0; r_loc < size_local; ++r_loc) //loop over all local rows
    {
      if(type_ == VankaType::CELL_JACOBI)
      {//in jacobi case, the residual has been calculated globally already
        rhs_local.at(r_loc) = residual_glo.at(dof_map.at(r_loc));
      }
      else
      {
        size_t r_glo = dof_map.at(r_loc); //the corresponding global dof
        int begin_r_glo = matrix_global_->get_row_ptr()[r_glo];
        int end_r_glo = matrix_global_->get_row_ptr()[r_glo + 1];

        double temp = rhs.get_entries()[r_glo];

        for( int i = begin_r_glo ; i < end_r_glo ; ++i )
        {
          int c_glo = matrix_global_->get_vector_columns()[i];
          double matrix_glo_entry = matrix_global_->GetEntries()[i];
          temp -= matrix_glo_entry * solution.get_entries()[c_glo];
        }

        rhs_local.at(r_loc) = temp;
      }
    }

    /* ******** Solve the local system with LAPACK. *********** */
    //copy rhs_local into solution_local
    solution_local = rhs_local;
    local_systems_.at(i)->solve(&solution_local.at(0));

    if (!store_systems_)
    {
      delete local_systems_[i]; local_systems_[i] = nullptr;
    }

    /* ******** Add damped local solution to global solution. *********** */
    for(size_t r_loc = 0; r_loc < size_local; ++r_loc) //loop over all local rows
    {
      double damp = damp_factor_;
      size_t r_glo = dof_map.at(r_loc); //the corresponding global dof
      if(type_ == VankaType::CELL_JACOBI)
      {
        // Store the updates in update_jacobi and count in update_counts.
        updates_jacobi.get_entries()[r_glo] += damp*solution_local.at(r_loc);
        update_counts.get_entries()[r_glo] += 1;
      }
      else
      {
        solution.get_entries()[r_glo] += damp*solution_local.at(r_loc);
      }
    }
  }

  if(type_ == VankaType::CELL_JACOBI)
  {
#ifdef _MPI
    // Updates and numbers are stored additive. Make the vectors consistent.
    for (size_t bl = 0; bl < comms_.size();++bl)
    {
      comms_.at(bl)->update_from_additive_to_consistent_storage(update_counts.block(bl),0);
      comms_.at(bl)->update_from_additive_to_consistent_storage(updates_jacobi.block(bl),0);
    }
#endif
    for(size_t i = 0; i < updates_jacobi.length(); ++i)
    {
      if(update_counts[i] != 0) //don't divide by zero...
      {
        solution[i] += updates_jacobi[i]/(double)update_counts[i];
      }
    }
  }

}

/*!
 * @brief Assembling of the Pressure Dofs.
 *
 * @param pressureSpacePtr A pointer to the pressure space the dofs refer to.
 *
 * This method determines the size, the setup and the order of the pressure dof patches.
 * Playing with it renders different Vanka smoothers!
 */
void VankaSmoother::set_up_pressure_patches(const TFESpace& pressureSpace){

  switch (type_) {
    case VankaType::NODAL:
    { // pressure node oriented Vanka
      //There are as many patches as nodes in this case.
      int nPatches = pressureSpace.get_n_dof();
      //Loop over cells.
      for (int i=0;i<nPatches;i++){
        //To every node belongs a pressure patch. Create that patch here.
        DofPatch currentPressurePatch{};
        currentPressurePatch.addDof(i);
        //Make the pressure patch nice and clean and copy it into the list.
        currentPressurePatch.tidyUp();
        press_dofs_local_.push_back(currentPressurePatch);
      }
      //End loop over cells
    } break;

    case VankaType::CELL:
    case VankaType::CELL_JACOBI:
    case VankaType::PATCH:
    { //cell, cell_jacobi and patch oriented Vanka treat the pressure dofs the same.

      //There are as many patches as cells in this case.
      int nPatches = pressureSpace.GetN_Cells();

      //Loop over cells.
      for (int i=0;i<nPatches;i++){
#ifdef _MPI
        // In MPI case choose only own cells to form local systems.
        if(pressureSpace.GetCollection()->GetCell(i)->IsHaloCell())
          continue;
#endif
        //To every cell belongs a pressure patch. Create that patch here.
        DofPatch currentPressurePatch{};
        //Loop over cell dofs.
        for(int j=0, n_local_p_dof = pressureSpace.get_n_local_dof(i);
            j < n_local_p_dof; ++j)
        {
          //Put the current dof into current pressure patch
          currentPressurePatch.addDof(pressureSpace.get_global_dof(i, j));
        }
        // End loop over cell dofs.

        //Make the pressure patch nice and clean and copy it into the list.
        currentPressurePatch.tidyUp();
        press_dofs_local_.push_back(currentPressurePatch);
      }
      //End loop over cells
    }
    break;

    default:
      ErrThrow("Unknown or unimplemented Vanka smoother type! ");

  }

}

/*!
 * @brief Assembling of the Velocity Dofs.
 *
 * @param[in]   pressureVelocityMatrix A matrix block from which the velo-pressure coupling
 *      can be determined (via non-zero entries).
 * @param[in] velocitySpace the velocity space which is "needed" for by-the-book assembling in CELL case.
 *
 * Look for non-zero entries in the matrix pressureVelocityMatrix.
 * To each row (crsp. pressure dof) put the columns (crsp. velo dofs)
 * into the corresponding velocity dof patch.
 *
 * This method is the same for all nodal vankas.
 */
void VankaSmoother::set_up_velocity_patches(const TMatrix& pressureVelocityMatrix,
                                            const TFESpace& velocitySpace){
  switch (type_) {
    case VankaType::NODAL:
    case VankaType::PATCH: //nodal and cellpatch Vanka treat assorting of velocity patches the same.
    {
      // Hold pointers for lookup in the coupling matrix.
      const int* columnsPointer = pressureVelocityMatrix.GetStructure().get_vector_columns();
      const int* rowPointer = pressureVelocityMatrix.GetStructure().get_row_ptr();

      // Loop over pressure patches.
      for(auto pPatchesIterator : press_dofs_local_){
        // Construct the new corresponding velocity patch.
        DofPatch currentVeloPatch{};
        //Loop over pressure dofs in current pressure patch
        for (auto pDof : pPatchesIterator){
          //Loop over all the places in columnsPointer, where entries are connected to row pDof.
          for (int i = rowPointer[pDof]; i != rowPointer[pDof+1]; i++){
            //Every column index that appears in the sparsity structure
            //is a velo dof connected to the pressure dof.
            currentVeloPatch.addDof(columnsPointer[i]);
          }

        }
        //End loop over pressure dofs in current pressure patch.

        //Make the velocity patch nice and clean and add a copy to the list.
        currentVeloPatch.tidyUp();
        velo_dofs_local_.push_back(currentVeloPatch);
      }
      // End loop over pressure patches
      break;
    }
    case VankaType::CELL:
    case VankaType::CELL_JACOBI:
    { //for cell vanka gather all velo dofs in one cell into one patch
      //There are as many patches as cells in this case.
      int nPatches = velocitySpace.GetN_Cells();

      //Loop over cells.
      for (int i=0;i<nPatches;i++){
#ifdef _MPI
        // In MPI case choose only own cells to form local systems.
        if(velocitySpace.GetCollection()->GetCell(i)->IsHaloCell())
            continue;
#endif
        //To every cell belongs a pressure patch. Create that patch here.
        DofPatch currentVeloPatch{};
        //Loop over cell dofs.
        for(int j=0, n_local_u_dof = velocitySpace.get_n_local_dof(i);
            j < n_local_u_dof; ++j)
        {
          //Put the current dof into current pressure patch
          currentVeloPatch.addDof(velocitySpace.get_global_dof(i, j));
        }
        // End loop over cell dofs.

        //Make the pressure patch nice and clean and copy it into the list.
        currentVeloPatch.tidyUp();
        velo_dofs_local_.push_back(currentVeloPatch);
      }
      //End loop over cells

      break;
    }
    default:
    {
      ErrThrow("Unknown or unimplemented Vanka smoother type! ");
      break;
    }

  }
}


VankaSmoother::~VankaSmoother()
{
  //Delete all local systems.
  for(auto sys : local_systems_)
  {
    delete sys; sys = nullptr;
  }
}

#undef TFESpaceXD
