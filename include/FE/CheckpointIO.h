#ifndef CHECKPOINTWRITER_H
#define CHECKPOINTWRITER_H

#include "ParameterDatabase.h"
#include "BlockVector.h"
#include "TimeDiscretizations.h"
#include <string>

/**
 * @brief enable writing of solution vectors to binary as a checkpoint.
 * 
 * Whenever longer simulations are interrupted for some unforeseeable reason, we
 * want to be able to restart it with the last solution as initial data.
 * Typically, this is done for a time loop, but could as well be extended for
 * any other iteration. This class provides the interface to writing a solution
 * to files which can be read in again, for example in a restarted simulation.
 * 
 * Possible/future extensions:
 *    - read/write HDF5 files instead of binary, leads to a single file even in
 *      MPI mode
 *    - read/write more than one BlockVector, useful for multistep time 
 *      discretizations and coupled problems with multiple solutions
 *    - move the reading/writing from the BlockVector class to this class
 * 
 * @ruleof0
 */
class CheckpointIO
{
  public:
    /// @brief constructor reading parameters from database to fill members
    CheckpointIO(const ParameterDatabase& input_db);
    
    CheckpointIO(const CheckpointIO&) = default;
    CheckpointIO(CheckpointIO&&) = default;
    CheckpointIO& operator=(const CheckpointIO&) = default;
    CheckpointIO& operator=(CheckpointIO&&) = default;
    ~CheckpointIO() = default;
    
    /// @brief A database to control reading/writing of checkpoints.
    static ParameterDatabase default_checkpoint_io_database();
    
    /// @brief write the 'data' to a binary file.
    /// This does nothing if CheckpointIO::write_solution is false or if
    /// the number of completed time steps (TimeDiscretization::current_step_)
    /// is not a multiple of CheckpointIO::write_solution_n_steps.
    void write(const BlockVector& data,
               const TimeDiscretization& time_stepping_scheme) const;
    /// @brief write the 'data' to a binary file.
    /// This does nothing if CheckpointIO::write_solution is false. This is
    /// intended for stationary problems, in the time-dependent case, use 
    /// CheckpointIO::write(data, time_stepping_scheme).
    void write(const BlockVector& data) const;
    /// @brief read data from file.
    /// The file is determined during construction by the parameter 
    /// 'initial_solution_file'. The input 'data' is only changed if the 
    /// parameter 'read_initial_solution' is true.
    /// @return return true if data has been read, otherwise nothing is done.
    bool read(BlockVector& data) const;
    
  protected:
    /// @brief write binary solutions if set to true, otherwise the method 
    /// CheckpointIO::write() does nothing.
    bool write_solution;
    /// @brief read solution from file if set to true, otherwise the method
    /// CheckpointIO::read() does nothing.
    bool read_solution;
    /// @brief write a solution only in every n-th step
    unsigned int write_solution_n_steps;
    /// @brief overwrite a previously written file, otherwise new files will be
    /// created (with a filename suffix).
    bool overwrite_solution;
    /// @brief prefix for all files written by this class object
    std::string base_name;
    /// @brief the base name used when reading from a file to a BlockVector
    std::string base_name_for_reading;
};

#endif // CHECKPOINTWRITER_H
