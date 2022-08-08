#include "CheckpointIO.h"
#include "MooNMD_Io.h"

/// @brief only in MPI mode. This makes sure we read the right file on each
/// processor.
constexpr auto parallel_extension = ".proc";
void modify_file_name(std::string& base_name);

CheckpointIO::CheckpointIO (const ParameterDatabase& input_db)
{
  ParameterDatabase db = CheckpointIO::default_checkpoint_io_database();
  db.merge(input_db);
  
  write_solution = db["write_solution_binary"];
  read_solution = db["read_initial_solution"];
  write_solution_n_steps = db["write_solution_binary_all_n_steps"];
  overwrite_solution = db["overwrite_solution_binary"];
  base_name = db["write_solution_binary_file"].get<std::string>();
  base_name_for_reading = db["initial_solution_file"].get<std::string>();
}

ParameterDatabase CheckpointIO::default_checkpoint_io_database()
{
  ParameterDatabase db("default ParMooN checkpoint io database");

  // to read or not to read
  db.add("read_initial_solution", false, "Choose true if the initial "
         "solution is given in a binary file. Do not forget to specify "
         "'initial_solution_file' in that case, too.", {true, false});

  // from which file to read
  db.add("initial_solution_file", "my_solution_in.txt",
         "If 'read_initial_solution' is set to 'true', this parameter "
         "determines from which binary file to read the initial solution. That "
         "file should be prodcued by a former run of the same program, using "
         "the same finite element space.");

  // to write or not to write
  db.add("write_solution_binary", false, "Choose true if the computed solution "
         "should be written out to a file in a binary format. This is helpful "
         "if you plan to read it in as initial solution later. Do not forget "
         "to specify 'write_solution_binary_all_n_steps' for the output "
         "interval and 'write_solution_binary_file', the file path and name.",
         {true,false});

  // write all n steps
  db.add("write_solution_binary_all_n_steps", 1u, "Determine at which interval "
         "a backup solution should be written to the file "
         "'write_solution_binary_file'. The number refers to the number of "
         "timesteps, i.e., at a constant time step length of 0.01s and this "
         "paramter set to 10, the current solution will be written into the "
         "file each 0.1s.", 1u, 1000000u);

  // into which file to write
  db.add("write_solution_binary_file", "parmoon_solution_binary",
         "If 'write_solution_binary' is set to 'true', this parameter sets the "
         "name (and path) of the file to write the solution into.");

  // whether the file should be overwritten every time
  db.add("overwrite_solution_binary", false, " Choose true, if the solution "
         "binary should be overwritten every time. If it stays false, a new "
         "file will be written every time, with the current time added to the "
         "filename.", {true, false} );

  return db;
}

void CheckpointIO::write(const BlockVector& data,
                         const TimeDiscretization& time_stepping_scheme)
const
{
  if(!write_solution)
    return;
  
  const unsigned int step = time_stepping_scheme.current_step_;
  const bool initial_time = (step == 0);
  const bool write_this_data = (step % write_solution_n_steps == 0);
  if(write_this_data && !initial_time) // skip initial time
  {
    std::string file_name = base_name;
    if(!overwrite_solution)
    {
      file_name += "." + std::to_string(time_stepping_scheme.current_time_);
    }
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if(my_rank == 0)
      Output::print<3>("writing checkpoint files ", file_name,
                       parallel_extension, "*");
#endif
    modify_file_name(file_name);
#ifndef _MPI
    Output::print<3>("writing checkpoint file ", file_name);
#endif
    data.write_to_file(file_name);
  }
}

void CheckpointIO::write(const BlockVector& data) const
{
  if(!write_solution)
    return;
  std::string file_name = base_name;
#ifdef _MPI
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if(my_rank == 0)
      Output::print<3>("writing checkpoint files ", file_name,
                       parallel_extension, "*");
#endif
  modify_file_name(file_name);
#ifndef _MPI
  Output::print<3>("writing checkpoint file ", file_name);
#endif
  data.write_to_file(file_name);
}

bool CheckpointIO::read(BlockVector& data) const
{
  if(read_solution)
  {
    std::string file_name = base_name_for_reading;
    modify_file_name(file_name);
    data.read_from_file(file_name);
  }
  return read_solution;
}

void modify_file_name(std::string& base_name)
{
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  base_name += parallel_extension + std::to_string(my_rank);
#endif
  (void) base_name; // avoid compiler warning
  (void) parallel_extension; // avoid compiler warning
}
