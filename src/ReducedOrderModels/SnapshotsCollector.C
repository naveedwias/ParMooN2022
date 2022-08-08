/** ****************************************************************************
*
* @name       SnapshotsCollector
* @brief      write/read snapshot data into/from a file
*
* @author     Swetlana Giere & Alfonso Caiazzo
* @date       08.05.2012 (start of implementaion). Restarted on 15.1.2019
*
*******************************************************************************/

#include <SnapshotsCollector.h>
#include "CheckpointIO.h"
#include <sys/stat.h>

#ifdef _MPI
#include <mpi.h>
#endif

/** ***************************************************************************/
ParameterDatabase SnapshotsCollector::default_snapshots_database()
{
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("Default ParMooN parameter database for storing snapshots");

  db.add("snaps_directory", ".",
         "This directory is where the snapshots will be written. This "
         "directory will be created, if it does not exist already. Files in "
         "this directory will be overwritten without any warning.");

  db.add("snaps_basename", "parmoon_snapshots",
         "Basename for file where the snapshots are stored.");

  db.add("steps_per_snap", (size_t) 5,
         "This integer specifies how many time steps are performed "
         "before a snapshot has to be written into file.");

  db.add("snaps_time_derivative", false,
         "Flag to save the time derivative in the snapshots. If true, for each "
         "output time step but the first one, the snapshots will contain the "
         "solution and its time derivative", { true, false });

  db.add("snaps_time_derivative_projection", false,
         "Set this to true if your time derivative snapshots are stored as "
         "linear functionals and need to be projected back into the FE space.",
         { true, false });

  db.add("snaps_time_derivative_tau", 1.0,
         "Scaling factor for time derivative snapshots.",
         0.0, 1.e+6);

  // Merge used to check if the solution has to be written in binary format
  db.merge(CheckpointIO::default_checkpoint_io_database(), true);

  return db;
}

/** ***************************************************************************/
SnapshotsCollector::SnapshotsCollector(const ParameterDatabase& param_db,
                                       std::string filename_suffix)
                   :  db(SnapshotsCollector::default_snapshots_database()),
                      snap_count(0)
{
  this->db.merge(param_db, false);
  this->snapshot_filename = SnapshotsCollector::get_snapshot_filename(this->db, filename_suffix);

  if (db["write_snaps"])
  {
    // create directory db["snaps_directory"]
    std::string directory_name = this->db["snaps_directory"].get<std::string>();
    mkdir(directory_name.c_str(), 0777);

    if (db["write_solution_binary"])
    {
      Output::print<1>("Filename for storing snapshots: ", snapshot_filename);
      this->datastream.open(snapshot_filename.c_str(),
                            std::ios::binary | std::ios::trunc | std::ios::out);
    }
    else
    {
      Output::print<1>("Filename for storing snapshots: ", snapshot_filename);
      this->datastream.open(snapshot_filename.c_str(),
                            std::ios::trunc | std::ios::out);
      this->datastream << setprecision(12);
    }

    if (!this->datastream.good())
    {
      ErrThrow("Error: File ", snapshot_filename,
               " could not be opened in SnapshotsCollector.\n"
               "(No read access to file or file already open)");
    }
  }
}

std::string SnapshotsCollector::get_snapshot_filename(const ParameterDatabase& db, std::string filename_suffix)
{
  std::string snapshot_filename = db["snaps_directory"].get<std::string>() + "/";
  snapshot_filename += db["snaps_basename"].get<std::string>();

#ifdef _MPI
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  snapshot_filename += "_r" + std::to_string(mpi_rank);
#endif

  if (!filename_suffix.empty())
  {
    snapshot_filename += "." + filename_suffix;
  }

  if (db["write_solution_binary"])
  {
    snapshot_filename += ".bin";
  }

  return snapshot_filename;
}

/** ***************************************************************************/
SnapshotsCollector::~SnapshotsCollector()
{
  if (this->datastream.is_open())
  {
    this->datastream.close();
  }
}

/** ***************************************************************************/
bool SnapshotsCollector::write_data(const double* const solution,
                                    const unsigned int length,
                                    size_t time_step)
{
  if ((time_step % (int)db["steps_per_snap"]) == 0)
  {
    Output::print<1>("Time step ", time_step, ", writing snapshot n. ",
                     snap_count + 1, " on ", snapshot_filename);
    snap_count++;

    if (!this->datastream.is_open())
    {
      ErrThrow("Error: Snapfile is not open.");
    }

    if (db["write_solution_binary"])
    {
      try
      {
        this->datastream.write((char*) &(solution[0]), sizeof(double) * length);

        if (!this->datastream.good())
        {
          throw std::runtime_error("Snapshot was not correctly written.");
        }
      }
      catch (std::exception& e)
      {
        ErrThrow(e.what());
      }
    }
    else
    {
      for (size_t i = 0; i < length; ++i)
      {
        this->datastream << solution[i] << " ";
      }

      this->datastream << endl;
    }

    this->datastream.flush();
    return true;
  }
  else
  {
    return false;
  }
}

/** ***************************************************************************/
void SnapshotsCollector::write_data(const double* const solution,
                                    const unsigned int length,
                                    size_t time_step,
                                    const double* const sol_m1,
                                    double tau)
{
  bool is_to_be_written = write_data(solution, length, time_step);

  if (is_to_be_written && db["snaps_time_derivative"])
  {
    if (tau == 0.0)
    {
      ErrThrow("Error write_data(): time step length tau cannot be 0");
    }

    if (db["write_solution_binary"])
    {
      std::vector<double> sol_dt(length);
      for (size_t i = 0; i < length; ++i)
      {
        sol_dt[i] = (solution[i] - sol_m1[i]) / tau;
      }
      this->datastream.write((char*) &(sol_dt[0]),
                             sizeof(double) * length);
    }
    else
    {
      double sol_dt;
      for (size_t i = 0; i < length; ++i)
      {
        sol_dt = (solution[i] - sol_m1[i]) / tau;
        this->datastream << sol_dt << " ";
      }

      this->datastream << endl;
    }

    this->datastream.flush();
  }
}

void SnapshotsCollector::write_data(const double* const solution,
                                    const unsigned int length,
                                    size_t time_step,
                                    const double* const dt_solution)
{
  bool is_to_be_written = write_data(solution, length, time_step);

  if (is_to_be_written && db["snaps_time_derivative"])
  {
    if (db["write_solution_binary"])
    {
      this->datastream.write((char*) dt_solution,
                             sizeof(double) * length);
    }
    else
    {
      for (size_t i = 0; i < length; ++i)
      {
        this->datastream << dt_solution[i] << " ";
      }

      this->datastream << endl;
    }

    this->datastream.flush();
  }
}

/** ***************************************************************************/
void SnapshotsCollector::read_data(std::string                        filename,
                                   bool                               bin,
                                   int                                length,
                                   std::vector<std::vector<double>>& mat)
{
  std::ifstream file;

  if (bin)
  {
    file.open(filename.c_str(), std::ios::in | std::ios::binary);
  }
  else
  {
    file.open(filename.c_str(), std::ios::in);
  }

  if (!file.good())
  {
    ErrThrow("Error: File ", filename, " could not be opened");
  }

  if (bin)
  {
    // get length of file:
    file.seekg(0, file.end);
    uint_fast64_t file_length = file.tellg();

    file.seekg(0, file.beg);
    const uint_fast64_t nb_values = file_length / sizeof(double);

    if ((nb_values % length) != 0)
    {
      ErrThrow("Error: in file: ", filename, " unexpected size to be read. "
               "Expected: ", length);
    }

    unsigned int number_snaps = nb_values / length;
    mat.resize(number_snaps);

    for (unsigned int i = 0; i < number_snaps; ++i)
    {
      mat[i].resize(length);
      file.read((char*)&(mat[i][0]), sizeof(double) * length);
    }
  }
  else
  {
    std::string line;

    while (getline(file, line))
    {
      std::vector<double> data;
      double value;
      std::istringstream iss(line);
      while (iss >> value)
      {
        data.push_back(value);
      }
      mat.push_back(data);
    }
  }

  file.close();
}
