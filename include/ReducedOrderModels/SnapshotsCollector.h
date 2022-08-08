/** ****************************************************************************
*
* @name       SnapshotsCollector
* @brief      write snapshots (finite elements coefficients of the solution)
*             into a file
*
* @author     Swetlana Giere & Alfonso Caiazzo
* @date       08.03.2017. Restarted on 15.1.2019
*
*******************************************************************************/

#ifndef SNAPSHOTSCOLLECTOR_H
#define SNAPSHOTSCOLLECTOR_H

#include <MooNMD_Io.h>
#include <ParameterDatabase.h>
#include <BlockVector.h>

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>


class SnapshotsCollector
{
  public:
    /**
    * @brief constructor
    *
    * The filename for storing snapshots is constructed as
    * filename = db["snaps_directory"] + "/" + db["snaps_basename"].
    * Then, the stream is created and stored as a member of the class.
    */
    SnapshotsCollector(const ParameterDatabase& param_db,
                       std::string filename_suffix="");

    /**
    * @brief destructor closes the datafile
    */
    ~SnapshotsCollector();

    /**
     * Creates a database filled with default parameters. This database will
     * contain all necessary parameters for the SnapshotsCollector.
     */
    static ParameterDatabase default_snapshots_database();

    static std::string get_snapshot_filename(const ParameterDatabase& db, std::string filename_suffix);

    /**
    * @brief Write snapshots into file
    *
    * Write finite elemenent coefficients of the solution row-wise into the file
    * represented by the class member 'datastream'. Moreover, the snapshots
    * will be written only if the database parameter db["write_snaps"] is set
    * to true and the condition time_step % db["steps_per_snap"] == 0 is
    * satisfied.
    *
    * Note: For 3D problems a binary format could be more appropriate to avoid
    * huge storage volume requirements, for this purpose set the parameter
    * 'write_solution_binary' to true
    *
    * @param[in] solution  pointer to the finite elemenent coefficients
    *                      of the solution (snapshot)
    * @param[in] length    number of values in solution
    * @param[in] time_step the time step from the main program time loop
    *
    * @return flag to signal if the conditions to be written were satisfied
    */
    bool write_data(const double* const solution,
                    const unsigned int length,
                    size_t time_step=0);

    /**
    * @brief Write snapshots containing solution and time derivative into file
    *
    * This function call the function write_data(const double* const, const int,
    * size_t) and use its return value as written condition. The time derivative
    * will be written if the parameter db["snaps_time_derivative"] is set to
    * true, but not for the first time step
    *
    * @param[in] solution  pointer to the finite elemenent coefficients
    *                      of the solution (snapshot)
    * @param[in] length    number of values in solution
    * @param[in] time_step the time step from the main program time loop
    * @param[in] sol_m1    pointer to the finite elemenent coefficients
    *                      of the solution at the last time step
    * @param[in] tau       the current time step length
    */
    void write_data(const double* const solution,
                    const unsigned int length,
                    size_t time_step,
                    const double* const sol_m1,
                    double tau);

    /**
    * @brief Write snapshots containing solution and time derivative into file
    *
    * This function call the function write_data(const double* const, const int,
    * size_t) and use its return value as written condition. The time derivative
    * will be written if the parameter db["snaps_time_derivative"] is set to
    * true, but not for the first time step
    *
    * @param[in] solution     pointer to the finite elemenent coefficients
    *                         of the solution (snapshot)
    * @param[in] dt_solution  pointer to the finite elemenent coefficients
    *                         of the time derivative (snapshot)
    * @param[in] length       number of values in solution
    * @param[in] time_step    the time step from the main program time loop
    */
    void write_data(const double* const solution,
                    const unsigned int length,
                    size_t time_step,
                    const double* const dt_solution);

    /**
    * @brief Read snapshots from file
    *
    * Read data from the file filename and save it into the matrix mat
    *
    * @param[in]  filename full name of the file
    * @param[in]  bin      flag to indicate if file filename is a binary file
    * @param[in]  length   length of a snapshot
    * @param[out] mat      matrix
    */
    static void read_data(std::string                        filename,
                          bool                               bin,
                          int                                length,
                          std::vector< std::vector<double>>& mat);

  private:
    /* stream for writing snapshots */
    std::ofstream datastream;

    ParameterDatabase db;

    int snap_count;

    std::string snapshot_filename;
};

#endif // SNAPSHOTSCOLLECTOR_H
