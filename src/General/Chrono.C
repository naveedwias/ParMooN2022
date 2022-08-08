/**
 * @file Chrono.C
 * Implementation of class Chrono declared in Chrono.h.
 *
 * @date 2016/04/18
 */
#include <Chrono.h>
#include <MooNMD_Io.h>
#include <time.h>
#include <array>

#ifdef _MPI
#include <mpi.h>
#endif

#ifdef _OMP
#include <omp.h>
#endif

/**
 * Evaluates system time spent in user mode and in kernel mode so far.
 * Returns sum of both.
 *
 * @return added up time spent in user and kernel mode since start of the 
 * program
 */
double get_rusage_time()
{
  rusage usage;

  if(getrusage(RUSAGE_SELF, &usage) == -1)
    ErrThrow("Error in GetTime!");

  // user mode time and system mode time in s
  double user_mode_time = usage.ru_utime.tv_sec;
  double kernel_mode_time =  usage.ru_stime.tv_sec;

  //add microseconds
  user_mode_time += ((double) usage.ru_utime.tv_usec)/1000000;
  kernel_mode_time += ((double) usage.ru_stime.tv_usec)/1000000;

  //return sum of both times
  return user_mode_time + kernel_mode_time;
}

/// @return Current wall time.
double get_wall_time()
{
#ifdef _OMP
  return omp_get_wtime();
#endif
#ifdef _MPI
  return MPI_Wtime();
#endif
  //wall time
  struct timeval wall_time;
  if (gettimeofday(&wall_time,nullptr))
    ErrThrow("Error in gettimeofday!");
  return wall_time.tv_sec + wall_time.tv_usec/(double)1000000;
}


Chrono::Chrono()
{
  reset();
  start();
}


void Chrono::reset()
{
  this->cumulative_time_rusage = 0;
  this->cumulative_time_wall = 0;
  this->running = false; // to avoid a warning in the start() method
}

void Chrono::start()
{
  if(running)
  {
    Output::warn<1>("Chrono class", "starting a Chrono object which is "
                    "already running. Did you want to do a reset()?");
    this->stop();
  }
  this->start_time_rusage = get_rusage_time();
  this->start_time_wall = get_wall_time();
  this->running = true;
}


double Chrono::stop()
{
  if(running)
  {
    this->cumulative_time_rusage += this->time_since_last_start();
    double t = this->wall_time_since_last_start(); // will be returned
    this->cumulative_time_wall += t;
    this->running = false;
    return t;
  }
  else
    return 0.0;
}

// a method for convenience, called from both print_total_time and
// print_time_since_last_start
void print_times(std::array<double, 2> times, const std::string& program_part)
{
#ifndef _MPI
  Output::print("--- TIME: time for ", program_part,": ", times[0], ", ", 
                times[1], " s (CPU time, wall time)");
#else
  // MPI
  MPI_Comm comm = MPI_COMM_WORLD;
  int my_rank;
  int size;
  MPI_Comm_rank(comm, &my_rank);
  MPI_Comm_size(comm, &size);

  //rusage time in root
  std::array<double, 2> t_min;
  std::array<double, 2> t_max;
  std::array<double, 2> t_avg;
  MPI_Reduce(&times[0], &t_max[0], 2, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&times[0], &t_min[0], 2, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&times[0], &t_avg[0], 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  t_avg[0] /= size;
  t_avg[1] /= size;
  if(my_rank == 0)
  {
    Output::print("--- TIME: time for ", program_part,": ", t_max[0], ", ",
                  t_min[0], ", ", t_avg[0],
                  " s (max, min, average CPU time over all processes)");
    Output::print("--- TIME: time for ", program_part,": ", t_max[1], ", ",
                  t_min[1], ", ", t_avg[1],
                  " s (max, min, average wall time over all processes)");
  }
#endif
}

void Chrono::print_total_time(const std::string& program_part) const
{
  std::array<double, 2> times = {{this->elapsed_time(),        // CPU time
                                  this->elapsed_wall_time()}}; // wall time
  print_times(times, program_part);
}

void Chrono::restart_and_print(const std::string& program_part)
{
  stop_and_print(program_part);
  start();
}

void Chrono::stop_and_print(const std::string& program_part)
{
  std::array<double, 2> times = {{time_since_last_start(),        // CPU
                                  wall_time_since_last_start()}}; // wall
  stop();
  print_times(times, program_part);
}

double Chrono::time_since_last_start() const
{
  if(this->running)
    return get_rusage_time() - this->start_time_rusage;
  else
    return 0.0;
}

double Chrono::wall_time_since_last_start() const
{
  if(this->running)
    return get_wall_time() - start_time_wall;
  else
    return 0.0;
}


double Chrono::elapsed_time() const
{
  return this->cumulative_time_rusage + this->time_since_last_start();
}

double Chrono::elapsed_wall_time() const
{
  return this->cumulative_time_wall + this->wall_time_since_last_start();
}
