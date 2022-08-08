#include <LoopInfo.h>
#include <MooNMD_Io.h>
#include "Residuals.h"
#include <limits>
#include <cmath>

#ifdef _MPI
#include <mpi.h>
#endif

template <typename T>
LoopInfo<T>::LoopInfo(const std::string& name)
 : name(name),
   initial_residual(std::is_arithmetic<T>::value ?
                    T{std::numeric_limits<T>::max()} : T{}),
   old_residual(std::is_arithmetic<T>::value ?
                T{std::numeric_limits<T>::max()} : T{}),
   timer(), n_previous_iterations(0)
{

}

template <typename T>
LoopInfo<T>::LoopInfo(const char* name)
 : LoopInfo(std::string(name))
{
  
}

template <typename T>
LoopInfo<T>::LoopInfo(const std::string& name,
                   bool print_time_every_step_in,
                   bool print_reduction_rates_in,
                   size_t verbosity_threshold_in)
  : LoopInfo(name)
{
  print_time_every_step = print_time_every_step_in;
  print_reduction_rates = print_reduction_rates_in;
  verbosity_threshold   = verbosity_threshold_in;
}

template <typename T>
void LoopInfo<T>::restart(const std::string& name, T initial_residual)
{
  this->name = name;
  this->initial_residual = initial_residual;
  this->old_residual = std::is_arithmetic<T>::value 
                       ? T{std::numeric_limits<T>::max()} : T{};
  this->timer.reset();
}

template <typename T>
void LoopInfo<T>::print(unsigned int loop_index, T current_residual)
{
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#else
  int my_rank = 0;
#endif

  // print only if verbosity is high enough.
  if(this->verbosity_threshold <= Output::getVerbosity() && my_rank == 0)
  {
    using namespace std;
    std::stringstream s;

    s << this->name << " iteration: " << setw(3) << loop_index << "  ";
    if(std::is_arithmetic<T>::value)
    {
      s << "residual: " << left << setprecision(10) << setw(15);
    }
    s << current_residual;
    if(loop_index > 0)
    {
      if(this->print_reduction_rates)
      {
        s << "  red/step: " << setprecision(10) << setw(15) 
          << double(current_residual)/double(old_residual)
          << "  red: " << setprecision(10) << setw(15)
          << double(current_residual)/double(initial_residual);
      }
      if(this->print_time_every_step)
      {
        s << "  t[s]/step: " << setprecision(4) << setw(9)
          << timer.time_since_last_start()
          << "  t[s]: " << setprecision(4) << setw(9) 
          << timer.elapsed_time();
      }
    }
    else
    {
      // should we reset the initial_time here or rather not?
      this->restart(this->name, current_residual);
    }
    Output::print<3>(s.str());
  }
  old_residual = current_residual;
  timer.stop();
  timer.start();
}

template <typename T>
void LoopInfo<T>::finish(unsigned int loop_index, T current_residual)
{
  this->n_previous_iterations += loop_index;
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#else
  int my_rank = 0;
#endif
  
  timer.stop();
  using namespace std;
  std::stringstream s;
  
  s << this->name << " iteration finished: " << setw(3) << loop_index
    << " residual: " << left << setprecision(10) << setw(15)
    << current_residual;
  if(loop_index > 0)
  {
    s << "  red: " << setprecision(10) << setw(15) 
      << current_residual/this->initial_residual 
      << "  red/step: " << setprecision(10) << setw(15)
      << std::pow(current_residual/this->initial_residual, 1./loop_index);
    s << "  t[s]: " << setprecision(4) << setw(9) 
      << timer.elapsed_time()
      << "  t[s]/step: " << setprecision(4) << setw(9)
      << this->timer.elapsed_time()/loop_index;
  }
  if(loop_index != this->n_previous_iterations)
  {
    s << "  total iterations: " << setw(3) <<this-> n_previous_iterations;
  }
  if(my_rank == 0)
    Output::print<2>(s.str());
}

template <typename T>
T LoopInfo<T>::get_initial_residual() const
{
  return this->initial_residual;
}

template <typename T>
T LoopInfo<T>::get_previous_residual() const
{
  return this->old_residual;
}

template <typename T>
unsigned int LoopInfo<T>::get_n_previous_iterations() const
{
  return n_previous_iterations;
}

template class LoopInfo<double>;
template class LoopInfo<Residuals>;
