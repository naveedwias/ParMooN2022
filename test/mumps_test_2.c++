/**
 * @brief A test program for actuating the MUMPS solver from ParMooN.
 *
 * This test program will test, whether the Mumps libraries we have included
 * can be used at all to solve a small system, distributed on four processes.
 *
 * Here we test the same solution strategy which ParMooN will employ when
 * making use MUMPS: the matrix is distributed among the processes (ICNTL(18)=3).
 *
 * Solves the system
 *   (0 0 1 1) (w)    (7)
 *   (0 0 1 1) (x)    (7)
 *   (2 2 3 3) (y)  = (27)
 *   (2 2 3 3) (z)    (27),
 * where A is distributed among four processes (each one holds the entries which
 * equal its own rank number). The solution is (1 2 3 4)^t.
 *
 * This has nothing to do yet with interfacing the ParMooN parallel structure
 * and distributed BlockFEMatrices with the MUMPS solver.
 *
 * @date 2016/03/09
 * @author Clemens Bartsch
 */
#ifdef _MPI
#include "all_defines_external_libraries.h"
#ifdef PARMOON_WITH_MUMPS
#include "dmumps_c.h"
#include <mpi.h>
#include <array>

//macros for mumps integer paramters - makes the code better readable
#define JOB_INIT -1
#define JOB_ANALYZE_FACTORIZE_SOLVE 6
#define JOB_END -2

int main(int argc, char* argv[])
{
  //set up and init a World MPI communicator
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);

  // Hold mpi rank and size ready
  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // build the mumps solver and set comm, par and sym
  DMUMPS_STRUC_C id; //a MUMPS structure
  id.comm_fortran = (MUMPS_INT) MPI_Comm_c2f(comm);
  id.par = 1; //host takes part in computations
  id.sym = 0; //unsymmetric matrix

  // let id initialize itself
  id.job = JOB_INIT;
  dmumps_c(&id);

  // macro to do the c-to-fortran index shift for the control parameters
  #define ICNTL(I) icntl[(I)-1]
  // set parameters which determine the solving process
  id.ICNTL(5) = 0;  // matrix in "assembled" format
  id.ICNTL(18) = 3; // matrix distributed among processes
  id.ICNTL(20) = 0; // rhs in dense format (default value)
  id.ICNTL(21) = 0; // solution centralized on host (default value)
  // explicitly put output control to default values (chose better?)
  id.ICNTL(1)= 6; // error message output stream
  id.ICNTL(2)= 0; // diagnostic message otuput stream
  id.ICNTL(3)= 0; // host output stream
  id.ICNTL(4)= 2; //verbosity level

  int nz_loc;
  std::array<int, 4> irn_loc;
  std::array<int, 4> jcn_loc;
  std::array<double, 4> A_loc;
  std::array<double, 4> rhs; // starts its life on host as right hand side,
                             // comes out as solution on host, workers just need the memory

  // each process fills the arrays with its local information
  switch (mpi_rank)
  {
    case 0: //this is the host
      nz_loc = 0;
      rhs = {1,7,23,27};
      break;
    case 1:
      nz_loc = 4;
      irn_loc = {1,1,2,2};
      jcn_loc = {3,4,3,4};
      A_loc = {-1,1,1,1};
      break;
    case 2:
      nz_loc = 4;
      irn_loc = {3,3,4,4};
      jcn_loc = {1,2,1,2};
      A_loc = {-2,2,2,2};
      break;
    case 3:
      nz_loc = 4;
      irn_loc = {3,3,4,4};
      jcn_loc = {3,4,3,4};
      A_loc = {3,3,3,3};
      break;
    default:
      throw std::runtime_error("That's an unexpected process number!");
  }
  //fill the generated data into the mumps structure id
  id.nz_loc = nz_loc;
  id.irn_loc = &irn_loc.front();
  id.jcn_loc = &jcn_loc.front();
  id.a_loc = &A_loc.front();

  //data filling which is only done by host
  if (mpi_rank ==0)
  {
    //global matrix size
    id.n = 4;
    //centralized rhs
    id.rhs = &rhs.front();
    id.nrhs = 1;
    id.lrhs = 4;
  }

  //do work
  id.job = JOB_ANALYZE_FACTORIZE_SOLVE;
  dmumps_c(&id);

  //check work
  if (mpi_rank == 0) {
    printf("Solution is : (%8.2f  %8.2f %8.2f %8.2f)\n", rhs[0],rhs[1],rhs[2],rhs[3]);
    for (int i = 1; i<5 ;++i)
    {
      if( abs(rhs[i-1] - i) > 1e-10)
      {
        throw std::runtime_error("That's the wrong result!");
      }
    }
  }

  //clean up
  id.job=JOB_END;
  dmumps_c(&id);

  // finalize MPI
  MPI_Finalize();

}

#else
int main(){}
#endif // PARMOON_WITH_MUMPS
#endif // _MPI
