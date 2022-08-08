/**
 * @brief A test program for actuating the MUMPS solver from ParMooN.
 *
 * This test program will test, whether the Mumps libraries we have included
 * can be used at all to solve a tiny system, distributed on two processes.
 * It is our re-implementation of the C-interface example program which ships
 * with the MUMPS package.
 * Solves the system
 *   (1 0) (x)    (1)
 *   (0 2) (y)  = (4)
 *
 * in ... way.
 * Solution is (1 2)^t.
 *
 * It is only done in MPI compile and run mode so far.
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
#include <cmath>
#include <stdexcept>

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

  //arrays needed to describe matrix and rhs
  int irn[2] = {1,2};
  int jcn[2] = {1,2};
  double A[2] = {1.0, 2.0};
  double rhs[2] = {1.0, 4.0};

  // host process communicates them to its mumps object
  if (mpi_rank == 0) {
    id.n = 2;
    id.nz = 2;
    id.irn = irn;
    id.jcn = jcn;
    id.a = A;   // matrix A entries
    id.rhs = rhs; // rhs vector
  }

  // macro to do the c-to-fortran index shift for the control parameters
  #define ICNTL(I) icntl[(I)-1]
  // explicitly put output control to default values
  id.ICNTL(1)= 6; // error message output stream
  id.ICNTL(2)= 0; // diagnostic message otuput stream
  id.ICNTL(3)= 0; // host output stream
  id.ICNTL(4)= 2; //verbosity level

  // do the work
  id.job = JOB_ANALYZE_FACTORIZE_SOLVE; //6 means: do 1 (analyze), 2 (factorize) and 3 (solve) one after the other
  dmumps_c(&id);

  // terminate the work
  id.job=JOB_END;
  dmumps_c(&id);

  if (mpi_rank == 0) {
    printf("Solution is : (%8.2f  %8.2f)\n", rhs[0],rhs[1]);
    //check the result
    if (std::abs(rhs[0] - 1.0) < 1e-10 && std::abs(rhs[1] - 2.0) < 1e-10)
    {
      printf("I mumpsed it!\n");
    }
    else
    {
      throw std::runtime_error("That's the wrong result!");
    }
  }

  // finalize the mpi communicator
  MPI_Finalize();

}

#else // PARMOON_WITH_MUMPS
int main(){}
#endif // PARMOON_WITH_MUMPS
#endif // _MPI
