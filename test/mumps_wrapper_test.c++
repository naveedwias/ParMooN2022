/**
 * Unit test for the MUMPS solver wrapper.
 *
 * The big problem of this test is, that we have to run a good chunk of
 * ParMooN code in order to generate the objects which we need to initialize
 * the actual Mumps Solver. The parallelization is badly modularized.
 * A particular problem is that the mesh partitioning (achieved with Metis)
 * might produce different partitioning on different systems, and thus the
 * hard-coded values might be useless for testing.
 *
 * @date 2016/03/14
 * @author Clemens Bartsch
 */
#ifdef _MPI
#include "all_defines_external_libraries.h"
#ifdef PARMOON_WITH_MUMPS

#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <Database.h>
#include <FESpace3D.h>
#include <Domain.h>
#include "MumpsWrapper.h"
#include <MeshPartition.h>
#include <MainUtilities.h>
#include <ParFECommunicator3D.h>
#include <ParameterDatabase.h>
#include "ParMooN.h"

#include <mpi.h>
#include <string>

int main(int argc, char* argv[])
{
  ParameterDatabase db = parmoon::parmoon_initialize(argc, argv);
  // generate the needed BlockFEMatrices by running a parallel program
  // up to the point where the whole par infrastructure is available
  //  - at the moment, modularizing the process is out of reach

  MPI_Comm comm = MPI_COMM_WORLD; //mpi setup
  TDatabase::ParamDB->Comm = comm;
  int mpiRank, mpiSize;
  MPI_Comm_rank(comm, &mpiRank);
  MPI_Comm_size(comm, &mpiSize);


  db.add("boundary_file", "Default_UnitCube", "");
  db.add("geo_file", "Default_UnitCube_Hexa", "", 
         {"Default_UnitCube_Hexa", "Default_UnitCube_Tetra"});
  // default construct a domain object
  TDomain domain(db);


  // the domain is initialised with default description and default
  // initial mesh
  //domain.Init((char*)"Default_UnitCube", (char*)"Default_UnitCube_Hexa");

  // refine grid just once
  domain.RegRefineAll();

  //partition grid and determine maxSubDomainPerDof
  domain.GenerateEdgeInfo();
  int maxCellsPerVertex;
  Partition_Mesh3D(comm, &domain, 0, maxCellsPerVertex);
  domain.GenerateEdgeInfo();

  Output::print(" Rank ", mpiRank, ": domain.GetN_OwnCells() = ", domain.GetN_OwnCells());
  //  // this is hardcoded for the 2 processes case - we expect rank 0
  //  // to have 5 own cells and rank 1 to have 3.
  //  if (mpiRank == 0 && domain.GetN_OwnCells() != 5)
  //    ErrThrow("The domain partitioning did not give 5 own cells to process 0.");
  //  if (mpiRank == 1 && domain.GetN_OwnCells() != 3)
  //    ErrThrow("The domain partitioning did not give 3 own cells to process 1.");

  //collection from finest cell level
  TCollection *coll = domain.GetCollection(It_EQ, 0);

  // Create two FESpace3D to fiddle around with
  size_t first_ansatz_order = 1;
  size_t second_ansatz_order = 2;

  std::shared_ptr<const TFESpace3D> fe_space_1(
    new TFESpace3D(coll, "first_fe_space", BoundCondition_FEM_FCT,
                   first_ansatz_order)); //with dirichlet boundaries

  std::shared_ptr<const TFESpace3D> fe_space_2(
    new TFESpace3D(coll, "second_fe_space", BoundaryConditionNewton,
                   second_ansatz_order)); //actives only
  //end blackbox

  //build ParFECommunicators, -mappers and a BlockFEMatrix for testing
  TParFEMapper3D mapper_1(1, fe_space_1.get());
  TParFECommunicator3D comm_1(&mapper_1);
  TParFEMapper3D mapper_2(1, fe_space_2.get());
  TParFECommunicator3D comm_2(&mapper_2);

  BlockFEMatrix bfem({fe_space_1, fe_space_2});

  FEMatrix fe_matrix_1(fe_space_1);
  FEMatrix fe_matrix_2(fe_space_2, fe_space_1);
  FEMatrix fe_matrix_2_t(fe_space_1, fe_space_2);
  FEMatrix fe_matrix_3(fe_space_2);
  fe_matrix_1.setEntries(std::vector<double>(fe_matrix_1.get_n_entries(),1.0));
  fe_matrix_2.setEntries(std::vector<double>(fe_matrix_2.get_n_entries(),2.0));
  fe_matrix_2_t.setEntries(std::vector<double>(fe_matrix_2_t.get_n_entries(),2.0));
  fe_matrix_3.setEntries(std::vector<double>(fe_matrix_3.get_n_entries(),3.0));

  bfem.replace_blocks(fe_matrix_1, {{0,0}}, {false});
  //bfem.replace_blocks(fe_matrix_2, {{0,1},{1,0}}, {false, true});
  bfem.replace_blocks(fe_matrix_2_t, {{0,1}},{false});
  bfem.replace_blocks(fe_matrix_2, {{1,0}}, {false});
  bfem.replace_blocks(fe_matrix_3, {{1,1}}, {false});
  bfem.print_and_check("test block fe matrix");

  // actually the combined matrix should be the same (resp. permutations)
  // on both processes, for all dofs are known to both processes here.
  bfem.get_combined_matrix()->write((std::string("1s") +
      std::to_string(mpiRank)).c_str());

  // first goal: setting up the master-only matrices on both processes
  // and re-combining (printout, MATLAB) should lead to a matrix of the same
  // norms as the combined matrix has

  {//start new scope, so that Mumps solver is destructed before MPI_FINALIZE
    //set up a mumps wrapper for the block matrix
    MumpsWrapper wrap(bfem);

    // //write matrices to file
    //wrap.write_matrix_distributed(std::string("1d"));
  }

  {
    //for fiddling around with: nse-style matrix with one off-diagonal zero block
    BlockFEMatrix nsemat({fe_space_1, fe_space_1, fe_space_2});
    nsemat.replace_blocks(fe_matrix_1,
                          {{0,0},{0,1},{1,1}},
                          {false, false, false});
    nsemat.replace_blocks(fe_matrix_2_t,
                          {{0,2},{1,2},{2,0},{2,1}}, {false, false, true, true});
    nsemat.replace_blocks(fe_matrix_3,
                          {{2,2}}, {false});

    for(int i = 0; i < mpiSize; i++) {
      MPI_Barrier(MPI_COMM_WORLD);
      if (i == mpiRank) {
        nsemat.print_and_check("nsemat");
      }
    }
    nsemat.get_combined_matrix()->write((std::string("2s") +
        std::to_string(mpiRank)).c_str());
    MumpsWrapper wrap(nsemat);
    wrap.write_matrix_distributed(std::string("2d"));

    // Here I used MATLAB to read in the matrices - those gained by get_combined_matrix
    // as well as those gained by write_matrix_distributed and compared their Frobenius
    // norms. Since for the examples all dofs are known to all processes, the matrices
    // gained with get_combined_matrix are globally complete on all processes. Thus,
    // if the matrix gained by adding the matrices from write_matrix_distributed
    // (in MATLAB) matches its Frobenius norm, that is a strong hint for having
    // produced the correct matrices by the mumps solver wrapper.

    //now create BlockVectors to try the solver with
    BlockVector rhs(nsemat ,true);
    BlockVector sol(nsemat, false);
    for(unsigned int i =0; i<rhs.length(); ++i)
    {//fill rhs with ones
      rhs.at(i) = 1;
    }
    wrap.solve(rhs, sol );

    //hard-coded test, I don't know anything better at the moment...
    Output::print(sol.norm(nsemat.get_communicators()));
    if(std::abs(sol.norm(nsemat.get_communicators()) - 24.6171) > 1e-4)
    {// solution is same (resp. permutations)
     // on all processes, cause all dofs are known to both processes,
     // which makes this easy test possible. For bigger matrix- check residual
     // of the local system!
        ErrThrow("Wrong solution norm!");
    }

  }
  // say goodbye to poppa
  parmoon::parmoon_finalize();
}

#else // PARMOON_WITH_MUMPS
int main(){}
#endif // PARMOON_WITH_MUMPS
#endif // _MPI


