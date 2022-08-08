/*
 * Unit testing of a few iterative solver classes.
 */
#include <BlockMatrix.h>
#include <Solver.h>
#include <MooNMD_Io.h>
#include <BlockVector.h>

int main(int, char**)
{
  Output::increaseVerbosity(1);
  /**
   * A = [ 5.5, 0,    0,  0  ; B = [ 4, 0, 0, 1; C = B^T
   *       0,   7.53, 0,  0  ;       0, 0, 5, 0;
   *       0,   0,    10, 0  ;       0, 6, 0, 1 ]
   *       0,   0,    0,  15 ]
   *
   * D = [ 50 0     0  ;
   *       0   19   0.1;
   *       0   0.1  85.6 ]
   */
  
  //Matrix A
  int * RowA = new int[5];
  int * ColA = new int[4];
  ColA[0] = 0; ColA[1] = 1; ColA[2] =  2; ColA[3] = 3;
  RowA[0] = 0; RowA[1] = 1; RowA[2] = 2; RowA[3] = 3; RowA[4] = 4;
  std::shared_ptr<TStructure> structureA(new TStructure(4, 4, ColA, RowA));
  auto matA = std::make_shared<TMatrix>(structureA);
  matA->GetEntries()[0] = 5.5;
  matA->GetEntries()[1] = 7.53;
  matA->GetEntries()[2] = 10;
  matA->GetEntries()[3] = 15;
  
  //Matrix B
  int * RowB = new int[4];
  int * ColB = new int[5];
  RowB[0] = 0; RowB[1] = 2; RowB[2] = 3; RowB[3] = 5;
  ColB[0] = 0; ColB[1] = 3; ColB[2] = 2; ColB[3] = 1; ColB[4] = 3;
  std::shared_ptr<TStructure> structureB(new TStructure(3, 4, 5, ColB, RowB));
  auto matB = std::make_shared<TMatrix>(structureB);
  matB->GetEntries()[0] = 4;
  matB->GetEntries()[1] = 1;
  matB->GetEntries()[2] = 5;
  matB->GetEntries()[3] = 6;
  matB->GetEntries()[4] = 1;
  
  //Matrix C
  auto matC = std::make_shared<TMatrix>(*matB->get_transposed());
  
  //Matrix D
  int * RowD = new int [4];
  int * ColD = new int [5];
  RowD[0] = 0; RowD[1] = 1; RowD[2] = 3; RowD[3] = 5;
  ColD[0] = 0; ColD[1] = 1; ColD[2] = 2; ColD[3] = 1; ColD[4] = 2;
  std::shared_ptr<TStructure> structureD(new TStructure(3, 3, 5, ColD, RowD));
  auto matD = std::make_shared<TMatrix>(structureD);
  matD->GetEntries()[0] = 50;
  matD->GetEntries()[1] = 19;
  matD->GetEntries()[2] = 0.1;
  matD->GetEntries()[3] = 0.1;
  matD->GetEntries()[4] = 85.6;

  // write the individual matrices
  //matA->PrintFull("matA", 4);
  //matB->PrintFull("matB", 4);
  //matC->PrintFull("matC", 4);
  //matD->PrintFull("matD", 4);
  
  // ##########################################################################
  // ##########################################################################
  // test the solving for this very small example
  
  //BlockMatrix mat(2, 2, {matA, matC, matB, matD});
  BlockMatrix mat({size_t(matA->get_n_rows()), size_t(matB->get_n_rows())},
                         {size_t(matA->get_n_columns()), 
                          size_t(matC->get_n_columns()) } );
  mat.replace_blocks(*matA, {{0,0}}, {false});
  mat.replace_blocks(*matC, {{0,1}}, {false});
  mat.replace_blocks(*matB, {{1,0}}, {false});
  mat.replace_blocks(*matD, {{1,1}}, {false});
  
  //BlockMatrix mat({size_t(matA->get_n_rows())}, {size_t(matA->get_n_columns())});
  //mat.replace_blocks(*matA, {{0,0}}, {false});
  
  BlockVector sol(mat, false);
  BlockVector exact_sol(sol);
  BlockVector rhs(mat, true);
  
  // fill exact solution with some nonzero numbers
  for(size_t i = 0; i < sol.length(); ++i)
    exact_sol.at(i) = double(i);
  // compute correct rhs = mat*exact_sol
  mat.apply(exact_sol, rhs);
  
  // write entire matrix, solution and right hand side
  mat.get_combined_matrix()->PrintFull("M", 4);
  exact_sol.print("exact_sol");
  rhs.print("rhs");
  
  // create a database which is modified a little bit for the different solvers
  Solver<BlockMatrix, BlockVector> dummy(ParameterDatabase("dummy database"));
  // copy the (default) database, `dummy` no longer needed then
  ParameterDatabase db(dummy.get_db());
  db["solver_type"] = "iterative"; // set to iterative solvers
  db["residual_tolerance"] = 1.0e-11;
  db["max_n_iterations"] = 1000;
  // cg requires a symmetric preconditioner
  db["preconditioner"] = "ssor"; // "jacobi" also ok
  db["sor_omega"] = 1.25;
  
  // lambda function to do the solving and to check the error
  auto solve_and_check
    = [&mat, &rhs, &sol, &exact_sol]
      (Solver<BlockMatrix, BlockVector>& solver) -> bool
      {
        sol.reset(); // start with same solution for every method
        solver.solve(mat, rhs, sol);
        // substract exact solution sol =  A^{-1}*rhs - exact_sol
        sol.add_scaled(exact_sol, -1.);
        // this is non-negative and should be small
        const double error = sol.norm();
        double desired_error = solver.get_db()["residual_tolerance"];
        // the following number is larger because we check the error, not the 
        // residual (which is A multiplied by the error).
        desired_error*= 100.;
        // check if error is not too large, the second condition is true for 
        // nan (not a number)
        if(error > desired_error || error != error)
        {
          Output::print("norm of error is too big ", error);
          return false;
        }
        return true;
      };
  
  // a map to store all solvers to test here
  std::map<std::string, std::string> solvers;
  solvers["Jacobi"] = "jacobi";
  solvers["SOR"] = "sor";
  solvers["SSOR"] = "ssor";
  solvers["Richardson"] = "richardson";
  solvers["cg"] = "cg";
  solvers["cgs"] = "cgs";
  solvers["bi-cgstab"] = "bi_cgstab";
  solvers["left gmres"] = "left_gmres";
  solvers["right gmres"] = "right_gmres";
  solvers["flexible gmres"] = "fgmres";
  
  for(auto & solver_info : solvers) // s is a std::pair
  {
    Output::print<1>("\ntesting ", solver_info.first, " iteration");
    db["iterative_solver_type"] = solver_info.second;
    
    Solver<BlockMatrix, BlockVector> s(db);
    if(!solve_and_check(s))
    {
      Output::print<1>(solver_info.first, " iteration failed");
      return 1;
    }
  }
  
  Output::print("\nTest program finished.");

}


