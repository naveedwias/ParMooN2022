#include <DirectSolver.h>
#include <BlockMatrix.h>
#include <BlockVector.h>
#include <MooNMD_Io.h>
#include <cmath>

bool equal(const double a, const double b)
{
  if(std::abs(a) > 1e-12)
    return std::abs((a-b)/a) < 1e-14;
  else
    return std::abs(a-b) < 1e-14;
}

int main(int, char**) 
{
  { // create local scope
  Output::increaseVerbosity(5);
  /**
   * A = [ 0.5, 0,    0, 0  ; B = [ 4, 0, 0, 1; C = [ 7, 0, 8;
   *       0,   0.75, 0, 0  ;       0, 0, 5, 0;       9,10, 0;
   *       0,   0,    1, 0  ;       0, 6, 0, 1 ]      0,11,12;
   *       0,   0,    0, 1.5 ]                        13,0,14 ]
   *
   * D = [ 0.5 0   0  ;
   *       0   1   0.1;
   *       0   0.1 1.5 ]
   */
  
  //Matrix A
  int * RowA = new int[5];
  int * ColA = new int[4];
  ColA[0] = 0; ColA[1] = 1; ColA[2] =  2; ColA[3] = 3;
  RowA[0] = 0; RowA[1] = 1; RowA[2] = 2; RowA[3] = 3; RowA[4] = 4;
  std::shared_ptr<TStructure> structureA(new TStructure(4, 4, ColA, RowA));
  auto matA = std::make_shared<TMatrix>(structureA);
  matA->GetEntries()[0] = 0.5;
  matA->GetEntries()[1] = 0.75;
  matA->GetEntries()[2] = 1;
  matA->GetEntries()[3] = 1.5;

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
  int * RowC = new int[5];
  int * ColC = new int[8];
  RowC[0] = 0; RowC[1] = 2; RowC[2] = 4; RowC[3] = 6; RowC[4] = 8;
  ColC[0] = 0; ColC[1] = 2; ColC[2] = 0; ColC[3] = 1; ColC[4] = 1; ColC[5] = 2;
  ColC[6] = 0; ColC[7] = 2;
  std::shared_ptr<TStructure> structureC(new TStructure(4, 3, 8, ColC, RowC));
  auto matC = std::make_shared<TMatrix>(structureC);
  matC->GetEntries()[0] = 7;
  matC->GetEntries()[1] = 8;
  matC->GetEntries()[2] = 9;
  matC->GetEntries()[3] = 10;
  matC->GetEntries()[4] = 11;
  matC->GetEntries()[5] = 12;
  matC->GetEntries()[6] = 13;
  matC->GetEntries()[7] = 14;

  
  //Matrix D
  int * RowD = new int [4];
  int * ColD = new int [5];
  RowD[0] = 0; RowD[1] = 1; RowD[2] = 3; RowD[3] = 5;
  ColD[0] = 0; ColD[1] = 1; ColD[2] = 2; ColD[3] = 1; ColD[4] = 2;
  std::shared_ptr<TStructure> structureD(new TStructure(3, 3, 5, ColD, RowD));
  auto matD = std::make_shared<TMatrix>(structureD);
  matD->GetEntries()[0] = 0.5;
  matD->GetEntries()[1] = 1;
  matD->GetEntries()[2] = 0.1;
  matD->GetEntries()[3] = 0.1;
  matD->GetEntries()[4] = 1.5;

  matA->PrintFull("matA", 4);
  matB->PrintFull("matB", 4);
  matC->PrintFull("matC", 4);
  matD->PrintFull("matD", 4);
  
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
  BlockVector sol(mat, false);
  BlockVector exact_sol(sol);
  BlockVector rhs(mat, true);
  //sol.info();
  //rhs.info();
  
  // fill exact solution with some nonzero numbers
  for(size_t i = 0; i < sol.length(); ++i)
    exact_sol.at(i) = double(i);
  // compute correct rhs = mat*exact_sol
  mat.apply(exact_sol, rhs);
  
  // compute solution
  auto type = DirectSolver::DirectSolverTypes::umfpack;
  class DirectSolver ds(mat, type);
  ds.solve(rhs, sol); // sol = A^{-1}*rhs   which  should be exact_sol
  
  //rhs.print("rhs");
  sol.print("solution");
  
  // substract exact solution sol =  A^{-1}*rhs - exact_sol
  sol.add_scaled(exact_sol, -1.);
  
  const double error = sol.norm(); // this is non-negative and should be zero
  if(error > 1.1e-14)
  {
    Output::print("norm of error is too big ", error);
    return 1;
  }
  
  // ##########################################################################
  // ##########################################################################
  // test the class DirectSolver a little bit
  //class DirectSolver ds2(ds); // should not compile
  //class DirectSolver ds2 = ds; // should not compile
  // move constructor
  class DirectSolver ds2(std::move(ds));
  // move assignment
  ds2 = DirectSolver(mat, type);
  } // end local scope
  
  Output::print("test successful");
  return 0;
}
