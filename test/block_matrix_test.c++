/*
 * Unit testing of the block matrix.
 *
 *  Created on: Dec 10, 2015
 *      Author: bartsch
 */


#include <vector>
#include <tuple>


#include <BlockMatrix.h>
#include <BlockVector.h>

#include <MooNMD_Io.h>

int main(int, char**)
{
  {
  //check correct splitting and merging of colors when replacing blocks
  BlockMatrix myMatrix({5 , 5 , 5}, {5, 5 , 5});

  TMatrix placeHolderA(5,5);
  auto t_1 = std::make_tuple(0,0,false);
  auto t_2 = std::make_tuple(0,1,false);
  auto t_3 = std::make_tuple(1,0,false);
  auto t_4 = std::make_tuple(1,1,false);
  std::vector<std::tuple<size_t, size_t, bool>> vecA = {t_1, t_2, t_3, t_4};

  TMatrix placeHolderB(5,5);
  std::vector<std::tuple<size_t, size_t, bool>> vecB = {t_1, t_4};


  myMatrix.replace_blocks(placeHolderA,
                          {{0,0}, {0,1}, {1,0}, {1,1} },
                          {false,false,false,false});
  myMatrix.check_coloring();

  myMatrix.replace_blocks(placeHolderB,
                          { {0,0}, {1,1} },
                          { false, false } );
  myMatrix.check_coloring();

  myMatrix.replace_blocks(placeHolderA,
                          {{0,0}, {0,1}, {1,0}, {1,1} },
                          {false,false,false,false});
  myMatrix.check_coloring();

  }

  {
  //check if it works with tranposition
  BlockMatrix myBlockMatrix({3,7},{3,7});
  myBlockMatrix.check_coloring();

  TMatrix placeHolderA(7,3);

  myBlockMatrix.replace_blocks(placeHolderA,
                               {{0,1}, {1,0}},
                               {true, false}  );
  myBlockMatrix.print_and_check("insert of transp");
  myBlockMatrix.print_coloring_pattern("insert of transp", true);

  }

  //check blockwise addition
  {
    /**
     * A = [ 1,  0,  0,  0  ;
     *       0,  1,  0,  0  ;
     *       0,  0,  1,  0  ;
     *       2,  0,  0,  1 ]
     */

    //Matrix A
    int * RowA = new int[5];
    int * ColA = new int[5];
    ColA[0] = 0; ColA[1] = 1; ColA[2] =  2; ColA[3] = 0; ColA[4] = 3;
    RowA[0] = 0; RowA[1] = 1; RowA[2] = 2; RowA[3] = 3; RowA[4] = 5;
    std::shared_ptr<TStructure> structureA(new TStructure(4, 5, ColA, RowA));
    TMatrix matA(structureA);
    matA.GetEntries()[0] = 1;
    matA.GetEntries()[1] = 1;
    matA.GetEntries()[2] = 1;
    matA.GetEntries()[3] = 1;
    matA.GetEntries()[4] = 1;

    BlockMatrix myMatrix({4 , 2 , 4}, {4 , 2, 4});

    myMatrix.replace_blocks(
        matA,
        {{0,0}, {0,2}, {2,0}, {2,2}},
        {false, false, true, false} );

    myMatrix.print_and_check("after replace");
    myMatrix.print_coloring_pattern("myMatrix", true);
    //myMatrix.get_combined_matrix()->PrintFull();

    //check copy constructor
    BlockMatrix yourMatrix(myMatrix);
    yourMatrix.print_and_check("yourMatrix");
    yourMatrix.print_coloring_pattern("yourMatrix", true);

    TMatrix summand(matA);
    summand.GetEntries()[0] = 1;
    summand.GetEntries()[1] = 1;
    summand.GetEntries()[2] = 1;
    summand.GetEntries()[3] = 1;
    summand.GetEntries()[4] = 1;

    myMatrix.add_unscaled_matrix(summand,
                                  {{2,0}, {2,2}},
                                  {true, false}  );

    myMatrix.print_and_check("addition to replace");
    myMatrix.get_combined_matrix()->PrintFull();

    myMatrix.add_unscaled_matrix(
        summand ,
        { {2,0}, {2,2} },
        {true, false} );

    myMatrix.print_and_check("addition transpose");
    myMatrix.get_combined_matrix()->PrintFull();

    std::vector<std::vector<size_t>> positions = {{0,0},{2,0}};
    std::vector<bool> transp_states = {false,true};

    myMatrix.add_matrix(
        summand, 1.0, positions, transp_states);


    myMatrix.print_and_check("addition 2");
    myMatrix.get_combined_matrix()->PrintFull();

    //check copy assignment
    BlockMatrix yourAssignedMatrix({0},{0});
    yourAssignedMatrix = myMatrix;
    yourAssignedMatrix.print_and_check("yourAssignedMatrix");
    yourAssignedMatrix.print_coloring_pattern("yourAssignedMatrix", true);


    delete[] RowA;
    delete[] ColA;
  }

  // check constructor taking existing TMatrix blocks directly
  {
    Output::print("\n\nstarting to test the BlockMatrix constructor taking ",
                  "existing blocks\n");
    /**
     * A = [ 1,   0,    0, 0  ; B = [ 4, 0, 0, 1; C = [ 7, 0, 8;
     *       0,   0.75, 0, 0  ;       0, 0, 5, 0;       9,10, 0;
     *       0,   0,    1, 0  ;       0, 6, 0, 1 ]      0,11,12;
     *       0,   0,    0, 1.5 ]                        13,0,14 ]
     *
     * D = [ 0.5 0   0  ;
     *       0   1   0.1;
     *       0   0.1 1.5 ]
     * 
     * !! the entries are not written, they are not important here
     */
    // create a few TMatrix objects
    //Matrix A
    int * RowA = new int[5];
    int * ColA = new int[4];
    ColA[0] = 0; ColA[1] = 1; ColA[2] =  2; ColA[3] = 3;
    RowA[0] = 0; RowA[1] = 1; RowA[2] = 2; RowA[3] = 3; RowA[4] = 4;
    std::shared_ptr<TStructure> structureA(new TStructure(4, 4, ColA, RowA));
    std::shared_ptr<TMatrix> matA = std::make_shared<TMatrix>(structureA);
    matA->GetEntries()[0] = 10.; // set just one entry for testing later
    matA->GetEntries()[1] = 75;
    matA->GetEntries()[2] = 100;
    matA->GetEntries()[3] = 150;
    
    //Matrix B
    int * RowB = new int[4];
    int * ColB = new int[5];
    RowB[0] = 0; RowB[1] = 2; RowB[2] = 3; RowB[3] = 5;
    ColB[0] = 0; ColB[1] = 3; ColB[2] = 2; ColB[3] = 1; ColB[4] = 3;
    std::shared_ptr<TStructure> structureB(new TStructure(3, 4, 5, ColB, RowB));
    std::shared_ptr<TMatrix> matB = std::make_shared<TMatrix>(structureB);
    matB->GetEntries()[0] = 4;
    matB->GetEntries()[1] = 1;
    matB->GetEntries()[2] = 5;
    matB->GetEntries()[3] = 6;
    matB->GetEntries()[4] = 1;
    
    //Matrix C
    int * RowC = new int[5];
    int * ColC = new int[8];
    RowC[0] = 0; RowC[1] = 2; RowC[2] = 4; RowC[3] = 6; RowC[4] = 8;
    ColC[0] = 0; ColC[1] = 2; ColC[2] = 0; ColC[3] = 1; ColC[4] = 1; ColC[5] =2;
    ColC[6] = 0; ColC[7] = 2;
    std::shared_ptr<TStructure> structureC(new TStructure(4, 3, 8, ColC, RowC));
    std::shared_ptr<TMatrix> matC = std::make_shared<TMatrix>(structureC);
    matC->GetEntries()[0] = 7;
    matC->GetEntries()[1] = 8;
    matC->GetEntries()[2] = 9;
    matC->GetEntries()[3] = 10;
    matC->GetEntries()[4] = 11;
    matC->GetEntries()[5] = 12;
    matC->GetEntries()[6] = 13;
    matC->GetEntries()[7] = 14;
    
    std::shared_ptr<TMatrix> matCT(matC->get_transposed());

    //Matrix D
    int * RowD = new int [4];
    int * ColD = new int [5];
    RowD[0] = 0; RowD[1] = 1; RowD[2] = 3; RowD[3] = 5;
    ColD[0] = 0; ColD[1] = 1; ColD[2] = 2; ColD[3] = 1; ColD[4] = 2;
    std::shared_ptr<TStructure> structureD(new TStructure(3, 3, 5, ColD, RowD));
    std::shared_ptr<TMatrix> matD = std::make_shared<TMatrix>(structureD);
    matD->GetEntries()[0] = 50.;
    matD->GetEntries()[1] = 10.;
    matD->GetEntries()[2] = 0.1;
    matD->GetEntries()[3] = 0.1;
    matD->GetEntries()[4] = 15.;
    
    // create a BlockMatrix
    std::vector<std::shared_ptr<TMatrix>> blocks = {matA, matC, matB, matD};
    BlockMatrix bm(2, 2, blocks);
    
    // do some simple tests
    if(bm.get_n_total_rows() != 7)
      ErrThrow("total number of rows in BlockMatrix is incorrect");
    if(bm.get_n_total_columns() != 7)
      ErrThrow("total number of columns in BlockMatrix is incorrect");
    if(bm.get_n_total_entries() != 22)
      ErrThrow("total number of entries in BlockMatrix is incorrect");
    // check an entry in the sparsity pattern
    if(bm.get(0, 0) != 10.)
      ErrThrow("wrong entry at position (0,0) ", bm.get(0, 0));
    // check an entry which is not in the sparsity pattern
    if(bm.get(4, 1) != 0.)
      ErrThrow("wrong entry at position (0,0) ", bm.get(0, 0));
    // check if exceptions are thrown when accessing a value which is too large
    try
    {
      double a = bm.get(7, 7); // the size is 6 by 6
      Output::print("I was able to get an entry where there should not even "
                    "exist one: ", a);
      return 1;
    }
    catch(...)
    {
      // correct behavior, nothing more to do
    }

    // check if exceptions are thrown when appropriate
    try
    {
      blocks = {matA, matCT, matB, matD};
      BlockMatrix bm_wrong(2, 2, blocks);
      Output::print("It was possible to create a BlockMatrix where the blocks ",
                    "did not match in size. This should have thrown an ",
                    "exception");
      return 1; // failure
    }
    catch(...)
    {
      // correct behavior, nothing more to do
    }
    
    try
    {
      BlockMatrix bm_wrong(3, 2, blocks);
      Output::print("It was possible to create a BlockMatrix where not enough ",
                    "blocks were given. This should have thrown an exception");
      return 1; // failure
    }
    catch(...)
    {
      // correct behavior, nothing more to do
    }
    
    try
    {
      
      blocks = {matA, matC, matB, std::shared_ptr<TMatrix>(nullptr), matD};
      BlockMatrix bm_wrong(2, 2, blocks);
      Output::print("It was possible to create a BlockMatrix where one block ",
                    "is a nullptr. This should have thrown an exception");
      return 1; // failure
    }
    catch(...)
    {
      // correct behavior, nothing more to do
    }

    {//test the sub-blockmatrix creation routine
      BlockMatrix bm({2,2,4,6},{2,4,6});

      //fill in some matrices
      bm.replace_blocks(TMatrix(2,2) , {{1,0}},{false});
      bm.replace_blocks(TMatrix(4,2) , {{2,0},{1,1}},{false, true});
      bm.replace_blocks(TMatrix(4,4) , {{2,1}},{false});

      BlockMatrix sub_bm;
      sub_bm = bm.get_sub_blockmatrix({1,0}, {2,1});

      sub_bm.print_and_check("sub_mb");

      bm.get_combined_submatrix({0,0},{1,2})->info(3);

    }

    {
      // test BlockMatrix::sor_sweep
      BlockVector rhs(bm, true);
      rhs[0] = 1.0;
      BlockVector sol(bm, false);
      BlockVector exact_sol(bm, false);
      exact_sol[0] = 0.1063356097695721;
      exact_sol[1] = 0.001025959742357445;
      exact_sol[2] = 5.817149251470894e-05;
      exact_sol[3] = 0.0007817668937139531;
      exact_sol[4] = -0.008522484119439989;
      exact_sol[5] = -2.446236018275369e-05;
      exact_sol[6] = -0.0004623386074561225;
      bm.get_combined_matrix()->PrintFull();
      for(size_t i = 0; i < 20; ++i)
      {
        bm.sor_sweep(rhs, sol, 1.1, 2);
        BlockVector diff(sol);
        diff.add_scaled(exact_sol, -1.);
        Output::print("norm ", i, ": ", diff.norm());
      }
    }
    
    delete[] RowA;
    delete[] ColA;
    delete[] RowB;
    delete[] ColB;
    delete[] RowC;
    delete[] ColC;
    delete[] RowD;
    delete[] ColD;
  }
  
  Output::print("Test program finished.");

}


