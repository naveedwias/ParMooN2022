/*
 * Unit testing of the block matrix. Since this requires
 * the presence of FEspaces, this starts as a usual ParMooN program,
 * which is needed to build up FESpaces which can than be used for
 * the BlockFEMatrix.
 *
 *
 *  Created on: Dec 17, 2015
 *      Author: bartsch
 */
#include <cmath>
#include <vector>
#include <tuple>
#include <type_traits>


#include <BlockFEMatrix.h>
#include <FEMatrix.h>

#include <MooNMD_Io.h>
#include <Database.h>
#include <FESpace2D.h>
#include <MainUtilities.h>
#include <BlockVector.h>
#include "ParMooN.h"


int main(int, char**)
{
  ParameterDatabase db = parmoon::parmoon_initialize();
  db.merge(ParameterDatabase::default_nonlinit_database());

  db.add("refinement_n_initial_steps", (size_t) 2, "");
  db.add("multigrid_n_levels", (size_t) 0, "");

  // default construct a domain object
  db.add("boundary_file", "Default_UnitSquare", "");
  db.add("geo_file", "UnitSquare", "", {"UnitSquare", "TwoTriangles"});
  TDomain domain(db);

  Output::setVerbosity(1);

  // refine grid up to the coarsest level
  size_t n_ref = domain.get_n_initial_refinement_steps();
  for(unsigned int i=0; i < n_ref; i++)
  {
    domain.RegRefineAll();
  }

  //collection from finest cell level
  TCollection *coll = domain.GetCollection(It_EQ, 0);
  // Create FeSpaces2D to fiddle around with

  size_t first_ansatz_order = 2;
  size_t second_ansatz_order = 1;
  size_t third_ansatz_order = 3;

  std::shared_ptr<const TFESpace2D> first_fe_space(
    new TFESpace2D(coll, "first_fe_space", //may act as velo space dummy
                  BoundConditionNSE, first_ansatz_order));
  
  std::shared_ptr<const TFESpace2D> second_fe_space(
    new TFESpace2D(coll, "second_fe_space", //may act as pressure space dummy
                   BoundCondition_FEM_FCT, second_ansatz_order));

  std::shared_ptr<const TFESpace2D> third_fe_space(
    new TFESpace2D(coll, "third_fe_space", //yet another space
                   BoundConditionNSE, third_ansatz_order));
  Output::print("grid has ", coll->GetN_Cells(), " cells");
  Output::print("first_fe_space has ", first_fe_space->get_n_dof(), " dofs");
  Output::print("second_fe_space has ", second_fe_space->get_n_dof(), " dofs");
  Output::print("third_fe_space has ", third_fe_space->get_n_dof(), " dofs");
  
  {
    //test default constructor
    BlockFEMatrix zero_matrix;
    zero_matrix.check_coloring();
    zero_matrix.check_pointer_types();
  }
  { // test standard methods with custom-made 2x2 FEMatrix, including
    // one transposed storage and one transposed-storage memory hack

    //create four FE Matrices to fiddle around with
    FEMatrix fe_matrix_1(first_fe_space);
    FEMatrix fe_matrix_2(first_fe_space, second_fe_space);
    FEMatrix fe_matrix_3(second_fe_space, first_fe_space);
    FEMatrix fe_matrix_4(second_fe_space);

    // one TMatrix, too (copied from one FEMatrix for the sake of a non-emtpy structure)
    TMatrix t_matrix_1(fe_matrix_3);
    t_matrix_1.setEntries(std::vector<double>(t_matrix_1.get_n_entries(),1.0));

   //custom construct and check
   BlockFEMatrix myMatrix({first_fe_space, second_fe_space});
   myMatrix.check_pointer_types(); //casts to FEMatrix work?
   myMatrix.check_coloring(); //coloring is unbroken?
   
   FEMatrix fe_matrix_4_copy(second_fe_space, static_cast<TMatrix&>(fe_matrix_4));

   //replace method and check
   myMatrix.replace_blocks(fe_matrix_1, {{0,0}}, { false } );
   myMatrix.replace_blocks(fe_matrix_4_copy, {{1,1}}, { false });

   myMatrix.replace_blocks(fe_matrix_3, {{0,1}, {1,0}}, { true, false }); //replace which leads to color merge
   myMatrix.replace_blocks(fe_matrix_2, {{0,1}}, {false} ); //replace which leads to color split

   myMatrix.check_pointer_types(); //casts to FEMatrix work?
   myMatrix.check_coloring(); //coloring is unbroken?

   myMatrix.replace_blocks(fe_matrix_3, {{0,1}, {1,0}}, { true, false }); //replace which leads to color merge
   myMatrix.check_pointer_types(); //casts to FEMatrix work?
   myMatrix.check_coloring(); //coloring is unbroken?

   // try combination method of base
   {
     std::shared_ptr<TMatrix> combined_base
     = myMatrix.BlockMatrix::get_combined_matrix();
//     combined_base->PrintFull("combined base");
     if(combined_base->GetNorm(-2) != 0)
     {
       ErrThrow("Wrong norm of base-class combined matrix!");
     }
   }

   //try combination method of derived
   {
     std::shared_ptr<TMatrix> combined_derived
     = myMatrix.BlockFEMatrix::get_combined_matrix();
     //combined_derived->PrintFull("combined derived");
     if( std::abs(combined_derived->GetNorm(-2) - 3.) > 1e-5) // sqrt(9)
     {
       ErrThrow("Wrong norm of derived-class combined matrix!");
     }
   }

   // test the space getter methods
   if(myMatrix.get_ansatz_space(0,0) != first_fe_space )
   {
     ErrThrow("get_ansatz_space not working correctly.")
   }
   if(myMatrix.get_test_space(0,1) != first_fe_space )
   {
     ErrThrow("get_ansatz_space not working correctly.")
   }
   if(myMatrix.get_row_space(1) != second_fe_space )
   {
     ErrThrow("get_row_space not working correctly.")
   }
   if(myMatrix.get_column_space(1) != second_fe_space )
   {
     ErrThrow("get_column_space not working correctly.")
   }

   // check standard adding (no color split intended), here with a simple TMatrix
   myMatrix.add_matrix(t_matrix_1, 1.0 ,{{1,0}, {0,1}}, {false, true});
   myMatrix.check_pointer_types(); //casts to FEMatrix work?
   myMatrix.check_coloring(); //coloring is unbroken?
   if( std::abs( myMatrix.get_combined_matrix()->GetNorm(-2) - std::sqrt(40.) ) > 1e-10)
   {//check norm of dirichlet-row corrected matrix
     ErrThrow("Wrong matrix norm after adding of a TMatrix!")
   }
   if( std::abs( myMatrix.BlockMatrix::get_combined_matrix()->GetNorm(-2) - std::sqrt(72) ) > 1e-10 )
   {//check norm of rough, algebraic matrix
     ErrThrow("Wrong matrix norm after adding of a TMatrix!")
   }

   //check actives adding (color split intended) with an FEMatrix
   fe_matrix_1.setEntries(std::vector<double>(fe_matrix_1.get_n_entries(), 2.0));
   myMatrix.add_matrix_actives(fe_matrix_1, -2.0 ,{{0,0}}, { false });
   myMatrix.check_pointer_types(); //casts to FEMatrix work?
   myMatrix.check_coloring(); //coloring is unbroken?

   if( std::abs( myMatrix.get_combined_matrix()->GetNorm(-2) - std::sqrt(184.) ) > 1e-10 )
   {//check norm of dirichlet-row corrected matrix
     ErrThrow("Wrong matrix norm after actives adding!")
   }
   if( std::abs( myMatrix.BlockMatrix::get_combined_matrix()->GetNorm(-2) - std::sqrt(216.) ) > 1e-10 )
   {//check norm of rough, algebraic matrix
     ErrThrow("Wrong matrix norm after actives adding!")
   }


   //check scaling of entries
   const std::vector<std::vector<size_t>> cell_positions = {{0,0},{1,1}};
     myMatrix.scale_blocks(0.5, cell_positions);
   myMatrix.check_pointer_types(); //casts to FEMatrix work?
   myMatrix.check_coloring(); //coloring is unbroken?

   if( std::abs( myMatrix.get_combined_matrix()->GetNorm(-2) - std::sqrt(76.) ) > 1e-10 )
   {//check norm of dirichlet-row corrected matrix
     ErrThrow("Wrong matrix norm after block scaling!")
   }
   if( std::abs( myMatrix.BlockMatrix::get_combined_matrix()->GetNorm(-2) - std::sqrt(108.) ) > 1e-10 )
   {//check norm of rough, algebraic matrix
     ErrThrow("Wrong matrix norm after block scaling!")
   }

   //check scaling of active entries
   myMatrix.replace_blocks(fe_matrix_1, {{0,0}}, { false });//give us a new block in {0,0}
   const std::vector<std::vector<size_t>> cell_positions_2 = {{0,0}};
   myMatrix.scale_blocks_actives(-2, cell_positions_2);
   myMatrix.check_pointer_types(); //casts to FEMatrix work?
   myMatrix.check_coloring(); //coloring is unbroken?

   if( std::abs( myMatrix.get_combined_matrix()->GetNorm(-2) - std::sqrt(184.) ) > 1e-10 )
   {//check norm of dirichlet-row corrected matrix
     ErrThrow("Wrong matrix norm after block scaling!")
   }
   if( std::abs( myMatrix.BlockMatrix::get_combined_matrix()->GetNorm(-2) - std::sqrt(504.) ) > 1e-10 )
   {//check norm of rough, algebraic matrix
     ErrThrow("Wrong matrix norm after block scaling!")
   }


  }

  {
    // test standard methods with custom-made 2x2 FEMatrix, including
    // one transposed storage and one transposed-storage memory hack

    //create four FE Matrices to fiddle around with
    auto fe_mat1 = std::make_shared<FEMatrix>(first_fe_space);
    auto fe_mat2 = std::make_shared<FEMatrix>(first_fe_space, second_fe_space);
    auto fe_mat3 = std::make_shared<FEMatrix>(second_fe_space, first_fe_space);
    auto fe_mat4 = std::make_shared<FEMatrix>(second_fe_space);
    auto fe_mat4_copy = std::make_shared<FEMatrix>(second_fe_space,
                                                   static_cast<TMatrix&>(
                                                     *fe_mat4));
    // custom construct and check
    BlockFEMatrix myMatrix(2, 2, {fe_mat1, fe_mat2, fe_mat3, fe_mat4_copy});
    myMatrix.check_pointer_types(); //casts to FEMatrix work?
    myMatrix.check_coloring(); //coloring is unbroken?
    
  }
  
  { //make a default NSE Matrix and two vectors
    BlockFEMatrix blockmat=
            BlockFEMatrix::NSE2D_Type1(first_fe_space, second_fe_space);

    //check matrix-vector multiplication (incl. actives)
    BlockVector preimage_act(blockmat, false);
    BlockVector image_act(blockmat,true);

    for(size_t i = 0 ; i< preimage_act.length() ; ++i)
    { //fill one vector with with ones
      preimage_act.at(i) = 1;
    }

    blockmat.apply(preimage_act, image_act);

    if(image_act.norm() != 4)
    {
      ErrThrow("Norm of BlockVector from multiplication is not correct!");
    }

    //check usage in std::vector
    std::vector<BlockFEMatrix> myMatrices;
    myMatrices.push_back(blockmat);
    myMatrices.push_back(blockmat);
    myMatrices.push_back(blockmat); //three pushes
    myMatrices.pop_back(); //one pop
    myMatrices.at(0).check_pointer_types(); //random access and check element
    myMatrices.at(0).check_coloring();


    //check copying
    Output::setVerbosity(5);
    Output::print("Check copying.");
    //copy construction
    BlockFEMatrix hisMatrix(blockmat);
    //hisMatrix.print_coloring_pattern("copy constructed fe matrix",true);
    hisMatrix.check_pointer_types();
    hisMatrix.check_coloring();

    //copy assignment
    BlockFEMatrix herMatrix({first_fe_space});
    herMatrix = blockmat;
    //herMatrix.print_coloring_pattern("copy assigned fe matrix",true);
    herMatrix.check_pointer_types();
    herMatrix.check_coloring();

    //check moving


    Output::print("Test moving.");

    //move constructor
    //BlockFEMatrix moveConstructedMatrix(BlockFEMatrix::NSE2D_Type14(first_fe_space, second_fe_space));
    BlockFEMatrix moveConstructedMatrix(
      std::move(BlockFEMatrix(
        std::vector<std::shared_ptr<const TFESpace2D>>({{first_fe_space, second_fe_space}}))));
    moveConstructedMatrix.check_pointer_types();
    moveConstructedMatrix.check_coloring();

    Output::print("Test move assigning.");

    BlockFEMatrix moveAssignedMatrix({first_fe_space});
    moveAssignedMatrix = BlockFEMatrix::Darcy2D(first_fe_space, second_fe_space);
    moveAssignedMatrix.check_pointer_types();
    moveAssignedMatrix.check_coloring();

    Output::setVerbosity(0);
  }


  // try out some named constructors, plus the block getter
  // methods for solvers and assemblers
  {//CD2D
    BlockFEMatrix blockmat = BlockFEMatrix::CD2D(third_fe_space);
    //blockmat.print_and_check("matrix for cd2d");
    blockmat.check_pointer_types(); //casts to FEMatrix work?
    blockmat.check_coloring(); //coloring is unbroken?
    // assemble blocks getter
    std::vector<std::shared_ptr<FEMatrix>> blocks_for_assembler
    = blockmat.get_blocks_uniquely();
    if (blocks_for_assembler.size() != 1)
    {
      ErrThrow("Incorrect blocks_for_assembler.size() !");
    }
    //solver blocks getter
    std::vector<std::shared_ptr<const FEMatrix>> blocks_for_solver
    = blockmat.get_blocks();
    if (blocks_for_assembler.size() != 1)
    {
      ErrThrow("Incorrect blocks_for_solver.size() !");
    }

  }
   {//NSE2D Typ 1
     BlockFEMatrix blockmat=
         BlockFEMatrix::NSE2D_Type1(first_fe_space, second_fe_space);
     //blockmat.print_and_check("matrix for nse2d, nstype 1");
     blockmat.check_pointer_types(); //casts to FEMatrix work?
     blockmat.check_coloring(); //coloring is unbroken?
     // assemble blocks getter
     std::vector<std::shared_ptr<FEMatrix>> blocks_for_assembler
     = blockmat.get_blocks_uniquely();
     if (blocks_for_assembler.size() != 3)
     {
       ErrThrow("Incorrect blocks_for_assembler.size() !");
     }
     //solver blocks getter
     std::vector<std::shared_ptr<const FEMatrix>> blocks_for_solver
     = blockmat.get_blocks();
     if (blocks_for_solver.size() != 9)
     {
       ErrThrow("Incorrect blocks_for_solver.size() !");
     }
   }
   {//NSE2D Typ 2
     BlockFEMatrix blockmat=
         BlockFEMatrix::NSE2D_Type2(first_fe_space, second_fe_space);
     //blockmat.print_and_check("matrix for nse2d, nstype 2");
     blockmat.check_pointer_types(); //casts to FEMatrix work?
     blockmat.check_coloring(); //coloring is unbroken?
     // assemble blocks getter
     std::vector<std::shared_ptr<FEMatrix>> blocks_for_assembler
     = blockmat.get_blocks_uniquely();
     if (blocks_for_assembler.size() != 5)
     {
       ErrThrow("Incorrect blocks_for_assembler.size() !");
     }
     //solver blocks getter
     std::vector<std::shared_ptr<const FEMatrix>> blocks_for_solver
     = blockmat.get_blocks();
     if (blocks_for_solver.size() != 9)
     {
       ErrThrow("Incorrect blocks_for_solver.size() !");
     }
   }
   {//NSE2D Typ 3
     BlockFEMatrix blockmat=
         BlockFEMatrix::NSE2D_Type3(first_fe_space, second_fe_space);
     //blockmat.print_and_check("matrix for nse2d, nstype 3");
     blockmat.check_pointer_types(); //casts to FEMatrix work?
     blockmat.check_coloring(); //coloring is unbroken?
     // assemble blocks getter
     std::vector<std::shared_ptr<FEMatrix>> blocks_for_assembler
     = blockmat.get_blocks_uniquely();
     if (blocks_for_assembler.size() != 6)
     {
       ErrThrow("Incorrect blocks_for_assembler.size() !");
     }
     //solver blocks getter
     std::vector<std::shared_ptr<const FEMatrix>> blocks_for_solver
     = blockmat.get_blocks();
     if (blocks_for_solver.size() != 9)
     {
       ErrThrow("Incorrect blocks_for_solver.size() !");
     }
   }
   {//NSE2D Typ 4
     BlockFEMatrix blockmat=
         BlockFEMatrix::NSE2D_Type4(first_fe_space, second_fe_space);
     //blockmat.print_and_check("matrix for nse2d, nstype 4");
     blockmat.check_pointer_types(); //casts to FEMatrix work?
     blockmat.check_coloring(); //coloring is unbroken?
     // assemble blocks getter
     std::vector<std::shared_ptr<FEMatrix>> blocks_for_assembler
     = blockmat.get_blocks_uniquely();
     if (blocks_for_assembler.size() != 8)
     {
       ErrThrow("Incorrect blocks_for_assembler.size() !");
     }
     //solver blocks getter
     std::vector<std::shared_ptr<const FEMatrix>> blocks_for_solver
     = blockmat.get_blocks();
     if (blocks_for_solver.size() != 9)
     {
       ErrThrow("Incorrect blocks_for_solver.size() !");
     }
   }
   {//NSE2D Typ 14
     BlockFEMatrix blockmat=
         BlockFEMatrix::NSE2D_Type14(first_fe_space, second_fe_space);
     //blockmat.print_and_check("matrix for nse2d, nstype 14");
     blockmat.check_pointer_types(); //casts to FEMatrix work?
     blockmat.check_coloring(); //coloring is unbroken?
     // assemble blocks getter
     std::vector<std::shared_ptr<FEMatrix>> blocks_for_assembler
     = blockmat.get_blocks_uniquely();
     if (blocks_for_assembler.size() != 9)
     {
       ErrThrow("Incorrect blocks_for_assembler.size() !");
     }
     //solver blocks getter
     std::vector<std::shared_ptr<const FEMatrix>> blocks_for_solver
     = blockmat.get_blocks();
     if (blocks_for_solver.size() != 9)
     {
       ErrThrow("Incorrect blocks_for_solver.size() !");
     }
   }
   {//Darcy2D
     BlockFEMatrix blockmat=
         BlockFEMatrix::Darcy2D(first_fe_space, second_fe_space);
     //blockmat.print_and_check("matrix for darcy 2d");
     blockmat.check_pointer_types(); //casts to FEMatrix work?
     blockmat.check_coloring(); //coloring is unbroken?
     // assemble blocks getter
     std::vector<std::shared_ptr<FEMatrix>> blocks_for_assembler
     = blockmat.get_blocks_uniquely();
     if (blocks_for_assembler.size() != 4)
     {
       ErrThrow("Incorrect blocks_for_assembler.size() !");
     }
     //solver blocks getter
     std::vector<std::shared_ptr<const FEMatrix>> blocks_for_solver
     = blockmat.get_blocks();
     if (blocks_for_solver.size() != 4)
     {
       ErrThrow("Incorrect blocks_for_solver.size() !");
     }
   }

   {//more tests on block getter methods (assemble parts getter, "true" for zero blocks)
     BlockFEMatrix blockmat= //use NSE Type 2 for the checks
         BlockFEMatrix::NSE2D_Type2(first_fe_space, second_fe_space);

     std::vector<std::vector<size_t>> positions_A = {{0,0},{0,1},{1,0},{1,1}};
     std::vector<std::vector<size_t>> positions_B = {{0,2},{1,2},{2,0},{2,1}};
     // adapt and use that if you want to check the handling of not-so clever input
     //std::vector<std::vector<size_t>> stupid_input = {{1,1}};

     // positionwise assemble getter
     std::vector<std::shared_ptr<FEMatrix>> A_blocks_for_assembler
       = blockmat.get_blocks_uniquely(positions_A);
     if(A_blocks_for_assembler.size() != 1)
     {
       ErrThrow("Incorrect A_blocks_for_assembler.size() !");
     }
     std::vector<std::shared_ptr<FEMatrix>> A_blocks_for_assembler_with_zeroes
       = blockmat.get_blocks_uniquely(positions_A, true); //with "true" for zeroes
     if(A_blocks_for_assembler_with_zeroes.size() != 2)
     {
       ErrThrow("Incorrect A_blocks_for_assembler_with_zeroes.size() !");
     }
     std::vector<std::shared_ptr<FEMatrix>> B_blocks_for_assembler
       = blockmat.get_blocks_uniquely(positions_B);
     if(B_blocks_for_assembler.size() != 4)
     {
       ErrThrow("Incorrect B_blocks_for_assembler.size() !");
     }
   }

   {
     Output::setVerbosity(1);
    // this is a test on the apply_scaled_add method which is suspected to produce
    // wrong results when transposition is involved
     BlockFEMatrix with_transp =
         BlockFEMatrix::NSE2D_Type1(first_fe_space, second_fe_space);
     BlockFEMatrix no_transp =
         BlockFEMatrix::NSE2D_Type4(first_fe_space, second_fe_space);
     //fill some blocks in
     FEMatrix A(first_fe_space, first_fe_space);
     FEMatrix B(second_fe_space, first_fe_space);
     FEMatrix BT(first_fe_space, second_fe_space);
     A.setEntries(std::vector<double>(A.get_n_entries(), 1.0));
     B.setEntries(std::vector<double>(B.get_n_entries(), 1.0));
     BT.setEntries(std::vector<double>(BT.get_n_entries(), 1.0));

     //do some careful replacing, maintaining the coloring
     with_transp.replace_blocks(A, {{0,0},{1,1}}, {false, false});
     with_transp.replace_blocks(B,{{0,2},{2,0}},{true,false});
     with_transp.replace_blocks(B,{{1,2},{2,1}},{true,false});

     no_transp.replace_blocks(A, {{0,0}}, {false});
     //no_transp.replace_blocks(A, {{0,1}}, {false});
     //no_transp.replace_blocks(A, {{1,0}}, {false});
     no_transp.replace_blocks(A, {{1,1}}, {false});
     no_transp.replace_blocks(BT,{{0,2}}, {false});
     no_transp.replace_blocks(BT,{{1,2}}, {false});
     no_transp.replace_blocks(B,{{2,0}}, {false});
     no_transp.replace_blocks(B,{{2,1}}, {false});

     //put up BlockVectors which fit the spaces
     BlockVector preimage(with_transp, false);
     for(size_t i = 0 ; i< preimage.length() ; ++i)
     { //fill the preimage vector with with ones
       preimage.at(i) = 1;
     }
     BlockVector image_with(with_transp, true);
     BlockVector image_no(no_transp, true);

     with_transp.apply_scaled_add(preimage, image_with, -1.0);
     //image_with.print("image_with");

     no_transp.apply_scaled_add(preimage, image_no, -1.0);
     //image_no.print("image_no");
     image_with.add_scaled(image_no,-1.0);
     if( image_with.norm() != 0)
       ErrThrow("Error in test of apply_scaled_add!");

   }

   {
     //test the method which returns a sub-blockfematrix
     BlockFEMatrix blockmat=
         BlockFEMatrix::NSE2D_Type1(first_fe_space, second_fe_space);

     BlockFEMatrix sub_velocity = blockmat.get_sub_blockfematrix(0,1);
     sub_velocity.check_pointer_types();
     sub_velocity.check_coloring();

     BlockFEMatrix sub_1D = blockmat.get_sub_blockfematrix(1,2);
     sub_1D.check_pointer_types();
     sub_1D.check_coloring();


   }

  parmoon::parmoon_finalize();
}

