#include <PETScSolver.h>
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <MooNMD_Io.h>
#ifdef _MPI
#include <mpi.h>
#include <ParFECommunicator3D.h>
#endif // _MPI
#include <algorithm>

#ifdef PARMOON_WITH_PETSC
/**
 * Adding PETSc parameters from the database to the PETSc arguments string.
 *
 * Supported parameters are
 *  - maximum number of iterations
 *  - absolute residual tolerance
 *  - residual reduction
 *  - gmres restart
 *  - damping factor for Richardson iteration
 *
 * If a parameters is set in this function it will take precedence over
 * the parameter in the db["petsc_arguments"] string, e.g. if
 * your .dat file reads
 *
 * max_n_iterations: 1000
 *
 * [...]
 *
 * petsc_arguments: -pc_type lu -ksp_max_it 5
 *
 * The 1000 will be set in this function and your code will run with
 * 1000 iterations maximum instead of 5.
 *
 * @param db[in] database from the solver
 * @return PETSc arguments conforming string
 */
std::string addParameters2str(const ParameterDatabase& db)
{
  std::string petsc_args = "";
  if(db.contains("petsc_arguments"))
  {
    petsc_args = db["petsc_arguments"].value_as_string();
  }
  // adding maximum number of iterations
  petsc_args.append(" -ksp_max_it "+db["max_n_iterations"].value_as_string());
  // adding absolute tolerance
  petsc_args.append(" -ksp_atol "+db["residual_tolerance"].value_as_string());
  // adding reduction tolerance
  petsc_args.append(" -ksp_rtol "+db["residual_reduction"].value_as_string());
  // adding gmres restart
  petsc_args.append(" -ksp_gmres_restart "+db["gmres_restart"].value_as_string());
  // adding sor omega
  // For me (Ulrich) petsc did not do a single iteration if omega was not set
  // to one. I have no clue why, but it could be fixed by adding the extra 
  // option "-ksp_error_if_not_converged 1".
    
  petsc_args.append(" -pc_sor_omega "+db["sor_omega"].value_as_string()
                    + " -ksp_error_if_not_converged 1");
  // adding damping for Richardson iteration
  petsc_args.append(" -ksp_richardson_scale "+db["damping_factor"].value_as_string());
  // tell PETSc to start with a nonzero solution. This is as the solution 
  // BlockVector in the PETScSolver::solve method. Sometimes this solution is
  // zero, but this does not hurt.
  petsc_args.append(" -ksp_initial_guess_nonzero");

    
//// AENDERUNG
//   petsc_args.append(" -pc_type "+db["preconditioning"].value_as_string());
//   petsc_args.append(" -ksp_type "+db["Krylov Subspace Method"].value_as_string());
//     Output::print('Hier Gucken');
//    Output::print(db["Krylov Subspace Method"].value_as_string());
//----------
                     
  return petsc_args;
}

/**
 * @brief Converts a std::string to a vector of char* by splitting words at 
 * white spaces.
 * 
 * The idea is that the returned vector mimics the standard arguments of a C/C++
 * main function: int main(int argc, char** argv). These arguments are needed to
 * properly call PetscInitialize.
 * 
 * @param[in] s input string
 * @return vector of char* containing all words + nullptr
 */
std::vector<char*> str_to_vector_char_p(const std::string &s)
{
  // number of arguments, the first is the program name by default
  size_t n_arg = 1;
  // a position in the string s
  auto current_position = s.find(" ");
  // count the white spaces and therefore the number of arguments
  while(current_position != std::string::npos)
  {
    n_arg++;
    current_position = s.find(" ", current_position + 1);
  }
  // the number of white spaces is smaller than the number arguments by 1
  n_arg++;
  
  // store where words (separated by words) are
  std::vector<size_t> word_positions;
  word_positions.reserve(n_arg+1);
  // zero is always the first position
  word_positions.push_back(0u);

  // fill the vector word_positions
  current_position = s.find(" ");
  while(current_position != std::string::npos)
  {
    word_positions.push_back(current_position + 1);
    current_position = s.find(" ", current_position + 1);
  }

  // push the length of the string as last "position"
  // now word_positions[i+1]-word_positions[i] is valid for i = 0, ... n_arg-1
  word_positions.push_back(s.size()+1);

  // plus 1 here for the terminating nullptr, this will be returned
  std::vector<char*> argv(n_arg+1);

  // PETSc expects the program name at first.
  char* progname = new char[8];
  strcpy(progname, "./dummy");
  argv[0] = &progname[0];

  for(size_t i = 1; i < n_arg; ++i)
  {
    size_t wordLength = word_positions[i] - word_positions[i-1] - 1;
    argv[i] = new char[wordLength+1];
    strncpy( argv[i], &s.c_str()[word_positions[i-1]], wordLength );
    // add terminating \0
    argv[i][wordLength] = '\0';
  }

  // set the nullptr at the end
  // demanded by the C++ Standard for main function argument argv
  argv[n_arg] = nullptr;
  return argv;
}

// *****************************************************************************
// Return a vector n_masters whose size is the number of mpi processes. For 
// each process p the vector n_master[p] has the size of the input vector. The
// element n_master[p][b] is the number of master dofs in the b-th block in the
// p-th process.
// This requires communication!
#ifdef _MPI
std::vector<std::vector<size_t>> get_n_masters_per_block_and_process(
  std::vector<const TParFECommunicator3D*> comms)
{
  int my_rank;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  // the number of communicators is the number of blocks
  size_t n_comms = comms.size();
  // the number of master dofs for each process and each block
  int n_masters_send[size][n_comms];
  // set all elements to 0
  for(size_t r = 0; (int)r < size; ++r)
    for(size_t c = 0; c < n_comms; ++c)
      n_masters_send[r][c] = 0;
  
  for(size_t c = 0; c < n_comms; ++c)
  {
    // each process only writes one row, others are left zero
    n_masters_send[my_rank][c] = comms[c]->GetN_Master();
  }
  // n_masters needs to be communicated
  int n_masters[size][n_comms];
   // set all elements to 0
  for(size_t r = 0; (int)r < size; ++r)
    for(size_t c = 0; c < n_comms; ++c)
      n_masters[r][c] = 0;
  // we only add all entries, each process only has non-zeros in exactly one row
  MPI_Allreduce(n_masters_send, n_masters, size*n_comms, MPI_INT, MPI_SUM, 
                MPI_COMM_WORLD);
  
  std::vector<std::vector<size_t>> ret(size, std::vector<size_t>(n_comms, 0));
  for(size_t c = 0; c < n_comms; ++c)
  {
    for(size_t r = 0; (int)r < size; ++r)
    {
      ret[r][c] = n_masters[r][c];
    }
  }
  return ret;
}
#endif // _MPI

// *****************************************************************************
// return a vector whose size is the total number of (process) local dofs. For 
// each entry it holds the index of the master process of this dof.
#ifdef _MPI
std::vector<int> get_block_masters(
  std::vector<const TParFECommunicator3D*> comms)
{
  size_t n_comms = comms.size();
  // number of (process) local dofs, sum over all blocks
  size_t n_dofs_local = 0;
  for(size_t c = 0; c < n_comms; ++c)
    n_dofs_local += comms[c]->GetNDof();
  std::vector<int> block_masters(n_dofs_local);
  auto it = block_masters.begin();
  for(size_t c = 0; c < n_comms; ++c)
  {
    auto masters = comms[c]->GetMaster();
    auto n_dof = comms[c]->GetNDof(); // local within this process and block
    it = std::copy(masters, masters + n_dof, it);
  }
  return block_masters;
}
#endif // _MPI

// *****************************************************************************
// return a vector whose size is the total number of (process) local dofs. For
// each entry it holds the global index in a vector whose entries are ordered
// according to the processes.
// This is a vector of int (rather than size_t) because PETSc likes ints better.
#ifdef _MPI
std::vector<int> local_to_global_block(
  std::vector<const TParFECommunicator3D*> comms)
{
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  // the number of communicators is the number of blocks
  size_t n_comms = comms.size();
  auto n_masters = get_n_masters_per_block_and_process(comms);
  auto offset = [&n_masters, &n_comms](const int rank, const int block)
  {
    int offset = 0;
    for(int r = 0; r < rank; ++r)
    {
      offset -= n_masters.at(r).at(block);
    }
    for(int r = 0; r < rank; ++r)
    {
      for(size_t b = 0; b < n_comms; ++b)
      {
        offset += n_masters[r][b];
      }
    }
    for(int b = 0; b < block; ++b)
    {
      offset += n_masters[rank][b];
    }
    return offset;
  };
  
  // number of (process) local dofs, sum over all blocks
  size_t n_dofs_local = 0;
  for(size_t c = 0; c < n_comms; ++c)
    n_dofs_local += comms[c]->GetNDof();
  // the vector to be returned
  std::vector<int> ret(n_dofs_local, 0);
  auto it = ret.begin(); // iterator for the returned vector
  for(size_t c = 0; c < n_comms; ++c)
  {
    // this is the local_to_global array local within this process and within
    // this block
    auto local_to_global = comms[c]->Get_Local2Global();
    auto n_dof = comms[c]->GetNDof(); // local within this process and block
    std::copy(local_to_global, local_to_global + n_dof, it);
    // each of the copied entries must now be adjusted so that the PETSc 
    // ordering is accomplished
    auto masters = comms[c]->GetMaster();
    std::vector<int> offset_per_process(size);
    for(size_t p = 0; (int)p < size; ++p)
    {
      offset_per_process[p] = offset(p, c);
    }
    for(size_t i = 0; (int)i < n_dof; ++i, ++it)
    {
      *it += offset_per_process[masters[i]];
    }
  }
  return ret;
}
#endif // _MPI

// *****************************************************************************
void create_sub_matrix(const BlockFEMatrix& matrix, 
                       const std::pair<size_t, size_t>& start, 
                       const std::pair<size_t, size_t>& end, Mat& petsc_mat)
{
  // sanity checks
  if(start.first > end.first || start.second > end.second)
  {
    ErrThrow("extracting blocks from a BlockFEMatrix into a single PETSc "
             "matrix. Starting block is after ending block: (", start.first, 
             ",", start.second, ") and (", end.first, ",", end.second, ")");
  }
  if(end.first >= matrix.get_n_cell_rows() 
     || end.second >= matrix.get_n_cell_columns())
  {
    ErrThrow("extracting blocks from a BlockFEMatrix into a single PETSc "
             "matrix. Indicated subblocks are out of range (", start.first, 
             ",", start.second, ") and (", end.first, ",", end.second, ")");
  }
  
  auto sub_block = matrix.get_combined_submatrix(start, end);
  
#ifndef _MPI 
  // sequential case
  
  // 1) create matrix
  // create a matrix with known non-zero distribution among the rows
  size_t n_rows = sub_block->get_n_rows();
  size_t n_cols = sub_block->get_n_columns();
  // number of entries for each row
  std::vector<int> nnz(n_rows, 0);
  // loop over all rows to get the number of entries in each row
  for(size_t r = 0; r < n_rows; ++r)
  {
    nnz[r] = sub_block->get_n_entries_in_row(r);
  }
  MatCreateSeqAIJ(PETSC_COMM_WORLD, (int)n_rows, (int)n_cols, 0, &nnz[0], 
                  &petsc_mat);
  
  // 2) fill matrix
  // copy the entries
  const int * row_ptr = sub_block->get_row_ptr();
  const int * col_ptr = sub_block->get_vector_columns();
  const double * entries = sub_block->GetEntries();
  // loop over all rows, add entries for each row
  for(size_t r = 0; r < n_rows; ++r)
  {
    size_t n_entries_in_row = nnz[r];
    if(n_entries_in_row == 0)
      continue;
    int row_index = r; // conversion: size_t to int
    MatSetValues(petsc_mat, 1, &row_index, n_entries_in_row, 
                 &col_ptr[row_ptr[r]], entries+row_ptr[r], INSERT_VALUES);
  }
#else
  // MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  // collecting some useful vectors/arrays from the involved sub blocks
  size_t n_block_rows = end.first - start.first + 1;
  size_t n_block_cols = end.second - start.second + 1;
  // communicators for the row spaces (test spaces) of this matrix
  std::vector<const TParFECommunicator3D*> row_communicators(n_block_rows);
  // for PETSc only the master rows belong to a (process local) matrix
  size_t n_master_rows = 0;
  size_t n_global_rows = 0;
  for(size_t r = 0; r < n_block_rows; ++r)
  {
    auto comm = &matrix.get_row_space(r + start.first)->get_communicator();
    row_communicators[r] = comm;
    n_master_rows += comm->GetN_Master();
    n_global_rows += comm->get_n_global_dof();
  }
  // communicators for the column spaces (ansatz spaces) of this matrix
  std::vector<const TParFECommunicator3D*> col_communicators(n_block_cols);
  size_t n_master_cols = 0;
  size_t n_global_cols = 0;
  for(size_t c = 0; c < n_block_cols; ++c)
  {
    auto comm = &matrix.get_column_space(c + start.second)->get_communicator();
    col_communicators[c] = comm;
    n_master_cols += comm->GetN_Master();
    n_global_cols += comm->get_n_global_dof();
  }
  
  size_t n_rows = sub_block->get_n_rows();
  const auto row_ptr = sub_block->get_row_array();
  const auto col_ptr = sub_block->get_vector_columns();
  auto row_block_masters = get_block_masters(row_communicators);
  auto col_block_masters = get_block_masters(col_communicators);
  auto row_block_local2global = local_to_global_block(row_communicators);
  auto col_block_local2global = local_to_global_block(col_communicators);
  
  
  // 1) create matrix
  // compute the number of entries in each row 
  std::vector<int> nnz_diagonal(n_master_rows, 0);
  std::vector<int> nnz_off_diagonal(n_master_rows, 0);
  
  // fill vectors nnz_diagonal and nnz_off_diagonal
  
  // loop over the matrix rows, for each row we find the number of master 
  // entries (whose column index corresponds to a master dof)
  for(size_t row = 0, master_index = 0; row < n_rows; ++row)
  {
    if(row_block_masters[row] != my_rank)
      continue;
    // number of nonzero entries in a master column in this row
    size_t n_masters_in_row = 0;
    int begin_row = row_ptr[row];
    int end_row = row_ptr[row+1];
    // loop over all entries in this row
    for(size_t index = begin_row; (int)index < end_row; ++index)
    {
      // column index within the combined matrix
      size_t column = col_ptr[index];
      if(col_block_masters[column] == my_rank)
        n_masters_in_row++;
    }
    nnz_diagonal.at(master_index) = n_masters_in_row;
    // end_row - begin_row = number of entries in this row
    nnz_off_diagonal.at(master_index) = end_row - begin_row - n_masters_in_row;
    master_index++;
  }
  
  MatCreateAIJ(PETSC_COMM_WORLD, (int)n_master_rows, (int)n_master_cols, 
               (int)n_global_rows, (int)n_global_cols, 0, &nnz_diagonal[0], 0, 
               &nnz_off_diagonal[0], &petsc_mat);
  
  // 2) fill matrix
  // get raw data of the matrix
  const std::vector<double>& entries = sub_block->get_entries();
  for(size_t row = 0; row < n_rows; ++row)
  {
    if(row_block_masters[row] != my_rank)
      continue; // skip slave rows
    int global_row = row_block_local2global[row];
    int begin_row = row_ptr[row];
    int end_row = row_ptr[row+1];
    for(size_t col_index = begin_row; (int)col_index < end_row; ++col_index)
    {
      int col = col_ptr[col_index];
      double entry = entries[col_index];
      int global_col = col_block_local2global[col];
      MatSetValues(petsc_mat, 1, &global_row, 1, &global_col, &entry, 
                   INSERT_VALUES);
    }
  }
#endif
  
  // the following two functions must be called after MatSetValues (says PETSc 
  // documentation)
  MatAssemblyBegin(petsc_mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(petsc_mat, MAT_FINAL_ASSEMBLY);
}

// *****************************************************************************
void create_one_big_petsc_matrix(const BlockFEMatrix& matrix, Mat& petsc_mat)
{
  if(matrix.get_n_total_rows() != matrix.get_n_total_columns())
    ErrThrow("PETScSolver: non-square matrix");
  size_t n = matrix.get_n_cell_rows();
  size_t m = matrix.get_n_cell_columns();
  create_sub_matrix(matrix, {0,0}, {n-1,m-1}, petsc_mat);
}

// *****************************************************************************
PETScSolver::PETScSolver(const BlockFEMatrix& matrix,
                         const ParameterDatabase& db)
 : petsc_mat(nullptr)
#ifdef _MPI
  , comms_(matrix.get_communicators())
#endif // _MPI
{
  // we will create petsc_mat (type: Mat). All entries of 'matrix' will be
  // copied to it. That copying is really slow if petsc_mat has not been created
  // with a known sparsity structure. This would cause many memory allocations 
  // and copies. Instead all that PETSc needs to know at first is the number of 
  // non-zero entries in each row. This is therefore computed first here.
#if (defined _MPI && defined __2D__)
  static_assert(false, "you can not use MPI in 2D (yet)");
#endif
  // sanity check:
  if(matrix.get_n_total_rows() != matrix.get_n_total_columns())
    ErrThrow("PETScSolver: non-square matrix ", matrix.get_n_total_rows(),
             ",", matrix.get_n_total_columns());
  size_t n_block_rows = matrix.get_n_cell_rows();
  size_t n_block_cols = matrix.get_n_cell_columns();
  
  // indicate if matrix is stored as a transposed, this is set when calling
  // matrix.get_block(...)
  bool transposed;
  auto block = matrix.get_block(n_block_rows-1, n_block_cols-1, transposed);
  // indicate if we have a saddle point problem
  // this is considered a saddle point problem if the last block is all zero
  /// @todo save this information in the BlockMatrix
  bool saddlepoint_problem = (block->GetNorm() == 0.0);
  
  // Is pc_type LU, ILU or CHOLESKY? I.e. do we want to use a direct solver
  bool use_direct_petsc = false;
  
  // call PetscInitialize
  {
    // the following is necessary to properly call PetscInitialize
    std::string petsc_args = addParameters2str(db);
    if(saddlepoint_problem)
    {
      // PETSc can figure out which block is zero. Also this is the only way to
      // use specialized saddle point solvers
      petsc_args.append(" -pc_fieldsplit_detect_saddle_point");
    }
    auto char_vector = str_to_vector_char_p(petsc_args);
    char ** params_pointer = char_vector.data();
    int argc = char_vector.size();

    PetscInitialize(&argc, &params_pointer, (char*)0, nullptr);

    // properly delete all the char* which where created
    // in str_to_vector_char_p(...)
    for(auto cp : char_vector)
      delete [] cp;
    
    auto npos = std::string::npos;
    if ( petsc_args.find("-pc_type lu", 0) != npos
      || petsc_args.find("-pc_type ilu", 0) != npos
      || petsc_args.find("-pc_type cholesky", 0) != npos
      || petsc_args.find("-pc_type hypre", 0) != npos)
    {
      use_direct_petsc = true;
    }
  }

  bool single_block = (n_block_rows*n_block_cols == 1);
  if(use_direct_petsc || single_block)
  {
    // the nested petsc matrix type will not work, so we write a single big 
    // petsc matrix. This is similar to calling matrix.get_combined_matrix
    create_one_big_petsc_matrix(matrix, this->petsc_mat);
  }
  else if(saddlepoint_problem && !single_block)
  {
    Output::print("PETScSolver class saddlepoint_problem");
    // reserve space for enough sub matrices for PETSc. You have to have a 
    // two-by-two block system in a saddle point problem.
    // It is assumed here, that the pressure block is last, i.e., has the row 
    // and colums index `n_block_rows-1` and `n_block_cols-1` respectively.
    size_t last_block_row = n_block_rows-1;
    size_t last_block_col = n_block_cols-1;
    sub_petsc_mats.resize(2 * 2);
    
    // velocity-velocity coupling:
    create_sub_matrix(matrix, {0,0}, {last_block_row-1, last_block_col-1},
                      sub_petsc_mats[0]);
    // velocity-pressure coupling:
    create_sub_matrix(matrix, {0,last_block_col}, {last_block_row-1, last_block_col},
                      sub_petsc_mats[1]);
    // pressure-velocity coupling:
    create_sub_matrix(matrix, {last_block_row, 0}, {last_block_row, last_block_col-1},
                      sub_petsc_mats[2]);
    // pressure-pressure coupling (should be a zero matrix)
      create_sub_matrix(matrix, {last_block_row, last_block_col},
                        {last_block_row, last_block_col},
                        sub_petsc_mats.at(3));
    
    // currently the last block should be all zero, we check this here
    PetscReal norm = 0.0;
    MatNorm(sub_petsc_mats.at(3), NORM_FROBENIUS, &norm);
    if(norm != 0.0)
    {
      // We are here usually in case of an all-Dirichlet Navier-Stokes problem,
      // where the pressure-pressure block is modified in one row. usually in 
      // saddle point problems this is zero. This is what PETSc expects.
      // The problem is that PETSc can not automatically identify the correct 
      // blocks if they are not the ones given as submatrices. We would have to
      // somehow define index sets for PETSc. I don't know how to do that.
      ErrThrow("Currently one can use PETSc for a saddle point problem only if "
               "the lower right block is zero. Here it has norm ", norm);
    }
    
    // This is the PETSc way of having a BlockMatrix, direct solvers wont work
    MatCreateNest(PETSC_COMM_WORLD, 2, nullptr, 2, nullptr, &sub_petsc_mats[0],
                  &petsc_mat);
  }
  else
  {
    create_one_big_petsc_matrix(matrix, this->petsc_mat);
  }
  
  // create solver and preconditioner objects
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, petsc_mat, petsc_mat);
  
  // PETSc preconditioner object
  PC  pc;
  // set preconditioner context (PC)
  KSPGetPC(ksp, &pc);
  
  // set the options from the string passed to PetscInitialize which includes
  // the values in the parameter "petsc_arguments".
  KSPSetFromOptions(ksp);
  PCSetFromOptions(pc);
  KSPSetUp(ksp);
  // some PETSc information about the solver
  /// @todo nicer output for the PETSc solver
  if(Output::getVerbosity() > 3)
    KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD);
}

/* ************************************************************************** */
PETScSolver::PETScSolver(const BlockMatrix&, const ParameterDatabase&)
{
  ErrThrow("can not construct a PETScSolver using a BlockMatrix. I need a "
           "BlockFEMatrix instead. It would be possible to implement this "
           "here, in case it is really needed.");
}

PETScSolver::PETScSolver(const CompositeOperator<BlockVector>&,
      const ParameterDatabase&)
{
  ErrThrow("Cannot use PETScSolver for a composite operator!");
}

/* ************************************************************************** */
PETScSolver::PETScSolver(PETScSolver && other)
{
  // the type of this->petsc_mat is 'Mat' which is in fact a pointer
  // destroy the old petsc_mat
  MatDestroy(&petsc_mat);
  // copy the pointer
  this->petsc_mat = other.petsc_mat;
  // reset to some zero matrix in the other object
  int n_rows = 0;
  int n_cols = 0;
  int n_nonzero_per_row = 0;
  MatCreateSeqAIJ(PETSC_COMM_WORLD, n_rows, n_cols, n_nonzero_per_row, nullptr,
                  &other.petsc_mat);
}

/* ************************************************************************** */
PETScSolver& PETScSolver::operator=(PETScSolver && other)
{
  // the type of this->petsc_mat is 'Mat' which is in fact a pointer
  // destroy the old petsc_mat
  MatDestroy(&petsc_mat);
  // copy the pointer
  this->petsc_mat = other.petsc_mat;
  // reset to some zero matrix in the other object
  int n_rows = 0;
  int n_cols = 0;
  int n_nonzero_per_row = 0;
  MatCreateSeqAIJ(PETSC_COMM_WORLD, n_rows, n_cols, n_nonzero_per_row, nullptr,
                  &other.petsc_mat);
  return *this;
}

/* ************************************************************************** */
PETScSolver::~PETScSolver()
{
  KSPDestroy(&ksp);
  // if we have just one sub matrix
  // it is the petsc_mat
  if (sub_petsc_mats.size() > 1)
  {
    for (size_t i = 0; i < sub_petsc_mats.size(); ++i)
    {
      MatDestroy(&sub_petsc_mats[i]);
    }
  }
  MatDestroy(&petsc_mat);
  // remember you have to call 'PetscFinalize();' at the very end of the program
}

/* ************************************************************************** */
void PETScSolver::solve(const BlockVector& rhs, BlockVector& solution)
{
  /// @todo check if matrix sizes are ok with solution and rhs sizes
  /// MatGetSize() does not work for MatNest mat type, which is used for
  /// saddle point problems.
  Output::info<5>("PETScSolver::solve");
  size_t n_local = solution.length();
  if(n_local != rhs.length())
  {
    ErrThrow("PETScSolver::solve: size of the rhs and the solution are not "
             "equal, ", rhs.length(), " != ", n_local);
  }
  
  // length of vector
  int n_global = n_local;
  int n_local_masters = n_local;
#ifdef _MPI
  n_local_masters = 0;
  n_global = 0;
  for(size_t bl = 0; bl < comms_.size(); ++bl) // loop over all blocks
  {
    n_local_masters += comms_[bl]->GetN_Master();
    n_global += comms_[bl]->get_n_global_dof();
  }
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
  auto block_local2global = local_to_global_block(comms_);
  auto block_masters = get_block_masters(comms_);
  if(block_local2global.size() != n_local)
  {
    ErrThrow("block_local2global has wrong size ", block_local2global.size(),
             " != ", n_local);
  }
  if(block_masters.size() != n_local)
  {
    ErrThrow("block_masters has wrong size ", block_masters.size(), " != ", 
             n_local);
  }
#endif
  
  
  // create PETSc vectors
  Vec x, b;
  VecCreate(PETSC_COMM_WORLD, &x);
  PetscObjectSetName((PetscObject) x, "Solution");
  VecSetSizes(x, n_local_masters, n_global);
  VecSetFromOptions(x);
  VecDuplicate(x, &b);
  PetscObjectSetName((PetscObject) b, "Rhs");
  
  // copy ParMooN BlockVectors to PETSc vectors
  for(size_t i = 0; i < n_local; ++i)
  {
    int index = i; // conversion to int
#ifdef _MPI
    index = block_local2global[i];
    if(block_masters[i] == my_rank) // only consider master dofs
#endif // _MPI
    {
      VecSetValue(b, index, rhs[i], INSERT_VALUES);
      VecSetValue(x, index, solution[i], INSERT_VALUES);
    }
  }
  
  VecAssemblyBegin(x);
  VecAssemblyEnd(x);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);
  
  // solve the linear system using PETSc
  KSPSolve(ksp, b, x);
  
  PetscInt its;
  KSPGetIterationNumber(ksp, &its);
#ifdef _MPI
  if(my_rank == 0)
#endif
  Output::print(prefix, " PETSc solver: number of iterations: ", its);
  
  
  // copy back to solution:
  for(size_t i = 0; i < n_local; ++i)
  {
    int index = i; // conversion to int
#ifdef _MPI
    index = block_local2global[i];
    if(block_masters[i] == my_rank) // only consider master dofs
#endif //_MPI
    {
      VecGetValues(x, 1, &index, &solution[i]);
    }
  }
  // delete petsc vectors
  VecDestroy(&b);
  VecDestroy(&x);
}

#else // !PARMOON_WITH_PETSC
PETScSolver::PETScSolver(const BlockFEMatrix&, const ParameterDatabase&)
{
  ErrThrow("ParMooN has been compiled without PETSc, therefore a PETScSolver "
           "can not be created.");
}
PETScSolver::PETScSolver(const BlockMatrix&, const ParameterDatabase&)
{
  ErrThrow("ParMooN has been compiled without PETSc, therefore a PETScSolver "
           "can not be created.");
}
PETScSolver::~PETScSolver()
{
  
}
void PETScSolver::solve(const BlockVector&, BlockVector&)
{
  ErrThrow("ParMooN has been compiled without PETSc, therefore "
           "PETScSolver::solve can not be called.");
}
#endif // PARMOON_WITH_PETSC
