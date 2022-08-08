// =======================================================================
// @(#)DirectSolver.h
//
// Purpose:     solve equation system by direct solver
//
// Author:      Gunar Matthies (06.09.05)
//
// History:     start of implementation 06.09.05 (Gunar Matthies)
//
// =======================================================================

#include "all_defines_external_libraries.h"
#include <DirectSolver.h>
#ifdef PARMOON_WITH_UMFPACK
#include "umfpack.h"
#endif
#include <BlockFEMatrix.h>
#include <BlockVector.h>
#include <Matrix.h>
#include "MooNMD_Io.h"
#include <string.h>

#ifdef PARMOON_WITH_PARDISO
#include <omp.h>
/* Pardiso Fortran Function */
extern "C" void pardiso_(
    void * handle,
    int *max_factorizations,
    int *matrix_num,
    int *matrix_type,
    int *ido,
    int *neqns,
    double *a,
    int *ia,
    int *ja,
    int *perm_user,
    int *nb,
    int iparam[pardiso_options_array_length],
    int *msglvl,
    double *b,
    double *x,
    int *error
);
#else
void pardiso_(void * , int *, int *, int *, int *, int *, double *, int *,
              int *, int *, int *, int [pardiso_options_array_length], int *,
              double *, double *, int *)
{
  ErrThrow("you compiled without pardiso, therefore a call to 'pardiso_' does "
           "not work.");
}
#endif


void handle_error_umfpack(int ierror)
{
#ifdef PARMOON_WITH_UMFPACK
  if (ierror == UMFPACK_OK)
    return;

  switch(ierror)
  {
    //WARNINGS
    case UMFPACK_WARNING_singular_matrix:
      Output::warn("UMFPACK", "Matrix is singular!");
      break;
    case UMFPACK_WARNING_determinant_underflow:
      Output::warn("UMFPACK", "Determinant smaller than eps");
      break;
    case UMFPACK_WARNING_determinant_overflow:
      Output::warn("UMFPACK", "Determinant is larger than IEEE Inf");
      break;
    //ERRORS
    case UMFPACK_ERROR_out_of_memory:
      ErrThrow("umfpack: Out of Memory");
      break;
    case UMFPACK_ERROR_invalid_Numeric_object:
      ErrThrow("umfpack: Invalid numeric factorization object");
      break;
    case UMFPACK_ERROR_invalid_Symbolic_object:
      ErrThrow("umfpack: Invalid symbolic factorization object");
      break;
    case UMFPACK_ERROR_argument_missing:
      ErrThrow("umfpack: Argument Missing.");
      break;
    case UMFPACK_ERROR_n_nonpositive:
      ErrThrow("umfpack: Matrix dimensions not positive.");
      break;
    case UMFPACK_ERROR_invalid_matrix:
      ErrThrow("umfpack: Invalid Matrix Structure.");
      break;
    case UMFPACK_ERROR_different_pattern:
      ErrThrow("umfpack: Different sparse pattern.");
      break;
    case UMFPACK_ERROR_invalid_system:
      ErrThrow("umfpack: Invalid system provided with sys.");
      break;
    case UMFPACK_ERROR_invalid_permutation:
      ErrThrow("umfpack: Invalid permutation vector.");
      break;
    case UMFPACK_ERROR_file_IO:
      ErrThrow("umfpack: Fille IO error.");
      break;
    case UMFPACK_ERROR_internal_error:
      ErrThrow("umfpack: Internal error.");
      break;
    
    default:
      ErrThrow("umfpack: unkown error. Error number ", ierror); break;
    break;
  }
#endif // PARMOON_WITH_UMFPACK
}
void handle_error_pardiso(int ierror)
{
  if (ierror==0) return;

  switch(ierror)
  {
  case -1:
    ErrThrow("Input inconsistent."); break;
  case -2:
    ErrThrow("Not enough memory."); break;
  case -3:
    ErrThrow("Reordering problem."); break;
  case -4:
    ErrThrow("Zero pivot, numerical factorization or iterative refinement problem."); break;
  case -5:
    ErrThrow("Unclassified (internal) error."); break;
  case -6:
    ErrThrow("Preordering failed"); break;
  case -7:
    ErrThrow("Diagonal matrix problem."); break;
  case -8:
    ErrThrow("32-bit integer overflow problem."); break;
  case -10:
    ErrThrow("No license file pardiso.lic found."); break;
  case -11:
    ErrThrow("License is expired."); break;
  case -12:
    ErrThrow("Wrong user name or host name."); break;
  case -100:
    ErrThrow("Reached maximum number of Krylov-subspace iteration in iterative solver."); break;
  case -101:
    ErrThrow("No sufficient convergence in Krylov-subspace iteration within 25 iterations"); break;
  case -102:
    ErrThrow("Error in Krylov-subspace iteration."); break;
  case -103:
    ErrThrow("Break-Down in Krylov-subspace iteration."); break;
  default:
    ErrThrow("pardiso: unkown error. Error number ", ierror); break;
  break;
  }
}

/** ************************************************************************ */
DirectSolver::DirectSolver(std::shared_ptr<TMatrix> matrix, 
                           DirectSolver::DirectSolverTypes type)
 : type(type), matrix(matrix), cols(), rows(),
   symbolic(nullptr), numeric(nullptr), pt(), maxfct(10), mnum(1), 
   mtype(11), perm(0), nrhs(1), iparm(), msglvl(0)
{
#ifndef PARMOON_WITH_UMFPACK
  if(type == DirectSolverTypes::umfpack)
    ErrThrow("ParMooN has been configured without umfpack, therefore a "
             "DirectSolver object with type 'umfpack' can not be created.");
#endif // not PARMOON_WITH_UMFPACK
#ifndef PARMOON_WITH_PARDISO
  if(type == DirectSolverTypes::pardiso)
    ErrThrow("ParMooN has been configured without pardiso, therefore a "
             "DirectSolver object with type 'pardiso' can not be created.");
#endif
  
  Output::print<5>("constructing a DirectSolver object");
  if(!matrix->is_square())
  {
    ErrThrow("unable to factorize a non-square matrix ", matrix->get_n_rows(),
             "  ", matrix->get_n_columns());
  }
  
  if(type == DirectSolverTypes::pardiso)
  {
    this->matrix->fortran_shift();
    // initialize all values to zero/nullptr
    for(size_t i = 0; i < pardiso_options_array_length; ++i)
    {
      this->pt[i] = nullptr;
      this->iparm[i] = 0;
    }
    this->iparm[0] = 0; /* override defaults */
    this->iparm[1] = 2; /* Metis reordering, default: minimize fill-in; 3D */
#ifdef _OMP
    Output::print<2>("FYI: omp_get_max_threads when setting up ParDiso is ",
                     omp_get_max_threads());

    this->iparm[2] = omp_get_max_threads(); //number of threads set to OMP_NUM_THREADS as recommended in the doc
#endif
    this->iparm[3] = 0; /* precond cgs, solver type */
    this->iparm[4] = 0; /* user permute */
  }

  // the threshold is rather small here, it should furthermore depend on the 
  // dimension (2 or 3) and the polynomial degree (and possibly more).
  if(this->matrix->get_n_rows() > 1e4 && type == DirectSolverTypes::umfpack)
  {
    this->cols.resize(this->matrix->get_n_entries(), 0);
    this->rows.resize(this->matrix->get_n_rows()+1, 0);
    for(int i = 0; i < this->matrix->get_n_entries(); ++i)
      this->cols[i] = this->matrix->get_vector_columns()[i];
    for(int i = 0; i < this->matrix->get_n_rows()+1; ++i)
      this->rows[i] = this->matrix->get_row_ptr()[i];
  }
  
  this->symbolic_factorize();
  this->numeric_factorize();
}

/** ************************************************************************ */
DirectSolver::DirectSolver(const BlockMatrix& matrix, 
                           DirectSolver::DirectSolverTypes type)
 : DirectSolver(matrix.get_combined_matrix(), type)
{
}

/** ************************************************************************ */
DirectSolver::DirectSolver(const BlockFEMatrix& matrix, 
                           DirectSolver::DirectSolverTypes type)
 : DirectSolver(matrix.get_combined_matrix(), type)
{
}

/** ************************************************************************ */
DirectSolver::DirectSolver(const TMatrix& matrix,
                           DirectSolver::DirectSolverTypes type)
 : DirectSolver(std::make_shared<TMatrix>(matrix), type)
{
}

DirectSolver::DirectSolver(const CompositeOperator<BlockVector>&,
      DirectSolverTypes)
{
  ErrThrow("Cannot use DirectSolver for a composite operator!");
}


/** ************************************************************************ */
DirectSolver::DirectSolver(DirectSolver&& other)
 : type(other.type), matrix(other.matrix),
   cols(std::move(other.cols)), rows(std::move(other.rows)),
   symbolic(other.symbolic), numeric(other.numeric)
{
  other.symbolic = nullptr;
  other.numeric = nullptr;
  Output::print<5>("DirectSolver::DirectSolver(DirectSolver&&)");
}

/** ************************************************************************ */
class DirectSolver& DirectSolver::operator=(DirectSolver&& other)
{
  this->type = other.type;
  this->matrix = other.matrix;
  this->symbolic = other.symbolic;
  this->numeric = other.numeric;
  this->cols = std::move(other.cols);
  this->rows = std::move(other.rows);
  other.symbolic = nullptr;
  other.numeric = nullptr;
  Output::print<5>("DirectSolver::operator=(DirectSolver&&)");
  return *this;
}

/** ************************************************************************ */
DirectSolver::~DirectSolver()
{
  switch(type)
  {
    case DirectSolver::DirectSolverTypes::umfpack:
#ifdef PARMOON_WITH_UMFPACK
      if(this->cols.size() == 0)
      {
        // using int for indices
        umfpack_di_free_symbolic(&symbolic);
        umfpack_di_free_numeric(&numeric);
      }
      else
      {
        // using long for indices
        umfpack_dl_free_symbolic(&symbolic);
        umfpack_dl_free_numeric(&numeric);
      }
#endif // PARMOON_WITH_UMFPACK
      break;
    case DirectSolver::DirectSolverTypes::pardiso:
    {
      int phase = -1; /* Release internal memory. */
      int n_eq = this->matrix->get_n_rows();
      double * entries = this->matrix->GetEntries();
      int * rows = this->matrix->get_row_ptr();
      int * cols = this->matrix->get_vector_columns();
      double dzero; // ??
      int error;
      pardiso_(this->pt, &this->maxfct, &this->mnum, &this->mtype, &phase,
               &n_eq, entries , rows, cols, &this->perm, &this->nrhs,
               this->iparm, &this->msglvl, &dzero, &dzero, &error);
      handle_error_pardiso(error);
      this->matrix->fortran_shift();
      break;
    }
    default:
      ErrThrow("unknown DirectSolverTypes ",
               static_cast<typename 
                 std::underlying_type<DirectSolverTypes>::type> (type));
  }
  Output::print<5>("destructed a DirectSolver object");
}

/** ************************************************************************ */
void DirectSolver::symbolic_factorize()
{
  int n_eq = matrix->get_n_rows();
  switch(type)
  {
    case DirectSolverTypes::umfpack:
    {
#ifdef PARMOON_WITH_UMFPACK
      // symbolic factorization
      if(this->cols.size() == 0)
      {
        // using int for indices
        int error = umfpack_di_symbolic(n_eq, n_eq, matrix->get_row_ptr(),
                                        matrix->get_vector_columns(), matrix->GetEntries(),
                                        &symbolic, nullptr, nullptr);
        handle_error_umfpack(error);
      }
      else
      {
        // using long for indices
        int error = umfpack_dl_symbolic(n_eq, n_eq, &this->rows[0], 
                                        &this->cols[0], matrix->GetEntries(), 
                                        &symbolic, nullptr, nullptr);
        handle_error_umfpack(error);
      }
#endif // PARMOON_WITH_UMFPACK
      break;
    }
    case DirectSolverTypes::pardiso:
    {
      int phase = 11; // analysing phase
      int n_eq = this->matrix->get_n_rows();
      double * entries = this->matrix->GetEntries();
      int * rows = this->matrix->get_row_ptr();
      int * cols = this->matrix->get_vector_columns();
      double dzero; // ??
      int error;
      pardiso_(this->pt, &this->maxfct, &this->mnum, &this->mtype, &phase,
               &n_eq, entries , rows, cols, &this->perm, &this->nrhs,
               this->iparm, &this->msglvl, &dzero, &dzero, &error);
      handle_error_pardiso(error);
      break;
    }
    default:
      ErrThrow("unknown DirectSolverTypes ",
               static_cast<typename 
                 std::underlying_type<DirectSolverTypes>::type> (type));
      break;
  }
}

/** ************************************************************************ */
void DirectSolver::numeric_factorize()
{
  switch(type)
  {
    case DirectSolverTypes::umfpack:
    {
#ifdef PARMOON_WITH_UMFPACK
      double Info[UMFPACK_INFO];
      double Control[UMFPACK_CONTROL];
      umfpack_di_defaults(Control);
      
      if(this->cols.size() == 0)
      {
        // using int for indices
        int error = umfpack_di_numeric(matrix->get_row_ptr(), matrix->get_vector_columns(),
                                       matrix->GetEntries(), symbolic, &numeric,
                                       Control, Info);
        handle_error_umfpack(error);
      }
      else
      {
        // using long for indices
        int error = umfpack_dl_numeric(&this->rows[0], &this->cols[0],
                                       matrix->GetEntries(), symbolic, &numeric,
                                       Control, Info);
        handle_error_umfpack(error);
      }
      Output::print<4>("umfpack: Peak memory ", Info[UMFPACK_PEAK_MEMORY],
                       ", estimated condition number ", 1./Info[UMFPACK_RCOND]);
#endif // PARMOON_WITH_UMFPACK
      break;
    }
    case DirectSolverTypes::pardiso:
    {
      int phase = 22; // numerical factorization
      int n_eq = this->matrix->get_n_rows();
      double * entries = this->matrix->GetEntries();
      int * rows = this->matrix->get_row_ptr();
      int * cols = this->matrix->get_vector_columns();
      double dzero; // ??
      int error;
      pardiso_(this->pt, &this->maxfct, &this->mnum, &this->mtype, &phase,
               &n_eq, entries, rows, cols, &this->perm, &this->nrhs,
               this->iparm, &this->msglvl, &dzero, &dzero, &error);
      handle_error_pardiso(error);
      break;
    }
    default:
      ErrThrow("unknown DirectSolverTypes ",
               static_cast<typename 
                 std::underlying_type<DirectSolverTypes>::type> (type));
      break;
  }
}

/** ************************************************************************ */
void DirectSolver::solve(const double* rhs, double* solution)
{
  Output::print<5>("solving using a direct solver");
  switch(type)
  {
    case DirectSolverTypes::umfpack:
    {
#ifdef PARMOON_WITH_UMFPACK
      // symbolic factorization
      if(this->cols.size() == 0)
      {
        // using int for indices
        int error = umfpack_di_solve(UMFPACK_At, matrix->get_row_ptr(), 
                                     matrix->get_vector_columns(), matrix->GetEntries(),
                                     solution, rhs, numeric, nullptr, nullptr);
        handle_error_umfpack(error);
      }
      else
      {
        // using long for indices
        int error = umfpack_dl_solve(UMFPACK_At, &this->rows[0], 
                                     &this->cols[0], matrix->GetEntries(),  
                                     solution, rhs, numeric, nullptr, nullptr);
        handle_error_umfpack(error);
      }
#endif // PARMOON_WITH_UMFPACK
      break;
    }
    case DirectSolverTypes::pardiso:
    {
      //make a copy of rhs since pardiso wants it non-const
      size_t length_rhs = matrix->get_n_rows(); //this was checked in the calling method
      double* rhs_nonconst = new double[length_rhs];
      memcpy(rhs_nonconst, rhs, length_rhs*sizeof(double));

      int phase = 33; //PHASE: solve
      int n_eq = this->matrix->get_n_rows();
      double * entries = this->matrix->GetEntries();
      int * rows = this->matrix->get_row_ptr();
      int * cols = this->matrix->get_vector_columns();
      int error;
      pardiso_(this->pt, &this->maxfct, &this->mnum, &this->mtype, &phase,
               &n_eq, entries, rows, cols, &this->perm, &this->nrhs,
               this->iparm, &this->msglvl, rhs_nonconst, solution, &error);
      handle_error_pardiso(error);
      delete [] rhs_nonconst;
      break;
    }
    default:
      ErrThrow("unknown DirectSolverTypes ",
               static_cast<typename 
                 std::underlying_type<DirectSolverTypes>::type> (type));
      break;
  }
}

/** ************************************************************************ */
void DirectSolver::solve(const BlockVector& rhs, BlockVector& solution)
{
  if(  (int) rhs.length() != this->matrix->get_n_rows()
    || (int) solution.length() != this->matrix->get_n_columns())
    ErrThrow("solution or right hand side vector has wrong size. ",
             "Size of the matrix: ", this->matrix->get_n_rows(), " x ", 
             this->matrix->get_n_columns(),"\t rhs size: ", rhs.length(), 
             "\tsolution size: ", solution.length());


  solve(rhs.get_entries(), solution.get_entries());

}
