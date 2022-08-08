// =======================================================================
// @(#)Matrix.C        1.2 11/20/98
//
// Class:       TMatrix
//
// Purpose:     store a  matrix (ansatz != test space)
//
// Author:      Gunar Matthies
//
// History:     26.08.1998 start implementation
//
// =======================================================================

#include <Matrix.h>
#include <string.h>
#include <LinAlg.h>

#include <MooNMD_Io.h>

#ifdef _MPI
#include<ParFECommunicator3D.h>
#endif

#include <fstream>
#include <algorithm>
#include <cmath>

TMatrix::TMatrix(std::shared_ptr<TStructure> structure
#ifdef _MPI
            , const SparsityType& sparse_type_in
#endif
)
 : structure(structure), entries(this->structure->get_n_entries(), 0.)
#ifdef _MPI
            ,sparse_type(sparse_type_in)
#endif
{
}

TMatrix::TMatrix(int nRows, int nCols)
 : TMatrix(std::make_shared<TStructure>(nRows,nCols))
{
  ;
}
/* *********************************************************** */

TMatrix::TMatrix(const std::vector<double>& diag,
                 std::shared_ptr<TStructure> Structure)
 : TMatrix(Structure)
{
	this->reset();

  // loop over all diagonal entries
  for(size_t i = 0; i < diag.size(); i++)
  {
	  this->set(i, i, diag[i]);
  }
}

/* ****************************************************************** */
void TMatrix::reset()
{
  memset(this->GetEntries(), 0., this->structure->get_n_entries()*sizeof(double));
}

int TMatrix::clean_denormals()
{
  int c = 0;

  size_t n = entries.size();
  for (size_t i = 0; i < n; i++)
  {
    if (!std::isnormal(entries[i]) && entries[i] != 0.0)
    {
      ++c;
      entries[i] = 0.0;
    }
  }

#ifdef _MPI
  MPI_Allreduce(MPI_IN_PLACE, &c, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

  return c;
}

void TMatrix::setEntries(const std::vector<double>& entries)
{
  if(this->entries.size() != entries.size())
  {
    ErrThrow("resetting the entries of a matrix with the wrong number of ",
             "entries ", this->entries.size(), " != ", entries.size());
  }
  this->entries = entries;
}

void TMatrix::write(const std::string& filename) const
{
  std::ofstream matrixfile;
  matrixfile.open(filename.c_str());

  //write the header line - coordinate format, real values, no symmetry used
  matrixfile << "%%MatrixMarket matrix coordinate real general \n";

  //write general matrix information
  matrixfile << get_n_rows() << "\t" << get_n_columns() << "\t" << get_n_entries() << "\n";

  //loop through matrix and print row - column - entry
  // for each entry in the sparsity structure
  const int* RowPtr = structure->get_row_ptr();
  const int* KCol = structure->get_vector_columns();
  int begin, end, pos=0;

  for (unsigned int i=0; i < structure->get_n_rows(); ++i)
  {
    begin = RowPtr[i];
    end   = RowPtr[i+1];

    for (int j=begin; j<end; ++j)
    {
      // shift row and col by +1 (fortran style)
      matrixfile << setprecision(16) << i +1 << "\t" << KCol[pos] + 1 << "\t" << entries[pos] << "\n";
      ++pos;
    }
  }


  matrixfile.close();
}

void TMatrix::Print(const char *name) const
{
  const int* RowPtr = structure->get_row_ptr();
  const int* KCol = structure->get_vector_columns();
  int begin, end, pos=0;

  for (unsigned int i=0; i < structure->get_n_rows(); ++i)
  {
    begin = RowPtr[i];
    end   = RowPtr[i+1];

    for (int j=begin; j<end; ++j)
    {
      Output::print<1>(name, "(", i, ",", KCol[pos], ") = ",
		       entries[pos],";");
      ++pos;
    }
  }
}

void TMatrix::PrintFull(const std::string& name, int fieldWidth) const
{
  const int* rowPtr = structure->get_row_ptr();
  const int* KCol = structure->get_vector_columns();

  cout << endl << name << " = " << endl;
  for (unsigned int curRow = 0; curRow < structure->get_n_rows(); curRow++)
  {
    int rowEnd = rowPtr[curRow+1];
    int posKCol = rowPtr[curRow];
    for (unsigned int curCol = 0; curCol < structure->get_n_columns(); curCol++)
    {
      if(static_cast<int>(curCol) == KCol[posKCol] && posKCol < rowEnd)
      {
        cout << setw(fieldWidth) << entries[posKCol] << ", ";
        posKCol++;
      }
      else
      {
        cout << setw(fieldWidth) << 0.0 << ", ";
      }
    }
    cout << endl;
  }
  cout << endl;
}

// add val to a matrix element
// return an error if the entry is not in the sparse structure
void TMatrix::add(int i,int j, double val)
{
  if(val != 0.0)
    this->get(i, j) += val;
}

void TMatrix::add(int i, const std::map<int,double>& vals, double factor)
{
  if(i < 0 || i > this->get_n_rows())
  {
    ErrThrow("This matrix does not have a row ", i,
             ".\nThe dimension of this matrix is ", this->get_n_rows(), " x ",
             this->get_n_columns());
  }
  const int* RowPtr = structure->get_row_ptr();
  const int* KCol = structure->get_vector_columns();
  auto it = vals.begin();
  for (int m=RowPtr[i];m < RowPtr[i+1] && it != vals.end(); m++)
  {
    if (KCol[m] == it->first)
    {
      entries[m] += factor*it->second;
      ++it;
    }
  }
  if(it != vals.end())
  {
    ErrThrow("Error in TMatrix::add. There are entries in 'vals' which are ",
             "not in the sparse structure. row ", i, ", column ", it->first,
             ".\nExit\n");
  }
}

void TMatrix::add(const std::map<int, std::map<int,double>>& vals,
                  double factor)
{
  // add every row to the matrix
  for(auto it = vals.begin(); it != vals.end(); ++it)
    add(it->first, it->second, factor);
}



// set val of a matrix element
// return an error if the entry is not in the sparse structure
void TMatrix::set(int i,int j, double val)
{
  int get_index_of_entry = this->structure->get_index_of_entry(i, j);
  if(get_index_of_entry == -1)
  {
    ErrThrow("Entry ", i,",",j," not in sparsity structure.");
  }
  this->entries[get_index_of_entry] = val;
}

// get val of a matrix element
// return an error if the entry is not in the sparse structure
const double& TMatrix::get(int i,int j) const
{
  // index of the entry (i,j) within the vector Entries
  int index = this->structure->get_index_of_entry(i, j);
  if(index >= 0 )
    return this->entries[index];
  ErrThrow("could not find the entry (", i, ",", j,
           ") in the sparsity structure");
}

double& TMatrix::get(int i,int j)
{
  // index of the entry (i,j) within the vector Entries
  int index = this->structure->get_index_of_entry(i, j);
  if(index >= 0 )
    return this->entries[index];
  ErrThrow("could not find the entry (", i, ",", j,
           ") in the sparsity structure");
}


std::vector<double> TMatrix::get_matrix_column(size_t j) const
{
  if((int) j > this->get_n_columns())
    ErrThrow("Requested column number greater than number of columns.");

  std::vector<double> col(get_n_rows(), 0.0);
  for(int i =0; i < get_n_rows();++i)
  {
    for(int l = get_row_ptr()[i]; l < get_row_ptr()[i+1] ; ++l)
    {
      if(get_vector_columns()[l] == (int) j)
      {
        col.at(i) = GetEntries()[l];
      }
    }
  }
  return col;
}

std::vector<double> TMatrix::get_matrix_row(size_t i) const
{
  if((int) i > this->get_n_rows())
    ErrThrow("Requested row number greater than number of rows.");

  std::vector<double> row(get_n_columns() , 0.0);
    for(int l = get_row_ptr()[i]; l < get_row_ptr()[i+1] ; ++l)
    {
      int j = get_vector_columns()[l];
      row.at(j) = GetEntries()[l];
    }
  return row;
}

void TMatrix::set_matrix_column(size_t j, const std::vector<double>& col)
{
  if((int) col.size() != get_n_rows())
    ErrThrow("Wrong number of entries in the given matrix column!");

  for(int i =0; i < get_n_rows();++i)
  {
    for(int l = get_row_ptr()[i]; l < get_row_ptr()[i+1] ; ++l)
    {
      if(get_vector_columns()[l] == (int) j)
      {
        set(i,j,col.at(i));
      }
    }
  }
}

double & TMatrix::operator()(const int i, const int j)
{
  return this->get(i, j);
}

const double & TMatrix::operator()(const int i, const int j) const
{
  return this->get(i, j);
}

std::vector<double> TMatrix::get_diagonal() const
{
  size_t n_diag = std::min(this->get_n_rows(), this->get_n_columns());
  std::vector<double> ret(n_diag, 0.0);
  // loop over all diagonal entries
  for(size_t d = 0; d < n_diag; ++d)
  {
    // this is basically this->get but we write a 0 into the vector in case the
    // diagonal is not in the sparsity structure
    int index = this->structure->get_index_of_entry(d, d);
    if(index >= 0 )
      ret[d] = this->entries[index];
    // else just leave the value to be zero
  }
  return ret;
}

std::vector<double> TMatrix::get_row_sums() const
{
  size_t n_rows = get_n_rows();
  std::vector<double> ret(n_rows, 0.0);

  for (size_t i = 0; i < n_rows; i++)
  {
    size_t n_in_row = structure->get_n_entries_in_row(i);
    size_t row_start = structure->get_n_entries_up_to_row(i);

    for (size_t j = 0; j < n_in_row; j++)
    {
      ret[i] += entries[row_start + j];
    }
  }

  return ret;
}

std::vector<double> TMatrix::get_col_sums() const
{
  const int *columns = structure->get_vector_columns();

  size_t n_rows = get_n_rows();
  size_t n_cols = get_n_columns();
  std::vector<double> ret(n_cols, 0.0);

  for (size_t i = 0; i < n_rows; i++)
  {
    size_t n_in_row = structure->get_n_entries_in_row(i);
    size_t row_start = structure->get_n_entries_up_to_row(i);

    for (size_t j = 0; j < n_in_row; j++)
    {
      size_t index = row_start + j;
      ret[columns[index]] += entries[index];
    }
  }

  return ret;
}

/* *********************************************************** */
double TMatrix::GetNorm(int p) const
{
  double result = 0.0;
  switch(p)
  {
    case -2:
      result = Dnorm(this->get_n_entries(), this->GetEntries());
      break;
    case -1:
    {
      const int * rows = this->get_row_ptr();
      for(int row=0; row<this->get_n_rows(); row++)
      {
        double row_sum = 0.0;
        //#pragma omp parallel for
        for(int i=rows[row]; i<rows[row+1]; i++)
        {
          row_sum += std::abs(entries[i]);
        }
        if(row_sum>result && row_sum!=1.0)
          result = row_sum;
      }
      break;
    }
    case 0:
      for(int i=0; i<this->get_n_entries(); i++)
      {
        double a = std::abs(entries[i]);
        if(a > result)
          result = a;
      }
      break;
    case 1:
      ErrThrow("maximum absolute column sum norm of a matrix not yet ",
               "implemented!");
      break;
    case 2:
      ErrThrow("spectral norm of a matrix not yet implemented!");
      break;
    default:
      ErrThrow("undefined norm of a matrix!");
      break;
  }
  return result;
}


double* operator*(const TMatrix & A,const double* x)
{
  const double *AEntries = A.GetEntries();
  const int *ARowPtr = A.get_row_ptr();
  const int *AColIndex = A.get_vector_columns();

  int nrows = A.get_n_rows();

  double *y = new double[nrows];
  double value;
  int index;

  for(int i=0;i<nrows;i++)
  {
    value = 0;
//#pragma omp parallel for
    for (int j=ARowPtr[i]; j<ARowPtr[i+1]; j++)
    {

      index = AColIndex[j];
      value += AEntries[j] * x[index];
    }
    y[i] = value;
  }
  return y;
}

TMatrix & TMatrix::operator+=(const TMatrix* A)
{
  if(this->GetStructure() != A->GetStructure()) // compare objects
  {
    this->GetStructure().info();
    A->GetStructure().info();
    Output::print("diff in norms ", this->GetNorm(-2) - A->GetNorm(-2));
    ErrThrow("TMatrix::operator+= : the two matrices do not match.");
  }

  int n_entries = this->get_n_entries();
  const double *AEntries = A->GetEntries();
  for(int i = 0; i < n_entries; ++i)
  {
    this->entries[i] += AEntries[i];
  }
  return *this;
}



TMatrix & TMatrix::operator-=(const TMatrix* A)
{
  if(this->GetStructure() != A->GetStructure()) // compare objects
  {
    ErrThrow("TMatrix::operator-= : the two matrices do not match.");
  }

  int n_entries = this->get_n_entries();
  const double *AEntries = A->GetEntries();
  for(int i = 0; i < n_entries; ++i)
  {
    this->entries[i] -= AEntries[i];
  }
  return *this;
}


void TMatrix::multiply(const double * const x, double * const y, double a) const
{
  if(a == 0.0)
    return;

  const int *rowPtr = get_row_ptr();
  const int *colIndex = get_vector_columns();

  int nrows = get_n_rows();

  double value;
  int i, j, end;
  for(i = 0; i < nrows; i++)
  {
    value = 0;
    end = rowPtr[i+1];
    for (j = rowPtr[i]; j < end; j++)
    {
      value += entries[j] * x[colIndex[j]];
    }
    y[i] += a * value;
  }
}

void TMatrix::transpose_multiply(const double * const x, double * const y, double a)
      const
{
  if(a == 0.0)
    return;

  const int *rowPtr = get_row_ptr();
  const int *colIndex = get_vector_columns();

  int nrows = get_n_rows();

  double value = 0.;
  int i, j, end;
  for(i = 0; i < nrows; i++)
  {
    value = a*x[i];
    end = rowPtr[i + 1];
    for(j = rowPtr[i]; j < end; j++)
    {
      // Entries[j] is the (i,colIndex[j])-th entry of this matrix
      y[colIndex[j]] += entries[j] * value;
    }
  }
}

TMatrix* TMatrix::multiply(const TMatrix * const B, double a) const
{
  const int n_A_rows = this->get_n_rows();   // = n_C_rows
  const int n_A_cols = this->get_n_columns();
  const int n_B_rows = B->get_n_rows();

  if(n_A_cols != n_B_rows)
  {
    ErrThrow("Dimension mismatch during matrix-matrix multiplication: n_rows in this->matrix:", n_A_cols,", n_rows in new matrix:", n_B_rows );
  }
  const int * const a_rows = this->get_row_ptr();
  const int * const a_cols = this->get_vector_columns();

  const TStructure & strucB = B->GetStructure();
  const double * const b_entries = B->GetEntries();

  std::shared_ptr<TStructure> struc_c = get_product_structure(this->GetStructure(), strucB);
  const int * c_rows = struc_c->get_row_ptr();
  const int * c_cols = struc_c->get_vector_columns();
  TMatrix * c = new TMatrix(struc_c);
  double * c_entries = c->GetEntries();

  // fill the entries
  // loop over all rows in C
//#pragma omp parallel for
  for(int row = 0; row < n_A_rows; row++)
  {

    // loop over all entries in this row in C
    for(int col = c_rows[row]; col < c_rows[row + 1]; col++)
    {
      // multiply 'this row of A' x 'this column of B'
      // loop over all entries in this row in A
      for(int i = a_rows[row]; i < a_rows[row+1]; i++)
      {
        int ib = strucB.get_index_of_entry(a_cols[i], c_cols[col]);
        if(ib != -1)
        {
          c_entries[col] += entries[i] * a * b_entries[ib];
        }
      }
    }
  }
  return c;
}

std::shared_ptr<DenseMatrix> TMatrix::multiply(const DenseMatrix* const B,
                                               bool transpose,
                                               double a) const
{
  double value;
  const int n_A_rows = this->get_n_rows();
  const int n_A_cols = this->get_n_columns();
  const int n_B_rows = transpose ? B->getNColumns() : B->getNRows();
  const int n_B_cols = transpose ? B->getNRows() : B->getNColumns();
  const int n_C_rows = n_A_rows;
  const int n_C_cols = n_B_cols;
  const int LeadingDimension_B = B->getLeadingDimension();
  const int LeadingDimension_C = n_C_rows;

  if(n_A_cols != n_B_rows)
  {
    ErrThrow("Dimension mismatch during matrix-matrix multiplication: ",
             n_A_cols, " n_cols in this->matrix and ",
             n_B_rows, " n_rows in the other one.");
  }

  auto C = std::make_shared<DenseMatrix>(n_C_rows, n_C_cols);

  const int * const a_rows = this->get_row_ptr();
  const int * const a_cols = this->get_vector_columns();

  int idx_b, idx_c;
  const double * const b_entries = B->get_entries();
  double * c_entries = C->get_entries();

  // fill the entries
  // loop over all rows in C
//#pragma omp parallel for
  for(int row = 0; row < n_A_rows; row++)
  {
    // loop over all entries in this row in C
    for(int col = 0; col < n_C_cols; col++)
    {
      value = 0.;

      // multiply 'this row of A' x 'this column of B'
      // loop over all entries in this row in A
      for(int i = a_rows[row]; i < a_rows[row+1]; i++)
      {
        idx_b = transpose ? a_cols[i] * LeadingDimension_B + col
                          : col * LeadingDimension_B + a_cols[i];
        value += a * entries[i] * b_entries[idx_b];
      }
      idx_c = col * LeadingDimension_C + row;
      c_entries[idx_c] = value;
    }
  }
  return C;
}

TMatrix* TMatrix::multiply(const TMatrix* const B,
                           const std::vector<double>& d) const
{
  const int n_A_rows = this->get_n_rows();   // = n_C_rows
  const int n_A_cols = this->get_n_columns();
  const int n_B_rows = B->get_n_rows();

  if(n_A_cols != n_B_rows)
  {
    ErrThrow("dimension mismatch during matrix-matrix multiplication", n_A_cols,"    ", n_B_rows );
  }
  const int * const a_rows = this->get_row_ptr();
  const int * const a_cols = this->get_vector_columns();

  const TStructure & strucB = B->GetStructure();
  const double * const b_entries = B->GetEntries();

  std::shared_ptr<TStructure> struc_c = get_product_structure(this->GetStructure(), strucB);
  const int * c_rows = struc_c->get_row_ptr();
  const int * c_cols = struc_c->get_vector_columns();
  TMatrix * c = new TMatrix(struc_c);
  double * c_entries = c->GetEntries();

  // fill the entries
  // loop over all rows in C
//#pragma omp parallel for
  for(int row = 0; row < n_A_rows; row++)
  {
    // loop over all entries in this row in C
    for(int col = c_rows[row]; col < c_rows[row + 1]; col++)
    {
      // multiply 'this row of A' x 'this column of B'
      // loop over all entries in this row in A
      for(int i = a_rows[row]; i < a_rows[row+1]; i++)
      {
        int ib = strucB.get_index_of_entry(a_cols[i], c_cols[col]);
        if(ib != -1)
        {
          c_entries[col] += entries[i]  * d[a_cols[i]] * b_entries[ib];
        }
      }
    }
  }
  return c;
}

TMatrix* TMatrix::multiply_with_transpose_from_right(
#ifdef _MPI
  const std::vector<const TParFECommunicator3D*>& test_comms
  , const std::vector<const TParFECommunicator3D*>& ansatz_comms
  , bool additive_storage
#endif
) const{

  // put up a unity diagonal matrix as a vector
  std::vector<double> unityScaling(structure->get_n_columns(), 1.0);

  // let the other implementations do the work.
  return multiply_with_transpose_from_right(unityScaling
#ifdef _MPI
      ,test_comms, ansatz_comms, additive_storage
#endif
  );
}

TMatrix* TMatrix::multiply_with_transpose_from_right(
  const std::vector<double>& diagonalScaling
#ifdef _MPI
  , const std::vector<const TParFECommunicator3D*>& test_comms
  , const std::vector<const TParFECommunicator3D*>& ansatz_comms
  , bool additive_storage
#endif
) const
{
  // put up the structure by call to the specific structure generating method
  TStructure * productStructure =
    structure->get_structure_of_product_with_transpose_from_right();

  // construct the matrix and return a pointer
  TMatrix * ret =  multiply_with_transpose_from_right(diagonalScaling,
                                                      *productStructure
#ifdef _MPI
      ,test_comms, ansatz_comms, additive_storage
#endif
  );

  delete productStructure; // this has been (deep) copied in the above method
  return ret;
}

TMatrix* TMatrix::multiply_with_transpose_from_right(
  const std::vector<double>& diagonalScaling, const TStructure& knownStructure
#ifdef _MPI
  , const std::vector<const TParFECommunicator3D*>& test_comms
  , const std::vector<const TParFECommunicator3D*>& ansatz_comms
  , bool additive_storage
#endif
)
const
{
  //check if the dimensions match
  if(diagonalScaling.size() != structure->get_n_columns())
  {
    ErrThrow("Dimension mismatch! ", diagonalScaling.size(), "  ",
             structure->get_n_columns());
  }

#ifdef _MPI
  //check the sparsity type
  if(this->sparse_type == SparsityType::B_TIMES_BT)
    ErrThrow("Multiplying a matrix of SparsityType::B_TIMES_BT with "
             "its transposed breaks parallel functionality.");
  //check if number of dofs fit
  int n_col_dof_comm = 0;
  int n_row_dof_comm = 0;
  for(auto tc : test_comms)
    n_col_dof_comm += tc->GetNDof();
  if(n_col_dof_comm != this->get_n_rows())
    ErrThrow("Test space communicators hold wrong number of dof.");
  for(auto ac : ansatz_comms)
    n_row_dof_comm += ac->GetNDof();
  if(n_row_dof_comm != this->get_n_columns())
    ErrThrow("Ansatz space communicators hold wrong number of dof.");
  //fill a vector of ansatz master processes
  std::vector<int> ansatz_masters(n_row_dof_comm,0);
  size_t shift = 0;
  for(auto ac : ansatz_comms)
  {
    const int* copy_start = ac->GetMaster();
    int n_to_copy = ac->GetNDof();
    std::copy(copy_start, copy_start + n_to_copy,
              ansatz_masters.begin() + shift);
    shift += n_to_copy;
  }
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  const size_t nProductRows = structure->get_n_rows();

  // copy construct the product's TStructure
  std::shared_ptr<TStructure> productStructure =
    std::make_shared<TStructure>(knownStructure);
  // create new matrix with this structure
#ifndef _MPI
  TMatrix* product = new TMatrix(productStructure);
#endif
#ifdef _MPI
  TMatrix* product = new TMatrix(productStructure, SparsityType::B_TIMES_BT);
#endif

  int * productRowPtr = productStructure->get_row_ptr();
  double * productEntries = product->GetEntries();

  // fill the entries
  // loop over all rows in the product
  for(size_t row = 0; row < nProductRows; row++)
  {
    // loop over all entries in this row in the product
    for(int iEntries = productRowPtr[row]; iEntries < productRowPtr[row + 1];
        ++iEntries)
    {
      // hold row1 and row2 to fix ideas
      int row1 = row;
      int row2 = productStructure->get_vector_columns()[iEntries];

      //store begin and end indices for the columns and entries segments
      int beginRow1 = structure->get_row_ptr()[row1];
      int endRow1 = structure->get_row_ptr()[row1+1];
      int beginRow2 = structure->get_row_ptr()[row2];
      int endRow2 = structure->get_row_ptr()[row2+1];

      // work on two segments of column array
      const int* row1ColBegin = &structure->get_vector_columns()[beginRow1];
      const double* row1EntriesBegin = &entries[beginRow1];
      const int row1SegmentSize = endRow1 - beginRow1;
      const int* row2ColBegin = &structure->get_vector_columns()[beginRow2];
      const double* row2EntriesBegin = &entries[beginRow2];
      const int row2SegmentSize = endRow2 - beginRow2;

      //initialize the value to be written later on
      double temp = 0;
      //better use index for control of the loop!
      size_t index1 = 0;
      size_t index2 = 0;
      while ((int)index1 < row1SegmentSize && (int)index2 < row2SegmentSize)
      {
        if (row1ColBegin[index1] > row2ColBegin[index2])
        {
          index2++;
        }
        else if (row2ColBegin[index2] > row1ColBegin[index1])
        {
          index1++;
        }
        else
        {
          // we found a pair of indices with equal entry in KCol.
#ifdef _MPI
          int col = row1ColBegin[index1];
          if(additive_storage && ansatz_masters.at(col) == rank)
          {//the coupling column is master? take this entry!
#endif
            temp +=  row1EntriesBegin[index1]
                   * row2EntriesBegin[index2]
                   * diagonalScaling[row1ColBegin[index1]];
#ifdef _MPI
          }
#endif
          index1++;
          index2++;
        }
      }
      //set the current entry to be the freshly calculated vector product
      productEntries[iEntries]= temp;
    } //end loop over all entries in this row
  } //end loop over all rows of the product

  // return the pointer
  return product;
}

std::shared_ptr< TMatrix >
  TMatrix::multiply_with_transpose_from_right(const TMatrix& B
#ifdef _MPI
  , const std::vector<const TParFECommunicator3D*>&
  , const std::vector<const TParFECommunicator3D*>&
#endif
  ) const
{
#ifdef _MPI
  ErrThrow("multiply_with_transpose_from_right(const TMatrix& B) "
           "is not finished yet for MPI.");
#endif
  // dimension check
  if(B.get_n_rows() != this->get_n_columns()
    || B.get_n_columns() !=this->get_n_columns())
  {
    ErrThrow("Dimension mismatch  ", B.get_n_rows(), "  ", this->get_n_columns());
  }

  // construct a product structure
  std::shared_ptr<TStructure> productStructure(
    structure->get_structure_of_product_with_transpose_from_right(B.GetStructure()));

  // lambda function which returns the entry in the BA^T matrix
  auto return_BAT_entry = [this, B](int i, int j)
  {
    const int* row1ColBegin = &B.get_vector_columns()[B.get_row_ptr()[i]];
    const int row1ColSize = B.get_row_ptr()[i+1] - B.get_row_ptr()[i];
    const double *entriesB = B.GetEntries();

    const int* row2ColBegin = &this->get_vector_columns()[this->get_row_ptr()[j]];
    const int row2ColSize = this->get_row_ptr()[j+1]
                               - this->get_row_ptr()[j];

    size_t indexB = 0;
    size_t indexA = 0;

    double entry_in_BAT_product = 0;
    while ((int)indexB < row1ColSize && (int)indexA < row2ColSize)
    {
      if (row1ColBegin[indexB] > row2ColBegin[indexA])
      {
        indexA++;
      }
      else if (row2ColBegin[indexA] > row1ColBegin[indexB])
      {
        indexB++;
      }
      else
      {
        entry_in_BAT_product +=  entriesB[indexB+B.get_row_ptr()[i]]
                               * entries[indexA+this->get_row_ptr()[j]];
        indexB++;
        indexA++;
      }
    }
    return entry_in_BAT_product;
  };

  // number of rows in product structure
  const size_t nProductRows = productStructure->get_n_rows();
  int * productRowPtr = productStructure->get_row_ptr();

  // create product matrix
  std::shared_ptr<TMatrix> productMatrix
          = std::make_shared<TMatrix>(productStructure);
  // entries in the product matrix
  double * productEntries = productMatrix->GetEntries();

  // fill the entries
  // loop over all rows in the product
  for(unsigned int row=0; row<nProductRows; row++)
  {
    unsigned int begin = productRowPtr[row];
    unsigned int end = productRowPtr[row+1];
    // loop over entries in "this row" in the product
    for(unsigned int iEntries = begin; iEntries<end; iEntries++)
    {
      int rowABAT = row;
      int colABAT = productStructure->get_vector_columns()[iEntries];

      // store begin and end indices for columns and entries segments
      int beginRowA = this->GetStructure().get_row_ptr()[rowABAT];
      int endRowA   = this->GetStructure().get_row_ptr()[rowABAT+1];

      double entry_ABAT = 0;
      for(int k=beginRowA; k<endRowA; ++k)
      {
        // compute the entry in the BA^T
        double BAT_entry = return_BAT_entry(this->get_vector_columns()[k], colABAT);
        double a_entry = entries[k];
        entry_ABAT += a_entry * BAT_entry;
      }
      // set the current entry
      productEntries[iEntries] = entry_ABAT;
    } //endfor loop over all entries in this row
  }// endfor loop over all entries in the product
  // return point of this matrix
  return productMatrix;
}

TMatrix* TMatrix::multiply_with_pseudotranspose_from_right(
  const TMatrix& BT
#ifdef _MPI
  , const std::vector<const TParFECommunicator3D*>& test_comms
  , const std::vector<const TParFECommunicator3D*>& ansatz_comms
  , bool additive_storage
#endif
) const
{
  // put up a unity diagonal matrix as a vector
  std::vector<double> unityScaling(structure->get_n_columns(), 1.0);

  // let the other implementations do the work.
  return multiply_with_pseudotranspose_from_right(BT, unityScaling
#ifdef _MPI
    , test_comms, ansatz_comms, additive_storage
#endif
  );
}

TMatrix* TMatrix::multiply_with_pseudotranspose_from_right(
  const TMatrix& BT, const std::vector<double>& diagonalScaling
#ifdef _MPI
  , const std::vector<const TParFECommunicator3D*>& test_comms
  , const std::vector<const TParFECommunicator3D*>& ansatz_comms
  , bool additive_storage
#endif
) const
{
  // put up the structure by call to the specific structure generating method
  TStructure * productStructure =
    structure->get_structure_of_product_with_transpose_from_right();

  // construct the matrix and return a pointer
  TMatrix * ret =  multiply_with_pseudotranspose_from_right(BT,
    diagonalScaling, *productStructure
#ifdef _MPI
    ,test_comms, ansatz_comms, additive_storage
#endif
  );

  delete productStructure; // this has been (deep) copied in the above method
  return ret;
}

TMatrix* TMatrix::multiply_with_pseudotranspose_from_right(
  const TMatrix& BT,
  const std::vector<double>& diagonalScaling, const TStructure& knownStructure
#ifdef _MPI
  , const std::vector<const TParFECommunicator3D*>& test_comms
  , const std::vector<const TParFECommunicator3D*>& ansatz_comms
  , bool additive_storage
#endif
)
const
{
  // check if the dimensions match
  if (diagonalScaling.size() != structure->get_n_columns())
  {
    ErrThrow("Dimension mismatch! ", diagonalScaling.size(), "  ",
             structure->get_n_columns());
  }

  if (structure->get_n_columns() != BT.structure->get_n_rows()
    || structure->get_n_rows() != BT.structure->get_n_columns())
  {
    ErrThrow("Dimension mismatch! B^T and B' are different sizes.");
  }

#ifdef _MPI

  if (!additive_storage)
  {
    ErrThrow("multiply_with_pseudotranspose_from_right is not implemented "
      "for non-additive storage!");
  }

  // check the sparsity type

  if (sparse_type == SparsityType::B_TIMES_BT
    || BT.sparse_type == SparsityType::B_TIMES_BT)
  {
    ErrThrow("Multiplying a matrix of SparsityType::B_TIMES_BT with "
             "its transposed breaks parallel functionality.");
  }

  // check if number of dofs fit
  int n_col_dof_comm = 0;
  int n_row_dof_comm = 0;

  for (auto tc: test_comms)
  {
    n_col_dof_comm += tc->GetNDof();
  }

  if (n_col_dof_comm != this->get_n_rows())
  {
    ErrThrow("Test space communicators hold wrong number of dof.");
  }

  for (auto ac: ansatz_comms)
  {
    n_row_dof_comm += ac->GetNDof();
  }

  if (n_row_dof_comm != this->get_n_columns())
  {
    ErrThrow("Ansatz space communicators hold wrong number of dof.");
  }

  // fill a vector of ansatz master processes
  std::vector<int> ansatz_masters(n_row_dof_comm, 0);
  size_t shift = 0;
  for (auto ac: ansatz_comms)
  {
    const int* copy_start = ac->GetMaster();
    int n_to_copy = ac->GetNDof();

    std::copy(copy_start, copy_start + n_to_copy,
              ansatz_masters.begin() + shift);

    shift += n_to_copy;
  }

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  const size_t n_product_rows = structure->get_n_rows();
  const size_t n_product_cols = BT.structure->get_n_rows();

  // copy construct the product's TStructure
  std::shared_ptr<TStructure> productStructure =
    std::make_shared<TStructure>(knownStructure);

  // create new matrix with this structure
#ifndef _MPI
  TMatrix* product = new TMatrix(productStructure);
#else
  TMatrix* product = new TMatrix(productStructure, SparsityType::B_TIMES_BT);
#endif

  const int* product_rows = productStructure->get_row_ptr();
  const int* product_cols = productStructure->get_vector_columns();

  const int* B_rows = structure->get_row_ptr();
  const int* B_cols = structure->get_vector_columns();

  const int* BT_rows = BT.structure->get_row_ptr();
  const int* BT_cols = BT.structure->get_vector_columns();

  std::vector<std::vector<int>> BT_col_indices(n_product_rows);
  std::vector<std::vector<int>> BT_col_rows(n_product_rows);

  for (size_t row = 0; row < n_product_cols; row++)
  {
#ifdef _MPI
    if (additive_storage)
    {
      if (ansatz_masters[row] != rank)
      {
        // skip B' rows that belong to other processes
        continue;
      }
    }
    else
    {
      // non-additive storage is not implemented
      break;
    }
#endif

    int BT_0 = BT_rows[row];
    int BT_1 = BT_rows[row + 1];

    for (int j = BT_0; j < BT_1; j++)
    {
      BT_col_rows[BT_cols[j]].push_back((int)row);
      BT_col_indices[BT_cols[j]].push_back(j);
    }
  }

  double* productEntries = product->GetEntries();

  for (size_t prod_row = 0; prod_row < n_product_rows; prod_row++)
  {
    int prod_0 = product_rows[prod_row];
    int prod_1 = product_rows[prod_row + 1];

    int B_0 = B_rows[prod_row];
    int B_1 = B_rows[prod_row + 1];

    for (int prod_i = prod_0; prod_i < prod_1; prod_i++)
    {
      int col_BT = product_cols[prod_i];
      auto& BT_rows = BT_col_rows[col_BT];
      auto& BT_indices = BT_col_indices[col_BT];
      int BT_n = (int)BT_indices.size();

      double sum = 0.0;

      int BT_k = 0;
      for (int B_i = B_0; B_i < B_1; B_i++)
      {
        int col_B = B_cols[B_i];

        while (BT_k < BT_n)
        {
          int row_BT = BT_rows[BT_k];

          if (row_BT == col_B)
          {
            sum += entries[B_i] * BT.entries[BT_indices[BT_k]]
              * diagonalScaling[col_B];
            break;
          }
          else if (row_BT > col_B)
          {
            if (BT_k > 0)
            {
              --BT_k;
            }

            break;
          }

          ++BT_k;
        }

        if (BT_k >= BT_n)
        {
          break;
        }
      }

      productEntries[prod_i] = sum;
    }
  }

  return product;
}

std::shared_ptr<TMatrix>
  TMatrix::multiply_with_pseudotranspose_from_right(const TMatrix&,
    const TMatrix&
#ifdef _MPI
  , const std::vector<const TParFECommunicator3D*>&
  , const std::vector<const TParFECommunicator3D*>&
#endif
  ) const
{
  ErrThrow("multiply_with_pseudotranspose_from_right with a central matrix is "
    "not yet implemented.");
}

TMatrix* TMatrix::get_transposed() const
{
  // get transposed structure
  std::shared_ptr<TStructure> structureT = structure->get_transposed();
  const int * rowsT= structureT->get_row_ptr();
  const int * colsT= structureT->get_vector_columns();
  const int * rows = this->get_row_ptr();
  const int * cols = this->get_vector_columns();

  // transpose the entries:
  TMatrix * mT = new TMatrix(structureT);
  double *entriesT = mT->GetEntries();
  // loop over all rows of the original matrix
  for(int i=0; i<this->get_n_rows(); i++)
  {
    // loop over all entries in this row of the original matrix
    for(int j=rows[i]; j<rows[i+1]; j++)
    {
      // cols[j] is the column of this entry in the original matrix,
      // it corresponds to a row of the transposed matrix
      // look for the column index in that row in the transposed matrix
      // which equals this (non-transposed) row index
      for(int k=rowsT[cols[j]]; k<rowsT[cols[j]+1] ;k++)
      {
        if(i==colsT[k])
        {
          entriesT[k] = entries[j];
          continue; // entry found
        }
      }
    }
  }

  return mT;
}

void TMatrix::remove_zeros(double tol)
{
  if(!this->structure.unique()) //could be adapted in time: work on a copy
    ErrThrow("Cannot remove zeroes in a structure shared by multiple matrices!");

  if(tol < 0)
    tol = this->GetNorm(0) * 1e-15; // largest (in magnitude) entry

  int *row_ptr = structure->get_row_ptr();
  int n_rows = structure->get_n_rows();
  int* kcol = structure->get_vector_columns();
  //"entries" will be treated as the std::vector it is

  int row_begin_old = row_ptr[0];

  for(int i = 0 ; i<n_rows ; ++i)
  {
    int row_end_old = row_ptr[i+1];
    std::vector<int> kcol_part_new;
    std::vector<double> entries_part_new;
    kcol_part_new.reserve( row_end_old -row_begin_old );
    entries_part_new.reserve( row_end_old -row_begin_old );

    for(int j = row_begin_old; j < row_end_old; ++j)
    {
      if(std::abs(entries.at(j) - 0) > tol) //entry does not count as zero
      {
        kcol_part_new.push_back(kcol[j]);
        entries_part_new.push_back(entries.at(j));
      }
    } //new kcol and entries part for row i are filled

    //copy new parts into old arrays
    int row_begin_new = row_ptr[i]; //was updated in loop step i-1
    size_t n_entries_row_new = kcol_part_new.size();
    if(n_entries_row_new != 0) //this is not a zero row
    {
      memcpy(&kcol[row_begin_new], &kcol_part_new.at(0), n_entries_row_new*sizeof(int));
      memcpy(&entries.at(row_begin_new), &entries_part_new.at(0), n_entries_row_new*sizeof(double)); //stl here!!
    }
    //update row_ptr
    row_begin_old = row_ptr[i+1];
    row_ptr[i+1] = row_begin_new + n_entries_row_new;
  }

  //Re-fit the TStructure (nEntries and columns array)
  structure->reset_n_entries();
  //Re-fit the entries array (throw out all trailing (nearly) zeroes)
  entries.resize(structure->get_n_entries());


}


void TMatrix::sor_sweep(const double* b, double* x, double omega, size_t flag
#ifdef _MPI
    , const std::string& par_strat, const TParFECommunicator3D& comm
#endif
    )
const
{
  if (flag > 2)
  {
    ErrThrow("TMatrix::sor_sweep with flag not 0,1, or 2.");
  }

  if (!this->is_square())
  {
    ErrThrow("TMatrix::sor_sweep for non-square matrix is not tested");
  }

  // make sure all diagonal entries are non-zero
  const std::vector<double> diagonal = this->get_diagonal();

  if (std::find_if(diagonal.begin(), diagonal.end(),
    [] (const double& d)
    {
      return d == 0.0;
    })
    != diagonal.end())
  {
    ErrThrow("There is a zero on the diagonal. You can not use `sor` in this "
             "case");
  }

#ifdef _MPI //decide which parallelization strategy to use
  bool average_at_interface = false;
  bool skip_halos = false;

  if (par_strat == std::string("own_cells"))
  {
    skip_halos = true;
    average_at_interface = true;
  }
  else if (par_strat != std::string("all_cells"))
  {
    ErrThrow("Parallelization strategy ", par_strat, " is not implemented.");
  }
#endif

  size_t n_rows = this->get_n_rows();
  int * row_ptr = this->structure->get_row_ptr();
  int * col_ptr = this->structure->get_vector_columns();

  // a lambda function doing a forward or backward solve. Note that in case
  // both steps are done, one could find a better implementation which is not
  // twice as long (in flops) as two individual solves. It would require more
  // memory though.
  auto do_sweep = [&](bool backward)
  {
#ifdef _MPI
    const char* markers = comm.get_dof_markers();
#endif

    for (size_t r = 0; r < n_rows; ++r)
    {
      size_t row = r;

      if (backward)
      {
        row = n_rows - 1 - r;
      }

#ifdef _MPI
      if (skip_halos)
      {
        bool is_halo = markers[row] == 'h' || markers[row] == 'H';

        if (is_halo) //the halos won't be updated
        {
          //Output::print("Halo skipped!");
          continue;
        }
      }
#endif

      // multiply sol with the current row of this matrix
      double sol_x_row = b[row];
      size_t row_begin = row_ptr[row];
      size_t row_end   = row_ptr[row + 1];

      // loop over all entries in this row (index is the index in the entries
      // vector)

      for (size_t index = row_begin; index < row_end; ++index)
      {
        sol_x_row -= this->entries[index] * x[col_ptr[index]];
      }

      x[row] = x[row] + omega * sol_x_row / diagonal[row];
    }
  };

  // do the actual solving step(s)

  if (flag == 0 || flag == 2) // forward_sweep
  {
    do_sweep(false);
  }

  if (flag == 1 || flag == 2) // backward_sweep
  {
    do_sweep(true);
  }

#ifdef _MPI
  if (average_at_interface)
  {
    // compute an all-processors average at interface masters and slaves
    comm.average_at_interface(x);
  }
#endif
}

/* ******************************************************************************** */

void TMatrix::add_scaled(const TMatrix& m, double factor)
{
	if (this->GetStructure() != m.GetStructure()) // compare objects
	{
		int pos = 0;
		for (int i = 0; i < m.get_n_rows(); ++i)
		{
			int begin = m.get_row_ptr()[i];
			int end   = m.get_row_ptr()[i + 1];

			for (int j = begin; j < end; ++j)
			{
				double val = factor * m.entries[pos];
				this->add(i, m.get_vector_columns()[pos], val);
				++pos;
			}
		}
	}
	else
	{
		Daxpy(this->get_n_entries(), factor, m.GetEntries(), this->GetEntries());
	}
}

#ifdef _MPI

void TMatrix::add_scaled(const TMatrix& m, double factor,
  const TParFECommunicator3D &comm, bool masters_only)
{
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  const int* masters = comm.GetMaster();

  const int* my_rows = get_row_ptr();
  const int* other_rows = m.get_row_ptr();

  const int* my_cols = get_vector_columns();
  const int* other_cols = m.get_vector_columns();

  for (int row = 0; row < m.get_n_rows(); row++)
  {
    int my_0 = my_rows[row];
    int my_1 = my_rows[row + 1];

    int other_0 = other_rows[row];
    int other_1 = other_rows[row + 1];

    int other_i = other_0;

    for (int my_i = my_0; my_i < my_1 && other_i < other_1; my_i++)
    {
      int my_col = my_cols[my_i];

      if (masters_only && masters[my_col] != mpi_rank)
      {
        continue;
      }

      do
      {
        int other_col = other_cols[other_i];

        if (other_col == my_col)
        {
          entries[my_i] += factor * m.entries[other_i];
          break;
        }
        else if (other_col > my_col)
        {
          if (other_i > other_0)
          {
            --other_i;
          }
          break;
        }

        ++other_i;
      }
      while (other_i < other_1);
    }
  }
}

#endif

void TMatrix::scale(double factor)
{
  Dscal(this->get_n_entries(), factor, this->GetEntries());
}

void TMatrix::scale(const double * const factor, bool from_left)
{
  const int *rowPtr = get_row_ptr();
  const int *colIndex = get_vector_columns();

  if (from_left)
  {
    for (int i = 0, nrows = get_n_rows(); i < nrows; i++)
    {
      int end = rowPtr[i + 1];

      // scale entire row with the same factor
      for (int j = rowPtr[i]; j < end; j++)
      {
        entries[j] *= factor[i];
      }
    }
  }
  else
  {
    for (int i = 0, nrows = get_n_rows(); i < nrows; i++)
    {
      int end = rowPtr[i + 1];
      for(int j = rowPtr[i]; j < end; j++)
      {
        // scale columnwise
        entries[j] *= factor[colIndex[j]];
      }
    }
  }
}

TMatrix & TMatrix::operator*=(const double a)
{
  this->scale(a);
  return *this;
}


void TMatrix::reorderMatrix()
{
  // make a deep copy of the structure in case of other matrices sharing it
  this->copyOwnStructure();
  int l, begin, end;
  double value;

  const int* Row = structure->get_row_ptr();
  int* KCol = structure->get_vector_columns();

  //Output::print("\tPARDISO: reordering of the columns will be performed");
  //Output::print("\tPARDISO:   no back ordering implemented !!!");

  for(unsigned int i=0;i<structure->get_n_rows();i++)
  {
    begin=Row[i];
    end=Row[i+1];
    for(int j=begin;j<end;j++)
    {
      for(int k=j+1;k<end;k++)
      {
        if(KCol[j] > KCol[k])
        {
          l = KCol[j];       value = entries[j];
          KCol[j] = KCol[k]; entries[j] = entries[k];
          KCol[k] = l;       entries[k] = value;
        }
      }
    }
  }
}

/////////////// Routines for periodic boundary conditions /////////////////
/** ************************************************************************* */
void TMatrix::changeRows(const std::map<int,std::map<int,double>>& entries)
{
	if(entries.size() == 0)
		return; // nothing needs to be done

	const int *oldRows = structure->get_row_ptr();
	const int *oldCols = structure->get_vector_columns();

	// find out how many entries there are after all changes are applied, i.e.
	// how many entries are deleted/created
	int offset = 0;
	for(auto it = entries.begin(); it != entries.end(); ++it)
	{
		int row = it->first;
		offset -= oldRows[row+1]-oldRows[row];// number of entries in old structure
		offset += (it->second).size();        // number of entries in new structure
	}

	int n_rows = structure->get_n_rows();// new number of rows = old number of rows
	// new number of columns = old number of columns
	int n_cols = structure->get_n_columns();
	int n_entries = structure->get_n_entries() + offset; // new number of entries

	int *columns = new int[n_entries];  // new pointer to columns
	int *rows = new int[n_rows+1];      // new row pointer
	rows[0] = 0;

	// create new array to store the entries
	std::vector<double> new_entries(n_entries);

	// fill the arrays 'rows', 'columns' and 'new_entries'
	for(int row=0; row<n_rows; row++)
	{
		auto it = entries.find(row);
		if(it == entries.end())
		{
			// this row stays unchanged
			// number of (old) entries in this row
			unsigned int n_old_entries = oldRows[row+1] - oldRows[row];
			// copy pointer to columns in this row
			memcpy(columns+rows[row], oldCols+oldRows[row], n_old_entries*sizeof(int));
			// update row pointer
			rows[row+1] = rows[row] + n_old_entries;
			// copy entries
			memcpy(&new_entries[0]+rows[row], this->GetEntries()+oldRows[row],
					n_old_entries*sizeof(double));
		}
		else
		{
			// this row will be replaced
			std::map<int,double> newRow = it->second;
			// loop over all new entries in this row
			int columnIndex=0;
			for(std::map<int,double>::iterator it2 = newRow.begin();
					it2 != newRow.end(); ++it2)
			{
				int colInd = it2->first; // column index of new entry
				double entry = it2->second; // value of new entry
				columns[columnIndex+rows[row]] = colInd;
				new_entries[columnIndex+rows[row]] = entry;
				columnIndex++;
			}
			rows[row+1] = rows[row] + newRow.size();
			//if(newRow.size() != columnIndex)
			//  Output::print("ERROR: wrong number of columns in this row ",
      //                newRow.size(), "\t", columnIndex, "\t", row);
		}
	}

	// change Structure of this matrix
	this->structure = std::make_shared<TStructure>(n_rows, n_cols, n_entries,
                                                 columns, rows);
	delete [] columns;
	delete [] rows;

	this->entries = new_entries;
}

///////////////// ///////////////// ///////////////// /////////////////

/** ************************************************************************* */
void TMatrix::copyOwnStructure()
{
	if(this->structure.unique())
	{
		// no one else knows this->structure
		return;
  }
  // create a deep copy of the own structure
  std::shared_ptr<TStructure> new_structure(new TStructure(*this->structure));
  this->structure = new_structure;
}


/** ************************************************************************* */
void TMatrix::info(size_t verbose) const
{
  Output::print(" TMatrix M with ", this->get_n_rows(), " rows and ",
                this->get_n_columns(), " columns");
  if(verbose > 2)
  {
    this->Print("M");
  }
}
