#include <FEMatrix.h>
#include "HNDesc.h"
#include "HangingNode.h"
#include <MooNMD_Io.h>
#include <algorithm>


FEMatrix::FEMatrix(std::shared_ptr<const TFESpace1D> space)
: FEMatrix(space, std::make_shared<TStructure>(space))
{
  
}

FEMatrix::FEMatrix(std::shared_ptr<const TFESpace2D> space)
: FEMatrix(space, std::make_shared<TStructure>(space))
{
  
}

#ifdef __3D__
FEMatrix::FEMatrix(std::shared_ptr<const TFESpace3D> space)
: FEMatrix(space, std::make_shared<TStructure>(space))
{
  
}
#endif // 3D

FEMatrix::FEMatrix(std::shared_ptr<const TFESpace2D> space, const TMatrix& m)
  : TMatrix(m), AnsatzSpace1D(nullptr), AnsatzSpace2D(space),
    AnsatzSpace3D(nullptr), TestSpace1D(nullptr), TestSpace2D(space),
    TestSpace3D(nullptr)
{
  int n_dof = space->get_n_dof();
  if(n_dof != this->get_n_rows() || n_dof != this->get_n_columns())
  {
    ErrThrow("unable to construct a FEMatrix using a TMatrix and an ",
             "incompatible space. ", n_dof, " ", this->get_n_rows(), "  ", 
             this->get_n_columns());
  }
}
#ifdef __3D__
FEMatrix::FEMatrix(std::shared_ptr<const TFESpace3D> space, const TMatrix& m)
  : TMatrix(m), AnsatzSpace1D(nullptr), AnsatzSpace2D(nullptr),
    AnsatzSpace3D(space), TestSpace1D(nullptr), TestSpace2D(nullptr),
    TestSpace3D(space)
{
  int n_dof = space->get_n_dof();
  if(n_dof != this->get_n_rows() || n_dof != this->get_n_columns())
  {
    ErrThrow("unable to construct a FEMatrix using a TMatrix and an ",
             "incompatible space. ", n_dof, " ", this->get_n_rows(), "  ", 
             this->get_n_columns());
  }
}
#endif // 3D

FEMatrix::FEMatrix(std::shared_ptr<const TFESpace2D> testspace,
                   std::shared_ptr<const TFESpace2D> ansatzspace, bool is_empty)
 : TMatrix(std::make_shared<TStructure>(testspace, ansatzspace, is_empty)),
   AnsatzSpace1D(nullptr), AnsatzSpace2D(ansatzspace), AnsatzSpace3D(nullptr),
   TestSpace1D(nullptr), TestSpace2D(testspace), TestSpace3D(nullptr)
{
  
}

#ifdef __3D__
FEMatrix::FEMatrix(std::shared_ptr<const TFESpace3D> testspace,
                   std::shared_ptr<const TFESpace3D> ansatzspace, bool is_empty)
: TMatrix(std::make_shared<TStructure>(testspace, ansatzspace, is_empty)),
  AnsatzSpace1D(nullptr), AnsatzSpace2D(nullptr), AnsatzSpace3D(ansatzspace),
  TestSpace1D(nullptr), TestSpace2D(nullptr), TestSpace3D(testspace)
{
 
}
#endif // 3D

FEMatrix::FEMatrix(std::shared_ptr<const TFESpace1D> space,
                   std::shared_ptr<TStructure> structure)
: TMatrix(structure),
  AnsatzSpace1D(space), AnsatzSpace2D(nullptr), AnsatzSpace3D(nullptr),
  TestSpace1D(space), TestSpace2D(nullptr), TestSpace3D(nullptr)
{
  if(!structure->is_square())
  {
    ErrThrow("The structure must be square for this FEMatrix constructor");
  }
  if(space->get_n_dof() != static_cast<int>(structure->get_n_rows()))
  {
    ErrThrow("The given matrix and structure to not properly match");
  }
}

FEMatrix::FEMatrix(std::shared_ptr<const TFESpace2D> space,
                   std::shared_ptr<TStructure> structure)
: TMatrix(structure),
  AnsatzSpace1D(nullptr), AnsatzSpace2D(space), AnsatzSpace3D(nullptr),
  TestSpace1D(nullptr), TestSpace2D(space), TestSpace3D(nullptr)
{
  if(!structure->is_square())
  {
    ErrThrow("The structure must be square for this FEMatrix constructor");
  }
  if(space->get_n_dof() != static_cast<int>(structure->get_n_rows()))
  {
    ErrThrow("The given matrix and structure to not properly match");
  }
}

#ifdef __3D__
FEMatrix::FEMatrix(std::shared_ptr<const TFESpace3D> space,
                   std::shared_ptr<TStructure> structure)
: TMatrix(structure),
  AnsatzSpace1D(nullptr), AnsatzSpace2D(nullptr), AnsatzSpace3D(space),
  TestSpace1D(nullptr), TestSpace2D(nullptr), TestSpace3D(space)
{
  if(!structure->is_square())
  {
    ErrThrow("The structure must be square for this FEMatrix constructor");
  }
  if(space->get_n_dof() != static_cast<int>(structure->get_n_rows()))
  {
    ErrThrow("The given matrix and structure to not properly match");
  }
}
#endif // 3D

void FEMatrix::resetActive()
{
  // numer of entries in active rows
  int nActive = this->get_n_active_rows();
  size_t n_active_entries = this->structure->get_n_entries_up_to_row(nActive);
  std::fill(this->entries.begin(), this->entries.begin()+n_active_entries, 0.0);
}

void FEMatrix::resetNonActive()
{
  // number of entries in active rows (including hanging)
  int nActive = this->GetTestSpace()->get_n_active_non_hanging();
  size_t n_active_entries = this->structure->get_n_entries_up_to_row(nActive);
  std::fill(this->entries.begin() + n_active_entries, this->entries.end(), 0.0);
}

void FEMatrix::scaleActive(double factor)
{
  if(factor == 1.0)
    return; // no scaling
  if(factor == 0.0)
    this->resetActive();
  
  // number of entries in active rows
  int nActive = this->get_n_active_rows();
  size_t n_active_entries = this->structure->get_n_entries_up_to_row(nActive);
  std::for_each(this->entries.begin(), this->entries.begin() + n_active_entries,
                 [factor](double & a){ a = a*factor; } );
}

void FEMatrix::scale_non_active_diagonals(double factor)
{
  const int n_rows = this->get_n_rows();
  const int active_bound = this->GetTestSpace()->get_n_active();
  const int * rowPtr = this->get_row_ptr();
  const int * colIndex = this->get_vector_columns();
  for(int i = active_bound; i < n_rows; ++i)
  {
    for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j)
    {
      if(colIndex[j] == i)
      {
        this->entries[j] *= factor;
      }
    }
  }
}

void FEMatrix::set_dirichlet_diagonals()
{
  const int n_rows = this->get_n_rows();
  const int active_bound = this->GetTestSpace()->get_n_active();
  const int * rowPtr = this->get_row_ptr();
  const int * colIndex = this->get_vector_columns();
  for(int i = active_bound; i < n_rows; ++i)
  {
    for(int j = rowPtr[i]; j < rowPtr[i+1]; ++j)
    {
      if(colIndex[j] == i)
      {
        this->entries[j] = 1.0;
      }
    }
  }
}

void FEMatrix::addActive(const FEMatrix& m, double factor)
{
  if(this->GetStructure() != m.GetStructure()) // compare objects
  {
    ErrThrow("FEMatrix::add : the two matrices do not match.");
  }
  
  // numer of entries in active rows
  int nActive = this->get_n_active_rows();
  size_t n_active_entries = this->structure->get_n_entries_up_to_row(nActive);
  std::transform(this->entries.begin(), this->entries.begin() + n_active_entries,
                 m.entries.begin(), this->entries.begin(), 
                 [factor](const double & a, const double & b)
                 { return a + factor * b; } );
}

void FEMatrix::multiplyActive(const double* x, double* y, double factor) const
{
  int nActive= this->get_n_active_rows();
  const int * rowPtr = this->get_row_ptr();
  const int * colIndex = this->get_vector_columns();
  
  
  for(int i=0; i<nActive; ++i)
  {
    double val=0.;
    for(int j=rowPtr[i]; j<rowPtr[i+1]; ++j)
      val += this->entries[j]*x[colIndex[j]];
    y[i] += factor*val;
  }  
}

void FEMatrix::multiplyTransposedActive(const double *x, double *y, double factor) const
{
  //be careful! we have to rely on y's actives being as many as this ansatz spaces
  //FIXME this can be sped up of course, but for the moment do it like this
  //assume that y is as long as this has columns
  int n_actives = this->GetAnsatzSpace()->get_n_active();
  std::vector<double> y_non_actives(this->GetAnsatzSpace()->get_n_dirichlet());
  for (int i= n_actives; i<this->get_n_columns() ; ++i )
  {//store non-actives
    y_non_actives[i-n_actives]=y[i];
  }
  this->TMatrix::transpose_multiply(x,y,factor);
  for (int i= n_actives; i<this->get_n_columns() ; ++i )
  {//put back non-actives
    y[i] = y_non_actives[i-n_actives];
  }

}

void FEMatrix::ModifyMatrixAccordingToCoupling(bool assemble_dirichlet_rows)
{
  double *Entries = this->GetEntries();
  const int *RowPtr = this->get_row_ptr();
  const int *ColInd = this->get_vector_columns();
  
  auto fespace = this->GetTestSpace();
  auto hanging_nodes = fespace->get_sorted_hanging_nodes();
  int hanging_bound = 0;
  //For AFC schemes we need to assemble the Dirichlet rows and hence the hanging bound increases
  //to all the DOFs
  if(assemble_dirichlet_rows)
    hanging_bound = fespace->get_n_dof();
  else
    hanging_bound = fespace->get_n_active();
  for(auto hn_dof_pair : hanging_nodes)
  {
    int k = hn_dof_pair.first->GetN_Nodes();
    auto Coupling = hn_dof_pair.first->GetCoeff();
    auto DOF = hn_dof_pair.first->GetDOF();
    
    int end = RowPtr[hn_dof_pair.second+1];
    for(int n = RowPtr[hn_dof_pair.second]; n < end; n++)
    {
      double v = Entries[n];
      int m = ColInd[n];
      for(int l=0;l<k;l++)
      {
        int l1 = DOF[l];
        if(l1<hanging_bound)
        {
          bool found = false;
          int last=RowPtr[l1+1];
          for(int l2=RowPtr[l1];l2<last;l2++)
          {
            if(ColInd[l2] == m)
            {
              Entries[l2] += Coupling[l] * v;
              found = true;
              break;
            }
          }                                     // endfor l2  
          if(!found)
          {
            ErrThrow("ModifyMatrixAccordingToCoupling: hanging row ",
                     hn_dof_pair.second, " column ", m,
                     " could not be added to row ", l1,
                     " because the structure does not include an entry here");
          }
        }                                       // endif  
      }                                         // endfor l
    }                                           // endfor n  
  }                                             // endfor i
}

void FEMatrix::ModifyMatrixAccordingToCouplingAFC()
{
  double *Entries = this->GetEntries();
  const int *RowPtr = this->get_row_ptr();
  const int *ColInd = this->get_vector_columns();
  int ndofs = this->get_n_rows();
  
  auto fespace = this->GetTestSpace();
  auto hanging_nodes = fespace->get_sorted_hanging_nodes();
  //hanging_dof : DOF of the hanging node
  //coupling_dof_1 : DOF of the first coupling DOF
  //coupling_dof_2 : DOF of the second coupling DOF
  int hanging_dof = 0, coupling_dof_1 = 0, coupling_dof_2 = 0;
  int hanging_entry = 0, coupling_entry_1 = 0, coupling_entry_2 = 0;
  int check_break = 0;

  for(auto hn_dof_pair : hanging_nodes)
  {
    int k = hn_dof_pair.first->GetN_Nodes();
    auto Coupling = hn_dof_pair.first->GetCoeff();
    auto DOF = hn_dof_pair.first->GetDOF();
    //Store the DOF for hanging node
    hanging_dof = hn_dof_pair.second;
    /** 
    *  NOTE: This is very specific to P1 elements. As AFC is only applicable to
    *  P1, hence this works. One needs modification if this idea needs to be
    *  implemented for higher order elements
    **/
    coupling_dof_1 = DOF[0];
    coupling_dof_2 = DOF[1];
    //Loop over all rows
    for(int i = 0; i<ndofs; i++)
    {
      hanging_entry = 0;
      coupling_entry_1 = 0;
      coupling_entry_2 = 0;
      check_break = 0;
      for(int j = RowPtr[i]; j<RowPtr[i+1]; j++)
      {
        int m = ColInd[j];
        //Output::print("Row : ", i, " ColInd : ", m);
        //Check if hanging dof found in the ith row
        if(hanging_dof == m)
        {
          hanging_entry = j;
          check_break++;
        }
        //Check if coupling dof 1 found in the ith row
        else if(coupling_dof_1 == m)
        {
          coupling_entry_1 = j;
          check_break++;
        }
        //Check if coupling dof 2 found in the ith row
        else if(coupling_dof_2 == m)
        {
          coupling_entry_2 = j;
          check_break++;
        }
        //If all them are found then exit the loop
        if(check_break == 3)
          break;
      }
      if(check_break == 3)
      {
        Entries[coupling_entry_1] += Coupling[0]*Entries[hanging_entry];
        Entries[coupling_entry_2] += Coupling[1]*Entries[hanging_entry];
        //Set the hanging column to unit vector
        if(i == hanging_dof)
          Entries[hanging_entry] = 1.0;
        else
          Entries[hanging_entry] = 0.0;
      }
    }
  }
}

void FEMatrix::correct_hanging_rows()
{
  if(!this->is_square())
    return;
  int colindex;
  auto fespace = this->GetTestSpace();
  int N_ = fespace->get_n_hanging();
  int ActiveBound = fespace->get_n_active_non_hanging();

  const int *RowPtr = this->get_row_ptr();
  const int *ColInd = this->get_vector_columns();
  //Get the Entry position of Diagonal in Hanging Row

  for(int i=0;i<N_;i++)
  {
    const THangingNode *hn = fespace->get_hanging_node(i);
    int k = hn->GetN_Nodes();
    auto Coupling = hn->GetCoeff();
    auto DOF = hn->GetDOF();

    for(int m = RowPtr[ActiveBound+i] ; m < RowPtr[ActiveBound+i+1]; m++)
    {
      colindex = ColInd[m];
      double new_entry = 0.;
      if(colindex == ActiveBound + i)
      {
        // Diagonal Entry
        new_entry = 1.0;
      }
      else
      {
        for(int l = 0; l<k ;l++)
        {
          if(colindex == DOF[l])
          {
            new_entry = -Coupling[l];
            break;
          }
        }
      }
      entries[m] = new_entry;
    }
  }                                             // endfor i
}


int FEMatrix::get_n_active_rows() const
{
  return this->GetTestSpace()->get_n_active();
}

std::shared_ptr<const TFESpace1D> FEMatrix::GetTestSpace1D() const
{
  return TestSpace1D;
}

std::shared_ptr<const TFESpace1D> FEMatrix::GetAnsatzSpace1D() const
{
  return AnsatzSpace1D;
}

std::shared_ptr<const TFESpace2D> FEMatrix::GetTestSpace2D() const
{
  return TestSpace2D;
}

std::shared_ptr<const TFESpace2D> FEMatrix::GetAnsatzSpace2D() const
{
  return AnsatzSpace2D;
}

#ifdef __3D__
std::shared_ptr<const TFESpace3D> FEMatrix::GetTestSpace3D() const
{
  return TestSpace3D;
}

std::shared_ptr<const TFESpace3D> FEMatrix::GetAnsatzSpace3D() const
{
  return AnsatzSpace3D;
}
#endif // 3D

std::shared_ptr<const TFESpace> FEMatrix::GetTestSpace() const
{
  if(TestSpace1D)
  {
    return TestSpace1D;
  }
  else if(TestSpace2D)
  {
    return TestSpace2D;
  }
  else
  {
    return TestSpace3D;
  }
}

std::shared_ptr<const TFESpace> FEMatrix::GetAnsatzSpace() const
{
  if(AnsatzSpace1D)
  {
    return AnsatzSpace1D;
  }
  else if(AnsatzSpace2D)
  {
    return AnsatzSpace2D;
  }
  else
  {
    return AnsatzSpace3D;
  }
}

std::shared_ptr<const TFESpace1D> FEMatrix::GetFESpace1D() const
{
  if(this->structure->is_square())
    return TestSpace1D;
  else
    ErrThrow("accessing FESpace for non-square matrix, but which one?");
}

std::shared_ptr<const TFESpace2D> FEMatrix::GetFESpace2D() const
{
  if(this->structure->is_square())
    return TestSpace2D;
  else
    ErrThrow("accessing FESpace for non-square matrix, but which one?");
}

#ifdef __3D__
std::shared_ptr<const TFESpace3D> FEMatrix::GetFESpace3D() const
{
  if(this->structure->is_square())
    return TestSpace3D;
  else
    ErrThrow("accessing FESpace for non-square matrix, but which one?");
}
#endif // 3D

#ifdef _MPI

#include "mpi.h"
#include <ParFECommunicator3D.h>

/* ************************************************************************ */
// The following are two functions which were used for debug reasons. They do
// not belong to the class FEMatrix.
// Yet they should not be lost entirely, so they were 'parked' here.

/**
 * This method can be used to check a matrix A in distributed storage.
 * A certain process ('sending_ps') will extract one master column after the other
 * from A, and ask all other processes for a consistency update of that vector.
 *
 * After that it will compare the updated column with the original one and inform
 * the programmer via Output::print on any differences.
 * It is then in the responsibility of the programmer to interpret the output.
 *
 * @param[in] A The matrix to check.
 * @param[in] sending_ps The process which checks its master columns.
 * @param[in] cons_level The consistency level to be updated to.
 */
void check_column_consistency(const FEMatrix& A, int sending_ps, int cons_level)
{

  const TParFECommunicator3D& comm = A.GetFESpace3D()->get_communicator();
  const int* masters = comm.GetMaster();
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int nDof = A.get_n_rows();

  const char* markers = comm.get_dof_markers();

  //SENDER process executes the following code.
  if(rank == sending_ps)
  {
    Output::info("CHECK", "Checking master columns of process ", sending_ps, " on consistency level ", cons_level);
    int n_vecs_send = comm.GetN_Master();
    MPI_Bcast(&n_vecs_send, 1, MPI_INT, sending_ps, MPI_COMM_WORLD);
    for(int s = 0; s < nDof; ++s)
    {
      if(masters[s] == rank) //master col found, ping it
      {
        Output::suppressAll();
        comm.dof_ping(sending_ps, s);
        Output::setVerbosity(1);
        int dummy_sb = -1;
        int dof_remote = -1;
        MPI_Allreduce(&dummy_sb, &dof_remote, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if(dof_remote != -1)
        {//ping was received on at least one other ps, now set up the matrix row
          std::vector<double> col_master = A.get_matrix_column(s);
          std::vector<double> col_master_cpy = col_master;
          comm.consistency_update(col_master.data(), cons_level); //CONSIST UPDATE
          //now check the differences between updated and original column
          for(size_t i = 0 ; i < col_master.size(); ++i)
          {
            if(col_master[i] != col_master_cpy[i])
            {
              char type_col = markers[s];
              char type_row = markers[i];
              Output::print("Local entry (", i, "[", type_row ,"] , ", s, "[", type_col,"]) was changed by an update.");
            }

          }

        }
      }
    }
    Output::print("Test on process ", rank, " complete.");
  }
  //RECEIVER processes execute the following code.
  else
  {
    int n_vecs_recv;
    MPI_Bcast(&n_vecs_recv, 1, MPI_INT, sending_ps, MPI_COMM_WORLD);
    for(int r = 0; r< n_vecs_recv; ++r)
    {
      int dof_local = comm.dof_ping(sending_ps, -1);
      int dof_remote = -1;
      MPI_Allreduce(&dof_local, &dof_remote, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if(dof_remote != -1)
      {
        std::vector<double> col_slave = A.get_matrix_column(dof_local);
        comm.consistency_update(col_slave.data(), cons_level); //CONSIST UPDATE
      }
    }
  }
}

// Restores consistency of master columns.
void update_column_consistency(FEMatrix& B, int sending_ps, int cons_level)
{
  const TParFECommunicator3D& comm = B.GetFESpace3D()->get_communicator();
  const int* masters = comm.GetMaster();
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int nDof = B.get_n_rows();

  const char* markers = comm.get_dof_markers();

  //SENDER process executes the following code.
  if(rank == sending_ps)
  {
    int n_vecs_send = comm.get_n_interface_master();
    MPI_Bcast(&n_vecs_send, 1, MPI_INT, sending_ps, MPI_COMM_WORLD);
    Output::info("MPI UPDATE", n_vecs_send, " interface master columns receive update.");
    for(int s = 0; s < nDof; ++s)
    {
      if(masters[s] == rank && markers[s] == 'm') //interface master col found, ping it
      {
        Output::suppressAll();
        comm.dof_ping(sending_ps, s);
        Output::setVerbosity(1);
        int dummy_sb = -1;
        int dof_remote = -1;
        MPI_Allreduce(&dummy_sb, &dof_remote, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if(dof_remote != -1)
        {//ping was received on at least one other ps, now set up the matrix row
          std::vector<double> col_master = B.get_matrix_column(s);
          comm.consistency_update(col_master.data(), cons_level); //CONSIST UPDATE
          B.set_matrix_column(s,col_master);

        }
      }
    }
  }
  //RECEIVER processes execute the following code.
  else
  {
    int n_vecs_recv;
    MPI_Bcast(&n_vecs_recv, 1, MPI_INT, sending_ps, MPI_COMM_WORLD);
    for(int r = 0; r< n_vecs_recv; ++r)
    {
      int dof_local = comm.dof_ping(sending_ps, -1);
      int dof_remote = -1;
      MPI_Allreduce(&dof_local, &dof_remote, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if(dof_remote != -1)
      {
        std::vector<double> col_slave = B.get_matrix_column(dof_local);
        comm.consistency_update(col_slave.data(), cons_level); //CONSIST UPDATE
      }
    }
  }
}
#endif //_MPI
