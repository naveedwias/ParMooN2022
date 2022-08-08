#include "Structure.h"
#include "MooNMD_Io.h"
#include "FESpace.h"
#include "BaseCell.h"
#include "HangingNode.h"
#include <algorithm>


TStructure::TStructure(std::shared_ptr<const TFESpace> space, bool is_empty,
                       bool face_integrals)
  : TStructure(space, space, is_empty, face_integrals)
{ // map on other constructor with twice the same space
}

/// count the number of local dofs in all cells neighboring the given `cell`
/// this is only meaningful if you want a 'far coupling', i.e., two dofs couple
/// if they both couple to some dof in the usual sense.
unsigned int get_n_far_dof(const TBaseCell* cell, const TFESpace& ansatzspace)
{
  unsigned int n_far_dof = 0;
  auto n_faces = cell->GetN_Joints();   // # faces
  for(int joint_k = 0; joint_k < n_faces; joint_k++)
  { // for all faces / neighbours
    auto neigh = cell->GetJoint(joint_k)->GetNeighbour(cell); // neighbour
    if(neigh)
    {
      int neigh_id = neigh->GetClipBoard();
      if (neigh_id < 0)
      {
#ifdef _MPI
        continue;
#else
        ErrThrow("Neighbour has wrong clipboard number. Internal face ",
            "integrals in combination with hanging nodes may not work.");
#endif
      }
      int n_neigh_dof_a = ansatzspace.get_n_local_dof(neigh_id);
      n_far_dof += n_neigh_dof_a;
    } //endif
  } //endfor joint_k
  return n_far_dof;
}

// find all dofs in this cell. If this cell is not part of the fe space's
// collection, this method searches through all children. It returns true if the
// `hanging_dof` was in the cell or one of its children.
bool loop_over_children(const TBaseCell* cell, const TFESpace& ansatzspace,
                        const TFESpace& testspace,
                        std::vector<int>& hanging_dofs,
                        std::function<void(int cell_index, const TFESpace& space)> f)
{
  const TCollection* coll = ansatzspace.GetCollection();
  int cell_index = -1;
  try { cell_index = coll->get_cell_index(cell); }
  catch(const std::out_of_range&){}
  bool found = false;
  if(cell_index != -1)
  {
    f(cell_index, ansatzspace);
    for(int hn_dof : hanging_dofs)
      found |= testspace.is_dof_in_cell(hn_dof, cell_index);
  }
  else
  {
    int n_children = cell->GetN_Children();
    for(int c = 0; c < n_children; ++c)
    {
      found |= loop_over_children(cell->GetChild(c), ansatzspace, testspace,
                                  hanging_dofs, f);
    }
  }
  return found;
}

void add_far_hanging_coupling(int cell_i, int hanging_dof, 
                              const TFESpace& testspace,
                              const TFESpace& ansatzspace,
                              std::vector<int>& upper_bounds)
{
  const TCollection* coll = testspace.GetCollection();
  int ActiveBound = testspace.get_n_active_non_hanging();
  const THangingNode *hn = testspace.get_hanging_node(hanging_dof-ActiveBound);
  std::vector<int> hn_partners = hn->get_partners();
  hn_partners.push_back(hanging_dof);
  int n_coupling_dof = hn->GetN_Nodes();
  const TBaseCell * cell = coll->GetCell(cell_i);
  const TBaseCell * parent = cell->GetParent();
  int n_children = parent->GetN_Children();
  // loop over all dofs which couple with this hanging dof
  for(int i = 0; i < n_coupling_dof; ++i)
  {
    int coupling_dof = hn->GetDOF()[i];
    if(!testspace.is_dof_in_cell(coupling_dof, cell_i))
    {
      // the coupling dof is in a different cell than the hanging_dof, so we 
      // add their coupling explicitly here (does not occur for P1 but for P2 or
      // higher orders)
      upper_bounds[hanging_dof] += 1;
      continue; // this is done from another cell
    }
    // we want to loop over all cells which share the hanging_dof, because this
    // is not so easy to do, we loop over all siblings (children of parent). If
    // a child is not in the collection, we look at its children recursively
    for(int c = 0; c < n_children; ++c)
    {
      const TBaseCell * child = parent->GetChild(c);
      if(cell == child)
        continue;

      int upper_bound_extra = 0;
      bool found = loop_over_children(
        child, ansatzspace, testspace, hn_partners,
        [&upper_bound_extra](int cell_index, const TFESpace& ansatzspace)
        { upper_bound_extra += ansatzspace.get_n_local_dof(cell_index); });
      if(!found)
        continue;
      upper_bounds[coupling_dof] += upper_bound_extra;
    }
  }
}

void add_entry(const std::vector<int>& upper_bounds,
               std::vector<int>& columns_aux, int glob_num_t, int glob_num_a,
               unsigned int& nEntries)
{
  auto index = upper_bounds[glob_num_t];
  auto l = columns_aux[index];

  // check whether this column is already in this row
  while(l!=-1 && l!= glob_num_a && index < upper_bounds[glob_num_t+1])
  {
    index++;
    l = columns_aux[index];
  }
  if(l==-1)
  { // this is a new column for this row
    columns_aux[index] = glob_num_a;
    nEntries++;
  }
}

void add_columns_far_hanging_coupling(int cell_i, int hanging_dof,
                                      const TFESpace& testspace,
                                      const TFESpace& ansatzspace,
                                      const std::vector<int>& upper_bounds,
                                      std::vector<int>& columns_aux,
                                      unsigned int& nEntries)
{
  const TCollection* coll = testspace.GetCollection();
  int ActiveBound = testspace.get_n_active_non_hanging();
  int HangingBound = testspace.get_n_active();
  const THangingNode *hn = testspace.get_hanging_node(hanging_dof-ActiveBound);
  std::vector<int> hn_partners = hn->get_partners();
  hn_partners.push_back(hanging_dof);
  int n_coupling_dof = hn->GetN_Nodes();
  const TBaseCell * cell = coll->GetCell(cell_i);
  const TBaseCell * parent = cell->GetParent();
  int n_children = parent->GetN_Children();
  std::vector<int> considered_cells;
  bool is_square = (&testspace == &ansatzspace);

  auto f = [&considered_cells](int cell_index, const TFESpace& )
           { considered_cells.push_back(cell_index); };

  // loop over all dofs which couple with this hanging dof
  for(int i = 0; i < n_coupling_dof; ++i)
  {
    int coupling_dof = hn->GetDOF()[i];
    if(is_square && !testspace.is_dof_in_cell(coupling_dof, cell_i))
    {
      add_entry(upper_bounds, columns_aux, hanging_dof, coupling_dof, nEntries);
      continue; // this is done from another cell
    }
    // we want to loop over all cells which share the hanging_dof, because this
    // is not so easy to do, we loop over all siblings (children of parent)
    for(int c = 0; c < n_children; ++c)
    {
      const TBaseCell * child = parent->GetChild(c);
      if(cell == child)
        continue;

      considered_cells.clear();
      bool found = loop_over_children(child, ansatzspace, testspace, hn_partners, f);
      if(!found)
        continue;

      for(int child_index : considered_cells)
      {
        unsigned int n_dofs_a = ansatzspace.get_n_local_dof(child_index);
        for(unsigned int dof_k = 0; dof_k < n_dofs_a; dof_k++)
        {
          int ansatz_dof = ansatzspace.GetGlobalDOF(child_index)[dof_k];
          add_entry(upper_bounds, columns_aux, coupling_dof, ansatz_dof,
                    nEntries);
          
          auto used_test_element = testspace.get_fe_type(0);
          auto used_ansatz_element = ansatzspace.get_fe_type(0);
          /** Reason for checking both the spaces is because of the test adaptive_hanging_2d.test
           * This can be done for all elements but one needs to think about it.
           * Currently, this is set for AFC schemes specifically. The database has 
           * no flag for checking AFC scheme and hence using this approach
           **/
          if(used_test_element == C_P1_2D_T_A && used_ansatz_element == C_P1_2D_T_A)
          {
            //Add the transposed entry. For AFC schemes and hanging nodes
            add_entry(upper_bounds, columns_aux, ansatz_dof, coupling_dof,
                   nEntries);
            //Check if the Ansatz dof is a hanging node
            if(ansatz_dof >= ActiveBound && ansatz_dof <HangingBound)
            {
              const THangingNode *hn_new = testspace.get_hanging_node(ansatz_dof-ActiveBound);
              int coupling_dof_new_1 = hn_new->GetDOF()[0];
              int coupling_dof_new_2 = hn_new->GetDOF()[1];
              add_entry(upper_bounds, columns_aux, coupling_dof, coupling_dof_new_1,
                      nEntries);  
              add_entry(upper_bounds, columns_aux, coupling_dof, coupling_dof_new_2,
                      nEntries);
            }
          }
        }
      }
    }
  }
}

TStructure::TStructure(std::shared_ptr<const TFESpace> testspace,
                       std::shared_ptr<const TFESpace> ansatzspace,
                       bool is_empty, bool face_integrals)
  : TStructure()
{
  // exclude cases that are not (correctly) implemented
  if(testspace->GetCollection() != ansatzspace->GetCollection())
  {
    ErrThrow("Structure2D.C : grid for test and ansatz space is not the same!");
  }

  if (face_integrals &&
      (testspace->get_n_hanging() > 0 || ansatzspace->get_n_hanging() >0 ))
  {
    ErrThrow("Structure.C : Internal face integrals in combination with hanging"
        " nodes are not implemented yet.");
  }

#ifdef _MPI
  if ((testspace->get_n_hanging() > 0 || ansatzspace->get_n_hanging() >0 ))
  {
    ErrThrow("Structure.C: Hanging nodes together with MPI does not work.");
  }
#endif

  // internal full matrix structure, false means a reduced structure in the
  // Dirichlet rows, so that there is only one entry on the diagonal in these
  // rows.
  bool ifms = true;

  // ###########################################################################
  // INITIALIZATION
  // -> Determine basic values of structure and geometry
  // ###########################################################################

  // determine basic numbers for the structure
  nRows = testspace->get_n_dof();
  nColumns = ansatzspace->get_n_dof();
  int n_hanging_rows = testspace->get_n_hanging();
  rows = std::vector<int>(nRows+1, 0);
  nEntries = 0;

  if (is_empty)
  { //no need to create anything else...
    return;
  }

  // get collection of mesh cells
  auto coll = testspace->GetCollection();

  // get information about the spaces from the possibly different ansatz and
  // test space
  int n_dirichlet_test = testspace->get_n_dirichlet();  // number Dirichlet nodes
  unsigned int n_active_rows = testspace->get_n_active_non_hanging();
  auto hanging_bound = n_active_rows + n_hanging_rows;

  // global numbers of degrees of freedom and information about hanging nodes
  auto n_hanging_test = testspace->get_n_hanging();

  // upper_bounds[i] will contain the maximal number of couplings counted with
  // repetition for test dof i, i.e. the number of all ansatz dofs of cells to
  // which the particular dof i belongs to. It is therefore an upper bound for
  // the number of entries in row i.
  // If additional_face_integrals also the number of ansatz dofs of
  // neighbour cells of cells to which dof belongs to.
  // This array is needed to determine the indices in the following array
  std::vector<int> upper_bounds(nRows + 1, 0);

  // get information about the boundary nodes
  auto n_active_rows_orig = n_active_rows;
  if (ifms)
  { // internal full matrix structure
    // if ifms it can be seen as if there were no Dirichlet or hanging nodes,
    // since they now couple with all of their coupling partners, even if most
    // of the entries are 0.
    n_dirichlet_test = 0;
    n_hanging_test = 0;
    n_active_rows = nRows;
  }

  // ###########################################################################
  // Fill upper_bounds
  // -> Find out maximal possible number of couplings
  // ###########################################################################

  // Now follows a loop over all cells in which upper_bounds is filled. After
  // the loop it will contain the maximal number of entries per row. This is
  // done as follows: On each cell there is a loop over all test dofs in which
  // the number of ansatz dofs are added to the test dofs' position in
  // upper_bounds. This is done due to the fact that at least every ansatz dof
  // couples with all test dofs on this cell.
  // For internal_face_integrals an additional loop over the neighbour cells is
  // done in which also the number of the test dofs of the neighbour cells are
  // added.
  auto n_cells = coll->GetN_Cells();
  for(int cell_i = 0; cell_i < n_cells; cell_i++)
  {
    // global numbers of ansatz and test dofs
    auto global_dof_numbers_t = testspace->GetGlobalDOF(cell_i);
    // number local ansatz and dofs
    auto n_dofs_a = ansatzspace->get_n_local_dof(cell_i);
    auto n_dofs_t = testspace->get_n_local_dof(cell_i);

    auto n_far_dof
     = face_integrals ? get_n_far_dof(coll->GetCell(cell_i), *ansatzspace) : 0u;

    // loop over test dofs
    for(unsigned int t_dof_j = 0; t_dof_j < n_dofs_t; t_dof_j++)
    { // each row corresponds to a test dof
      unsigned int glob_num_t = global_dof_numbers_t[t_dof_j];

      // at least all the dofs on a particular cell couple, and therefore
      // produce an entry
      if(glob_num_t < n_active_rows)
      { // "real" dof not Dirichlet / Hanging
        upper_bounds[glob_num_t] += n_dofs_a;
      }
      if(glob_num_t >= n_active_rows_orig && glob_num_t < hanging_bound)
      { // hanging dof
        add_far_hanging_coupling(cell_i, glob_num_t, *testspace, *ansatzspace,
                                 upper_bounds);
      }
      // additional entries necessary if integrals over the edges
      // appear in the discretization. If so, there is also a coupling with
      // the ansatz dofs of the neighbour cell
      upper_bounds[glob_num_t] += n_far_dof;
    } // endfor t_dof_j
  } // endfor cell_i


  // If not Internal Full Matrix, up to now the hanging nodes and the
  // Dirichlet nodes were neglected, since they are completely determined due
  // to the boundary data or the conformity of the solution. Therefore, they
  // have to be considered now.
  if (!ifms)
  {
    // add couplings for hanging nodes of space
    // each hanging node produces one entry in the diagonal + couplings to still
    // have a conformal solution
    for(size_t hn_i = 0; hn_i < n_hanging_test; hn_i++)
    {
      auto hn = testspace->get_hanging_node(hn_i);  // current hanging node
      int index = n_active_rows + hn_i; // hanging nodes after Neumann nodes
      auto n = hn->GetN_Nodes() + 1;  // number of couplings + diagonal
      upper_bounds[index] = n;
      nEntries += n;
    }

    // add rows for Dirichlet nodes in  space
    // each Dirichlet node produces exactly one entry on the diagonal
    nEntries += n_dirichlet_test;
    for( int i = 0; i < n_dirichlet_test; i++)
    {
      // Dirichlet nodes come after active dof and hanging dofs
      int index = n_active_rows + n_hanging_test + i;
      upper_bounds[index] = 1;
    }
  }

  // Cumulatively sum up the array upper_bounds. upper_bounds[i] will then
  // contain the index for the upcoming array that specifies the start of the
  // interval in the array related to test dof i.
  auto old_value = upper_bounds[0];
  upper_bounds[0] = 0;
  for(unsigned int dof_i = 0; dof_i < nRows; dof_i++)
  {
    auto new_value = upper_bounds[dof_i+1];
    upper_bounds[dof_i+1] = upper_bounds[dof_i] + old_value;
    old_value = new_value;
  }

  // ###########################################################################
  // Fill columns_aux
  // -> Find out correct number of couplings and total number of entries
  // ###########################################################################

  // columns_aux will contain the correct values for columns and possibly many
  // more -1 that will originate from possibly multiple counting of couplings.
  // The entries for row/test dof i are in the interval
  // [columns_aux[upper_bounds[i]], columns_aux[upper_bounds[i+1]])
  // With the help of this array the actual arrays of interest namely rows and
  // columns, and nEntries will be constructed.
  std::vector<int> columns_aux(upper_bounds[nRows], -1);// initialize it with -1

  // In whats following the array columns_aux is filled and the total number of
  // entries is  counted. In columns_aux there will be packages/intervals of
  // size upper_bounds[1]-upper_bounds[0], upper_bounds[2]-upper_bounds[1], ...,
  // upper_bounds[end]-upper_bounds[end-1], filled with all column  numbers for
  // the particular row without repetitions. The possible repetitions stay at
  // -1.
  // Inside a loop over all cells is an inner loop over the test
  // dofs of the cell, and inside this loop there is a loop over the
  // ansatz dofs of the cell. Altogether on each cell every combination between
  // ansatz and test dof is checked. In other words, all possible columns are
  // checked for each row. If the column is new, there is a new entry leading to
  // an increase of nEntries and to the fact that columns_aux is set to the
  // column number at the correct position. Else columns_aux stays at -1.
  // For additional face integrals also a loop over the neighbours is performed.
  // This will be used afterwords to construct the column and the row array.

  for(int cell_i = 0; cell_i < n_cells; cell_i++)
  {
    // global numbers of ansatz and test dofs
    auto global_dof_numbers_a = ansatzspace->GetGlobalDOF(cell_i);
    auto global_dof_numbers_t = testspace->GetGlobalDOF(cell_i);
    // number local ansatz and dofs
    auto n_dofs_a = ansatzspace->get_n_local_dof(cell_i);
    auto n_dofs_t = testspace->get_n_local_dof(cell_i);

    for(unsigned int dof_j = 0; dof_j < n_dofs_t; dof_j++)
    { // test dofs
      unsigned int glob_num_t = global_dof_numbers_t[dof_j]; //global dof number

      for(unsigned int dof_k = 0; dof_k < n_dofs_a; dof_k++)
      { // ansatz dofs
        auto glob_num_a = global_dof_numbers_a[dof_k];

        if(glob_num_t < n_active_rows)
        { // Neumann or inner node
          add_entry(upper_bounds, columns_aux, glob_num_t, glob_num_a,
                    nEntries);
        }
      } // endfor dof_k

      if(glob_num_t >= n_active_rows_orig && glob_num_t < hanging_bound)
      {
        add_columns_far_hanging_coupling(cell_i, glob_num_t, *testspace,
                                         *ansatzspace, upper_bounds,
                                         columns_aux, nEntries);
      }
      // additional entries necessary if integrals over the edges
      // appear in the discretization
      if (face_integrals)
      {
        auto cell = coll->GetCell(cell_i); // current cell with number cell_i
        int n_faces = cell->GetN_Joints(); // number of neighbours
        for(int joint_l = 0; joint_l < n_faces; joint_l++)
        { // neighbour cells
          auto neigh = cell->GetJoint(joint_l)->GetNeighbour(cell);
          if(neigh)
          {
            auto neigh_id = neigh->GetClipBoard();
            if (neigh_id < 0)
            {
              continue;
            }
            auto NumbersNeighbour = ansatzspace->GetGlobalDOF(neigh_id);

            auto n_neigh_dof_a = ansatzspace->get_n_local_dof(neigh_id);
            for(unsigned int neigh_local_dof_id = 0; neigh_local_dof_id <
                n_neigh_dof_a; neigh_local_dof_id++)
            {
              auto current_neigh_nr = NumbersNeighbour[neigh_local_dof_id];
              if(glob_num_t < n_active_rows)
              { // this dof is a real node (inner or Neumann)
                add_entry(upper_bounds, columns_aux, glob_num_t,
                          current_neigh_nr, nEntries);
              } // endif real dof
            } // endfor neigh_local_dof_id
          } // endif neigh
        } // endfor p
      } // endif ifms

    } // endfor dof_j
  } // endfor cell_i


  // If not Internal Full Matrix the hanging nodes are not yet
  // included in columns_aux. Hence, this has to be treated now
  if (!ifms)
  {
    // add the additional columns from hanging nodes to other nodes
    for(size_t hn_i = 0; hn_i < n_hanging_test; hn_i++)
    {
      auto hn = testspace->get_hanging_node(hn_i);
      auto coupling_dofs = hn->GetDOF();
      auto index = upper_bounds[n_active_rows + hn_i];
      auto n_couplings = hn->GetN_Nodes();
      columns_aux[index] = n_active_rows_orig + hn_i ;

      for(int dof_j = 0; dof_j < n_couplings; dof_j++)
      { // loop over all nodes in coupling
        auto glob_dof_num = coupling_dofs[dof_j];
        columns_aux[index + dof_j + 1] = glob_dof_num;
      } // endfor dof_j
    } // endfor hn_i
  } // endif ifms

  // ###########################################################################
  // Put everything together
  // -> construct rows, columns, hanging rows, hanging columns
  // ###########################################################################

  // compress columns_aux array to KCol by deleting all -1's
  // build the rows array
  columns.resize(nEntries);
  rows = upper_bounds;

  // Construct columns by a loop over all elements in columns_aux and if they
  // are not -1 add them to rows. While doing so, count the adding and
  // cumulative add the numbers to rows.
  // Analogously for hanging columns afterwards.
  int index_entry = 0; // index of current entry, counts the number of entries
  for(unsigned int t_dof_i = 0; t_dof_i < n_active_rows + n_hanging_test;
      t_dof_i++)
  {
    for(int j = upper_bounds[t_dof_i];
        j < upper_bounds[t_dof_i+1] && columns_aux[j] != -1;
        j++)
    {
      columns[index_entry] = columns_aux[j];
      index_entry++;
    } // endfor j
    rows[t_dof_i + 1] = index_entry;
  } // endfor t_dof_i

  // If not Internal Full Matrix Structure, the Dirichlet nodes has to be
  // treated as well, since they were neglected up to now.
  if (!ifms)
  {
    auto offset = rows[hanging_bound];
    for(int diri_i = 0; diri_i < n_dirichlet_test; diri_i++)
    {
      int index = hanging_bound + diri_i + 1;
      rows[index + 1] = rows[index] + 1 ;
      columns[offset + diri_i] = hanging_bound + diri_i ;
    }
  }

  this->Sort();
}

TStructure::TStructure()
: nRows(0), nColumns(0), nEntries(0), columns(), rows()
{
}

TStructure::TStructure(int n, int nEntries, int *col_ptr, int *row_ptr)
  : TStructure(n, n, nEntries, col_ptr, row_ptr)
{
}

TStructure::TStructure(int nRows, int nCols, int nEntries, int *col_ptr,
                       int *row_ptr)
 : nRows(nRows), nColumns(nCols), nEntries(nEntries), columns(nEntries, 0),
   rows(nRows+1)
{
  std::copy(col_ptr, col_ptr + this->get_n_entries(), this->columns.begin());
  std::copy(row_ptr, row_ptr + this->get_n_rows() + 1, this->rows.begin());
  this->Sort();
}

TStructure::TStructure(int nRows, int nCols) : nRows(nRows), nColumns(nCols),
  nEntries(0), columns(), rows(nRows+1, 0)
{
}

TStructure::TStructure(int n) : TStructure(n, n)
{
}

/* sort one row [BeginPtr, AfterEndPtr) */
void TStructure::SortRow(int *BeginPtr, int *AfterEndPtr)
{
  int *IPtr, *JPtr, T;

  for(IPtr=BeginPtr;IPtr<AfterEndPtr;IPtr++)
  {
    for(JPtr=IPtr+1;JPtr<AfterEndPtr;JPtr++)
    {
      if( *IPtr > *JPtr )
      {
        T = *IPtr;
        *IPtr = *JPtr;
        *JPtr = T;
      }
    } // endfor JPtr
  } // endfor IPtr
}

/* sort numbers within each row */
void TStructure::Sort()
{
  int end, begin;

  end = 0;
  for(unsigned int i=0; i<nRows; i++)
  {
    begin = end;
    end = rows[i+1];
    SortRow(&columns[0]+begin, &columns[0]+end);
  } // endfor i
}

size_t TStructure::get_n_entries_in_row(size_t row_index) const
{
  if((unsigned int) row_index > nRows)
  {
    ErrThrow("unable to return the number of entries in row ", row_index,
        ". There are only ", this->nRows, " rows.");
  }
  auto n = this->rows[row_index+1]-this->rows[row_index];
  if(n < 0)
    ErrThrow("it seems this structure is in an invalid state");
  return (size_t)n;
}

size_t TStructure::get_n_entries_up_to_row(size_t row_index) const
{
  if((unsigned int) row_index > nRows)
  {
    ErrThrow("unable to return the number of entries up to row ", row_index,
        ". There are only ", this->nRows, " rows.");
  }
  auto n = this->rows[row_index];
  if(n < 0)
    ErrThrow("it seems this structure is in an invalid state");
  return (size_t)n;
}

void TStructure::reset_n_entries()
{
  //throw if number of rows changed
  if (this->rows.size() - 1 != this->nRows)
    ErrThrow("TStructure: nRows != rows.size() - 1 ",nRows, " != ", rows.size() + 1);

  //reset number of entries to last entry of row ptr
  this->nEntries = this->rows.back();

  //columns array must be resized  to nEntries - everything behind that is erased
  this->columns.resize(this->nEntries);
}

int TStructure::get_index_of_entry(const int i, const int j) const
{
  if(i < 0 || i >= (int) this->get_n_rows())
  {
    ErrThrow("row index is too large or too small");
  }
  if(j < 0 || j > (int) this->get_n_columns())
  {
    ErrThrow("column index is too large or too small");
  }

  for (int m=rows[i];m < rows[i+1]; m++)
  {
    if (columns[m]== j)
    {
      // index found in sparsity pattern
      return m;
    }
  }
  // index not in the sparsity pattern
  return -1;
}

/* return a new structure for a transposed matrix */
std::shared_ptr<TStructure> TStructure::get_transposed() const
{
  // variables for transposed structure:
  int nRowsT = nColumns;
  int nColsT = nRows;
  // number of entries does not change
  std::vector<int> rowsT(nRowsT+1, 0);
  std::vector<int> colsT(nEntries, 0);

  std::vector<int> ColB_count(nColumns, 0);

  // count number of entries per column (will be number of entries in each row)
  for(int i=0;i<rows[nRows];i++)
    rowsT[columns[i]]++;
  // change to increasing numbering as in rows
  for(int i=0,k=0;i<=nRowsT;i++)
  {
    int j = rowsT[i];
    rowsT[i] = k;
    k += j;
  }

  // fill 'colsT'
  // loop over (non-transposed) rows
  for(unsigned int i=0; i<nRows; i++)
  {
    // loop over all entries in this row
    for(int j=rows[i]; j<rows[i+1]; j++)
    {
      int col = columns[j]; // (non-transposed) column = transposed row
      int l  = rowsT[col];
      int offset = ColB_count[col];
      colsT[l+offset] = i;
      ColB_count[col]++;
    }
  }
  std::shared_ptr<TStructure> structureT(new TStructure(nRowsT, nColsT,
        nEntries, colsT.data(), rowsT.data()));
  return structureT;
}

std::shared_ptr<TStructure> get_product_structure(
    TStructure const & strucA, TStructure const & strucB)
{
  const int n_A_rows = strucA.get_n_rows();   // = n_C_rows
  const int n_A_cols = strucA.get_n_columns();
  const int n_B_rows = strucB.get_n_rows();
  const int n_B_cols = strucB.get_n_columns();   // = n_C_cols

  if(n_A_cols != n_B_rows)
  {
    ErrThrow("dimension mismatch during matrix-matrix multiplication");
  }
  const int * const a_rows = strucA.get_row_ptr();
  const int * const a_cols = strucA.get_vector_columns();

  // everything needed to call the constructor of TStructure later on:
  int n_c_entries = 0; // number of entries in product structure C
  std::vector<int> c_rows(n_A_rows + 1, 0.0); // row pointer
  std::vector<std::vector<int> > dofs(n_A_rows);

  // loop over all rows of C
  for(int row = 0; row < n_A_rows; row++)
  {
    // loop over all columns of C
    for(int col = 0; col < n_B_cols; col++)
    {
      // check whether 'this row of A' x 'this column of B' would give an entry
      int n_a_entries_in_row = a_rows[row+1] - a_rows[row];
      // loop over all entries in this row in A
      for(int i = 0; i < n_a_entries_in_row; i++)
      {
        if(strucB.get_index_of_entry(a_cols[i+a_rows[row]], col) != -1)
        {
          dofs[row].push_back(col);
          break;
        }
      }
    }
    n_c_entries += dofs[row].size();
    c_rows[row+1] = n_c_entries; // c_rows[0] has been set to 0 already
  }

  std::vector<int> c_cols(n_c_entries); // columns pointer
  // now fill the array c_cols
  // loop over all rows of C
  for(int row = 0; row < n_A_rows; row++)
  {
    // loop over all columns of C
    int nEntries_in_this_row = c_rows[row+1] - c_rows[row];
    for(int col = 0; col < nEntries_in_this_row; col++)
    {
      c_cols[c_rows[row] + col] = dofs[row].at(col);
    }
  }  
  return std::make_shared<TStructure>(n_A_rows, n_B_cols, n_c_entries,
                                      c_cols.data(), c_rows.data());
}


TStructure* TStructure::get_structure_of_product_with_transpose_from_right() 
  const
{
  int nProductEntries = 0; // number of entries in product structure
  int nProductRows = nRows;
  int nProductColumns = nRows;

  std::vector<int> productRowPtr(nProductRows + 1, 0); // row pointer

  // this is a temporary storing structure "gridPlaces"
  std::vector<std::vector<int> > gridPlaces(nRows); 
  // gridPlaces[row] stores, in which columns in this row there are non-zero 
  // entries

  // loop over all rows of the product
  for(int row = 0; row < nProductRows; row++)
  {
    // loop over all columns of the product
    for(int col = 0; col < nProductColumns; col++)
    {
      // check whether 'this row of A' x 'this column of A^T' would give an 
      // entry
      // this boils down to checking 'row1 of A' x 'row2 of A' with row1 = 
      // row, row2 = col
      int row1 = row;
      int row2 = col; //two definitions to fix ideas

      // work on two segments of column array
      const int* row1ColBegin = &columns[rows[row1]];
      const int row1ColSize = rows[row1+1] - rows[row1];
      const int* row2ColBegin = &columns[rows[row2]];
      const int row2ColSize = rows[row2+1] - rows[row2];

      // find out if the two arrays contain a common value
      // exploit the fact that both are sorted
      size_t index1 = 0;
      size_t index2 = 0;
      while ((int)index1 < row1ColSize && (int)index2 < row2ColSize)
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
          //we found a pair of indices with equal entry in KCol
          gridPlaces[row].push_back(col);
          break;
        }
      }

    } //end for columns of the product

    // update the total number of entries in the product
    nProductEntries += gridPlaces[row].size();
    // productRowPtr[0] has been set to 0 already
    productRowPtr[row+1] = nProductEntries;
  }//end for rows of the product

  // FROM HERE IT'S ONLY TRANSFORMING THE TEMPORARY DATA STRUCTURE TO THE 
  // REQUIRED ONE
  // now fill the array productColumnsPtr
  //initialise the product's columns pointer
  std::vector<int> productColumnPtr(nProductEntries, 0);

  // loop over all rows of C
  for(int row = 0; row < nProductRows; row++)
  {
    // loop over all columns of C
    int n_entries_in_this_row = productRowPtr[row+1] - productRowPtr[row];
    for(int col = 0; col < n_entries_in_this_row; col++)
    {
      productColumnPtr[productRowPtr[row] + col] = gridPlaces[row].at(col);
    }
  }
  // hand over a pointer to the product's structure.
  return new TStructure(nProductRows, nProductColumns, nProductEntries, 
      &productColumnPtr[0], &productRowPtr[0]);
}


TStructure* TStructure::get_structure_of_product_with_transpose_from_right(
    const TStructure& B) const
{
  if(this->nColumns != B.get_n_rows() || B.get_n_columns() != this->nColumns)
    ErrThrow("dimension mismatch, inner matrix has wrong dimensions, ",
        B.get_n_rows(), "  ", B.get_n_columns());

  int nProductEntries = 0; // number of entries in product structure
  int nProductRows = this->nRows;
  int nProductColumns = this->nRows;

  // lambda function to check if the product of 'B*A^T' has an entry (i,j)
  auto check_entry = [this, B](int i, int j) -> bool
  {
    // work on two segments of column array
    const int* row1ColBegin = &B.get_vector_columns()[B.get_row_ptr()[i]];
    const int row1ColSize = B.get_row_ptr()[i+1] - B.get_row_ptr()[i];
    const int* row2ColBegin=&this->get_vector_columns()[this->get_row_ptr()[j]];
    const int row2ColSize = this->get_row_ptr()[j+1]
      - this->get_row_ptr()[j];

    // find out if the two arrays contain a common value
    // exploit the fact that both are sorted
    size_t index1 = 0;
    size_t index2 = 0;
    while ((int)index1 < row1ColSize && (int)index2 < row2ColSize)
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
        return true;
      }
    }
    return false;
  };

  std::vector<int> productRowPtr(nProductRows + 1, 0); // row pointer
  // this is a temporary storing structure "gridPlaces"
  std::vector<std::vector<int>> gridPlaces(this->nRows);
  // gridPlaces[row] stores, in which columns in this row there are non-zero 
  // entries
  // loop over all rows of the product
  for(int row = 0; row < nProductRows; row++)
  {
    // loop over all columns of the product
    for(int col = 0; col < nProductColumns; col++)
    {
      //Output::print("row ", row, "  column ", col);
      // check whether 'this row of A' x 'this column of B*A^T' would give an 
      // entry
      // this boils down to checking 'row1 of A' x 'row2 of B*A' with row1 = 
      // row, row2 = col
      int row1 = row;
      //int row2 = col; //two definitions to fix ideas

      // work on two segments of column array
      const int* row1ColBegin = &columns[rows[row1]];
      const int row1ColSize = rows[row1+1] - rows[row1];
      //const int* row2ColBegin = &columns[rows[row2]];
      //const int row2ColSize = rows[row2+1] - rows[row2];

      // find out if the two arrays contain a common value
      // exploit the fact that both are sorted
      size_t index1 = 0;
      size_t index2 = 0;
      while ((int)index1 < row1ColSize && (unsigned int)index2 < this->nColumns)
      {
        //Output::print("  index1 ", index1, "  index2 ", index2);
        if (row1ColBegin[index1] > (int)index2)
        {
          index2++;
        }
        else if ((int)index2 > row1ColBegin[index1])
        {
          index1++;
          // I think we should never be here (Ulrich)
        }
        else
        {
          // now we need to check if 'B*A^T' has an entry at (index2,col)
          if(check_entry(index2, col))
          {
            //we found a pair of indices with equal entry in KCol
            gridPlaces[row].push_back(col);
            break;
          }
          ++index1;
          ++index2;
        }
      }

    } //end for columns of the product

    // update the total number of entries in the product
    nProductEntries += gridPlaces[row].size();
    // productRowPtr[0] has been set to 0 already
    productRowPtr[row+1] = nProductEntries;
  }//end for rows of the product

  // FROM HERE IT'S ONLY TRANSFORMING THE TEMPORARY DATA STRUCTURE TO THE 
  // REQUIRED ONE
  // now fill the array productColumnsPtr
  //initialise the product's columns pointer
  std::vector<int> productColumnPtr(nProductEntries, 0);

  // loop over all rows of C
  for(int row = 0; row < nProductRows; row++)
  {
    // loop over all columns of C
    int n_entries_in_this_row = productRowPtr[row+1] - productRowPtr[row];
    for(int col = 0; col < n_entries_in_this_row; col++)
    {
      productColumnPtr[productRowPtr[row] + col] = gridPlaces[row].at(col);
    }
  }
  // hand over a pointer to the product's structure.
  return new TStructure(nProductRows, nProductColumns, nProductEntries, 
      &productColumnPtr[0], &productRowPtr[0]);
}


void TStructure::fortran_shift()
{
  // the first entry in the this->rows is always zero (when indices start with 
  // zero), so it is used as an indicator if the entries in this->rows and 
  // this->columns are in Fortran or c++ style.
  int s = 1; // forward shift
  if(this->is_fortran_shifted())
  {
    // backward shift (indices start with 1)
    s = -1;
  }
  // else // forward shift (indices start with 0)
  std::for_each(this->rows.begin(), this->rows.end(),       [s](int& i){i+=s;});
  std::for_each(this->columns.begin(), this->columns.end(), [s](int& i){i+=s;});
}

bool TStructure::is_fortran_shifted() const
{
  if(this->rows[0] == 1)
  {
    return true;
  }
  else if(this->rows[0] == 0)
  {
    return false;
  }
  else
    ErrThrow("broken matrix structure, first index in row vector is neither ",
        "0 nor 1, but ", this->rows[0]);
}

void TStructure::info() const
{
  Output::info<3>("TStructure","Information on the stored matrix structure");
  Output::dash<3>("Number of rows: ", nRows);
  Output::dash<3>("Number of columns: ", nColumns);
  Output::dash<3>("Number of matrix entries: ", nEntries);
}

void TStructure::draw(const std::string& filename) const
{
  std::ofstream out_stream(filename);
  if(!out_stream)
  {
    Output::warn("TStructure::draw", "cannot open '", filename, "' for output");
    return;
  }

  double scale = 3;
  int BX = (int) (scale * this->nColumns); // width of picture
  int BY = (int) (scale * this->nRows);    // height of picture
  int offset = 2 * scale; // picture is a bit away from picture boundary

  out_stream << "%!PS-Adobe-3.0\n";
  out_stream << "%%Creator: ParMooN (Ulrich Wilbrandt)\n";
  out_stream << "%%DocumentFonts: Helvetica\n";
  out_stream << "%%BoundingBox: 0 0 " << 2*offset + BX << " " << 2*offset + BY;
  out_stream << endl;
  out_stream << "%%Pages: 1\n";
  out_stream << "%%EndComments\n";
  out_stream << "%%EndProlog\n";
  out_stream << "%%Page: 1 1\n";
  out_stream << "%% n_rows " << this->nRows << ",  n_columns " 
    << this->nColumns << ",  n_entries " << this->nEntries << endl;
  out_stream << "/M {" << scale-1 << " " << scale-1 << " rectfill} def\n";
  for(unsigned int row = 0; row < this->nRows; ++row)
  {
    for(int col = this->rows[row]; col < this->rows[row+1]; ++col)
    {
      out_stream << (offset + scale * this->columns[col]) << " " 
        << (offset + scale * (this->nRows - row)) << " M\n";
    }
  }
  out_stream << "stroke" << endl;
  out_stream << "showpage" << endl;
  out_stream << "%%Trailer" << endl;
  out_stream << "%%Pages: 1" << endl;
  out_stream.close();
  Output::print<2>("postscript picture of matrix structure drawn in file ",
      filename);
}


bool operator==(const TStructure &lhs, const TStructure &rhs)
{
  return lhs.nRows == rhs.nRows
    && lhs.nColumns == rhs.nColumns
    && lhs.nEntries == rhs.nEntries
    && lhs.rows == rhs.rows
    && lhs.columns == rhs.columns;
}

bool operator!=(const TStructure &lhs, const TStructure &rhs)
{
  return !(rhs == lhs);
}

