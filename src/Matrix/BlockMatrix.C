#include <BlockMatrix.h>
#include <BlockVector.h>
#include <LinAlg.h>
#include "MooNMD_Io.h"
#include <algorithm>
#include <limits>
#include <list>

#ifdef _MPI
#include <BlockFEMatrix.h> //needed in sor_sweep for a dynamic type downcast
#endif

/* ************************************************************************* */
BlockMatrix::CellInfo::CellInfo()
: BlockMatrix::CellInfo::CellInfo(0, 0) //delegate construction
{
}

/* ************************************************************************* */
BlockMatrix::CellInfo::CellInfo(size_t nRows, size_t nColumns)
: n_rows_{nRows}, n_columns_{nColumns},
  //block_in_cell gets default initialised to null
  color_{std::numeric_limits<size_t>::max()}, // no colour
  is_transposed_{false}, //non-transposed
  re_color_flag_{ReColoringFlag::KEEP}

  {
  }

  /* ************************************************************************* */
  // IMPLEMENTATION OF PUBLIC METHODS
  /* ************************************************************************* */

  /* ************************************************************************* */

  BlockMatrix::BlockMatrix()
  : BlockMatrix{{},{}}
  {
  }

  // throws an error
  BlockMatrix::BlockMatrix(CompositeOperator<BlockVector>)
  : BlockMatrix()
  {
    ErrThrow("Cannot convert a CompositeOperator to a BlockMatrix!");
  }

  BlockMatrix::BlockMatrix(const std::vector<size_t>& cellRowNumbers,
                           const std::vector<size_t>& cellColumnNumbers)
  : n_cell_rows_{cellRowNumbers.size()},
    n_cell_columns_{cellColumnNumbers.size()},
    cell_grid_(n_cell_rows_, std::vector<CellInfo>(n_cell_columns_))
    // color_count_ gets default initialized as vector of size zero
    // combined_matrix gets default initialized to smart nullptr
    {
      //traverse the cell info grid from left to right, top to bottom
      // and fill in newly constructed, correctly dimensioned zero matrices as blocks
      for(size_t i = 0; i < n_cell_rows_ ; ++i )
      {
        //hold the number of rows each cell in this row will have
        size_t nRowsOfCell = cellRowNumbers[i];

        for(size_t j = 0; j < n_cell_columns_ ; ++j )
        {
          //hold the number of columns each cell in this column will have
          size_t nColumnsOfCell = cellColumnNumbers[j];

          // construct the new cell info and the zero matrix it will hold
          CellInfo newInfo(nRowsOfCell, nColumnsOfCell);
          newInfo.block_ = std::make_shared<TMatrix>(nRowsOfCell, nColumnsOfCell);

          // the next color is get_n_colors() (because the colors 0 to
          // get_n_colors() -1 are already assigned)
          newInfo.color_ = get_n_colors();

          //start as non-transposed
          newInfo.is_transposed_ = false;

          // put the new cell info to the correct place by copy assignment
          cell_grid_ [i][j] = newInfo;

          // one new block has been added and colored
          color_count_.push_back(1);

        }
      }
    }

/* ************************************************************************* */
BlockMatrix::BlockMatrix(int nRows, int nCols,
                         const std::vector<std::shared_ptr<TMatrix>>& blocks)
 : n_cell_rows_(nRows), n_cell_columns_(nCols),
   cell_grid_(nRows, std::vector<CellInfo>(nCols)), color_count_(nRows*nCols, 1)
{
  // make sure enough blocks are provided
  if(blocks.size() < n_cell_rows_ * n_cell_columns_)
  {
    ErrThrow("unable to create a ", n_cell_rows_, "x", n_cell_columns_,
             " BlockMatrix out of only ", blocks.size(), " blocks");
  }
  // check if all given blocks exist (none of the pointers are to nullptr)
  for(auto& m : blocks)
  {
    if(!m)
    {
      // we might be able to create a zero matrix of the right size, but most
      // likely this is an error
      ErrThrow("unable to create BlockMatrix with one block given as nullptr");
    }
  }

  for(size_t i = 0; i < n_cell_rows_ ; ++i )
  {
    //hold the number of rows each cell in this row will have
    size_t nRowsOfCell = blocks.at(i * this->n_cell_columns_)->get_n_rows();
    for(size_t j = 0; j < n_cell_columns_ ; ++j )
    {
      //hold the number of columns each cell in this column will have
      size_t nColumnsOfCell = blocks.at(j)->get_n_columns();

      // construct the new cell info
      CellInfo newInfo(nRowsOfCell, nColumnsOfCell);

      // add the matrix
      newInfo.block_ = blocks.at(i * this->n_cell_columns_ + j);

      // check if the block has the correct dimensions, basically this checks if
      // all blocks in one row/column have the same number of rows/columns.
      if(  newInfo.block_->get_n_rows() != (int) nRowsOfCell
        || newInfo.block_->get_n_columns() != (int) nColumnsOfCell)
      {
        ErrThrow("unable to create a block matrix at entry ", i, " ", j,
                 ". The matrix block should have dimensions ", nRowsOfCell,
                 "x", nColumnsOfCell, ", but has ", newInfo.block_->get_n_rows(),
                 "x", newInfo.block_->get_n_columns());
      }

      // each block has its own color, ordered row wise
      newInfo.color_ = i*n_cell_columns_ + j;

      //start as non-transposed
      newInfo.is_transposed_ = false;

      // put the new cell info to the correct place by copy assignment
      cell_grid_[i][j] = newInfo;
    }
  }

  // perform a few checks to make sure this matrix is properly defined
  this->check_coloring();
}


    /* ************************************************************************* */
    void BlockMatrix::add_matrix_to_blocks(const TMatrix& summand,
                                                  std::vector<grid_place_and_mode> row_column_transpose_tuples)
    {
      add_scaled_matrix_to_blocks(summand, 1.0, row_column_transpose_tuples);
    }

    /* ************************************************************************* */
    void BlockMatrix::add_unscaled_matrix(
        const TMatrix& summand,
        const std::vector<std::vector<size_t>>& cell_positions,
        const std::vector<bool>& transposed_states)
    {
      add_matrix(summand, 1.0, cell_positions, transposed_states);
    }



    /* ************************************************************************* */
    void BlockMatrix::add_matrix(
        const TMatrix& summand, double scaling_factor,
        const std::vector<std::vector<size_t>>& cell_positions,
        const std::vector<bool>& transposed_states)
    {
      std::vector<std::tuple<size_t, size_t, bool>> input_tuples(
          check_and_tupelize_vector_input(cell_positions, transposed_states));


      add_scaled_matrix_to_blocks(summand, scaling_factor, input_tuples);

    }

    /* ************************************************************************* */
    void BlockMatrix::apply(const BlockVector & x, BlockVector & y) const
    {
      // reset all values in 'y' to 0 and delegate to apply_scaled_add
      y.reset();
      apply_scaled_add(x, y, 1.0);
    }

    /* ************************************************************************* */
    void BlockMatrix::apply_scaled_add(const BlockVector & x, BlockVector & y,
                                              double a) const
    {
      // check if the vectors fit, if not so the program throws an error
      check_vector_fits_pre_image(x);
      check_vector_fits_image(y);

      const double * xv = x.get_entries(); // array of values in x
      double * yv = y.get_entries(); // array of values in y

      size_t row_offset = 0;

      // n_rows, n_cols are the number of cell rows/columns
      for (size_t i = 0; i < n_cell_rows_; ++i)
      {
        int col_offset = 0;

        for (size_t j = 0; j < n_cell_columns_; j++)
        {
          std::shared_ptr<TMatrix> current_block = cell_grid_[i][j].block_;
          bool transp_state = cell_grid_[i][j].is_transposed_;

          if (transp_state)
          {
            current_block->transpose_multiply(xv + col_offset, yv + row_offset, a);
          }
          else
          {
            current_block->multiply(xv + col_offset, yv + row_offset, a);
          }

          col_offset += cell_grid_[i][j].n_columns_;
        }

        row_offset += cell_grid_[i][0].n_rows_;
      }
    }

/* ************************************************************************* */
void BlockMatrix::apply_transpose(const BlockVector & x, BlockVector & y) const
{
  // reset all values in 'y' to 0 and delegate to apply_scaled_add
  y.reset();
  apply_transpose_scaled_add(x, y, 1.0);
}

void BlockMatrix::apply_transpose_scaled_add(const BlockVector & x, BlockVector & y,
                                              double a) const
{
  // check if the vectors fit, if not so the program throws an error
  check_vector_fits_image(x);
  check_vector_fits_pre_image(y);

  const double * xv = x.get_entries(); // array of values in x
  double * yv = y.get_entries(); // array of values in y

  size_t row_offset = 0;

  // n_rows, n_cols are the number of cell rows/columns
  for (size_t i = 0; i < n_cell_rows_; ++i)
  {
    int col_offset = 0;

    for (size_t j = 0; j < n_cell_columns_; j++)
    {
      std::shared_ptr<TMatrix> current_block = cell_grid_[i][j].block_;
      bool transp_state = cell_grid_[i][j].is_transposed_;

      if (!transp_state)
      {
        current_block->transpose_multiply(xv + row_offset, yv + col_offset, a);
      }
      else
      {
        current_block->multiply(xv + row_offset, yv + col_offset, a);
      }

      col_offset += cell_grid_[i][j].n_columns_;
    }

    row_offset += cell_grid_[i][0].n_rows_;
  }
}

/* ************************************************************************* */
#ifdef _MPI
void BlockMatrix::sor_sweep(const BlockVector& b, BlockVector& x, double omega,
                            size_t flag, const std::string& par_strat) const
#else
void BlockMatrix::sor_sweep(const BlockVector& b, BlockVector& x, double omega,
                            size_t flag) const
#endif
{
  if (flag > 2)
  {
    ErrThrow("BlockMatrix::sor_sweep with flag not 0,1, or 2.");
  }

  size_t n_diag_blocks = this->get_n_cell_rows();

  if (this->get_n_cell_columns() != n_diag_blocks)
  {
    ErrThrow("BlockMatrix::sor_sweep not tested for non square block "
             "structure");
  }

  // only for ssor we do two sweeps, first forward then backward.

  size_t n_sweeps = (flag == 2) ? 2 : 1;

  for (size_t sweep = 0; sweep < n_sweeps; ++sweep)
  {
    // do either a forward or a backward sweep

    bool backward_sweep = flag == 1 || (flag == 2 && sweep == 1);

    // loop over all diagonal blocks
    for (size_t i = 0; i < n_diag_blocks; ++i)
    {
      size_t d = i; // block on the diagonal

      if(backward_sweep)
      {
        d = n_diag_blocks - 1 - i; // reverse count
      }

      bool transposed;
      auto diag_block = this->get_block(d, d, transposed);

      if (transposed)
      {
        ErrThrow("Can't handle transposed blocks on the diagonal "
          "in BlockMatrix::sor_sweep");
      }

      if (!diag_block->is_square())
      {
        ErrThrow("unable to handle non-square diagonal block "
          "in BlockMatrix::sor_sweep");
      }

      size_t block_size = diag_block->get_n_rows();
      std::vector<double> modified_rhs(block_size, 0.0);

      std::copy(b.block(d), b.block(d) + block_size, modified_rhs.begin());

      // multiply solution with non-diagonal blocks in this block row.
      for (size_t row_block = 0; row_block < n_diag_blocks; ++row_block)
      {
        if(row_block == d)
        {
          continue; // skip diagonal block
        }

        auto non_diagonal_block = this->get_block(d, row_block, transposed);

        if (!transposed)
        {
          non_diagonal_block->multiply(
            x.block(row_block), &modified_rhs[0], -1.0);
        }
        else
        {
          non_diagonal_block->transpose_multiply(
            x.block(row_block), &modified_rhs[0], -1.0);
        }
      }

      // do the sor sweep

#ifdef _MPI
      const auto* block_fe_matrix = dynamic_cast<const BlockFEMatrix*>(this);

      if (block_fe_matrix == nullptr)
      {
        ErrThrow("Invoked BlockMatrix::sor_sweep on a non-BlockFEMatrix in "
          "MPI mode!");
      }

      auto comms = block_fe_matrix->get_communicators();

      if (backward_sweep)
      {
        diag_block->sor_sweep(&modified_rhs[0], x.block(d), omega, 1,
          par_strat, *comms.at(d));
      }
      else
      {
        diag_block->sor_sweep(&modified_rhs[0], x.block(d), omega, 0,
          par_strat, *comms.at(d));
      }

#else
      if (backward_sweep)
      {
        diag_block->sor_sweep(&modified_rhs[0], x.block(d), omega, 1);
      }
      else
      {
        diag_block->sor_sweep(&modified_rhs[0], x.block(d), omega, 0);
      }
#endif
    }
  }
}

    /* ************************************************************************* */
    void BlockMatrix::check_coloring() const
    {
      // check if the color_count_ array holds
      // the correct information
      check_color_count();

      // check the ascending order of the first appearing colors
      check_coloring_order();

      // check the equivalence of the relations "has the same color as"
      // and "holds a pointer to the same matrix as"
      check_equivalence_of_relations();

      // check if the recoloring flags are all in neutral state
      check_re_coloring_flags();

    }

    ///! Check whether a BlockVector b is fit to be the rhs b of the equation Ax=b.
    void BlockMatrix::check_vector_fits_image(const BlockVector& b) const
    {
      size_t n_vec_blocks = b.n_blocks();

      //check if number of blocks fits number of block rows
      if(n_vec_blocks != n_cell_rows_)
      {
        ErrThrow("Vector blocks number does not fit cell row number. ",
                 n_vec_blocks, "!=", n_cell_rows_);
      }

      //check if each vector block's length fits the number of rows in the corresponding row
      for(size_t i = 0; i<n_vec_blocks; ++i)
      {
        if (b.length(i) != get_n_rows_in_cell(i, 0))
        {
          ErrThrow("Length of Vector Block ", i, " is ", b.length(i),
                   "which does not fit n_rows_in_cell ", get_n_rows_in_cell(i, 0));
        }
        // this method is (indirectly) called from BlockFEMatrix as well, so we
        // should not print a warning here.
        //handle_discovery_of_vector_non_actives(b.length(i)-b.active(i), i);
      }
    }

    ///! Check whether a BlockVector x is fit to be to be factor x in the equation Ax=b.
    void BlockMatrix::check_vector_fits_pre_image(const BlockVector& x) const
    {
      size_t n_vec_blocks = x.n_blocks();

      //check if number of blocks fits number of block columns
      if(n_vec_blocks != n_cell_columns_)
      {
        ErrThrow("Vector blocks number does not fit cell column number. ",
                 n_vec_blocks, "!=", n_cell_columns_);
      }

      //check if each vector block's length fits the number of columns in the corresponding column
      for(size_t j = 0; j<n_vec_blocks; ++j)
      {
        if (x.length(j) != get_n_columns_in_cell(0,j))
        {
          ErrThrow("Length of Vector Block ", j, " is ", x.length(j),
                   "which does not fit n_columns_in_cell ", get_n_columns_in_cell(0,j));
        }
        // this method is (indirectly) called from BlockFEMatrix as well, so we
        // should not print a warning here.
        //handle_discovery_of_vector_non_actives(x.length(j)-x.active(j), j);
      }
    }

    /* ************************************************************************* */
    std::shared_ptr<TMatrix> BlockMatrix::get_combined_matrix() const
    {
      // fill an array with smart pointers to blocks,
      // such thate those stored as transposed really get transposed

      std::vector<std::vector<std::shared_ptr<TMatrix>>> temp_block_grid
      (n_cell_rows_,
       std::vector<std::shared_ptr<TMatrix>>(n_cell_columns_,nullptr));

      // store smart pointers to the already treated transposed colors
      std::vector<std::shared_ptr<TMatrix>> treated_transp_colors{color_count_.size(), nullptr};

      for ( size_t i = 0; i < n_cell_rows_ ; ++i)
      {
        for ( size_t j = 0; j < n_cell_columns_; ++j)
        {
          if (! cell_grid_[i][j].is_transposed_)
          { // non-transposed - just store the cells smart pointer
            temp_block_grid[i][j] = cell_grid_[i][j].block_;
          }
          else
          { //transposed case - check in the treated_transp_colors vector
            size_t color = cell_grid_[i][j].color_;

            if(!treated_transp_colors[color])
            { // this is our sign to make a transposed copy
              treated_transp_colors[color].reset(cell_grid_[i][j].block_->get_transposed());
            }

            temp_block_grid[i][j] = treated_transp_colors[color];

          } // end transp -non-transp

        }// end cell_columns
      }// end cell_rows


      // number of entries of the combined matrix
      size_t n_comb_entries = get_n_total_entries();
      size_t n_comb_rows = get_n_total_rows();
      size_t n_comb_cols = get_n_total_columns();

      // we will create a sparsity structure for the combined matrix. The
      // following two vectors are needed for the constructor
      int * column_of_entry = new int[n_comb_entries];
      int * entries_in_rows = new int[n_comb_rows+1];
      std::vector<double> comb_entries(n_comb_entries, 0.0);
      entries_in_rows[0] = 0;

      // filling the vectors:
      size_t row_offset = 0;
      // position of current entry in combined matrix
      size_t pos = 0;
      // loop over all block rows of this BlockMatrix
      for(size_t block_row = 0; block_row < n_cell_rows_; ++block_row)
      {
        // number of matrix rows in this cell row
        size_t n_rows = cell_grid_[block_row][0].n_rows_;
        // loop over all rows in this cell row
        for(size_t row = 0; row < n_rows; ++row)
        {
          size_t column_offset = 0;
          // loop over all block columns of this (block) row
          for(size_t block_col = 0; block_col < n_cell_columns_; ++block_col)
          {
            // current matrix block
            const TMatrix& cm = *temp_block_grid[block_row][block_col];
            const int * row_ptr = cm.get_row_ptr();
            const int * col_ptr = cm.get_vector_columns();
            const double * entries = cm.GetEntries();
            // loop over entire row in this block
            for(int e = row_ptr[row]; e < row_ptr[row+1]; ++e)
            {
              comb_entries[pos] = entries[e];
              column_of_entry[pos] = col_ptr[e] + column_offset;
              ++pos;
            }
            column_offset += cm.get_n_columns();
          }
          entries_in_rows[row_offset + row + 1] = pos;
        }
        row_offset += n_rows;
      }

      // create sparsity structure
      std::shared_ptr<TStructure> sp(
          new TStructure(n_comb_rows, n_comb_cols, n_comb_entries,
                         column_of_entry, entries_in_rows));

      //Structure copies these object
      delete[] column_of_entry;
      delete[] entries_in_rows;

      // create Matrix
      std::shared_ptr<TMatrix> combined (new TMatrix(sp));
      combined->setEntries(comb_entries);

      return combined;
    }

    std::shared_ptr<TMatrix> BlockMatrix::get_combined_submatrix(
        std::pair<size_t,size_t> upper_left,
        std::pair<size_t,size_t> lower_right) const
    {
      //delegate the work
      BlockMatrix temp = get_sub_blockmatrix(upper_left, lower_right);
      return temp.get_combined_matrix();
    }

    /// @brief total number of columns(> n_block_columns)
    size_t BlockMatrix::get_n_total_columns() const
    { // sum the number of columns in the first block row
      size_t n_total_columns = 0;
      for (size_t j = 0; j < n_cell_columns_ ;++j)
      {
        n_total_columns += cell_grid_[0][j].n_columns_;
      }
      return n_total_columns;
    }

    /// @brief total number of entries
    size_t BlockMatrix::get_n_total_entries() const
    {
      size_t n_total_entries = 0;
      //sum the number of entries from all cells
      for(size_t i= 0; i < n_cell_rows_ ; ++i)
      {
        for (size_t j = 0; j < n_cell_columns_ ;++j)
        {
          n_total_entries += cell_grid_[i][j].block_->get_n_entries();
        }
      }
      return n_total_entries;
    }

    /// @brief total number of rows (> n_block_rows)
    size_t  BlockMatrix::get_n_total_rows() const
    { // sum the number of rows in the first block column
      size_t n_total_rows = 0;
      for (size_t i = 0; i < n_cell_rows_ ;++i)
      {
        n_total_rows += cell_grid_[i][0].n_rows_;
      }
      return n_total_rows;
    }

    size_t BlockMatrix::get_n_square_blocks() const
    {
      size_t n_square_blocks = 0;
      for(size_t i= 0; i < n_cell_rows_ ; ++i)
      {
        for (size_t j = 0; j < n_cell_columns_ ;++j)
        {
          if(cell_grid_[i][j].n_rows_ == cell_grid_[i][j].n_columns_)
            n_square_blocks++;
        }
      }
      return n_square_blocks;
    }


    BlockMatrix BlockMatrix::get_sub_blockmatrix(
        const std::pair<size_t,size_t>& upper_left,
        const std::pair<size_t,size_t>& lower_right) const
    {
      size_t r_first = upper_left.first;
      size_t r_last  = lower_right.first;
      size_t c_first = upper_left.second;
      size_t c_last  = lower_right.second;

      //input check
      if(r_first > r_last)
        ErrThrow("upper_left.first > lower_right.first");
      if(c_first > c_last)
        ErrThrow("upper_left.second > lower_right.second");
      if(r_last >= n_cell_rows_)
        ErrThrow("lower_right.first >= this->n_cell_rows_");
      if(c_last >= n_cell_columns_)
        ErrThrow("lower_right.second >= n_cell_columns_");

      //step 1: construct a blanco sumatrix of correct dimensions
      std::vector<size_t> r_sizes;
      std::vector<size_t> c_sizes;
      for(size_t r = upper_left.first; r<=lower_right.first ; ++r)
      {//fill r_sizes
        r_sizes.push_back(this->cell_grid_[r][0].n_rows_);
      }
      for(size_t c = upper_left.second; c <= lower_right.second; ++c)
      {
        c_sizes.push_back(this->cell_grid_[0][c].n_columns_);
      }
      BlockMatrix sub_matrix(r_sizes, c_sizes);

      //step 2: fill in the correctly colored blocks
      // step 2: fill in the blocks - maintaining the correct coloring
      std::vector< int > known_colors;
      std::vector< std::shared_ptr<TMatrix> > known_mats; //actually FEMatrices...

      for( size_t r = r_first; r < r_last + 1; ++r)
      {
        for( size_t c = c_first; c < c_last + 1; ++c)
        {

          int old_color = this->cell_grid_[r][c].color_;
          auto known = std::find(known_colors.begin(), known_colors.end(), old_color);

          if(known == known_colors.end())
          {//case: this color appears for the first time
            known_colors.push_back(old_color);

            //this is actually a shared ptr to FEMatrix...
            std::shared_ptr<TMatrix> new_mat =
                create_block_shared_pointer(*cell_grid_[r][c].block_.get());

            known_mats.push_back(new_mat);
            //reset known
            known = known_colors.end()-1;
          }

          size_t new_color = std::distance(known_colors.begin(), known); //position in vector
          size_t transp = this->cell_grid_[r][c].is_transposed_;
          size_t new_r = r - r_first;       // force element and color
          size_t new_c = c - c_first;       // into the new sub matrix' cell grid
          size_t old_color_sub =  sub_matrix.cell_grid_[new_r][new_c].color_;
          sub_matrix.cell_grid_[new_r][new_c].color_ = new_color;
          sub_matrix.cell_grid_[new_r][new_c].is_transposed_ = transp;
          sub_matrix.cell_grid_[new_r][new_c].block_ = known_mats.at(new_color);
          //color count
          ++sub_matrix.color_count_.at(new_color);
          --sub_matrix.color_count_.at(old_color_sub);
        }
      }

      //tidy up the color count vector
      sub_matrix.color_count_.erase(
          std::remove( sub_matrix.color_count_.begin(),
                       sub_matrix.color_count_.end(), 0),
                       sub_matrix.color_count_.end()
      );

      return sub_matrix;
    }


    /* ************************************************************************* */
    void BlockMatrix::print_and_check(const std::string& matrix_name) const
    {
      // do both prints
      print_coloring_pattern(matrix_name);
      print_coloring_count(matrix_name);

      //...and perform the consistency check
      check_coloring();
    }

    /* ************************************************************************* */
    void BlockMatrix::print_coloring_count(const std::string& matrix_name) const
    {
      Output::print("----------------------");
      Output::print(" Color count: ", matrix_name);
      Output::print("----------------------");
      size_t color = 0;
      for (auto count : color_count_)
      {
        Output::print(color, "\t : \t", count );
        ++color;
      }
      Output::print("----------------------");
    }

    /* ************************************************************************* */
    void BlockMatrix::print_coloring_pattern(const std::string& matrix_name,
                                                    bool print_adress) const
    {
      Output::print("-------------------------");
      Output::print(" Coloring pattern: ", matrix_name);
      Output::print("-------------------------");
      for(size_t i = 0; i < n_cell_rows_; ++i)
      {
        std::stringstream out_row;
        out_row << "( ";
        for(size_t j = 0; j < n_cell_columns_; ++j)
        {
          out_row << "\t";
          if(cell_grid_[i][j].block_->get_n_entries()==0)
          {//for a zero-map block the number is printed in brackets
            out_row << "(";
          }
          out_row << cell_grid_[i][j].color_;
          if(cell_grid_[i][j].is_transposed_)
          {
            out_row << "^T";
          }
          if(cell_grid_[i][j].block_->get_n_entries()==0)
          {
            out_row << ")";
          }
          if(print_adress)
          {
            out_row << " ";
            out_row << cell_grid_[i][j].block_.get();
          }
        }
        out_row << "\t )";
        Output::print(out_row.str());
      }
      Output::print("-------------------------");
    }

    /* ************************************************************************* */
    void BlockMatrix::replace_blocks(
        const TMatrix& new_block,
        const std::vector<std::vector<size_t>>& cell_positions,
        const std::vector<bool>& transposed_states)
    {
      std::vector<std::tuple<size_t, size_t, bool>> input_tuples(
          check_and_tupelize_vector_input(cell_positions, transposed_states));

      replace_blocks(new_block, input_tuples );
    }

    /* ************************************************************************* */
    void BlockMatrix::scale(double factor)
    {
      std::vector<std::vector<size_t>> positions;
      for ( size_t i = 0; i < n_cell_rows_ ; ++i)
      {
        for ( size_t j = 0; j < n_cell_columns_; ++j)
        {
          positions.push_back({i,j});
        }
      }
      // scale it!
      scale_blocks(factor, positions);
    }

    /* ************************************************************************* */
    void BlockMatrix::scale_blocks(
        double scaling_factor,
        const std::vector<std::vector<size_t>>& cell_positions )
    {
      std::vector<bool> transposed_states (cell_positions.size(), false);
      std::vector<std::tuple<size_t, size_t, bool>> input_tuples(
          check_and_tupelize_vector_input(cell_positions, transposed_states));

      scale_blocks(scaling_factor, input_tuples);
    }

/* ************************************************************************** */
double BlockMatrix::get(unsigned int i, unsigned int j) const
{
  // first find the block in which the indices are in
  size_t block_row = 0;
  size_t block_col = 0;

  unsigned int row = i; // local copy for a better error message
  unsigned int col = j; // local copy for a better error message

  for(; block_row < this->n_cell_columns_ ; ++block_row)
  {
    if(row >= this->cell_grid_[block_row][0].n_columns_)
      row -= this->cell_grid_[block_row][0].n_columns_;
    else
      // found the block
      break;
  }
  if(block_row == this->n_cell_columns_)
    ErrThrow("could not find an entry in row ", i, ". There are ",
             this->get_n_total_rows(), " rows in this BlockMatrix");
  for(; block_col < this->n_cell_columns_ ; ++block_col)
  {
    if(col >= this->cell_grid_[block_row][block_col].n_columns_)
      col -= this->cell_grid_[block_row][block_col].n_columns_;
    else
      // found the block
      break;
  }
  if(block_col == this->n_cell_columns_)
    ErrThrow("could not find an entry in column ", j, ". There are ",
             this->get_n_total_columns(), " columns in this BlockMatrix");

  auto block = this->cell_grid_[block_row][block_col].block_;
  try
  {
    double ret = block->get(row, col);
    return ret;
  }
  catch(...) // entry is not in the sparsity structure
  {
    return 0.;
  }
}

std::vector<double> BlockMatrix::get_diagonal() const
{
  size_t n_diag_entries = std::min(this->get_n_total_rows(),
                                   this->get_n_total_columns());
  std::vector<double> ret(n_diag_entries); // to be returned
  size_t n_diag_blocks = std::min(this->get_n_cell_rows(),
                                  this->get_n_cell_columns());
  size_t offset = 0;
  // loop over all diagonal blocks
  for(size_t d_block = 0; d_block < n_diag_blocks; ++d_block)
  {
    bool transposed; // flag saying if the block is stored as the transposed
    auto diag_block = this->get_block(d_block, d_block, transposed);
    if(!diag_block->is_square())
      ErrThrow("getting the diagonal of a BlockMatrix with non-square block "
               "on the diagonal. I am not sure if this works as expected");
    // note that this is not a great way to do this, because it creates these
    // additional vectors, instead of writing the diagonal entries directly into
    // the vector 'ret'.
    std::vector<double> block_diagonal = diag_block->get_diagonal();
    std::copy(block_diagonal.begin(), block_diagonal.end(), ret.begin()+offset);
    offset += block_diagonal.size();
  }
  return ret;
}

std::vector<double> BlockMatrix::get_row_sums() const
{
  size_t n_cell_rows = get_n_cell_rows();
  size_t n_cell_cols = get_n_cell_columns();

  size_t n_rows = this->get_n_total_rows();
  size_t row_offset = 0;

  std::vector<double> ret(n_rows, 0.0);

  for (size_t block_row = 0; block_row < n_cell_rows; block_row++)
  {
    size_t n_block_rows = 0;

    for (size_t block_col = 0; block_col < n_cell_cols; block_col++)
    {
      bool transposed;
      auto block = get_block(block_row, block_col, transposed);

      auto block_sums = transposed ?
        block->get_col_sums() : block->get_row_sums();

      size_t n = block_sums.size();

      for (size_t i = 0; i < n; i++)
      {
        ret[row_offset + i] += block_sums[i];
      }

      n_block_rows = n;
    }

    row_offset += n_block_rows;
  }

  return ret;
}

std::vector<double> BlockMatrix::get_col_sums() const
{
  size_t n_cell_rows = get_n_cell_rows();
  size_t n_cell_cols = get_n_cell_columns();

  size_t n_cols = this->get_n_total_columns();

  std::vector<double> ret(n_cols, 0.0);

  for (size_t block_row = 0; block_row < n_cell_rows; block_row++)
  {
    size_t col_offset = 0;

    for (size_t block_col = 0; block_col < n_cell_cols; block_col++)
    {
      bool transposed;
      auto block = get_block(block_row, block_col, transposed);

      auto block_sums = transposed ?
        block->get_row_sums() : block->get_col_sums();

      size_t n = block_sums.size();

      for (size_t i = 0; i < n; i++)
      {
        ret[col_offset + i] += block_sums[i];
      }

      col_offset += n;
    }
  }

  return ret;
}


/* ************************************************************************** */
std::shared_ptr<const TMatrix> BlockMatrix::get_block(size_t cell_row,
                                                      size_t cell_col,
                                                      bool& is_transposed) const
{
  //find out the transposed state
  const BlockMatrix::CellInfo& cell = cell_grid_.at(cell_row).at(cell_col);
  is_transposed = cell.is_transposed_;

  //cast const and FEMatrix (range check is done via "at")
  //std::shared_ptr<const TMatrix> shared
  //= std::dynamic_pointer_cast<const TMatrix>(cell.block_);
  return cell.block_;
}

    /* ************************************************************************* */
    // IMPLEMENTATION OF SPECIAL MEMBER FUNCTIONS
    /* ************************************************************************* */
    BlockMatrix::BlockMatrix(const BlockMatrix& other)
    : n_cell_rows_(other.n_cell_rows_), n_cell_columns_(other.n_cell_columns_),
      cell_grid_(other.cell_grid_), color_count_(other.color_count_)

    {
      // each block instance has to be copied once and all the shared pointers of
      // the same color have to be set pointing to the new block
      std::vector<std::shared_ptr<TMatrix>> treated_colors{color_count_.size(), nullptr};
      for ( size_t i = 0; i < n_cell_rows_ ; ++i)
      {
        for ( size_t j = 0; j < n_cell_columns_; ++j)
        {
          size_t color = cell_grid_[i][j].color_;
          if(!treated_colors[color])
          { // this is our sign to make a copy
            treated_colors[color].reset(new TMatrix(*other.cell_grid_[i][j].block_));

            cell_grid_[i][j].block_ = treated_colors[color];
          }
          else
          { // a pointer is stored already
            cell_grid_[i][j].block_ = treated_colors[color];
          }
        }
      }
    }


    /* ************************************************************************* */

    void swap(BlockMatrix& first, BlockMatrix& second)
    {
      std::swap(first.n_cell_columns_, second.n_cell_columns_);
      std::swap(first.n_cell_rows_, second.n_cell_rows_);
      std::swap(first.cell_grid_, second.cell_grid_);
      std::swap(first.color_count_, second.color_count_);
    }
    /* ************************************************************************* */

    BlockMatrix& BlockMatrix::operator=(BlockMatrix other)
    {
      //do a swap with the copy constructed object "other"
      swap(*this, other);

      return *this;
    }

    /* ************************************************************************* */
    // IMPLEMENTATION OF PRIVATE METHODS
    /* ************************************************************************* */

    /* ************************************************************************* */
    void BlockMatrix::add_scaled_matrix_to_blocks(
        const TMatrix& summand, double scaling_factor,
        std::vector<grid_place_and_mode> row_column_transpose_tuples)
    {
      // first of all check the input, modify if reparable or throw if not so.
      check_grid_fit(summand, row_column_transpose_tuples);

      // check if the replacement requires color splits and if so perform them.
      size_t colorToSplit = std::numeric_limits<size_t>::max();
      while (does_modification_require_color_split(colorToSplit, row_column_transpose_tuples))
      {
        mark_for_color_split(colorToSplit, row_column_transpose_tuples);
        split_color(colorToSplit);
      }

      // after the input is order nicely and the colors split,
      // check if all transposed states of the input match
      // the transposed states of the cells
      compare_transposed_mode(row_column_transpose_tuples);

      size_t searched_color = 0;

      //      // list of colors whose first instance is not in the same transposed
      //      // state as the summand - they have to be treated specially
      //      std::list<size_t> treat_special_colors;

      // delegate the additions to the TMatrices
      for (auto it: row_column_transpose_tuples)
      {
        size_t cell_row = std::get<0>(it);
        size_t cell_column = std::get<1>(it);
        CellInfo& current_cell = cell_grid_[cell_row][cell_column];
        size_t cell_color = current_cell.color_;


        if (cell_color >= searched_color)
        { // we found an untreated color
          //          bool addend_transp = current_cell.is_transposed_;
          //          bool summand_transp = std::get<2>(it);

          //          if(summand_transp == addend_transp) // transposed states match ?
          //          {
          current_cell.block_.get()->add_scaled(summand, scaling_factor);
          //          }
          //          else
          //          { //the color has to be kept in mind for later treatment
          //            treat_special_colors.push_back(cell_color);
          //            Output::print("Special treatment color found, transp not fitting"
          //                          "at ", cell_row ,",",cell_column);
          //          }

          searched_color = cell_color + 1;
          continue;
        }

        //        // check if the color is among those which need a special treatment
        //        auto hit = std::find(treat_special_colors.begin(),
        //                             treat_special_colors.end(), cell_color);
        //        if (hit != treat_special_colors.end())
        //        {
        //          // color needs special treatment
        //          bool addend_transp = current_cell.is_transposed_;
        //          bool summand_transp = std::get<2>(it);
        //
        //          if(summand_transp == addend_transp) // transposed states match ?
        //          {
        //            current_cell.block_.get()->add_scaled(summand, scaling_factor);
        //            //color does not need special treatment anymore
        //            treat_special_colors.erase(hit);
        //            Output::print("Special treatment color removed, found fitting transp"
        //                          "at ", cell_row ,",",cell_column);
        //          }
        //        }
      }

      //      //treat those colors which are still in the treat_special_colors list
      //      // - those where the transposed state of all appearances does not fit
      //      for (auto it : treat_special_colors)
      //      {
      //        //if this is nowhere the case, we have to add the transposed somewhere
      //        {
      //          size_t cell_row;
      //          size_t cell_column;
      //          find_first_appearance_of_color( it, cell_row, cell_column );
      //
      //          CellInfo& current_cell = cell_grid_[cell_row][cell_column];
      //
      //          //make a transposed version of the summand - this is expensive
      //          TMatrix* summand_transposed = summand.get_transposed();
      //
      //          Output::print("Warning! Must transpose a TMatrix to perform addition. \n"
      //              "Consider different way to use BlockMatrix::add_scaled_matrix_to_blocks");
      //
      //          current_cell.block_.get()->add_scaled(*summand_transposed, scaling_factor);
      //
      //          delete summand_transposed;
      //        }
      //      }
    }

    /* ************************************************************************* */
    void BlockMatrix::check_and_edit_input(
        std::vector<grid_place_and_mode>& row_column_transpose_tuples)
    {
      // check if the input is
      //  a) sorted w.r.t. two first indices

      // a lambda function for alphanumerical sort
      auto sort_alphnum = [] (grid_place_and_mode a, grid_place_and_mode b) -> bool {
        if (std::get<0>(a) < std::get<0>(b))
        {
          return true;
        }
        else if (std::get<0>(a) > std::get<0>(b))
        {
          return false;
        }
        else
        {
          return (std::get<1>(a) < std::get<1>(b));
        }
      };
      std::sort(row_column_transpose_tuples.begin(),
                row_column_transpose_tuples.end(),
                sort_alphnum);

      //  b) unique w.r.t. two first indices - throw if there is an ambiguity of the
      //     type "do for transposed and non-transposed"

      // lambda function for that purpose
      auto unique_wrt_position = [] (grid_place_and_mode a, grid_place_and_mode b) -> bool {
        // catch the case of one block being modified in transp and non-tranps mode
        if( std::get<0>(a) == std::get<0>(b)
            && std::get<1>(a) == std::get<1>(b)
            && std::get<2>(a) != std::get<2>(b))
        {
          ErrThrow("Block can not be"
              " modified in transposed AND non-transposed state!"
              " [", std::get<0>(a) ," , ", std::get<1>(a) ,"]" );
        }
        return (std::get<0>(a) == std::get<0>(b)
            && std::get<1>(a) == std::get<1>(b));
      };

      // move duplicates to the end of the vector and hold an iterator
      // to the first duplicate-position
      auto unique_end = std::unique(row_column_transpose_tuples.begin(),
                                    row_column_transpose_tuples.end(),
                                    unique_wrt_position);

      // give a warning if there was duplicate input
      if(unique_end != row_column_transpose_tuples.end())
      {
        Output::print("Warning! Duplicates in the input vector removed.");
      }

      //crop the vector to range of uniques
      row_column_transpose_tuples.erase(
          unique_end, row_column_transpose_tuples.end());

    }

    /* ************************************************************************* */
    std::vector<std::tuple<size_t, size_t, bool>> BlockMatrix::
    check_and_tupelize_vector_input(
        const std::vector<std::vector<size_t>>& cell_positions,
        const std::vector<bool>& transposed_states )
    {
      if (cell_positions.size() != transposed_states.size())
      {
        ErrThrow("Different input vector lenghts!");
      }

      std::vector<grid_place_and_mode> tuples;

      for (size_t i = 0; i < cell_positions.size(); ++i)
      {
        if (cell_positions[i].size() != 2)
        {
          ErrThrow("Cell_positions at ", i, " is not an index pair!");
        }
        size_t cell_row = cell_positions[i][0];
        size_t cell_column = cell_positions[i][1];
        bool transp = transposed_states[i];

        tuples.push_back(std::make_tuple(cell_row, cell_column, transp));
      }

      return tuples;

    }

    /* ************************************************************************* */
    void BlockMatrix::check_color_count() const
    {
      //a list of the cell colors from left to right, top to bottom
      std::list<size_t> foundColors;
      for(size_t i = 0 ; i < n_cell_rows_ ; ++i)
      {
        for(size_t j = 0; j < n_cell_columns_; ++j)
        {
          foundColors.push_back(cell_grid_[i][j].color_);
        }
      }

      // check for every color if it is assigned and counted as often as
      // color_count_ states
      for(size_t color = 0; color < color_count_.size() ; ++color )
      {

        // count how often color appears in foundColors
        size_t n_finds = std::count(foundColors.begin(), foundColors.end(), color);

        if(n_finds == 0) //throw if something's wrong
        {
          ErrThrow("Here is an unassigned color: ", color);
        }
        if( n_finds != color_count_[color] )
        {
          ErrThrow("Number of found colors and stored color_count_"
              "do not match for color: ", color);
        }

      }
    }

    /* ************************************************************************* */
    void BlockMatrix::check_coloring_order() const
    {
      //we rely on condition 1 here
      size_t searched_for = 0;

      for(size_t i = 0 ; i < n_cell_rows_ ; ++i)
      {
        for(size_t j = 0; j < n_cell_columns_; ++j)
        {
          if (cell_grid_[i][j].color_ > searched_for)
          {
            ErrThrow("The color ordering of this BlockMatrix is incorrect.")
          }
          else if (cell_grid_[i][j].color_ == searched_for)
          {
            ++searched_for;
          }
        }
      }
    }



    /* ************************************************************************* */
    void BlockMatrix::check_equivalence_of_relations() const
    {
      //traverse the matrix to get the first element for the comparison
      for(size_t i_first = 0 ; i_first < n_cell_rows_ ; ++i_first)
      {
        for(size_t j_first = 0; j_first < n_cell_columns_; ++j_first)
        {
          if( is_last_index_pair(i_first, j_first) )
          {
            //last index pair reached
            return;
          }

          //traverse the matrix to get the second element for the comparison

          // start in the same row as first element, unless that lies in the last column:
          // then go one row further down
          size_t i_second{i_first};
          size_t j_second{j_first};

          get_next_cell_grid_index(i_second, j_second);

          for( ; i_second < n_cell_rows_ ; ++i_second)
          {
            // start one column right from first element, unless we are already in
            // a further donw row. then start at zero
            for( ; j_second < n_cell_columns_; ++j_second)
            {
              const CellInfo first = cell_grid_[i_first][j_first];
              const CellInfo& second = cell_grid_[i_second][j_second];
              if (! does_color_match_block(first, second) )
              {
                ErrThrow(" The equaivalence of relations is broken in this BlockMatrix.");
              }
            }
          }
        }
      }
    }

    /* ************************************************************************* */

    void BlockMatrix::check_grid_fit(
        const TMatrix& matrix,
        std::vector<grid_place_and_mode>& row_column_transpose_tuples
    ) const
    {
      // sort input and remove duplicates if need be
      check_and_edit_input(row_column_transpose_tuples);

      // check for index out of bounds
      check_indices(row_column_transpose_tuples);

      //  check if the matrix fits into the given cells
      for (auto it: row_column_transpose_tuples)
      {
        size_t cell_row = std::get<0>(it);
        size_t cell_column = std::get<1>(it);

        const CellInfo& currentCell = cell_grid_[cell_row][cell_column];
        bool transposed = std::get<2>(it);

        if(!does_block_fit_cell(matrix, currentCell, transposed))
        {
          ErrThrow("The given matrix will not fit into cell "
              "[", std::get<0>(it) ," , ", std::get<1>(it) ,"]");
        }
      }
    }

    /* ************************************************************************* */
    void BlockMatrix::check_indices(
        const std::vector<grid_place_and_mode>& row_column_transpose_tuples
    ) const
    {
      for (auto it: row_column_transpose_tuples)
      {
        size_t cell_row = std::get<0>(it);
        size_t cell_column = std::get<1>(it);
        if (cell_row >= n_cell_rows_ || cell_column >= n_cell_columns_)
        {
          ErrThrow("Cell index pair out of block matrix bounds! "
              "[", std::get<0>(it) ," , ", std::get<1>(it) ,"]");
        }
      }
    }

    /* ************************************************************************* */
    void BlockMatrix::check_re_coloring_flags() const
    {
      for(size_t i = 0 ; i < n_cell_rows_ ; ++i)
      {
        for(size_t j = 0; j < n_cell_columns_; ++j)
        {
          if (cell_grid_[i][j].re_color_flag_ != CellInfo::ReColoringFlag::KEEP)
          {
            ErrThrow("Oops! Somebody forgot to clean up a"
                "recoloring flag in cell [",i,",",j,"].");
          }
        }
      }
    }

    /* ************************************************************************* */
    void BlockMatrix::compare_transposed_mode(
        std::vector<grid_place_and_mode>& row_column_transpose_tuples) const
    {
      for ( auto it : row_column_transpose_tuples)
      {
        size_t i = std::get<0>(it);
        size_t j = std::get<1>(it);
        bool you_are_transposed = std::get<2>(it);
        bool i_am_transposed=cell_grid_[i][j].is_transposed_;
        if(you_are_transposed != i_am_transposed)
        {
          ErrThrow("Transposed state of input does not fit transposed state of cell."
              "Rearrange input!");
        }
      }
    }

    /* ************************************************************************* */
    std::shared_ptr<TMatrix> BlockMatrix::create_block_shared_pointer(const TMatrix& block) const
    {

      //Output::print("Called base class copy and store");
      return std::make_shared<TMatrix>(block);
    }

    /* ************************************************************************* */
    bool BlockMatrix::does_block_fit_cell(
        const TMatrix& block,
        const CellInfo& cell,
        bool transposed
    )
    {
      size_t n_rows_cell = cell.n_rows_;
      size_t n_columns_cell = cell.n_columns_;

      size_t n_rows_block = transposed ? block.get_n_columns() : block.get_n_rows();
      size_t n_columns_block = transposed ? block.get_n_rows() : block.get_n_columns();

      // if both dimension fit return true, otherwise false
      return (n_rows_cell == n_rows_block && n_columns_cell ==  n_columns_block);
    }

    /* ************************************************************************* */
    bool BlockMatrix::does_color_match_block(const BlockMatrix::CellInfo& first,
                                                    const BlockMatrix::CellInfo& second)
    {
      return (first.color_ == second.color_) == (first.block_ == second.block_);
    }

    /* ************************************************************************* */
    bool BlockMatrix::does_modification_require_color_split(
        size_t& color_to_split,
        std::vector<grid_place_and_mode> row_column_transposed_tuples) const
    {
      // the part from here performs in O(m log m), when m is size of row_column_tuples

      // fill a list with the colors which will be affected
      // by the modification and sort it
      std::list<size_t> color_touches;
      for(auto it : row_column_transposed_tuples)
      {
        size_t color = cell_grid_[std::get<0>(it)][std::get<1>(it)].color_;
        color_touches.push_back( color );
      }
      color_touches.sort();

      //push back one more element as "stop" element
      color_touches.push_back(std::numeric_limits<size_t>::max());

      // count out how many times each affected color gets affected
      // and compare to the number of its overall appearances
      size_t color = color_touches.front();
      size_t n_touches = 0;

      for (auto it : color_touches)
      {
        if(it != color) // finished counting of a certain color
        {
          if (n_touches - color_count_[color] != 0)
          { // we found a color class which is not entirely affected by
            // the modification - this class must undergo the splitting procedure

            // update the output and return true
            color_to_split = color;
            return true;
          }
          else
          {
            //reset number of touches and color
            color = it;
            n_touches = 0;
          }

        }
        ++n_touches;
      }


      // great, no color split needed - we finished in O(m log m)
      return false;
    }

    /* ************************************************************************* */
    bool BlockMatrix::does_replace_require_color_merge(
        size_t& color_a, size_t& color_b,
        std::vector<grid_place_and_mode> row_column_transposed_tuples) const
    {
      // grab the first color as "fixed_color"
      size_t first_cell_row = std::get<0>(*row_column_transposed_tuples.begin());
      size_t first_cell_column = std::get<1>(*row_column_transposed_tuples.begin());
      size_t fixed_color = cell_grid_[first_cell_row][first_cell_column].color_;

      // find out if any other color is affected by the replacement
      for(auto it : row_column_transposed_tuples)
      {
        size_t cell_row = std::get<0>(it);
        size_t cell_column = std::get<1>(it);
        size_t cell_color = cell_grid_[cell_row][cell_column].color_;

        if (fixed_color != cell_color)
        { // found two different participating colors!
          color_a = fixed_color;
          color_b = cell_color;
          return true;
        }

      }
      // all participating cells belong to the same color class already
      color_a = fixed_color;
      color_b = fixed_color;
      return false;

    }


    /* ************************************************************************* */
    void BlockMatrix::find_first_appearance_of_color(
        size_t color_to_find, size_t& block_row , size_t& block_column) const
    {
      if( color_to_find >= get_n_colors() )
      {
        ErrThrow("That color does not exist in this BlockMatrix.");
      }

      for(size_t i = 0; i < n_cell_rows_ ; ++i)
      {
        for(size_t j = 0; j < n_cell_columns_ ; ++j)
        {
          if (cell_grid_[i][j].color_ == color_to_find)
          {
            // color was found
            block_row = i;
            block_column = j;
            return;
          }
        }
      }
    }

    /* ************************************************************************* */
    void BlockMatrix::find_first_appearance_of_color_and_mode(
        size_t color_to_find, bool find_transposed,
        size_t& block_row , size_t& block_column ) const
    {
      if( color_to_find >= get_n_colors() )
      {
        ErrThrow("That color does not exist in this BlockMatrix.");
      }

      for(size_t i = 0; i < n_cell_rows_ ; ++i)
      {
        for(size_t j = 0; j < n_cell_columns_ ; ++j)
        {
          if (cell_grid_[i][j].color_ == color_to_find
              && cell_grid_[i][j].is_transposed_ == find_transposed)
          {
            // color was found
            block_row = i;
            block_column = j;
            return;
          }
        }
      }
      throw std::runtime_error("Could not find a cell of that color and mode.");
    }

    /* ************************************************************************* */
    void BlockMatrix::get_next_cell_grid_index (
        size_t& block_row, size_t& block_column ) const
    {

      if( block_row > n_cell_rows_
          || block_column > n_cell_columns_)
      {
        throw std::runtime_error("Index is out of bounds!");
      }

      if( is_last_index_pair ( block_row, block_column ) )
      {
        throw std::logic_error("Next index pair after last required!");
      }

      size_t new_row = (block_column < n_cell_columns_ - 1) ? block_row : block_row + 1;
      size_t new_column = (block_row == new_row) ? block_column + 1 : 0;

      block_row = new_row;
      block_column = new_column;

    }

    /* ************************************************************************* */
    void BlockMatrix::handle_discovery_of_vector_non_actives(
      const int n_nonActive, const int) const
    {
      if(n_nonActive != 0)
      {
        //maybe put to virtual method: handle_discovery_of_vector_non_actives
        // give a warning if the vector has actives - the matrix has definitely not!
        Output::print<2>("Warning! The BlockVector has actives, but BlockMatrix does not."
                         " Did you want to use a BlockFEMatrix instead?");
      }
    }
    /* ************************************************************************* */
    bool BlockMatrix::is_last_index_pair(size_t block_row, size_t block_column) const
    {
      if( block_row > n_cell_rows_
          || block_column > n_cell_columns_)
      {
        throw std::runtime_error("Index is out of bounds!");
      }

      if( block_row == n_cell_rows_ -1
          && block_column == n_cell_columns_ - 1)
      {
        return true;
      }
      return false;
    }

    /* ************************************************************************* */
    void BlockMatrix::mark_for_color_split(
      size_t color_to_mark,
      const std::vector<grid_place_and_mode>& row_column_transposed_tuples)
    {
      for (auto it : row_column_transposed_tuples)
      {
        if (cell_grid_[std::get<0>(it)][std::get<1>(it)].color_ == color_to_mark)
        {
          //the cell is of the color to split and belongs to the given subset
          // - so mark it for split
          cell_grid_[std::get<0>(it)][std::get<1>(it)].re_color_flag_ =
              CellInfo::ReColoringFlag::SPLIT;
        }
      }
    }

    /* ************************************************************************* */
    void BlockMatrix::merge_colors(size_t color_a, size_t color_b)
    {
      //hold the color of the merged class
      size_t new_color = std::min(color_a, color_b);

      // hold the color which will get replaced
      size_t replaced_color = std::max(color_a, color_b);

      // a shared pointer to the matrix which the merged color class will share
      std::shared_ptr<TMatrix> shared_matrix;
      bool searching_shared_matrix = true;

      //traverse the entire cell grid
      for (size_t i = 0; i < n_cell_rows_; ++i)
      {
        for (size_t j = 0; j < n_cell_columns_; ++j)
        {
          CellInfo& currentCell = cell_grid_[i][j];


          if (currentCell.color_ == color_a || currentCell.color_ == color_b)
          {
            // found  matrix to be merged into new class
            if(searching_shared_matrix)
            {
              //found the first appeareance of the lower class, hold its block
              shared_matrix = currentCell.block_;
              searching_shared_matrix = false;
            } // end if searching shared matrix

            // replace matrix and color and update color count (it's even done when unnecessary)
            color_count_[currentCell.color_] -= 1;
            currentCell.color_ = new_color;
            currentCell.block_ = shared_matrix;
            color_count_[currentCell.color_] += 1;

          } //end: found cell of new merged class

          else if (currentCell.color_ > replaced_color)
          {
            // we found a class whose color has to be reduced by one
            color_count_[currentCell.color_] -= 1;
            currentCell.color_ -= 1;
            color_count_[currentCell.color_] += 1;
          }
        }
      } //end traversing entire matrix

      // remove the last color class - should be empty (0)!
      if(color_count_.back() != 0)
      {
        ErrThrow("Something's terribly wrong here - last color class is not empty!")
      }

      color_count_.pop_back();
    }

    /* ************************************************************************* */
    void BlockMatrix::replace_blocks(
        const TMatrix& new_block,
        std::vector<grid_place_and_mode> row_column_transpose_tuples)
    {
      // first of all check the input, modify if reparable or throw if not so.
      check_grid_fit(new_block, row_column_transpose_tuples);

      // check if the replacement requires color splits and if so perform them.
      size_t colorToSplit = std::numeric_limits<size_t>::max();
      while (does_modification_require_color_split(colorToSplit, row_column_transpose_tuples))
      {
        mark_for_color_split(colorToSplit, row_column_transpose_tuples);
        split_color(colorToSplit);
      }

      // check if the replacement requires color merges and if so perform them.
      size_t colorA = std::numeric_limits<size_t>::max();
      size_t colorB = std::numeric_limits<size_t>::max();
      while(does_replace_require_color_merge(colorA, colorB, row_column_transpose_tuples))
      {
        merge_colors(colorA, colorB);
      }

      // now everything works out - do the actual replacement!

      //wrap shared ptr around (correctly typed) copy of new_block
      std::shared_ptr<TMatrix> new_block_shared = create_block_shared_pointer(new_block);

      for (auto it : row_column_transpose_tuples)
      {
        size_t cell_row = std::get<0>(it);
        size_t cell_column = std::get<1>(it);
        bool transposed = std::get<2>(it);
        //copy assign new shared pointer
        cell_grid_[cell_row][cell_column].block_ = new_block_shared;
        // set the transposed state correctly
        cell_grid_[cell_row][cell_column].is_transposed_ = transposed;
      }
    }

    /* ************************************************************************* */
    void BlockMatrix::scale_blocks(
        double scaling_factor,
        std::vector<grid_place_and_mode> row_column_transpose_tuples)
    {
      // first of all check the input, modify if reparable or throw if not so.
      check_and_edit_input( row_column_transpose_tuples );

      //check for index-out-of-bound
      check_indices( row_column_transpose_tuples );

      // check if the replacement requires color splits and if so perform them.
      size_t colorToSplit = std::numeric_limits<size_t>::max();
      while (does_modification_require_color_split(colorToSplit, row_column_transpose_tuples))
      {
        mark_for_color_split(colorToSplit, row_column_transpose_tuples);
        split_color(colorToSplit);
      }

      size_t searched_color = 0;

      // delegate the scaling to the TMatrices
      for (auto it: row_column_transpose_tuples)
      {
        size_t cell_row = std::get<0>(it);
        size_t cell_column = std::get<1>(it);
        CellInfo& current_cell = cell_grid_[cell_row][cell_column];
        size_t cell_color = current_cell.color_;

        if (cell_color >= searched_color)
        {
          // we found an untreated color

          current_cell.block_.get()->scale( scaling_factor );
          searched_color = cell_color + 1;
        }
      }
    }

    /* ************************************************************************* */
    void BlockMatrix::split_color(size_t color_to_split)
    {
      Output::root_warn("BlockMatrix", "A color of this BlockMatrix had to be split."
          " Is that what you intended?");

      // a new color is going to appear - start with 0 count
      color_count_.push_back(0);

      size_t i = 0;
      size_t j = 0;

      find_first_appearance_of_color(color_to_split, i , j);

      //deep copy the matrix and make the shared_ptr temporarily responsible
      std::shared_ptr<TMatrix> matrix_copy = create_block_shared_pointer(*cell_grid_[i][j].block_.get());

      // Whether the first found cell of the color to split
      // has flag SPLIT or flag KEEP determines, which part of the color
      // class gets the new numbering - the one marked with value
      // equal to first_split maintains its current color number
      CellInfo::ReColoringFlag first_split(cell_grid_[i][j].re_color_flag_);

      //reset recoloring flag in this first found cell
      cell_grid_[i][j].re_color_flag_ = CellInfo::ReColoringFlag::KEEP;

      if(is_last_index_pair(i,j))
      {
        ErrThrow("It should not happen that the first found element of a split "
            "color is the last cell grid entry. Something's very wrong!");
      }

      size_t new_color = color_to_split + 1;
      bool searching_first_cell_of_new_class{true};

      // traverse through the remainder of the matrix (after the first appearance)
      while (! is_last_index_pair(i,j) )
      {
        get_next_cell_grid_index(i,j); //iterate up

        CellInfo& current_cell = cell_grid_[i][j];


        if(current_cell.color_ >= new_color) // a cell of a higher color
        {//found a cell of higher color
          if (searching_first_cell_of_new_class)
          {//this means we have not yet found the first element
            // of the second part of the color class to split
            // and thus still have to determine its number

            // update new_color, because we have found another class
            // which lies between our two split parts
            ++new_color;
          }
          else
          { //we are not searching anymore, the new_color is determined
            //so recolor and update the number of colored matrices array
            color_count_[current_cell.color_] -= 1;
            ++current_cell.color_;
            color_count_[current_cell.color_] += 1;

          }
        } //end: cell of higher color class
        else if (current_cell.color_ == color_to_split)
        {//found another element of the split color

          if (current_cell.re_color_flag_ == first_split)
          {//the element belongs to the first set and thus keeps its color
            //just reset the re color flags
            current_cell.re_color_flag_ = CellInfo::ReColoringFlag::KEEP;
          }
          else
          {
            //end the search, by now we found the new color!
            // (now this is done every time, but nevermind)
            searching_first_cell_of_new_class = false;

            //The element will get assigned the new color
            color_count_[current_cell.color_] -= 1;
            current_cell.color_ = new_color;
            color_count_[current_cell.color_] += 1;

            current_cell.block_ = matrix_copy; //the new shared_ptr

            //..and resets the re color flag
            current_cell.re_color_flag_ = CellInfo::ReColoringFlag::KEEP;
          }
        } //end: cell of split color class
      } //endwhile

    }

/* ************************************************************************* */
void BlockMatrix::reset()
{
  this->scale(0.0);
}

int BlockMatrix::clean_denormals()
{
  int c = 0;

  for (size_t i = 0; i < n_cell_rows_; ++i)
  {
    for (size_t j = 0; j < n_cell_columns_; ++j)
    {
      int cij = cell_grid_[i][j].block_->clean_denormals();

      if (cij > 0)
      {
        Output::root_warn("BlockMatrix", "Block ", i, " x ", j, " had ",
          cij, " non-normal components");
      }

      c += cij;
    }
  }

  // MPI: no allreduce here because TMatrix::clean_denormals already does that

  return c;
}