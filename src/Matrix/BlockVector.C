#include <BlockVector.h>
#include <BlockMatrix.h>
#include <BlockFEMatrix.h>
#include <LinAlg.h>
#include "MooNMD_Io.h"
#include <cmath>

#ifdef _MPI
#include <ParFECommunicator3D.h>
#include <mpi.h>
#endif

#include <algorithm>

/** ************************************************************************ */
BlockVector::BlockVector() : entries(), lengths(), actives(), inners()
{
  Output::print<5>("Constructor of BlockVector with no arguments");
}

/** ************************************************************************ */
BlockVector::BlockVector(const std::vector<unsigned int>& lengths)
{
  unsigned int n_entries = std::accumulate(lengths.begin(), lengths.end(), 0);
  entries = std::vector<double>(n_entries , 0.0);

  this->lengths = lengths;

  actives = lengths; //assume only active entries
  inners = lengths;
}

/** ************************************************************************ */
BlockVector::BlockVector(unsigned int length)
 : entries(length, 0.0), lengths(1, length),
 actives(1, length), inners(1, length)
{
  Output::print<5>("Constructor of BlockVector with length ", length);
}

/** ************************************************************************ */
BlockVector::BlockVector(int length)
 : BlockVector(length > 0 ? (unsigned int)length : 0)
{
  if (length < 0)
  {
    ErrThrow("cannot construct BlockVector with negative number of entries");
  }
}

BlockVector::BlockVector(const BlockMatrix& mat, bool result)
{
  Output::print<5>("Constructor of BlockVector using a BlockMatrix");

  // the total length of this vector
  size_t total_length = result ? mat.get_n_total_rows() : mat.get_n_total_columns();

  // number of blocks in this BlockVector
  size_t n_blocks = result ? mat.get_n_cell_rows() : mat.get_n_cell_columns();

  entries.resize(total_length, 0.0);

  // set all entries to zero
  std::fill(entries.begin(), entries.end(), 0.0);
  lengths.resize(n_blocks, 0);
  actives.resize(n_blocks, 0);
  inners.resize(n_blocks, 0);

  // determine the length of each vector block
  for (size_t b = 0; b < n_blocks; b++)
  {
    if (result)
    {
      //vector on rhs
      lengths[b] = mat.get_n_rows_in_cell(b,0);
    }
    else
    {
      //vector as factor
      lengths[b] = mat.get_n_columns_in_cell(0,b);
    }

    // all entries are active
    actives[b] = lengths[b];
    inners[b] = lengths[b];
  }

  Output::print<5>("(end) Constructor of BlockVector using a BlockMatrix");
  Output::print<5>("entries: ", this->entries.size(), " ",this->length());
}

BlockVector::BlockVector(const BlockFEMatrix& mat, bool result)
: BlockVector(static_cast<const BlockMatrix&>(mat), result)
{
  // now set actives correctly

  for (size_t b = 0; b < n_blocks(); ++b)
  {
    if (result)
    {
      // take actives from row space
      actives[b] = mat.get_n_row_actives(b);
      inners[b] = mat.get_n_row_inner(b);
    }
    else
    {
      // take actives from column space
      actives[b] = mat.get_n_column_actives(b);
      inners[b] = mat.get_n_column_inner(b);
    }
  }
}

/** ************************************************************************ */
unsigned int BlockVector::offset(unsigned int b) const
{
  if (b >= this->n_blocks())
  {
    ErrThrow("trying to access block ", b, ", but there are only ",
             this->n_blocks(), " blocks in this BlockVector");
  }

  return std::accumulate(lengths.cbegin(), lengths.cbegin() + b, 0);
}

/** ************************************************************************ */
void BlockVector::reset()
{
  std::fill(this->entries.begin(), this->entries.end(), 0.0);
}

/** ************************************************************************ */
void BlockVector::ResetActive()
{
  unsigned int n_blocks = this->n_blocks();

  auto it = this->entries.begin();

  for(unsigned int i = 0; i < n_blocks; i++)
  {
    std::fill(it, it + this->active(i), 0.0);
    std::advance(it, this->length(i)); // it += this->length(i);
  }
}

/** ************************************************************************ */
void BlockVector::ResetNonActive()
{
  unsigned int n_blocks = this->n_blocks();

  auto it = this->entries.begin();

  for(unsigned int i = 0; i < n_blocks; i++)
  {
    std::fill(it + this->active(i), it + this->length(i), 0.0);
    std::advance(it, this->length(i)); // it += this->length(i);
  }
}

/** ************************************************************************ */
void BlockVector::ResetInner()
{
  unsigned int n_blocks = this->n_blocks();

  auto it = this->entries.begin();

  for(unsigned int i = 0; i < n_blocks; i++)
  {
    std::fill(it, it + this->inner(i), 0.0);
    std::advance(it, this->length(i)); // it += this->length(i);
  }
}

/** ************************************************************************ */
void BlockVector::ResetBoundary()
{
  unsigned int n_blocks = this->n_blocks();

  auto it = this->entries.begin();

  for(unsigned int i = 0; i < n_blocks; i++)
  {
    std::fill(it + this->inner(i), it + this->length(i), 0.0);
    std::advance(it, this->length(i)); // it += this->length(i);
  }
}

/** ************************************************************************ */
void BlockVector::scale(const double a, const unsigned int i)
{
  if ((unsigned int)i < lengths.size())
  {
    Dscal(lengths.at(i), a, this->block(i)); // scale i-th subvector
  }
  else
  {
    ErrThrow("trying to scale subvector ", i, " which does not exist");
  }
}

/** ************************************************************************ */
void BlockVector::scaleActive(const double a)
{
  if (a == 0)
  {
    this->reset(); /// @todo only reset active or all dofs?
  }
  else if (a != 1.0)
  {
    for (unsigned int i = 0; i < this->n_blocks(); ++i)
    {
      Dscal(actives.at(i), a, this->block(i));
    }
  }
}
/** ************************************************************************* */
void BlockVector::scaleNonActive(const double a)
{
  if (a == 0)
  {
    this->ResetNonActive();
  }
  else if (a != 1.0)
  {
    for (unsigned int i = 0; i < this->n_blocks(); ++i)
    {
      // do the dscal for the nonactives only (includes pointer arithmetic...)
      int n_non_actives = this->lengths.at(i)- this->actives.at(i);
      Dscal(n_non_actives, a, this->block(i) + this->actives.at(i));
    }
  }
}

/** ************************************************************************ */
void BlockVector::scale(const double a)
{
  if (a == 0.0)
  {
    this->reset();
  }
  else if (a != 1.0)
  {
    Dscal(this->length(), a, this->get_entries());
  }
}


/** ************************************************************************ */
void BlockVector::add_scaled(const BlockVector& r, double factor)
{
  const unsigned int l = this->length(); // length of this BlockVector
  if(r.length() != l)
  {
    ErrThrow("unable to add two BlockVectors of different lengths\t",
             l, "\t", r.length());
  }
  else if(this->n_blocks() != r.n_blocks())
  {
    Output::print("WARNING: BlockVector::operator+=\n adding to BlockVectors "
                  "with the same length but different numbers of blocks");
  }
  Daxpy(l, factor, r.get_entries(), this->get_entries());
}

/** ************************************************************************ */
void BlockVector::addScaledActive(const BlockVector& r, double factor)
{
  if(this->n_blocks() != r.n_blocks())
    ErrThrow("number of blocks in two vectors must be same");

  for(unsigned int i=0; i<this->n_blocks(); ++i)
  {
    Daxpy(this->actives.at(i), factor, r.block(i), this->block(i));
  }
}

/** ************************************************************************ */
void BlockVector::addScaledNonActive(const BlockVector& r, double factor)
{
  if(this->n_blocks() != r.n_blocks())
    ErrThrow("number of blocks in two vectors must be same");

  for(unsigned int i=0; i<this->n_blocks(); ++i)
  {
    if(this->actives.at(i) != r.active(i))
    {
      ErrThrow("Number of actives in block ",i," do not match!");
    }
    // do the daxpy for the nonactives only (includes pointer arithmetic...)
    int n_non_actives = this->lengths.at(i)- this->actives.at(i);
    Daxpy(n_non_actives , factor, r.block(i) + this->actives.at(i), this->block(i) + this->actives.at(i));
  }
}

/** ************************************************************************ */
void BlockVector::copy(const double * x, const int i)
{
  if(i < 0)
    *this = x; // copy entire vector
  else if((unsigned int)i < this->n_blocks())
    std::copy(x, x + this->length(i), this->block(i));
    //memcpy(this->block(i), x, this->length(i)*sizeof(double)); // in string.h
  else
  {
    ErrThrow("trying to copy to subvector ", i, " which does not exist.");
  }
}

/** ************************************************************************ */
void BlockVector::copy_nonactive(const BlockVector& r)
{
  const unsigned int l = this->length(); // length of this BlockVector
  if(r.length() != l)
  {
    ErrThrow("unable to copy from one BlockVector to another one of different ",
             "length\t", l, "\t", r.length());
  }
  else if(this->n_blocks() != r.n_blocks())
  {
    Output::print("WARNING: BlockVector::operator+=\n adding to BlockVectors "
                  "with the same length but different numbers of blocks\n");
  }

  for(unsigned int b = 0, n_b = this->n_blocks(); b < n_b; ++b)
  {
    if(this->lengths[b] > this->active(b))
      std::copy(r.block(b) + r.active(b), r.block(b) + r.length(b),
                this->block(b) + this->active(b));
      // in string.h
      //memcpy(this->block(b) + this->active(b), r.block(b) + this->active(b),
      //       (this->lengths[b] - this->active(b))*sizeof(double));
  }
}

/** ************************************************************************ */
void BlockVector::add(const double* x, const int i, double a)
{
  if(i < 0)
  {
    const unsigned int l = length(); // length of this BlockVector
    Daxpy(l, a, x, this->get_entries());
  }
  else if((unsigned int)i < this->lengths.size())
  {
    double * b = this->block(i);
    for(unsigned int j = 0; j < this->length(i); j++)
    {
      b[j] += a*x[j];
    }
  }
  else
  {
    ErrThrow("trying to add to subvector ", i, " which does not exist.");
  }
}


/** ************************************************************************ */
double BlockVector::norm(const std::vector<unsigned int>& blocks
#ifdef _MPI
   , std::vector<const TParFECommunicator3D*> comms
#endif
) const
{
#ifdef _MPI
  if (comms.size() == 0)
  {
    //If comms is not provided, a default comm with size 0 is given, and
    // the norm is calculated as it is in the seq case on each process,
    // but it won't give the correct results. Once, all iterative solvers
    // are correctly parallelised, the following warning can be converted
    // to ErrThrow.
    Output::warn<1>("BlockVector", "You are calculating the norm of a vector "
        "in the parallel case, without providing the Communicators. This will "
        "definitely give you a WRONG norm value in terms of consistency.");
#endif

  if (blocks.empty())
  {
    return Dnorm(this->length(), this->get_entries());
  }
  else
  {
    auto n_considered_blocks = blocks.size();
    double norm = 0;

    for (auto i = 0u; i < n_considered_blocks; ++i)
    {
      if (blocks[i] >= this->n_blocks())
      {
        ErrThrow("unable to compute norm for block ", blocks[i],
                 ". There are only ", this->n_blocks(),
                 " blocks in this BlockVector");
      }

      norm += Ddot(this->length(blocks[i]),
                   this->block(blocks[i]),
                   this->block(blocks[i]));
    }

    return std::sqrt(norm);
  }

#ifdef _MPI
  }
  else
  {
    auto blocks_copy = blocks; // work on a copy
    if (blocks_copy.empty())
    {
      blocks_copy.resize(this->n_blocks(), 0);
      std::iota (std::begin(blocks_copy), std::end(blocks_copy), 0); // {0, 1, 2, ...}
    }

    // This MPI method makes only use of values of master dofs, therefore
    // "this" does not have to be updated, consistency level 0 is enough.

    auto n_considered_blocks = blocks_copy.size();

    /// First check if vector and communicators fit
    if (comms.size() != n_considered_blocks)
    {
      ErrThrow("Number of blocks does not equal number of communicators.",
               n_considered_blocks, " ", comms.size());
    }

    for (size_t i = 0; i < n_considered_blocks; ++i)
    {
      if (comms[i]->GetNDof() != (int)length(blocks_copy[i]))
      {
        ErrThrow("Length of Block ", blocks_copy[i], " and comms ", i,
                 " do not match. ", comms[i]->GetNDof(), " ",
                 length(blocks_copy[i]));
      }
    }

    int my_rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    double sum_local = 0;
    for (size_t i = 0; i < n_considered_blocks; ++i)
    {
      const int* masters = comms[i]->GetMaster();

      size_t offset = this->offset(blocks_copy[i]);
      size_t length_of_block = length(blocks_copy[i]);

      for (size_t j = 0; j < length_of_block; ++j)
      {
        if (masters[j] == my_rank)
        {
          double val = entries[offset + j];
          sum_local += val * val;
        }
      }
    }

    // Now add up all local sums via MPI_Allreduce.
    MPI_Allreduce(MPI_IN_PLACE, &sum_local, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return std::sqrt(sum_local); //square root of the global sum
  }
#endif
}

double BlockVector::norm_infty(const std::vector<unsigned int>& blocks
#ifdef _MPI
   , std::vector<const TParFECommunicator3D*> comms
#endif
) const
{
#ifdef _MPI
  if (comms.size() == 0)
  {
    ErrThrow("No communicators provided.");
#endif

  if (blocks.empty())
  {
    auto n = this->length();
    const double* entries = this->get_entries();

    double norm = 0.0;

    for (auto j = 0u; j < n; j++)
    {
      norm = std::max(norm, std::abs(entries[j]));
    }

    return norm;
  }
  else
  {
    auto n_considered_blocks = blocks.size();
    double norm = 0;

    for (auto i = 0u; i < n_considered_blocks; ++i)
    {
      if (blocks[i] >= this->n_blocks())
      {
        ErrThrow("unable to compute norm for block ", blocks[i],
                 ". There are only ", this->n_blocks(),
                 " blocks in this BlockVector");
      }

      const double* b = this->block(blocks[i]);
      unsigned int n = this->length(blocks[i]);

      for (auto j = 0u; j < n; j++)
      {
        norm = std::max(norm, std::abs(b[j]));
      }
    }

    return norm;
  }

#ifdef _MPI
  }
  else
  {
    auto blocks_copy = blocks; // work on a copy
    if (blocks_copy.empty())
    {
      blocks_copy.resize(this->n_blocks(), 0);
      std::iota (std::begin(blocks_copy), std::end(blocks_copy), 0); // {0, 1, 2, ...}
    }

    // This MPI method makes only use of values of master dofs, therefore
    // "this" does not have to be updated, consistency level 0 is enough.

    auto n_considered_blocks = blocks_copy.size();

    /// First check if vector and communicators fit
    if (comms.size() != n_considered_blocks)
    {
      ErrThrow("Number of blocks does not equal number of communicators.",
               n_considered_blocks, " ", comms.size());
    }

    for (size_t i = 0; i < n_considered_blocks; ++i)
    {
      if (comms[i]->GetNDof() != (int)length(blocks_copy[i]))
      {
        ErrThrow("Length of Block ", blocks_copy[i], " and comms ", i,
                 " do not match. ", comms[i]->GetNDof(), " ",
                 length(blocks_copy[i]));
      }
    }

    int my_rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    double max_local = 0.0;
    for (size_t i = 0; i < n_considered_blocks; ++i)
    {
      const int* masters = comms[i]->GetMaster();

      size_t offset = this->offset(blocks_copy[i]);
      size_t length_of_block = length(blocks_copy[i]);

      for (size_t j = 0; j < length_of_block; ++j)
      {
        if (masters[j] == my_rank)
        {
          max_local = std::max(max_local, std::abs(entries[offset + j]));
        }
      }
    }

    // Now take the global maximum
    MPI_Allreduce(MPI_IN_PLACE, &max_local, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    return max_local;
  }
#endif
}

/** ************************************************************************ */
void BlockVector::print(const std::string& name, const int iB) const
{
  if(iB < 0)
  { // print full BlockVector
    for (unsigned int i = 0, l = this->length(); i < l; i++)
      Output::print(name , "(" , i , ")= " , this->at(i) , ";");
  }
  else if((unsigned int)iB < lengths.size())
  {
    const double *myBlock = this->block(iB);
    for (unsigned int i = 0; i < this->lengths[iB]; i++)
      Output::print(name ,"(" , i , ")= " , myBlock[i] ,";");
  }
  else
    ErrThrow("trying to print subvector ", iB, " which does not exist");
}

/** ************************************************************************ */
void BlockVector::write(const std::string& filename) const
{
  std::ofstream vectorfile;
  vectorfile.open(filename.c_str());

  //write the header line - array format, real values
  vectorfile << "%%MatrixMarket matrix array real general \n";

  //write general matrix information
  vectorfile << length() << "\t" << 1 << "\t" << "\n";

  //loop through matrix and print each entry
  for(auto it : entries)
  {
    vectorfile << setprecision(16) <<  it << "\n";
  }


  vectorfile.close();
}

/** ************************************************************************ */
void BlockVector::info() const
{
  using namespace Output;
  Output::print(" | info to BlockVector :");
  unsigned int nb = n_blocks();
  Output::print(" | -- total length ", length(), "\tnumber of blocks ", nb);
  if(nb > 1)
  {
    for(unsigned int i = 0; i < nb; i++)
    {
      if(this->lengths[i] != this->actives[i])
      {
        Output::print(" | ---- block ", i, " has length ", this->lengths[i],
                      " (", this->actives[i], " active)");
      }
      else
      {
        Output::print(" | ---- block ", i, " has length ", this->lengths[i]);
      }
    }
  }
}

/** ************************************************************************ */
BlockVector& BlockVector::operator=(const double *r)
{
  if(this->get_entries() == r)
    // both are the same, no copying necessary
    return *this;
  std::copy(r, r + this->length(), this->entries.begin());
  //memcpy(this->get_entries(), r, this->length()*sizeof(double));// in string.h
  return *this;
}

/** ************************************************************************ */
BlockVector& BlockVector::operator=(const double a)
{
  std::fill(this->entries.begin(), this->entries.end(), a);
  return *this;
}

/** ************************************************************************ */
BlockVector& BlockVector::operator*=(const double a)
{
  this->scale(a);
  return *this;
}

/** ************************************************************************ */
BlockVector& BlockVector::operator+=(const BlockVector& r)
{
  this->add_scaled(r, 1.0);
  return *this;
}

/** ************************************************************************ */
BlockVector& BlockVector::operator-=(const BlockVector& r)
{
  this->add_scaled(r, -1.0);
  return *this;
}

/** ************************************************************************ */
double dot(const BlockVector& a, const BlockVector& b
#ifdef _MPI
           , std::vector<const TParFECommunicator3D*> comms
#endif
)
{
  unsigned int l = a.length();
  // Check that the 2 vectors are compatible (same length and block number)
  if(b.length() != l)
  {
    ErrThrow("BlockVector class: Unable to compute dot product of two "
        "BlockVectors of different lengths\t", l, "\t", b.length());
  }
  else if(a.n_blocks() != b.n_blocks())
  {
    Output::warn<1>("BlockVector", "computing the dot product of two "
                    "vectors with the same length but "
                    "different numbers of blocks");
  }
  // Compute the Dot-product in SEQ and MPI cases
#ifndef _MPI
  return Ddot(l, a.get_entries(), b.get_entries());
#elif _MPI
  if (comms.size()==0)
  {
    // If comms is not provided, a default comm with size 0 is given, and
    // the dotproduct is calculated as it is in the seq case on each
    // process,but it won't give the correct results. Once, all
    // iterative solvers are correctly parallelised, the following
    // warning can be converted to ErrThrow.
    Output::warn<1>("BlockVector", "You are calculating the dot-product of 2 vectors "
        "in the parallel case, without providing the Communicators. This will "
        "definitely give you a WRONG value in terms of consistency.");
    return Ddot(l, a.get_entries(), b.get_entries());
  }
  else
    return dot_global(a, b, comms);
#endif
}

/** ************************************************************************ */
#ifdef _MPI
double dot_global(const BlockVector& a, const BlockVector& b,
                  std::vector<const TParFECommunicator3D*> comms)
{
  // This MPI method makes only use of values of master dofs, therefore
  // a and b do not have to be updated, consistency level 0 is enough.

  /* First check if vector and communicators fit.
   * before calling this method, it was already checked that the length
   * and block numbers of the 2 vectors match. So it is enough to make
   * the following check only for one of the 2 vectors, for example a.
   */
  if(comms.size() != a.n_blocks() || comms.size() != b.n_blocks() )
  {
    ErrThrow("Number of blocks does not equal number of communicators.",
             a.n_blocks(), " ", b.n_blocks(), " ", comms.size() );
  }
  for(size_t n =0; n < a.n_blocks(); ++n)
  {
    if(comms[n]->GetNDof() != (int) a.length(n)
        || comms[n]->GetNDof() != (int) b.length(n) )
    {
      ErrThrow("Length of Block number", n, " and comms", n, " do not match. ",
               comms[n]->GetNDof(), " ", a.length(n), " ", b.length(n));
    }
  }

  int my_rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  double sum_local = 0;
  for(size_t n =0; n < a.n_blocks(); ++n) // n is the block number
  {
    const TParFECommunicator3D* comm = comms[n]; //for convenience
    const int* masters = comm->GetMaster();

    size_t offset_a = a.offset(n);
    size_t offset_b = b.offset(n);
    if (offset_a != offset_b)
      ErrThrow("BlockVector : 2 Vectors have different block offsets!");

    for(size_t i = 0; i < a.length(n); ++i )
    {
      if(masters[i] == my_rank)
      {
        double val_a = a.entries[offset_a + i];
        double val_b = b.entries[offset_b + i];
        sum_local += val_a * val_b;
      }
    }
  }

  // Now add up all local sums via MPI_Allreduce.
  double sendbf[1] = {sum_local};
  double recvbf[1];
  MPI_Allreduce(sendbf,recvbf,1,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);

  double product_result = recvbf[0]; // global sum (no square root)
  return product_result;
}

class DotInfo
{
  public:

    enum class OperationType { dot, norm_2, norm_infty };

    static void Flush()
    {
      MPI_Waitall(pendingRequests.size(), pendingRequests.data(), MPI_STATUSES_IGNORE);

      pendingRequests.clear();

      for (DotInfo* info: pendingInfo)
      {
        info->Finalize();
      }

      pendingInfo.clear();
    }

    static void Queue(const double source, double& target, OperationType op = OperationType::dot)
    {
      pendingInfo.push_back(new DotInfo(source, target, op));
    }
  private:
    DotInfo(const double source, double& target, OperationType op)
      : op(op), source(source), target(target)
    {
      Send();
    }

    void Send()
    {
      MPI_Request req;

      switch (op)
      {
        case OperationType::norm_infty:
          MPI_Iallreduce(&source, &target, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, &req);
          break;

        default:
          MPI_Iallreduce(&source, &target, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &req);
          break;
      }

      pendingRequests.push_back(req);
    }

    void Finalize()
    {
      switch (op)
      {
        case OperationType::norm_2:
          target = std::sqrt(target);
          break;

        default:
          break;
      }

      delete this;
    }

    int size;
    const OperationType op;
    const double source;
    double &target;

    static std::vector<DotInfo*> pendingInfo;
    static std::vector<MPI_Request> pendingRequests;
};

std::vector<DotInfo*> DotInfo::pendingInfo;
std::vector<MPI_Request> DotInfo::pendingRequests;

void queue_dot_global(const BlockVector& a, const BlockVector& b,
                  std::vector<const TParFECommunicator3D*> comms,
                  double& recv)
{
  // This MPI method makes only use of values of master dofs, therefore
  // a and b do not have to be updated, consistency level 0 is enough.

  /* First check if vector and communicators fit.
   * before calling this method, it was already checked that the length
   * and block numbers of the 2 vectors match. So it is enough to make
   * the following check only for one of the 2 vectors, for example a.
   */
  if (comms.size() != a.n_blocks() || comms.size() != b.n_blocks())
  {
    ErrThrow("Number of blocks does not equal number of communicators.",
             a.n_blocks(), " ", b.n_blocks(), " ", comms.size() );
  }

  for (size_t n = 0; n < a.n_blocks(); ++n)
  {
    if (comms[n]->GetNDof() != (int) a.length(n)
      || comms[n]->GetNDof() != (int) b.length(n))
    {
      ErrThrow("Length of Block number", n, " and comms", n, " do not match. ",
               comms[n]->GetNDof(), " ", a.length(n), " ", b.length(n));
    }
  }

  int my_rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  double sum_local = 0;
  for (size_t n = 0; n < a.n_blocks(); ++n) // n is the block number
  {
    const TParFECommunicator3D* comm = comms[n]; // for convenience
    const int* masters = comm->GetMaster();

    size_t offset_a = a.offset(n);
    size_t offset_b = b.offset(n);

    if (offset_a != offset_b)
    {
      ErrThrow("BlockVector : 2 Vectors have different block offsets!");
    }

    size_t len = a.length(n);
    for (size_t i = 0; i < len; ++i )
    {
      if (masters[i] == my_rank)
      {
        int o = offset_a + i;
        sum_local += a.entries[o] * b.entries[o];
      }
    }
  }

  DotInfo::Queue(sum_local, recv, DotInfo::OperationType::dot);
}

void queue_norm_global(const BlockVector& v,
                             std::vector<const TParFECommunicator3D*> comms,
                             double& recv)
{
  if (comms.size() != v.n_blocks())
  {
    ErrThrow("Number of blocks does not equal number of communicators.",
             v.n_blocks(), " ", comms.size());
  }

  for (size_t n = 0; n < v.n_blocks(); ++n)
  {
    if (comms[n]->GetNDof() != (int) v.length(n))
    {
      ErrThrow("Length of Block number", n, " and comms", n, " do not match. ",
               comms[n]->GetNDof(), " ", v.length(n));
    }
  }

  int my_rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  double sum_local = 0;
  for (size_t n = 0; n < v.n_blocks(); ++n) // n is the block number
  {
    const TParFECommunicator3D* comm = comms[n]; // for convenience
    const int* masters = comm->GetMaster();

    size_t offset = v.offset(n);

    size_t len = v.length(n);
    for (size_t i = 0; i < len; ++i )
    {
      if (masters[i] == my_rank)
      {
        double val = v.entries[offset + i];
        sum_local += val * val;
      }
    }
  }

  DotInfo::Queue(sum_local, recv, DotInfo::OperationType::norm_2);
}

void queue_norm_infty_global(const BlockVector& v,
                             std::vector<const TParFECommunicator3D*> comms,
                             double& recv)
{
  if (comms.size() != v.n_blocks())
  {
    ErrThrow("Number of blocks does not equal number of communicators.",
             v.n_blocks(), " ", comms.size());
  }

  for (size_t n = 0; n < v.n_blocks(); ++n)
  {
    if (comms[n]->GetNDof() != (int) v.length(n))
    {
      ErrThrow("Length of Block number", n, " and comms", n, " do not match. ",
               comms[n]->GetNDof(), " ", v.length(n));
    }
  }

  int my_rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  double max_local = 0;
  for (size_t n = 0; n < v.n_blocks(); ++n) // n is the block number
  {
    const TParFECommunicator3D* comm = comms[n]; // for convenience
    const int* masters = comm->GetMaster();

    size_t offset = v.offset(n);

    size_t len = v.length(n);
    for (size_t i = 0; i < len; ++i )
    {
      if (masters[i] == my_rank)
      {
        max_local = std::max(max_local, std::abs(v.entries[offset + i]));
      }
    }
  }

  DotInfo::Queue(max_local, recv, DotInfo::OperationType::norm_infty);
}

void flush_dot_queue()
{
  DotInfo::Flush();
}

void flush_norm_queue()
{
  DotInfo::Flush();
}

#endif

int BlockVector::clean_denormals()
{
  int c = 0;
  unsigned int n = length();

  for (unsigned int i = 0; i < n; i++)
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

/** ************************************************************************ */
void BlockVector::copy_structure(const BlockVector& r)
{
  entries.resize(r.length(), 0.0);

  lengths.resize(r.n_blocks());
  actives.resize(r.n_blocks());
  inners.resize(r.n_blocks());

  for (unsigned int i = 0; i < r.n_blocks(); i++)
  {
    lengths[i] = r.length(i);
    actives[i] = r.active(i);
    inners[i] = r.inner(i);
  }
}

/** ************************************************************************ */
void BlockVector::write_to_stream(std::ostream& os) const
{
  if(!os.good())
  {
    ErrThrow("cannot read stream to save data. return");
  }
  os << this->length() << endl;
  os.write((char *) this->get_entries(), sizeof(double) * this->length());
}

/** ************************************************************************ */
void BlockVector::read_from_stream(std::istream& is)
{
  if(!is.good())
  {
    ErrThrow("stream in bad state, cannot read data");
  }

  // check if the size of this vector and the number of doubles in the file
  // coincide:
  unsigned int length_read;
  std::string line;
  std::getline(is, line);
  std::istringstream parser( std::string( line.begin(), line.end() ) );
  parser >> length_read >> std::ws;

  if(parser.good())
  {
    ErrThrow("formatting error, the first line in the file should contain ",
             "only a number indicating the number of entries to be read");
  }
  else if(this->length() == 0)
  {
    // the vector is not yet filled with any data, allocate the arrays:
    this->entries.resize(length_read, 0.0);
    lengths.resize(1, length_read);
    actives.resize(1, length_read);
    inners.resize(1, length_read);
  }
  else if(this->length() != length_read)
  {
    ErrThrow("unexpected size to be read. Expected: ", this->length(),
             "\tin file: ", length_read);
  }
  // now we are sure that (this->length == length_read)

  is.read((char *)this->get_entries(), sizeof(double) * length_read);
}


/** ************************************************************************ */
void BlockVector::write_to_file(const std::string& filename) const
{
  std::ofstream dat(filename);
  if(!dat)
  {
    Output::warn("BlockVector::write_to_file", "cannot open file '", filename,
                 "' to save data. return");
    return;
  }
  this->write_to_stream(dat);
  dat.close();
}

/** ************************************************************************ */
void BlockVector::read_from_file(const std::string& filename)
{
  std::ifstream dat(filename);
  if(!dat)
  {
    ErrThrow("cannot open file '", filename, "' to read data");
  }
  this->read_from_stream(dat);
  dat.close();
}

/** ************************************************************************ */
double* BlockVector::block(const unsigned int i)
{
  return this->get_entries() + this->offset(i);
}

/** ************************************************************************ */
const double* BlockVector::block(const unsigned int i) const
{
  return this->get_entries() + this->offset(i);
}

/** ************************************************************************ */
double& BlockVector::at(const unsigned int i)
{
  try
  {
    return entries.at(i);
  }
  catch(...)
    ErrThrow("index out of bounds");
}

/** ************************************************************************ */
const double& BlockVector::at(const unsigned int i) const
{
  try
  {
    return entries.at(i);
  }
  catch(...)
  {
    ErrThrow("index out of bounds");
  }
}
