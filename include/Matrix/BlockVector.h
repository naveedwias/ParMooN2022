/** ************************************************************************
*
* @class       BlockVector
* @brief      Store a block-vector
*
* Store block-vector which consists of several subvectors
* (diffenrent length is possible). Some algebraic operations on the
* subvectors are implemented
*
* @author     Ulrich Wilbrandt
* @date       01.11.2013
*
* @ruleof0
*
****************************************************************************/

#ifndef __BLOCKVECTOR__
#define __BLOCKVECTOR__

class BlockMatrix; // forward declaration
class BlockFEMatrix;

#include <numeric>
#include <vector>
#include <string>

#ifdef _MPI
class TParFECommunicator3D;
#endif

class BlockVector;

// This is a forward declaration with default argument for the friend
// method dot. We need to split this function this way to make it work
// with both compilers Clang and Gnu. If we just define friend dot in the
// class with a default argument, Clang doesn't compile.
double dot(const BlockVector& a, const BlockVector& b
#ifdef _MPI
    , std::vector<const TParFECommunicator3D*> comms={}
#endif
);

class BlockVector
{
  protected:
    /** @brief vector of (pointers to) blocks. blocks[0] is the entire vector */
    std::vector<double> entries;

    /** @brief the lengths of each block
     *
     * The size of this vector also determines the number of blocks.
     */
    std::vector<unsigned int> lengths;

    /** @brief the number of active degrees of freedom for each block
     *
     * The size of this vector is the same as that of 'lengths'.
     */
    std::vector<unsigned int> actives;

    /** @brief the number of inner degrees of freedom for each block
     *
     * The size of this vector is the same as that of 'lengths'.
     */
    std::vector<unsigned int> inners;

    /** @brief return the accumulated length of all blocks i with i < b
     *
     * This tells you the start index of a given block in the entries vector.
     */
    unsigned int offset(unsigned int b) const;

  public:

    /** standard constructor */
    BlockVector();

    /** Construct a BlockVector of length length.size(), where block i has
     * length[i] entries, filled with zeroes.*/
    explicit BlockVector(const std::vector<unsigned int>& lengths);

    /** constructor for a BlockVector consisting of a single block of length
     * 'l' filled with zeros.
     */
    explicit BlockVector(unsigned int l);

    /** constructor for a BlockVector consisting of a single block of length
     * 'l' filled with zeros.
     *
     * @note This constructor is for backwards compatibility, to avoid
     * explicit casts of integers to unsigned int. Use carefully.
     */
    explicit BlockVector(int l);


    /// Construct a BlockVector which is suitable to serve as factor ("false")
    /// or result ("true") in multiplication with a BlockMatrix.
    /// Vector is filled with zeroes and all entries are non-active.
    BlockVector(const BlockMatrix& mat, bool result = false);

    BlockVector(const BlockFEMatrix& mat, bool result = false);

    //Declaration of special member functions - rule of zero

    //! Default copy constructor. Performs deep copy.
    BlockVector(const BlockVector&) = default;

    //! Default move constructor.
    BlockVector(BlockVector&&) = default;

    //! Default copy assignment operator. Performs deep copy.
    BlockVector& operator=(const BlockVector&) = default;

    //! Default move assignment operator
    BlockVector& operator=(BlockVector&&) = default;

    //! Default destructor.
    ~BlockVector() = default;

    /**
     * @brief Set all entries to zero
     *
     * Reset all entries of all subvectors to zero.
     */
    void reset();

    /**
     * @brief Set all active entries to zero
     *
     * Reset all active entries of all subvectors to zero leaving only Dirichlet
     * values unchanged.
     *
     */
    void ResetActive();

    /**
     * @brief Set all non-active entries to zero
     *
     * Reset all non-active entries of all subvectors to zero leaving all
     * non-Dirichlet values unchanged.
     */
    void ResetNonActive();

    /**
     * @brief Set all inner entries to zero
     *
     * Reset all inner entries of all subvectors to zero leaving only boundary
     * values unchanged.
     *
     */
    void ResetInner();

    /**
     * @brief Set all boundary entries to zero
     *
     * Reset all boundary (not just Dirichlet) entries of all subvectors to zero
     * leaving all inner values unchanged.
     */
    void ResetBoundary();

    /**
     * @brief Scale a subvector
     *
     * The subvector with index i is scaled by a.
     *
     * @param a factor by which the subvector is scaled
     * @param i index of target subvector
     *
     */
    void scale(const double a, const unsigned int i);

    /**
     * @brief scale a subvector
     * This subvector with index i is scaled by a, only the actives.
     */
    void scaleActive(const double a);

    /**
     * @brief scale a subvector
     * This subvector with index i is scaled by a, only the non actives.
     */
    void scaleNonActive( const double a);

    /**
     * @brief Scale the entire vector
     */
    void scale(const double a);

    /**
     * @brief add scaled vector to this
     *
     * @param r BlockVector which is added to this
     * @param factor factor with which r is multiplied
     *
     */
    void add_scaled(const BlockVector& r, double factor);

    /**
     * @brief add scaled vector to this, only actives
     *
     * @param r BlockVector which is added to this
     * @param factor factor with which r is multiplied
     *
     */
    void addScaledActive(const BlockVector& r, double factor);

    /**
     * @brief add scaled vector to this, only non-actives
     *
     * @param r BlockVector which is added to this
     * @param factor factor with which r is multiplied
     *
     */
    void addScaledNonActive(const BlockVector& r, double factor);

    /** @brief copy the structure of another BlockVector,
     *
     * No values are copied. If there are old values, they are deleted! New
     * values are set to 0.
     *
     * @param r BlockVector from which the structure is copied to this one
     */
    void copy_structure(const BlockVector & r);

    /**
     * @brief add (scaled) values into a subvector, "this += a*x"
     *
     * The values given by 'x' are added to the subvector with index 'i'. If
     * 'i' is negative, the entire BlockVector is added. The scalar 'a' is a
     * scaling factor.
     *
     * @param[in] x values which are copied to the target subvector
     * @param i index of target subvector, -1 means entire vector
     * @param a factor by which x is multiplied
     *
     */
    void add(const double * x, const int i = -1, double a = 1.0);

    /**
     * @brief copy values into a subvector
     *
     * The values given by x are added to the subvector with index 'i'.
     * If 'i' is negative, the entire BlockVector is added (as for operator=)
     *
     * @param x values which are copied to the target subvector
     * @param i index of target subvector
     *
     */
    void copy(const double * x, const int i = -1);

    /** @brief copy the nonactive values of each block of r to this */
    void copy_nonactive(const BlockVector& r);

    /**
     * @brief compute norm of this BlockVector
     *
     * compute the 2-norm (square root of sum of squares)
     * Note: possibly implement other norms as well
     * @param[in] blocks compute norm only for the specified subblocks, defaults
     * to all blocks.
     */
    double norm(const std::vector<unsigned int>& blocks
#ifdef _MPI
       , std::vector<const TParFECommunicator3D*> comms={}
#endif
    ) const;

    double norm_infty(const std::vector<unsigned int>& blocks
#ifdef _MPI
       , std::vector<const TParFECommunicator3D*> comms={}
#endif
    ) const;

#ifdef _MPI
    double norm(std::vector<const TParFECommunicator3D*> comms) const
    { return norm(std::vector<unsigned int>{}, comms); }

    void queue_norm(double& recv, std::vector<const TParFECommunicator3D*> comms) const
    { queue_norm_infty_global(*this, comms, recv); }

    double norm_infty(std::vector<const TParFECommunicator3D*> comms) const
    { return norm_infty(std::vector<unsigned int>{}, comms); }

    void queue_norm_infty(double& recv, std::vector<const TParFECommunicator3D*> comms) const
    { queue_norm_global(*this, comms, recv); }
#endif

    double norm() const
    { return norm(std::vector<unsigned int>{}); }

    double norm_infty() const
    { return norm_infty(std::vector<unsigned int>{}); }

    int clean_denormals();

    /**
     * @brief Print subvector iB to console in Matlab format
     *
     * If iB < 0, the entire block vector is printed to console. Otherwise the
     * subvector of index iB is printed.
     * Format: name(i)= val       (i>0)
     *
     * @param name name for output
     * @param iB index of subvector
     *
     */
    void print(const std::string& name = "rhs", const int iB = -1) const;

    /**
     * Write entire Vector into an Outfile in MatrixMarket array format.
     * @param filename Desired filename (and path).
     */
    void write(const std::string& filename) const;

    /**
     * @brief Print some information without explicitly printing values
     */
    void info() const;

    /**
     * @brief write the entire BlockVector to a stream
     *
     * The stream will start with a line containing only one number, indicating
     * the length of the vector. The following lines are just written and are
     * then not human readable.
     *
     * A stream created with this method can be read into a BlockVector using
     * BlockVector::read_from_stream(std::istream);
     */
    void write_to_stream(std::ostream& os) const;

    /**
     * @brief read data into this BlockVector from a stream
     *
     * The stream has to start with a line containing only one number, indicating
     * the length of the vector. The following lines are just read in and are
     * not human readable.
     *
     * A file created with BlockVector::write_fo_strean(std::ostream) can be
     * read using this method.
     */
    void read_from_stream(std::istream& is);

    /**
     * @brief write the entire BlockVector to a file
     *
     * The file will start with a line containing only one number, indicating
     * the length of the vector. The following lines are just written and are
     * then not human readable.
     *
     * A file created with this method can be read into a BlockVector using
     * BlockVector::read_from_file(std::string);
     */
    void write_to_file(const std::string& filename) const;

    /**
     * @brief read data into this BlockVector from a file
     *
     * The file has to start with a line containing only one number, indicating
     * the length of the vector. The following lines are just read in and are
     * not human readable.
     *
     * A file created with BlockVector::write_fo_file(std::string) can be read
     * using this method.
     */
    void read_from_file(const std::string& filename);

    /** getters */

    unsigned int length() const
    { return std::accumulate(lengths.begin(),lengths.end(),0); }

    unsigned int length(const int i) const
    { return lengths.at(i); }

    /** @brief return the number of active entries for a given block i */
    unsigned int active(const int i) const
    { return actives.at(i); }

    /** @brief return the number of inner entries for a given block i */
    unsigned int inner(const int i) const
    { return inners.at(i); }

    /** @brief return the number of non-active entries for a given block i */
    size_t n_non_actives(const int i) const
    { return lengths.at(i) - actives.at(i); }

    unsigned int n_blocks() const
    { return lengths.size(); }

    const double* block(const unsigned int i) const;

    double* block(const unsigned int i);

    const std::vector<double>* get_vector_entries() const
    { return &(this->entries); }

    const double* get_entries() const
    { return &(this->entries.at(0)); }

    double* get_entries()
    { return &(this->entries.at(0)); }

    /**
     * @return A const reference to the entire entries vector (all blocks).
     */
    const std::vector<double>& get_entries_vector() const
    {
      return entries;
    }

    /**
     * @return A non-const reference to the entire entries vector (all blocks).
     */
    std::vector<double>& get_entries_vector()
    {
      return entries;
    }

    /* *********************************************************************
     * overloading standard operators
     * we avoid operators +,-,*,/ because they would require to alocate new
     * memory, which is often undesirable and unnecessary.
     */

    /** array acces:
     * For a BlockVector v and unsigned int i use either 'v.at(i)', 'v[i]' or
     * 'v(i)' to acces the i-th element. All three methods provide the same
     * functionality. The function 'at()' also checks if the given index is not
     * too large. Furthermore for each method two versions are implemented. One
     * returns an lvalue so that a value can be set at the given index
     * (e.g. a[i]=5;). The other leaves this BlockVector constant
     * (e.g. in cout << a[i]);
     */
    double & at(const unsigned int i);
    const double & at(const unsigned int i) const;
    double& operator[](const unsigned int i)
    { return entries[i]; }
    const double& operator[](const unsigned int i) const
    { return entries[i]; }
    double& operator()(const unsigned int i)
    { return entries[i]; }
    const double& operator()(const unsigned int i) const
    { return entries[i]; }


    // copy, r should be as long as entire BlockVector
    BlockVector& operator=(const double *r);
    BlockVector& operator=(const double a); // set all values to a
    BlockVector& operator*=(const double a); // multiply by a scalar
    BlockVector& operator+=(const BlockVector& r); // add other BlockVector
    BlockVector& operator-=(const BlockVector& r);// substract other BlockVector

    /** ******************************************************************** */

    /** friends */

    /**
     * @brief compute dot product of two BlockVectors
     *
     * it is checked if the two BlockVectors are of equal length
     *
     * @note we put 'friend' here, so that we don't need a static function call
     *       in e.g. a template iterative solver
     * @param a,b the two BlockVectors
     */
    friend double dot(const BlockVector& a, const BlockVector& b
#ifdef _MPI
        , std::vector<const TParFECommunicator3D*> comms
#endif
    );

#ifdef _MPI
    /// Compute global dot-product of 2 vectors distributed among the processes. All
    /// processes will return the same results.
    friend double dot_global(const BlockVector& a, const BlockVector& b,
                             std::vector<const TParFECommunicator3D*> comms);

    friend void queue_dot_global(const BlockVector& a, const BlockVector& b,
                             std::vector<const TParFECommunicator3D*> comms,
                             double& recv);

    friend void queue_norm_global(const BlockVector& v,
                             std::vector<const TParFECommunicator3D*> comms,
                             double& recv);

    friend void queue_norm_infty_global(const BlockVector& v,
                             std::vector<const TParFECommunicator3D*> comms,
                             double& recv);
#endif

};

#ifdef _MPI
void flush_dot_queue();
void flush_norm_queue();
#endif

#endif