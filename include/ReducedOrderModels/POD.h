/** ****************************************************************************
*
* @name       POD
* @brief      Computation/manipulation of POD basis
*             independently of problem dimension
*
* @author     Swetlana Giere & Alfonso Caiazzo
* @date       08.03.2017 (start of implementation), 15.1.2019 (restart)
*
*******************************************************************************/

#ifndef POD_H
#define POD_H

#include <BlockFEMatrix.h>
#include <DenseMatrix.h>
#include <MooNMD_Io.h>
#include <Matrix.h>
#include <ParameterDatabase.h>
#include <fstream>

class POD
{
  public:
    /**
    * @brief constructor
    *
    */
    POD(const ParameterDatabase& param_db, std::string filename_suffix="");

    /**
    * @brief default destructor
    *
    */
    ~POD();

    /**
     * Creates a database filled with default parameters. This database will
     * contain all necessary parameters for the POD.
     */
    static ParameterDatabase default_pod_database();

   /** @brief check parameters in database
    *
    * This function checks if the parameters of the input file are consistent
    * with the parameters stored in the database
    * If some parameters are inconsistent, an error occurs and throws
    * an exception.
    */
    void check_parameters();

    int get_rank() const
    { return rank; }

    int get_length() const
    { return length; }

    const DenseMatrix* get_basis() const
    { return &pod_basis; }

    /// @brief get the j-th vector of the POD basis
//     std::vector<double> get_basis(size_t j) const
//     { return pod_basis.get_matrix_column(j); }
//     const double* get_basis_vect(size_t j) const
//     { return const_cast<double*>(const_cast<DenseMatrix>(pod_basis).getEntry_ptr(0, j)); }
    double* get_basis_vect(size_t j)
    { return pod_basis.getEntry_ptr(0, j); }

//     const double* get_snaps_avr_ptr() const
//     { return *snaps_mean; }
    double* get_snaps_avr_ptr()
    { return &(this->snaps_mean.at(0)); }

    const std::vector<double> get_snaps_avr() const
    { return snaps_mean; }

    const double* get_pod_eigs(void) const
    { return eigs; }

    std::string get_basename(void) const
    { return basename; }


    /**
    * @brief set pointer to gramian matrix (matrix describing inner product for
    *        POD). Make sure that the TMatrix mat is still the gramian matrix
    *        when needed
    *
    * Set the gramian matrix describing the inner product with respect to which
    * POD basis will be computed. In the MPI case, the matrix must be stored such
    * that the global inner product can be computed by computing the rank-wise
    * inner product and then summing over the ranks. (For e.g. the L^2 inner
    * product, simply computing the L^2 product for each local pair of DOFs
    * regardless of ownership/tag will accomplish this.)
    *
    * NOTE: If no gramian matrix is set, then it is automatically set to
    * identity, i.e. POD basis will be computed with respect to the euclidean
    * inner product.
    *
    * @param[out] mat pointer to ParMooN matrix which should be the gramian
    *                 matrix
    */
    void set_gramian_ptr(std::shared_ptr<TMatrix> mat);

    void set_projection(const BlockFEMatrix& mat, const ParameterDatabase& solver_db);

    /**
    * @brief compute POD basis
    *
    * Compute POD basis from snapshots stored in member variable snaps_mat.
    * If this->db["pod_fluctuations_only"]==true then POD basis will be computed
    * from the fluctuating part of the snapshots.
    * The rank of the basis is given by the class member variable this->rank
    * (or db["pod_rank"]). If this->rank<=0, then all possible POD basis
    * functions will be computed and stored.
    */
    void compute_basis(int l);


  protected:
    /* parameter database */
    ParameterDatabase db;

    /* basename suffix (useful whenn different POD basis are computed in a model
    * for instance in TNSE: velocity and pressure)
    */
    std::string filename_suffix;

    /* basename for write/read files
    * this->db["pod_directory"] + "/" + this->db["pod_basename"]
    * + filename_suffix + "."
    * + this->db["pod_inner_product"] + "." + ("fluc.")
    */
    std::string basename;

    /* rank of pod basis */
    int rank;

    /* dof of pod basis functions */
    int length;

    /* dof of snapshots */
    int length_snaps;

    /* number of snapshots */
    int number_snaps;

    /* threshold value for eigenvalues of autocorrelation matrix */
    double eigen_threshold;

    /* number of eigenvalues greater than eigen_threshold */
    int valid_eigs;

    /* array for pod eigenvalues */
    double* eigs;

    /* pointer to inner product matrix for pod computation */
    std::shared_ptr<TMatrix> gramian_ptr;

    // inner product matrix for time derivative projection
    BlockFEMatrix projection_matrix;

    // solver database
    ParameterDatabase projection_solver_db;

    /* matrix with snapshots values (DO WE NEED IT? TOO MUCH OF STORAGE?) */
    /* one snapshot per row */
    DenseMatrix snaps_mat;

    /* time_mean of the snapshots */
    std::vector<double> snaps_mean;

    /* matrix for pod basis */
    /* one basis-vector per row */
    DenseMatrix pod_basis;

    /**
    * @brief read snapshots
    *
    * Read snapshots from file
    * this->db["snaps_directory"] + "/" + this->db["snaps_basename"] (+ ".bin").
    * This function is called by the compute_basis().
    * The snapshots will be stored into the class member this->snaps_mat, one
    * snapshot per row.
    */
    void read_snapshots();

    /**
    * @brief read basis from file
    *
    * Read POD basis from the file with the name basename + "pod"
    * and store it in the member variable pod_basis, one basis-vector per row.
    * If this->db["pod_fluctuations_only"] is set to true, then the function
    * read_averages() will be automatically called.
    *
    * @param rom_rank the number of POD-modes to be read, in order to build
    *                 a POD-basis of rank 'rom_rank' for the ROM computations.
    *                 If rom_rank <= 0 all possible POD modes will be used.
    */
    void read_basis(int rom_rank=0);

    /**
    * @brief change the basename
    *
    * Enable to change the name basename + filename_suffix, usefull in case of
    * model with many variables like TNSE.
    */
    void set_basename(std::string filename_suffix);

  private:

    /**
    * @brief Decompose snapshots by substracting their mean values
    *
    * Decompose snapshots into the time average of the snapshots and the
    * fluctuating part of the snapshots stored in this->snaps_mat.
    */
    void decompose_snaps();

    /**
    * @brief Write POD data into file
    *
    * Write POD basis functions(row-wise) into basename + "pod", and call
    * write_eigenvalues().
    * If this->db["pod_fluctuations_only"]==true call write_averages().
    */
    void write_pod();

    /**
    * @brief Write time average of snapshots into file
    *
    * Write time average of snapshots into basename + "mean".
    * It is automatically called from write_pod() if
    * db["pod_fluctuations_only"]==true.
    *
    * @param[in] basename basename of the file
    */
    void write_averages();

    /**
    * @brief Write POD eigenvalue data into file
    *
    * Write time steps, POD eigenvalues and missing energy ratios into
    * basename + "mean". It is automatically called from write_pod().
    *
    */
    void write_eigenvalues();

    /**
    * @brief Read eigenvalues from file and print the sum of resting eigenvalues
    *
    * Read eigenvalues from the file: basename + "eigs".
    * It is automatically called  at the end of read_basis() if
    * this->rank < this->valid_eigs.
    *
    */
    void read_eigenvalues();

    /**
    * @brief Read time average of snapshots from file
    *
    * Read snapshots mean from the file: basename + "mean".
    * It is automatically called from read_basis() if
    * db["pod_fluctuations_only"]==true.
    * The average of the snapshots is stored in the member variable snaps_mean.
    *
    */
    void read_averages();

    /**
    * @brief read data from file
    *
    * Read data from the file _filename and save it into the vector 'data'
    *
    * @param[in]  filename full name of the file
    * @param[out] data     vector
    * @param[in]  format   if 0 all the data from file will be read
    *                      if i only the data from the i-th column
    */
    void read_data(std::string filename, std::vector<double>& data,
                   int format=0);

    /**
    * @brief Compute auto-correlation matrix for eigenvalue problem
    *
    * Compute autocorrelation matrix for POD computation which is computed as
    * U^T*S*U, where U is the snapshot matrix (this->snaps_mat) and S is the
    * inner product matrix (gramian_mat). If S is not specified (i.e.
    * set_gramian() was not called), then the POD basis will be computed with
    * respect to the euclidian inner product with S = Id.
    *
    * @param[out] corr_mat resulting auto-correlation matrix
    */
    void compute_autocorr_mat( DenseMatrix& corr_mat );
};

#endif // POD_H
