/** ****************************************************************************
*
* @name       ROM
* @brief      This class provides different routines for basis transformation
*             (between the underlying finite element and POD space),
*             solving the ROM system.
*
* @date       14/03/2017 (start of implementation)
*
*******************************************************************************/

#ifndef ROM_H
#define ROM_H

#include <fstream>
#include <Matrix.h>
#include <memory>
#include <MooNMD_Io.h>
#include <ParameterDatabase.h>
#include <POD.h>


class ROM : public POD
{
  public:
    /**
    * @brief default constructor
    */
    ROM();

    /**
    * @brief constructor
    */
    ROM(const ParameterDatabase& param_db, int rom_rank=0,
        std::string filename_suffix="");

    /**
    * @brief default destructor
    */
    ~ROM()=default;

    /**
     * Creates a database filled with default parameters. This database will
     * contain all necessary parameters for the ROM.
     */
    static ParameterDatabase default_rom_database();

    /**
    * @brief Reduce full sparse matrix using POD basis
    *
    * red_mat = pod_basis^T * full_mat * pod_basis
    *
    * @param[in]  full_mat Sparse ParMooN matrix to reduce
    * @param[out] red_mat  Reduced matrix to be computed
    */
    void reduce(std::shared_ptr<TMatrix> full_mat,
                DenseMatrix&             red_mat);

    /**
    * @brief Reduce full sparse matrix using two different POD basis
    *        for Navier Stokes: P_basis is this basis and
    *                           U_basis is other_basis
    *
    * red_mat = pod_basis^T * full_mat * other_basis
    *
    * @param[in]  full_mat    Sparse ParMooN matrix to reduce
    * @param[out] red_mat     Reduced matrix to be computed
    * @param[in]  other_basis other POD basis
    */
    void reduce(std::shared_ptr<TMatrix> full_mat,
                DenseMatrix&             red_mat,
                const DenseMatrix* const other_basis);

    /**
    * @brief Reduce full-order sparse matrix by multiplying by the POD basis and
    *        the snapshots mean
    *
    * red_vec = pod_basis^T * full_mat * snaps_mean
    *
    * @param[in]  full_mat Sparse ParMooN matrix
    * @param[out] red_vec  Reduced array to be computed
    *
    */
    void reduce_mat_mean(std::shared_ptr<TMatrix> full_mat,
                         std::vector<double>&     red_vec);

    /**
    * @brief Reduce full-order sparse matrix by multiplying by the POD basis and
    *        an other snapshots mean
    *        for Navier Stokes: P_basis is this basis and
    *                           U_snaps_mean is other_snaps_mean
    *
    * red_vec = pod_basis^T * full_mat * other_snaps_mean
    *
    * @param[in]  full_mat Sparse ParMooN matrix
    * @param[out] red_vec  Reduced array to be computed
    *
    */
    void reduce_mat_mean(std::shared_ptr<TMatrix> full_mat,
                         std::vector<double>&     red_vec,
                         const double* const      other_snaps_avr);

    /**
    * @brief Reduce full-order vector using POD basis
    *
    * red_vec = pod_basis^T * full_vec
    *
    * @param[in]  first    iterator pointing to the fisrt element of the vector
    *                      full_vec
    * @param[in]  last     iterator pointing to the fisrt element of the vector
    *                      full_vec
    * @param[out] red_vec  Reduced-order vector
    *
    */
    void reduce(std::vector<double>::const_iterator first,
                std::vector<double>::const_iterator last,
                std::vector<double>&                red_vec);

    /**
    * @brief Reduce full-order vector using POD basis
    *
    * red_vec = pod_basis^T * full_vec
    *
    * @param[in]  full_vec Full-order vector
    * @param[out] red_vec  Reduced-order vector
    *
    */
    void reduce(const std::vector<double>* const full_vec,
                std::vector<double>&       red_vec);

    /**
    * @brief Reduce initial full-order solution
    *
    * Reduce initial full-order solution by first subtracting from it the
    * snapshots' mean (rom_db["pod_fluctuations_only"] == true) and then
    * projecting it onto POD basis
    *
    * red_sol = pod_basis^T * Gram_mat * (full_sol - mean_snaps),
    *
    * where Gram_mat represents the inner product used for the computation
    * of POD basis.
    * NOTE: For the pointer full_sol enough memory must be allocated.
    *
    * @param[in]  full_sol Finite element solution
    * @param[out] red_sol  ROM solution
    *
    */
    void reduce_init_solution(double* full_sol, std::vector<double>& red_sol);

    /**
    * @brief Compute full-order solution from the ROM solution
    *
    * Compute full-order solution form the ROM solution by first multiplying
    * the ROM solution with the POD matrix and then adding to it the snapshots'
    * mean if rom_db["pod_fluctuations_only"] == true.
    *
    * full_sol = mean_snaps + pod_basis * red_sol
    *
    * NOTE: For the pointer full_sol enough memory must be allocated.
    *
    * @param[in]  red_sol  ROM solution
    * @param[out] full_sol Finite element solution
    *
    */
    void get_full_solution(const std::vector<double>& red_sol,
                           double*                    full_sol);
};

#endif // ROM_H
