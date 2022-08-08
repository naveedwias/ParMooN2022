/** ****************************************************************************
*
* @name       ROM
* @brief      This class provides different routines for basis transformation
*             (between the underlying finite element and POD space),
*             solving the ROM system.
*
*******************************************************************************/

#include <ROM.h>

#include <ParameterDatabase.h>

/**************************************************************************** */
ParameterDatabase ROM::default_rom_database()
{
  ParameterDatabase db("Default ParMooN parameter database for "
                       "POD-based ROM problems");

  db.add("mat_time_dependent", false,
         "This is the flag indicating if stiffness and mass matrices are time "
         "dependent and then if they should be reassembled at every time step.",
         {true,false});

  db.add("rom_init_regularized", false,
         "This is the flag whether the ROM initial condition should be "
         "regularized.",
         {true,false});

  db.add("differential_filter_width", 1.0,
         "Filter width for the differential filter (Helmoltz equation) for the "
         "computation of the regularized ROM initial condition.",
         0., 10.);

  // Merge with POD database needed to get information about the POD basis
  db.merge(POD::default_pod_database(), true);

  return db;
}

/**************************************************************************** */
ROM::ROM(const ParameterDatabase& param_db, int rom_rank,
         std::string filename_suffix)
: POD(param_db, filename_suffix)
{
  POD::read_basis(rom_rank);
}

/**************************************************************************** */
void ROM::reduce(std::shared_ptr<TMatrix> full_mat,
                 DenseMatrix&             red_mat)
{
  bool transpose = true;
  std::shared_ptr<DenseMatrix> tmp_mat = full_mat->
                                         multiply(POD::get_basis(), transpose);

  red_mat = *(POD::get_basis()->multiply(tmp_mat.get()));
}

/**************************************************************************** */
void ROM::reduce(std::shared_ptr<TMatrix> full_mat,
                 DenseMatrix&             red_mat,
                 const DenseMatrix* const other_basis)
{
  bool transpose = true;
  std::shared_ptr<DenseMatrix> tmp_mat = full_mat->
                                         multiply(other_basis, transpose);

  red_mat = *(POD::get_basis()->multiply(tmp_mat.get()));
}

/**************************************************************************** */
void ROM::reduce_mat_mean(std::shared_ptr<TMatrix> full_mat,
                          std::vector<double>&     red_vec)
{
  if(!POD::db["pod_fluctuations_only"])
  {
    ErrThrow("ROM::reduce_mat_mean(..): Function is only available if the "
             "parameter 'pod_fluctuations_only' is set to true.");
  }

  std::vector<double> tmp_vec (full_mat->get_n_rows());
  full_mat->multiply( &POD::get_snaps_avr()[0], &tmp_vec[0] );

  red_vec = POD::pod_basis.multiply(&tmp_vec);
}


/**************************************************************************** */
void ROM::reduce_mat_mean(std::shared_ptr<TMatrix> full_mat,
                          std::vector<double>&     red_vec,
                          const double* const      other_snaps_avr)
{
  if(!POD::db["pod_fluctuations_only"])
  {
    ErrThrow("ROM::reduce_mat_mean(..): Function is only available if the "
             "parameter 'pod_fluctuations_only' is set to true.");
  }

  std::vector<double> tmp_vec (full_mat->get_n_rows());
  full_mat->multiply( other_snaps_avr, &tmp_vec[0] );

  red_vec = POD::pod_basis.multiply(&tmp_vec);
}

/**************************************************************************** */
void ROM::reduce(std::vector<double>::const_iterator first,
                 std::vector<double>::const_iterator last,
                 std::vector<double>&       red_vec)
{
  red_vec = POD::pod_basis.multiply(first, last);
}

/**************************************************************************** */
void ROM::reduce(const std::vector<double>* const full_vec,
                 std::vector<double>&       red_vec)
{
  red_vec = POD::pod_basis.multiply(full_vec);
}

/**************************************************************************** */
void ROM::reduce_init_solution(double* full_sol, std::vector<double>& red_sol)
{
  std::vector<double> tmp_vec(POD::length);
  for (int i=0; i<POD::length; ++i)
  {
    if(POD::db["pod_fluctuations_only"])
    {
      tmp_vec[i] = full_sol[i] - POD::snaps_mean[i];
    }
    else
    {
      tmp_vec[i] = full_sol[i];
    }
  }

  if(POD::db["pod_inner_product"].get<std::string>() != "euclidean")
  {
    if(POD::gramian_ptr == nullptr)
    {
      ErrThrow("Gramian matrix is not available! Before reducing initial full "
              "solution POD::set_gramian_ptr() must be called!");
    }
    std::vector<double> copy_vec(POD::length);
    copy_vec.swap(tmp_vec);
    POD::gramian_ptr->multiply(&copy_vec[0], &tmp_vec[0]);
  }

  red_sol = POD::pod_basis.multiply(&tmp_vec);
}

/**************************************************************************** */
void ROM::get_full_solution(const std::vector<double>& red_sol,
                            double*                    full_sol)
{
  if((int)red_sol.size() != POD::rank)
  {
    ErrThrow("ROM::get_full_solution(): Length of input reduced-order vector "
             "must coincide with the dimension of the POD basis. Length of "
             "input vector: ", red_sol.size(), ", dimension of POD basis: ",
             POD::rank);
  }

  bool transpose = true;
  std::vector<double> tmp_vec(POD::length);
  tmp_vec = POD::pod_basis.multiply(&red_sol, 1., transpose);

  for (int i=0; i<POD::length; ++i)
  {
    if(POD::db["pod_fluctuations_only"])
    {
      full_sol[i] = tmp_vec[i] + POD::snaps_mean[i];
    }
    else
    {
      full_sol[i] = tmp_vec[i];
    }
  }
}
