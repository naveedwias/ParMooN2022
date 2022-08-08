/**
 * @file This file holds routines used to compute different values of interess
 *       for the example TNSE_3D/ChannelTau.h
 */

#ifndef _TNS_ChannelTau_
#define _TNS_ChannelTau_

#include <TimeNavierStokes.h>

#include <deque>
#include <vector>

class TDomain;
class ParameterDatabase;

// TODO: - enable choice of streamwise direction and derivative direction
//       - test different finite elements than C_Q2_3D_H_A, C_Q2_3D_H_M,
//         C_Q1_3D_H_A, C_Q1_3D_H_M


/** @brief check if vector V is containing at least an element equals to x
 *         up to the tolerance */
bool is_containing(const std::vector<double>& V,
                   const double x,
                   double tolerance=1e-12);

/** @brief eliminates all except the first element from every consecutive group
*         of equivalent elements up to the tolerance from the range
*         [first, last) and returns a past-the-end iterator for the new
*         logical end of the range
*/
std::vector<double>::iterator unique_double(std::vector<double>::iterator first,
                                            std::vector<double>::iterator last,
                                            double tolerance=1e-12);


class TNS_ChannelTau : public TimeNavierStokes<3>
{
  public:
    /** @brief following values are used in ChannelTau.h to maintain the flow */
    // expected streamwise bulk velocity
    static double bulkVelocity_expected;
    // computed streamwise bulk velocity
    static double bulkVelocity_sim;

    /** @brief set-get methods to communicate with the example ChannelTau.h */
    static void set_bulkVelocityExpected(double value)
    { bulkVelocity_expected = value; }
    static void set_bulkVelocity_sim(double value)
    { bulkVelocity_sim = value; }
    static double get_bulkVelocity_sim()
    { return bulkVelocity_sim; }

    /** @brief The standard constructor
     *  @param[in] Domain   The computational domain providing the grid
     *  @param[in] param_db A parameter database with parameters concerning this
     *                      class or any of its members
     */
    TNS_ChannelTau(const TDomain&           domain,
                   const ParameterDatabase& param_db);

    /** @brief compute errors and write solution */
    void output();


  protected:
    // Reynolds number
    double reynolds_number;

    /** @brief coordinates x, y and z of the degrees of freedom */
    std::vector<double> xDofs;
    std::vector<double> yDofs;
    std::vector<double> zDofs;
    // array with coordinates of dof-layers in z-direction
    std::vector<double> zLayers;
    // number of dof-layers in z-direction
    size_t nZLayers = 0;

    // number of dof per layer, used for averaging
    std::vector<int> sum_layer_dofs;

    // tolerance used for the comparison between double
    const double tolerance = 1e-10;

    // as long as the value is <0 the time averaging is not started,
    // the time averaging starts when CURRENTTIME >= start_time_averaging_at
    // and t0_avg is set to the actual starting time for time averaging
    double t0_avg = -1.;

    /** @brief mean velocities on layers
     * only spatial mean for CURRENTTIME < start_time_averaging_at
     * time averaging of spatial mean for CURRENTTIME >= start_time_averaging_at
     */
    std::deque<std::vector<double>> MeanVelocity;

    /** @brief mean Reynold stress on layers
     * [0]:u1u1, [1]:u2u2, [2]:u3u3, [3]:u1u2, [4]:u1u3, [5]:u2u3
     * only spatial mean for CURRENTTIME < start_time_averaging_at
     * time averaging of spatial mean for CURRENTTIME >= start_time_averaging_at
     */
    std::deque<std::vector<double>> MeanReynoldsStress;

    /** @brief z-derivative of the mean velocity
     * only the first component of velocity (streamwise)
     */
    std::vector<double> DerivMeanVelocity;

    void check_and_set_parameters();

    /** @brief computing the coordinates of the d.o.f. */
    void GetCoordinatesOfDof();

    /** @brief set the memory */
    void set_up_memory();

    /** @brief compute the summation of layers for averaging */
    void count_dofs_per_layer(std::vector<int>& summ);

    /** @brief compute mean velocity and the corresponding mean Reynolds stress
     * if start_time_averaging_at>CURRENTTIME only the spatial mean is computed
     * otherwise the time averaging of the spatial mean is computed.
     */
    void computeMeanValues();

    /** @brief compute the summation of velocity components and Reynolds Stress
     *         on layers, used for the computation of the spatial mean values
     */
    void compute_sum_u_uiuj_on_layers(
                                  std::deque<std::vector<double>>& velocity_sum,
                                  std::deque<std::vector<double>>& uiuj_mean);

    /** @brief compute the bulk Velocity (total average streamwise velocity) */
    void compute_bulkVelocity(const std::vector<double>& spatialVMean);

    /** @brief compute the temporal mean,
     *         is only called when t0_avg>=start_time_averaging_at
     *                             (and CURRENTTIME>t0_avg)
     *  @param[in] spatial mean of the current time step
     *  @param[in/out] temporal mean of the previous time step will be
     *                 updated to the current time
     */
    void temporalMean(std::vector<double>  spatial_mean,
                      std::vector<double>& temporal_mean);

    /** compute the z-derivative of values */
    void compute_derivative(const std::vector<double>& values,
                            std::vector<double>&       derivative);

    /** @brief print out the interesting quantities */
    void print_quantity_of_interest();

    /** @brief computes the friction viscosity u_tau */
    double get_FrictionVelocity();

    /** @brief compute the root mean square velocities
     *        (formula 11 Volker's paper)
     *  @return root mean square velocity all components
     */
    std::deque<std::vector<double>> get_rms();
};

#endif
