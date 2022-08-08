#ifndef TIMEDISCRETIZATIONS_H
#define TIMEDISCRETIZATIONS_H


/** ************************************************************************
 *
 * @name         TimeDiscretization
 * @brief        This class is useful for time dependent problems. The 
 *               idea is to separate the time stepping schemes and parameters
 *               needed for time dependent problems. This class stores the 
 *               different time stepping schemes like backward-Euler, BDF
 *               schemes, and to include other schemes using Butcher tableau.
 *
 * @author       Naveed Ahmed
 * @History      20. 02. 2017
**************************************************************************/

#include <ParameterDatabase.h>
#include <vector>
#include <memory>
#include <BlockFEMatrix.h>
#include <BlockVector.h>


class TimeDiscretization
{
  protected:
    /// @brief the RungeKuttaTable object
    // std::shared_ptr<RungeKuttaTable> rk;

    /// @brief coeffcients of the BDF schemes
    std::vector<double> bdf_coefficients;

    /// @brief a local parameter database which controls this class
    ParameterDatabase db;

    /// @brief time step legth
    double current_time_step_length;

    /// @brief start time
    double start_time;

    /// @brief end time
    double end_time;

    double get_scale_factor(double tau);
    double last_scale_factor;
    double last_b_factor;

public:
    // 
    explicit TimeDiscretization(const ParameterDatabase& param_db);

    /// @brief This parameter is used for time step counter.
    /// It is set to be zero before the time iteration starts
    /// in the main code and increment is done within the 
    /// the time iteration. 
    unsigned int current_step_;
    
    /// @brief Parameter for counting nonlinear iterations.
    unsigned int nonlin_iteration;
    
    /// @brief current time
    double current_time_;
    
    /// @brief return a default TimeDiscretization parameter database
    /// Using the TimeDiscretization class requires these parameters.
    static ParameterDatabase default_TimeDiscretization_database();
    
    /// @brief This parameter is to be set for scaling and descaling 
    /// of the pressure-velocity blocks: this is set in the main class
    /// and indicates how many blocks have to be scaled during the 
    /// time step and how many during the nonlinear 
    unsigned int n_scale_block;
    
    /// @brief this will be used to indicate the matrices B, BT and C
    /// are linear or nonlinear or solution dependent and are setted 
    // in the class Time_NSE2D
    std::string b_bt_linear_nl;
    
    /// @brief This function scales all B, BT (or C) blocks at the 
    /// very first time and only if the nonlinear iteration counter
    /// is 0. For the bdf schemes,  the scaling will be done at first
    /// two, three ... time steps depending on the order of the scheme
    void scale_descale_all_b_blocks(BlockFEMatrix& matrix);
    
    /// @brief This scales the B, BT or C blocks that are nonlinear. This 
    /// is the special case in the SUPG or Residual based VMS schemes
    /// TODO: combine some of the functions but getting no clue at the moment
    void scale_nl_b_blocks(BlockFEMatrix& matrix);

    /// @brief Thinking to use this function which can be used to perform
    /// the complete time step ????
    void set_time_disc_parameters();

    // a function which takes the matrices, right-hand side and solution vector
    // and perform the ...
    void prepare_timestep(BlockFEMatrix& system_matrix, 
        const BlockFEMatrix& mass_matrix, std::vector<BlockVector>& rhs, 
        const std::vector<BlockVector>& old_solutions);
    //
    void prepare_rhs_from_time_disc(BlockFEMatrix& system_matrix, 
        const BlockFEMatrix& mass_matrix, std::vector<BlockVector>& rhs, 
        const std::vector<BlockVector>& old_solutions);
    
    /// @brief this function will prepare the BlockFEMatrix 
    /// which passes to the solver. In details, the system_matrix
    /// will be scaled depending on the time stepping scheme
    /// and the mass matrix will be added to it. The B-blocks
    /// will be only scaled with the corresponding factor:
    /// @param system_matrix stiffness matrix
    /// @param mass_matrix  the mass matrix
    /// @param nl_b_block the BT and B-blocks are nonlinear
    ///        nl_b_block = 0, 2 means B and BT's blocks have 
    ///        to be scaled 
    ///        nl_b_block = 1 means only the BT's have to be scaled
    void prepare_system_matrix(BlockFEMatrix& system_matrix, 
        const BlockFEMatrix& mass_matrix);
    
    /// @brief This function is used to reset the linear part of the 
    /// A blocks. The nonlinear blocks are re-assembled during the 
    /// nonlinear iteration and then the system matrix will be set
    /// in the prepare_system_matrix function 
    void reset_linear_matrices(BlockFEMatrix& matrix, 
                               const BlockFEMatrix& mass);
    
    /// @brief This is an important parameter for the BDF schemes
    /// which uses for the first step solution that is 
    /// obtained either by the backward-Euler or Crank-Nicolson scheme.
    bool pre_stage_bdf;
    
    /// @brief a method which returns how many old solutions I need
    unsigned int n_old_solutions() const
    {
      /*if(rk)
      {
         return 1;
      }
      else */if(bdf_coefficients.size() > 0)
      {
        return bdf_coefficients.size() - 1;
      }
      else
      {
        return 1;
      }
    }

    //getters
    std::vector<double> get_bdf_coefficients() const
    { return bdf_coefficients; }

    double get_step_length() const
    {
      return current_time_step_length;
    }

    void set_step_length(double tau);

    double get_start_time() const
    {
      return start_time;
    }

    double get_end_time() const
    {
      return end_time;
    }

    bool reached_final_time_step() const
    {
      return current_time_ > end_time - 1e-10;
    }
};

#endif // TIMEDISCRETIZATIONS_H
