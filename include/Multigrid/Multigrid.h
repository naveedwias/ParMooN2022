/**
 * @file A multigrid class, which holds the necessary
 * grid information for executing a geometric multigrid iteration.
 *
 * @note When using multigrid you must take care of the matrices being
 * correctly assembled on all levels before calling "cycle()". This might,
 * should some nonlinearity or time-dependency be included
 * include one or more calls of GridTransfer::RestrictFunctionRepeatedly,
 * which informs every grid about a current approximate solution.
 *
 * Here are the (bigger) tasks and functionalities to regain. Most of the work
 * necessary will not amass in this class but elsewhere.
 *
 * @todo TODO Reenable step length control (work in VankaSmoother).
 *
 *
 * @date 2016/05/10
 * @author Clemens Bartsch
 */

#ifndef INCLUDE_MULTIGRID_MULTIGRID_H_
#define INCLUDE_MULTIGRID_MULTIGRID_H_

#include <Chrono.h>
#include <CycleControl.h>
#include <MultigridLevel.h>
#include <ParameterDatabase.h>

#include <list>
#include <vector>

//Forward declarations.
class BlockVector;

enum class MultigridType{ STANDARD , MDML };
MultigridType string_to_multigrid_type(const std::string& code);

class Multigrid
{
  public:

    /** @brief construct a multigrid object.
     * 
     * To really use it, you have to call initialize as well.
     */
    explicit Multigrid(const ParameterDatabase& db);

    /** @brief Actual set up a multigrid object. 
     *
     * Note that there is ABSOLUTELY NO CHECKS performed whether the matrices 
     * found a reasonable Multigrid hierarchie. The crucial point are the 
     * restriction and prolongation methods - everything else will work 
     * for any set of Rubbish matrices.
     *
     * @param matrices A vector of the matrices per level, ordered from coarsest
     * (0) to finest level.
     */
    void initialize(const std::list<BlockFEMatrix*>& matrices);

    /// @brief Apply one complete multigrid cycle. 
    /// @details Which kind of cycle that is is determined at the time of 
    /// construction via a parameter database. So far V,W and F cycle are 
    /// implemented.
    void cycle();

    /// @brief Get the solution on the finest grid.
    /// @details Call this after a cycle is complete to get the result.
    const BlockVector& get_finest_sol();

    /// @brief get the type of this multigrid object
    MultigridType get_type() const { return type_; }
    
    /// @brief find out it mdml is used, otherwise standard multigrid is used
    bool is_using_mdml() const { return type_ == MultigridType::MDML; }

    /// @brief return the number of multigrid levels.
    /// @details This returns the number of geometric levels.
    size_t get_n_geometric_levels() const { return db["multigrid_n_levels"]; };

    /// @brief return the number of multigrid levels.
    /// @details This returns the number of algebraic levels.
    size_t get_n_algebraic_levels() const { return n_algebraic_levels_; };

    /// Set the right hand side on the finest grid. It must of course fit the
    /// matrix stored on the finest grid.
    /// When right hand side and (initial) solution are set on the finest grid,
    /// and a call of "update" was made since the last change of the matrices,
    /// the multigrid object is ready for a cycle.
    void set_finest_rhs(const BlockVector& bv);

    /// Set the (initial) solution on the finest grid. It must of course fit the
    /// matrix stored on the finest grid.
    /// When right hand side and (initial) solution are set on the finest grid,
    /// and a call of "update" was made since the last change of the matrices,
    /// the multigrid object is ready for a cycle.
    void set_finest_sol(const BlockVector& bv);

    /// Calling this method is the sign for updating the smoothers on all levels.
    /// It must be called before every cycle which was preceded by a change in
    /// the matrices (e.g. assembling).
    void update();

    /// Set up a database which contains default values of all control parameters
    /// necessary to control a multigrid cycle.
    static ParameterDatabase default_multigrid_database();

    /// Print the total time spent in solving on the coarsest grid.
    /// This feature is of interest for our current project (September 2016,
    /// ParMooN paper) and might be removed soon.
    void print_coarse_grid_time_total() const;

  private:

    /// @brief Parameter database to store all parameters related to multigrid
    ParameterDatabase db;
    
    size_t n_algebraic_levels_;

    /// A list of the participating levels, ordered from coarsest (0)
    /// to finest level.
    std::vector<MultigridLevel> levels_;

    /// A list of damping factors which are used when updating the solution
    /// on a finer level by adding a coarser level's solution.
    std::vector<double> damp_factors_;

    /// The number of pre-smoothing steps to perform per level before
    /// descending to the next coarser level.
    size_t n_pre_smooths_;

    /// Maximal number of smoother iteration on the coarsest level.
    size_t coarse_n_maxit;

    /// The residual to reach in the coarsest level's equation by smoothing
    /// on the coarsest level.
    double coarse_epsilon;

    /// The number of pre-smoothing steps to perform per level after
    /// ascending from the next coarser level.
    size_t n_post_smooths_;

    /// An object taking care of the order of ascends and descends between the levels.
    CycleControl control_;

    /// the object knows whether it is of standard or MDML type.
    MultigridType type_;

    /// A timer object used to measure time spent in coarse grid solve
    /// This object is of interest for our current project (ParMooN Paper, Sep 2016)
    /// and might be removed after that.
    Chrono coarse_grid_timer;


    /// Restrict defect on level lvl and store it as rhs in the next coarsest level.
    void update_rhs_in_coarser_grid(size_t lvl);

    /// Prolongate solution on level lvl and add it (damped!) to solution on the
    /// next finest level.
    void update_solution_in_finer_grid(size_t lvl);

    /// Nuke the solution on the next coarsest grid. This is done in order
    /// to ensure a zero start iterate when smoothing.
    void set_solution_in_coarser_grid_to_zero(size_t lvl);

    /// Perform the operations necessary on one grid.
    int cycle_step(size_t step, size_t level);


};




#endif /* INCLUDE_MULTIGRID_MULTIGRID_H_ */
