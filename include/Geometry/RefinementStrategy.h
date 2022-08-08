// =============================================================================
//
// Purpose: This class provides refinement strategies. Provided eta_K's from
//          an error estimator, cells are marked in a boolean vector which are
//          to be refined.
//
// Author: Moritz Hoffmann, 06/07/15.
//
// =============================================================================

#ifndef PARMOON_REFINEMENTSTRATEGY_H
#define PARMOON_REFINEMENTSTRATEGY_H


#include "ParameterDatabase.h"
#include "Enumerations_geometry.h"
#include "ErrorEstimator.h"
#include <array>
#include <functional>
#include <vector>

template <int d>
class RefinementStrategy 
{
  protected:
    RefinementStrategyType refine_strategy;
    double reftol;
    double coarsetol;
    double min_fraction_to_change;
    double decrease_reftol_factor;
    double increase_coarsetol_factor;
    double fraction_of_error;
    int max_cell_level;
    std::vector<bool> refinements;
    int current_estimator;
  public:
    
    RefinementStrategy(const ParameterDatabase& db);
    ~RefinementStrategy() = default;
    
    /// @brief default parameters used in this class
    static ParameterDatabase default_refinement_strategy_database();

    /// @brief With a given ErrorEstimator find out which cells should be 
    /// refined/coarsened.
    /// This fills the member vector 'refinements'
    void apply_estimator(ErrorEstimator<d> &estimator);
    
    /// @brief Find out which cells should be refined using an indicator.
    /// The indicator function maps coordinates (x,y,z) (z=0 in 2D) to double.
    /// All cells which have vertices with both positive and negative indicator
    /// are marked to be refined. This fills the member vector 'refinements'.
    void apply_indicator(const TCollection* coll,
                         std::function<double(std::array<double, 3>)> indicator);
    
    /// @brief specify explicitly which cells are to be refined. The vector
    /// cell_markers is copied to the member `refinements` and its length should
    /// be the number of cells. This method is mainly for debugging purposes.
    void apply_indicator(std::vector<bool> cell_markers)
    { refinements = cell_markers; }

    /// @brief find out if a given cell (by index) should be refined.
    /// @note You have to call 'apply_estimator' before this makes any sense.
    bool should_refine_cell(size_t cell_index) const;
};


#endif //PARMOON_REFINEMENTSTRATEGY_H
