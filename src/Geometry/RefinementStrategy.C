// =============================================================================
//
// Purpose: Implementation of RefinementStrategy.h.
//
// Author: Moritz Hoffmann, 06/07/15.
//
// =============================================================================

#include "RefinementStrategy.h"
#include "BaseCell.h"
#include "Collection.h"
#include "MooNMD_Io.h"

template <int d>
ParameterDatabase RefinementStrategy<d>::default_refinement_strategy_database()
{
  ParameterDatabase db("default refinement strategy database");
  
  db.add("refinement_strategy", 0u, "Choose the strategy for refinement, see "
         "the enum class 'RefinementStrategyType'.", 0u, 2u);
  db.add("refinement_tolerance", 0.5, 
         "Fraction to determine which cells to refine. Depending on the "
         "'refinement_strategy', typically, cells whose local error is larger "
         "than this tolerance times the maximal global error are refined. The "
         "smaller the value, the more cells are being refined.", 0., 1.);
  db.add("coarsening_tolerance", 0.0, 
         "Fraction to determine which cells to coarsen. Depending on the "
         "'refinement_strategy', typically, cells whose local error is smaller "
         "than this tolerance times the maximal global error are coarsened. "
         "The larger the value, the more cells are being coarsened.", 0., 1.);
  db.add("min_fraction_to_change", 0.1, 
         "The minimum number of cells to be refined/coarsened as a fraction of "
         "the total number of cells. If too few cells are refined/coarsened "
         "then the 'refinement_tolerance' ('coarsening_tolerance') are changed "
         "according to 'decrease_reftol_factor' ('increase_coarsetol_factor'). "
         "This avoids adaptive refinements with only very few changed cell "
         "which would be inefficient.");
  db.add("decrease_reftol_factor", 0.8, 
         "Whenever too few cells are to be refined/coarsened (see "
         "'min_fraction_to_change'), the 'refinement_tolerance' is multiplied "
         "by this factor.", 0., 1.);
  db.add("increase_coarsetol_factor", 1.1,
         "Whenever too few cells are to be refined/coarsened (see "
         "'min_fraction_to_change'), the 'coarsening_tolerance' is multiplied "
         "by this factor.", 1., 1.e10);
  db.add("fraction_of_error", 0.25,
         "The sum of local errors of refined cells should be smaller than this "
         "fraction of the global error. This only applies to the "
         "'refinement_strategy' PortionOfEstimatedGlobalError (1).", 0., 1.);
  db.add("max_cell_level", 1000, 
         "The maximal refinement leven for each cell, i.e., no cell is refined "
         "more often than this number.", 1, 1000000);
  // why is this not in the respective ErrorEstimator?
  db.add("current_estimator", 0, 
         "type of error estimator, depends on problem type "
         "(convection-diffusion, Navier-Stokes, ...).", 0, 10);
  return db;
}

template <int d>
RefinementStrategy<d>::RefinementStrategy(const ParameterDatabase& db_in)
 : refinements() // filled during 'RefinementStrategy::apply_estimator'
{
  auto db(RefinementStrategy::default_refinement_strategy_database());
  db.merge(db_in, false);
  refine_strategy = RefinementStrategyType((size_t)db["refinement_strategy"]);
  reftol = db["refinement_tolerance"];
  coarsetol = db["coarsening_tolerance"];
  min_fraction_to_change = db["min_fraction_to_change"];
  decrease_reftol_factor = db["decrease_reftol_factor"];
  increase_coarsetol_factor = db["increase_coarsetol_factor"];
  fraction_of_error = db["fraction_of_error"];
  max_cell_level = db["max_cell_level"];
  current_estimator = db["current_estimator"];
  
  if(coarsetol >= reftol)
    this->coarsetol = 0.001 * reftol;
  if(coarsetol != 0.)
  {
    ErrThrow("Coarsening during adaptive grid refinement is not yet supported. "
             "Please set 'coarsening_tolerance' to 0.");
  }
}

template <int d>
bool RefinementStrategy<d>::should_refine_cell(size_t cellIndex) const
{
  if(refinements.empty())
    ErrThrow("please call 'RefinementStrategy::apply_estimator' before using "
             "'RefinementStrategy::should_refine_cell'");
  return refinements.at(cellIndex);
}

template <int d>
void RefinementStrategy<d>::apply_estimator(ErrorEstimator<d> &estimator)
{
  // get the collection of the finest available mesh
  auto collection = estimator.get_collection();
  // number of cells in that finest mesh
  size_t n_cells = collection->GetN_Cells();
  // resize indicator vector accordingly
  refinements.resize((unsigned long) n_cells);
  // initialize it with false, i.e., do not refine anything
  std::fill(refinements.begin(), refinements.end(), false);

  // number of marked cells, current iteration and maximal number of iterations
  int changed = -1, it = 0, max_it = 100;
  // maximal local error
  double eta_max = estimator.get_maximal_local_error();
  Output::print<1>("Refinement strategy eta_max: ", eta_max);
  // reference tolerance as percentage value of the maximal local error
  double reftol = this->reftol * eta_max;
  // reference coarsening tolerance as percentage value of the maximal local error
  double coarsetol = this->coarsetol * eta_max;
  // minimal number of cells to change
  double min_changed = n_cells * min_fraction_to_change;
  // current cell
  TBaseCell *cell;

  Output::print("Refine strategy: ", refine_strategy, ", min cells to change: ",
         ((int) min_changed));
  switch (refine_strategy)
  {
    case RefinementStrategyType::MaximalEstimatedLocalError:
    {
      // loop while the number of changed cells is below the minimal number and we have not reached max_it
      while ((changed < min_changed) && (it < max_it))
      {
        // reset number of changed cells
        // if this is not the first pass, the conditions have been relaxed and the cells that have been marked
        // will now be marked in particular
        changed = 0;
        // loop over all cells
        for (size_t m = 0; m < n_cells; m++)
        {
          cell = collection->GetCell((int) m);
          // mark cell for refinement
          if ((estimator.get_eta_K(m) >= reftol) && (cell->GetGeoLevel() < max_cell_level))
          {
            refinements[m] = true;
            changed++;
          }
          if (estimator.get_eta_K(m) <= coarsetol)
          {
            ; //TODO: mark cell for coarsening
          }
        }
        // not enough cells marked, relax conditions
        if (changed < min_changed) {
            reftol *= decrease_reftol_factor;
            coarsetol *= increase_coarsetol_factor;
        }
        Output::print<4>("total ", n_cells, " changed ", changed);
        it++;
      }
      break;
    }
    case RefinementStrategyType::PortionOfEstimatedGlobalError:
    {
      // mark smallest set of cells whose sum of local errors is a prescribed
      // fraction of the global
      double sum_local_errors = 0.0;
      // sum of local errors of marked cells
      double eta = estimator.get_estimated_global_error();
      // because we compare to the sum over the eta_K^2
      eta *= eta;
      // if number of changed cells is smaller than min_changed and the sum of local errors is smaller than a
      // fraction of the global error, refine more cells
      while ((changed < min_changed) && (sum_local_errors <= fraction_of_error * eta) && (it < max_it))
      {
        // reset number of changed cells
        // if this is not the first pass, the conditions have been relaxed and the cells that have been marked
        // will now be marked in particular
        changed = 0;
        // loop over all cells
        for (size_t m = 0; m < n_cells; m++)
        {
          // consider m-th cell
          cell = collection->GetCell((int) m);
          // mark cell for refinement
          if ((estimator.get_eta_K(m) >= reftol) && (cell->GetGeoLevel() < max_cell_level))
          {
            refinements[m] = true;
            // accumulate sum_local_errors
            sum_local_errors += estimator.get_eta_K()[m];
            // number of changed cells increased
            changed++;
          }
          if (estimator.get_eta_K(m) <= coarsetol)
          {
            ; //TODO: mark cell for derefinement
          }
        }
        Output::print<2>("total ", n_cells, " changed ", changed,
                         " global error ", eta, " error of refined cells ",
                         sum_local_errors);
        // if the number of marked cells is still smaller than min_changed or the accumulated errors are still
        // smaller than a presribed fraction of the global error, relax conditions
        if ((changed < min_changed) || (sum_local_errors <= fraction_of_error * eta))
        {
          reftol *= decrease_reftol_factor;
          coarsetol *= increase_coarsetol_factor;
          sum_local_errors = 0.0;
        }
        it++;
      }
      break;
    }
    case RefinementStrategyType::EquilibrationStrategy:
    {
      double sum_local_errors = 0.0;
      double eta = estimator.get_estimated_global_error();
      eta *= eta;
      while(sum_local_errors < (1.-reftol)*eta)
      {
        // find maximal eta_K
        double eta_max_value = -1;
        {
          for (size_t m = 0; m < n_cells; ++m)
          {
            if(!refinements[m] && estimator.get_eta_K(m) > eta_max_value
              && (collection->GetCell((int) m)->GetGeoLevel() < max_cell_level))
            {
              eta_max_value = estimator.get_eta_K(m);
            }
          }
        }
        // prevent infinite loop
        if(eta_max_value == -1) break;
        // mark all cells with value of eta_K
        for (size_t m = 0; m < n_cells; ++m)
        {
          if(!refinements[m] && estimator.get_eta_K(m) == eta_max_value
            && (collection->GetCell((int) m)->GetGeoLevel() < max_cell_level))
          {
            sum_local_errors += eta_max_value;
            changed++;
            refinements[m] = true;
          }
        }
      }
      Output::print<2>("total ", n_cells, " changed ", changed,
                       " global error ", eta, " error of refined cells ",
                       sum_local_errors);
      break;
    }
  }
  Output::print<1>("Cells marked for refinement: ", changed, ", fraction ",
                   double(changed)/n_cells);
  if(it >= max_it && changed < min_changed)
  {
    Output::print<1>("Failed to mark ", int(min_changed),
                     " cells, since parameters were relaxed it=", it,
                     " >= max_it times. Consider decreasing ",
                     "decrease_reftol_factor.");
  }
}

template <int d>
void RefinementStrategy<d>::apply_indicator(
  const TCollection* collection,
  std::function<double(std::array<double, 3>)> indicator)
{
  size_t n_cells = collection->GetN_Cells();
  // resize indicator vector accordingly
  refinements.resize((unsigned long) n_cells);
  // initialize it with false, i.e., do not refine anything
  std::fill(refinements.begin(), refinements.end(), false);
  
  for (size_t i = 0; i < n_cells; ++i)
  {
    auto curr_cell = collection->GetCell(i);
    int n_inner = 0;
    int n_outer = 0;
    int k = curr_cell->GetN_Vertices();
    for(int j = 0; j < k; ++j)
    {
      auto vert = curr_cell->GetVertex(j);
      double val = indicator({vert->GetX(), vert->GetY(), vert->GetZ()});
      if(val <= 1e-10) n_inner++;
      if(val >= -1e-10) n_outer++;
    }
    if((n_inner > 0) && (n_outer > 0)) 
    {
      // there vertices on both sides
      refinements[i] = true;
    }
  }
}


#ifdef __3D__
template class RefinementStrategy<3>;
#else
template class RefinementStrategy<2>;
#endif
