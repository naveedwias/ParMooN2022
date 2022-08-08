/*
 * CycleControl.h
 *
 * @ruleof0
 *
 *  Created on: May 14, 2016
 *      Author: bartsch
 */

#ifndef INCLUDE_MULTIGRID_CYCLECONTROL_H_
#define INCLUDE_MULTIGRID_CYCLECONTROL_H_

#include <string>
#include <vector>

enum class MGDirection{Down, Up, End, Dummy};
enum class MGCycle{V, W, F};

class CycleControl
{
  public:
    CycleControl();

    CycleControl(const std::string& control, size_t n_levels_);

    void print_cycle_control() const;

    MGDirection get_direction(size_t step);

    MGCycle get_type()
    {
      return cycle_;
    }

    size_t get_n_steps()
    {
      return cycle_control_.size() - 1;//-1, for first is always a dummy object
    }

  private:

    MGCycle cycle_;

    size_t n_levels_;

    /// This vector visibly controls moving up and down in the grid hierarchy.
    std::vector<MGDirection> cycle_control_;

    MGCycle string_to_cycle_code(const std::string& code);

    void fill_recursively(std::vector<int>& mg_recursions, int level);
};



#endif /* INCLUDE_MULTIGRID_CYCLECONTROL_H_ */
