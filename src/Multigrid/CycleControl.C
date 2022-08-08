/*
 * CycleControl.C
 *
 *  Created on: May 14, 2016
 *      Author: bartsch
 */


#include <CycleControl.h>
#include <MooNMD_Io.h>

CycleControl:: CycleControl()
: n_levels_(0)
{
   ;
}
CycleControl::CycleControl(const std::string& type_str, size_t n_levels)
{
  cycle_ = string_to_cycle_code(type_str);
  n_levels_ = n_levels;

  std::vector<int> mg_recursions(n_levels);
  // coarsest grid
  if (cycle_ == MGCycle::V)
  {
    std::fill(mg_recursions.begin(), mg_recursions.end(), 1);
  }
  else if (cycle_ == MGCycle::W)
  {
    std::fill(mg_recursions.begin(), mg_recursions.end(), 2);
  }
  else if (cycle_ == MGCycle::F)
  {
    std::fill(mg_recursions.begin(), mg_recursions.end(), 2);
  }
  mg_recursions[n_levels-1] = 1;

  cycle_control_.push_back(MGDirection::Dummy); //start with a dummy

  fill_recursively(mg_recursions, n_levels-1);

}

MGDirection CycleControl::get_direction(size_t step)
{
  return cycle_control_.at(step);
}

void CycleControl::print_cycle_control() const
{
  Output::print("Look at that cycle_control_ !");
  for(auto cc : cycle_control_)
  {
    if(cc == MGDirection::Up)
      Output::print("Up");
    if(cc == MGDirection::Down)
      Output::print("Down");
    if(cc == MGDirection::End)
      Output::print("End");
    if(cc == MGDirection::Dummy)
      Output::print("Dummy");
  }
}

MGCycle CycleControl::string_to_cycle_code(const std::string& code)
{
  if(code == std::string("V"))
    return MGCycle::V;
  else if(code == std::string("W"))
    return MGCycle::W;
  else if(code == std::string("F"))
    return MGCycle::F;
  else
  {
    Output::warn("SmootherCode", "The string ", code,
                 " does not equal a MGCycle code. "
                 "Defaulting to V cycle.");
    return MGCycle::V;
  }
}

void CycleControl::fill_recursively(std::vector<int>& mg_recursions, int level)
{
  int coarsest_level = 0;
  int finest_level = n_levels_ - 1;

  if(level == coarsest_level)
  {
    cycle_control_.push_back(MGDirection::Up);
    return;
  }

  for(int j=0;j<mg_recursions[level];j++)
  {
    cycle_control_.push_back(MGDirection::Down);
    fill_recursively(mg_recursions, level-1);
  }

  if (cycle_ == MGCycle::F)
    mg_recursions[level] = 1;

  if(level == finest_level) //finest level
  {
    cycle_control_.push_back(MGDirection::End);
  }
  else //other levels
  {
    cycle_control_.push_back(MGDirection::Up);
  }

}
