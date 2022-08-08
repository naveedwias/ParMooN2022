#ifndef INCLUDE_GENERAL_PARMOON_H_
#define INCLUDE_GENERAL_PARMOON_H_

#include "ParameterDatabase.h"

namespace parmoon
{
  /**
   * @brief initialize parmoon, call this at the beginning of any main function
   * 
   * This generates global objects which needs to be done exactly once. Also
   * this will read the input file, if present, and return the respective 
   * ParameterDatabase.
   * 
   * Note that a main program can be well formed even without this, but you have
   * to know what you're doing.
   */
  ParameterDatabase parmoon_initialize(int argc = 0, char* argv[] = nullptr);
  
  void parmoon_finalize();
}

#endif // INCLUDE_GENERAL_PARMOON_H_
