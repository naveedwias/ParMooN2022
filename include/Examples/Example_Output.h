#ifndef __EXAMPLE_OUTPUT_VARIABLES__
#define __EXAMPLE_OUTPUT_VARIABLES__

#include <vector>
#include <string>

namespace ExampleOutput
{
  void DeclareOutputVariable(std::string name, double default_value = 0.0);
  void UpdateOutputVariable(std::string name, double value);

  bool HasOutputVariable(std::string name);
  double GetOutputVariable(std::string name);
  std::vector<std::string> ListOutputVariables();
}

#endif