#include <MooNMD_Io.h>
#include <Example_Output.h>
#include <unordered_map>

std::unordered_map<std::string, double> example_output_variables;
std::vector<std::string> declared_example_output_variable_names;

void ExampleOutput::DeclareOutputVariable(std::string name, double default_value)
{
  if (example_output_variables.count(name) > 0)
  {
    Output::root_warn("ExampleOutput", "Double declaration of output variable '", name, "'");
    return;
  }

  declared_example_output_variable_names.push_back(name);
  example_output_variables[name] = default_value;
}

void ExampleOutput::UpdateOutputVariable(std::string name, double value)
{
  if (example_output_variables.count(name) == 0)
  {
    ErrThrow("Tried to set undeclared output variable '", name, "'");
  }

  example_output_variables[name] = value;
}

bool ExampleOutput::HasOutputVariable(std::string name)
{
  return example_output_variables.count(name) > 0;
}

double ExampleOutput::GetOutputVariable(std::string name)
{
  if (example_output_variables.count(name) == 0)
  {
    ErrThrow("Tried to retrieve undeclared output variable '", name, "'");
  }

  return example_output_variables[name];
}

std::vector<std::string> ExampleOutput::ListOutputVariables()
{
  return declared_example_output_variable_names;
}