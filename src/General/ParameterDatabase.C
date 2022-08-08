#include <ParameterDatabase.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <memory>
#include <ParMooN_repository_info.h>
#include <Utilities.h>
#include <MooNMD_Io.h>


#ifdef _MPI
#include <mpi.h>
#endif

// helper functions to find a parameter or database by their name in a list 
template<class T>
typename std::list<T>::const_iterator 
find_in_list(const std::string& name, const std::list<T>& search_list)
{
  auto search_lambda = [=](const T& p){ return p.get_name() == name; };
  return std::find_if(search_list.begin(), search_list.end(), search_lambda);
}
template<class T>
typename std::list<T>::iterator 
find_in_list(const std::string& name, std::list<T>& search_list)
{
  auto search_lambda = [=](const T& p){ return p.get_name() == name; };
  return std::find_if(search_list.begin(), search_list.end(), search_lambda);
}

bool ParameterDatabase::contains(const std::string& param_name) const
{
  if(find_in_list(param_name, this->parameters) == this->parameters.end())
    return false;
  else
    return true;
}

/* ************************************************************************** */
ParameterDatabase::ParameterDatabase(std::string name) 
 : name(name), parameters(), databases()
{
  
}

/* ************************************************************************** */
ParameterDatabase::ParameterDatabase(const ParameterDatabase& db)
 : name(db.name), parameters(), databases()
{
  // iterate through the parameters in db
  for(const auto& p : db.parameters)
  {
    this->parameters.emplace_back(Parameter(p));
  }
  for(const auto& database : db.databases)
  {
    this->databases.emplace_back(database);
  }
}

/* ************************************************************************** */
ParameterDatabase& ParameterDatabase::operator=(const ParameterDatabase& db)
{
  this->name = db.get_name();
  // remove parameters and databases previously stored
  this->parameters.clear();
  this->databases.clear();
  // add parmeter copies from db to this
  // iterate through the parameters in db
  for(const auto& p : db.parameters)
  {
    this->parameters.emplace_back(Parameter(p));
  }
  for(const auto& database : db.databases)
  {
    this->databases.emplace_back(database);
  }
  return *this;
}

/* ************************************************************************** */
template <typename T>
void ParameterDatabase::add(std::string name, T value, std::string description)
{
  Parameter p(name, value, description);
  this->add(std::move(p));
}
template void ParameterDatabase::add(std::string n, bool v,      std::string d);
template void ParameterDatabase::add(std::string n, int v,       std::string d);
template void ParameterDatabase::add(std::string n, size_t v,    std::string d);
template<> void ParameterDatabase::add(std::string n, unsigned int v,
                                       std::string d)
{
  this->add(n, static_cast<size_t>(v), d);
}
template void ParameterDatabase::add(std::string n, double v,    std::string d);
template void ParameterDatabase::add(std::string n,std::string v,std::string d);
template<> void ParameterDatabase::add(std::string n,const char* v,
                                       std::string d)
{
  this->add(n, std::string(v), d);
}

/* ************************************************************************** */
template <typename T>
void ParameterDatabase::add(std::string name, std::vector<T> value,
                            std::string description)
{
  Parameter p(name, value, description);
  this->add(std::move(p));
}
template void ParameterDatabase::add(std::string n, std::vector<bool> v,
                                     std::string d);
template void ParameterDatabase::add(std::string n, std::vector<int> v,
                                     std::string d);
template void ParameterDatabase::add(std::string n, std::vector<size_t> v,
                                     std::string d);
template void ParameterDatabase::add(std::string n, std::vector<double> v,
                                     std::string d);
template void ParameterDatabase::add(std::string n, std::vector<std::string> v,
                                     std::string d);
template <typename T>
void ParameterDatabase::add(std::string name, std::initializer_list<T> value,
                            std::string description)
{
  Parameter p(name, value, description);
  this->add(std::move(p));
}
template void ParameterDatabase::add(std::string n,
                                     std::initializer_list<bool> v, 
                                     std::string d);
template void ParameterDatabase::add(std::string n,
                                     std::initializer_list<int> v,
                                     std::string d);
template void ParameterDatabase::add(std::string n,
                                     std::initializer_list<size_t> v,
                                     std::string d);
template<> void ParameterDatabase::add(std::string n,
                                     std::initializer_list<unsigned int> v,
                                     std::string d)
{
  std::vector<size_t> vals(v.begin(), v.end());
  this->add(n, vals, d);
}
template void ParameterDatabase::add(std::string n,
                                     std::initializer_list<double> v,
                                     std::string d);
template void ParameterDatabase::add(std::string n,
                                     std::initializer_list<std::string> v,
                                     std::string d);
/* ************************************************************************** */
template <typename T>
void ParameterDatabase::add(std::string name, T value, std::string description,
                            T min, T max)
{
  Parameter p(name, value, description);
  p.set_range(min, max);
  this->add(std::move(p));
}
template void ParameterDatabase::add(std::string n, int v,    std::string d,
                                     int min, int max);
template void ParameterDatabase::add(std::string n, size_t v, std::string d,
                                     size_t min, size_t max);
template <>
void ParameterDatabase::add(std::string n, unsigned int v, std::string d,
                            unsigned int min, unsigned int max)
{
  this->add<size_t>(n, v, d, min, max);
}
template void ParameterDatabase::add(std::string n, double v, std::string d,
                                     double min, double max);

/* ************************************************************************** */
template <typename T>
void ParameterDatabase::add(std::string name, std::vector<T> value, 
                            std::string description, T min, T max)
{
  Parameter p(name, value, description);
  p.set_range(min, max);
  this->add(std::move(p));
}
template void ParameterDatabase::add(std::string n, std::vector<int> v,   
                                     std::string d, int min, int max);
template void ParameterDatabase::add(std::string n, std::vector<size_t> v,
                                     std::string d, size_t min, size_t max);
template<> void ParameterDatabase::add(std::string n, std::vector<unsigned int> v,
                                       std::string d, unsigned int min,
                                       unsigned int max)
{
  std::vector<size_t> val(v.begin(), v.end());
  this->add(n, val, d, (size_t)min, (size_t)max);
}
template void ParameterDatabase::add(std::string n, std::vector<double> v,
                                     std::string d, double min, double max);

/* ************************************************************************** */
template <typename T>
void ParameterDatabase::add(std::string name, T value, std::string description, 
                            std::set<T> range)
{
  Parameter p(name, value, description);
  p.set_range(range);
  this->add(std::move(p));
}
template void ParameterDatabase::add(std::string n, bool v,      std::string d,
                                     std::set<bool> range);
template void ParameterDatabase::add(std::string n, int v,       std::string d,
                                     std::set<int> range);
template void ParameterDatabase::add(std::string n, size_t v,    std::string d,
                                     std::set<size_t> range);
template void ParameterDatabase::add(std::string n,std::string v,std::string d,
                                     std::set<std::string> range);
template<> void ParameterDatabase::add(std::string n,const char* v,
                                       std::string d, 
                                       std::set<const char*> range)
{
  std::set<std::string> range_of_strings;
  for(const auto& r : range) range_of_strings.insert(std::string(r));
  this->add(n, std::string(v), d, range_of_strings);
}

/* ************************************************************************** */
template <typename T>
void ParameterDatabase::add(std::string name, std::vector<T> value, 
                            std::string description, std::set<T> range)
{
  Parameter p(name, value, description);
  p.set_range(range);
  this->add(std::move(p));
}
template void ParameterDatabase::add(std::string n, std::vector<bool> v,
                                     std::string d, std::set<bool> range);
template void ParameterDatabase::add(std::string n, std::vector<int> v,
                                     std::string d, std::set<int> range);
template void ParameterDatabase::add(std::string n, std::vector<size_t> v,
                                     std::string d, std::set<size_t> range);
template void ParameterDatabase::add(std::string n, std::vector<std::string> v,
                                     std::string d, std::set<std::string> range);

/* ************************************************************************** */
void ParameterDatabase::add(Parameter && p)
{
  // check if such a parameter already exists in this database
  if(!this->contains(p.get_name()))
    this->parameters.emplace_back(p);
  else
    ErrThrow("parameter with this name already exists ", p.get_name());
}

/* ************************************************************************** */
void ParameterDatabase::add(const Parameter & p)
{
  // check if such a parameter already exists in this database
  if(!this->contains(p.get_name()))
    this->parameters.emplace_back(p);
  else
    ErrThrow("parameter with this name already exists ", p.get_name());
}

/* ************************************************************************** */
const Parameter& ParameterDatabase::operator[](std::string parameter_name) const
{
  auto it = find_in_list(parameter_name, this->parameters);
  if(it == this->parameters.end())
    ErrThrow("unknown parameter ", parameter_name);
  return *it;
}
Parameter& ParameterDatabase::operator[](std::string parameter_name)
{
  auto it = find_in_list(parameter_name, this->parameters);
  if(it == this->parameters.end())
    ErrThrow("unknown parameter ", parameter_name);
  return *it;
}

/* ************************************************************************** */

template <typename T>
T ParameterDatabase::try_get_value(std::string parameter_name, T default_value) const
{
  auto it = find_in_list(parameter_name, this->parameters);
  return it == this->parameters.end() ? default_value : it->get<T>();
}

template bool ParameterDatabase::try_get_value(std::string parameter_name, bool default_value) const;
template int ParameterDatabase::try_get_value(std::string parameter_name, int default_value) const;
template size_t ParameterDatabase::try_get_value(std::string parameter_name, size_t default_value) const;

template double ParameterDatabase::try_get_value(std::string parameter_name, double default_value) const;
template std::string ParameterDatabase::try_get_value(std::string parameter_name, std::string default_value) const;

/* ************************************************************************** */
std::string ParameterDatabase::get_name() const
{
  return this->name;
}

/* ************************************************************************** */
void ParameterDatabase::set_name(const std::string& new_name)
{
  this->name = new_name;
}

/* ************************************************************************** */
size_t ParameterDatabase::get_n_parameters() const
{
  return this->parameters.size();
}

/* ************************************************************************** */
void ParameterDatabase::add_nested_database(ParameterDatabase db)
{
  if(find_in_list(db.get_name(), this->databases) == this->databases.end())
  {
    this->databases.emplace_back(db);
  }
  else
  {
    ErrThrow("additional database with this name already exists ", 
             db.get_name());
  }
}

/* ************************************************************************** */
bool ParameterDatabase::has_nested_database(const std::string& name) const
{
  auto it = find_in_list(name, this->databases);
  if(it == this->databases.end())
    return false;
  else
    return true;
}

/* ************************************************************************** */
const ParameterDatabase & ParameterDatabase::get_nested_database(
  const std::string& name) const
{
  auto it = find_in_list(name, this->databases);
  if(it == this->databases.end())
    ErrThrow("unknown nested database ", name);
  return *it;
}

/* ************************************************************************** */
ParameterDatabase & ParameterDatabase::get_nested_database(
  const std::string& name)
{
  auto it = find_in_list(name, this->databases);
  if(it == this->databases.end())
    ErrThrow("unknown nested database ", name);
  return *it;
}

/* ************************************************************************** */
size_t ParameterDatabase::get_n_nested_databases() const
{
  return this->databases.size();
}

/* ************************************************************************** */
void ParameterDatabase::write(std::ostream& os, bool verbose,
                              bool include_general_info) const
{
  if(!os.good())
  {
    Output::print("Error in ParameterDatabase::write. stream not good");
    return;
  }
  
  if(include_general_info) // write some information only once
  {
    os << "# current date and time: " << utilities::get_date_and_time();
    os << "# ParMooN git revision: " << parmoon::git_revision << "\n";
    os << "# ParMooN git branch  : " << parmoon::git_branch << "\n";
    os << "# ParMooN, local changes: " << std::boolalpha 
       << parmoon::git_local_changes << "\n";
    os << "# ParMooN build information: " << parmoon::build_type << ", " 
       << parmoon::parallel_type << "\n";
#ifdef _MPI
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    os << "# number of ParMooN mpi processes: " << size << "\n";
#endif
    os << "# hostname: " << utilities::get_host_name() << "\n";
  }
  os << "# A ParMooN parameter database with "
     << this->get_n_parameters() << " parameters and "
     << this->get_n_nested_databases() << " nested databases\n\n";
  if(include_general_info)
  {
    os << "# The name of the database.\n";
  }
  os << "[ " << this->name << " ]\n\n";
  
  for(const auto& p : this->parameters)
  {
    if(verbose)
    {
      // print the description as well as the parameter with its value and range
      p.print_description(os, "## ", 60, "");
      os << p.get_name() << ": " << p << "   " 
         << p.range_as_string() << "\n\n";
    }
    else
    {
      os << p.get_name() << ": " << p << "\n";
    }
  }
  auto n_nested = this->databases.size();
  if(n_nested != 0)
  {
    os << "\n# Nested databases:\n\n";
    unsigned int n = 1; // count the nested databases
    for(const auto& db : this->databases)
    {
      os << "# A nested parameter database of '" << this->name << "'. Number "
         << n << " of " << n_nested << "\n";
      db.write(os, verbose, false);
      n++;
    }
  }
}

/* ************************************************************************** */
void ParameterDatabase::write(const std::string& filename, bool verbose,
                              bool include_general_info) const
{
  std::ofstream os(filename);
  this->write(os, verbose, include_general_info);
  os.close();
}

/* ************************************************************************** */
// helper function to remove whitespace at the beginning of a string
void remove_leading_whitespace(std::string& s)
{
  // position of first non white space character
  auto pos = s.find_first_not_of(" \t");
  if(pos != std::string::npos)
    s = s.substr(pos);
  else 
    s = std::string(); // empty string
}
// helper function to remove whitespace at the end of a string
void remove_trailing_whitespace(std::string& s)
{
  // position of first non white space character
  auto pos = s.find_last_not_of(" \t");
  if(pos != std::string::npos)
    s = s.substr(0, pos+1);
  else 
    s = std::string(); // empty string
}
// helper function to create a list of all lines for one database from a stream
// additionally get the name of the database
std::list<std::string> get_lines_of_database(std::istream& is,
                                             std::string& name)
{
  std::list<std::string> lines_read; // read all lines into local variable
  std::string line; // read each line into this string
  bool found_name = false; // indicate if a database name has been found
  auto position = is.tellg(); // remember position of stream before 'getline'
  // loop over all lines until a new database is found. This is indicated by a
  // pair of brackets surrounding the name of the database 
  // '[_name_of_database_]'
  while(std::getline(is, line))
  {
    if(!line.empty() && line.at(0) == '[')
    {
      // starting position of the new name in string 'line', is right after the
      // opening bracket. The name ends at the closing bracket
      auto end_position_for_name = line.find("]");
      if(end_position_for_name != std::string::npos) // found closing bracket
      {
        if(found_name) // we found a second database here
        {
          // reset the stream position such that the next line to be read is
          // 'line', so that the second database can be read as well
          is.seekg(position);
          break; // no more lines are read into the list
        }
        found_name = true;
        name = line.substr(1, end_position_for_name-1);
        remove_leading_whitespace(name);
        remove_trailing_whitespace(name);
      }
      // else // we do nothing, maybe give a warning?
    }
    else
      lines_read.push_back(line);
    position = is.tellg(); // update to new position
  }
  if(!is)
    Output::print<3>("read last line of input stream");
  
  if(!found_name)
  {
    // no parameter database name found
    name = "custom parameter database without a given name";
  }
  return lines_read;
}

// helper function to identify the type of a variable in `value_string`. Exactly
// one of the other parameters are updated accordingly. The value_string must 
// not be a vector-valued parameter. This is only for a single value (possibly 
// within a vector).
Parameter::types get_parameter_and_type(const std::string& value_string,
                                        bool & bool_ret,
                                        int & int_ret, size_t & size_t_ret, 
                                        double & double_ret,
                                        std::string & string_ret)
{
  if(value_string == "True" || value_string == "true")
  {
    bool_ret = true;
    return Parameter::types::_bool;
  }
  else if(value_string == "False" || value_string == "false")
  {
    bool_ret = false;
    return Parameter::types::_bool;
  }
  else if(value_string.find('.') != std::string::npos 
         || value_string.find(std::string("e-")) != std::string::npos)
  {
    try // this could be a double
    {
      double_ret = stod(value_string);
      return Parameter::types::_double;
    }
    catch(...)
    {
      // this could be a filename or a path
    }
  }
  else 
  {
    try
    {
      int value = stoi(value_string);
      // in the following, we check if input = scientific notation, e.g 1e5
      int exponent = 0 ;
      auto position = value_string.find(std::string("e"));
      if (position != std::string::npos)
      {
        exponent = stoi(value_string.substr(position + 1));
      }
      value = value * std::pow(10,exponent);
      if(value < 0) // this is an int not a size_t
      {
        int_ret = value;
        double_ret = (double)value;
        return Parameter::types::_int;
      }
      // this is an unsigned parameter
      size_t_ret = (size_t)value; // safe, because value >= 0
      // we set int_ret and double_ret as well because maybe this parameter 
      // turns out to be of those types really, then we need these numbers
      int_ret = value;
      double_ret = (double)value;
      return Parameter::types::_size_t;
    }
    catch(...)
    {}
  }
  
  // then it must be a string
  string_ret = value_string;
  return Parameter::types::_string;
}
// helper function to identify the type of a variable in `value_string`. 
// Depending on the returned type the respective other parameter is filled.
// The other parameters might be used as well, but should then be discarded.
Parameter::types get_parameter_and_type(const std::string& value_string,
                                        bool & bool_ret,
                                        int & int_ret, size_t & size_t_ret, 
                                        double & double_ret,
                                        std::string & string_ret,
                                        std::vector<bool> & bool_vec,
                                        std::vector<int> & int_vec,
                                        std::vector<size_t> & size_t_vec,
                                        std::vector<double> & double_vec,
                                        std::vector<std::string> & string_vec)
{
  // find out if this is a vector:
  std::vector<std::string> vector_entries;
  auto npos = std::string::npos;
  auto left = value_string.find("(");
  if(left != npos)
  {
    // possibly a vector
    auto right = value_string.find(")");
    if(right != npos)
    {
      // consider this to be a vector
      left += 1; // just after the "("
      auto middle = value_string.find(",", left);
      while(middle < right && middle != npos)
      {
        std::string value = value_string.substr(left, middle-left);
        remove_trailing_whitespace(value);
        remove_leading_whitespace(value);
        vector_entries.push_back(value);
        left = middle+1; // just after the comma
        middle = value_string.find(",", left); // next comma
      }
      // find last entry
      std::string value = value_string.substr(left, right-left);
      remove_trailing_whitespace(value);
      remove_leading_whitespace(value);
      if(value.empty())
        ErrThrow("empty vector, provide at least one value");
      vector_entries.push_back(value);
    }
    // else // this could be a filename or a path (or an error)
  }
  if(!vector_entries.empty())
  {
    // only to get the type
    auto type = get_parameter_and_type(vector_entries.at(0), bool_ret, int_ret,
                                       size_t_ret,  double_ret, string_ret);
    
    typedef Parameter::types tps; // just for shorter code
    for(const auto s : vector_entries)
    {
      auto type_local = get_parameter_and_type(s, bool_ret, int_ret, size_t_ret, 
                                               double_ret, string_ret);
      if(type != type_local)
      {
        if(type == tps::_bool || type_local == tps::_bool)
          ErrThrow("cannot read vector where entries have different types ",
                   type, ", ", type_local);
        if(type == tps::_string || type_local == tps::_string)
          ErrThrow("can not read vector with different types, ", type, " ",
                   type_local);
        if(type == tps::_int)
        {
          if(type_local == tps::_size_t)
            int_ret = (int)size_t_ret;
          if(type_local == tps::_double)
          {
            // switch type to double
            double_vec.clear();
            double_vec.insert(double_vec.begin(), int_vec.begin(), 
                              int_vec.end());
            type = tps::_double;
          }
        }
        else if(type == tps::_size_t)
        {
          if(type_local == tps::_int)
          {
            // switch type to int
            int_vec.clear();
            int_vec.insert(int_vec.begin(), size_t_vec.begin(), 
                           size_t_vec.end());
            type = tps::_int;
          }
          else if(type_local == tps::_double)
          {
            // switch type to double
            double_vec.clear();
            double_vec.insert(double_vec.begin(), size_t_vec.begin(),
                              size_t_vec.end());
            type = tps::_double;
          }
        }
        else if(type == tps::_double)
        {
          if(type_local == tps::_int)
            double_ret = (double)int_ret;
          if(type_local == tps::_size_t)
            double_ret = (double)size_t_ret;
        }
      }
      switch(type)
      {
        case tps::_bool:
          bool_vec.push_back(bool_ret);
          break;
        case tps::_int:
          int_vec.push_back(int_ret);
          break;
        case tps::_size_t:
          size_t_vec.push_back(size_t_ret);
          break;
        case tps::_double:
          double_vec.push_back(double_ret);
          break;
        case tps::_string:
          string_vec.push_back(string_ret);
          break;
        default:
          ErrThrow("unknown type.");
      }
    }
    switch(type)
    {
      case tps::_bool:
        return tps::_bool_vec;
        break;
      case tps::_int:
        return tps::_int_vec;
        break;
      case tps::_size_t:
        return tps::_size_t_vec;
        break;
      case tps::_double:
        return tps::_double_vec;
        break;
      case tps::_string:
        return tps::_string_vec;
        break;
      default:
        ErrThrow("unknown type.");
    }
  }
  // else // not a vector
  return get_parameter_and_type(value_string, bool_ret, int_ret, size_t_ret, 
                                double_ret, string_ret);
}

// if the parameter range consists of values of a different type than `type`,
// that type needs to be adjusted. For example if the type is size_t but the 
// range has negative values, then the type is set to int. Also if the range 
// has a double value in it and the type is an integer type, we change the type
// to double.
void adjust_type(Parameter::types& type,
                 const std::pair<bool,std::set<std::string>>& range_list)
{
  // we need the following values to call get_parameter_and_type
  bool bool_val;
  int int_val;
  size_t size_t_val;
  double double_val;
  std::string string_val;
  std::vector<bool> bool_vec;
  std::vector<int> int_vec;
  std::vector<size_t> size_t_vec;
  std::vector<double> double_vec;
  std::vector<std::string> string_vec;
  
  for(const std::string & s : range_list.second)
  {
    Parameter::types t = get_parameter_and_type(s, bool_val, int_val, 
                                                size_t_val, double_val, 
                                                string_val, bool_vec, int_vec,
                                                size_t_vec, double_vec, 
                                                string_vec);
    if(type == t)
      continue;
    typedef Parameter::types types; // just to have less code
    if( (type == types::_bool_vec && t == types::_bool)
     || (type == types::_int_vec && t == types::_int)
     || (type == types::_size_t_vec && t == types::_size_t)
     || (type == types::_double_vec && t == types::_double)
     || (type == types::_string_vec && t == types::_string)
    )
      continue;
    // type and t are different
    if(type == types::_bool || t == types::_bool)
      ErrThrow("value and range have different types, one is bool");
    if(type == types::_string || t == types::_string)
      ErrThrow("value and range have different types, one is string");
    if(type == types::_size_t && t == types::_int)
      type = types::_int;
    if( (type == types::_size_t || type == types::_int) && t == types::_double)
      type = types::_double;
    // in other cases we don't change the type 
  }
}

// the bool is true in case of an interval and false otherwise
std::pair<bool,std::set<std::string>> get_range_list(const std::string& range)
{
  std::pair<bool, std::set<std::string>> ret; // to be returned
  if(range.empty())
    return ret; // nothing to do here
  auto npos = std::string::npos;
  // expected format: [ min, max ] or { a, b, c, d, e, ... }
  auto begin = range.find("[");
  if(begin != npos) // interval range
  {
    begin += 1;
    auto middle = range.find(",", begin);
    auto end    = range.find("]", begin);
    if(middle == npos || end == npos)
    {
      ErrThrow("wrong range format for an interval.");
    }
    if(range.find(",", middle+1) < end)
    {
      // another comma between '[' and ']'
      Output::print("WARNING a parameter interval should only have two values. "
                    "Did you mean a set? Then use '{' and '}'.");
    }
    std::string min_string = range.substr(begin, middle-begin);
    std::string max_string = range.substr(middle+1, end-middle-1);
    remove_leading_whitespace(min_string);
    remove_leading_whitespace(max_string);
    remove_trailing_whitespace(min_string);
    remove_trailing_whitespace(max_string);
    ret.first = true;
    ret.second.insert(min_string);
    ret.second.insert(max_string);
  }
  else // set range
  {
    auto begin = range.find("{");
    if(begin == npos)
      return ret; // no meaningful range given
    begin += 1; // just after the "{"
    auto middle = range.find(",", begin); // find first comma
    auto end    = range.find("}", begin);
    if(end == npos)
    {
      ErrThrow("wrong range format for a set.");
    }
    ret.first = false;
    while(middle < end && middle != npos)
    {
      std::string value = range.substr(begin, middle-begin);
      remove_trailing_whitespace(value);
      remove_leading_whitespace(value);
      ret.second.insert(value);
      begin = middle+1; // just after the comma
      middle = range.find(",", begin); // next comma
    }
    // find last entry
    std::string value = range.substr(begin, end-begin);
    remove_trailing_whitespace(value);
    remove_leading_whitespace(value);
    ret.second.insert(value);
  }
  return ret;
}

// helper function to call set_range on the parameter with the right type
void add_range_to_parameter(
  Parameter& p, const std::pair<bool, std::set<std::string>>& range_list)
{
  if(range_list.second.empty())
    return; // nothing we can do here
  switch(p.get_type())
  {
    case Parameter::types::_bool:
    case Parameter::types::_bool_vec:
    {
      if(range_list.first)
        ErrThrow("unable to set an interval range for a boolean parameter");
      auto& rl = range_list.second;
      if(rl.count("true") == 1 || rl.count("True") == 1)
      {
        if(rl.count("false") == 1 || rl.count("False") == 1)
          p.set_range(std::set<bool>{true, false});
        else
          p.set_range(std::set<bool>{ true });
      }
      else if(rl.count("false") == 1 || rl.count("False") == 1)
        p.set_range(std::set<bool>{ false });
      // else not a valid range for a bool
      break;
    }
    case Parameter::types::_int:
    case Parameter::types::_size_t:
    case Parameter::types::_int_vec:
    case Parameter::types::_size_t_vec:
    {
      if(range_list.first)
      {
        // interval range
        std::string min_string = *range_list.second.begin();
        std::string max_string;
        if(range_list.second.size() == 1) // one-element range
          max_string = min_string;
        else
          max_string = *range_list.second.rbegin();
        try
        {
          long min = stol(min_string);
          long max = stol(max_string);
          if(min > max)
            std::swap(min, max);
          if(p.get_type() == Parameter::types::_size_t 
            || p.get_type() == Parameter::types::_size_t_vec)
          {
            if(min < 0 || max < 0)
              ErrThrow("reading negative value for size_t??");
            p.set_range((size_t)min, (size_t)max);
          }
          else
            p.set_range((int)min, (int)max);
        }
        catch(...)
        {
          ErrThrow("could not read interval range for int parameter");
        }
      }
      else
      {
        // set range (not necessarily interval)
        std::set<long> long_range;
        try
        {
          for(const std::string& s : range_list.second)
          {
            long_range.insert(stol(s));
          }
        }
        catch(...)
        {
          ErrThrow("could not read range for int parameter");
        }
        // convert to the right `set` and set_range
        if(p.get_type() == Parameter::types::_size_t
          || p.get_type() == Parameter::types::_size_t_vec)
        {
          auto is_negative = [](long i){ return i < 0; };
          if(std::any_of(long_range.begin(), long_range.end(), is_negative))
            ErrThrow("reading negative value for size_t??");
          std::set<size_t> size_t_range(long_range.begin(), long_range.end());
          p.set_range(size_t_range);
        }
        else
        {
          std::set<int> int_range(long_range.begin(), long_range.end());
          p.set_range(int_range);
        }
      }
      break;
    }
    case Parameter::types::_double:
    case Parameter::types::_double_vec:
    {
      if(!range_list.first)
        ErrThrow("unable to set a set range for a double parameter");
      std::string min_string = *range_list.second.begin();
      std::string max_string;
        if(range_list.second.size() == 1)
          max_string = min_string;
        else
          max_string = *range_list.second.rbegin();
      try
      {
        double min = stod(min_string);
        double max = stod(max_string);
        if(min > max)
          std::swap(min, max);
        p.set_range(min, max);
      }
      catch(...)
      {
        ErrThrow("could not read range for double parameter ", min_string,
                 ",", max_string);
      }
      break;
    }
    case Parameter::types::_string:
    case Parameter::types::_string_vec:
    {
      if(range_list.first)
        ErrThrow("unable to set an interval range for a string parameter");
      p.set_range(range_list.second);
      break;
    }
    default:
      ErrThrow("unknown parameter type");
      break;
  }
}

// helper function to create a new parameter where the type is to be determined
Parameter create_parameter(const std::string& name,
                           const std::string& value_string,
                           const std::string& description,
                           const std::string& range_string)
{
  bool bool_val;
  int int_val;
  size_t size_t_val;
  double double_val;
  std::string string_val;
  std::vector<bool> bool_vec;
  std::vector<int> int_vec;
  std::vector<size_t> size_t_vec;
  std::vector<double> double_vec;
  std::vector<std::string> string_vec;
  
  
  
  // p will be returned at the end. We use a pointer because the Parameter class
  // has no default constructor.
  std::unique_ptr<Parameter> p;
  // find out what type this is
  Parameter::types type = get_parameter_and_type(value_string, bool_val,
                                                 int_val, size_t_val,
                                                 double_val, string_val,
                                                 bool_vec, int_vec,
                                                 size_t_vec, double_vec,
                                                 string_vec);
  
  auto range_list = get_range_list(range_string);
  
  adjust_type(type, range_list);
  
  // create a new parameter according to `type`
  switch(type)
  {
    case Parameter::types::_bool:
      p.reset(new Parameter(name, bool_val, description));
      break;
    case Parameter::types::_int:
      p.reset(new Parameter(name, int_val, description));
      break;
    case Parameter::types::_size_t:
      p.reset(new Parameter(name, size_t_val, description));
      break;
    case Parameter::types::_double:
      p.reset(new Parameter(name, double_val, description));
      break;
    case Parameter::types::_string:
      p.reset(new Parameter(name, string_val, description));
      break;
    case Parameter::types::_bool_vec:
      p.reset(new Parameter(name, bool_vec, description));
      break;
    case Parameter::types::_int_vec:
      p.reset(new Parameter(name, int_vec, description));
      break;
    case Parameter::types::_size_t_vec:
      p.reset(new Parameter(name, size_t_vec, description));
      break;
    case Parameter::types::_double_vec:
      p.reset(new Parameter(name, double_vec, description));
      break;
    case Parameter::types::_string_vec:
      p.reset(new Parameter(name, string_vec, description));
      break;
    default:
      ErrThrow("unknown type ");
      break;
  }
  add_range_to_parameter(*p, range_list);
  return std::move(*p);
}

// helper function to read a single parameter from one line and return strings 
// for name, value and range
bool read_parameter(const std::string& line, std::string& name, 
                    std::string& value, std::string& range_string)
{
  auto npos = std::string::npos;
  // find ':'. Whatever is before that is the parameter name
  auto position_for_name = line.find(":");
  if(position_for_name == npos)
    return false;
  
  name = line.substr(0, position_for_name);
  remove_leading_whitespace(name);
  remove_trailing_whitespace(name);
  
  // whatever is after ':' is the value
  position_for_name += 1;
  // what remains is the value and possibly the range
  std::string remainder = line.substr(position_for_name);
  // remove leading whitespace
  remove_leading_whitespace(remainder);
  
  if(remainder.size() == 0)
    return false; // no value given, maybe write a warning
  
  auto position_for_value = remainder.find_first_of("{["); // may be npos
  value = remainder.substr(0, position_for_value);
  remove_trailing_whitespace(value);
  
  range_string.clear();
  if(position_for_value != npos)
  {
    // the now remaining part is the range which is not always given
    range_string = remainder.substr(position_for_value);
    remove_leading_whitespace(range_string);
    remove_trailing_whitespace(range_string);
  }
  return true;
}

/* ************************************************************************** */
void ParameterDatabase::read(std::istream& is, bool root)
{
  if(!is.good())
  {
    Output::print("Error in ParameterDatabase::read. stream not good");
    return;
  }
  
  ParameterDatabase tmp("");
  
  Output::print<2>("\nReading database from stream");
  
  // get all lines of the input stream as a list of strings
  std::list<std::string> lines_read = get_lines_of_database(is, tmp.name);
  
  
  std::string description; // accumulates from several lines
  // loop over all lines
  for(std::string & line : lines_read)
  {
    remove_leading_whitespace(line);
    // check if this line is empty
    if(line.length() == 0)
    {
      // no parameter on this line, reset the description
      description.clear();
      continue; // go to next line
    }
    if(line.at(0) == '#')
    {
      // documentation for a parameter on this line
      std::string des = line.substr(1);
      remove_leading_whitespace(des);
      // if there were two '#' at the beginning, then interpret this as a 
      // comment which should not become part of the description
      if(des.length() != 0 && des.at(0) == '#')
        description.clear();
      else
        // add this line to the description
        description += des + " ";
      continue; // go to next line
    }
    
    // check for a parameter
    std::string param_name, value_string, range_string;
    bool found_parameter = read_parameter(line, param_name, value_string,
                                          range_string);
    Output::print<4>("Reading parameter from stream: ", param_name, " \t ",
                     value_string, " \t ->", range_string, "<-     ",
                     found_parameter);
    if(!found_parameter)
    {
      description.clear();
      continue;
    }
    
    tmp.add(create_parameter(param_name, value_string, description,
                             range_string));
    description.clear();
  }
  Output::print<2>("Done reading database '", tmp.get_name(),
                   "' from stream. Read ", tmp.get_n_parameters(),
                   " parameters.");
  
  this->name = tmp.get_name();
  this->merge(tmp);
  // in case there are more databases in the stream, read those as well, but
  // only in the root database
  while(is.good() && root)
  {
    // additional database found
    ParameterDatabase nested_db("Nested parameter Database");
    nested_db.read(is, false);
    this->add_nested_database(nested_db);
  }
}

/* ************************************************************************** */
void ParameterDatabase::read(const std::string& filename)
{
  std::ifstream ifs(filename);
  this->read(ifs);
  ifs.close();
}

/* ************************************************************************** */
void ParameterDatabase::read(std::istream& is)
{
  this->read(is, true);
}

/* ************************************************************************** */
void ParameterDatabase::merge(const ParameterDatabase &other,
                              bool create_new_parameters, bool recursive)
{
  for(const Parameter & p : other.parameters)
  {
    if(this->contains(p.get_name()))
    {
      this->operator[](p.get_name()).impose(p);
    }
    else if(create_new_parameters)
      this->add(Parameter(p)); // add a copy of the parameter p
  }
  if(recursive)
  {
    for(const ParameterDatabase& db : other.databases)
    {
      auto it = find_in_list(db.get_name(), this->databases);
      if(it == this->databases.end())
      {
        this->databases.emplace_back(db);
      }
      else
      {
        it->merge(db, create_new_parameters);
      }
    }
  }
}

/* ************************************************************************** */
void ParameterDatabase::info(bool only_names) const
{
  bool do_print = true;
#ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  do_print = (rank == 0);
#endif
  if(do_print)
  {
    Output::print("Parameter database: ", this->name);
    Output::print("  number of parameters: ", this->parameters.size());
    auto n_nested = this->databases.size();
    if(n_nested != 0)
    {
      Output::print("  number of nested parameter databases: ", n_nested);
    }
    for(const auto& p : this->parameters)
    {
      if(only_names)
        Output::print("    ", p.get_name(), ": ", p.value_as_string());
      else
        p.info();
    }
    unsigned int n = 1; // count the nested databases
    for(const auto& db : this->databases)
    {
      Output::print("Info on a nested parameter database of '", this->name,
                    "'. Number ", n, " of ", n_nested);
      db.info(only_names);
      n++;
    }
  }
}

/* ************************************************************************** */
// set default parameters
ParameterDatabase ParameterDatabase::parmoon_default_database()
{
  ParameterDatabase db("default ParMooN parameter database");

  // add parameters which are needed by all ParMooN programs and don't belong
  // anywhere else.

  db.add("outfile", "default_parmoon_outfile.out",
         "This is the file where all output of ParMooN is (usually) written "
         "to. In general ParMooN produces text output on the console as well "
         "as in this file. For this to properly work, you should call "
         "`Output::set_outfile(db[\"outfile\"]);` in your main program.");

  db.add("problem_type", 0u, 
         "Determine which kind of problem you want to solve. A value of 0 "
         "means not set. Other values have the following meanings: "
         "1: stationary convection-diffusion,  2: time-dependent "
         "convection-diffusion,  3: stationary Stokes,  4: time-dependent "
         "Stokes,  5: stationary Navier-Stokes,  6: time-dependent "
         "Navier-Stokes 7: Brinkman.",
         0u, 7u);

  db.add("output_write_ps", false,
         "Draw a postscript file of the domain. This only works in two space "
         "dimensions. Usually this is used in the main program.",
         {true,false});

  db.add("verbosity", 1u,
         "Set the verbosity of ParMooN. The higher the number, the more will "
         "output you will get. Such output will be written to console and the "
         "'outfile'.", 1u, 5u);

  db.add("script_mode", false, "Set ParMooN into script mode. This means all "
         "output is written to the outfile and not to console.");

  db.add("write_snaps", false,
         "Write Snapshots. If set to true, snapshots of the numerical solution "
         "will be written into a file.",
         {true,false});

  db.add("compute_POD_basis", false,
         "Compute the POD basis from snapshots. Snapshots have to be already "
         "available in 'db[snaps_directory]/db[snaps_basename]', or can be "
         "acquire by activating the parameter 'db[write_snaps]'.",
         {true,false});

  db.add("ROM_method", false,
         "Solve the PDE using Reduced Order Modelling. The POD basis has "
         "to be already available in 'db[pod_directory]/db[pod_basename].pod' "
         "or can be computed by activating the parameter "
         "'db[compute_POD_basis]'.",
         {true,false});

  return db;
}

/* ************************************************************************** */
ParameterDatabase ParameterDatabase::default_nonlinit_database()
{
  ParameterDatabase db("default ParMooN nonlinear iteration parameters database");

  //TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE, TDatabase::ParamDB->SC_NONLIN_MAXIT_SCALAR;
  db.add("nonlinloop_maxit", 100u,
         "The maximum number of iterations to perform in a non-linear loop.",
         0u, 1000u);

  db.add("nonlinloop_minit", 1u,
         "The minimum number of iterations to perform in a non-linear loop.",
         0u, 1000u);

  // TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE, TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALAR,
  db.add("nonlinloop_epsilon", 1e-10,
         "At which absolute residual 2-norm to break the nonlinear loop.",
         0. , 1. );

  db.add("nonlinloop_max_epsilon", 1e-10,
         "At which absolute residual infty-norm to break the nonlinear loop.",
         0. , 1. );

  db.add("nonlinloop_hard_epsilon", false,
         "Should the iteration be considered a failure if nonlinloop_epsilon "
         "is not reached?",
         {true,false});

  db.add("nonlinloop_damping_factor", 1.0,
         "Damping factor 'w' for the nonlinear iteration. The solution of the "
         "k-th iterate will be scaled by 'w'. Then The previous solution, "
         "scaled by '1-w', will be added. Setting to it to zero makes no "
         "sense.",
         0., 1.);

  db.add("nonlinloop_damping_auto", false,
         "If true, the best damping factor for the nonlinear iteration will "
         "be selected automatically. Currently only implemented for TNSE.",
         {true,false});

  db.add("nonlinloop_damping_auto_precision", 10,
         "If auto damping is on, sets the selection's precision in the fully "
         "nonlinear case. Note that in the fully nonlinear case, large "
         "values will be expensive due to repeated reassembly.",
         {1, 5, 10, 20, 50, 100});

  db.add("nonlinloop_damping_auto_binary_search", false,
         "If auto-damping is on, increases precision with a binary search.",
         {true,false});

  db.add("nonlinloop_damping_auto_skip", 2,
         "Number of nonlinear iterations to skip before applying auto-damping.",
         0, 100);

  db.add("nonlinloop_diagnostics", false,
         "If true, diagnostic information will be written to a CSV file after "
         "each step of the nonlinear loop.",
         {true,false});

  //TDatabase::ParamDB->SC_NONLIN_DIV_FACTOR
  db.add("nonlinloop_slowfactor", 1.e10,
         "Determines at which reduction rate over x iterations"
         "(usually x = 10, see system classes) a convergence is interpreted"
         "as too slow and therefore the iteration is stopped.",
         0., 1.e10);

  //TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALE_SADDLE, TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALE_SCALAR
  db.add("nonlinloop_scale_epsilon_with_size", false,
         "Whether or not to scale the absolute residual breaking criterion"
         "with the square root of the problem size.",
         {true,false});

  db.add("nonlinloop_residual_relative_to_rhs", false,
         "If true, not the pure residual is compared to nonlinloop_epsilon ,"
         "but the residual is divided by the norm of the initial right hand "
         "side vector. Note that this is not yet enabled in each problem class, "
         "you must check that, if you plan to use it.", {true, false});

  return db;
}

ParameterDatabase ParameterDatabase::default_output_database()
{
  ParameterDatabase db("default ParMooN output control parameters database");

  db.add("output_write_vtk", false,
         "This parameter can control, whether an output method"
         "of a system class will produce VTK output or not.",
         {true,false});

  db.add("output_write_xdmf", false,
         "This parameter can control, whether an output method "
         "of a system class will produce xdmf output or not.",
         {true,false});

  db.add("output_xdmf_format", 2u,
         "This parameter controls the output format of the heavy data for "
         "the xdmf file output. 0 means plaintext, 1 means binary, and 2 "
         "means hdf5.",
         0u, 2u);

  db.add("output_xdmf_timetype", 1u,
         "This parameter controls the time type of the xdmf file output. "
         "0 means list of discrete values, 1 means a collection of values "
         "with (Start, Stride, Count).",
         0u, 1u);

  db.add("separate_h5_files", true,
         "This parameter can control, whether the xdmf output method stores "
         "the solution data in separate h5 files.",
         {true,false});
  
  db.add("output_compress_h5_files", true,
         "This parameter controls whether the h5 files written by the "
         "DataWriter class should be compressed by zlib or not.",
         {true,false});

  db.add("output_collective_h5_files", true,
         "Controls whether the h5 files written by the DataWriter "
         "class should be written using collective I/O or collated "
         "and written on the root process only.",
         {true,false});

  db.add("output_write_csv", false,
         "Controls whether the DataWriter class will write an additional "
         "CSV file tabulating the global output variables.",
         {true,false});

  db.add("output_xdmf_globals", true,
         "Controls whether the xdmf files written by the DataWriter "
         "class should include global output variables.",
         {true,false});

  db.add("output_write_subdomain_info", false,
         "This parameter controls whether the output should contain the "
         "subdomains for each processor in mpi computations. This is often "
         "useful for debugging.",
         {true,false});
  
  db.add("output_compute_errors", true,
         "Do or do not compute errors after computing a solution. This makes "
         "much sense if an analytical solution is known. If not then it is "
         "often simply set to zero and computing errors then means computing "
         "norms, e.g. the L^2-norm of the solution.",
         {true,false});
  
  db.add("output_compute_minmax", true,
         "Do or do not compute the minimum and maximum of a solution. Note "
         "that this computation is rather time consuming.");

  db.add("output_directory", ".",
         "This directory is where the output is written. This "
         "directory will be created, if it does not exist already. Files in "
         "this directory will be overwritten without any warning.");

  db.add("output_basename", "parmoon",
         "This string is prepended to most files written by ParMooN. "
         "This includes also vtk- and case-files");

  db.add("steps_per_output", (size_t)1,
         "This integer specifies how many (time) steps are performed "
         "before writing the results ");

  db.add("continue_output_after_restart", false,
         "This parameter, when true, allows the output numbering "
         "to continue after a simulation restart, which is often a "
         "desirable behavior. Otherwise, it behaves as if the "
         "initial time was 0, and starts the numbering at 0, leading "
         "sometimes to overwritting existing files.", 
         {true,false});

  db.add("output_compute_time_average", false,
         "Do or do not compute time average of the solution.",
         {true,false});

  db.add("output_along_line", false,
         "Do or do not write the solution along lines, "
         "according to the lines defined in the nested database used "
         "by the class LinesEval, see EvalTools.h",
         {true,false});

  db.add("start_time_averaging_at", 0.,
         "Time at which the time averaging will start.",
         0., 1.0e10);

  db.add("output_write_exact_solution", false,
         "If set to true, this parameter allows to write the exact solution "
         "into the vtk or case output file. "
         "Note: the exact solution must be specified in the example file",
         {true,false});

  db.add("compute_vorticity_divergence", false,
         "To Compute vorticity and Divergence, vorticity is used to compute"
         " the thickness, Note: the routine is implemented in the example file",
         {true,false});

  db.add("compute_wall_shear_stress", false,
         "Compute the solution's wall shear stress?",
         {true,false});

  db.add("wall_shear_stress_mode", "face",
         "Compute the solution's wall shear stress by face or by vertex?",
         {"face", "vertex"});

  db.add("compute_turbulent_kinetic_energy", false,
         "Compute the solution's turbulent kinetic energy? "
         "This requires an LES discretization (\"smagorinsky\") and uses the "
         "approximation \\nu_t ~= 2 C_e \\delta \\sqrt k, with a constant C_e "
         "set using the parameter \"turbulent_kinetic_energy_scale\".",
         {true,false});

  db.add("turbulent_kinetic_energy_scale", 1.0,
         "Constant C_e to use in computing the turbulent kinetic energy.",
         0.0, 1.e+12);

  db.add("turbulent_kinetic_energy_mode", "cell",
         "Compute the solution's wall shear stress by cell or by vertex?",
         {"cell", "vertex"});

  db.add("compute_effective_viscosity", false,
         "Compute the solution's effective viscosity, in the "
         "non-Newtonian case?",
         {true,false});

  db.add("compute_time_derivative", false,
         "Compute the solution's time derivative?",
         {true,false});

  return db;
}

