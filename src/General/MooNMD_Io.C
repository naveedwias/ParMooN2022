#include <MooNMD_Io.h>

namespace Output
{
  // the verbosity threshold is usually between 1 (low verbosity) and
  // maxVerbosity(high verbosity). A value of zero is possible and means no 
  // output at all. This value is only possible if you called
  // Output::suppressAll().
  unsigned int verbosityThreshold = 1;
  
  void increaseVerbosity(unsigned int i)
  {
    verbosityThreshold = std::min(verbosityThreshold + i, maxVerbosity);
  }
  
  void decreaseVerbosity(unsigned int i) // increase the verbosityThreshold
  {
    if(i >= verbosityThreshold)
    {
      verbosityThreshold = 1;
    }
    else
    {
      verbosityThreshold -= i;
    }
  }
  
  void setVerbosity(unsigned int i)
  {
    verbosityThreshold = std::max(1u, std::min(5u,i));
  }
  
  unsigned int getVerbosity()
  {
    return verbosityThreshold;
  }
  
  void suppressAll()
  {
    verbosityThreshold = 0;
  }
  
  // 'verbosity' is greater than 0 (not equal)
  bool verbose(unsigned int verbosity)
  {
    return verbosity <= verbosityThreshold;
  }
  
  
  // the outfile into which all output is written
  std::ofstream outfile;
  // another outfile to be used to write some local output. This enables to
  // write a file with the contents of for example all output of a function
  // 
  // Output::redirect("filename.txt");
  // Output::increaseVerbosity(5); // set some high verbosity
  // f(); // call function which uses print a lot
  // Output::decreaseVerbosity(5); // set some low verbosity
  // Output::resetOutfile(); // set to original outfile
  //
  std::ofstream local_outfile;
  // the streambuf for std::cout. it is stored here so that after calling
  // redirect (any calls to std::cout will be written to some file), one can
  // restore the original behavior.
  std::streambuf * coutBuffer = nullptr;
  // no output to console, respect verbosity settings
  bool script_mode = false;
  // store all warnings here, to print them again at the end of each program
  std::ostringstream collected_warnings;
  
  void set_outfile(const std::string& filename, bool sm)
  {
    if(outfile.is_open())
    {
      outfile.close();
    }
    outfile.open(filename, std::ofstream::out);
    outfile.setf(std::ios::scientific);
    script_mode = sm;
  }
  
  void redirect(const std::string& filename, bool overwrite_file)
  {
    if(local_outfile.is_open())
    {
      local_outfile.close();
    }
    (overwrite_file) ? local_outfile.open(filename, std::ofstream::out) :
      local_outfile.open(filename, std::ofstream::app);
    // store cout buffer to be able to restore later again
    coutBuffer = std::cout.rdbuf();
    std::cout.rdbuf(local_outfile.rdbuf());
  }
  
  void resetOutfile()
  {
    if(local_outfile.is_open())
    {
      local_outfile.close();
    }
    if(coutBuffer != nullptr)
    {
      std::cout.rdbuf(coutBuffer);
    }
  }
  
  std::ofstream& get_outfile()
  {
    if(local_outfile.is_open())
    {
      return local_outfile;
    }
    else
    {
      return outfile;
    }
  }
  
  bool writeOnlyToFile()
  {
    return local_outfile.is_open();
  }
  
  void set_script_mode(bool sm)
  {
    script_mode = sm;
  }
  
  bool in_script_mode()
  {
    return script_mode;
  }
  
  std::ostringstream& get_warnings()
  {
    return collected_warnings;
  }
  
  void print_warnings()
  {
    std::string w = collected_warnings.str();
    if(!w.empty())
    {
      // note that warnings were stored starting with a new line ("\n")
      print<1>("##### Summary of issued warnings #####", w);
    }
  }

  
  void close_file()
  {
    if(outfile.is_open())
    {
      outfile.close();
    }
    if(local_outfile.is_open())
    {
      local_outfile.close();
    }
  }
}
