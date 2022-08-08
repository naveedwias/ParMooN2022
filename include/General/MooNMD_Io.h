// =======================================================================
// 
// Purpose:     Collection of routine and methods for IO
//
// Author:      Gunar Matthies
//
// Date:        2004/02/19
// 
// =======================================================================

#ifndef __MOONMD_IO__
#define __MOONMD_IO__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdexcept>

#ifdef _MPI
#include <mpi.h>
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::ends;
using std::setw;
using std::setprecision;


/**
@brief everything related to writing to console and files
 
Especially the methods print and printToFile are defined in this namespace. The
proposed usage within ParMooN is as follows:

At the beginning of any ParMooN program make a call to 
`parmoon::parmoon_initialize()` will call  Output::set_outfile().
This will make sure that all calls to print will write its arguments to
`std::cout` as well as to the outfile. Before calling 
Output::set_outfile the print methods only write to `std::cout`.

Use print with any number of arguments. An additional template parameter 
(`unsigned int`) is recommended to indicate at which verbosity level you 
really want to print. Higher numbers mean smaller likelihood of being 
printed. There is a threshold value stored in this namespace 
which determines whether or not to print and is in the range of 1 and 5 
(maxVerbosity). For example

    print<2>("print only if verbosity level is 2 or larger")
    print<3>("print only if verbosity level is 3 or larger")
    print<4>("print only if verbosity level is 4 or larger")

Calling Output::print without a template argument, is equivalent to calling 
Output::print<1>. This means it will be printed for any verbosity level.
There is only one exception: If Output::suppressAll() is called, then nothing
is ever written to `std::cout` or the outfile. However in that case it is still
possible to write to a second outfile through calling redirect.


In summary the Output::print method writes to

- only `std::cout` if Output::set_outfile has not yet been called (discouraged)
- `std::cout` and the outfile if Output::set_outfile has been called
  (recommended)
- only to an extra outfile if Output::redirect has been called
- nowhere if Output::suppressAll has been called

where only in the first two cases the verbosity matters at all.

You can control the verbosity level using the methods 
Output::setVerbosity(unsigned int i),
Output::increaseVerbosity() and Output::decreaseVerbosity().

At the end of the program parmoon::parmoon_finalize() calls 
Output::close_file().

The user should only call the following methods from this namespace and ignore
the rest of it:

- Output::increaseVerbosity
- Output::decreaseVerbosity
- Output::setVerbosity
- Output::print<>
- Output::printToFile<>
- Output::set_outfile
- Output::close_file
- Output::suppressAll
- Output::redirect
- Output::resetOutfile
- Output::set_script_mode
- Output::warn, Output::stat, Output::info, Output::root_warn,
  Output::root_info, Output::dash
- Output::print_warnings
- #ErrThrow
 */
namespace Output
{
  /// @brief increase the verbosity by \p i
  void increaseVerbosity(unsigned int i = 1);
  /// @brief decrease the verbosity by \p i
  void decreaseVerbosity(unsigned int i = 1);
  /// @brief set the verbosity to a given level (must be between 1 and 5)
  void setVerbosity(unsigned int i);
  /// @brief get the verbosity, typically to reset it to this value later again
  unsigned int getVerbosity();
  /// @brief suppress all output to std::cout and the outfile
  ///
  /// Only calls to Output::print are suppressed. Direct calls to std::cout are
  /// still possible. This lasts until Output::increaseVerbosity, 
  /// Output::decreaseVerbosity, or Output::setVerbosity is called again.
  void suppressAll();
  
  /// @class Temporary_verbosity
  /// @brief create an object of this class to temporarily set the verbosity to
  /// a specific value.
  /// The destructor resets the verbosity to the previous value.
  /// @warning do not create multiple objects of this type in the same scope.
  class Temporary_verbosity;
  
  /// @brief set the outfile
  ///
  /// The print methods print to both std::cout and the outfile. This method 
  /// opens the file and it should be closed again calling Output::close_file().
  void set_outfile(const std::string& filename, bool script_mode = false);
  
  /// @brief close files which have been previously opened
  ///
  /// Files can be opened through Output::set_outfile and Output::redirect.
  void close_file();
  
  /// @brief print all warnings which were issued until now through Output::warn
  ///
  /// This makes most sense at the end of a program and should help not missing
  /// issued warnings in longer program output.
  void print_warnings();
  
  /** \brief redirect all output to this file, nowhere else
   *
   * This method opens the file and all subsequent calls to print and std::cout
   * will write into this file. To return to the normal behavior of the print
   * methods call Output::resetOutfile().
   * \note The default behaviour overwrites the file if it already exists. If
   * you want to append information instead of overwriting the file, set \p
   * overwrite_file to false.
   *
   * Typically this is used to write all
   * output of a block of code into a separate file:
   * \code{.cpp}
   *      Output::redirect("file_for_output_of_f");
   *      f(); // call some function which calls print a lot
   *      Output::resetOutfile();
   * \endcode
   * \note Using `printf` will still give you output on the console
   *
   * \note As it's the usual behaviour of print, the string is only written to
   * the file if the verbosity is high enough. So you might want to increase the
   * verbosity before using redirect and decrease it after
   * Output::resetOutfile().
   */
  void redirect(const std::string& filename, bool overwrite_file = true);
  /// @brief close the file after a call to Output::redirect.
  ///
  /// This restores the normal behavior of the print methods.
  void resetOutfile();
  
  /// @brief do not write anything to console, only to the outfile
  ///
  /// In contrast to the Output::redirect method, this respects the verbosity 
  /// setting. If you do use Output::redirect the script mode does not play a
  /// role.
  void set_script_mode(bool script_mode);
  
  /// @brief Write something to std::cout and the outfile 
  ///
  /// This is printed only if the template parameter verbosity is not too large.
  /// That means verbosity 1 (the default) will always be printed, while with 
  /// larger numbers it is only printed if the verbosity level is high enough.
  /// Note that verbosity must be within 1 and 5 (maxVerbosity), other values
  /// (including 0) give a compile error. 
  /// 
  /// The template parameters `Arguments` can be of any type which can be 
  /// inserted into a std::ostream via `<<`.
  ///
  /// This is a variadic template to enable variable number of arguments to 
  /// print.
  template<unsigned int verbosity = 1, typename ... Arguments>
  void print(Arguments const& ... rest);
  
  template<unsigned int verbosity = 1, typename ... Arguments>
  void warn(const std::string& context, Arguments const& ... rest);

  template<unsigned int verbosity = 1, typename ... Arguments>
  void stat(const std::string& context, Arguments const& ... rest);

  template<unsigned int verbosity = 1, typename ... Arguments>
  void info(const std::string& context, Arguments const& ... rest);

  template<unsigned int verbosity = 1, typename ... Arguments>
  void root_info(const std::string& context, Arguments const& ... rest);

  template<unsigned int verbosity = 1, typename ... Arguments>
  void root_warn(const std::string& context, Arguments const& ... rest);

  template<unsigned int verbosity = 1, typename ... Arguments>
  void dash(Arguments const& ... rest);

  /// @brief print only to the outfile depending on verbosity
  ///
  /// This method behaves like the corresponding print, but does not print to
  /// std::cout, only to the outfile.
  ///
  /// This is a variadic template to enable variable number of arguments to 
  /// print.
  template<unsigned int verbosity = 1, typename ... Arguments>
  void printToFile(Arguments const& ... rest);
  
  /// @brief write an error message to the outfile and throw an exception
  ///
  /// This is the preferred way to terminate the program if necessary. The 
  /// exception thrown is of type `std::runtime_error`. For convenience there 
  /// is a macro #ErrThrow defined which calls this method with the correct 
  /// \p file and \p line.
  template<typename ... Arguments>
  [[ noreturn ]] void errThrow(std::string file, int line,
                               Arguments const& ... args);
  
  /// @brief macro to call Output::errThrow with the correct file and line
  ///
  /// You have to provide a meaningful error message here.
  #define ErrThrow(...) {Output::errThrow(__FILE__, __LINE__, __VA_ARGS__);}
}


// this is an attempt to separate declaration and implementation. The template
// methods unfortunately can not be put in the source file, so they are 
// implemented here together with some helping methods. But in general you only
// have to worry about the above declarations.
namespace Output
{
  /// @brief maximum allowed verbosity level
  constexpr unsigned int maxVerbosity = 5;
  /// @brief a helper function the user should not call
  ///
  /// Return true if is v smaller or equal to the current verbosity level.
  /// v should be between 1 and maxVerbosity, but this is not enforced here.
  bool verbose(unsigned int v);
  
  class Temporary_verbosity
  {
    public:
      explicit Temporary_verbosity(unsigned int new_verbosity)
      {
        previous_verbosity = Output::getVerbosity();
        if(new_verbosity != 0)
          Output::setVerbosity(new_verbosity);
      }
      ~Temporary_verbosity()
      {
        Output::setVerbosity(previous_verbosity);
      }
    private:
      unsigned int previous_verbosity;
  };
  
  /// @brief a helper function the user should not call
  ///
  /// Returns true after a call to Output::redirect, otherwise false.
  bool writeOnlyToFile();
  /// @brief a helper function the user should not call
  /// 
  /// Get the stream of the outfile. If Output::redirect has been called this 
  /// returns the stream of that file.
  std::ofstream& get_outfile();
  /// @brief a helper function the user should not call
  /// 
  /// Get the stream which stores all the warnings generated in the program.
  std::ostringstream& get_warnings();
  
  /// @brief determine if scrip mode is turned on.
  bool in_script_mode();
  
  template<typename ... Arguments>
  std::ostringstream combine(Arguments const& ... args)
  {
    // the following writes all the arguments into the stream. With c++17 we can
    // also use 'fold expressions':
    // (std::cout << ... << args) << "\n";
    // (get_outfile() << ... << args) << "\n";
    // http://stackoverflow.com/a/21812549
    // http://stackoverflow.com/a/33621132
    std::ostringstream stream;
#ifdef _MPI
    // In MPI case, tell us which process is talking to us.
    // This feature is in testing stage.
    int finalized = 0;
    MPI_Finalized(&finalized);
    int initialized = 0;
    MPI_Initialized(&initialized);
    if(!finalized && initialized)
    {
      int my_rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
      stream << my_rank << ": ";
    }
    else
    {
      stream << "(?): ";
    }
#endif
    using List= int[];
    (void)List{0, ( (void)(stream << args), 0 ) ... };
    return stream;
  }
  
  // implementation of the print method
  template<unsigned int verbosity, typename ... Arguments>
  void print(Arguments const& ... args)
  {
    static_assert(verbosity > 0,
        "It makes no sense to call print with 0 verbosity"
        "Use a value greater than 0");
    static_assert(verbosity <= maxVerbosity, 
        "calling Output::print with verbosity too large");

    if(!verbose(verbosity))
      return; // nothing to write, skip the rest

    std::ostringstream stream = combine(args ...);

    if(!writeOnlyToFile() && !in_script_mode())
    {
      std::cout << stream.str() << "\n";
    }
    get_outfile() << stream.str() << "\n";
  }
  
  template<unsigned int verbosity, typename ... Arguments>
  void warn(const std::string& context, Arguments const& ... rest)
  {
    // is it a good idea to check the verbosity here, or should a warning always
    // be written?
    if(!verbose(verbosity))
      return;
    print<verbosity>("WARNING (",context,"): ", rest ...);
    // (get_warnings() << ... << args) << "\n"; // c++17
    std::ostringstream stream = combine("WARNING (",context,"): ", rest ...);
    get_warnings() << "\n" << stream.str();
  }

  template<unsigned int verbosity, typename ... Arguments>
  void stat(const std::string& context, Arguments const& ... rest)
  {
	  print<verbosity>("STATS (",context,"): ", rest ...);
  }

  template<unsigned int verbosity, typename ... Arguments>
  void info(const std::string& context, Arguments const& ... rest)
  {
	  print<verbosity>("INFO (",context,"): ", rest ...);
  }

  /// Use this method whenever only the root process should give
  /// an info. In sequential, this performs just the same as info.
  template<unsigned int verbosity, typename ... Arguments>
  void root_info(const std::string& context, Arguments const& ... rest)
  {
    int my_rank = 0;

#ifdef _MPI
    int finalized = 0;
    MPI_Finalized(&finalized);

    int initialized = 0;
    MPI_Initialized(&initialized);

    if(!finalized && initialized)
    {
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    }
#endif

    if(my_rank == 0)
    {
      print<verbosity>("INFO (",context,"): ", rest ...);
    }
  }

  /// Use this method whenever only the root process should give
  /// an info. In sequential, this performs just the same as warn.
  template<unsigned int verbosity, typename ... Arguments>
  void root_warn(const std::string& context, Arguments const& ... rest)
  {
    int my_rank = 0;

#ifdef _MPI
    int finalized = 0;
    MPI_Finalized(&finalized);

    int initialized = 0;
    MPI_Initialized(&initialized);

    if(!finalized && initialized)
    {
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    }
#endif

    if(my_rank == 0)
    {
      warn<verbosity>(context, rest ...);
    }
  }

  template<unsigned int verbosity, typename ... Arguments>
  void dash(Arguments const& ... rest)
  {
	  print<verbosity>(" > ", rest ...);
  }

  // implementation of the printToFile method
  template<unsigned int verbosity, typename ... Arguments>
  void printToFile(Arguments const& ... args)
  {
    bool sm = in_script_mode();
    set_script_mode(true); // temporarily turn off output to console
    print<verbosity>(args ...);
    set_script_mode(sm);
  }
  
  // implementation of the errThrow method
  /// @todo adapt to MPI (write rank and call MPI_Finalize())
  template<typename ... Arguments>
  [[ noreturn ]] void errThrow(std::string file, int line,
                               Arguments const& ... args)
  {
    // see the implementation of Output::print
    std::ostringstream stream;
    stream << std::string(60, '*') << "\nError in file " << file << ", line " 
           << line << ":\n";
    using List= int[];
    (void)List{0, ( (void)(stream << args), 0 ) ... };
    
    get_outfile() << stream.str() << "\n";
    throw std::runtime_error(stream.str());
  }
}

#endif
