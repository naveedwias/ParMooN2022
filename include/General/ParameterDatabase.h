#ifndef __PARAMETERDATABASE__
#define __PARAMETERDATABASE__

#include <Parameter.h>
#include <list>
#include <iostream>


/** @brief store a number of Parameter objects
 * 
 * This class stores `Parameter` objects where no two Parameters have the same
 * name. One can add parameters via some `add` methods and access them via
 * `operator[]` given the name of the desired parameter as a string.
 * 
 * Furthermore this class can write all its contents to an output stream (e.g. a
 * file stream) and also read from an input stream, see methods `read` and 
 * `write` for a description of the format.
 * 
 * If one has two ParameterDatabase objects, one can merge the one into the 
 * other, see the method `merge` for details.
 * 
 * To see all stored Parameters call `ParameterDatabase::info()`.
 * 
 * Furthermore, it is possible to store additional parameter databases. This 
 * can be nested up to an arbitrary order.
 */
class ParameterDatabase
{
  public:
    /// @brief construct an empty parameter database with a given name
    explicit ParameterDatabase(std::string name);
    
    /// @brief construct a database filled with parameters of general interest
    ///
    /// These parameters include "outfile", "boundary_file", "geo_file", 
    /// "problem_type", "base_name", ...
    static ParameterDatabase parmoon_default_database();
    
    /// @brief construct a database filled with parameters which holds
    /// controls and stopping criteria for a nonlinear iteration loop
    static ParameterDatabase default_nonlinit_database();

    /// @brief construct a database which holds all those parameters
    /// that are typically needed in the output methods of the system
    /// classes
    /// TODO CB The more of these default databases I add here, the more I am
    /// convinced: they are a symptom for classes which we ought to have but don't.
    static ParameterDatabase default_output_database();
    
    /// A database to control the writing of snapshots into a file
    /// This feature is used to build a POD and a reduced-order model
    static ParameterDatabase get_default_snapshots_database();

    /// A database to control the computation, writing, and reading
    /// of POD basis
    static ParameterDatabase get_default_pod_database();

    /// @brief delete all parameters from this database
    ~ParameterDatabase() = default;
    
    /// @brief copy constructor, this is a deep copy
    ParameterDatabase(const ParameterDatabase&);
    
    /// @brief default move constructor
    ParameterDatabase(ParameterDatabase&&) = default;
    
    /// @brief copy assignment, this is a deep copy
    ParameterDatabase& operator=(const ParameterDatabase&);
    
    /// @brief default move assignemt
    ParameterDatabase& operator=(ParameterDatabase&&) = default;
    
    
    /// @brief construct and add one parameter without given range
    ///
    /// T can be bool, int, size_t, double, or std::string.
    template <typename T>
    void add(std::string name, T value, std::string description);
    
    /// @brief construct and add a vector-valued parameter without given range.
    ///
    /// T can be bool, int, size_t, double, or std::string.
    template <typename T>
    void add(std::string name, std::vector<T> value, std::string description);
    
    /// @brief construct and add a vector-valued parameter without given range.
    ///
    /// T can be bool, int, size_t, double, or std::string.
    /// This exists only for convenience, because without it often template 
    /// argument deduction/substitution fails.
    template <typename T>
    void add(std::string name, std::initializer_list<T> value, 
             std::string description);
    
    /// @brief construct and add one parameter with a given interval as range
    ///
    /// T can be int, size_t, or double. The range is `[min, max]`.
    template <typename T>
    void add(std::string name, T value, std::string description, T min, T max);
    
    /// @brief construct and add one vector-valued parameter with a given 
    /// interval as range.
    ///
    /// T can be int, size_t, or double. The range is `[min, max]`.
    template <typename T>
    void add(std::string name, std::vector<T> value, std::string description, 
             T min, T max);
    
    /// @brief construct and add one parameter with a given set as range
    ///
    /// T can be bool, int, size_t, or string.
    template <typename T>
    void add(std::string name, T value, std::string description, 
             std::set<T> range);
    
    /// @brief construct and add one vector-valued parameter with a given set 
    /// as range.
    ///
    /// T can be bool, int, size_t, or string.
    template <typename T>
    void add(std::string name, std::vector<T> value, std::string description, 
             std::set<T> range);
    
    /// @brief add a parameter (move)
    void add(Parameter&& p);
    /// @brief add a parameter (copy)
    void add(const Parameter & p);
    
    /// @brief return parameter with a given name
    const Parameter& operator[](std::string parameter_name) const;
    Parameter& operator[](std::string parameter_name);
    
    /// @brief return a parameter's value if present, or a default value otherwise
    template <typename T>
    T try_get_value(std::string parameter_name, T default_value) const;
    
    /// @brief return the name of this database
    std::string get_name() const;
    
    /// @brief change the name of this database
    void set_name(const std::string&);
    
    /// @brief get the number of parameters in this database
    size_t get_n_parameters() const;
    
    /// @brief find out if a parameter with a given name exists in this database
    /// 
    /// @note This does not search for such a parameter in nested databases.
    /// better name? one would write e.g.: if(db.contains("param_name"))
    bool contains(const std::string& name) const;
    
    /// @brief add a nested parameter database
    void add_nested_database(ParameterDatabase db);
    
    /// @brief Kindly ask, whether a nested database of the name
    ///        'name' is contained.
    bool has_nested_database(const std::string& name) const;

    /// @brief return additional parameter database with a given name.
    const ParameterDatabase& get_nested_database(const std::string& name) const;
    ParameterDatabase& get_nested_database(const std::string& name);
    
    /// @brief get the number of nested databases in this database
    size_t get_n_nested_databases() const;
    
    /// @brief write database to a stream, which can be read again
    ///
    /// The idea is to use this function to write the database into a file which
    /// then can be used to run the program again.
    ///
    /// Set the parameter verbose to 
    /// - false, if you only want the parameters,
    /// - true, if you want the parameters, their description (documentation), 
    ///   and range.
    ///
    /// Additionally to the above, the date, the database name, some git revision
    /// information and the host name is printed whenever 'include_general_info'
    /// is set to true.
    ///
    /// Nested databases are printed at the end. This means that nesting with 
    /// two or more levels is flattened as if all nested datases were stored
    /// directly in this one. This is where read-write will produce different 
    /// results.
    void write(std::ostream& stream, bool verbose = false,
               bool include_general_info = true) const;
    /// @brief convenience function which calls write(std::ostream&, bool, bool)
    void write(const std::string& filename, bool verbose = false,
               bool include_general_info = true) const;
    
    
    /// @brief read parameters from a stream
    ///
    /// There are formatting restrictions. In general try to use 
    /// ParameterDatabase::write to get a conforming stream. 
    ///
    /// The name of the database is surrounded by '[' and ']' where the first
    /// character of that line must be '['. In one stream there can be multiple
    /// databases to be read. They must then be seperated by lines indicating 
    /// names, e.g.
    /// \code
    /// [name of first database]
    /// parameter_name_1: value
    /// parameter_name_2: value
    /// 
    /// [name of second database]
    /// parameter_name_3: value
    /// parameter_name_4: value
    /// \endcode
    /// The parameters for the database are all between the first occurence of 
    /// '[_some_name_]' and the second (or the end of the stream respectively).
    /// That means any parameter before the first '[_some_name_]' is ignored!
    /// If no name is found at all, then a default name is given and all lines 
    /// are considered to be read. In case there are more than one database, 
    /// all remaining ones will be stored as nested databases.
    ///
    /// Each parameter must be on one line. Empty lines and lines without a 
    /// colon (':') are ignored. Parameters must not have spaces in their names.
    /// The general form is one of the following
    ///     parameter_name: value
    ///     parameter_name: value [ range_min, range_max ]
    ///     parameter_name: value { range_1, range_2, range_3 }
    ///     parameter_name: (v1, v2, v3)
    ///     parameter_name: (v1, v2, v3) [ range_min, range_max ]
    ///     parameter_name: (v1, v2, v3) { range_1, range_2, range_3 }
    ///
    /// There must be a colon (':') after the name. For booleans the value can 
    /// be 'true' or 'false'. The vectors (indicated by '(' and ')') are 
    /// translated to std::vector, so can have arbitrary length.
    ///
    /// All parameters can have a range. That means the parameter will never
    /// be outside some specified values (e.g. an interval or some set of 
    /// values). Ranges for boolean and string parameters should be determined 
    /// using braces  (e.g. '{true}' or '{ true, false }'). For parameters of 
    /// type int or size_t either braces ('{}') or brackets ('[]') are ok, one
    /// indicating a set of individual entries and the other an interval. For 
    /// double parameters only brackets (e.g. '[-2.5, 5.7]') are possible, 
    /// indicating an interval. There has to be a space (' ') between the value
    /// and the range. For vector valued parameters each entry must be in the 
    /// specified range.
    ///
    /// A parameter object is always in a consistent state, meaning that the
    /// value(s) are always within the specified range.
    ///
    /// It is possible to read documentation for each parameter which will then
    /// become its Parameter::description. All lines directly before the line 
    /// with the parameter which start with a '#' are considered. The line 
    /// before the documentation should be either empty, some other parameter, 
    /// or consist of anything but a colon (':'). If you want to document the
    /// Parameters in your input stream, but don't want that to be read, use 
    /// two '#' instead of just one. This is how the default documentation in
    /// ParMooN can be kept inside the Parameter objects.
    void read(std::istream& is);
    /// @brief convenience function which calls read(std::istream&)
    void read(const std::string& filename);
    
    /// @brief merge another database into this one
    ///
    /// Parameters which exist in this one get the values from the other
    /// database. All parameters in `other` which do not exist in this database,
    /// are copied into this one only if `create_new_parameters` is true.
    /// If 'recursive' is true, nested databases are copied into this one. If a
    /// nested database with the same name already exists and 'recursive' is 
    /// true, 'merge' is called.
    void merge(const ParameterDatabase &other,
               bool create_new_parameters = true,
               bool recursive = false);
    
    /// @brief out some information on the parameters
    ///
    /// Setting `only_names_and_values` to true prints only little information
    /// on each parameter. Setting it to false prints all details.
    void info(bool only_names_and_values = true) const;
  
  private:
    /// @brief name of this parameter database
    std::string name;
    
    // We need a container of Parameter objects with the following properties:
    // - easy insertion (easy deleting is not so much of importance)
    // - unique elements, no two Parameter objects shall have the same name
    // - full access to the Parameters (not just const)
    // A std::set<Parameter> with a custom compare function seems to be a good 
    // candidate to store all the Parameter objects. However a set will not let
    // you access the items as non-const references, as this might break the 
    // ordering. Another idea would be std::map<std::string, Parameter>, but
    // then the name of each parameter would be safed as the key and in the 
    // Parameter class, not so nice. These would have to be kept in sync. So it
    // seems std::list<Parameter> is ok as long as one makes sure no two entries
    // have the same name.
    // Note that finding a parameter in a list is slower than in a set::set or
    // std::map. However we expect this list to be rather short, so there 
    // should be no performance problem.
    
    /// @brief the set of all parameters stored in this parameter database
    std::list<Parameter> parameters;
    
    /// @brief additional parameter databases (identified by their unique name)
    std::list<ParameterDatabase> databases;
    
    /// @brief a helper function which is called by the other 'read' method.
    /// 
    /// This allows for a simpler implementation of nested database (where root
    /// is false).
    void read(std::istream& is, bool root);
};

#endif // __PARAMETERDATABASE__
