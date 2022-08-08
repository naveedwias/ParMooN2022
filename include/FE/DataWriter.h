#ifndef __DATAWRITER__
#define __DATAWRITER__

#include "all_defines_external_libraries.h"
#include "ParameterDatabase.h"
#include "templateNames.h"
#include <vector>
#ifdef _MPI
#include "mpi.h"
#endif
#ifdef PARMOON_WITH_HDF5
#include "hdf5.h"
#else
typedef int64_t hid_t;
#endif

class TCollection;

enum class xdmf_data_formats{xml, binary, hdf5};
enum class xdmf_time_types{list, hyperslab};

/**
 * @brief store given data and write output, Supported formats: xdmf, vtk
 *
 * The recommended format is xdmf, while vtk can be considered obsolete and will
 * be removed in some time.
 *
 * Using xdmf will produce an xml-file containing the meta (light) data. In 
 * addition the heavy data is written to hdf5 or binary files (depending on
 * DataWriter<d>::xdmf_format). It is possible to write all data in text form
 * into the xdmf file (DataWriter<d>::xdmf_format = xdmf_data_formats::xml), but
 * this is only meaningful for very small grids and few time steps.
 */
template <int d>
class DataWriter
{
  public:
    ///@brief default constructor: parameters are copied from Database
    explicit DataWriter(const ParameterDatabase& param_db);

    using FEFunction = typename Template_names<d>::FEFunction;
    using FEVectFunct = typename Template_names<d>::FEVectFunct;

    /// @brief Add a finite element function to this object which will be
    /// written when write() or write(double) is called.
    /// @note Additional finite element functions must live on the same grid
    /// (TCollection). This class can only write one grid.
    /// @note In parallel mode, the underlying vector should be in consistency 
    /// level 1
    void add_fe_function(const FEFunction* fefunction);

    /// @brief Add a finite element vector function to this object which will be
    /// written when write() or write(double) is called.
    /// @note Additional finite element vector functions must live on the same
    /// grid (TCollection). This class can only write one grid.
    /// @note In parallel mode, the underlying vector should be in consistency 
    /// level 1
    void add_fe_vector_function(const FEVectFunct* fevectfunct);

    /// @brief Add a global variable to this object which will be written
    /// when write() or write(double) is called. Use this for data such
    /// as flow statistics etc.
    /// @note This is a function pointer, evaluated at write time. Null pointers
    /// default to zero.
    /// @note In MPI mode, this is called on every process, but only the return
    /// value from the root process is written.
    /// @note This currently only works in XDMF/HDF5 mode. The value of each global
    /// output variable is written as a "Grid"-centered attribute.
    void add_global_output_variable(const std::string& name,
      std::function<double()> func);

    /** @brief write data to files during time-dependent problems.
     *
     * You need to have added at least one finite element function or vector
     * function, see add_fe_function() and add_fe_vector_function(),
     * respectively.
     * @param[in] t the physical time
     */
    void write(double t);

    /** @brief write data files (mainly for stationary problems).
     * You need to have added at least one finite element function or vector
     * function, see add_fe_function() and add_fe_vector_function(),
     * respectively.
     */
    void write();

    /** @brief write extra data to files.
     *
     * You need to have added at least one finite element function or vector
     * function, see add_fe_function() and add_fe_vector_function(),
     * respectively.
     *
     *  This function does not check if a file with the same index already
     *  exists, in this case the old file will be replaced.
     *
     * @param[in] index in the file naming.
     */
    void write(int index);

    std::string get_base_name() { return base_name; };

  // -------------------------------------------------------------------------
  protected:

    /// @name parameters taken from the database during construction.
    //@{
    /// @brief output directory. All files written by this class go here.
    std::string output_dir;

    /// @brief base name of output files.
    /// The file name of every file written by this class starts with this
    /// string.
    std::string base_name;

    /// @brief write vtk file whenever write() or write(double) is called.
    bool writeVTK;

    /// @brief write xdmf file together with heavy data whenever write() or
    /// write(double) is called.
    bool write_xdmf;

    /// @brief write csv file tabulating the global output variables
    /// whenever write() or write(double) is called.
    bool write_csv;

    /// @brief write separate h5 files for the solutions in each time step.
    /// This is not used for stationary problems.
    bool separate_h5;

    /// @brief enable DEFLATE compression for all created h5 files.
    /// This sets the maximum compression level 9.
    bool compress_h5;

    /// @brief Use collective I/O when writing h5 files.
    /// On some filesystems, notably NFS, collective I/O can be unreliable,
    /// even unstable. If false, collates and writes all data on the root,
    /// for a gain in safety at the expense of I/O performance.
    bool collective_h5;

    /// @brief Write global output variables to XDMF files.
    /// With lots of global output variables, XDMF files may become large.
    /// Depending on use case, it may be more sensible to write them to CSV
    /// files only.
    bool globals_xdmf;

    /// @brief in mpi context, write the subdomain as a scalar function.
    /// This helps to see the subdomain of each processor in mpi mode. Together
    /// with the grid this is only written once and referenced in each time
    /// step.
    bool write_subdomain_info;

    /// @brief the format in which heavy data is written
    xdmf_data_formats xdmf_format;

    /// @brief the format in which time steps are written to file.
    /// In case of a constant time step length 'xdmf_time_types::hyperslab'
    /// makes much sense, because it stores only the initial time, the time step
    /// length, and the number of time steps. Otherwise 'xdmf_time_types::list'
    /// has to be used, which simply writes each time step to the xdmf file as a
    /// list.
    xdmf_time_types xdmf_timetype;

    /// @brief initial time of simulation for time-dependent problems, used only
    /// if DataWriter<d>::xdmf_time_type is xdmf_time_types::hyperslab.
    double time_start = 0;

    /// @brief time step length of simulation for time-dependent problems, used
    /// only if DataWriter<d>::xdmf_time_type is xdmf_time_types::hyperslab.
    double time_step_length = 0;

    /// @brief number of iterations between two outputs.
    size_t n_steps_per_output;
    //@}

    /** @brief representation of the xdmf-file.
     * It is updated and rewritten in every time step if 
     * DataWriter<d>::write_xdmf is true.
     */
    std::vector<std::string> xdmf_buffer;
    /// @brief index in the DataWriter<d>::xdmf_buffer where the time steps are
    /// declared. The corresponding entry in the DataWriter<d>::xdmf_buffer is
    /// updated in every time step.
    size_t buffer_ind;
    /// @brief counting the steps until next output is written
    int n_steps_until_next_output;

    /// @brief array of stored finite element functions
    std::vector<const FEFunction*> fe_functions;
    /// @brief array of stored vector-valued finite element functions
    std::vector<const FEVectFunct*> fe_vector_functions;

    /// @brief array of global variable names
    std::vector<std::string> global_output_variable_names;

    /// @brief array of global variable function pointers
    std::vector<std::function<double()>> global_output_variable_funcs;

    bool has_example_output = false;

    /// @brief do or do not write discontinuous data.
    /// This is set to true whenever one of the involved finite element spaces
    /// is discontinuous. Each vertex is listed as often as it appears in
    /// different cells. So the number of written vertices is
    /// 'numper_of_cells * number_of_vertices_per_cell', which means that all
    /// solutions are written as vectors with this length. Therefore, this makes
    /// the files much larger. If this member variable is false, a P1/Q1 
    /// projection is written.
    bool write_discontinuous_functions;

    ///@brief geometry (collection of cells)
    const TCollection* coll;

    /// @brief number of vertices, equals the size of solution arrays.
    /// Note that every vertex is counted multiple times if
    /// DataWriter<d>::write_discontinuous_functions is true. This is computed
    /// exactly once. This number is meant per process in mpi mode.
    unsigned int local_n_vertices = 0;

    /// @brief number of vertices, equals the size of solution arrays.
    /// This number is global over all processors in mpi mode and otherwise 
    /// equal to DataWriter<d>::local_n_vertices.
    unsigned int global_n_vertices = 0;

    /// @brief the number of vertices on processors with a smaller rank.
    /// This is only used in mpi mode.
    unsigned int global_vertex_offset = 0;

    /// @brief part of the xdmf (xml) file for each time step. Only used if
    /// DataWriter<d>::write_subdomain_info is true, i.e., in mpi mode.
    std::string subdomain_info_xml = "";

    /** @brief vector containing the (physical) output times
     * This array allows to display the physical time when looking at results
     * and is only filled for time-dependent problems, i.e., when calling 
     * DataWriter<d>::write(double) rather than DataWriter<d>::write().
     */
    std::vector<double> timeValues;

    ///@brief a boolean to manage output when restarting a simulation
    bool continue_output_after_restart;

    ///@brief the restart time (=initial time of the restarted simulation)
    double restart_time;

    // -------------------------------------------------------------------------

    /// @brief write vtk file with the given name.
    void writeVtk(const std::string& name);

    /// @brief helper methods used in DataWriter<d>::writeVtk()
    //@{
    /// @brief Write mesh in vtk format
    /// @note see also TCollection::writeMesh() for a different file format.
    void writeMesh(const std::string& name);

    /** @brief writes the coordinates of the vertices */
    void writeCoord(std::ofstream& f);

    /** @brief writes x-component for all vertices and then y- component for all
     * vertices*/
    void printVectCompwise(std::ofstream& dat, const std::string& name,
                           unsigned int N_Vertices, unsigned int N_Comp,
                           const std::vector<double>& solutionAtNode);

    /** @brief writes the vector (x- and y-component) for all vertices*/
    void printVectPointwise(std::ofstream& dat, const std::string& name,
                            unsigned int N_Vertices, unsigned int N_Comp,
                            const std::vector<double>& solutionAtNode);

    /** @brief writes the absolute value for all vertices*/
    void printVectAbsValue(std::ofstream& dat, const std::string& name,
                           unsigned int N_Vertices, unsigned int N_Comp,
                           const std::vector<double>& solutionAtNode);
    //@}

    void retrieve_example_output();

#ifdef _MPI
    /** @brief write stored data into a pvtu and vtu files (XML files for
     * paraview)
     */
    void Write_ParVTK(int img, std::string directory, std::string basename);
#endif

    void write_csv_file(const std::string& file_name);

    /// @brief write the xdmf file and all files (heavy data) it relies on.
    /// After this method ran the xdmf file is fully functional even during a 
    /// time-dependent simulation. It is overwritten in each time step.
    void write_xdmf_file(const std::string& file_name);

    /// @name helper methods used in DataWriter<d>::write_xdmf_file()
    //@{
    /// @brief write the heavy data for the grid (vertices+cells).
    /// The returned string is the part of xdmf file realated to the grid. This
    /// method is called exactly once, even in a time-dependent problem.
    std::string write_grid();

    /// @brief write the vertices to the h5/binary/xml file depending on the
    /// 'xdmf_format'. This sets the member variable 'n_vertices' and should 
    /// therefore be called before 'write_topology' and 'write_solution'.
    std::string write_geometry(hid_t h5gridfile);
    std::string write_topology(hid_t h5gridfile);

    /// @brief write the heavy data for all solutions.
    /// This is called once for every time step in a time-dependent setting. The 
    /// returned string is the part of the xdmf file related to the solution.
    /// @param[in] only_xml if true no heavy data is written, useful to write 
    /// xdmf file after restart.
    std::string write_solution(bool only_xml = false) const;
    std::string write_global_solution(hid_t& h5file, std::string h5_file_name,
                                      bool only_xml) const;
    std::string write_scalar_solution(hid_t& h5file, std::string h5_file_name,
                                      bool only_xml) const;
    std::string write_vector_solution(hid_t& h5file, std::string h5_file_name,
                                      bool only_xml) const;
    std::string write_solutions_before_restart();

    /** @brief method computes the value in of a given function in each vertex.
     * @param[in] function FEFunction[23]D or FEVectFunction[23]D function,
     * @param[out] solutionAtNode double vector with function values,
     * @param[out] dimension 1 for scalar functions, d for vector-valued ones
     */
    template <class T>
    void computeNodeValues(const T* function, std::vector<double>& solutionAtNode,
                           unsigned int& dimension) const;

    /// @brief helper function to return the vector of vertex coordinates.
    /// The length of this vector depends on 'write_discontinuous_functions' and
    /// the space dimension.
    std::vector<double> get_vertex_coord() const;
    /// @brief helper function to return the vector of cells.
    /// Each cell has an identifier and the vertex indices.
    /// @param[in] offset used in mpi mode to get the globally correct vertices
    std::vector<unsigned int> get_topology(unsigned int n_vertices,
                                           unsigned int offset = 0) const;
    //@}

};

typedef DataWriter<3> DataWriter3D;
typedef DataWriter<2> DataWriter2D;

#endif // __DATAWRITER__
