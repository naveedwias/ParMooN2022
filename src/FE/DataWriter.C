#include <array>
#include <vector>
#include <DataWriter.h>
#include <Example_Output.h>
#ifdef __2D__
#include <FEVectFunct2D.h>
#else
#include <FEVectFunct3D.h>
#endif
#include <Vertex.h>
#include "BaseCell.h"
#include <algorithm>
#include <sys/stat.h>
#include <type_traits>
#include <assert.h>
#include <TimeDiscretizations.h>
#include <Chrono.h>
#include <cmath>
#include <cstring> // memset
#include <stdexcept>
#include <limits>
#include <unordered_set>
/// @brief check if deprecated parameters have been provided and issue a warning
/// in such a case. Additionally if the 'output_directory' is not set and the
/// deprecated 'output_vtk_directory' is set, then save the latter as the former.
void check_for_old_parameters(const ParameterDatabase & param_db,
                              ParameterDatabase db);
unsigned int get_n_vertices_multiply_counted(const TCollection* coll);

template <int d>
DataWriter<d>::DataWriter(const ParameterDatabase& param_db)
  : output_dir(), base_name(), n_steps_per_output(1), xdmf_buffer(),
    buffer_ind(), n_steps_until_next_output(0), fe_functions(),
    fe_vector_functions(), write_discontinuous_functions(false), coll(nullptr),
    timeValues(), continue_output_after_restart(false), restart_time(0)
{
  // set the variables depending on input parameters
  ParameterDatabase db = ParameterDatabase::default_output_database();
  db.merge(param_db, false);

  check_for_old_parameters(param_db, db);

  writeVTK = db["output_write_vtk"];
  write_xdmf = db["output_write_xdmf"];
  write_csv = db["output_write_csv"];

  separate_h5 = db["separate_h5_files"];
  compress_h5 = db["output_compress_h5_files"];
  collective_h5 = db["output_collective_h5_files"];
  globals_xdmf = db["output_xdmf_globals"];
  write_subdomain_info = db["output_write_subdomain_info"];

  if (write_xdmf)
  {
    size_t format = db["output_xdmf_format"];
#ifdef _MPI
    if (format != 2)
    {
      Output::info("DataWriter", "writing xdmf-files in parallel requires the "
                   "use of the hdf5 format for the heavy data, i.e., "
                   "'output_xdmf_format' is changed from ", format,
                   " to 2 here.");
      format = 2;
    }
#ifdef PARMOON_WITH_HDF5
    unsigned int h5maj, h5min, h5rel;
    H5get_libversion(&h5maj, &h5min, &h5rel);

    if (compress_h5)
    {
      if (h5maj < 1 ||
         (h5maj == 1 && (h5min < 10 || (h5min == 10 && h5rel < 2))))
      {
        Output::root_info("DataWriter",
            "compression of h5 files in mpi mode is only supported for newer "
            "versions of the hdf5 library (at least 1.10.2).");
        compress_h5 = false;
      }
    }

    Output::root_info<4>("DataWriter", "HDF5 library version: ", h5maj, ".",
                         h5min, ".", h5rel);
#ifndef H5_HAVE_PARALLEL
    // note: this should be checked during configure time by cmake already
    ErrThrow("We need the parallel version of hdf5 when using mpi.");
#endif // not H5_HAVE_PARALLEL
#endif // PARMOON_WITH_HDF5
#endif // MPI
    if (format == 0)
    {
      xdmf_format = xdmf_data_formats::xml;
    }
    else if (format == 1)
    {
      xdmf_format = xdmf_data_formats::binary;
    }
    else if (format == 2)
    {
#ifndef PARMOON_WITH_HDF5
      ErrThrow("Cannot use HDF5 format, because HDF5 library was not found. "
               "Set 'output_write_vtk' to 1 to write plain text files.");
#endif
      xdmf_format = xdmf_data_formats::hdf5;
    }
    else
    {
      ErrThrow("unknown output_xdmf_format ", format);
    }
  }

  size_t timetypes = db["output_xdmf_timetype"];
  if (timetypes == 0)
  {
    xdmf_timetype =  xdmf_time_types::list;
  }
  else if (timetypes == 1)
  {
    if (param_db.contains("time_start") && param_db.contains("time_step_length"))
    {
      xdmf_timetype =  xdmf_time_types::hyperslab;
      time_start = param_db["time_start"];
      time_step_length = param_db["time_step_length"];
    }
    else
    {
      xdmf_timetype =  xdmf_time_types::list;
    }
  }
  else
  {
    ErrThrow("unknown output_xdmf_timetype ", timetypes);
  }

  base_name = db["output_basename"].get<std::string>();
  output_dir = db["output_directory"].get<std::string>();
  n_steps_per_output = std::max<size_t>(1, db["steps_per_output"]); // avoid 0

  continue_output_after_restart = db["continue_output_after_restart"];
  if (continue_output_after_restart)
  {
    // some consistency check
    if ( !param_db.contains("time_start")
      || !param_db.contains("time_step_length"))
    {
      ErrThrow("You want to restart a simulation, then you have to provide the "
               "parameters 'time_start' and 'time_step_length' to the "
               "DataWriter constructor."
      );
    }

    restart_time = param_db["time_start"];
    if (restart_time == 0.)
    {
      ErrThrow("The parameter continue_output_after_restart is TRUE, but the "
      "start time is 0! Did you forget to set time_start to the actual restart time? "
      "It should be equal to the time of the binary you are starting from. If you "
      "really want to restart from 0, put the boolean continue_output_after_restart "
      "to FALSE.");
    }

    if (xdmf_format == xdmf_data_formats::xml && write_xdmf)
    {
      ErrThrow("restarting simulation with xdmf output in xml (plain text) "
               "format is not supported. Use 'output_xdmf_format' 1 (binary) "
               "or 2 (hdf5).");
    }

    if (writeVTK && !write_xdmf)
    {
      write_solutions_before_restart(); // resets timeValues
    }
  }

  // create directory where all the output files from this class instance go
  int my_rank = 0;

#ifdef _MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

  if (my_rank == 0 && (writeVTK || write_xdmf || write_csv))
  {
    mkdir(output_dir.c_str(), 0777);
  }
}

template <int d>
void DataWriter<d>::add_fe_function(const FEFunction* fefunction)
{
  // check that FE functions have the same collection
  if(this->coll == nullptr)
  {
    this->coll = fefunction->GetFESpace()->GetCollection();
  }
  else
  {
    if(this->coll != fefunction->GetFESpace()->GetCollection())
    {
      // we could also just refuse to add this fe function and return without an
      // exception.
      ErrThrow("new FE function has a different collection");
    }
  }
  bool already_known = std::any_of(
                         fe_functions.begin(), fe_functions.end(),
      [fefunction](const FEFunction* f) { return f == fefunction; });
  if(!already_known)
  {
    fe_functions.push_back(fefunction);
    write_discontinuous_functions |= fefunction->GetFESpace()->
      is_discontinuous();
  }
}

template <int d>
void DataWriter<d>::add_fe_vector_function(const FEVectFunct* fevectfunction)
{
  // check that FE functions have the same collection
  if(this->coll == nullptr)
  {
    this->coll = fevectfunction->GetFESpace()->GetCollection();
  }
  else
  {
    if(this->coll != fevectfunction->GetFESpace()->GetCollection())
    {
      // we could also just refuse to add this fe function and return without an
      // exception.
      ErrThrow("new FE vector function has a different collection");
    }
  }
  bool already_known = std::any_of(
                         fe_vector_functions.begin(), fe_vector_functions.end(),
      [fevectfunction](const FEVectFunct* f) { return f == fevectfunction; });
  if(!already_known)
  {
    fe_vector_functions.push_back(fevectfunction);
    write_discontinuous_functions |= fevectfunction->GetFESpace()->
      is_discontinuous();
  }
}

template <int d>
void DataWriter<d>::add_global_output_variable(const std::string& name,
  std::function<double()> func)
{
  global_output_variable_names.push_back(std::string(name));
  global_output_variable_funcs.push_back(func);
}

template <int d>
void DataWriter<d>::write()
{
  if (!writeVTK && !write_xdmf && !write_csv)
  {
    return; // return early to avoid unnecessary output
  }

  retrieve_example_output();

  Chrono timer;
  std::string name = output_dir + "/" + base_name;

  if (writeVTK)
  {
#ifdef _MPI
    Write_ParVTK(0, output_dir, base_name);
#else
    Output::root_info<2>("DataWriter", "writing ", name + ".vtk");
    writeVtk(name + ".vtk");
#endif
  }

  if (write_xdmf)
  {
    Output::root_info<2>("DataWriter", "writing ", name + ".xdmf");
    write_xdmf_file(name + ".xdmf");
  }

  if (write_csv)
  {
    Output::root_info<2>("DataWriter", "writing ", name + ".csv");
    write_csv_file(name + ".csv");
  }

  timer.stop_and_print("writing vtk/xdmf files");
}

/// @brief create a string representing 'n' with prepended zeros, so that the
/// length of the string is (at least) 'n_digits'. If 'n' is too large to fit
/// within 'n_digits', the returned string is longer.
// see https://stackoverflow.com/a/225435
std::string num2str(int n, int n_digits)
{
  std::stringstream ss;
  ss << std::setw(n_digits) << std::setfill('0') << n;
  return ss.str();
}

template <int d>
void DataWriter<d>::write(double current_time)
{
  if (n_steps_until_next_output != 0)
  {
    n_steps_until_next_output--;
  }
  else
  {
    retrieve_example_output();

    n_steps_until_next_output = n_steps_per_output - 1;
    timeValues.push_back(current_time);
    // index of current time step in 'timeValues'
    auto index = timeValues.size() - 1;
    std::string name = output_dir + "/" + base_name;
    if (write_xdmf)
    {
      std::string xdmf_file_name = name;
      xdmf_file_name += ".xdmf";

      Output::root_info<2>("DataWriter", "writing ", xdmf_file_name);

      write_xdmf_file(xdmf_file_name);
    }

    if (writeVTK)
    {
#ifdef _MPI
      Write_ParVTK(index * n_steps_per_output, output_dir, base_name);
#else
      std::string vtk_name = name;
      vtk_name += num2str(index * n_steps_per_output, 8);
      vtk_name += ".vtk";

      Output::root_info<2>("DataWriter", "writing ", vtk_name);

      writeVtk(vtk_name);
#endif
    }

    if (write_csv)
    {
      Output::root_info<2>("DataWriter", "writing ", name + ".csv");
      write_csv_file(name + ".csv");
    }
  }
}

template <int d>
void DataWriter<d>::write(int index)
{
  retrieve_example_output();

  timeValues.push_back(index);
  std::string name = output_dir + "/" + base_name;

  if (write_xdmf)
  {
    std::string xdmf_file_name = name;
    xdmf_file_name += ".xdmf";
    Output::root_info<2>("DataWriter", "writing ", xdmf_file_name);
    write_xdmf_file(xdmf_file_name);
  }

  if (writeVTK)
  {
#ifdef _MPI
    Write_ParVTK(index, output_dir, base_name);
#else
    std::string vtk_name = name;
    //std::string idx = static_cast<std::ostringstream*>
    //                                ( &(std::ostringstream() << index) )->str();
    std::string idx = std::to_string(index);
    vtk_name += idx;
    vtk_name += ".vtk";
    Output::root_info<2>("DataWriter", "writing ", vtk_name);
    writeVtk(vtk_name);
#endif
  }

  if (write_csv)
  {
    Output::root_info<2>("DataWriter", "writing ", name + ".csv");
    write_csv_file(name + ".csv");
  }
}

template <int d>
int n_local_vertices_to_type(int n_loc_vert)
{
  switch(d)
  {
    case 2:
      switch(n_loc_vert)
      {
        case 4: return 9; break;
        case 3: return 5; break;
      }
      break;
    case 3:
      switch(n_loc_vert)
      {
        case 4: return 10; break;
        case 8: return 12; break;
      }
      break;
  }
  ErrThrow("a ", d, "D cell with ", n_loc_vert, " vertices is not supported");
}

template <int d>
void DataWriter<d>::writeMesh(const std::string& name)
{
  std::ofstream dat(name);
  if(!dat)
  {
    ErrThrow("cannot open file for output. ", name);
  }
  dat << setprecision(12);

  dat << "# vtk DataFile Version 4.0\n";
  dat << "file created by ParMooN"
      << " Time < " << (timeValues.empty() ? 0. : timeValues[0]) << " >\n";
  dat << "ASCII\n";
  dat << "DATASET UNSTRUCTURED_GRID\n";
  dat << "POINTS " << coll->GetN_Vertices() << " double\n";

  writeCoord(dat);

  unsigned int N_LocVertices = coll->GetNLocVertices();
  unsigned int N_Elements = coll->GetN_Cells();

  dat << "CELLS " << N_Elements << " " << N_Elements + N_LocVertices << "\n";
  for(unsigned int i = 0; i < N_Elements; i++)
  {
    int N_CellVertices = coll->GetCell(i)->GetN_Vertices();
    dat << N_CellVertices << " ";
    for(int j = 0; j < N_CellVertices; j++)
    {
      dat << this->coll->GetGlobalVerNo(i, j) << " ";
    }
    dat << "\n";
  }
  dat << "\n";

  dat << "CELL_TYPES " << N_Elements << "\n";
  for(unsigned int i = 0; i < N_Elements; i++)
  {
    int N_CellVertices = coll->GetCell(i)->GetN_Vertices();
    dat << n_local_vertices_to_type<d>(N_CellVertices) << " ";
  }
}

template <int d>
void DataWriter<d>::writeVtk(const std::string& name)
{
  bool tmp_write_discontinuous_functions = write_discontinuous_functions;
  if(write_discontinuous_functions)
  {
    Output::warn("DataWriter",
                 "Writing discontinuous functions into a vtk file will project "
                 "that function onto a P1/Q1 space, which is continuous. That "
                 "means you get average values at vertices. If you want to see "
                 "a discontinuous function, set 'output_write_xdmf' to true to "
                 "generate an xdmf-file."
                );
    // temporarily write continuous functions
    write_discontinuous_functions = false;
  }

  writeMesh(name);
  std::ofstream dat(name, std::ios::out | std::ios::app);
  if(!dat)
  {
    ErrThrow("cannot open file for output. ", name);
  }
  dat.setf(std::ios::fixed);
  dat << setprecision(12);

  local_n_vertices = coll->GetN_Vertices();
  unsigned int dimension;

  dat << "\n\n";
  dat << "POINT_DATA " << local_n_vertices << "\n";

  std::vector<double> solutionAtNode;
  for(unsigned int i = 0; i < fe_functions.size(); ++i)
  {
    computeNodeValues(fe_functions.at(i), solutionAtNode, dimension);
    std::string name = fe_functions.at(i)->GetName();
    printVectCompwise(dat, name, local_n_vertices, dimension, solutionAtNode);
    if(dimension > 1)
    {
      printVectAbsValue(dat, name, local_n_vertices, dimension, solutionAtNode);
      printVectPointwise(dat, name, local_n_vertices, dimension, solutionAtNode);
    }
  }
  for(unsigned int i = 0; i < fe_vector_functions.size(); ++i)
  {
    computeNodeValues(fe_vector_functions.at(i), solutionAtNode, dimension);
    std::string name = fe_vector_functions.at(i)->GetName();
    printVectCompwise(dat, name, local_n_vertices, dimension, solutionAtNode);
    printVectAbsValue(dat, name, local_n_vertices, dimension, solutionAtNode);
    printVectPointwise(dat, name, local_n_vertices, dimension, solutionAtNode);
  }
  dat << "\n";
  dat.close();
  write_discontinuous_functions = tmp_write_discontinuous_functions;
}

template <typename T>
std::string convert_to_string(T x)
{
  // std::to_string does not work for std::string and more importantly for
  // T=double it truncates to 6 digits, which is not what we want here
  std::ostringstream os;
  os << std::setprecision(std::numeric_limits< double >::digits10+1)
     << std::boolalpha << std::showpoint << x;
  std::string s = os.str();
  // remove trailing zeros (behind a decimal point). We leave one zero if it
  // is directly after the decimal point (to yield e.g. "1.0").
  std::size_t found = s.find_last_not_of(std::string("0"));
  auto npos = std::string::npos;
  bool has_decimal_point = (s.find_last_of(std::string(".")) != npos);
  bool has_exp_notation = (s.find_last_of(std::string("e")) != npos);
  if(has_decimal_point && found != npos && !has_exp_notation)
  {
    if(s[found] == '.' && s.size() >= found + 2)
      s.erase(found + 2);
    else
      s.erase(found + 1);
  }
  return s;
}

// convert xdmf_data_formats to the string which valid in xdmf-files
std::string as_string(xdmf_data_formats f)
{
  switch(f)
  {
    case xdmf_data_formats::hdf5:
    {
      return "HDF"; break;
    }
    case xdmf_data_formats::binary:
    {
      return "Binary"; break;
    }
    case xdmf_data_formats::xml:
    {
      return "XML"; break;
    }
    default:
      ErrThrow("unknown xdfm data format");
      break;
  }
}

std::string get_time_element(xdmf_time_types xdmf_timetype,
                             std::vector<double> timeValues, double time_start,
                             double time_step_length, double n_steps_per_output)
{
  std::string timetype;
  switch(xdmf_timetype)
  {
    case xdmf_time_types::list:
    {
      timetype = "List"; break;
    }
    case xdmf_time_types::hyperslab:
    {
      timetype = "HyperSlab"; break;
    }
  }
  std::string timeValuesString;
  std::string timeValuesDim;
  const unsigned int n_timesteps = timeValues.size();
  switch(xdmf_timetype)
  {
    case xdmf_time_types::list:
    {
      for(unsigned int i = 0; i < n_timesteps; i++)
      {
        timeValuesString += convert_to_string(timeValues[i]) + " ";
        if(i%20 == 19) // Linebreak every 20 lines.
        {
          timeValuesString += "\n";
        }
      }
      timeValuesDim = std::to_string(n_timesteps);
      break;
    }
    case xdmf_time_types::hyperslab:
    {
      timeValuesString = convert_to_string(time_start) + " "
          + convert_to_string(time_step_length * n_steps_per_output) + " "
          + std::to_string(n_timesteps);
      timeValuesDim = "3";
      break;
    }
  }

  std::string time_string =
    "      <Time TimeType=\"" + timetype + "\">\n"
    "        <DataItem Format=\"XML\" NumberType=\"Float\""
    " Dimensions=\"" + timeValuesDim + "\">\n"
    + timeValuesString + "\n"
    "        </DataItem>\n"
    "      </Time>\n";
  return time_string;
}

#ifdef PARMOON_WITH_HDF5

template <unsigned int filedim, typename T>
void write_h5_file(unsigned int d, unsigned int dim, unsigned int global_dim,
                   unsigned int
#ifdef _MPI
                   global_offset
#endif
                   ,
                   bool compressed, hid_t file,
                   const char* ds_name, const T* data)
{
  // Write and expand the datasets in the h5 files, in which the heavy data
  // of the grid and the solutions is stored.
  bool is_int = std::is_same<T, unsigned int>::value;

  // Create the h5 dataset with current and maximal possible dimensions of the
  // data to be stored. Then write the data into the file.
  hsize_t local_dims[filedim];
  hsize_t global_dims[filedim];
  hsize_t max_dims[filedim];

  local_dims[0] = dim;
  global_dims[0] = global_dim;
  max_dims[0] = global_dim;

  if (filedim == 2)
  {
    local_dims[1] = d;
    global_dims[1] = d;
    max_dims[1] = d;
  }

  hid_t prop = H5Pcreate(H5P_DATASET_CREATE);
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);

#ifdef _MPI
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Pset_dxpl_mpio_collective_opt(plist_id, H5FD_MPIO_COLLECTIVE_IO);
  hid_t node_memory_dataspace = H5Screate_simple(filedim, local_dims, nullptr);
#else
  hid_t node_memory_dataspace = H5S_ALL;
#endif

  if (compressed)
  {
    // compression only works for chunked data

#ifdef _MPI
    hsize_t chunk_dims[filedim];

    // parallel compression fails unless all chunks are the same size
    // use the largest local dim, for simplicity

    unsigned int max_dim = 0;

    MPI_Allreduce(&dim, &max_dim, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);

    chunk_dims[0] = max_dim;

    if (filedim == 2)
    {
      chunk_dims[1] = d;
    }

    H5Pset_chunk(prop, filedim, chunk_dims);

    // max dim has to be unlimited, since global_dim is not
    // usually a multiple of max_dim

    max_dims[0] = H5S_UNLIMITED;
#else
    H5Pset_chunk(prop, filedim, local_dims);
#endif
    H5Pset_deflate(prop, 9); // enables zlib compression.
  }

  hid_t space = H5Screate_simple(filedim, global_dims, max_dims);
  hid_t dataset = H5Dcreate2(file, ds_name,
                             (is_int ? H5T_NATIVE_UINT : H5T_NATIVE_DOUBLE),
                             space, H5P_DEFAULT, prop, H5P_DEFAULT);
  if (dataset < 0)
  {
    ErrThrow("unable to create a dataset '", ds_name, "' in the h5 file");
  }

  herr_t status;
#ifdef _MPI
  hsize_t offset[filedim] = {global_offset}; // if filedim==2 => second entry=0
  hid_t node_file_dataspace = H5Dget_space(dataset);
  H5Sselect_hyperslab(node_file_dataspace, H5S_SELECT_SET, offset,
                                 nullptr, local_dims, nullptr);
#else
  hid_t node_file_dataspace = H5S_ALL;
#endif

  status = H5Dwrite(dataset, (is_int ? H5T_NATIVE_UINT : H5T_NATIVE_DOUBLE),
                    node_memory_dataspace, node_file_dataspace, plist_id, data);
  if (status < 0)
  {
    ErrThrow("unable to write dataset ", ds_name);
  }

  status = H5Pclose(plist_id);
  if (status < 0)
  {
    ErrThrow("unable to close data xfer property list");
  }

  status = H5Dclose(dataset);
  if (status < 0)
  {
    ErrThrow("unable to close data set");
  }

  status = H5Sclose(space);
  if (status < 0)
  {
    ErrThrow("unable to close data space ");
  }

  status = H5Pclose (prop);
  if (status < 0)
  {
    ErrThrow("unable to close property list for data set");
  }
}

template <unsigned int filedim, typename T>
void write_h5_file_noncollective(unsigned int d, unsigned int dim,
                   bool compressed, hid_t file,
                   const char* ds_name, const T* data)
{
  // Write and expand the datasets in the h5 files, in which the heavy data
  // of the grid and the solutions is stored.
  bool is_int = std::is_same<T, unsigned int>::value;

  // Create the h5 dataset with current and maximal possible dimensions of the
  // data to be stored. Then write the data into the file.
  hsize_t dims[filedim];
  hsize_t max_dims[filedim];

  dims[0] = dim;
  max_dims[0] = dim;

  if (filedim == 2)
  {
    dims[1] = d;
    max_dims[1] = d;
  }

  hid_t prop = H5Pcreate(H5P_DATASET_CREATE);

  if (compressed)
  {
    // compression only works for chunked data

    hsize_t chunk_dims[filedim];

    chunk_dims[0] = 16384;

    if (filedim == 2)
    {
      chunk_dims[1] = d;
    }

    H5Pset_chunk(prop, filedim, chunk_dims);

    max_dims[0] = H5S_UNLIMITED;

    H5Pset_deflate(prop, 9); // enables zlib compression.
  }

  hid_t space = H5Screate_simple(filedim, dims, max_dims);
  hid_t dataset = H5Dcreate2(file, ds_name,
                             (is_int ? H5T_NATIVE_UINT : H5T_NATIVE_DOUBLE),
                             space, H5P_DEFAULT, prop, H5P_DEFAULT);
  if (dataset < 0)
  {
    ErrThrow("unable to create a dataset '", ds_name, "' in the h5 file");
  }

  herr_t status = H5Dwrite(dataset, (is_int ? H5T_NATIVE_UINT : H5T_NATIVE_DOUBLE),
                    H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  if (status < 0)
  {
    ErrThrow("unable to write dataset ", ds_name);
  }

  status = H5Dclose(dataset);
  if (status < 0)
  {
    ErrThrow("unable to close data set");
  }

  status = H5Sclose(space);
  if (status < 0)
  {
    ErrThrow("unable to close data space ");
  }

  status = H5Pclose(prop);
  if (status < 0)
  {
    ErrThrow("unable to close property list for data set");
  }
}

#endif // PARMOON_WITH_HDF5

void write_buffer_to_xdmf(const std::string& file_name,
                          const std::vector<std::string>& xdmf_buffer)
{
  // Create the XDMF file and write the lines from the buffer into the file.
  // The buffer is filled with the lines in the write_xdmf_file fuction. The
  // created xdmf-file can be opened in ParaView.
  std::ofstream f(file_name);
  f << std::fixed << std::scientific << std::setprecision(16);
  if(!f)
  {
    ErrThrow("cannot open file for output. ", file_name);
  }

  for(auto& line : xdmf_buffer)
  {
    f << line << "\n";
  }
}

template<int d>
std::string DataWriter<d>::write_grid()
{
  bool rank_has_file = true;

#ifdef _MPI
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  rank_has_file = mpi_rank == 0 || collective_h5;
#endif

  // write the actual h5 file containing the grid information (vertices+cells)
  // if xdmf_format is not hdf5, the 'h5gridfile' remains unused
  hid_t h5gridfile = 0;
#ifdef PARMOON_WITH_HDF5
  hid_t plist_id = H5P_DEFAULT;
#endif // #ifdef PARMOON_WITH_HDF5
  if (xdmf_format == xdmf_data_formats::hdf5 && !continue_output_after_restart)
  {
#ifdef PARMOON_WITH_HDF5
    std::string h5grid = output_dir + "/" + base_name + "_grid.h5";

#ifdef _MPI
    if (collective_h5)
    {
      plist_id = H5Pcreate(H5P_FILE_ACCESS);
      H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    }
#endif // _MPI

    if (rank_has_file)
    {
      Output::print<4>("creating h5 file ", h5grid);
      h5gridfile = H5Fcreate(h5grid.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, plist_id);
      if (h5gridfile < 0)
      {
        ErrThrow("Error creating the h5 file ", h5grid);
      }
    }
#endif // #ifdef PARMOON_WITH_HDF5
  }
  // Topology specifies the connectivity of the grid,
  // while Geometry specifies the location of the grid nodes.
  std::string grid_string;
  grid_string += write_geometry(h5gridfile);
  grid_string += "\n";
  grid_string += write_topology(h5gridfile);

  if(xdmf_format == xdmf_data_formats::hdf5 && !continue_output_after_restart)
  {
#ifdef PARMOON_WITH_HDF5
    herr_t status;
    if (plist_id != H5P_DEFAULT)
    {
      status = H5Pclose(plist_id);
      if (status < 0)
      {
        ErrThrow("Error closing property list");
      }
    }

    if (rank_has_file)
    {
      status = H5Fclose(h5gridfile);
      if (status < 0)
      {
        ErrThrow("Error closing h5 file for the grid");
      }
    }
#endif // #ifdef PARMOON_WITH_HDF5
  }
  return grid_string;
}

template<int d>
std::string DataWriter<d>::write_geometry(hid_t h5gridfile)
{
  auto vertex_coords = this->get_vertex_coord();
  assert(vertex_coords.size() % d == 0);
  local_n_vertices = vertex_coords.size() / d;
  global_n_vertices = local_n_vertices;
  global_vertex_offset = 0;

#ifdef _MPI
  MPI_Allreduce(&local_n_vertices, &global_n_vertices, 1, MPI_UNSIGNED, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Exscan(&local_n_vertices, &global_vertex_offset, 1, MPI_UNSIGNED, MPI_SUM,
             MPI_COMM_WORLD);
#endif // MPI
  // Header for the Geometry element of the XDMF file.
  std::string geostring = (d == 3 ? "Z" : "");
  std::string geo_string = "    <Geometry GeometryType=\"XY" + geostring + "\">\n"
    "      <DataItem Format=\"" + as_string(xdmf_format) +
    "\" DataType=\"Float\" Dimensions=\"" + std::to_string(global_n_vertices) + " "
    + std::to_string(d) + "\" Precision=\"8\">\n";

  // Write the array of stored data to .h5, .bin or directly to .xdmf files.
  switch(xdmf_format)
  {
    case xdmf_data_formats::hdf5:
    {
#ifdef PARMOON_WITH_HDF5
      if (!continue_output_after_restart)
      {
#ifdef _MPI
        if (collective_h5)
        {
#endif
          write_h5_file<2>(d, local_n_vertices, global_n_vertices,
                           global_vertex_offset, compress_h5, h5gridfile,
                           "Geometry", vertex_coords.data());
#ifdef _MPI
        }
        else
        {
          // gather and write vertex coordinates at root

          int mpi_rank, mpi_size;
          MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
          MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

          if (mpi_rank == 0)
          {
            std::vector<int> n_vertices(mpi_size);
            std::vector<double> vertex_coords_all(d * global_n_vertices);

            std::vector<int> n_coords(mpi_size);
            std::vector<int> n_displs(mpi_size);

            MPI_Gather(&local_n_vertices, 1, MPI_UNSIGNED,
              n_vertices.data(), 1, MPI_UNSIGNED,
              0, MPI_COMM_WORLD);

            int displ = 0;
            for (int i = 0; i < mpi_size; i++)
            {
              n_coords[i] = n_vertices[i] * d;
              n_displs[i] = displ;
              displ += n_coords[i];
            }

            MPI_Gatherv(vertex_coords.data(), vertex_coords.size(), MPI_DOUBLE,
              vertex_coords_all.data(), n_coords.data(), n_displs.data(), MPI_DOUBLE,
              0, MPI_COMM_WORLD);

            write_h5_file_noncollective<2>(d, global_n_vertices,
              compress_h5, h5gridfile,
              "Geometry", vertex_coords_all.data());
          }
          else
          {
            MPI_Gather(&local_n_vertices, 1, MPI_UNSIGNED,
              nullptr, 1, MPI_UNSIGNED,
              0, MPI_COMM_WORLD);

            MPI_Gatherv(vertex_coords.data(), vertex_coords.size(), MPI_DOUBLE,
              nullptr, nullptr, nullptr, MPI_DOUBLE,
              0, MPI_COMM_WORLD);
          }

          MPI_Barrier(MPI_COMM_WORLD);
        }
#endif
      }

      geo_string += "        " + base_name + "_grid.h5:/Geometry\n";
#else
      (void) h5gridfile; // suppress compiler warning
#endif // #ifdef PARMOON_WITH_HDF5
      break;
    }
    case xdmf_data_formats::binary:
    {
      if(!continue_output_after_restart)
      {
        std::ofstream g(output_dir + "/" + base_name + "_geometry.bin",
                        std::ios::binary);
        g.write((char*) &vertex_coords[0],
                sizeof(double) * vertex_coords.size());
      }
      geo_string += "        " + base_name + "_geometry.bin\n";
      break;
    }
    case xdmf_data_formats::xml:
    {
      for(unsigned int k = 0; k < local_n_vertices; k++)
      {
        geo_string += convert_to_string(vertex_coords[k*d]) + " "
                    + convert_to_string(vertex_coords[k*d+1]) + " "
                    + (d == 3 ? convert_to_string(vertex_coords[k*d+2]) : "")
                    + "\n";
      }
      break;
    }
  }
  geo_string += "      </DataItem>\n"
                "    </Geometry>\n";
  return geo_string;
}

template<int d>
std::string DataWriter<d>::write_topology(hid_t h5gridfile)
{
  auto topology = get_topology(local_n_vertices, global_vertex_offset);

  std::array<unsigned int, 2> local_n_data =
     {static_cast<unsigned int>(coll->GetN_Cells()),
      static_cast<unsigned int>(topology.size())};
  std::array<unsigned int, 2> global_n_data = local_n_data;
  std::array<unsigned int, 2> global_data_offset = {0, 0};
#ifdef _MPI
  local_n_data[0] = coll->GetN_OwnCells();
  MPI_Allreduce(local_n_data.data(), global_n_data.data(), 2, MPI_UNSIGNED,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Exscan(local_n_data.data(), global_data_offset.data(), 2, MPI_UNSIGNED,
             MPI_SUM, MPI_COMM_WORLD);
#endif // MPI

  // Header for the Topology element of the XDMF file.
  std::string topo_string = "    <Topology TopologyType=\"Mixed\""
    " NumberOfElements=\"" + std::to_string(global_n_data[0]) + "\">\n"
    "      <DataItem Format=\"" + as_string(xdmf_format) + "\" DataType=\"Int\""
    " Dimensions=\""+ std::to_string(global_n_data[1]) + "\">\n";

  // Write the array of stored data to .h5, .bin or directly to .xdmf files.
  switch(xdmf_format)
  {
    case xdmf_data_formats::hdf5:
    {
#ifdef PARMOON_WITH_HDF5
      if(!continue_output_after_restart)
      {
#ifdef _MPI
        if (collective_h5)
        {
#endif
          write_h5_file<1>(d, local_n_data[1], global_n_data[1], global_data_offset[1],
                           compress_h5, h5gridfile, "Topology", topology.data());
#ifdef _MPI
        }
        else
        {
          // gather and write tetrahedra at root

          int mpi_rank, mpi_size;
          MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
          MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

          if (mpi_rank == 0)
          {
            std::vector<int> n_data(mpi_size);
            std::vector<unsigned int> topology_all(global_n_data[1]);

            std::vector<int> n_displs(mpi_size);

            MPI_Gather(&(local_n_data[1]), 1, MPI_UNSIGNED,
              n_data.data(), 1, MPI_UNSIGNED,
              0, MPI_COMM_WORLD);

            int displ = 0;
            for (int i = 0; i < mpi_size; i++)
            {
              n_displs[i] = displ;
              displ += n_data[i];
            }

            MPI_Gatherv(topology.data(), topology.size(), MPI_UNSIGNED,
              topology_all.data(), n_data.data(), n_displs.data(), MPI_UNSIGNED,
              0, MPI_COMM_WORLD);

            write_h5_file_noncollective<1>(d, topology_all.size(),
              compress_h5, h5gridfile,
              "Topology", topology_all.data());
          }
          else
          {
            MPI_Gather(&(local_n_data[1]), 1, MPI_UNSIGNED,
              nullptr, 1, MPI_UNSIGNED,
              0, MPI_COMM_WORLD);

            MPI_Gatherv(topology.data(), topology.size(), MPI_UNSIGNED,
              nullptr, nullptr, nullptr, MPI_UNSIGNED,
              0, MPI_COMM_WORLD);
          }

          MPI_Barrier(MPI_COMM_WORLD);
        }
#endif
      }
      topo_string += "        " + base_name + "_grid.h5:/Topology\n";
#else
      (void) h5gridfile; // suppress compiler warning
#endif // #ifdef PARMOON_WITH_HDF5
      break;
    }
    case xdmf_data_formats::binary:
    {
      if(!continue_output_after_restart)
      {
        std::ofstream t(output_dir + "/" + base_name + "_topology.bin",
                        std::ios::binary);
        t.write((char*) &topology[0], sizeof(unsigned int) * local_n_data[1]);
      }
      topo_string += "        " + base_name + "_topology.bin\n";
      break;
    }
    case xdmf_data_formats::xml:
    {
      for(unsigned int icell = 0, loc_vert = 0; icell < global_n_data[0]; icell++)
      {
        const TBaseCell* cell = coll->GetCell(icell);
        unsigned int n_Loc_vert = cell->GetN_Vertices();
        std::string topo_cell = std::to_string(topology[loc_vert]);
        for(unsigned int ivert = 0; ivert < n_Loc_vert; ++ivert)
        {
          topo_cell += " " + std::to_string(topology[loc_vert + ivert + 1]);
        }
        topo_string += topo_cell + "\n";
        loc_vert += n_Loc_vert+1;
      }
      break;
    }
  }

  topo_string += "      </DataItem>\n"
                 "    </Topology>\n";
#if defined _MPI && defined PARMOON_WITH_HDF5
  if (write_subdomain_info)
  {
    // write subdomain information also to the xdmf file, this is done only once
    subdomain_info_xml = "        <DataItem Name=\"subdomain_data\" Format=\""
        + as_string(xdmf_format) + "\" Dimensions=\""
        + std::to_string(global_n_data[0]) + "\">\n";

    if (!continue_output_after_restart)
    {
      int mpi_rank, mpi_size;
      MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
      MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

      if (collective_h5)
      {
        std::vector<unsigned int> rank_data(local_n_data[0], mpi_rank);

        write_h5_file<1>(d, local_n_data[0], global_n_data[0],
                         global_data_offset[0], compress_h5, h5gridfile,
                         "Subdomain", rank_data.data());
      }
      else
      {
        // gather and write subdomain data at root

        if (mpi_rank == 0)
        {
          std::vector<int> n_data(mpi_size);
          std::vector<unsigned int> rank_data(global_n_data[0]);

          std::vector<int> n_displs(mpi_size);

          MPI_Gather(&(local_n_data[0]), 1, MPI_UNSIGNED,
            n_data.data(), 1, MPI_UNSIGNED,
            0, MPI_COMM_WORLD);

          int displ = 0;
          for (int i = 0; i < mpi_size; i++)
          {
            for (int j = 0; j < n_data[i]; j++)
            {
              rank_data[displ + j] = i;
            }

            n_displs[i] = displ;
            displ += n_data[i];
          }

          write_h5_file_noncollective<1>(d, rank_data.size(),
              compress_h5, h5gridfile,
              "Subdomain", rank_data.data());
        }
        else
        {
          MPI_Gather(&(local_n_data[0]), 1, MPI_UNSIGNED,
            nullptr, 1, MPI_UNSIGNED,
            0, MPI_COMM_WORLD);
        }

        MPI_Barrier(MPI_COMM_WORLD);
      }
    }

    subdomain_info_xml += "          " + base_name + "_grid.h5:/Subdomain\n";
    subdomain_info_xml += "        </DataItem>\n";
  }
#endif // MPI && PARMOON_WITH_HDF5
  return topo_string;
}

template<int d>
std::string DataWriter<d>::write_solution(bool only_xml) const
{
  bool rank_has_file = true;
  bool rank_has_global_file = true;
  bool separate_global_file = false;

#ifdef _MPI
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  rank_has_file = mpi_rank == 0 || collective_h5;

  separate_global_file = collective_h5;
  rank_has_global_file = mpi_rank == 0 && separate_global_file && global_output_variable_funcs.size() > 0;
#endif

  hid_t h5file = 0;
  std::string h5file_name;

  hid_t h5file_global = 0;
  std::string h5file_name_global;

#ifdef PARMOON_WITH_HDF5
  hid_t plist_id = H5P_DEFAULT;
#endif // #ifdef PARMOON_WITH_HDF5

  if (xdmf_format == xdmf_data_formats::hdf5)
  {
#ifdef PARMOON_WITH_HDF5
    const unsigned int n_timesteps = timeValues.size();
    bool first_call = (n_timesteps <= 1);

    std::string time_string;

    if (separate_h5 && n_timesteps > 0)
    {
      time_string += "_" + std::to_string(n_timesteps-1);
    }
    else if (continue_output_after_restart)
    {
      first_call = false; // only open solution file from previous simulation
    }

    h5file_name = base_name + "_solution" + time_string + ".h5";

#ifdef _MPI
    if (separate_global_file)
    {
      h5file_name_global = base_name + "_globals" + time_string + ".h5";
    }
    else
    {
      h5file_name_global = h5file_name;
    }
#else
    h5file_name_global = h5file_name;
#endif

    std::string h5filepath = output_dir + "/" + h5file_name;

    if (!only_xml)
    {
#ifdef _MPI
      if (collective_h5)
      {
        plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
      }
#endif // MPI

      if (rank_has_file)
      {
        if (first_call || separate_h5)
        {
          Output::print<4>("creating h5 file ", h5filepath);
          h5file = H5Fcreate(h5filepath.c_str(), H5F_ACC_EXCL, H5P_DEFAULT,
                             plist_id);
        }
        else
        {
          h5file = H5Fopen(h5filepath.c_str(), H5F_ACC_RDWR, plist_id);
        }

        if (h5file < 0)
        {
          ErrThrow("Error opening the h5 file ", h5filepath);
        }
      }

      if (rank_has_global_file && separate_global_file)
      {
        std::string h5filepath_global = output_dir + "/" + h5file_name_global;

        if (first_call || separate_h5)
        {
          Output::print<4>("creating h5 file ", h5filepath_global);
          h5file_global = H5Fcreate(h5filepath_global.c_str(),
            H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
        }
        else
        {
          h5file_global = H5Fopen(h5filepath_global.c_str(),
            H5F_ACC_RDWR, H5P_DEFAULT);
        }

        if (h5file_global < 0)
        {
          ErrThrow("Error opening the h5 file ", h5filepath_global);
        }
      }
      else
      {
        h5file_global = h5file;
      }
    }
#else
    (void) only_xml; // suppress compiler warning
#endif // PARMOON_WITH_HDF5
  }

  std::string sol_string;

  if (globals_xdmf)
  {
    sol_string += write_global_solution(h5file_global, h5file_name_global, only_xml);
  }

  sol_string += write_scalar_solution(h5file, h5file_name, only_xml);
  sol_string += write_vector_solution(h5file, h5file_name, only_xml);

#ifdef _MPI
  if(write_subdomain_info)
  {
    sol_string += "      <Attribute AttributeType=\"Scalar\" "
                  "Name=\"subdomain\" Center=\"Cell\">\n" + subdomain_info_xml +
                  "      </Attribute>";
  }
#endif

  if(xdmf_format == xdmf_data_formats::hdf5 && !only_xml)
  {
#ifdef PARMOON_WITH_HDF5
    herr_t status;
    if (plist_id != H5P_DEFAULT)
    {
      status = H5Pclose(plist_id);
      if (status < 0)
      {
        ErrThrow("Error closing property list");
      }
    }

    if (rank_has_file)
    {
      status = H5Fclose(h5file);
      if (status < 0)
      {
        ErrThrow("Error closing the h5 file for the solution");
      }
    }

    if (rank_has_global_file && separate_global_file)
    {
      status = H5Fclose(h5file_global);
      if (status < 0)
      {
        ErrThrow("Error closing the h5 file for the globals");
      }
    }
#endif // PARMOON_WITH_HDF5
  }

  return sol_string;
}

template<int d>
std::string DataWriter<d>::write_global_solution(hid_t& h5file,
                                                 std::string h5_file_name,
                                                 bool only_xml)
const
{
  std::string global_string;

  for (unsigned int i = 0; i < global_output_variable_funcs.size(); ++i)
  {
    double value = 0.0;

    if (global_output_variable_funcs[i] != nullptr)
    {
      value = (global_output_variable_funcs[i])();

      if (std::isnan(value))
      {
        Output::root_info("DataWriter", "Value of global output variable '", global_output_variable_names[i], "' is NaN; resetting to zero.");
        value = 0.0;
      }
    }

    global_string += "      <Attribute Center=\"Grid\" AttributeType=\"Scalar\" Name=\""
    + global_output_variable_names[i] + "\">\n        <DataItem Format=\""
     + as_string(xdmf_format) + "\" Dimensions=\"1\" Precision=\"8\">\n";

    switch(xdmf_format)
    {
      case xdmf_data_formats::hdf5:
      {
#ifdef PARMOON_WITH_HDF5
        std::string ds_name = "DS_Global_" + std::to_string(i);
        if (timeValues.size() >= 1 && !separate_h5)
        {
          ds_name += "_";
          ds_name += num2str(timeValues.size() - 1, 5);
        }

        if (!only_xml)
        {
          int mpi_rank = 0;

#ifdef _MPI
          MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

          if (mpi_rank == 0)
          {
            // global data is written to the HDF5 file only at root

            // note that we get a separate individual file just for root
            // if solution writing is collective

            write_h5_file_noncollective<1>(1, 1,
                  false, h5file,
                  ds_name.c_str(), &value);
          }
        }

        global_string += "          " + h5_file_name + ":/" + ds_name + "\n";
#else
       (void) h5file; // suppress compiler warning
       (void) h5_file_name; // suppress compiler warning
#endif // PARMOON_WITH_HDF5
        break;
      }
      case xdmf_data_formats::xml:
      {
        global_string += convert_to_string(value) + "\n";
        break;
      }
      case xdmf_data_formats::binary:
      {
        ErrThrow("Global output variables are currently not supported in XDMF/binary mode.");
        break;
      }
    }

    global_string += "        </DataItem>\n"
                     "      </Attribute>\n";
  }

  return global_string;
}

template<int d>
std::string DataWriter<d>::write_scalar_solution(hid_t& h5file,
                                                 std::string h5_file_name,
                                                 bool only_xml)
const
{
  std::string scalar_string;
  for(unsigned int i = 0; i < fe_functions.size(); ++i)
  {
    const auto* fefunction = fe_functions[i];
    std::vector<double> solutionAtNode;
    unsigned int dimension;
    bool scalar;
    if(!only_xml)
    {
      computeNodeValues(fefunction, solutionAtNode, dimension);
      assert(solutionAtNode.size() % dimension == 0);
      // indicate if this fe function is scalar or vector-valued (eg,
      // Raviart-Thomas)
      scalar = (dimension == 1);
    }
    else
    {
      scalar = (1 == fefunction->GetFESpace()->GetBaseVectDim());
    }
    //assert(scalar || (dimension == d));

    // Header for the scalar solution attribute of the XDMF file.
    std::string sca_vec_string= (scalar ? "Scalar" : "Vector");
    scalar_string += "      <Attribute AttributeType=\"" + sca_vec_string
     + "\" Name=\"" + fefunction->GetName() +"\">\n        <DataItem Format=\""
     + as_string(xdmf_format) + "\" Dimensions=\""+ std::to_string(global_n_vertices)
     + (scalar ? "" : " 3") + "\" Precision=\"8\">\n";

    switch(xdmf_format)
    {
      case xdmf_data_formats::hdf5:
      {
#ifdef PARMOON_WITH_HDF5
        std::string sds_name = "DS_Scalar_" + std::to_string(i);
        if(timeValues.size() >= 1 && !separate_h5)
        {
          sds_name += "_";
          sds_name += num2str(timeValues.size() - 1, 5);
        }

        // avoid rewriting initial solution after restart
        if (!only_xml)
        {
#ifdef _MPI
          if (collective_h5)
          {
#endif
            if (scalar)
            {
              write_h5_file<1>(d, local_n_vertices, global_n_vertices,
                               global_vertex_offset, compress_h5, h5file,
                               sds_name.c_str(), solutionAtNode.data());
            }
            else
            {
              write_h5_file<2>(dimension, local_n_vertices, global_n_vertices,
                               global_vertex_offset, compress_h5, h5file,
                               sds_name.c_str(), solutionAtNode.data());
            }
#ifdef _MPI
          }
          else
          {
            // gather and write data at root

            int mpi_rank, mpi_size;
            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
            MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

            if (mpi_rank == 0)
            {
              if (scalar)
              {
                dimension = 1;
              }

              std::vector<int> n_vertices(mpi_size);
              std::vector<double> solution_all(dimension * global_n_vertices);

              std::vector<int> n_values(mpi_size);
              std::vector<int> n_displs(mpi_size);

              MPI_Gather(&local_n_vertices, 1, MPI_UNSIGNED,
                n_vertices.data(), 1, MPI_UNSIGNED,
                0, MPI_COMM_WORLD);

              int displ = 0;
              for (int i = 0; i < mpi_size; i++)
              {
                n_values[i] = n_vertices[i] * dimension;
                n_displs[i] = displ;
                displ += n_values[i];
              }

              MPI_Gatherv(solutionAtNode.data(), solutionAtNode.size(), MPI_DOUBLE,
                solution_all.data(), n_values.data(), n_displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

              if (scalar)
              {
                write_h5_file_noncollective<1>(dimension, global_n_vertices,
                  compress_h5, h5file,
                  sds_name.c_str(), solution_all.data());
              }
              else
              {
                write_h5_file_noncollective<2>(dimension, global_n_vertices,
                  compress_h5, h5file,
                  sds_name.c_str(), solution_all.data());
              }
            }
            else
            {
              MPI_Gather(&local_n_vertices, 1, MPI_UNSIGNED,
                nullptr, 1, MPI_UNSIGNED,
                0, MPI_COMM_WORLD);

              MPI_Gatherv(solutionAtNode.data(), solutionAtNode.size(), MPI_DOUBLE,
                nullptr, nullptr, nullptr, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
            }

            MPI_Barrier(MPI_COMM_WORLD);
          }
#endif
        }
        scalar_string += "          " + h5_file_name + ":/" + sds_name + "\n";
#else
       (void) h5file; // suppress compiler warning
       (void) h5_file_name; // suppress compiler warning
#endif // PARMOON_WITH_HDF5
        break;
      }
      case xdmf_data_formats::binary:
      {
        std::string file_name = base_name + "_scalar_solution"
                                + std::to_string(i);
        if(!timeValues.empty())
          file_name += "_t" + std::to_string(timeValues.size() - 1);
        file_name += ".bin";
        // avoid rewriting initial solution after restart
        if(!only_xml)
        {
          // Create a new file for every solution step
          std::ofstream s(output_dir + "/" + file_name, std::ios::binary);
          s.write((char*) &solutionAtNode[0],
                  sizeof(double)*solutionAtNode.size());
        }
        scalar_string += "          " + file_name + "\n";
        break;
      }
      case xdmf_data_formats::xml:
      {
        if (scalar)
        {
          for (unsigned int k = 0; k < local_n_vertices; k++)
          {
            scalar_string += convert_to_string(solutionAtNode[k]) + "\n";
          }
        }
        else
        {
          for (unsigned int k = 0; k < local_n_vertices; k++)
          {
            scalar_string += convert_to_string(solutionAtNode[dimension*k]);
            for (unsigned int i = 1; i < dimension; ++i)
              scalar_string += " " + convert_to_string(solutionAtNode[dimension*k+i]);
            scalar_string += "\n";
          }
        }
        break;
      }
    }
    scalar_string += "        </DataItem>\n"
                    "      </Attribute>\n";
  }
  return scalar_string;
}

template<int d>
std::string DataWriter<d>::write_vector_solution(hid_t& h5file,
                                                 std::string h5_file_name,
                                                 bool only_xml)
const
{
  std::string vector_string;
  for (unsigned int i = 0; i < fe_vector_functions.size(); i++)
  {
    const auto* fe_vect = fe_vector_functions[i];
    std::vector<double> solutionAtNode;
    unsigned int dimension;
    if (!only_xml)
    {
      computeNodeValues(fe_vect, solutionAtNode, dimension);
      assert(solutionAtNode.size() % dimension == 0);
    }

    // Header for the vector solution attribute of the XDMF file.
    vector_string += "      <Attribute AttributeType=\"Vector\" Name=\""
      + fe_vect->GetName() + "\">\n        <DataItem Format=\"" + as_string(xdmf_format)
      + "\" Dimensions=\""+std::to_string(global_n_vertices)+" 3\" Precision=\"8\">\n";

    switch (xdmf_format)
    {
      case xdmf_data_formats::hdf5:
      {
#ifdef PARMOON_WITH_HDF5
        std::string vds_name = "DS_Vector_" + std::to_string(i);
        if (timeValues.size() >= 1 && !separate_h5)
        {
          vds_name += "_";
          vds_name += num2str(timeValues.size() - 1, 5);
        }

        // avoid rewriting initial solution after restart
        if (!only_xml)
        {
#ifdef _MPI
          if (collective_h5)
          {
#endif
            write_h5_file<2>(3, local_n_vertices, global_n_vertices,
                             global_vertex_offset, compress_h5, h5file,
                             vds_name.c_str(), solutionAtNode.data());
#ifdef _MPI
          }
          else
          {
            // gather and write data at root

            int mpi_rank, mpi_size;
            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
            MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

            if (mpi_rank == 0)
            {
              std::vector<int> n_vertices(mpi_size);
              std::vector<double> solution_all(d * global_n_vertices);

              std::vector<int> n_values(mpi_size);
              std::vector<int> n_displs(mpi_size);

              MPI_Gather(&local_n_vertices, 1, MPI_UNSIGNED,
                n_vertices.data(), 1, MPI_UNSIGNED,
                0, MPI_COMM_WORLD);

              int displ = 0;
              for (int i = 0; i < mpi_size; i++)
              {
                n_values[i] = n_vertices[i] * d;
                n_displs[i] = displ;
                displ += n_values[i];
              }

              MPI_Gatherv(solutionAtNode.data(), solutionAtNode.size(), MPI_DOUBLE,
                solution_all.data(), n_values.data(), n_displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

              write_h5_file_noncollective<2>(d, global_n_vertices,
                compress_h5, h5file,
                vds_name.c_str(), solution_all.data());
            }
            else
            {
              MPI_Gather(&local_n_vertices, 1, MPI_UNSIGNED,
                nullptr, 1, MPI_UNSIGNED,
                0, MPI_COMM_WORLD);

              MPI_Gatherv(solutionAtNode.data(), solutionAtNode.size(), MPI_DOUBLE,
                nullptr, nullptr, nullptr, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
            }

            MPI_Barrier(MPI_COMM_WORLD);
          }
#endif
        }
        vector_string += "          " + h5_file_name + ":/" + vds_name + "\n";
#else
        (void) h5file; // suppress compiler warning
        (void) h5_file_name; // suppress compiler warning
#endif // PARMOON_WITH_HDF5
        break;
      }
      case xdmf_data_formats::binary:
      {
        std::string file_name = base_name + "_vector_solution"
                               + std::to_string(i);
        if(!timeValues.empty())
          file_name += "_t" + std::to_string(timeValues.size() - 1);
        file_name += ".bin";
        // avoid rewriting initial solution after restart
        if(!only_xml)
        {
          // Create a new file for every solution step
          std::ofstream v( output_dir + "/" + file_name, std::ios::binary);
          v.write((char*) &solutionAtNode[0],
                  sizeof(double)*solutionAtNode.size());
        }
        vector_string += "          " + file_name + "\n";
        break;
      }
      case xdmf_data_formats::xml:
      {
        for(unsigned int k = 0; k < local_n_vertices; k++)
        {
          vector_string += convert_to_string(solutionAtNode[k*3]) + " "
              + convert_to_string(solutionAtNode[k*3+1]) + " "
              + convert_to_string(solutionAtNode[k*3+2]) + "\n";
        }
        break;
      }
    }
    vector_string += "        </DataItem>\n"
                     "      </Attribute>\n";
  }
  return vector_string;
}

template<int d>
std::string DataWriter<d>::write_solutions_before_restart()
{
  if(write_xdmf && (timeValues.size() != 1 || timeValues[0] != time_start))
  {
    ErrThrow("DataWriter<d>::write_solutions_before_restart, something is "
             "wrong");
  }
  /// @attention we assume the previous simulation started at t=0
  timeValues.clear();
  n_steps_until_next_output = 0;
  std::string previous_solutions;
  for(int i = 0, n_previous_time_steps = restart_time / time_step_length;
      i <= n_previous_time_steps; ++i)
  {
    if(n_steps_until_next_output != 0)
    {
      n_steps_until_next_output--;
    }
    else
    {
      n_steps_until_next_output = n_steps_per_output - 1;
      timeValues.push_back(i*time_step_length);
      // the solutions (scalar, vector) are always stored in
      // <Grid Name=...> ... </Grid>, here we set a name
      std::string data_name = "data";
      data_name = "TimeStep_" + std::to_string(i);

      previous_solutions += "\n"
          "    <Grid Name=\""+ data_name + "\">\n"
          "      <Topology Reference=\"/Xdmf/Domain/Topology[1]\"/>\n"
          "      <Geometry Reference=\"/Xdmf/Domain/Geometry[1]\"/>";

      previous_solutions += write_solution(true);
      previous_solutions += "    </Grid>\n"; // Name="data_name"
    }
  }
  time_start = 0.;
  return previous_solutions;
}

template<int d>
void DataWriter<d>::write_csv_file(const std::string& file_name)
{
  if (global_output_variable_names.empty())
  {
    // no global outputs, write nothing
    return;
  }

  bool write_file = true;

#ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  write_file = rank == 0;
#endif

  std::ofstream f = write_file ? std::ofstream(file_name, std::ios::app) : std::ofstream();

  if (write_file && !f)
  {
    ErrThrow("Cannot open file for output: ", file_name);
  }

  const bool isTemporal = !timeValues.empty();
  unsigned int n_timesteps = timeValues.size();

  // the following is only true during the first call to this method
  const bool first_call = (n_timesteps <= 1);

  if (first_call)
  {
    if (isTemporal)
    {
      f << "Time";
    }

    for (unsigned int i = 0; i < global_output_variable_names.size(); ++i)
    {
      if (isTemporal || i > 0)
      {
        f << "," << global_output_variable_names[i];
      }
      else
      {
        f << global_output_variable_names[i];
      }
    }

    f << "\n";
  }

  f << std::fixed << std::scientific << std::setprecision(16);

  if (isTemporal)
  {
    f << timeValues.back();
  }

  for (unsigned int i = 0; i < global_output_variable_funcs.size(); ++i)
  {
    double value = 0.0;

    if (global_output_variable_funcs[i] != nullptr)
    {
      value = (global_output_variable_funcs[i])();

      if (std::isnan(value))
      {
        Output::root_info("DataWriter", "Value of global output variable '", global_output_variable_names[i], "' is NaN; resetting to zero.");
        value = 0.0;
      }
    }

    if (isTemporal || i > 0)
    {
      f << ", " << value;
    }
    else
    {
      f << value;
    }
  }

  f << "\n";
}

/**
 * This function writes the xdmf file (light data) and the heavy data, i.e.,
 * the grid information together with the solution. The 'xdmf_buffer' is filled
 * only in this method and in no other method called here.
 *
 * Grid: If the 'xdmf_format' is xml, everything is written directly into the
 * xdmf file. If the 'xdmf_format' is 'binary', separate files for the vertices,
 * topology (cell->vertex - map), and each scalar and vector solution are
 * written. If the 'xdmf_format' is hdf5, the vertices and topology are written
 * into a single file and also all solutions (per time step) are also written
 * into a single file.
 *
 */
template<int d>
void DataWriter<d>::write_xdmf_file(const std::string& file_name)
{
  const bool isTemporal = !timeValues.empty();
  unsigned int n_timesteps = timeValues.size();
  // the following is only true during the first call to this method
  const bool first_call = (n_timesteps <= 1);

  // Fill the buffer with the lines to create a XDMF file. Geometry and
  // topology are created and written once in all cases, while the solutions
  // are created and written to buffer every timestep in case of a
  // time-dependent problem.
  // The following 'if' is true only during the first invocation of this method.
  if(first_call)
  {
    // Header of the XDMF file.
    xdmf_buffer.push_back("<?xml version=\"1.0\" ?>\n"
      "<!-- # This file was generated by ParMooN -->\n"
      "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"3.0\">\n"
      "  <Domain>");
    xdmf_buffer.push_back(write_grid());
    // Save the buffer position to renew the time step values after each step.
    if(isTemporal)
    {
      xdmf_buffer.push_back(
        "    <Grid GridType=\"Collection\" CollectionType=\"Temporal\">");
      buffer_ind = xdmf_buffer.size();
      xdmf_buffer.push_back("\n");
      if(continue_output_after_restart)
      {
        xdmf_buffer.push_back(write_solutions_before_restart());
        n_timesteps = timeValues.size();
      }
    }
  }
  else
  {
    xdmf_buffer.pop_back(); // "  </Domain>\n</Xdmf>"
    xdmf_buffer.pop_back(); // "    </Grid>\n"
  }

  // write solution, skip initial solution after restart, because it was written
  // as the final solution in the previous simulation.
  if(!continue_output_after_restart || !first_call)
  {
    // the solutions (scalar, vector) are always stored in
    // <Grid Name=...> ... </Grid>, here we set a name
    std::string data_name = "data";
    if(isTemporal)
      data_name = "TimeStep_" + std::to_string(n_timesteps-1);

    xdmf_buffer.push_back("\n"
        "    <Grid Name=\""+ data_name + "\">\n"
        "      <Topology Reference=\"/Xdmf/Domain/Topology[1]\"/>\n"
        "      <Geometry Reference=\"/Xdmf/Domain/Geometry[1]\"/>");

    xdmf_buffer.push_back(write_solution());
    xdmf_buffer.push_back("    </Grid>"); // Name="data_name"
  }

  if(isTemporal)
  {
    // Rewrite old time values and the time dimension once every timestep.
    xdmf_buffer[buffer_ind] = get_time_element(xdmf_timetype, timeValues,
                                               time_start, time_step_length,
                                               n_steps_per_output);
    xdmf_buffer.push_back("    </Grid>"); // CollectionType="Temporal"
  }
  xdmf_buffer.push_back("  </Domain>\n</Xdmf>");
#ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0)
#endif
  write_buffer_to_xdmf(file_name, xdmf_buffer);
}


// TODO: alter Write_ParVTK
/** write stored PARALLEL data into a pvtu and vtu files (XML files for
 *paraview)
 *(Sashikumaar Ganesan) */
#ifdef _MPI

template <int d>
void DataWriter<d>::Write_ParVTK(int img, std::string directory,
                                 std::string basename)
{
  int rank, size;
  double xi = 0, eta = 0, zeta = 0;
  double *WArray = nullptr, *DoubleArray = nullptr;
  double BFValues[MaxN_BaseFunctions3D];
  double* Coords = nullptr;
  static double HexaCoords[] = {-1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1,
                                -1, -1, 1,  1, -1, 1,  1, 1, 1,  -1, 1, 1};
  static double TetraCoords[] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};

  int N_ScalarVar = fe_functions.size();
  int N_VectorVar = fe_vector_functions.size();

  auto comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  std::string vtu("VTU/");
  std::string vtudir = directory + std::string("/") + vtu;

  time_t rawtime;
  struct tm* timeinfo;

  std::vector<TVertex*> Vertices;
  TVertex *Last, *Current;

  if(rank == 0)
  {
    // create the folder to store SubDomain vtu files
    mkdir(vtudir.c_str(), 0777);
  }
  MPI_Barrier(comm);

  time(&rawtime);
  timeinfo = localtime(&rawtime);


  // write the master pvtu file
  if(rank == 0)
  {
    std::ostringstream os;
    os << directory << "/" << basename;
    os << "." << num2str(img, 5) << ".pvtu";
    Output::info("DataWriter", "writing output into file: ", os.str());

    std::ofstream dat(os.str().c_str());
    if(!dat)
    {
      cerr << "cannot open file for output\n";
      MPI_Abort(MPI_COMM_WORLD, 0);
    }

    dat << "<?xml version=\"1.0\"?>\n";
    dat << "\n";
    dat << "<!--\n";
    dat << "      Title: Master file for parallel vtk data\n";
    dat << "    Program: ParMooN\n";
    dat << "    Version: v1.0.0\n";
    dat << "Date & Time: " << asctime(timeinfo) << "\n";
    if(!timeValues.empty())
      dat << "Problem Current Time " << timeValues.back() << "\n";
    dat << "  -->\n";
    dat << "\n";


    dat << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1"
        << "\" byte_order=\"LittleEndian\">\n";
    dat << "<PUnstructuredGrid GhostLevel=\"0\">\n";
    dat << "\n";

    dat << " <PPoints>\n";
    dat << "   <PDataArray type=\"Float32\" Name=\""
        << "Position\" NumberOfComponents=\"3\"/>\n";
    dat << "</PPoints>\n";
    dat << "\n";

    dat << "<PCells>\n";
    dat << "  <PDataArray type=\"Int32\" Name=\"connectivity\""
        << " NumberOfComponents=\"1\"/>\n";
    dat << "  <PDataArray type=\"Int32\" Name=\"offsets\""
        << "      NumberOfComponents=\"1\"/>\n";
    dat << "  <PDataArray type=\"UInt8\" Name=\"types\""
        << "        NumberOfComponents=\"1\"/>\n";
    dat << "</PCells>\n";
    dat << "\n";

    dat << "<PPointData Vectors=\"Vectors\" Scalars=\"Scalars\">\n";

    for(int i = 0; i < N_VectorVar; i++)
    {
      dat << "  <PDataArray type=\"Float32\" Name=\""
          << fe_vector_functions[i]->GetName() << "\""
          << " NumberOfComponents=\"3\" format=\"ascii\"/>\n";
      for(int j = 0; j < fe_vector_functions[i]->GetN_Components(); j++)
      {
        dat << "  <DataArray type=\"Float32\" Name=\""
            << fe_vector_functions[i]->GetName() << j << "\" NumberOfComponents=\""
            << "1\" format=\"ascii\"/>\n";
      }
      dat << "\n";
    }

    for(int i = 0; i < N_ScalarVar; i++)
      dat << "  <PDataArray type=\"Float32\" Name=\""
          << fe_functions[i]->GetName() << "\""
          << " NumberOfComponents=\"1\" format=\"ascii\"/>\n";
    dat << "</PPointData>\n";
    dat << "\n";

    dat << "<PCellData Scalars=\"SubDomainAndRegionID\">\n";
    dat << "  <PDataArray type=\"Int32\"   Name=\"SubDomain\""
        << "  NumberOfComponents=\"1\"/>\n";
    dat << "  <PDataArray type=\"Int32\"   Name=\"RegionID\""
        << "  NumberOfComponents=\"1\"/>\n";

    dat << "  <PDataArray type=\"Int32\"   Name=\"CellID\""
        << "  NumberOfComponents=\"1\"/>\n";
    dat << "  <PDataArray type=\"Int32\"   Name=\"LocalCellID\""
        << "  NumberOfComponents=\"1\"/>\n";

    dat << "</PCellData>\n";
    dat << "\n";


    // int begin = 1; // root does not take part in computation
    int begin = 0; // root take part in computation

    for(int i = begin; i < size; i++)
    {
      dat << "  <Piece Source=\"" << vtu << basename << "."
          << num2str(i, 4) << "." << num2str(img, 5) << ".pvtu\"/>\n";
    }
    dat << "\n";

    dat << "</PUnstructuredGrid>\n";
    dat << "</VTKFile>\n";
    dat.close();
  } // if(rank==0)

  int N_Elements = coll->GetN_OwnCells(); // possibly zero
  int N_LocVertices = 0;
  for(int i = 0; i < N_Elements; i++)
  {
    auto cell = coll->GetCell(i);
    N_LocVertices += cell->GetN_Vertices();
  }

  if(N_LocVertices)
    Vertices.resize(N_LocVertices);
  for(int i = 0, N_ = 0; i < N_Elements; i++)
  {
    auto cell = coll->GetCell(i);
    int k = cell->GetN_Vertices();
    for(int j = 0; j < k; j++)
    {
      Vertices[N_] = cell->GetVertex(j);
      N_++;
    }
  }

  if(!Vertices.empty())
    std::sort(Vertices.begin(), Vertices.end());

  Last = nullptr;
  int N_Vertices = 0;
  for(int i = 0; i < N_LocVertices; i++)
  {
    TVertex *Current = Vertices[i];
    if(Current != Last)
    {
      N_Vertices++;
      Last = Current;
    }
  }

  int *VertexNumbers = nullptr, *NumberVertex = nullptr;
  if(N_LocVertices)
  {
    Coords = new double[3 * N_Vertices];
    VertexNumbers = new int[N_LocVertices];
    NumberVertex = new int[N_LocVertices];
  }
  Last = nullptr;
  for(int i = 0, N_ = 0, k = -1; i < N_LocVertices; i++)
  {
    TVertex *Current = Vertices[i];
    if(Current != Last)
    {
      Vertices[i]->GetCoords(Coords[3 * N_], Coords[3 * N_ + 1],
                             Coords[3 * N_ + 2]);
      k++;
      N_++;
      Last = Current;
    }
    NumberVertex[i] = k;
  }

  for(int i = 0, m = 0; i < N_Elements; i++)
  {
    auto cell = coll->GetCell(i);
    int k = cell->GetN_Vertices();
    for(int j = 0; j < k; j++)
    {
      Current = cell->GetVertex(j);
      int l = std::distance(
          Vertices.begin(),
          std::lower_bound(Vertices.begin(), Vertices.end(), Current));
      VertexNumbers[m] = NumberVertex[l];
      m++;
    } // endfor j
  }   // endfor i

  if(N_LocVertices)
  {
    DoubleArray = new double[3 * N_Vertices];
    WArray = new double[N_Vertices];
  }

  std::ostringstream os;
  os << vtudir << "/" << basename;
  os << "." << num2str(rank, 4) << "." << num2str(img, 5) << ".pvtu" << ends;

  std::ofstream dat(os.str().c_str());
  if(!dat)
  {
    Output::warn("DataWriter<d>::Write_ParVTK", "cannot open file for output");
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  dat << setprecision(8);

  dat << "<?xml version=\"1.0\"?>\n";
  dat << "\n";
  dat << "<!--\n";
  dat << "      Title: SubDomain data for master ptvu file\n";
  dat << "    Program: ParMooN\n";
  dat << "    Version: v1.0.0\n";
  dat << "Date & Time: " << asctime(timeinfo) << "\n";
  dat << "  -->\n";
  dat << "\n";

  dat << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1"
      << "\" byte_order=\"LittleEndian\">\n";
  dat << "<UnstructuredGrid>\n";
  dat << "\n";

  dat << "<Piece NumberOfPoints=\"" << N_Vertices << "\" NumberOfCells=\""
      << N_Elements << "\">\n";
  dat << "<Points>\n";
  dat << "  <DataArray type=\"Float32\" Name=\"Position\""
      << " NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(int i = 0, N_ = 0; i < N_Vertices; i++)
  {
    dat << "   " << Coords[N_] << " " << Coords[N_ + 1] << " " << Coords[N_ + 2]
        << "\n";
    N_ += 3;
  }
  dat << "  </DataArray>\n";
  dat << "</Points>\n";


  dat << "<Cells>\n";
  dat << "  <DataArray type=\"Int32\" Name=\"connectivity\""
      << " NumberOfComponents=\"1\" format=\"ascii\">\n";

  for(int i = 0, l = 0; i < N_Elements; i++)
  {
    int N_CellVertices = coll->GetCell(i)->GetN_Vertices();
    dat << "       ";
    for(int j = 0; j < N_CellVertices; j++)
    {
      dat << VertexNumbers[l] << " ";
      l++;
    }
    dat << "\n";
  }
  dat << "  </DataArray>\n";
  dat << "\n";

  dat << "  <DataArray type=\"Int32\" Name=\"offsets\""
      << " NumberOfComponents=\"1\" format=\"ascii\">";
  for(int i = 1; i <= N_Elements; i++)
  {
    if(i%10 == 1) dat << "\n";
    int N_CellVertices = coll->GetCell(i - 1)->GetN_Vertices();
    dat << i * N_CellVertices << "  ";
  }
  dat << "\n   </DataArray>\n";
  dat << "\n";

  dat << "  <DataArray type=\"UInt8\"  Name=\"types\""
      << " NumberOfComponents=\"1\" format=\"ascii\">";
  for(int i = 0; i < N_Elements; i++)
  {
    if(i%10 == 0) dat << "\n";
    int N_CellVertices = coll->GetCell(i)->GetN_Vertices();
    dat << n_local_vertices_to_type<d>(N_CellVertices) << " ";
  }
  dat << "\n  </DataArray>\n";
  dat << "</Cells>\n";
  dat << "\n";
  dat << "<PointData Vectors=\"Velocity\" Scalars=\"Scalars\">\n";
  dat << "\n";

  // write vector variables into file
  if(N_LocVertices)
  {
    for(int k = 0; k < N_VectorVar; k++)
    {
      auto fespace = fe_vector_functions[k]->GetFESpace3D();
      int N_Comp = fe_vector_functions[k]->GetN_Components();
      int Length = fe_vector_functions[k]->GetLength();
      const double* Coeffs = fe_vector_functions[k]->GetValues();

      memset(DoubleArray, 0, sizeof(double) * N_Vertices * N_Comp);
      memset(WArray, 0, sizeof(double) * N_Vertices);

      for(int i = 0, m = 0; i < N_Elements; i++)
      {
        auto cell = coll->GetCell(i);
        int N_ = cell->GetN_Vertices();

        // find FE data for this element
        const BaseFunctions* bf = fespace->get_fe(i).GetBaseFunct();
        const int * DOF = fespace->GetGlobalDOF(i);
        int N_LocDOF = bf->GetDimension();
        for(int j = 0; j < N_; j++)
        {
          switch(cell->GetN_Vertices())
          {
            case 4:
              xi = TetraCoords[3 * j];
              eta = TetraCoords[3 * j + 1];
              zeta = TetraCoords[3 * j + 2];
              break;

            case 8:
              xi = HexaCoords[3 * j];
              eta = HexaCoords[3 * j + 1];
              zeta = HexaCoords[3 * j + 2];
              break;
          }
          bf->GetDerivatives(MultiIndex3D::D000, xi, eta, zeta, BFValues);
          bf->ChangeBF(coll, cell, BFValues);

          for(int n = 0; n < N_Comp; n++)
          {
            double value = 0;
            for(int l = 0; l < N_LocDOF; l++)
              value += BFValues[l] * Coeffs[DOF[l] + n * Length];
            DoubleArray[N_Comp * VertexNumbers[m] + n] += value;
          }
          WArray[VertexNumbers[m]] += 1.;
          m++;
        } // endfor j
      }   // endfor i

      // midle value
      for(int i = 0, l = 0; i < N_Vertices; i++)
      {
        for(int j = 0; j < N_Comp; j++)
        {
          if(WArray[i] > 1.)
            DoubleArray[l] /= WArray[i];
          l++;
        }
      } // endfor l

      //      for(i=0;i<2*N_Vertices;i++)
      //       cout << "Do[" << i << "]" << DoubleArray[i] << "\n";

      dat << "  <DataArray type=\"Float32\" Name=\""
          << fe_vector_functions[k]->GetName() << "\" NumberOfComponents=\""
          << "3\" format=\"ascii\">\n";
      for(int i = 0; i < N_Vertices; i++)
      {
        for(int j = 0; j < N_Comp; j++)
          dat << DoubleArray[N_Comp * i + j] << " ";
        dat << "\n";
      }
      dat << "  </DataArray>\n";
      dat << "\n";

      for(int j = 0; j < N_Comp; j++)
      {
        dat << "  <DataArray type=\"Float32\" Name=\""
            << fe_vector_functions[k]->GetName() << j << "\" NumberOfComponents=\""
            << "1\" format=\"ascii\">";
        for(int i = 0; i < N_Vertices; i++)
        {
          if(i%10 == 0) dat << "\n";
          dat << DoubleArray[i * N_Comp + j] << " ";
        }
        dat << "\n  </DataArray>\n";
        dat << "\n";
      }
    } // for(k=0;k<N_Vec
  }

  // write scalar variables into file
  if(N_LocVertices)
  {
    for(int k = 0; k < N_ScalarVar; k++)
    {
      auto fespace = fe_functions[k]->GetFESpace3D();
      const double* Coeffs = fe_functions[k]->GetValues();

      memset(DoubleArray, 0, sizeof(double) * N_Vertices);
      memset(WArray, 0, sizeof(double) * N_Vertices);

      for(int i = 0, m = 0; i < N_Elements; i++)
      {
        auto cell = coll->GetCell(i);
        int N_ = cell->GetN_Vertices();

        // find FE data for this element
        const BaseFunctions* bf = fespace->get_fe(i).GetBaseFunct();
        const int * DOF = fespace->GetGlobalDOF(i);
        int N_LocDOF = bf->GetDimension();
        for(int j = 0; j < N_; j++)
        {
          switch(cell->GetN_Vertices())
          {
            case 4: // Tetrahedron
              xi = TetraCoords[3 * j];
              eta = TetraCoords[3 * j + 1];
              zeta = TetraCoords[3 * j + 2];
              break;
            case 8: // Hexahedron
              xi = HexaCoords[3 * j];
              eta = HexaCoords[3 * j + 1];
              zeta = HexaCoords[3 * j + 2];
              break;
          }
          bf->GetDerivatives(MultiIndex3D::D000, xi, eta, zeta, BFValues);
          bf->ChangeBF(coll, cell, BFValues);
          double value = 0;
          for(int l = 0; l < N_LocDOF; l++)
            value += BFValues[l] * Coeffs[DOF[l]];
          DoubleArray[VertexNumbers[m]] += value;
          WArray[VertexNumbers[m]] += 1.;
          m++;
        } // endfor j
      }   // endfor i

      // non conforming
      for(int i = 0; i < N_Vertices; i++)
      {
        if(WArray[i] > 1.)
          DoubleArray[i] /= WArray[i];
      }

      dat << "  <DataArray type=\"Float32\" Name=\""
          << fe_functions[k]->GetName() << "\" NumberOfComponents=\""
          << "1\" format=\"ascii\">";
      for(int j = 0; j < N_Vertices; j++)
      {
        if(j%10 == 0) dat << "\n";
        dat << DoubleArray[j] << " ";
      }
      dat << "\n  </DataArray>\n";
      dat << "\n";
    } // for(k=0;k<N_ScalarVar
  }

  dat << "</PointData>\n";
  dat << "\n";

  dat << "<CellData Scalars=\"SubDomain\">\n";
  dat << "  <DataArray type=\"Int32\" Name=\"SubDomain\""
      << " NumberOfComponents=\"1\" format=\"ascii\">";
  for(int i = 0; i < N_Elements; i++)
  {
    if(i%10 == 0) dat << "\n";
    dat << (coll->GetCell(i))->GetSubDomainNo() << " ";
  }
  dat << "\n  </DataArray>\n";

  dat << "  <DataArray type=\"Int32\" Name=\"RegionID\""
      << " NumberOfComponents=\"1\" format=\"ascii\">";
  for(int i = 0; i < N_Elements; i++)
  {
    if(i%10 == 0) dat << "\n";
    dat << (coll->GetCell(i))->GetRegionID() << " ";
  }
  dat << "\n  </DataArray>\n";

  dat << "  <DataArray type=\"Int32\" Name=\"CellID\""
      << " NumberOfComponents=\"1\" format=\"ascii\">";
  for(int i = 0; i < N_Elements; i++)
  {
    if(i%10 == 0) dat << "\n";
    dat << (coll->GetCell(i))->GetGlobalCellNo() << " ";
  }
  dat << "\n  </DataArray>\n";

  dat << "  <DataArray type=\"Int32\" Name=\"LocalCellID\""
      << " NumberOfComponents=\"1\" format=\"ascii\">";
  for(int i = 0; i < N_Elements; i++)
  {
    if(i%10 == 0) dat << "\n";
    dat << i << " ";
  }
  dat << "\n  </DataArray>\n";

  dat << "</CellData>\n";

  dat << "</Piece>\n";
  dat << "</UnstructuredGrid>\n";
  dat << "</VTKFile>\n";

  dat.close();

  if(N_LocVertices)
  {
    delete[] NumberVertex;
    delete[] VertexNumbers;
    delete[] DoubleArray;
    delete[] WArray;
    delete[] Coords;
  }
} // DataWriter::Write_ParVTK
#endif // _MPI

template <int d>
void DataWriter<d>::writeCoord(std::ofstream& f)
{
  f.precision(8);
  unsigned int N_Vertices = coll->GetN_Vertices();
  for(unsigned int i = 0; i < N_Vertices; i++)
  {
    ///@attention case output works only with this format
    f << coll->GetCoord(i * d);
    f << " ";
     f << coll->GetCoord(i * d + 1);
    f << " ";
    if(d == 2)
    {
      f << 0.;
    }
    else
    {
      f << coll->GetCoord(i * d + 2);
    }
    f << "\n";
  }
  f << "\n";
}

template <int d>
template <class T>
void DataWriter<d>::computeNodeValues(const T* function,
                                      std::vector<double>& solutionAtNode,
                                      unsigned int& dimension) const
{
  auto FESpace = function->GetFESpace();
  // get function type of T
  bool IsVect = std::is_same<T, TFEVectFunct3D>::value
                || std::is_same<T, TFEVectFunct2D>::value;
  dimension = 1;
  if(IsVect || FESpace->GetBaseVectDim() > 1)
  {
    // writing vector data requires 3 components even in 2D
    dimension = 3;
  }

  solutionAtNode.assign(dimension * local_n_vertices, 0.0);
  std::vector<int> WArray(local_n_vertices, 0);
  unsigned int n_cells = coll->GetN_Cells();

  for(unsigned int icell = 0, loc_vert = 0; icell < n_cells; icell++)
  {
    const TBaseCell* cell = coll->GetCell(icell);
#ifdef _MPI
    if(cell->IsHaloCell())
      continue;
#endif // MPI
    unsigned int nLocalVertices = cell->GetN_Vertices();

    for(unsigned int nvert = 0; nvert < nLocalVertices; nvert++)
    {
      unsigned int globalVert_index = coll->GetGlobalVerNo(icell, nvert);
      if(write_discontinuous_functions)
        globalVert_index = loc_vert + nvert;

      double function_value[3] = {0.};
      double x, y, z;
      cell->GetVertex(nvert)->GetCoords(x, y, z);
#ifdef __2D__
      function->FindValueLocal(cell, icell, x, y, function_value);
#else
      function->FindValueLocal(cell, icell, x, y, z, function_value);
#endif
      for(unsigned int dim = 0; dim < dimension; dim++)
      {
        solutionAtNode[dimension * globalVert_index + dim]
            += function_value[dim];
      }
      WArray.at(globalVert_index) += 1;
    } // endfor nvert
    loc_vert += nLocalVertices;
  }   // endfor icell

  if(!write_discontinuous_functions)
  {
    for(unsigned int nvert = 0; nvert < local_n_vertices; nvert++)
    {
      if(WArray[nvert] != 0)
      {
        for(unsigned int dim = 0; dim < dimension; dim++)
        {
          solutionAtNode[dimension * nvert + dim] /= WArray[nvert];
        }
      }
    }
  }
}

// ****************
// scalar components of vector variables
// ****************
template <int d>
void DataWriter<d>::printVectCompwise(std::ofstream& dat,
                                      const std::string& name,
                                      unsigned int N_Vertices,
                                      unsigned int N_Comp,
                                      const std::vector<double>& solutionAtNode)
{
  for(unsigned int ncomp = 0; ncomp < N_Comp; ncomp++)
  {
    dat << "SCALARS " << name;
    if(N_Comp != 1)
      dat << ncomp;
    dat << " double\n";
    dat << "LOOKUP_TABLE default\n";
    for(unsigned int nvert = 0; nvert < N_Vertices; nvert++)
    {
      double value = solutionAtNode[nvert * N_Comp + ncomp];
      dat << std::scientific << value << "\n";
    }
    dat << "\n\n";
  }
}

// ***************
// absolute value of vector variables
// ***************
template <int d>
void DataWriter<d>::printVectAbsValue(
  std::ofstream& dat, const std::string& name, unsigned int N_Vertices,
  unsigned int N_Comp, const std::vector<double>& solutionAtNode)
{
  dat << "SCALARS |" << name << "|";
  dat << " double\n";
  dat << "LOOKUP_TABLE default\n";
  for(unsigned int nvert = 0, l = 0; nvert < N_Vertices; nvert++)
  {
    double t = 0;
    for(unsigned int ncomp = 0; ncomp < N_Comp; ncomp++)
    {
      t += solutionAtNode[l] * solutionAtNode[l];
      l++;
    }
    dat << std::sqrt(t) << "\n";
  }
  dat << "\n\n";
}

// ***************
// VECTORS
// ***************
template <int d>
void DataWriter<d>::printVectPointwise(
  std::ofstream& dat, const std::string& name, unsigned int N_Vertices,
  unsigned int N_Comp, const std::vector<double>& solutionAtNode)
{
  dat << "VECTORS " << name;
  dat << " double\n";
  for(unsigned int nvert = 0; nvert < N_Vertices; nvert++)
  {
    for(unsigned int ncomp = 0; ncomp < N_Comp; ncomp++)
    {
      double value = solutionAtNode[nvert * N_Comp + ncomp];
      dat << std::scientific << value << " ";
    }
    if(N_Comp == 2)
      dat << double(0);
    dat << "\n";
  }
  dat << "\n";
}

void check_for_old_parameters(const ParameterDatabase & param_db,
                              ParameterDatabase db)
{
  if(param_db.contains("output_write_vtu"))
  {
    Output::warn("DataWriter", "you provided the parameter 'output_write_vtu' "
                 "which is no longer supported. Please switch to xdmf format "
                 "setting 'output_write_xdmf'.");
  }
  if(param_db.contains("output_write_case"))
  {
    Output::warn("DataWriter", "you provided the parameter 'output_write_case' "
                 "which is no longer supported. Please switch to xdmf format "
                 "setting 'output_write_xdmf'.");
  }
  if(param_db.contains("output_vtk_directory"))
  {
    Output::warn("DataWriter", "you provided the parameter "
                 "'output_vtk_directory' which is no longer supported. Please "
                 "set 'output_directory' instead.");
    if(!param_db.contains("output_directory"))
    {
      Parameter p{param_db["output_vtk_directory"]};
      db["output_directory"].set(p.get<std::string>(), false);
    }
  }
}

template <int d>
std::vector<double> DataWriter<d>::get_vertex_coord() const
{
  if(write_discontinuous_functions)
  {
    // counting vertices (multiple times if appearing in multiple cells)
    unsigned int n_vertices = get_n_vertices_multiply_counted(coll);
    std::vector<double> vertex_coords(d*n_vertices);
    unsigned int n_cells = coll->GetN_Cells();
    for(unsigned int icell = 0, n_total_locvert = 0; icell < n_cells; icell++)
    {
      const TBaseCell* cell = coll->GetCell(icell);
#ifdef _MPI
      if(cell->IsHaloCell())
        continue;
#endif // MPI
      unsigned int n_Loc_vert = cell->GetN_Vertices();
      // Compute and store the coordinates as an array for every cell.
      for(unsigned int ivert = 0; ivert < n_Loc_vert; ++ivert)
      {
        const TVertex * v = cell->GetVertex(ivert);
        double x, y, z;
        v->GetCoords(x, y, z);
        vertex_coords[d*n_total_locvert + d*ivert] = x;
        vertex_coords[d*n_total_locvert + d*ivert + 1] = y;
        if(d==3)
        {
          vertex_coords[d*n_total_locvert + d*ivert + 2] = z;
        }
      }
      n_total_locvert += n_Loc_vert;
    }
    return vertex_coords;
  }
  else
    return coll->get_nodes_coords();
}

template <int d>
std::vector<unsigned int> DataWriter<d>::get_topology(unsigned int n_vertices,
                                                      unsigned int offset)
const
{
  if(!write_discontinuous_functions)
  {
    // count the number of vertices, each multiple times. If we write
    // discontinuous functions, this is already computed and passed into this
    // method as an argument, otherwise we have to recalculate it here:
    n_vertices = get_n_vertices_multiply_counted(coll);
  }

  unsigned int n_cells = coll->GetN_Cells();
#ifdef _MPI
  n_cells = coll->GetN_OwnCells();
#endif // MPI
  // for each cell he write: one number indicating the type and the vertices
  std::vector<unsigned int> topology(n_cells + n_vertices);
  // note that in sequential i==icell, but in parallel this it might be i>icell
  for(unsigned int i = 0, icell = 0, loc_vert = 0; i < n_cells; i++)
  {
    const TBaseCell* cell = coll->GetCell(icell);
#ifdef _MPI
    if(cell->IsHaloCell())
      continue;
#endif // MPI
    unsigned int n_Loc_vert = cell->GetN_Vertices();
    int cell_type;
    // Compute the cell types and the connections of the nodes for every cell.
    switch(n_Loc_vert)
    {
      case 3: cell_type = 4; break; // Triangle
      case 8: cell_type = 9; break; // Hexahedron
      case 4: cell_type = (d == 2 ? 5 : 6); break; // Quadrilateral/Tetrahedron
      default:
        ErrThrow("Unknown number of local vertices ",n_Loc_vert, " in ",d, "D.");
        break;
    }
    topology[loc_vert + icell] = cell_type;
    if(write_discontinuous_functions)
    {
      for(unsigned int ivert = 0; ivert < n_Loc_vert; ++ivert)
      {
        topology[loc_vert + icell + ivert + 1] = loc_vert + ivert + offset;
      }
    }
    else
    {
      for(unsigned int ivert = 0; ivert < n_Loc_vert; ++ivert)
      {
        topology[loc_vert + icell + ivert + 1]
          = this->coll->GetGlobalVerNo(icell, ivert) + offset;
      }
    }
    loc_vert += n_Loc_vert;
    icell++;
  }
  return topology;
}

unsigned int get_n_vertices_multiply_counted(const TCollection* coll)
{
  unsigned int n_vertices = 0;
  for(unsigned int i = 0, n_cells = coll->GetN_Cells(); i < n_cells; i++)
  {
    auto cell = coll->GetCell(i);
#ifdef _MPI
    if(cell->IsHaloCell())
      continue;
#endif // MPI
    n_vertices += cell->GetN_Vertices();
  }
  return n_vertices;
}

template <int d>
void DataWriter<d>::retrieve_example_output()
{
  if (has_example_output)
  {
    return;
  }

  for (std::string name: ExampleOutput::ListOutputVariables())
  {
    add_global_output_variable(name, [&, name] () -> double
    {
      return ExampleOutput::GetOutputVariable(name);
    });
  }

  has_example_output = true;
}

#ifdef __3D__
template class DataWriter<3>;
#else
template class DataWriter<2>;
#endif
