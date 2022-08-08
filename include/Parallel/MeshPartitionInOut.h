#ifndef INCLUDE_PARALLEL_MESH_PARTITION_IN_OUT_INCLUDED
#define INCLUDE_PARALLEL_MESH_PARTITION_IN_OUT_INCLUDED

#include "all_defines_external_libraries.h"
#ifdef _MPI
#ifdef PARMOON_WITH_METIS
#include <metis.h>  //for the id_x typedef
#else // PARMOON_WITH_METIS
typedef int32_t idx_t;
#endif // PARMOON_WITH_METIS

//forward declaration
class TDomain;

/**
 * A namespace holding to methods. They are used for reading/writing the
 * information, which cell of a domain belongs to which processor in a parallel
 * run of the program to/from a text file. The respective filename must be known
 * to the database of the domain.
 * Both methods need be called by root only, this is a sequential part of the
 * code (still).
 *
 *
 * TODO: Describe the used file format.
 */
namespace MeshPartitionInOut
{
  /// Read domain partitioning information from a file.
	void read_file(const TDomain& Domain, int size,
	               std::vector<idx_t>& Cell_Rank);

	/// Write domain partitioning information to a file.
	void write_file(const TDomain& Domain, int size,
	                const std::vector<idx_t>& Cell_Rank);
}
#endif // _MPI

#endif //INCLUDE_PARALLEL_MESH_PARTITION_IN_OUT_INCLUDED
