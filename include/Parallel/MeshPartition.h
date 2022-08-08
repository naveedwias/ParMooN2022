// =======================================================================
// @(#)MeshPartition.h
// 
// Purpose:     partition the domain into "npart" parts for parallel computing
// 
// Author:      Sashikumaar Ganesan
// History:      start of implementation  07/09/09 (Sashikumaar Ganesan)
// =======================================================================
#ifndef __MESH_PARTITION__
#define __MESH_PARTITION__


#ifdef  _MPI
#include <Collection.h>

class TDomain;
class TVertex;

/**@brief Calls one of two Metis methods to distribute the mesh among
 * participating processes. Which one is determined by the control parameter
 * metis_type. Setting that to 0 means call of METIS_PartMeshNodal, 1 a call of
 * METIS_PartMeshDual. The method further sets the datastructure of cells,
 * domain and collection to store the information on own and halo cells.
 * Works on finest refinement level only.
 * Can so far deal with homogeneous tetrahedral or hexahedral meshes only.
 *
 * @note One should keep in mind that after calling this function the Domain
 * is reduced to the process' subdomain - no process knows the whole domain anymore.
 *
 * @param[in] comm The used MPI Communicator.
 * @param[in,out] Domain The domain whose finest collection is to be distributed.
 * @param[in] metis_type must be 0 or 1, partitioning is done by either
              METIS_PartMeshNodal (0) or METIS_PartMeshDual (1).
 * @param[out] MaxCpV The globally unique maximum number of cells per vertex.
 *
 * @return 0 if partitioning was successful, 1 if unsuccessful.
 *           At the moment no strong guarantee can be given, because
 *           we don't know enough about the method yet.
 */
int Partition_Mesh3D(MPI_Comm comm, TDomain *Domain, int metis_type, int &MaxCpV);

/**
 * @brief Removes unneeded halo cells which get introduced by refining a process' subdomain
 * (so after the initial mesh has been partitioned).
 *
 * This method changes the Cell Tree of the Domain it works on in such a way,
 * that all remaining cells are set as root cells. This makes the GetCollection method
 * relatively useless, for after the call only one cell level (0) is known to the Domains tree.
 * \todo Rework Domain_Crop as to modify the Domain's cell tree less brutally!
 *
 * @note The method does not physically delete any cells, edges or nodes at the moment.
 * That code was commented out for some reason.
 *
 *  @param[in] comm The MPI Communicator. In ParMooN this is almost always MPI_COMM_WORLD.
 *  @param[in,out] Domain The domain on whose finest level superfluous halo cells have to be removed.
 */
void Domain_Crop(MPI_Comm comm, TDomain *Domain);

/**
 * @brief return all vertices ordered by the address in memory
 */
std::vector<const TVertex*> get_sorted_vertices(const TCollection* coll);

/**
 * @brief return the number of vertices
 */
int get_number_vertex(const std::vector<const TVertex*>& all_vertices,
                      const TCollection* coll, std::vector<int>& VertexNumbers);
#endif // _MPI


#endif // __MESH_PARTITION__

