/** ***************************************************************************
 *
 * @name       TParFEMapper3D
 * @brief      Class containing communication routines between processes which all
 * 			   administer their own subdomain and corresponding degrees of freedom of
 * 			   an FE system.
 *
 *             Objects of this type are set up after the domain is split, so
 *             there exists a particular one in each process. The data members
 *             refer to the data of the processes' Domain and FESpace.
 *
 *			   All the information the class needs is known to an object of class
 *			   TParFEMapper3D.
 *
 * @author     Sashikumaar Ganesan
 * @date       2015/04/24
 *
 * @ruleof0
 *
 ******************************************************************************/

#ifdef _MPI

#ifndef __PARFECOMMUNICATOR3D__
#define __PARFECOMMUNICATOR3D__

#include "mpi.h"

#include <ParFEMapper3D.h>
#include <vector>

class TParFECommunicator3D
{
  public:
    /**
     * Set up a ParFECommunicator belonging to a certain FESpace.
     *
     * TODO This should be replaced by a ctor taking an FESpace - ParFEMapper
     * should only work as the underlying data structure and not be visible from
     * the outside.
     */
    TParFECommunicator3D(TParFEMapper3D *mapper);

    /// Default constructor. Creates an almost empty object, which should not be used.
    TParFECommunicator3D();

    /// return the number of global degrees of freedom (across all processes)
    size_t get_n_global_dof() const;

    /// Gather information about this communicator in root and print it
    /// to console and output file.
    void print_info() const;

    /**
     * Restores the required level of consistency of a vector.
     * For the notion of consistency levels refer to our ParMooN paper:
     * TODO cite the ParMooN paper correctly...
     *
     * Not that this method is not range checked - you must make sure that
     * "vector" whcih is supposed to be updated is as long as N_Dof.
     */
    void consistency_update(double* vector, size_t level) const;

    /// @brief Queues an asynchronous consistency update.
    /// @param vector The vector to make consistent
    /// @param level The consistency level to achieve
    void queue_consistency_update(double* vector, size_t level) const;

    /// @brief Queues an asynchronous additive update, adding up all
    /// entries corresponding to each global DOF across all ranks. The vector
    /// will be left in consistency level 0.
    /// @param The vector to update
    void additive_update(double* vector, bool with_halo = true) const;

    /// @brief Executes a synchronous additive update.
    void queue_additive_update(double* vector, bool with_halo = true) const;

    /// @brief Waits for all outstanding asynchronous consistency updates
    /// to finish and writes received data into the associated vectors.
    static void flush_consistency_updates();

    /**
     * This concerns only interface masters and slaves. For each interface dof,
     * calculate the average of all its values on all processors where it is
     * present and store that on its master value. The vector will leave this
     * method in level 0 consistency.
     *
     * @param[in, out] vector The vector whose interface dofs should be averaged.
     *
     * This is a relatively common technique and therefore we offer a
     * specialized method for it.
     */
    void average_at_interface(double* vector) const;

    /**
     * Bring a vector stored in additive ("inconsistent") format to a certain
     * consistency level.
     *
     * Therefore we "invert" the usual communication channels: all interface
     * slave, halo 1 and halo 2 d.o.f. send their value to the responsible
     * master, which then adds up all received values.
     *
     * With the parameter "consistency_level" one can determine, whether the
     * masters should report their new values to their slaves, and to which
     * group of them.
     *
     * 0 - do not report at all (only masters have correct value)
     * 1 - report to interface slaves only
     * 2 - level 1 plus report to halo1 d.o.f.
     * 3 - level 2 plus report to halo2 d.o.f. (full consistency)
     *
     */
    void update_from_additive_to_consistent_storage(
        double* vector, size_t consistency_level) const;

    // Die Methode addiert die Werte von Interface-Slaves und Interface-
    // Mastern auf und stellt mit diesen addierten Werten Level-1-Konsistenz her.
    // Dies scheint nach Multigrid-Gridtransfers noetig zu sein ( zumindest wird
    // die Methode dort aufgerufen).
    void CommUpdateReduce(double *rhs) const;

    /**
     * This was just a quick idea to identify a dof across processors and print
     * some information about it. Might be worked up some time.
     *
     * It is costly, for a vector of size N_Dof has to be communicated.
     *
     * @param process The process which sends the ping.
     * (Must be same on all processes).
     * @param dof The dof on that process which is pinged.
     * (And should be master on 'process' otherwise this method will just give
     *  a warning and -1 on all processes.)
     *
     * @return The local number of the pinged dof, -1 if the dof is not present
     * locally.
     */
    int dof_ping(size_t process, size_t dof) const;

    /// Returns an array which knows the master process of each locally known dof.
    const int *GetMaster() const
    {return Master;}
    
    /// Get the dof marker array from the fe mapper, which holds encoded
    /// information about the type of a certain d.o.f.
    /// The encoding is as follows:
    ///   i - independent (inner) master dof
    ///   D - dependent 1 master dof
    ///   d - dependent 2 master dof
    ///   m - interface master dof
    ///   s - interface slave dof
    ///   H - halo 1 slave dof
    ///   h - halo 2 slave dof
    const char* get_dof_markers() const
    {
      return Mapper->Get_DofMarker();
    }

    /// Give the total number of dof present on this communicator.
    int GetNDof() const
    {return N_Dof;}
    
    /// Return number of dimensions the communicator is used for.
    int get_n_dim() const
    {
      return Mapper->get_n_dim();
    }

    /// Return an array which contains an ad-hoc (=decomposition dependent)
    /// global numbering of the dof. Its i'th entry is the (ad hoc!) global
    /// dof number of local dof i. This is used for the interface to external
    /// solvers (MUMPS, soon the PETSc family).
    const int* Get_Local2Global() const
    { return Mapper->Get_Local2Global();}
    
    /// Returns the number of master d.o.f. on this process.
    int GetN_Master() const
    { return Mapper->GetN_Master();}
    
    /// Returns the number of interface master d.o.f. on this process
    int get_n_interface_master() const
    { return Mapper->GetN_InterfaceM();}

    /// Return a pointer to the fe space this communicator belongs to.
    const TFESpace3D* get_fe_space() const
    {
      return Mapper->get_fe_space();
    }

    //Special member functions. Rule of zero applied (class does not manage ressources).

    //Declaration of special member functions - rule of zero

    //! Default copy constructor.
    TParFECommunicator3D(const TParFECommunicator3D&) = default;

    //! Default move constructor.
    TParFECommunicator3D(TParFECommunicator3D&&) = default;

    //! Default copy assignment operator. Performs deep copy.
    TParFECommunicator3D& operator=(const TParFECommunicator3D&) = default;

    //! Default move assignment operator
    TParFECommunicator3D& operator=(TParFECommunicator3D&&) = default;

    //! Default destructor.
    ~TParFECommunicator3D() = default;

  private:
    /// The ParFEMapper, thje underlying communicative data structure.
    /// \todo Better separation of tasks of mapper and communicator!
    TParFEMapper3D *Mapper;

    /// The number of FESpaces this mapper is used for at once. Should always be 1.
    /// \todo I found this concept confusing, and therefore dropped it (CB)
    /// The concept might be revived when some day progressing to "BlockFESpaces".
    int N_Dim;

    /// The total number of d.o.f. known to this communicator - must equal the
    /// corresponding number in the FESpace.
    int N_Dof;

    /** There is three different lines of communication, each updating a different
     * class of slave d.o.f. These are:
     *  - interface master-slave communication ("MS")
     *  - halo1 communication ("H1")
     *  - halo2 communication ("H2")
     * The usual communication direction is "from master to slave", this motivated
     * the naming of the following variables. Communication follows this standard
     * direction in the methods wrapped up in "consistency_update".
     * The other direction "from slave to master" occurs in
     * "update_from_additive_to_consistent_storage", where all masters gather
     * the values of their slaves and store the sum of these values.
     *      *
     * Of the following values, there is one for each line of communication xx
     *  - N_SendDofxx: The total number of d.o.f. this process sends in xx.
     *  - Send_Infoxx: A pre-allocated send buffer of size N_SendDofxx
     *  - Recv_Infoxx: A pre-allocted receive buffer of size N_xxxxxxx
     *  - N_DofSendxx: An array of size 'mpi_size', listing how many values to
     *                 send to each process (its entries should sum up to N_SendDofxx)
     *  - N_DofRecvxx: An array of size 'mpi_size', listing how manz values to
     *                 receive from each process (its entries should sum up to N_xxxxxxx)
     *  - sdisplxx:    An array of size 'mpi_size', containing send displacements.
     *                 No holes or overlap should occur.
     *  - rdisplxx:    An array of size 'mpi_size', containing receive displacements.
     *                 No holes or overlap should occur.
     *  - DofSendxx:   This array helps to fill the send buffer - it contains the local
     *                 d.o.f. indices of those d.o.f. to be send in communication xx.
     *  - DofRecvxx:   This array helps to interpret the receive buffer - it contains
     *                 the local d.o.f. indices of those d.o.f. to be received in communication xx.
     *  - N_xxxxxxx:   The total number of d.o.f. to be received in communication xx.
     *
     *  \todo The first values in each row are from the outdated former communication
     *        structure and should be removed!
     */
    int N_SendDof, N_SendDofMS, N_SendDofH1, N_SendDofH2;

    double *Send_Info, *Send_InfoMS, *Send_InfoH1, *Send_InfoH2;

    double *Recv_Info, *Recv_InfoMS, *Recv_InfoH1, *Recv_InfoH2;

    int *N_DofSend, *N_DofSendMS, *N_DofSendH1, *N_DofSendH2;

    int *N_DofRecv, *N_DofRecvMS, *N_DofRecvH1, *N_DofRecvH2;

    int *sdispl, *sdisplMS, *sdisplH1, *sdisplH2;

    int *rdispl, *rdisplMS, *rdisplH1, *rdisplH2;

    int *DofSend, *DofSendMS, *DofSendH1, *DofSendH2;

    int *DofRecv, *DofRecvMS, *DofRecvH1, *DofRecvH2;

    int N_Slave, N_InterfaceS, N_Halo1, N_Halo2;

    /// The MPI Communicator controlling the communication. Should usually be
    /// MPI_COMM_WORLD (since we do not use Comm forcefully in the implementation,
    /// but mostly MPI_COMM_WORLD, only that will work at the moment anyway).
    MPI_Comm Comm;

    const int *Master;

    std::vector<double> averaging_multiplier;
    friend class BufferInformation;
};

#endif
#endif