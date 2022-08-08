/** ***************************************************************************
 *
 * @name       TParFEMapper3D
 * @brief      Class containing all info needed for communication between subdomains.
 *
 *             Objects of this type are set up after the domain is split, so
 *             there exists a particular one in each process. The data members
 *             refer to the data of the processes' Domain and FESpace.
 *             The class holds the information which was formerly stored in
 *             TParVector3D.
 *
 * @author     Sashikumaar Ganesan
 * @date       2015/04/24
 *
 * @ruleof0
 *
 ******************************************************************************/

#ifdef _MPI

#ifndef __TParFEMapper3D__
#define __TParFEMapper3D__

#include "mpi.h"

#include <FESpace3D.h>
#include <Structure.h>

class FESpace3D;

class TParFEMapper3D
{
  protected:
    //! The used MPI communicator.
    MPI_Comm Comm;
     
    /**
     * Spatial dimension of the  range of the FE functions the Mapper belongs to.
     *
     * Put it to 1 for convection-diffusion, 3 for velocity in NSE, 1 for pressure in NSE.
     */
    int N_Dim;
    
    //! Number of degrees of freedom in the process.
    int N_Dof;
    
#ifdef _OMP
    int* RowPtr;

    //! THIS IS UNUSED EXCEPT FOR ONE PLACE IN HYBRID - REMOVE!
    int* KCol;
#endif
    
    //! The underlying fe space for which the communications are needed.
    const TFESpace3D *FESpace;
    
    /** array containing the master of the dofs **/ 
    int* Master;

    int* OwnDofs;
    
    /** An array which contains chars for the type of the dof. This is internally used
     * when constructing the dof map, and the meaning of the stored chars changes between program parts.
     **/
    char *DofMarker;
    
    /** array which carries info across procs**/
    double *Send_Info, *Recv_Info; 
    
    /**array to store globalDofNo of all its dofs**/
    int *Local2Global;
    
    /** info of number of dofs to be send/recv from other procs**/
    int *sdispl, *rdispl;
    
    /* The following data members have been added when implementing the master/slave/halo
     * mapping strategy and are used in there.
     */

    int N_InterfaceM, N_InterfaceS, N_Halo1, N_Halo2, N_Dept1, N_Dept2;
    
    int N_Slave, N_Halo, N_Master, N_Dept, N_Int, N_OwnDof;
    
    int *N_DofSend, *N_DofRecv, *N_DofSendMS, *N_DofSendH1, *N_DofSendH2, *N_DofRecvMS, *N_DofRecvH1, *N_DofRecvH2;
     
    int *sdisplMS, *sdisplH1, *sdisplH2, *rdisplMS, *rdisplH1, *rdisplH2;
     
    int N_SendDof, N_SendDofMS, N_SendDofH1, N_SendDofH2;
     
    int *DofSend, *DofSendMS, *DofSendH1, *DofSendH2, *DofRecv, *DofRecvMS, *DofRecvH1, *DofRecvH2;
     
    double *Send_InfoMS, *Send_InfoH1, *Send_InfoH2, *Recv_InfoMS, *Recv_InfoH1, *Recv_InfoH2; 
    
    int *NewGN, *Reorder, *Reorder_M, *Reorder_I, *Reorder_D1, *Reorder_D2;
     
    int N_CMaster, N_CDept1, N_CDept2, N_CInt, *ptrCMaster, *ptrCDept1, *ptrCDept2, *ptrCInt;
     
  public:
    /** Default constructor.
     *
     * This is intended to be used only when default initializing ParFEMapper
     * as a class member. It sets a lot of zeroes and allocates zero/small arrays -
     * in order for the Destructor not to throw when called upon a default constructed
     * object.
     *
     * @note Try not to call any methods on a default constructed object,
     * because it holds a nullptrs to FESpace.
     *
     */
    TParFEMapper3D();

    /**
     * The standard constructor. Idea: Add a constructor which takes a TStructure/TMatrix?
     * @param[in] N_dim The number of dimensions associated with a certain dof (range), e.g. "3" if
     * the ParFEMapper is for the velocity of a 3D NSE problem or "1" for a 3D CDR problem.
     * @param[in] fespace the FE space the dofs belong to.
     */
#ifndef _OMP
    TParFEMapper3D(int N_dim, const TFESpace3D *fespace);
#else
    TParFEMapper3D(int N_dim, const TFESpace3D *fespace, int *rowptr, int *kcol);
#endif
    
    //! A getter method which gives all information needed by the ParFECommunicator.
    //! All parameters are out-parameters.
    void GetCommInfo(int &n_Dim, int &n_Dof,
		     int &n_SendDof, int &n_SendDofMS, int &n_SendDofH1, int &n_SendDofH2,
		     double *&send_Info, double *&send_InfoMS, double *&send_InfoH1, double *&send_InfoH2,
		     double *&recv_Info, double *&recv_InfoMS, double *&recv_InfoH1, double *&recv_InfoH2,
		     int *&n_DofSend, int *&n_DofSendMS, int *&n_DofSendH1, int *&n_DofSendH2,
		     int *&n_DofRecv, int *&n_DofRecvMS, int *&n_DofRecvH1, int *&n_DofRecvH2,
		     int *&Sdispl, int *&SdisplMS, int *&SdisplH1, int *&SdisplH2,
		     int *&Rdispl, int *&RdisplMS, int *&RdisplH1, int *&RdisplH2,
		     int *&dofSend, int *&dofSendMS, int *&dofSendH1, int *&dofSendH2,
		     int *&dofRecv, int *&dofRecvMS, int *&dofRecvH1, int *&dofRecvH2,
		     int &n_Slave, int &n_InterfaceS, int &n_Halo1, int &n_Halo2) const
    {
      n_Dim        = N_Dim; 
      n_Dof        = N_Dof;
      
      n_SendDof    = N_SendDof;
      n_SendDofMS  = N_SendDofMS;
      n_SendDofH1  = N_SendDofH1;
      n_SendDofH2  = N_SendDofH2;
      
      send_Info    = Send_Info;
      send_InfoMS  = Send_InfoMS;
      send_InfoH1  = Send_InfoH1;
      send_InfoH2  = Send_InfoH2;
      
      recv_Info    = Recv_Info;
      recv_InfoMS  = Recv_InfoMS;
      recv_InfoH1  = Recv_InfoH1;
      recv_InfoH2  = Recv_InfoH2;
      
      n_DofSend    = N_DofSend;
      n_DofSendMS  = N_DofSendMS;
      n_DofSendH1  = N_DofSendH1;
      n_DofSendH2  = N_DofSendH2;
    
      n_DofRecv    = N_DofRecv;
      n_DofRecvMS  = N_DofRecvMS;
      n_DofRecvH1  = N_DofRecvH1;
      n_DofRecvH2  = N_DofRecvH2;
    
      Sdispl       = sdispl;
      SdisplMS     = sdisplMS;
      SdisplH1     = sdisplH1;
      SdisplH2     = sdisplH2;
      
      Rdispl       = rdispl;
      RdisplMS     = rdisplMS;
      RdisplH1     = rdisplH1;
      RdisplH2     = rdisplH2;
    
      dofSend      = DofSend;
      dofSendMS    = DofSendMS;
      dofSendH1    = DofSendH1;
      dofSendH2    = DofSendH2;
      
      dofRecv      = DofRecv;
      dofRecvMS    = DofRecvMS;
      dofRecvH1    = DofRecvH1;
      dofRecvH2    = DofRecvH2;
    
      n_Slave      = N_Slave;
      n_InterfaceS = N_InterfaceS;
      n_Halo1      = N_Halo1;
      n_Halo2      = N_Halo2;
    }
    
    int GetN_Master() const
    {
      return N_Master;
    }
    
    const int *GetMaster() const
    {
      return Master;
    }
    
    const char *Get_DofMarker() const
    {
      return DofMarker;
    }
    
    const int* GetReorder_M() const
    {
      return Reorder_M;
    }
    
    const int* GetReorder_I() const
    {
      return Reorder_I;
    }
    
    const int* GetReorder_D1() const
    {
      return Reorder_D1;
    }
    
    const int* GetReorder_D2() const
    {
      return Reorder_D2;
    }
    
    int GetN_InterfaceM() const
    {
      return N_InterfaceM;
    }
    
    int GetN_Int_light() const
    {
      return N_Int;
    }
    
    int GetN_Dept1() const
    {
      return N_Dept1;
    }
    
    int GetN_Dept2() const
    {
      return N_Dept2;
    }
    
    const int* Get_Local2Global() const
    {
      return Local2Global;
    }

    /// Return number of dimensions the mapper is used for.
    int get_n_dim() const
    {
      return N_Dim;
    }

    /// Return a pointer to the fe space this communicator belongs to.
    const TFESpace3D* get_fe_space() const
    {
      return FESpace;
    }

#ifdef _OMP
    void Color(int &numColors, int *&ptrColors, char type);
    
    int GetN_CMaster()
    {return N_CMaster;}
    int* GetptrCMaster()
    {return ptrCMaster;}
    
    int GetN_CDept1()
    {return N_CDept1;}
    int* GetptrCDept1()
    {return ptrCDept1;}
    
    int GetN_CDept2()
    {return N_CDept2;}
    int* GetptrCDept2()
    {return ptrCDept2;}
    
    int GetN_CInt()
    {return N_CInt;}
    int* GetptrCInt()
    {return ptrCInt;}
    
#endif

    //Special member functions. Rule of three (and a half: swap) applied.
    /** All special member function rely on the switch TDatabase::ParamDB->MapperType
     * to be the same at the time of a ParFEMapper object's creation/copying/moving/destruction.
     * That switch decides, whether the map is constructed the old way (master-slave)
     * or the new way (master-slave-halo).
     * TODO This is a bit unsafe and could be changed.
     * TODO Implement additional things for HYBRID case.
     * TODO No move implemented, class will always copy instead. (Not performance critical)
     * TODO These special member functions are untested and there is a chance that they fail.
     */

    //! Copy constructor. Performs deep copy of all owned members.
    //! Does only shallow copy of: RowPtr, KCol, FESpace. Keep that in mind!
    TParFEMapper3D(const TParFEMapper3D& other);

    /** Swap function used for copy-and swap in copy assignment.
     * @param[in,out] first The object to be swapped with second.
     * @param[in,out] second The object to be swapped with first.
     */
    friend void swap(TParFEMapper3D& first, TParFEMapper3D& second);


    /** Copy assignment operator. Uses copy and swap idiom.
    * @note Note that the object to be copied is passed by value -
    * so there is a call to the copy constructor, in the function
    * body there happens a swap.
    * Does only shallow copy of: RowPtr, KCol, FESpace. Keep that in mind!
	*/
    TParFEMapper3D& operator=(TParFEMapper3D other);

    //! Destructor.
    ~TParFEMapper3D();

  protected:
    /**
     *  Allocates and fills the member Local2Global, which contains global identification
     *	numbers for the dofs known to this process. The method is only called by the
     *  constructor.
     *	@note The method deals with memory which is only allocated when the
     *	master/slave/halo mapping concept is used.
     */
    void Assign_GlobalDofNo();

    /** Constructs the entire mapping between processes, using the master/slave/halo mapping concept.
     *  Called in the constructor, if Database::ParamDB->MapperType is not 2.
     */
    void ConstructDofMap_Master_Halo();

    /** @brief Helper function for ConstructDofMap_Master_Halo().
     *  Assembles consensus on which process is responsible for shared DOFs.
     */
    void NegotiateDOFOwnership(int* LocalIndex, int** MappingData);

    /** @brief Helper function for ConstructDofMap_Master_Halo().
     *  Optionally performs some sanity checks on DOF distribution.
     */
    void CheckDOFValidity(int* LocalIndex, int** MappingData);

    /** @brief Helper function for ConstructDofMap_Master_Halo().
     *  Finalizes DOF type marking.
     */
    void MarkDOFs();

    /** @brief Helper function for ConstructDofMap_Master_Halo().
     *  Sets up process-to-process mapping of MS DOFs.
     */
    void MapDOFsMS(int* LocalIndex, int** MappingData);

    /** @brief Helper function for ConstructDofMap_Master_Halo().
     *  Sets up process-to-process mapping of halo DOFs.
     */
    void MapDOFsHalo(int* LocalIndex, int** MappingData);

    /** @brief Constructs the entire mapping between processes, using the master/slave mapping concept.
     *  Called in the constructor, if Database::ParamDB->MapperType equals 2.
     */
    void ConstructDofMap();
};

#endif
#endif






















