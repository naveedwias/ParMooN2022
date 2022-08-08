#ifdef _MPI

#include "mpi.h"
#include <ParFEMapper3D.h>
#include <Database.h>
#include <SubDomainJoint.h>
#include <Edge.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define GLOBAL_NO 0
#define DOF_NO 1

#define HALOCELL 0
#define NONHALO  1


#ifndef _OMP
TParFEMapper3D::TParFEMapper3D(int N_dim, const TFESpace3D *fespace)
#else
TParFEMapper3D::TParFEMapper3D(int N_dim, const TFESpace3D *fespace, int *rowptr, int *kcol)
#endif
{
  N_Dim       = N_dim; 
  Comm        = TDatabase::ParamDB->Comm;
  FESpace     = fespace;

#ifdef _OMP
  RowPtr      = rowptr;
  KCol        = kcol;
#endif
 
  N_Dof = FESpace->get_n_dof();
  
  if (TDatabase::ParamDB->MapperType != 2)
  {
    ConstructDofMap_Master_Halo();

    // CB Moved the following call here, because segfaults when calling
    // Assign_GlobalDofNo without ConstructDofMap_Master_Halo.
    // (ConstructDofMap does not allocate all the memory which
    // Assign_GlobalDofNo watns to access)

    Assign_GlobalDofNo();
#ifdef _OMP
    Color(N_CInt,ptrCInt,'i');
    Color(N_CMaster,ptrCMaster,'m');
    Color(N_CDept1,ptrCDept1,'D');
    Color(N_CDept2,ptrCDept2,'d');
#endif
  }
  else
  {
    ConstructDofMap();
  }
}

TParFEMapper3D::TParFEMapper3D()
{
    //Take the Database communicator
    Comm = TDatabase::ParamDB->Comm;

    FESpace = nullptr;

#ifdef _OMP
    RowPtr = nullptr;
    KCol = nullptr;
#endif

    //assign non-array built-in type data members
    N_Dim = 0;
    N_Dof = 0;

    N_InterfaceM = 0;
    N_InterfaceS = 0;
    N_Halo1 = 0;
    N_Halo2 = 0;
    N_Dept1 = 0;
    N_Dept2 = 0;

    N_Slave = 0;
    N_Halo = 0;
    N_Master = 0;
    N_Dept = 0;
    N_Int = 0;
    N_OwnDof = 0;

    N_SendDof = 0;
    N_SendDofMS = 0;
    N_SendDofH1 = 0;
    N_SendDofH2 = 0;

    N_CMaster = 0;
    N_CDept1 = 0;
    N_CDept2 = 0;
    N_CInt = 0;

    //determine rank and size
    int mpiRank, mpiSize;
    MPI_Comm_rank(Comm, &mpiRank);
    MPI_Comm_size(Comm, &mpiSize);

    //switch over the two MapperTypes
    if(TDatabase::ParamDB->MapperType != 2)
    { //master-slave-halo mapping
      Master = new int[N_Dof];
      DofMarker = new char[N_Dof];

      sdispl = new int[mpiSize];
      rdispl = new int[mpiSize];
      sdisplMS = new int[mpiSize];
      rdisplMS = new int[mpiSize];
      sdisplH1 = new int[mpiSize];
      rdisplH1 = new int[mpiSize];
      sdisplH2 = new int[mpiSize];
      rdisplH2 = new int[mpiSize];

      N_DofSend = new int[mpiSize];
      N_DofSendMS = new int[mpiSize];
      N_DofSendH2 = new int[mpiSize];
      N_DofSendH1 = new int[mpiSize];

      N_DofRecv = new int[mpiSize];
      N_DofRecvMS = new int[mpiSize];
      N_DofRecvH2 = new int[mpiSize];
      N_DofRecvH1 = new int[mpiSize];

      OwnDofs = new int[N_OwnDof];

      DofSend = new int[N_SendDof];
      DofSendMS  = DofSend;
      DofSendH1  = DofSend + N_SendDofMS;
      DofSendH2  = DofSend + N_SendDofMS + N_SendDofH1; //pointer arithmetics!


      DofRecv    = new int[N_InterfaceS+N_Halo1+N_Halo2];
      DofRecvMS  = DofRecv;
      DofRecvH1  = DofRecv + N_InterfaceS;
      DofRecvH2  = DofRecv + N_InterfaceS + N_Halo1; //pointer arithmetics!

      Reorder = new int[N_Dof];
      Reorder_M  = Reorder;
      Reorder_I  = Reorder + N_InterfaceM;
      Reorder_D1 = Reorder + N_InterfaceM + N_Int;
      Reorder_D2 = Reorder + N_InterfaceM + N_Int + N_Dept1; //pointer arithmetics!

      NewGN = new int[N_Dof];

      //TODO CB Maybe I do not like the following - conditional assignment mixed with deterministic.
      if(N_SendDof>0)
      {
        Send_Info   = new double[N_SendDof*N_Dim];
      }
      Send_InfoMS = Send_Info;
      Send_InfoH1 = Send_Info + N_SendDofMS*N_Dim;
      Send_InfoH2 = Send_Info + N_SendDofMS*N_Dim + N_SendDofH1*N_Dim;

      if(N_Slave>0){
        Recv_Info   = new double[N_Slave*N_Dim];
      }
      Recv_InfoMS = Recv_Info;
      Recv_InfoH1 = Recv_Info + N_InterfaceS*N_Dim;
      Recv_InfoH2 = Recv_Info + N_InterfaceS*N_Dim + N_Halo1*N_Dim;

      Local2Global = new int[N_Dof];

    }
    else
    { //master-slave mapping
      Master = new int[N_Dof];
      sdispl = new int[mpiSize];
      rdispl = new int[mpiSize];

      N_DofSend = new int[mpiSize];
      N_DofRecv = new int[mpiSize];

      OwnDofs = new int[N_OwnDof];

      DofSend = new int[N_SendDof];
      DofRecv    = new int[N_Slave];

        if(N_SendDof>0)
         Send_Info = new double[N_SendDof*N_Dim];
        if(N_Slave>0)
         Recv_Info = new double[N_Slave*N_Dim];
    }



}

int GetLocalIndex(int, int *array, int val)
{
  int m = 0;

  while (array[m] != val)
  {
    m++;
  }

  return m;
}

int find_min(int *arr, int N, char *temp_arr)
{
  int i, j;
  int min = -1;
  j = 0;

  for (i = 0; i < N; i++)
  {
    if (temp_arr[i] == 'x')
    {
      continue;
    }

    if (min < 0 || arr[i] < min)
    {
      min = arr[i];
      j = i;
    }
  }

  return j;
}

void TParFEMapper3D::NegotiateDOFOwnership(int* LocalIndex, int** MappingData)
{
  int rank, size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  auto Coll   = FESpace->GetCollection();
  int N_Cells = Coll->GetN_Cells();

  //------------------------------------------//
  /** Master DOF verification by other ranks **/
  //------------------------------------------//

  // NOTE: "Verification" here actually means that the ranks negotiate which
  // of them takes ownership over which DOFs
 
  int* N_DOFtobeverified_otherRank;        // Array containing how many DOF's need to be verified from each processor (rank) 
  int Total_DOFtobeverified_otherRank = 0; // Total no of DOF's that need to be verified from other ranks
  int* N_DOFtobeverified_thisRank;         // Array containing how many DOF's need to be verified by this rank for each other processor (rank)
  int Total_DOFtobeverified_thisRank = 0;  // Total no of DOF verified by this rank
  
  int *master_recvbuf;
  int *master_sendbuf;
 
  int** mapping_sendbuf = new int*[2];
  int** mapping_recvbuf = new int*[2];
  
  int* buffer_index_per_rank = new int[size]; // each rank's index into the send buffer (counts up as the buffer is filled)
  
  N_DOFtobeverified_otherRank = new int[size];
  N_DOFtobeverified_thisRank  = new int[size];
  
  memset(N_DOFtobeverified_otherRank, 0, size * sizeof(int));
  memset(N_DOFtobeverified_thisRank, 0, size * sizeof(int)); 
  
  // only the dofs in halo cells (not dependent cells) needs to be verified
  for (int i = 0; i < N_Dof; i++)
  {
    if (DofMarker[i] == 'h')
    {
      N_DOFtobeverified_otherRank[Master[i]]++;
      Total_DOFtobeverified_otherRank++;
    }
  }
  
  // communicate number of verification requests
  MPI_Alltoall(N_DOFtobeverified_otherRank, 1, MPI_INT, N_DOFtobeverified_thisRank, 1, MPI_INT, Comm);
  
  for (int i = 0; i < size; i++)
  {
    Total_DOFtobeverified_thisRank += N_DOFtobeverified_thisRank[i];
  }
  
  for (int i = 0; i < 2; i++)
  {
    mapping_sendbuf[i] = new int[Total_DOFtobeverified_otherRank];
    mapping_recvbuf[i] = new int[Total_DOFtobeverified_thisRank];
  }
  
  // compute send buffer offsets
  sdispl[0] = 0;
  for (int i = 1; i < size; i++)
  {
    sdispl[i] = N_DOFtobeverified_otherRank[i - 1] + sdispl[i - 1];
  }
  
  // compute receive buffer offsets
  rdispl[0] = 0;
  for (int i = 1; i < size; i++)
  {
    rdispl[i] = N_DOFtobeverified_thisRank[i - 1] + rdispl[i - 1];
  }
  
  // duplicate receive buffer offsets as starting indices
  memcpy(buffer_index_per_rank, sdispl, size * sizeof(int));
  
  // fill send buffer with DOF mapping data
  for (int i = 0; i < N_Dof; i++)
  {
    if (DofMarker[i] == 'h')
    {
      // we ask the lowest-ranked process we know is aware of each halo DOF

      mapping_sendbuf[GLOBAL_NO][buffer_index_per_rank[Master[i]]] = MappingData[GLOBAL_NO][i];
      mapping_sendbuf[DOF_NO][buffer_index_per_rank[Master[i]]]    = MappingData[DOF_NO][i];
      
      buffer_index_per_rank[Master[i]]++;
    }
  }
  
  // communicate DOF mapping to other ranks
  
  MPI_Alltoallv(mapping_sendbuf[GLOBAL_NO], N_DOFtobeverified_otherRank, sdispl, MPI_INT,
                mapping_recvbuf[GLOBAL_NO], N_DOFtobeverified_thisRank,  rdispl, MPI_INT, Comm);
  
  MPI_Alltoallv(mapping_sendbuf[DOF_NO], N_DOFtobeverified_otherRank, sdispl, MPI_INT,
                mapping_recvbuf[DOF_NO], N_DOFtobeverified_thisRank,  rdispl, MPI_INT, Comm);
  
  master_sendbuf = new int[Total_DOFtobeverified_thisRank];
  master_recvbuf = new int[Total_DOFtobeverified_otherRank];

  // assemble master info for each verification request
  for (int i = 0; i < Total_DOFtobeverified_thisRank; i++)
  {
    int temp_globalno = mapping_recvbuf[GLOBAL_NO][i];
    int temp = GetLocalIndex(N_Cells, LocalIndex, temp_globalno);
    temp = FESpace->get_global_dof(temp, mapping_recvbuf[DOF_NO][i]);
    
    master_sendbuf[i] = Master[temp];
  }
  
  // communicate master info to other ranks (note that the offsets are switched since we're sending data back)
  MPI_Alltoallv(master_sendbuf, N_DOFtobeverified_thisRank,  rdispl, MPI_INT,
                master_recvbuf, N_DOFtobeverified_otherRank, sdispl, MPI_INT, Comm);

  // retrieve master info. 
  for (int i = 0; i < Total_DOFtobeverified_otherRank; i++)
  {
    int temp_globalno = mapping_sendbuf[GLOBAL_NO][i];
    int temp_dofno    = mapping_sendbuf[DOF_NO][i];
    int temp = GetLocalIndex(N_Cells, LocalIndex, temp_globalno);
    temp = FESpace->get_global_dof(temp, temp_dofno);
    
    if (DofMarker[temp] != 'h')
    {
      ErrThrow("Error : This degree of Freedom (", temp, ") didn't require verification");
    }
    
    // if the other process we asked knows of a lower process adjacent to the halo DOF, we update it here
    // NB: this should not cause inconsistencies since the other process is the owner of an adjacent cell
    // (and as such should have every cell adjacent to this DOF as a halo cell)
    Master[temp] = master_recvbuf[i];
  }
  
  int numberOfMyMasterDOFs = 0;
  for (int i = 0; i < N_Dof; i++)
  {
    if (Master[i] == rank)
    {
      numberOfMyMasterDOFs++;
    }
  }
  
  if (numberOfMyMasterDOFs == 0)
  {
    ErrThrow(rank, ": No degrees of freedom assigned to this process.");
  }

  for (int i = 0; i < 2; i++)
  {
    delete[] mapping_sendbuf[i];
    delete[] mapping_recvbuf[i];
  }

  delete[] mapping_sendbuf;
  delete[] mapping_recvbuf;

  delete[] N_DOFtobeverified_otherRank;
  delete[] N_DOFtobeverified_thisRank;

  delete[] master_recvbuf;
  delete[] master_sendbuf;

  delete[] buffer_index_per_rank;
}

void TParFEMapper3D::CheckDOFValidity(int* LocalIndex, int** MappingData)
{
  int rank, size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  auto Coll      = FESpace->GetCollection();
  int N_Cells    = Coll->GetN_Cells();
  int N_OwnCells = Coll->GetN_OwnCells();
  
  int* buffer_index_per_rank = new int[size]; // each rank's index into the send buffer (counts up as the buffer is filled)

  //---------------------------------------------------------------------------------------------------------------------------------------------------------------
  int total_interface_dofs        = 0;                    // these are the dofs which lie on the interface of sub domains (total over all sub domains)
  int N_own_interface_dofs        = 0;                    // these are the own interface dofs
  int N_interface_dofs            = 0;                    // these are the interface dofs including shared ones
  int T_interface_dofs            = 0;                    // these are the total of N_interface_dofs over all sub domains
  int total_own_dofs              = 0;                    // these are the dofs which belongs to my rank (total implies total over all sub domains)
  int start = 0,           end    = 0;                    // temporary variables used for numbering the dofs
  int max_n_ranks_interface_dofs  = 0;                    // maximum number of ranks sharing any interface dof

  int *GlobalDofNo                = new int[N_Dof];       // this array is used to assign unique global dof no. over all sub domains 
      memset(GlobalDofNo, -1, N_Dof * sizeof(int));       // all is set to -1 for a default value
  int *N_Dof_Slave                = new int[size];        // array of number of slave interface dofs to be verified by other ranks
      memset(N_Dof_Slave, 0, size * sizeof(int));
  int *N_Dof_Master               = new int[size];        // array of interface dofs to be verified by this rank for other ranks 
  int **MasterPos                 = new int*[2];          // a double array for storing the position info for master interface dofs
  int **SlavePos                  = new int*[2];          // a double array for storing the position info for slave interface dofs 
  int *masterInfo, *slaveInfo;                            // these are used to update the global dof no for interface master-slave dofs
  int *GlobalDofNo_interface;                             // this array contains only the globaldof no. of interface dofs
  int *all_own_dofs_info          = new int[size];        // array of N_own dofs by each rank
  int *all_interface_dofs_info    = new int[size];        // array of N_own interface dofs by each rank
  int *all_T_interface_dofs_info  = new int[size];        // array of interface dofs(including shared ones) by each rank
  int *all_GlobalDofNo;                                   // array of global interface dof no from all ranks
  int *N_ranks_per_interface_dofs;                        // array of count of number of ranks sharing an interface dof
  int *N_allocated_masters        = new int[size];        // array containing number of masters(interface dofs) allocated to each processors
  char *tempc                     = new char[size];       // temporary array
  char **Master_Table             = new char*[size];      // a Table to mark the shared interface dofs
  //-----------------------------------------------------------------------------------------------------------------------------------------------------------------

  // mark all (interface) dofs as 'z'
  for (int i = N_OwnCells; i < N_Cells; i++)
  {
    auto DOF = FESpace->GetGlobalDOF(i);
    int N_LocDof = FESpace->get_n_local_dof(i);
    
    // z is both interface and dependent dof, while d is only dependent dof
    for (int j = 0; j < N_LocDof; j++)
    {
      int N = DOF[j];
      
      if (DofMarker[N] == 'd')
      {
         DofMarker[N] = 'z';
      }
    }
  }
  
  // count total number of OWN interface_dofs and own dofs
  N_OwnDof = 0;
  for (int i = 0; i < N_Dof; i++)
  {
    if (Master[i] == rank)
    {
      N_OwnDof++;
      
      if (DofMarker[i] == 'z')
      {
        N_own_interface_dofs++;
      }
    }
  }
  
  // Send the count info to all ranks
  MPI_Allgather(&N_own_interface_dofs, 1, MPI_INT, all_interface_dofs_info, 1, MPI_INT, Comm);
  MPI_Allgather(&N_OwnDof,             1, MPI_INT, all_own_dofs_info,       1, MPI_INT, Comm);
  
  // sum up all ranks' own and own interface dofs  
  for (int aa = 0; aa < size; aa++)
  {
    total_interface_dofs += all_interface_dofs_info[aa];
    total_own_dofs       += all_own_dofs_info[aa];
  }
  
  // "start" is the index of this rank's first interface dof
  for (int aa = 0; aa < rank; aa++)
  {
    start += all_interface_dofs_info[aa];
  }
  
  int Total_DOFtobeverified_otherRank = 0;
  int Total_DOFtobeverified_thisRank = 0;
  
  // number the own interface dofs
  // count the total dofs to be verified by this and other rank
  for (int i = 0; i < N_Dof; i++)
  {
    if (DofMarker[i] == 'z')
    {
      // this is an interface dof

      if (Master[i] == rank)
      {
        // this is ours - number it
        GlobalDofNo[i] = start;
        start++;
      }
      else
      {
        // this is not ours - count it up as a verification request
        N_Dof_Slave[Master[i]]++;
        Total_DOFtobeverified_otherRank++;
      }
    }
  }

  // communicate MS counts to other ranks
  MPI_Alltoall(N_Dof_Slave, 1, MPI_INT, N_Dof_Master, 1, MPI_INT, Comm);
  
  // count up other ranks' verification requests
  for (int i = 0; i < size; i++)
  {
    Total_DOFtobeverified_thisRank += N_Dof_Master[i];
  }
 
  for (int i = 0; i < 2; i++)
  {
    MasterPos[i] = new int[Total_DOFtobeverified_thisRank];
    SlavePos[i]  = new int[Total_DOFtobeverified_otherRank];
  }
  
  masterInfo = new int[Total_DOFtobeverified_thisRank];
  slaveInfo  = new int[Total_DOFtobeverified_otherRank];
  
  // count up send/receive offsets

  rdispl[0] = 0;
  sdispl[0] = 0;

  for (int i = 1; i < size; i++)
  {
    rdispl[i] = rdispl[i - 1] + N_Dof_Master[i - 1];
    sdispl[i] = sdispl[i - 1] + N_Dof_Slave[i - 1];
  }
 
  // duplicate send offsets to use as indices
  memcpy(buffer_index_per_rank, sdispl, size * sizeof(int));
 
  // store the position info of the slave dofs
  for (int i = 0; i < N_Dof; i++)
  {
    if (DofMarker[i] == 'z')
    {
      if (Master[i] != rank)
      {
        // this is an interface dof that isn't ours - write mapping data

        SlavePos[GLOBAL_NO][buffer_index_per_rank[Master[i]]] = MappingData[GLOBAL_NO][i];
        SlavePos[DOF_NO][buffer_index_per_rank[Master[i]]]    = MappingData[DOF_NO][i];
        
        buffer_index_per_rank[Master[i]]++;
      }
    }
  }

  // send the slave info to masters
  MPI_Alltoallv(SlavePos[GLOBAL_NO], N_Dof_Slave, sdispl, MPI_INT, MasterPos[GLOBAL_NO], N_Dof_Master, rdispl, MPI_INT, Comm);
  MPI_Alltoallv(SlavePos[DOF_NO],    N_Dof_Slave, sdispl, MPI_INT, MasterPos[DOF_NO],    N_Dof_Master, rdispl, MPI_INT, Comm);
 
  // gather response to verification requests
  start = 0; // rank 0's interface dofs start at 0
  for (int aa = 0; aa < size; aa++)
  {
    // add rank aa's interface dof count
    end = start + N_Dof_Master[aa];
    
    for (int i = start; i < end; i++)
    {
      // get the cell with global cell no = MasterPos[GLOBAL_NO][i]
      int temp = GetLocalIndex(N_Cells, LocalIndex, MasterPos[GLOBAL_NO][i]);
      
      // the master of this dof, i should be rank
      if (Master[FESpace->get_global_dof(temp, MasterPos[DOF_NO][i])] != rank)
      {
        ErrThrow("....................Wrong Master.....................\n");
      }
      
      // this dof should be a category z dof
      if (DofMarker[FESpace->get_global_dof(temp, MasterPos[DOF_NO][i])] != 'z')
      {
        ErrThrow("....................Wrong Dof.........................\n");
      }
      
      // this dof should have a GlobalDofNo
      if (GlobalDofNo[FESpace->get_global_dof(temp, MasterPos[DOF_NO][i])] == -1)
      {
        ErrThrow("....................Wrong GlobalDofNo.........................\n");
      }
      
      // everything checks out, write it down
      masterInfo[i] = GlobalDofNo[FESpace->get_global_dof(temp, MasterPos[DOF_NO][i])];
    }

    // the next rank's interface dofs start right after this one's
    start = end;
  }

  // distribute global numbering to other ranks
  MPI_Alltoallv(masterInfo, N_Dof_Master, rdispl, MPI_INT,
                slaveInfo, N_Dof_Slave, sdispl, MPI_INT, Comm);

  // apply received global numbering
  start = 0;
  for (int aa = 0; aa < size; aa++)
  {
    end = start + N_Dof_Slave[aa];

    for (int i = start; i < end; i++)
    {
      int temp_globalno = SlavePos[GLOBAL_NO][i];
      int temp_dofno    = SlavePos[DOF_NO][i];
      int temp = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
      temp = FESpace->get_global_dof(temp, temp_dofno);
      
      // the master of this dof, i should be aa
      if (Master[temp] != aa)
      {
        ErrThrow("...................2.Wrong Master.....master[%d]=%d....rank=%d............", temp, Master[temp], rank);
      }
      
      // this dof should be a category z dof
      if (DofMarker[temp] != 'z')
      {
        ErrThrow("...................2.Wrong Dof.........................");
      }
      
      // this dof should not have a GlobalDofNo
      if (GlobalDofNo[temp] != -1)
      {
        ErrThrow("...................2.Wrong GlobalDofNo.........................");
      }

      //assign the global dof no to the slave obtained from master
      GlobalDofNo[temp] = slaveInfo[i];
    }

    start = end;
  }
  
  // now we mark the own dofs apart from the interface ones

  start = total_interface_dofs;
  for (int aa = 0; aa < rank; aa++)
  {
    start += all_own_dofs_info[aa];
    aa++;
  }

  // just go ahead and number them, these are ours
  for (int i = 0; i < N_Dof; i++)
  {
    if (DofMarker[i] != 'z')
    {
      if (Master[i] == rank)
      {
        GlobalDofNo[i] = start;
        start++;
      }
    }
    
    if (Master[i] == rank)
    {
      if (GlobalDofNo[i] == -1)
      {
        ErrThrow("......................This shudnt happen l.679.................");
      }
    }
    
    if (DofMarker[i] == 'z')
    {
      if (GlobalDofNo[i] >= total_interface_dofs)
      {
        ErrThrow("......................This shudnt happen l.687.................");
      }
      
      // also count interface dofs
      N_interface_dofs++;
    }
  }
  
  // gather just the interface dofs into a separate array
  GlobalDofNo_interface = new int[N_interface_dofs];
  int idx = 0;
  for (int i = 0; i < N_Dof; i++)
  {
    if (DofMarker[i] == 'z')
    {
      GlobalDofNo_interface[idx] = GlobalDofNo[i];
      idx++;
    }
  }

  // gather other ranks' interface dof counts
  MPI_Allgather(&N_interface_dofs, 1, MPI_INT, all_T_interface_dofs_info, 1, MPI_INT, Comm);
  
  // count all (perspectives on) interface dofs
  T_interface_dofs = 0;
  for (int aa = 0; aa < size; aa++)
  {
    T_interface_dofs += all_T_interface_dofs_info[aa];
  }
  
  all_GlobalDofNo = new int[T_interface_dofs];
  
  // compute offsets into global interface dof array
  rdispl[0] = 0;
  for (int aa = 1; aa < size; aa++)
  {
    rdispl[aa] = rdispl[aa - 1] + all_T_interface_dofs_info[aa - 1];
  }

  // gather global interface dof numberings from other rnaks
  MPI_Allgatherv(GlobalDofNo_interface, N_interface_dofs, MPI_INT, all_GlobalDofNo, all_T_interface_dofs_info, rdispl, MPI_INT, Comm);

  // table is created only over the total own interface dofs over the sub domains
  for (int i = 0; i < size; i++)
  {
    Master_Table[i] = new char[total_interface_dofs];
    memset(Master_Table[i], 'x', total_interface_dofs * sizeof(char));
  }
  
  start = 0;
  end = 0;

  for (int i = 0; i < size; i++)
  {
    end += all_T_interface_dofs_info[i];

    for (int aa = start; aa < end; aa++)
    {
      Master_Table[i][all_GlobalDofNo[aa]] = 'y';
    }

    start = end; 
  }
  
  N_ranks_per_interface_dofs = new int[total_interface_dofs];
  memset(N_ranks_per_interface_dofs, 0, total_interface_dofs * sizeof(int));
  
  for (int j = 0; j < total_interface_dofs; j++)
  {
    // check how many ranks share the interface dof
    for (int i = 0; i < size; i++)
    {
      if (Master_Table[i][j] == 'y')
      {
        N_ranks_per_interface_dofs[j]++;
      }
    }

    // compute the max munber of rank sharing any interface dof

    if (max_n_ranks_interface_dofs < N_ranks_per_interface_dofs[j])
    {
      max_n_ranks_interface_dofs = N_ranks_per_interface_dofs[j];
    }
  }
  
  for (int j = 0; j < total_interface_dofs; j++)
  {
    // a interface dof should have at least one neighbour from other rank
    if (N_ranks_per_interface_dofs[j] == 1)
    {
      ErrThrow("......................This shudnt happen l.783.....j=%d............\n", j);
    }
  }
 
  memset(N_allocated_masters, 0, size * sizeof(int));
  for (int i = 0; i < N_Dof; i++)
  {
    if (DofMarker[i] == 'z')
    {
      N_allocated_masters[Master[i]]++;
    }
  }
  
  memset(N_allocated_masters, 0, size * sizeof(int));
  int min = 0;

  // a interface dof should have at least one neighbour from other rank
  for (int i = 2; i <= max_n_ranks_interface_dofs; i++)
  {
    for (int j = 0; j < total_interface_dofs; j++)
    {
      if (N_ranks_per_interface_dofs[j] == i)
      {
        for (int aa = 0; aa < size; aa++)
        {
          tempc[aa] = Master_Table[aa][j];
        }
  
        min = find_min(N_allocated_masters, size, tempc);
        N_allocated_masters[min]++;
  
        for (int aa = 0; aa < N_Dof; aa++)
        {
          if (GlobalDofNo[aa] == j)
          {
            Master[aa] = min;
          }
        }

        N_ranks_per_interface_dofs[j] = -1;
      }
    }
  }
  
  memset(N_allocated_masters, 0, size * sizeof(int));

  for (int i = 0; i < N_Dof; i++)
  {
    if (DofMarker[i] == 'z')
    {
      N_allocated_masters[Master[i]]++;
      
      if (GlobalDofNo[i] == -1)
      {
        ErrThrow("......................This shudnt happen l.846.................\n");
      }
    }
  }

  delete[] N_Dof_Master;
  delete[] N_Dof_Slave;
  
  for (int i = 0; i < 2; i++)
  {
    delete[] MasterPos[i];
    delete[] SlavePos[i];
  }

  delete[] MasterPos;
  delete[] SlavePos;
  
  delete[] masterInfo;
  delete[] slaveInfo;
  
  delete[] all_GlobalDofNo;
  delete[] GlobalDofNo;
  delete[] GlobalDofNo_interface;

  delete[] all_own_dofs_info;
  delete[] all_interface_dofs_info;
  delete[] all_T_interface_dofs_info;

  for (int i = 0; i < size; i++)
  {
    delete[] Master_Table[i];
  }

  delete[] Master_Table;

  delete[] N_ranks_per_interface_dofs;

  delete[] N_allocated_masters;

  delete[] buffer_index_per_rank;

  delete[] tempc;
}

void TParFEMapper3D::MarkDOFs()
{
  int rank, size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  auto Coll  = FESpace->GetCollection();
  int N_Cells    = Coll->GetN_Cells();
  int N_OwnCells = Coll->GetN_OwnCells();

  //------------------------------------------------------------------//
  /** Gather information about other DOFs ::                          */
  /**    ## dofs marked as 'i'             -->independent             */
  /**    ## dofs marked as 's'             -->slave                   */
  /**    ## dofs marked as 'D'             -->dependent_type1         */
  /**    ## dofs marked as 'd'             -->dependent_type2         */
  /**    ## dofs marked as 'H'             -->halo_type1              */
  /**    ## dofs marked as 'h'             -->halo_type2              */
  //------------------------------------------------------------------// 
  
  // all dofs marked as independent
  memset(DofMarker, 'i', N_Dof * sizeof(char));
  
  // all dofs in halo cell marked as halo_type2 (unused)
  for (int i = N_OwnCells; i < N_Cells; i++)
  {
    auto DOF = FESpace->GetGlobalDOF(i);
    int N_LocDof = FESpace->get_n_local_dof(i);
   
    for (int j = 0; j < N_LocDof; j++)
    {
      DofMarker[DOF[j]] = 'h';
    }
  }
  
  for (int i = 0; i < N_OwnCells; i++)
  {
    auto cell = Coll->GetCell(i);

    // now mark dofs in dependent cells
    if (cell->IsDependentCell())
    {
      auto DOF = FESpace->GetGlobalDOF(i);
      int N_LocDof = FESpace->get_n_local_dof(i);
     
      for (int j = 0; j < N_LocDof; j++)
      {
        int N = DOF[j];

        if (DofMarker[N] != 'h')  // if dof is not marked halo then it is dependent dof
        {
          DofMarker[N] = 'd';
        }
      }
    }
  }
  
  for (int i = 0; i < N_OwnCells; i++)
  {
    auto cell = Coll->GetCell(i);

    // now mark dofs in dependent cells
    if (cell->IsDependentCell())
    {
      auto DOF = FESpace->GetGlobalDOF(i);
      int N_LocDof = FESpace->get_n_local_dof(i);
     
      for (int j = 0; j < N_LocDof; j++)
      {
        int N = DOF[j];

        if (DofMarker[N] == 'h')
        {
          if (Master[N] != rank)
          {
            DofMarker[N] = 's'; // if dof is marked halo and master of it is some other proc then it is a slave
          }
          else
          {
            DofMarker[N] = 'm'; //else it is a master
          }
        }
      }
    }
  }
  
  //mark dependent type1 & type2

  for (int i = 0; i < N_OwnCells; i++)
  {
    auto cell = Coll->GetCell(i);

    // now mark dofs in dependent cells
    if (cell->IsDependentCell())
    {
      bool any_foreign = false;

      auto DOF = FESpace->GetGlobalDOF(i);
      int N_LocDof = FESpace->get_n_local_dof(i);
     
      for (int j = 0; j < N_LocDof; j++)
      {
        int N = DOF[j];

        if (Master[N] != rank)
        {
          any_foreign = true;
          break;
        }
      }
     
      // dependent dofs connected to slave dofs are marked as type1(D) else type2(d)
      // type1 dofs must be smoothed first as they are the halo dofs (type1 i.e. useful) to other procs
      if (any_foreign)
      {
        for (int j = 0; j < N_LocDof; j++)
        {
          int N = DOF[j];

          if (DofMarker[N] == 'd')
          {
            DofMarker[N] = 'D';
          }
        }
      }
    }
  }

#ifdef _OMP
  if (TDatabase::ParamDB->Par_P5 == 1)
  {  
    N_InterfaceM = 0;
    N_Halo1 = 0;
    N_Dept2 = 0;
    N_Int = 0;
  
    for (int i = 0; i < N_Dof; i++)
    {
      switch (DofMarker[i])
      {
        case 'H':
          N_Halo1++;
          break;
        case 'd':
          N_Dept2++;
          break;
        case 'i':
          N_Int++;
          break;
        case 'm':
          N_InterfaceM++;
          break;
        case 's':
          N_InterfaceS++;
          break;
        case 'D':
          N_Dept1++;
          break;
      }
    }
  
    double reqd_ratio = (N_InterfaceM + N_InterfaceS) / (N_Dept1 + N_Halo1);
    double ratio;
    int ctr = 0;
  
    N_Dept2 = 0;
    N_Int = 0;
    for (int i = 0; i < N_Dof; i++)
    {
      if (DofMarker[i] == 'd')
      {
        N_Dept2++;
      }
      else if (DofMarker[i] == 'i')
      {
        N_Int++;
      }
    }
  
    if (N_Dept2 != 0)
    {
      ratio = N_Int / N_Dept2;
    }
    else
    {
      ratio = 0;
    }
    
    while (ratio < reqd_ratio && ctr < 4)
    {
      ctr++;

      for (int i = 0; i < N_OwnCells; i++)
      {
        auto cell = Coll->GetCell(i);

        // now mark dofs in independent cells
        if (!(cell->IsDependentCell()))
        {
          bool any_dependent = false;

          auto DOF = FESpace->GetGlobalDOF(i);
          int N_LocDof = FESpace->get_n_local_dof(i);

          for (int j = 0; j < N_LocDof; j++)
          {
            int N = DOF[j];

            if (DofMarker[N] == 'd' || DofMarker[N] == 'D')
            {
              any_dependent = true;
              break;
            }
          }

          // dependent dofs connected to slave dofs are marked as type1(D) else type2(d)
          // here we increase the ratio of dependent_type2 dofs to independent dofs
          if (any_dependent)
          {
            for (int j = 0; j < N_LocDof; j++)
            {
              int N = DOF[j];

              if (DofMarker[N] == 'i')
              {
                DofMarker[N] = 't';
              }
            }
          }
        }
      }
      
      for (int i = 0; i < N_Dof; i++)
      {
        if (DofMarker[i] == 't')
        {
          DofMarker[i] = 'd';
        }
      }
    
      N_Dept2 = 0;
      N_Int = 0;

      for (int i = 0; i < N_Dof; i++)
      {
        if (DofMarker[i] == 'd')
        {
          N_Dept2++;
        }
        else if (DofMarker[i] == 'i')
        {
          N_Int++;
        }
      }

      if (N_Dept2 != 0)
      {
        ratio = N_Int / N_Dept2;
      }
      else
      {
        ratio = 0;
      }
    }
  }
#endif
  
  // mark halo type1(H)-->useful & type2(h)
  for (int i = N_OwnCells; i < N_Cells; i++)
  {
    bool any_own = false;

    auto DOF = FESpace->GetGlobalDOF(i);
    int N_LocDof = FESpace->get_n_local_dof(i);
   
    for (int j = 0; j < N_LocDof; j++)
    {
      int N = DOF[j];
      if (Master[N] == rank)
      {
        any_own = true;
        break;
      }
    }

    // halo dofs connected to master dofs are marked as type1(H) else type2(h)
    // type2 halo dofs are not required in smoothing operations
    if (any_own)
    {
      for (int j = 0; j < N_LocDof; j++)
      {
        int N = DOF[j];
        if (DofMarker[N] == 'h')
        {
          DofMarker[N] = 'H';  
        }
      }
    }
  }
}

void TParFEMapper3D::MapDOFsMS(int* LocalIndex, int** MappingData)
{
  int rank, size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  auto Coll  = FESpace->GetCollection();
  int N_Cells    = Coll->GetN_Cells();
 
  N_InterfaceM = 0;
  N_InterfaceS = 0;

  N_Slave = 0;
  N_OwnDof = 0; 
  N_Master = 0;
  N_Int = 0; 

  N_Dept = 0;
  N_Dept1 = 0;
  N_Dept2 = 0;

  N_Halo = 0;
  N_Halo1 = 0;
  N_Halo2 = 0;
 
  N_DofSend = new int[size];
  N_DofSendMS = new int[size];
  N_DofSendH2 = new int[size];
  N_DofSendH1 = new int[size];
 
  N_DofRecv = new int[size];
  N_DofRecvMS = new int[size];
  N_DofRecvH2 = new int[size];
  N_DofRecvH1 = new int[size];
  
  memset(N_DofRecvMS, 0, size * sizeof(int));
  memset(N_DofRecvH1, 0, size * sizeof(int));
  memset(N_DofRecvH2, 0, size * sizeof(int));
  
  int* buffer_index_per_rank = new int[size];
  
  int** SlaveBuf = new int*[2];
  int** MasterBuf = new int*[2];

  for (int N = 0; N < N_Dof; N++)
  {
    if (Master[N] != rank)
    {
      N_Slave++;
      N_DofRecv[Master[N]]++;
     
      if (DofMarker[N] == 's')
      {
        N_InterfaceS++;
        N_DofRecvMS[Master[N]]++;
      }
      else if (DofMarker[N] == 'H')
      {
        N_Halo1++;       
        N_DofRecvH1[Master[N]]++;
      }
      else
      {
        N_Halo2++;
        N_DofRecvH2[Master[N]]++;
      }
    }
    else
    {
      N_OwnDof++;
      N_Master++;

      if (DofMarker[N] == 'm')
      {
        N_InterfaceM++;
      }
      else if (DofMarker[N] == 'D')
      {
        N_Dept1++;
      }
      else if (DofMarker[N] == 'd')
      {
        N_Dept2++;
      }
      else
      {
        N_Int++;
      }
    }
  }
 
  N_Halo = N_Halo1 + N_Halo2;
  N_Dept = N_Dept1 + N_Dept2;

  sdispl[0] = 0;
  rdispl[0] = 0;

  sdisplMS = new int[size];
  sdisplMS[0] = 0;  
  rdisplMS = new int[size];
  rdisplMS[0] = 0;

  sdisplH1 = new int[size];
  sdisplH1[0] = 0;  
  rdisplH1 = new int[size];
  rdisplH1[0] = 0;

  sdisplH2 = new int[size];
  sdisplH2[0] = 0;
  rdisplH2 = new int[size];
  rdisplH2[0] = 0;
  
  OwnDofs = new int[N_OwnDof];
  
  N_SendDofMS = 0;
  N_SendDofH1 = 0;
  N_SendDofH2 = 0;

  MPI_Alltoall(N_DofRecv  , 1, MPI_INT, N_DofSend  , 1, MPI_INT, Comm);
  MPI_Alltoall(N_DofRecvMS, 1, MPI_INT, N_DofSendMS, 1, MPI_INT, Comm);
  MPI_Alltoall(N_DofRecvH1, 1, MPI_INT, N_DofSendH1, 1, MPI_INT, Comm);
  MPI_Alltoall(N_DofRecvH2, 1, MPI_INT, N_DofSendH2, 1, MPI_INT, Comm);
  
  for (int i = 1; i < size; i++)
  {
    rdispl[i]   = rdispl  [i - 1] + N_DofRecv  [i - 1];
    rdisplMS[i] = rdisplMS[i - 1] + N_DofRecvMS[i - 1];
    rdisplH1[i] = rdisplH1[i - 1] + N_DofRecvH1[i - 1];
    rdisplH2[i] = rdisplH2[i - 1] + N_DofRecvH2[i - 1];
   
    sdispl[i]   = sdispl  [i - 1] + N_DofSend  [i - 1];
    sdisplMS[i] = sdisplMS[i - 1] + N_DofSendMS[i - 1];
    sdisplH1[i] = sdisplH1[i - 1] + N_DofSendH1[i - 1];
    sdisplH2[i] = sdisplH2[i - 1] + N_DofSendH2[i - 1];
  }
 
  for (int i = 0; i < size; i++)
  {
    N_SendDofMS += N_DofSendMS[i];
    N_SendDofH1 += N_DofSendH1[i];
    N_SendDofH2 += N_DofSendH2[i];
  }

  N_SendDof = N_SendDofMS + N_SendDofH1 + N_SendDofH2;

  memcpy(buffer_index_per_rank, rdisplMS, size * sizeof(int));
 
  DofSend    = new int[N_SendDofMS + N_SendDofH1 + N_SendDofH2];
  DofSendMS  = DofSend;
  DofSendH1  = DofSend + N_SendDofMS;
  DofSendH2  = DofSend + N_SendDofMS + N_SendDofH1;
  
  DofRecv    = new int[N_InterfaceS + N_Halo1 + N_Halo2];
  DofRecvMS  = DofRecv;
  DofRecvH1  = DofRecv + N_InterfaceS;
  DofRecvH2  = DofRecv + N_InterfaceS + N_Halo1;
  
  for (int i = 0; i < 2; i++)
  {
    if (N_InterfaceS > 0)
    {
      SlaveBuf[i] = new int[N_InterfaceS];
    }

    if (N_SendDofMS > 0)
    {
      MasterBuf[i] = new int[N_SendDofMS];
    }
  }
  
  Reorder = new int[N_Dof];
  NewGN = new int[N_Dof];
  int m = 0;
 
  int Mstr  = 0;
  Reorder_M  = Reorder;

  int Indpt = N_InterfaceM;
  Reorder_I  = Reorder + Indpt;

  int Dept1 = Indpt + N_Int;
  Reorder_D1 = Reorder + Dept1;

  int Dept2 = Dept1 + N_Dept1;
  Reorder_D2 = Reorder + Dept2;

  int Slv = Dept2 + N_Dept2;
  int Hl1 = Slv + N_InterfaceS;
  int Hl2 = Hl1 + N_Halo1;

  int ts = 0, th1 = 0, th2 = 0, ti = 0, tm = 0, td1 = 0, td2 = 0;
  
  for (int i = 0; i < N_Dof; i++)
  {
    // slave dofs
    if (Master[i] != rank)
    {
      if (DofMarker[i] == 's')
      {
        SlaveBuf[GLOBAL_NO][buffer_index_per_rank[Master[i]]] = MappingData[GLOBAL_NO][i];
        SlaveBuf[DOF_NO]   [buffer_index_per_rank[Master[i]]] = MappingData[DOF_NO]   [i];
   
        DofRecvMS[buffer_index_per_rank[Master[i]]] = i;
        buffer_index_per_rank[Master[i]]++;
         
        Reorder[Slv] = i;
        NewGN[i]     = ts++;
        Slv++;
      }
      else if (DofMarker[i] == 'H')
      {
        Reorder[Hl1] = i;
        NewGN[i]     = th1++;
        Hl1++;
      }
      else
      {
        Reorder[Hl2] = i;
        NewGN[i]     = th2++;
        Hl2++;
      }
    }
    else
    {
      OwnDofs[m++] = i;

      if (DofMarker[i] == 'm')
      {
        Reorder[Mstr] = i;
        NewGN[i]      = tm++;
        Mstr++;
      }
      else if (DofMarker[i] == 'D')
      {
        Reorder[Dept1] = i;
        NewGN[i]       = td1++;
        Dept1++;
      }
      else if (DofMarker[i] == 'd')
      {
        Reorder[Dept2] = i;
        NewGN[i]       = td2++;
        Dept2++;
      }
      else
      {
        Reorder[Indpt] = i;
        NewGN[i]       = ti++;
        Indpt++;
      }
    }
  }
 
  MPI_Alltoallv(SlaveBuf[GLOBAL_NO], N_DofRecvMS, rdisplMS, MPI_INT, MasterBuf[GLOBAL_NO], N_DofSendMS, sdisplMS, MPI_INT, Comm); 
  MPI_Alltoallv(SlaveBuf[DOF_NO],    N_DofRecvMS, rdisplMS, MPI_INT, MasterBuf[DOF_NO],    N_DofSendMS, sdisplMS, MPI_INT, Comm);
 
  for (int i = 0; i < N_SendDofMS; i++)
  {
    int temp_globalno = MasterBuf[GLOBAL_NO][i];
    int temp_dofno    = MasterBuf[DOF_NO][i];
    int temp = GetLocalIndex(N_Cells, LocalIndex, temp_globalno);

    DofSendMS[i] = FESpace->get_global_dof(temp, temp_dofno); 
  }

  delete[] buffer_index_per_rank;

  for (int i = 0; i < 2; i++)
  {
    if (N_InterfaceS > 0)
    {
      delete[] SlaveBuf[i];
    }

    if (N_SendDofMS > 0)
    {
      delete[] MasterBuf[i];
    }
  }

  delete[] SlaveBuf;
  delete[] MasterBuf;
}

void TParFEMapper3D::MapDOFsHalo(int* LocalIndex, int** MappingData)
{
  int rank, size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  auto Coll   = FESpace->GetCollection();
  int N_Cells = Coll->GetN_Cells();

  int **SlaveBufH1, **MasterBufH1;
  int **SlaveBufH2, **MasterBufH2;

  int* buffer_index_per_rank_H1 = new int[size];
  int* buffer_index_per_rank_H2 = new int[size];
 
  SlaveBufH1  = new int*[2];
  MasterBufH1 = new int*[2];
 
  SlaveBufH2  = new int*[2];
  MasterBufH2 = new int*[2];
 
  for (int i = 0; i < 2; i++)
  {
    if (N_Halo1 > 0)
    {
      SlaveBufH1[i] = new int[N_Halo1];

      for (int indx = 0; indx < N_Halo1; ++indx)
      {
        SlaveBufH1[i][indx] = 0;
      }
    }

    if (N_SendDofH1 > 0)
    {
      MasterBufH1[i] = new int[N_SendDofH1];

      for (int indx = 0; indx < N_SendDofH1; ++indx)
      {
        MasterBufH1[i][indx] = 0;
      }
    }

    if (N_Halo2 > 0)
    {
      SlaveBufH2[i] = new int[N_Halo2];
      
      for (int indx = 0; indx < N_Halo2; ++indx)
      {
        SlaveBufH2[i][indx] = 0;
      }
    }

    if (N_SendDofH2 > 0)
    {
      MasterBufH2[i] = new int[N_SendDofH2];

      for (int indx = 0; indx < N_SendDofH2; ++indx)
      {
        MasterBufH2[i][indx] = 0;
      }
    }
  }

  memcpy(buffer_index_per_rank_H1, rdisplH1, size * sizeof(int));
  memcpy(buffer_index_per_rank_H2, rdisplH2, size * sizeof(int));
  
  for (int i = 0; i < N_Dof; i++)
  {
    if (DofMarker[i] == 'H')
    {
      if (Master[i] == rank)
      {
        ErrThrow("This should be a HALO_1 dof");
      }
     
      SlaveBufH1[GLOBAL_NO][buffer_index_per_rank_H1[Master[i]]] = MappingData[GLOBAL_NO][i];
      SlaveBufH1[DOF_NO]   [buffer_index_per_rank_H1[Master[i]]] = MappingData[DOF_NO]   [i];
    
      DofRecvH1[buffer_index_per_rank_H1[Master[i]]] = i;
      buffer_index_per_rank_H1[Master[i]]++;
    }
    else if(DofMarker[i] == 'h')
    {
      if (Master[i] == rank)
      {
        ErrThrow("This should be a HALO_2 dof");
      }
     
      SlaveBufH2[GLOBAL_NO][buffer_index_per_rank_H2[Master[i]]] = MappingData[GLOBAL_NO][i];
      SlaveBufH2[DOF_NO]   [buffer_index_per_rank_H2[Master[i]]] = MappingData[DOF_NO]   [i];
    
      DofRecvH2[buffer_index_per_rank_H2[Master[i]]] = i;
      buffer_index_per_rank_H2[Master[i]]++;
    }
  }
  
  MPI_Alltoallv(SlaveBufH1[GLOBAL_NO], N_DofRecvH1, rdisplH1, MPI_INT, MasterBufH1[GLOBAL_NO], N_DofSendH1, sdisplH1, MPI_INT, Comm); 
  MPI_Alltoallv(SlaveBufH1[DOF_NO],    N_DofRecvH1, rdisplH1, MPI_INT, MasterBufH1[DOF_NO],    N_DofSendH1, sdisplH1, MPI_INT, Comm);
 
  for (int i = 0; i < N_SendDofH1; i++)
  {
    int temp_globalno = MasterBufH1[GLOBAL_NO][i];
    int temp_dofno    = MasterBufH1[DOF_NO][i];
    int temp          = GetLocalIndex(N_Cells, LocalIndex, temp_globalno);

    DofSendH1[i] = FESpace->get_global_dof(temp, temp_dofno); 
  }
  
  MPI_Alltoallv(SlaveBufH2[GLOBAL_NO], N_DofRecvH2, rdisplH2, MPI_INT, MasterBufH2[GLOBAL_NO], N_DofSendH2, sdisplH2, MPI_INT, Comm);
  MPI_Alltoallv(SlaveBufH2[DOF_NO],    N_DofRecvH2, rdisplH2, MPI_INT, MasterBufH2[DOF_NO],    N_DofSendH2, sdisplH2, MPI_INT, Comm);

  for (int i = 0; i < N_SendDofH2; i++)
  {
    int temp_globalno = MasterBufH2[GLOBAL_NO][i];
    int temp_dofno    = MasterBufH2[DOF_NO][i];
    int temp          = GetLocalIndex(N_Cells, LocalIndex, temp_globalno);

    DofSendH2[i] = FESpace->get_global_dof(temp, temp_dofno); 
  }
 
  for (int i = 0; i < 2; i++)
  {
    if (N_Halo1 > 0)
    {
      delete[] SlaveBufH1[i];
    }

    if (N_SendDofH1 > 0)
    {
      delete[] MasterBufH1[i];
    }

    if (N_Halo2 > 0)
    {
      delete[] SlaveBufH2[i];
    }

    if (N_SendDofH2 > 0)
    {
      delete[] MasterBufH2[i];
    }
  }

  delete[] SlaveBufH1;
  delete[] MasterBufH1;

  delete[] SlaveBufH2;
  delete[] MasterBufH2;

  delete[] buffer_index_per_rank_H1;
  delete[] buffer_index_per_rank_H2;
}

void TParFEMapper3D::ConstructDofMap_Master_Halo()
{
  double start_time, end_time, temp_time;

  int rank, size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  start_time = MPI_Wtime();
  temp_time  = start_time;
  
  auto Coll      = FESpace->GetCollection();
  int N_Cells    = Coll->GetN_Cells();
  int N_OwnCells = Coll->GetN_OwnCells();

  N_Dof = FESpace->get_n_dof();

  // Check whether the Collection is in proper order: own cells first, then halo cells.
  for (int i = 0; i < N_OwnCells; i++)
  {
    auto cell = Coll->GetCell(i);

    if (cell->IsHaloCell())
    {
      Output::print("Halo cell at position ", i, " although position 0 to ",
                    N_OwnCells - 1, " are reserved for own cells.");
    }
  }

  for (int i = N_OwnCells; i < N_Cells; i++)
  {
    auto cell = Coll->GetCell(i);

    if (!cell->IsHaloCell())
    {
      Output::print("Non-halo cell at position ", i, " although position ", N_OwnCells, " to ",
                    N_Cells - 1, " are reserved for halo cells.");
    }
  }

  // Array containing global number of all Local cells
  int* LocalIndex = new int[N_Cells];
  
  // For each DOF of this space on this processor, MappingData stores
  // - the global cell number (array 0)
  // - the (cell-)local dof number (array 1)

  int **MappingData;
  MappingData = new int*[2];
  for (int i = 0; i < 2; i++)
  {
    MappingData[i] = new int[N_Dof];
    memset(MappingData[i], 0, N_Dof * sizeof(int));
  }
  
  /** *************************************************************************************/
  /** Master   :: Array containing the rank to which the Dof belongs                      */ 
  /** DofMarker:: Array to help finding Master, Slave, Dependent, Halo                    */
  /** DofMarker:: ALL the Dofs are marked 'i'                                             */
  /** *************************************************************************************/
  Master = new int[N_Dof];
  
  // reserve all DOFs for this process until we look more closely later
  for (int i = 0; i < N_Dof; i++) 
  {
    Master[i] = rank;
  }

  DofMarker = new char[N_Dof];
  memset(DofMarker, 'i', N_Dof*sizeof(char));
  
  /** *************************************************************************************/
  /** Local Index Array is filled with GlobalNumbers of each cell                         */ 
  /** DofMarker::Dofs of Dependent cells are marked 'd'                                   */
  /** *************************************************************************************/
  for (int i = 0; i < N_Cells; i++)
  {
    auto cell = Coll->GetCell(i);
    LocalIndex[i] = cell->GetGlobalCellNo();

    if (cell->IsDependentCell() && !cell->IsHaloCell())
    {
      auto DOF = FESpace->GetGlobalDOF(i);
      int N_LocDof = FESpace->get_n_local_dof(i);

      for (int j = 0; j < N_LocDof; j++)
      {
        DofMarker[DOF[j]] = 'd';
      }
    }
  }
    
  /** *************************************************************************************/
  /** MappingData for dofs to be updated across domains are stored                        */ 
  /** DofMarker::Dofs of Halo cells not a part of dependent cells are marked 'h'          */
  /** *************************************************************************************/  
  for (int i = N_OwnCells; i < N_Cells; i++)
  {
    auto cell = Coll->GetCell(i);
    int ID = cell->GetSubDomainNo();

    auto DOF = FESpace->GetGlobalDOF(i);
    int N_LocDof = FESpace->get_n_local_dof(i);
    
    for (int j = 0; j < N_LocDof; j++)
    {
      int N = DOF[j];

      // DOFs in halo cells are given the master rank of the lowest adjacent cell this process knows about
      
      if (DofMarker[N] == 'i')
      {
        DofMarker[N] = 'h';
        Master[N]    = ID;
      }
      
      if (ID < Master[N])
      {
        Master[N] = ID;
      }
      
      MappingData[GLOBAL_NO][N] = cell->GetGlobalCellNo();
      MappingData[DOF_NO][N]    = j; 
    }
  }

  // these are protected members for historical reasons, so we initialize them here
  sdispl = new int[size]; // offsets into verification data send buffer
  rdispl = new int[size]; // offsets into verification data receive buffer
  
  NegotiateDOFOwnership(LocalIndex, MappingData);
  
  end_time = MPI_Wtime();
  
  if(rank == 0)
  {
    Output::print("total time taken for master verification = ", end_time - temp_time);
  }
  
  temp_time = end_time;
  //#################################################################  Master verification  ##################################################################################//   

  //################################################################# Redistribution of interface dofs #######################################################################//

  // NOTE: it's unclear why this is labelled "redistribution" - what actually
  // happens here is a series of optional sanity checks regarding which process
  // (rank) is assigned ownership of which DOFs. this is purely a debug feature


  if (TDatabase::ParamDB->Par_P4)
  {
    CheckDOFValidity(LocalIndex, MappingData);
    
    end_time = MPI_Wtime();
    if (rank == 0)
    {
      Output::print("Total Time Taken for Redistribution of master dofs = ", end_time - temp_time);
    }
    temp_time = end_time;
  }

  //################################################################# Redistribution of interface dofs #######################################################################//


  
  //#################################################################  Marking The Dofs ######################################################################################//     
  
  MarkDOFs();

  end_time = MPI_Wtime();
  if (rank == 0)
  {
    Output::print("Total Time Taken for marking the dofs = ", end_time - temp_time);
  }
  temp_time = end_time;

  //#################################################################  Marking The Dofs ######################################################################################// 
  
  //#############################################################  Mapper for Master Dof are set ##############################################################################//  
  
  MapDOFsMS(LocalIndex, MappingData);
  
  //*********************************************************************************************************************************//
  //                                                  Mapper for Halo Dof are set                                                    //
  //*********************************************************************************************************************************//
  
  MapDOFsHalo(LocalIndex, MappingData);
  
  if (N_Dim > 1)
  {
    for (int i = 0; i < size; i++)
    {
      sdispl[i]   *= N_Dim;
      sdisplMS[i] *= N_Dim;
      sdisplH1[i] *= N_Dim;
      sdisplH2[i] *= N_Dim;

      rdispl[i]     *= N_Dim;
      rdisplMS[i]   *= N_Dim;
      rdisplH1[i]   *= N_Dim;
      rdisplH2[i]   *= N_Dim;

      N_DofSend[i]   *= N_Dim;
      N_DofSendMS[i] *= N_Dim;
      N_DofSendH1[i] *= N_Dim;
      N_DofSendH2[i] *= N_Dim;

      N_DofRecv[i]    *= N_Dim;
      N_DofRecvMS[i]  *= N_Dim;
      N_DofRecvH1[i]  *= N_Dim;
      N_DofRecvH2[i]  *= N_Dim;
    }
  }
  
  if (N_SendDof > 0)
  {
    Send_Info = new double[N_SendDof * N_Dim];
  }

  Send_InfoMS = Send_Info;
  Send_InfoH1 = Send_Info + N_SendDofMS * N_Dim;
  Send_InfoH2 = Send_Info + N_SendDofMS * N_Dim + N_SendDofH1 * N_Dim;
 
  if (N_Slave > 0)
  {
    Recv_Info = new double[N_Slave * N_Dim];
  }

  Recv_InfoMS = Recv_Info;
  Recv_InfoH1 = Recv_Info + N_InterfaceS * N_Dim;
  Recv_InfoH2 = Recv_Info + N_InterfaceS * N_Dim + N_Halo1 * N_Dim;
 
  end_time = MPI_Wtime();
  if (rank == 0)
  {
    Output::print("Total Time Taken for mapping the dofs = ", end_time - temp_time);
  }
  
  if (rank == TDatabase::ParamDB->Par_P0)
  {
    Output::print("################       Mapping for slave-master dofs and halo_1, halo_2 dofs done !!!    ################");
  }
  
  if (rank == 0)
  {
    Output::print("total time taken by the ConstructDofMap_light() = ", MPI_Wtime() - start_time);  
  }

  /*// only for checking purpose   
  if (0)
  {
    int Total_Own = 0;

    for (int i = 0; i < N_Dof; i++)
    {
      if (Master[i] != rank)
      {
        if (!(DofMarker[i] == 's' || DofMarker[i] == 'h' || DofMarker[i] == 'H'))
        {
          ErrThrow("................1.This shudnt happen l.1527...................");
        }
      }
      else
      {
        Total_Own++;
        if (!(DofMarker[i] == 'm' || DofMarker[i] == 'd' || DofMarker[i] == 'D' || DofMarker[i] == 'i' || DofMarker[i] == 'x'))
        {
          ErrThrow("................2.This shudnt happen l.1535...........DofMarker[%d]=%c........", i, DofMarker[i]);
        }
      }
      else
      {
        ErrThrow("................3.This shudnt happen l.1541...................");
      }
    }
  }*/

  delete[] LocalIndex;

  for (int i = 0; i < 2; i++)
  {
    delete[] MappingData[i];
  }

  delete[] MappingData;
}

void TParFEMapper3D::ConstructDofMap()
{
   
 int rank, size, i, j, m, N_Cells, N_U, N_LocDof{0}, ID;
 int N;
 int *LocalIndex,*Verify;
 int N_OwnCells;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  auto Coll = FESpace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  N_U = FESpace->get_n_dof();
  N_OwnCells = Coll->GetN_OwnCells();
  
 // printf("N_U is %d----Rank %d\n",N_U,rank);
  
  
  for(i=0;i<N_OwnCells;i++)
  {
    auto cell = Coll->GetCell(i);
    if(cell->IsHaloCell()) printf("This shudnt happen l.1587---------------------------------------------\n");
  }
  for(i=N_OwnCells;i<N_Cells;i++)
  {
    auto cell = Coll->GetCell(i);
    if(!cell->IsHaloCell()) printf("This shudnt happen l.1592---------------------------------------------\n");
  }
  /** Array containing number of ranks(process ID) associated (surrounding) with each DOF (including halo DOF) */ 
//   N_DofRankIndex = new int[N_U];
 
//   memset(N_DofRankIndex, 1, N_U*sizeof(int));
  /** Array containing global number of all Local cells */
  LocalIndex = new int[N_Cells]; 
  
  int **MappingData;
  MappingData = new int*[2];
  for(i=0;i<2;i++)
  {
    MappingData[i] = new int[N_U];
    memset(MappingData[i], 0, N_U*sizeof(int));
  }
  
  Master = new int[N_U];
  for(i=0;i<N_U;i++) Master[i]=rank;
    
  Verify = new int[N_U];
  memset(Verify, 0, N_U*sizeof(int));
  
 /** START ---> [ DofRankIndex ------- N_DofRankIndex -------- N_DependentCells ] **/ 
 for(i=0; i<N_Cells; i++)
  {
     auto cell = Coll->GetCell(i);
     LocalIndex[i] = cell->GetGlobalCellNo();
     
     if(cell->IsDependentCell() && !cell->IsHaloCell())
     {
      auto DOF = FESpace->GetGlobalDOF(i);
      N_LocDof = FESpace->get_n_local_dof(i);

      for(j=0; j<N_LocDof; j++)
      {
  N = DOF[j];
  Verify[N] = -1;
      }  // for(j=0; j<N_DOF; j++)
      
     }
  } //  for(i=0; i<N_Cells
    

 for(i=N_OwnCells;i<N_Cells;i++)
 {
      auto cell = Coll->GetCell(i);
      ID = cell->GetSubDomainNo();
      auto DOF = FESpace->GetGlobalDOF(i);
      N_LocDof = FESpace->get_n_local_dof(i);

      for(j=0; j<N_LocDof; j++)
      {
  N = DOF[j];
  //if(N>=N_Active) continue;
  if(Verify[N] == 0) 
  {
    Verify[N] = 1;
    Master[N] = ID;
  }
  if(ID < Master[N]) Master[N] = ID;
  MappingData[GLOBAL_NO][N]    = cell->GetGlobalCellNo();
  MappingData[DOF_NO][N]       = j; 
      }
      
  }
 
   
    int temp,temp_globalno,temp_dofno;
  
  
   /** END ---> [ Master ------------ Verify] **/ 
  
    int *N_DOFtobeverified;    /** Array containing how many DOF's need to be verified from other processors(ranks) */
    int *N_DOFverifiedbythisrank; /** Array containing how many DOF's need to be verified by this rank for other processors(ranks) */
    int VerifiedDOF = 0;    /** Total no of DOF's that need to be verified from other ranks */
    int Verifiedbythisrank =0;    /** Total no of DOF verified by this rank */
    int **sendbuf,**recvbuf,*verrecvbuf,*versendbuf;
     
    /** PROTECTED VARIABLES  * int *sdispl,*rdispl;*/
    int *temp_arr;
  
    sendbuf   = new int*[2];
    recvbuf   = new int*[2];
    sdispl    = new int[size];
    rdispl    = new int[size];
    temp_arr  = new int[size];
  
    N_DOFtobeverified = new int[size];
    memset(N_DOFtobeverified,0,size*sizeof(int));
    N_DOFverifiedbythisrank = new int[size];
    memset(N_DOFtobeverified,0,size*sizeof(int));
  
    for(i=0;i<N_U;i++)
    {
      if(Verify[i] == 1)
      {
  N_DOFtobeverified[Master[i]]++;
  VerifiedDOF++;
      }
    }
    
    MPI_Alltoall(N_DOFtobeverified,1, MPI_INT,N_DOFverifiedbythisrank ,1, MPI_INT, Comm);
          
    for(i=0;i<size;i++)
      Verifiedbythisrank +=  N_DOFverifiedbythisrank[i];
    
    
    for(i=0;i<2;i++)
    {
      sendbuf[i] = new int[VerifiedDOF];
      recvbuf[i] = new int[Verifiedbythisrank];
    }
    
    //printf("Rank F %d----------%d\n",rank,Verifiedbythisrank);
    sdispl[0] = 0;
    for(i=1;i<size;i++)
    {
      sdispl[i] = N_DOFtobeverified[i - 1] + sdispl[i - 1];
//       if(rank == 0)
      //printf(" %d \n",sdispl[i]);
    }
          
    rdispl[0] = 0;
    for(i=1;i<size;i++)
      rdispl[i] = N_DOFverifiedbythisrank[i - 1] + rdispl[i - 1];
          
          
    memcpy (temp_arr,sdispl, size*sizeof(int) );
          
    for(i=0;i<N_U;i++)
     {
      if(Verify[i] == 1)
      {
  //printf("%d ",temp_arr[Master[i]]);
  sendbuf[GLOBAL_NO][temp_arr[Master[i]]] = MappingData[GLOBAL_NO][i];
  sendbuf[DOF_NO][temp_arr[Master[i]]] = MappingData[DOF_NO][i];
  temp_arr[Master[i]]++;
      }
    }

    MPI_Alltoallv(sendbuf[GLOBAL_NO],N_DOFtobeverified, sdispl, MPI_INT,recvbuf[GLOBAL_NO], N_DOFverifiedbythisrank, rdispl, MPI_INT, Comm);
    MPI_Alltoallv(sendbuf[DOF_NO],N_DOFtobeverified, sdispl, MPI_INT,recvbuf[DOF_NO], N_DOFverifiedbythisrank, rdispl, MPI_INT, Comm);
        
    verrecvbuf    = new int[Verifiedbythisrank];
    versendbuf    = new int[VerifiedDOF];

    for(i=0;i<Verifiedbythisrank;i++)
    {
      temp = GetLocalIndex(N_Cells,LocalIndex,recvbuf[GLOBAL_NO][i]);
      verrecvbuf[i] = Master[FESpace->get_global_dof(temp, recvbuf[DOF_NO][i])];
    }
  
    MPI_Alltoallv(verrecvbuf, N_DOFverifiedbythisrank, rdispl, MPI_INT,versendbuf,N_DOFtobeverified, sdispl, MPI_INT, Comm);
          
    for(i=0;i<VerifiedDOF;i++)
    {
      temp_globalno = sendbuf[GLOBAL_NO][i];
      temp_dofno    = sendbuf[DOF_NO][i];
      temp = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
      temp = FESpace->get_global_dof(temp, temp_dofno);
      if(Verify[temp] != 1)
  printf("Error : This degree of Freedom (%d) didn't require verification\n",temp);
      Master[temp] = versendbuf[i];
    }
          
/**===================================== END : MASTER FOR ALL DEGREES OF FREEDOM HAS BEEN DECIDED ========================*/

    delete [] verrecvbuf;
    delete [] versendbuf;
    delete [] N_DOFtobeverified;
    delete [] N_DOFverifiedbythisrank;
    delete [] temp_arr;

    
    for(i=0;i<2;i++)
    {
      delete [] sendbuf[i];
      delete [] recvbuf[i];
    }

    delete [] sendbuf;
    delete [] recvbuf;
    
        
/**==================================START : NOW FIRST GATHERING INFORMATION REGARDING WHICH 
 *                                    DOF's have to be received from OTHER PROCESSORS==================*/
          
    int **SlaveBuf,**MasterBuf,SizeDofSend = 0;
    /** PROTECTED VARIABLES * int *N_DofSend,*N_DofRecv,*DofSend,*DofRecv;*/
    N_Slave = 0;
    N_OwnDof = 0;
          
    N_DofSend = new int[size];
    N_DofRecv = new int[size];
    temp_arr  = new int[size];
    SlaveBuf  = new int*[2];
    MasterBuf = new int*[2];
  
  //for(i=N_Active;i<N_U;i++) Master[i]=rank;
  
    memset(N_DofRecv,0,size*sizeof(int));
      
    for(i=0;i<N_U;i++)
    {
     if(Master[i] != rank)
     {  
      N_Slave++;
      N_DofRecv[Master[i]]++;
     } 
     else 
       N_OwnDof++;
    }
  //  printf("Rank %d      OwnDofs is %d---------\n",rank,N_OwnDof);

    OwnDofs = new int[N_OwnDof];
    rdispl[0] = 0; 

    MPI_Alltoall(N_DofRecv,1, MPI_INT,N_DofSend ,1, MPI_INT, Comm);

    for(i=1;i<size;i++)
      rdispl[i] = rdispl[i - 1] + N_DofRecv[i - 1];
  
    sdispl[0] = 0;
    for(i=0;i<size;i++)
    {
      SizeDofSend += N_DofSend[i];
      sdispl[i] = sdispl[i - 1] + N_DofSend[i - 1];
    }


    DofSend  = new int[SizeDofSend];
    DofRecv  = new int[N_Slave];
    N_SendDof = SizeDofSend;
  
    for(i=0;i<2;i++)
    {
      SlaveBuf[i]  = new int[N_Slave];
      MasterBuf[i] = new int[SizeDofSend];
    }
  
    memcpy (temp_arr,rdispl, size*sizeof(int) );
  
    m=0;
    for(i=0;i<N_U;i++)
    {
      if(Master[i] != rank)
      {
  SlaveBuf[GLOBAL_NO][temp_arr[Master[i]]] = MappingData[GLOBAL_NO][i];
  SlaveBuf[DOF_NO]   [temp_arr[Master[i]]] = MappingData[DOF_NO]   [i];
  DofRecv[temp_arr[Master[i]]] = i;
  temp_arr[Master[i]]++;
      }
      else 
      {
  OwnDofs[m]=i;
  m++;
      }
    }
      
    MPI_Alltoallv(SlaveBuf[GLOBAL_NO],N_DofRecv, rdispl, MPI_INT,MasterBuf[GLOBAL_NO], N_DofSend, sdispl, MPI_INT, Comm);
    MPI_Alltoallv(SlaveBuf[DOF_NO],N_DofRecv, rdispl, MPI_INT,MasterBuf[DOF_NO], N_DofSend, sdispl, MPI_INT, Comm);
    
    for(i=0;i<SizeDofSend;i++)
    {
      temp_globalno = MasterBuf[GLOBAL_NO][i];
      temp_dofno    = MasterBuf[DOF_NO][i];
      temp = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
      temp = FESpace->get_global_dof(temp, + temp_dofno);
      DofSend[i] = temp;  
    }
  
  //printf(" Rank %d ------ NUMBER OF DOF's to be sent = %d -------- NUMBER OF DOF's to be recv = %d\n",rank,N_SendDof,N_Slave); 
    for(i=0;i<2;i++)
    {
      delete [] SlaveBuf[i];
      delete [] MasterBuf[i];
    }

    delete [] SlaveBuf;
    delete [] MasterBuf;
    
   if(N_Dim>1)
   {
     for(i=0;i<size;i++)
     {
        sdispl[i]     *= N_Dim;
  rdispl[i]     *= N_Dim;
  
  N_DofSend[i]  *= N_Dim;
        N_DofRecv[i]  *= N_Dim; 
     }
   }
  
  if(N_SendDof>0)
   Send_Info = new double[N_SendDof*N_Dim];
  if(N_Slave>0)
   Recv_Info = new double[N_Slave*N_Dim];
    
  if(rank==TDatabase::ParamDB->Par_P0)
  printf("ConstructDofMap done !!!\n");
      
       //DEBUG
 int aa;
 if(size<9)
 for(aa=0;aa<size;aa++){
     if(rank==aa){
          printf("\nRank::%d\n",rank);
    printf("N_OwnDof     = %d\t Not_Own = %d\n",N_OwnDof,N_Slave);
    }
     MPI_Barrier(MPI_COMM_WORLD);
  }//verified
  
//    printf("N_SendDofMS = %d\n",N_SendDofMS);
}


void TParFEMapper3D::Assign_GlobalDofNo()
{ 
  int rank, size;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
  
  int *Send_data = new int[N_SendDof];
  int *Recv_data = new int[N_Slave];
  
  // these are offset pointers into Send_data and Recv_data
  int *Send_dataMS, *Send_dataH1, *Send_dataH2;
  int *Recv_dataMS, *Recv_dataH1, *Recv_dataH2;

  double t1 = MPI_Wtime();
  int *all_own_dofs_info = new int[size];

  // initialize local-global mapping with -1
  Local2Global = new int[N_Dof];
  for (int i = 0; i < N_Dof; i++)
  {
    Local2Global[i] = -1;
  }
  
  // communicate dof counts to each other rank
  MPI_Allgather(&N_OwnDof, 1, MPI_INT, all_own_dofs_info, 1, MPI_INT, Comm);
  
  // start is the global index of our first dof
  int start = 0;
  for (int k = 0; k < rank; k++)
  {
    start += all_own_dofs_info[k];
  }
  
  // globally number our dofs
  for (int i = 0; i < N_Dof; i++)
  {
    if (Master[i] == rank)
    {
      Local2Global[i] = start++;
    }
  }
  
  // these offsets and counts were originally multiplied by N_Dim
  if (N_Dim > 1)
  {
    for (int i = 0; i < size; i++)
    {
      sdispl[i]   /= N_Dim;
      sdisplMS[i] /= N_Dim;
      sdisplH1[i] /= N_Dim;
      sdisplH2[i] /= N_Dim;

      rdispl[i]   /= N_Dim;
      rdisplMS[i] /= N_Dim;
      rdisplH1[i] /= N_Dim;
      rdisplH2[i] /= N_Dim;

      N_DofSend[i]   /= N_Dim;
      N_DofSendMS[i] /= N_Dim;
      N_DofSendH1[i] /= N_Dim;
      N_DofSendH2[i] /= N_Dim;

      N_DofRecv[i]   /= N_Dim;
      N_DofRecvMS[i] /= N_Dim;
      N_DofRecvH1[i] /= N_Dim;
      N_DofRecvH2[i] /= N_Dim;
    }
  }
   
  if (TDatabase::ParamDB->MapperType != 2)
  {
    // "new" mapper with halo cells

    // compute offset send pointers
    Send_dataMS = Send_data;
    Send_dataH1 = Send_data + N_SendDofMS;
    Send_dataH2 = Send_data + N_SendDofMS + N_SendDofH1;
  
    // compute offset receive pointers
    Recv_dataMS = Recv_data;
    Recv_dataH1 = Recv_data + N_InterfaceS;
    Recv_dataH2 = Recv_data + N_InterfaceS + N_Halo1;

    // copy global numbering of our MS dofs to send
    for (int i = 0; i < N_SendDofMS; i++)
    {
      Send_dataMS[i] = Local2Global[DofSendMS[i]];
    }

    // communicate global numbering of MS dofs
    MPI_Alltoallv(Send_dataMS, N_DofSendMS, sdisplMS, MPI_INT,
                  Recv_dataMS, N_DofRecvMS, rdisplMS, MPI_INT, Comm);
  
    // write out received global numbering of foreign MS dofs
    for (int i = 0; i < N_InterfaceS; i++)
    {
      Local2Global[DofRecvMS[i]] = Recv_dataMS[i];
    }

    // copy global numbering of our Halo 1 dofs to send
    for (int i = 0;i < N_SendDofH1; i++)
    {
      Send_dataH1[i] = Local2Global[DofSendH1[i]];
    }

    // communicate global numbering of Halo 1 dofs
    MPI_Alltoallv(Send_dataH1, N_DofSendH1, sdisplH1, MPI_INT,
                  Recv_dataH1, N_DofRecvH1, rdisplH1, MPI_INT, Comm);
  
    // write out received global numbering of foreign Halo 1 dofs
    for (int i = 0; i < N_Halo1; i++)  
    {
      Local2Global[DofRecvH1[i]] = Recv_dataH1[i];
    }

    // copy global numbering of our Halo 2 dofs to send
    for (int i = 0; i < N_SendDofH2; i++)
    {
      Send_dataH2[i] = Local2Global[DofSendH2[i]];
    }
  
    // communicate global numbering of Halo 2 dofs
    MPI_Alltoallv(Send_dataH2, N_DofSendH2, sdisplH2, MPI_INT,
                  Recv_dataH2, N_DofRecvH2, rdisplH2, MPI_INT, Comm);
  
    // write out received global numbering of foreign Halo 2 dofs
    for (int i = 0; i < N_Halo2; i++)  
    {
      Local2Global[DofRecvH2[i]] = Recv_dataH2[i];
    }
  }
  else
  {
    // old mapper, unclear whether it's worth preserving

    int *Send_data = new int[N_SendDof];
    int *Recv_data = new int[N_Slave];
    
    for (int i = 0; i < N_SendDof; i++)
    {
      Send_data[i] = Local2Global[DofSend[i]];
    }
  
    MPI_Alltoallv(Send_data, N_DofSend, sdispl, MPI_INT,
                  Recv_data, N_DofRecv, rdispl, MPI_INT, Comm);
  
    for (int i = 0; i < N_Slave; i++)  
    {
      Local2Global[DofRecv[i]] = Recv_data[i];
    }
  }
  
  // re-multiply offsets and counts with N_Dim
  if (N_Dim > 1)
  {
    for (int i = 0; i < size; i++)
    {
      sdispl[i]   *= N_Dim;
      sdisplMS[i] *= N_Dim;
      sdisplH1[i] *= N_Dim;
      sdisplH2[i] *= N_Dim;

      rdispl[i]     *= N_Dim;
      rdisplMS[i]   *= N_Dim;
      rdisplH1[i]   *= N_Dim;
      rdisplH2[i]   *= N_Dim;

      N_DofSend[i]   *= N_Dim;
      N_DofSendMS[i] *= N_Dim;
      N_DofSendH1[i] *= N_Dim;
      N_DofSendH2[i] *= N_Dim;

      N_DofRecv[i]    *= N_Dim;
      N_DofRecvMS[i]  *= N_Dim;
      N_DofRecvH1[i]  *= N_Dim;
      N_DofRecvH2[i]  *= N_Dim;
    }
  }
  
  delete[] Send_data;
  delete[] Recv_data;

  delete[] all_own_dofs_info;
  
  // verification
  /*for (int i = 0; i < N_Dof; i++)
  {
    if(Local2Global[i] == -1)
    {
      ErrThrow("..........all global dof number not assigned..............");
    }
  }*/
  
  if (rank == 0)
  {
    Output::print("Time taken for Assign_GlobalDofNo :: ", MPI_Wtime() - t1);
  }
}


#ifdef _OMP
void TParFEMapper3D::Color(int &numColors, int *&ptrColors, char type)
{
  int i,j,k;
  
  int *myReorder;
  int ndof;   

  // find the type of dof that needs to be colored
  if(type == 'm')
  {
    myReorder = Reorder_M;
    dof = N_InterfaceM;
  }
  else if(type == 'i')
  {
    myReorder = Reorder_I;
    ndof = N_Int;
  }
  else if(type == 'D')
  {
    myReorder = Reorder_D1;
    ndof = N_Dept1;
  }
  else if(type == 'd')
  {
    myReorder = Reorder_D2;
    ndof = N_Dept2;
  }
  else
  {
    ErrThrow("Unsupported type '", type, "' specified!");
  }
  
  //this stores the total number of colors used for coloring
  numColors = 0;
  if (!ndof)
  {
    ptrColors = nullptr;
    return;
  }
  
  int max, temp, t;

  //this stores the color number for dofs
  int *allocatedColor = new int[ndof];

  //initialize all colors with default
  for (i = 0; i < ndof; i++)
  {
    allocatedColor[i] = -1;
  }
  
  //now we start coloring the dofs
  for (i = 0; i < ndof;i ++)
  {
    temp = -1;
    k = myReorder[i];
    for (j = RowPtr[k]; j < RowPtr[k + 1]; j++)
    {
      if (KCol[j] >= k || DofMarker[KCol[j]] != type)
      {
        continue;
      }
      
      t = NewGN[KCol[j]]; // locate the pos of the dof in myreorder to identify its color in allocatedColor

      if (temp < allocatedColor[t])
      {
        temp = allocatedColor[t];
      }
    }
    
    allocatedColor[i] = temp + 1;
    
    if (numColors < allocatedColor[i])
    {
      numColors = allocatedColor[i];
    }
  }

  // colors were numbered from 0, hence total will be 1 more
  numColors++;
  
  ptrColors = new int[numColors + 1];
  temp = 0;
  k = 0;

  // arrange the dofs, such that same color dofs are together
  for (i = 0; i < numColors; i++)
  {
    ptrColors[i] = k;

    for (j = 0; j < ndof; j++)
    {
      if (allocatedColor[j] == i)
      {
        temp = myReorder[j];
        myReorder[j] = myReorder[k];
        myReorder[k] = temp;
        k++;
  
        allocatedColor[j] = -1;
      }
    }
  }
  
  ptrColors[numColors] = k;
  printf("numcolors (type = %c):: %d\t total dof = %d\n",type,numColors,ndof);
}
#endif


TParFEMapper3D::TParFEMapper3D(const TParFEMapper3D& other)
{
  //shallow copies
  FESpace = other.FESpace;

#ifdef _OMP
  RowPtr = other.RowPtr; //shallow copy
  KCol = other.KCol; //shallow copy
#endif

  //copy assign mpi communicator
  Comm = other.Comm;

  //assign non-array built-in type data members
  N_Dim = other.N_Dim;
  N_Dof = other.N_Dof;

  N_InterfaceM = other.N_InterfaceM;
  N_InterfaceS = other.N_InterfaceS;
  N_Halo1 = other.N_Halo1;
  N_Halo2 = other.N_Halo2;
  N_Dept1 = other.N_Dept1;
  N_Dept2 = other.N_Dept2;

  N_Slave = other.N_Slave;
  N_Halo = other.N_Halo;
  N_Master = other.N_Master;
  N_Dept = other.N_Dept;
  N_Int = other.N_Int;
  N_OwnDof = other.N_OwnDof;

  N_SendDof = other.N_SendDof;
  N_SendDofMS = other.N_SendDofMS;
  N_SendDofH1 = other.N_SendDofH1;
  N_SendDofH2 = other.N_SendDofH2;

  N_CMaster = other.N_CMaster;
  N_CDept1 = other.N_CDept1;
  N_CDept2 = other.N_CDept2;
  N_CInt = other.N_CInt;

  //determine rank an size
  int mpiRank, mpiSize;
  MPI_Comm_rank(Comm, &mpiRank);
  MPI_Comm_size(Comm, &mpiSize);

  //switch over the two MapperTypes
  if(TDatabase::ParamDB->MapperType != 2)
  { //master-slave-halo mapping
    Master = new int[N_Dof]; memcpy(Master,other.Master, N_Dof*sizeof(int));
    DofMarker = new char[N_Dof];
    //this is copied in a loop, for we don't store a "SizeOfChar"
    for(int i = 0; i < N_Dof; ++i)
    {
      DofMarker[i]=other.DofMarker[i];
    }

    sdispl = new int[mpiSize]; memcpy(sdispl,other.sdispl, mpiSize*sizeof(int));
    rdispl = new int[mpiSize]; memcpy(rdispl,other.rdispl, mpiSize*sizeof(int));
    sdisplMS = new int[mpiSize]; memcpy(sdisplMS,other.sdisplMS, mpiSize*sizeof(int));
    rdisplMS = new int[mpiSize]; memcpy(rdisplMS,other.rdisplMS, mpiSize*sizeof(int));
    sdisplH1 = new int[mpiSize]; memcpy(sdisplH1,other.sdisplH1, mpiSize*sizeof(int));
    rdisplH1 = new int[mpiSize]; memcpy(rdisplH1,other.rdisplH1, mpiSize*sizeof(int));
    sdisplH2 = new int[mpiSize]; memcpy(sdisplH2,other.sdisplH2, mpiSize*sizeof(int));
    rdisplH2 = new int[mpiSize]; memcpy(rdisplH2,other.rdisplH2, mpiSize*sizeof(int));

    N_DofSend = new int[mpiSize]; memcpy(N_DofSend,other.N_DofSend, mpiSize*sizeof(int));
    N_DofSendMS = new int[mpiSize]; memcpy(N_DofSendMS,other.N_DofSendMS, mpiSize*sizeof(int));
    N_DofSendH2 = new int[mpiSize]; memcpy(N_DofSendH2,other.N_DofSendH2, mpiSize*sizeof(int));
    N_DofSendH1 = new int[mpiSize]; memcpy(N_DofSendH1,other.N_DofSendH1, mpiSize*sizeof(int));

    N_DofRecv = new int[mpiSize]; memcpy(N_DofRecv,other.N_DofRecv, mpiSize*sizeof(int));
    N_DofRecvMS = new int[mpiSize]; memcpy(N_DofRecvMS,other.N_DofRecvMS, mpiSize*sizeof(int));
    N_DofRecvH2 = new int[mpiSize]; memcpy(N_DofRecvH2,other.N_DofRecvH2, mpiSize*sizeof(int));
    N_DofRecvH1 = new int[mpiSize]; memcpy(N_DofRecvH1,other.N_DofRecvH1, mpiSize*sizeof(int));

    OwnDofs = new int[N_OwnDof]; memcpy(OwnDofs,other.OwnDofs, N_OwnDof*sizeof(int));

    DofSend = new int[N_SendDof]; memcpy(DofSend,other.DofSend, N_SendDof*sizeof(int));
    DofSendMS  = DofSend;
    DofSendH1  = DofSend + N_SendDofMS;
    DofSendH2  = DofSend + N_SendDofMS + N_SendDofH1; //pointer arithmetics!


    DofRecv    = new int[N_InterfaceS+N_Halo1+N_Halo2];
    memcpy(DofRecv,other.DofRecv, (N_InterfaceS+N_Halo1+N_Halo2)*sizeof(int));
    DofRecvMS  = DofRecv;
    DofRecvH1  = DofRecv + N_InterfaceS;
    DofRecvH2  = DofRecv + N_InterfaceS + N_Halo1; //pointer arithmetics!

    Reorder = new int[N_Dof]; memcpy(Reorder, other.Reorder, N_Dof*sizeof(int));
    Reorder_M  = Reorder;
    Reorder_I  = Reorder + N_InterfaceM;
    Reorder_D1 = Reorder + N_InterfaceM + N_Int;
    Reorder_D2 = Reorder + N_InterfaceM + N_Int + N_Dept1; //pointer arithmetics!

    NewGN = new int[N_Dof]; memcpy(NewGN, other.NewGN, N_Dof*sizeof(int));

    //TODO CB Maybe I do not like the following - conditional assignment mixed with deterministic.
    if(N_SendDof>0)
    {
      Send_Info   = new double[N_SendDof*N_Dim];
      memcpy(Send_Info, other.Send_Info, N_SendDof*N_Dim*sizeof(int));
    }
    Send_InfoMS = Send_Info;
    Send_InfoH1 = Send_Info + N_SendDofMS*N_Dim;
    Send_InfoH2 = Send_Info + N_SendDofMS*N_Dim + N_SendDofH1*N_Dim;

    if(N_Slave>0)
    {
      Recv_Info   = new double[N_Slave*N_Dim];
      memcpy(Recv_Info, other.Recv_Info, N_Slave*N_Dim*sizeof(int));
    }
    Recv_InfoMS = Recv_Info;
    Recv_InfoH1 = Recv_Info + N_InterfaceS*N_Dim;
    Recv_InfoH2 = Recv_Info + N_InterfaceS*N_Dim + N_Halo1*N_Dim;

    Local2Global = new int[N_Dof]; memcpy(Local2Global, other.Local2Global, N_Dof*sizeof(int));

  }
  else
  { //master-slave mapping
    Master = new int[N_Dof]; memcpy(Master,other.Master, N_Dof*sizeof(int));

    sdispl = new int[mpiSize]; memcpy(sdispl,other.sdispl, mpiSize*sizeof(int));
    rdispl = new int[mpiSize]; memcpy(rdispl,other.rdispl, mpiSize*sizeof(int));

    N_DofSend = new int[mpiSize]; memcpy(N_DofSend,other.N_DofSend, mpiSize*sizeof(int));
    N_DofRecv = new int[mpiSize]; memcpy(N_DofRecv,other.N_DofRecv, mpiSize*sizeof(int));

    OwnDofs = new int[N_OwnDof]; memcpy(OwnDofs,other.OwnDofs, N_OwnDof*sizeof(int));

    DofSend = new int[N_SendDof]; memcpy(DofSend,other.DofSend, N_SendDof*sizeof(int));
    DofRecv    = new int[N_Slave]; memcpy(DofRecv,other.DofRecv, N_Slave*sizeof(int));

      if(N_SendDof>0)
      {
        Send_Info = new double[N_SendDof*N_Dim];
        memcpy(Send_Info,other.Send_Info, (N_SendDof*N_Dim)*sizeof(int));
      }
      if(N_Slave>0)
      {
         Recv_Info = new double[N_Slave*N_Dim];
         memcpy(Recv_Info,other.Recv_Info, (N_Slave*N_Dim)*sizeof(int));
      }
  }

}

//Befriended swap function.
void swap(TParFEMapper3D& first, TParFEMapper3D& second)
{
  //swap values or pointers
  std::swap(first.Comm, second.Comm);
  std::swap(first.N_Dim , second.N_Dim );
  std::swap(first.N_Dof , second.N_Dof );
#ifdef _OMP
  std::swap(first.RowPtr , second.RowPtr );
  std::swap(first.KCol , second.KCol );
#endif
  std::swap(first.FESpace , second.FESpace );
  std::swap(first.Master , second.Master );
  std::swap(first.OwnDofs , second.OwnDofs );
  std::swap(first.DofMarker , second.DofMarker );

  std::swap(first.Send_Info , second.Send_Info );
  std::swap(first.Recv_Info , second.Recv_Info );

  std::swap(first.Local2Global , second.Local2Global );

  std::swap(first.sdispl , second.sdispl );
  std::swap(first.rdispl , second.rdispl );

  std::swap(first.N_InterfaceM , second.N_InterfaceM );
  std::swap(first.N_InterfaceS , second.N_InterfaceS );
  std::swap(first.N_Halo1 , second.N_Halo1 );
  std::swap(first.N_Halo2 , second.N_Halo2 );
  std::swap(first.N_Dept1 , second.N_Dept1 );
  std::swap(first.N_Dept2 , second.N_Dept2 );

  std::swap(first.N_Slave , second.N_Slave );
  std::swap(first.N_Halo , second.N_Halo );
  std::swap(first.N_Master , second.N_Master );
  std::swap(first.N_Dept , second.N_Dept );
  std::swap(first.N_Int , second.N_Int );
  std::swap(first.N_OwnDof , second.N_OwnDof );

  std::swap(first.N_DofSend , second.N_DofSend );
  std::swap(first.N_DofRecv , second.N_DofRecv );
  std::swap(first.N_DofSendMS , second.N_DofSendMS );
  std::swap(first.N_DofSendH1 , second.N_DofSendH1 );
  std::swap(first.N_DofSendH2 , second.N_DofSendH2 );
  std::swap(first.N_DofRecvMS , second.N_DofRecvMS );
  std::swap(first.N_DofRecvH1 , second.N_DofRecvH1 );
  std::swap(first.N_DofRecvH2 , second.N_DofRecvH2 );

  std::swap(first.sdisplMS , second.sdisplMS );
  std::swap(first.sdisplH1 , second.sdisplH1 );
  std::swap(first.sdisplH2 , second.sdisplH2 );
  std::swap(first.rdisplMS , second.rdisplMS );
  std::swap(first.rdisplH1 , second.rdisplH1 );
  std::swap(first.rdisplH2 , second.rdisplH2 );

  std::swap(first.N_SendDof , second.N_SendDof );
  std::swap(first.N_SendDofMS , second.N_SendDofMS );
  std::swap(first.N_SendDofH1 , second.N_SendDofH1 );
  std::swap(first.N_SendDofH2 , second.N_SendDofH2 );

  std::swap(first.DofSend , second.DofSend );
  std::swap(first.DofSendMS , second.DofSendMS );
  std::swap(first.DofSendH1 , second.DofSendH1 );
  std::swap(first.DofSendH2 , second.DofSendH2 );
  std::swap(first.DofRecv , second.DofRecv );
  std::swap(first.DofRecvMS , second.DofRecvMS );
  std::swap(first.DofRecvH1 , second.DofRecvH1 );
  std::swap(first.DofRecvH2 , second.DofRecvH2 );

  std::swap(first.Send_InfoMS , second.Send_InfoMS );
  std::swap(first.Send_InfoH1 , second.Send_InfoH1 );
  std::swap(first.Send_InfoH2 , second.Send_InfoH2 );
  std::swap(first.Recv_InfoMS , second.Recv_InfoMS);
  std::swap(first.Recv_InfoH1 , second.Recv_InfoH1 );
  std::swap(first.Recv_InfoH2 , second.Recv_InfoH2 );

  std::swap(first.NewGN , second.NewGN );
  std::swap(first.Reorder , second.Reorder );
  std::swap(first.Reorder_M , second.Reorder_M );
  std::swap(first.Reorder_I , second.Reorder_I );
  std::swap(first.Reorder_D1 , second.Reorder_D1 );
  std::swap(first.Reorder_D2 , second.Reorder_D2 );

  std::swap(first.N_CMaster , second.N_CMaster );
  std::swap(first.N_CDept1 , second.N_CDept1 );
  std::swap(first.N_CDept2 , second.N_CDept2 );
  std::swap(first.N_CInt , second.N_CInt );
  std::swap(first.ptrCMaster , second.ptrCMaster );
  std::swap(first.ptrCDept1 , second.ptrCDept1 );
  std::swap(first.ptrCDept2 , second.ptrCDept2 );
  std::swap(first.ptrCInt , second.ptrCInt );

}

TParFEMapper3D& TParFEMapper3D::operator=(TParFEMapper3D other)
{
  //do a swap with the copy constructed object "other"
  swap(*this, other);

  return *this;
}

TParFEMapper3D::~TParFEMapper3D()
{
  //does NOT delete the objects pointed to by FESpace
  //     - this belong to other objects
  if(TDatabase::ParamDB->MapperType != 2)
  {
    //dofs were mapped with master/slave/halo concept
    //clean away everything which was allocated in function ConstructDofMap_Master_Halo()
    delete[] Master;
    delete[] DofMarker;
    delete[] sdispl;
    delete[] rdispl;

    delete[] N_DofSend;
    delete[] N_DofSendMS;
    delete[] N_DofSendH2;
    delete[] N_DofSendH1;

    delete[] N_DofRecv;
    delete[] N_DofRecvMS;
    delete[] N_DofRecvH2;
    delete[] N_DofRecvH1;

    delete[] sdisplMS;
    delete[] rdisplMS;
    delete[] sdisplH1;
    delete[] rdisplH1;
    delete[] sdisplH2;
    delete[] rdisplH2;

    delete[] OwnDofs;

    delete[] DofSend;
    delete[] DofRecv;

    delete[] Reorder;
    delete[] NewGN;

    if(N_SendDof>0)
    {
      delete[] Send_Info;
    }
    if(N_Slave>0)
    {
      delete[] Recv_Info;
    }

    // clean away stuff allocated by Assign_GlobalDofNo() (just one array)
    delete[] Local2Global;

  }
  else
  {
    //dofs were mapped with old master/slave concept
    //clean away everything which was allocated in ConstructDofMap()
    delete[] Master;
    delete[] sdispl;
    delete[] rdispl;
    delete[] N_DofSend;
    delete[] N_DofRecv;
    delete[] OwnDofs;

    delete[] DofSend;
    delete[] DofRecv;

    if(N_SendDof>0)
    {
      delete[] Send_Info;
    }
    if(N_Slave>0)
    {
      delete[] Recv_Info;
    }
  }

}

#endif
