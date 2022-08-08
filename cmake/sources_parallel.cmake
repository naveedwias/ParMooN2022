# CMakeLists.txt for subdirectory Parallel of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2015/08/20 Clemens Bartsch: Rework to supply 2D and 3D library at once.
#

if(PARMOON_USING_MPI)
  # Include header files. 
  list(APPEND PARMOON_INCLUDE_DIRS "include/Parallel")
  
  list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/MeshPartition.C")
  list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/MeshPartitionInOut.C")  
  list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/MumpsWrapper.C")
  list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/ParFEMapper3D.C")
  list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/ParFECommunicator3D.C")
  
  list(APPEND PARMOON_SOURCES_2D ${PAR_SOURCES})
  list(APPEND PARMOON_SOURCES_3D ${PAR_SOURCES})
endif(PARMOON_USING_MPI)
