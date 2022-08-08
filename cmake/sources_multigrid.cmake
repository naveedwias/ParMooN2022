# CMakeLists.txt for subdirectory Multigrid of ParMooN project. 
# Use only as subproject of ParMooN.


# Include header files. 
list(APPEND PARMOON_INCLUDE_DIRS "include/Multigrid")

list(APPEND MULTIGRID_SOURCES "${PROJECT_SOURCE_DIR}/src/Multigrid/CycleControl.C")
list(APPEND MULTIGRID_SOURCES "${PROJECT_SOURCE_DIR}/src/Multigrid/DirectSmoother.C")
list(APPEND MULTIGRID_SOURCES "${PROJECT_SOURCE_DIR}/src/Multigrid/DofPatch.C") 
list(APPEND MULTIGRID_SOURCES "${PROJECT_SOURCE_DIR}/src/Multigrid/GridTransfer.C")
list(APPEND MULTIGRID_SOURCES "${PROJECT_SOURCE_DIR}/src/Multigrid/JacobiSmoother.C")
list(APPEND MULTIGRID_SOURCES "${PROJECT_SOURCE_DIR}/src/Multigrid/Multigrid.C")
list(APPEND MULTIGRID_SOURCES "${PROJECT_SOURCE_DIR}/src/Multigrid/MultigridLevel.C")
list(APPEND MULTIGRID_SOURCES "${PROJECT_SOURCE_DIR}/src/Multigrid/SORSmoother.C")
list(APPEND MULTIGRID_SOURCES "${PROJECT_SOURCE_DIR}/src/Multigrid/VankaSmoother.C")

list(APPEND PARMOON_SOURCES_2D ${MULTIGRID_SOURCES})
list(APPEND PARMOON_SOURCES_3D ${MULTIGRID_SOURCES})
