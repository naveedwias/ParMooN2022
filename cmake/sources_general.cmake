# CMakeLists.txt for subdirectory General of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2015/08/20 Clemens Bartsch: Rework to supply 2D and 3D library at once.
#

# Include header files. 
list(APPEND PARMOON_INCLUDE_DIRS "include/General")

# Source files used in 2D and 3D.
list(APPEND GENERAL_SOURCES "${PROJECT_SOURCE_DIR}/src/General/Chrono.C")
list(APPEND GENERAL_SOURCES "${PROJECT_SOURCE_DIR}/src/General/Database.C")
list(APPEND GENERAL_SOURCES "${PROJECT_SOURCE_DIR}/src/General/FunctionEvaluation.C")
list(APPEND GENERAL_SOURCES "${PROJECT_SOURCE_DIR}/src/General/LoopInfo.C")
list(APPEND GENERAL_SOURCES "${PROJECT_SOURCE_DIR}/src/General/MainUtilities.C")
list(APPEND GENERAL_SOURCES "${PROJECT_SOURCE_DIR}/src/General/MooNMD_Io.C")
list(APPEND GENERAL_SOURCES "${PROJECT_SOURCE_DIR}/src/General/Parameter.C")
list(APPEND GENERAL_SOURCES "${PROJECT_SOURCE_DIR}/src/General/ParMooN.C")
list(APPEND GENERAL_SOURCES "${PROJECT_SOURCE_DIR}/src/General/ParameterDatabase.C")
list(APPEND GENERAL_SOURCES "${PROJECT_SOURCE_DIR}/src/General/Utilities.C")
list(APPEND GENERAL_SOURCES "${PROJECT_SOURCE_DIR}/src/General/anderson.f90")

list(APPEND PARMOON_SOURCES_2D ${GENERAL_SOURCES})
list(APPEND PARMOON_SOURCES_3D ${GENERAL_SOURCES})
