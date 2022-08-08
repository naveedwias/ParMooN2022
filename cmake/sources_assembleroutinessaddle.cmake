# CMakeLists.txt for subdirectory AssembleRoutinesSaddle of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2017/03/22 Najib Alia: Rework to separate FE and AssembleRoutines files.
# 2017/03/22 Naveed Ahmed: splitting of the local assembling routnies for the saddle point problems
#

# Include header files. 
list(APPEND PARMOON_INCLUDE_DIRS "include/AssembleRoutinesSaddle")

# Source files used in 2D and 3D.
list(APPEND ASSEMBLESADDLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/NSE2DSUPG.C")
list(APPEND ASSEMBLESADDLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/CommonRoutineTNSE3D.C")
list(APPEND ASSEMBLESADDLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/NSE_local_assembling_routines.C")
list(APPEND ASSEMBLESADDLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/LPS_scott_zhang.C")
list(APPEND ASSEMBLESADDLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/Time_NSE_local_assembling_routines.C")
list(APPEND ASSEMBLESADDLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/Variational_multiscale.C")
list(APPEND ASSEMBLESADDLE_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/NonNewtonianViscosity.C")

list(APPEND ASSEMBLESADDLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/CommonRoutineTNSE3D.C")
list(APPEND ASSEMBLESADDLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/NSE_local_assembling_routines.C")
list(APPEND ASSEMBLESADDLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/Time_NSE_local_assembling_routines.C")
list(APPEND ASSEMBLESADDLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/Variational_multiscale.C")
list(APPEND ASSEMBLESADDLE_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/AssembleRoutinesSaddle/NonNewtonianViscosity.C")

list(APPEND PARMOON_SOURCES_2D ${ASSEMBLESADDLE_SOURCES_2D})
list(APPEND PARMOON_SOURCES_3D ${ASSEMBLESADDLE_SOURCES_3D})
