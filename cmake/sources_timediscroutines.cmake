# CMakeLists.txt for subdirectory System of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
#

# Include header files. 
list(APPEND PARMOON_INCLUDE_DIRS "include/TimeDiscRoutines")

list(APPEND TimeDisc_SOURCES "${PROJECT_SOURCE_DIR}/src/TimeDiscRoutines/TimeDiscretizations.C")
# Source files used in 2D and 3D.

list(APPEND PARMOON_SOURCES_2D ${TimeDisc_SOURCES} )
list(APPEND PARMOON_SOURCES_3D ${TimeDisc_SOURCES} )

list(APPEND PARMOON_SOURCES_2D ${TimeDisc_SOURCES_2D} )
# list(APPEND PARMOON_SOURCES_3D ${TD_SOURCES_3D} )
