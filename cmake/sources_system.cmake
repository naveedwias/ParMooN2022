# CMakeLists.txt for subdirectory System of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2015/08/20 Clemens Bartsch: Rework to supply 2D and 3D library at once.
#

# Include header files. 
list(APPEND PARMOON_INCLUDE_DIRS "include/System")

# Source files to be added to the 2D library.
list(APPEND SYST_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/System/ConvectionDiffusion.C")
list(APPEND SYST_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/System/ConvectionDiffusion_AFC.C")
list(APPEND SYST_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/System/NavierStokes.C")
list(APPEND SYST_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/System/NavierStokes_DG.C")
list(APPEND SYST_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/System/Darcy.C")
list(APPEND SYST_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/System/TimeConvectionDiffusion.C")
list(APPEND SYST_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/System/TimeNavierStokes.C")

# Source files to be added to the 3D library.
list(APPEND SYST_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/System/ConvectionDiffusion.C")
list(APPEND SYST_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/System/ConvectionDiffusion_AFC.C")
list(APPEND SYST_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/System/NavierStokes.C")
#list(APPEND SYST_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/System/NavierStokes_DG.C")
list(APPEND SYST_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/System/Darcy.C")
list(APPEND SYST_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/System/TimeConvectionDiffusion.C")
list(APPEND SYST_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/System/TimeNavierStokes.C")

# Source files to be added to both 2D and 3D library
list(APPEND SYST_SOURCES "${PROJECT_SOURCE_DIR}/src/System/Residuals.C")
list(APPEND SYST_SOURCES "${PROJECT_SOURCE_DIR}/src/System/WindkesselBoundary.C")

list(APPEND PARMOON_SOURCES_2D ${SYST_SOURCES} ${SYST_SOURCES_2D})
list(APPEND PARMOON_SOURCES_3D ${SYST_SOURCES} ${SYST_SOURCES_3D})
