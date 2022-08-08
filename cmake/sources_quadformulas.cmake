# CMakeLists.txt for subdirectory QuadFormulas of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2015/08/20 Clemens Bartsch: Rework to supply 2D and 3D library at once.
#

# Include header files. 
list(APPEND PARMOON_INCLUDE_DIRS "include/QuadFormulas")

# Source files used in 2D and 3D
list(APPEND QUAD_SOURCES "${PROJECT_SOURCE_DIR}/src/QuadFormulas/Enumerations_quadrature_formula.C")
list(APPEND QUAD_SOURCES "${PROJECT_SOURCE_DIR}/src/QuadFormulas/QuadFormula.C")
list(APPEND QUAD_SOURCES "${PROJECT_SOURCE_DIR}/src/QuadFormulas/QuadratureFormulaDatabase.C")

list(APPEND PARMOON_SOURCES_2D ${QUAD_SOURCES})
list(APPEND PARMOON_SOURCES_3D ${QUAD_SOURCES})
