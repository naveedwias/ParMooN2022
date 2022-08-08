# CMakeLists.txt for subdirectory Matrix of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Change history:
# 2015/08/25 Ulrich Wilbrandt: Adjust to Clemens changes to the cmake system
#

# Include header files.
list(APPEND PARMOON_INCLUDE_DIRS "include/Matrix")

# Source files to be added to the 2D and 3D library.
list(APPEND MATRIX_SOURCES "${PROJECT_SOURCE_DIR}/src/Matrix/BlockMatrix.C")
list(APPEND MATRIX_SOURCES "${PROJECT_SOURCE_DIR}/src/Matrix/BlockFEMatrix.C")
list(APPEND MATRIX_SOURCES "${PROJECT_SOURCE_DIR}/src/Matrix/BlockVector.C")
list(APPEND MATRIX_SOURCES "${PROJECT_SOURCE_DIR}/src/Matrix/DenseMatrix.C")
list(APPEND MATRIX_SOURCES "${PROJECT_SOURCE_DIR}/src/Matrix/FEMatrix.C")
list(APPEND MATRIX_SOURCES "${PROJECT_SOURCE_DIR}/src/Matrix/Matrix.C")
list(APPEND MATRIX_SOURCES "${PROJECT_SOURCE_DIR}/src/Matrix/SquareMatrix.C")
list(APPEND MATRIX_SOURCES "${PROJECT_SOURCE_DIR}/src/Matrix/Structure.C")
list(APPEND MATRIX_SOURCES "${PROJECT_SOURCE_DIR}/src/Matrix/CompositeOperator.C")

#2D - only temporary, until there is no Matrix2D/Matrix3D anymore
list(APPEND MATRIX_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Matrix/Matrix2D.C")
list(APPEND MATRIX_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Matrix/SquareMatrix1D.C")
list(APPEND MATRIX_SOURCES_2D "${PROJECT_SOURCE_DIR}/src/Matrix/SquareMatrix2D.C")

#3D - only temporary, until there is no Matrix2D/Matrix3D anymore
list(APPEND MATRIX_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Matrix/Matrix3D.C")
list(APPEND MATRIX_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/Matrix/SquareMatrix3D.C")

list(APPEND PARMOON_SOURCES_2D ${MATRIX_SOURCES} ${MATRIX_SOURCES_2D})
list(APPEND PARMOON_SOURCES_3D ${MATRIX_SOURCES} ${MATRIX_SOURCES_3D})
