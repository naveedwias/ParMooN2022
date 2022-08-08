# CMakeLists.txt for subdirectory Solver of ParMooN project. 
# Use only as subproject of ParMooN.
# 
# Include header files. 
list(APPEND PARMOON_INCLUDE_DIRS "include/Solver")

# Source files to be added to the 2D and 3D library.
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/DirectSolver.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Utilities_gmres.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Iteration_bicgstab.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Iteration_mixed.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Iteration_cg.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Iteration_cgs.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Iteration_gmres.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Iteration_jacobi.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Iteration_multigrid.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Iteration_richardson.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Iteration_sor.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/PETScSolver.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Preconditioner_vanka.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Preconditioner_iterative.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Preconditioner_optimized.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Saddle_point_preconditioner.C")
list(APPEND SOLVER_SOURCES "${PROJECT_SOURCE_DIR}/src/Solver/Solver.C")

list(APPEND PARMOON_SOURCES_2D ${SOLVER_SOURCES})
list(APPEND PARMOON_SOURCES_3D ${SOLVER_SOURCES})
