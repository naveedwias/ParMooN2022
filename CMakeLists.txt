# ========================================================================
#     This file is part of ParMooN.
# 
#     ParMooN is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     ParMooN is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with ParMooN.  If not, see <http://www.gnu.org/licenses/>.
# =======================================================================
# PURPOSE :  configuration file for the ParMooN Version 1.1
# AUTHOR  : Sashikumaar Ganesan, SERC, IISc Bangalore, India
#			Clemens Bartsch, WIAS Berlin
# DATE    : 05 June 2015
# ========================================================================

# normally no need to alter any thing in this file
# please DO NOT modify unless it is necessary
# beginners are advised NOT to alter anything of the below contents
# ========================================================================

cmake_minimum_required(VERSION 3.5)

project(ParMooN VERSION 2.0)
enable_language(Fortran)
include(FortranCInterface) 
#get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} gfortran)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Debug' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

# determine the parallel type (see https://stackoverflow.com/a/39671055)
list(APPEND PARMOON_PARALLEL_TYPE_LIST "MPI" "OMPONLY" "HYBRID" "SEQUENTIAL")
set(PARMOON_PARALLEL_TYPE "SEQUENTIAL" 
    CACHE STRING "The used parallelism. Choose between ${PARMOON_PARALLEL_TYPE_LIST}.")
# set the allowed strings as property. This way (graphical) editors for the 
# cache (CMakeCache.txt) can select only the allowed values
set_property(CACHE PARMOON_PARALLEL_TYPE PROPERTY STRINGS ${PARMOON_PARALLEL_TYPE_LIST})
if(NOT PARMOON_PARALLEL_TYPE IN_LIST PARMOON_PARALLEL_TYPE_LIST)
  message(FATAL_ERROR "Unknown PARMOON_PARALLEL_TYPE ${PARMOON_PARALLEL_TYPE}. "
                      "Choose between ${PARMOON_PARALLEL_TYPE_LIST}")
endif()
message(STATUS "Chosen parallel type is " ${PARMOON_PARALLEL_TYPE}.)

# Set two auxiliary variables, which make the file better readable imho
if("${PARMOON_PARALLEL_TYPE}" STREQUAL "MPI" OR "${PARMOON_PARALLEL_TYPE}" STREQUAL "HYBRID")
  set(PARMOON_USING_MPI TRUE)
endif()
if("${PARMOON_PARALLEL_TYPE}" STREQUAL "OMPONLY" OR "${PARMOON_PARALLEL_TYPE}" STREQUAL "HYBRID")
  set(PARMOON_USING_OMP TRUE)
endif()

# include all external libraries, this creates a list of libraries and a list of
# include directories: PARMOON_EXTERNAL_LIBRARIES, PARMOON_EXTERNAL_INCLUDES
include(cmake/external_libraries.cmake)

###############################################################################
# parmoon defines 
set(PARMOON_CONFIG_DIR ${CMAKE_CURRENT_BINARY_DIR}/config)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/all_defines_external_libraries.h.in
               ${PARMOON_CONFIG_DIR}/all_defines_external_libraries.h)
###############################################################################
# make sure the correct information on the repository (current commit, branch,
# ...) and build configuration (build type, parallel, path to directory) are
# available during compile and run time. This information will be updated during
# before building the ParMooN libraries.
# For this we define a target 'repository_info' on which the ParMooN libraries
# depend on. Unfortunately, it is not possible to call a cmake function in
# add_custom_target(...), so that we have to call another cmake script instead.
# This script in turn is created via configure_file(...).
find_package(Git)
configure_file(${CMAKE_SOURCE_DIR}/cmake/get_repository_information.cmake.in
               ${PARMOON_CONFIG_DIR}/get_repository_information.cmake @ONLY)
add_custom_target(repository_info
                  COMMAND ${CMAKE_COMMAND}
                  -D PARMOON_PARALLEL_TYPE=${PARMOON_PARALLEL_TYPE}
                  -D PARMOON_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                  -D PARMOON_SOURCE_DIR=${CMAKE_SOURCE_DIR}
                  -P ${PARMOON_CONFIG_DIR}/get_repository_information.cmake
                  BYPRODUCTS ${PARMOON_CONFIG_DIR}/ParMooN_repository_info.cpp)
#list(APPEND PARMOON_INCLUDE_DIRS $<BUILD_INTERFACE:${PARMOON_CONFIG_DIR}>)
list(APPEND PARMOON_SOURCES_2D ${PARMOON_CONFIG_DIR}/ParMooN_repository_info.cpp)
list(APPEND PARMOON_SOURCES_3D ${PARMOON_CONFIG_DIR}/ParMooN_repository_info.cpp)

###############################################################################
# turn on all warnings
add_compile_options(
  # clang/GCC warnings
  $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wall>
  $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wextra>
  $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-pedantic>
  $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-Wno-attributes>
  # MSVC warnings
  $<$<CXX_COMPILER_ID:MSVC>:"/W4">)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  add_compile_options(-D__MAC64__) # how can we get rid of this?
endif()


###############################################################################
# If we want global compile flags to reach the subdirectories, we have to 
# run the following includes AFTER setting global flags.
include(cmake/sources_assembleroutines.cmake)
include(cmake/sources_assembleroutinessaddle.cmake)
include(cmake/sources_examples.cmake)
include(cmake/sources_fe.cmake)
include(cmake/sources_general.cmake)
include(cmake/sources_geometry.cmake)
include(cmake/sources_matrix.cmake)
include(cmake/sources_multigrid.cmake)
include(cmake/sources_parallel.cmake)
include(cmake/sources_quadformulas.cmake)
include(cmake/sources_refinement.cmake)
include(cmake/sources_rom.cmake)
include(cmake/sources_solver.cmake)
include(cmake/sources_system.cmake)
include(cmake/sources_timediscroutines.cmake)
# include(cmake/sources_pbe.cmake)
include(cmake/sources_prePostprocessing.cmake)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

###############################################################################
# Define two static libraries. Each has its own precompiler flag (2D/3D).
if(NOT PARMOON_USING_MPI)
  add_library(parmoon_2d_${PARMOON_PARALLEL_TYPE} STATIC ${PARMOON_SOURCES_2D})
  add_library(parmoon::parmoon_2d ALIAS parmoon_2d_${PARMOON_PARALLEL_TYPE})
  target_compile_definitions(parmoon_2d_${PARMOON_PARALLEL_TYPE} PUBLIC __2D__)
  target_compile_definitions(parmoon_2d_${PARMOON_PARALLEL_TYPE} PUBLIC 
                             $<$<BOOL:${PARMOON_USING_MPI}>:_MPI>
                             $<$<BOOL:${PARMOON_USING_OMP}>:_OMP>)
  target_include_directories(parmoon_2d_${PARMOON_PARALLEL_TYPE}
      PUBLIC 
          "$<INSTALL_INTERFACE:${PARMOON_INCLUDE_DIRS}>"
          $<BUILD_INTERFACE:${PARMOON_INCLUDE_DIRS}>
          $<BUILD_INTERFACE:${PARMOON_CONFIG_DIR}>
          )
  target_link_libraries(parmoon_2d_${PARMOON_PARALLEL_TYPE} parmoon::external_libraries)
  #target_compile_features(parmoon_2d_${PARMOON_PARALLEL_TYPE} PUBLIC cxx_std_11) # cmake 3.8
  add_dependencies(parmoon_2d_${PARMOON_PARALLEL_TYPE} repository_info)
endif(NOT PARMOON_USING_MPI)

add_library(parmoon_3d_${PARMOON_PARALLEL_TYPE} STATIC ${PARMOON_SOURCES_3D})
add_library(parmoon::parmoon_3d ALIAS parmoon_3d_${PARMOON_PARALLEL_TYPE})
target_compile_definitions(parmoon_3d_${PARMOON_PARALLEL_TYPE} PUBLIC __3D__)
target_compile_definitions(parmoon_3d_${PARMOON_PARALLEL_TYPE} PUBLIC 
                           $<$<BOOL:${PARMOON_USING_MPI}>:_MPI>
                           $<$<BOOL:${PARMOON_USING_OMP}>:_OMP>)
target_include_directories(parmoon_3d_${PARMOON_PARALLEL_TYPE} 
    PUBLIC
        "$<INSTALL_INTERFACE:${PARMOON_INCLUDE_DIRS}>"
        $<BUILD_INTERFACE:${PARMOON_INCLUDE_DIRS}>
        $<BUILD_INTERFACE:${PARMOON_CONFIG_DIR}>
        )
target_link_libraries(parmoon_3d_${PARMOON_PARALLEL_TYPE} parmoon::external_libraries)
#target_compile_features(parmoon_3d_${PARMOON_PARALLEL_TYPE} PUBLIC cxx_std_11) # cmake 3.8
add_dependencies(parmoon_3d_${PARMOON_PARALLEL_TYPE} repository_info)


##############################################################
# Start adding targets.

include(cmake/main_programs.cmake)
include(cmake/documentation_target.cmake)
include(cmake/install.cmake)

###############################################################################
enable_testing()

###############################################################################
# tests 
add_subdirectory(test EXCLUDE_FROM_ALL)
####################################################
#                THE CHECK TARGET                  #
#   making this target results in building and     #
#   running all tests corresponding to the current #
#   parallel type. Use it excessively!             #
####################################################
add_custom_target(build_all_tests DEPENDS ${parmoon_tests})
add_custom_target(build_all_main DEPENDS ${parmoon_main_programs})
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure --exclude-regex all
                  DEPENDS build_all_tests build_all_main generate_standard_input_files)
add_custom_target(check_all COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure 
                  DEPENDS build_all_tests build_all_main generate_standard_input_files)

include(cmake/additional_helper_targets.cmake)

############
# End adding executable targets.
############

message("=====================================================================================================")        
message("Configuration for ParMooN complete, Version " ${ParMooN_VERSION})
message("=====================================================================================================") 
