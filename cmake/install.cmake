# This file includes everything needed to install ParMooN

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  # we want to avoid accidental installs in system directories (for now)
  message(STATUS "You did not provide a 'CMAKE_INSTALL_PREFIX', therefore you "
                 "will not be able to install this built of parmmon. "
                 "Everything else continues to work.")
else(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  # provide some commonly used directory names. These are relative paths
  # prefixed by 'CMAKE_INSTALL_PREFIX'.
  include(GNUInstallDirs)
  
  set(PARMOON_INSTALL_CONFIG_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/parmoon)
  # no mpi in 2D
  if(NOT PARMOON_USING_MPI)
    list(APPEND PARMMON_INSTALL_LIBS parmoon_2d_${PARMOON_PARALLEL_TYPE})
  endif(NOT PARMOON_USING_MPI)
  list(APPEND PARMMON_INSTALL_LIBS parmoon_3d_${PARMOON_PARALLEL_TYPE})
  
  # Install the parmoon libraries. For this to work also in other cmake 
  # projects, we also have to explicitly include the external libraries here.
  # the destinations are standard directory names for installing libraries 
  # (include, lib)
  install(TARGETS ${PARMMON_INSTALL_LIBS} external_libraries
          EXPORT parmoon-targets
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
          INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  # during install, copy all the header files into ${CMAKE_INSTALL_INCLUDEDIR}
  install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/ 
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  # also copy header files which were created by cmake
  install(FILES ${PARMOON_CONFIG_DIR}/all_defines_external_libraries.h
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  # save the parmoon libraries (and external libs) in cmake file.
  install(EXPORT parmoon-targets
          FILE parmoon-targets.cmake
          NAMESPACE parmoon::
          DESTINATION ${PARMOON_INSTALL_CONFIG_DIR})
  
  # to use these libraries via find_package in another cmake project, we still
  # need 'parmoon-config.cmake' and 'parmoon-config-version.cmake'. The 
  # following commands simplify their creation
  include(CMakePackageConfigHelpers)
  
  configure_package_config_file(
      ${CMAKE_SOURCE_DIR}/cmake/parmoon-config.cmake.in
      ${PARMOON_CONFIG_DIR}/parmoon-config.cmake
      INSTALL_DESTINATION ${PARMOON_INSTALL_CONFIG_DIR})
  
  write_basic_package_version_file(
      ${PARMOON_CONFIG_DIR}/parmoon-config-version.cmake
      VERSION ${ParMooN_VERSION}
      COMPATIBILITY AnyNewerVersion)

  # these files also need to be copied into the appropriate subdirectory
  install(
      FILES
          ${PARMOON_CONFIG_DIR}/parmoon-config.cmake
          ${PARMOON_CONFIG_DIR}/parmoon-config-version.cmake
      DESTINATION ${PARMOON_INSTALL_CONFIG_DIR})
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
