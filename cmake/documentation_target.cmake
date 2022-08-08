# building the documentation using doxygen (define target 'doc')
find_package(Doxygen)
if(DOXYGEN_FOUND)
  set(PARMOON_DOCUMENTATION_DIRECTORY "" CACHE PATH 
      "The directory for the Doxygen documentation.")
  if(NOT PARMOON_DOCUMENTATION_DIRECTORY)
    add_custom_target(doc ${CMAKE_COMMAND} -E cmake_echo_color --red 
                      "Please set PARMOON_DOCUMENTATION_DIRECTORY in your "
                      "CMakeCache.txt file. No documentation is created!")
  else()
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/documentation/Doxyfile.in
                   ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
                      ${DOXYGEN_EXECUTABLE} 
                      ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                      COMMENT "Generating documentation with Doxygen" VERBATIM)
  endif(NOT PARMOON_DOCUMENTATION_DIRECTORY)
else()
  add_custom_target(doc ${CMAKE_COMMAND} -E cmake_echo_color --red 
                    "It seems cmake could not find a doxygen executable. "
                    "Therefore, no documentation is created!")
endif(DOXYGEN_FOUND)
