# CMakeLists.txt for subdirectory PrePostprocessing of ParMooN project.
# Use only as subproject of ParMooN.

# Include header files.
list(APPEND PARMOON_INCLUDE_DIRS "include/PrePostProcessing")

list(APPEND PREPOSTPROCESSING_SOURCES_3D "${PROJECT_SOURCE_DIR}/src/PrePostProcessing/ChannelTauRoutines.C")
list(APPEND PREPOSTPROCESSING_SOURCES "${PROJECT_SOURCE_DIR}/src/PrePostProcessing/SlopeLimiter.C")
list(APPEND PARMOON_SOURCES_2D ${PREPOSTPROCESSING_SOURCES})
list(APPEND PARMOON_SOURCES_3D ${PREPOSTPROCESSING_SOURCES} ${PREPOSTPROCESSING_SOURCES_3D})
