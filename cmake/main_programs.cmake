# defining targets for all mein programs

macro(register_main_program name file_name space_dim)
  add_executable(${name} "${PROJECT_SOURCE_DIR}/${file_name}")
  # rename the compiled executable
  set_target_properties(${name} PROPERTIES OUTPUT_NAME 
                        parMooN_${name}_${PARMOON_PARALLEL_TYPE})
  # Link in the required libraries.
  target_link_libraries(${name} PUBLIC 
                        parmoon_${space_dim}d_${PARMOON_PARALLEL_TYPE})
  # create list of all main programs
  list(APPEND parmoon_main_programs ${name})
endmacro(register_main_program)

macro(register_2d_main_program name file_name)
  register_main_program(${name} ${file_name} 2)
endmacro(register_2d_main_program)

macro(register_3d_main_program name file_name)
  register_main_program(${name} ${file_name} 3)
endmacro(register_3d_main_program)

###############################################################################
### Standard 2D Programs (no MPI support in 2D) ###
if(NOT PARMOON_USING_MPI)
  register_2d_main_program(cd2d 2DPrograms/CD2D_ParMooN.C)
  register_2d_main_program(cd2d_afc 2DPrograms/CD2D_AFC_ParMooN.C)
  register_2d_main_program(cd2d_afc_adaptive 2DPrograms/CD2D_AFC_Adaptive_ParMooN.C)
  register_2d_main_program(cd2d_adaptive 2DPrograms/CD2D_Adaptive_ParMooN.C)
  register_2d_main_program(darcy2d 2DPrograms/Darcy2D_ParMooN.C)
  register_2d_main_program(nse2d 2DPrograms/NSE2D_ParMooN.C)
  register_2d_main_program(nse2d_adaptive 2DPrograms/NSE2D_Adaptive_ParMooN.C)
  register_2d_main_program(tcd2d 2DPrograms/TCD2D_ParMooN.C)
  register_2d_main_program(tnse2d 2DPrograms/TNSE2D_ParMooN.C)
  register_2d_main_program(mesh 2DPrograms/mesh_ParMooN.C)
  register_2d_main_program(geo2mesh2d 2DPrograms/geo2mesh2d_ParMooN.C)
#  register_2d_main_program(tcd2dpod 2DPrograms/TCD2D_POD_ParMooN.C)
#  register_2d_main_program(tcd2drom 2DPrograms/TCD2D_ROM_ParMooN.C)
endif(NOT PARMOON_USING_MPI)

###############################################################################
### Standard 3D Programs ###
register_3d_main_program(cd3d 3DPrograms/CD3D_ParMooN.C)
register_3d_main_program(cd3d_afc 3DPrograms/CD3D_AFC_ParMooN.C)
register_3d_main_program(darcy3d 3DPrograms/Darcy3D_ParMooN.C)
register_3d_main_program(tcd3d 3DPrograms/TCD3D_ParMooN.C)
register_3d_main_program(nse3d 3DPrograms/NSE3D_ParMooN.C)
register_3d_main_program(tnse3d 3DPrograms/TNSE3D_ParMooN.C)
register_3d_main_program(mesh3d 3DPrograms/mesh3d_ParMooN.C)

###############################################################################
### Special 3D Programs for applications ###
register_3d_main_program(tnse3d_channeltau 3DPrograms/TNSE3D_ChannelTau.C)
