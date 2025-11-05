set(OPENSN_ENV_DIR "${CMAKE_INSTALL_PREFIX}/bin")
set(ENV_SCRIPT "${OPENSN_ENV_DIR}/set_opensn_env.sh")
file(MAKE_DIRECTORY "${OPENSN_ENV_DIR}")
file(WRITE "${ENV_SCRIPT}" "export CMAKE_PREFIX_PATH=\"${CMAKE_INSTALL_PREFIX}\"\${CMAKE_PREFIX_PATH:+:\${CMAKE_PREFIX_PATH}}\n")
file(CHMOD "${ENV_SCRIPT}" FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)

message(STATUS "
*******************************************************************************
*******************************************************************************

 OpenSn dependency install complete.

 To update your CMAKE_PREFIX_PATH for the current shell, run:

   source ${ENV_SCRIPT}

 (This adjusts only the current shell session. Add the line above to your
  shell init file - e.g., ~/.bashrc or ~/.zshrc - to make it permanent.)

*******************************************************************************
*******************************************************************************")
