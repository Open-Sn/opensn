# FindPETSc.cmake
#
# Once done this will define
#  PETSC_FOUND - System has PETSc
#  PETSC_INCLUDE_DIR - The PETSc include directory
#  PETSC_LIBRARY - The PETSc library
#  PETSC_VERSION - The PETSc version
#  SLEPC_FOUND - PETSc has SLEPc
#  SLEPC_INCLUDE_DIR - The SLEPc include directory
#  SLEPC_LIBRARY - The SLEPc library
#  SLEPC_VERSION

include(FindPackageHandleStandardArgs)

find_package(PETSc CONFIG QUIET
  PATHS
    $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/cmake/PETSc
    $ENV{PETSC_ROOT}/$ENV{PETSC_ARCH}/lib/cmake/PETSc
)

if (NOT PETSc_FOUND)
  find_path(PETSC_INCLUDE_DIR
    NAMES petsc.h
    PATHS
      $ENV{PETSC_DIR}/include
      $ENV{PETSC_ROOT}/include
    DOC "PETSc include directory"
  )

  find_library(PETSC_LIBRARY
    NAMES petsc
    PATHS
      $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
      $ENV{PETSC_DIR}/lib
      $ENV{PETSC_ROOT}/$ENV{PETSC_ARCH}/lib
      $ENV{PETSC_ROOT}/lib
    DOC "PETSc library"
  )

  set(PETSc_FOUND
    ${PETSC_INCLUDE_DIR}
    AND
    ${PETSC_LIBRARY}
  )
endif()

if (TARGET PETSc::PETSc)
  get_target_property(_petsc_version PETSc::PETSc INTERFACE_VERSION)
  set(PETSC_VERSION ${_petsc_version})
endif()

if(NOT DEFINED PETSC_VERSION)
  set(PETSC_VERSION "unknown")
  find_file(PETSCVERSION_H
    NAMES petscversion.h
    PATHS
      $ENV{PETSC_DIR}/include
      $ENV{PETSC_ROOT}/include
  )
  if(PETSCVERSION_H)
    file(READ ${PETSCVERSION_H} _petscver)
    string(REGEX MATCH "define[ ]+PETSC_VERSION_MAJOR[ ]+([0-9]+)" _ "${_petscver}")
    set(PETSC_VERSION_MAJOR ${CMAKE_MATCH_1})
    string(REGEX MATCH "define[ ]+PETSC_VERSION_MINOR[ ]+([0-9]+)" _ "${_petscver}")
    set(PETSC_VERSION_MINOR ${CMAKE_MATCH_1})
    string(REGEX MATCH "define[ ]+PETSC_VERSION_SUBMINOR[ ]+([0-9]+)" _ "${_petscver}")
    set(PETSC_VERSION_PATCH ${CMAKE_MATCH_1})
    set(PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_PATCH}")
  endif()
endif()

find_path(SLEPC_INCLUDE_DIR
  NAMES slepc.h
  PATHS
    $ENV{PETSC_DIR}/include/slepc
    $ENV{PETSC_ROOT}/include/slepc
    $ENV{SLEPC_DIR}/include
    $ENV{SLEPC_ROOT}/include
  DOC "SLEPc include directory"
)

find_library(SLEPC_LIBRARY
  NAMES slepc
  PATHS
    $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
    $ENV{PETSC_DIR}/lib
    $ENV{PETSC_ROOT}/$ENV{PETSC_ARCH}/lib
    $ENV{PETSC_ROOT}/lib
    $ENV{SLEPC_DIR}/$ENV{SLEPC_ARCH}/lib
    $ENV{SLEPC_ROOT}/$ENV{SLEPC_ARCH}/lib
    $ENV{SLEPC_DIR}/lib
    $ENV{SLEPC_ROOT}/lib
  DOC "SLEPc library"
)

set(SLEPC_VERSION "unknown")
find_file(SLEPCVERSION_H
  NAMES slepcversion.h
  PATHS
    $ENV{PETSC_DIR}/include/slepc
    $ENV{PETSC_ROOT}/include/slepc
    $ENV{SLEPC_DIR}/include
    $ENV{SLEPC_ROOT}/include
)
if(SLEPCVERSION_H)
  file(READ ${SLEPCVERSION_H} _slepcver)
  string(REGEX MATCH "define[ ]+SLEPC_VERSION_MAJOR[ ]+([0-9]+)" _ "${_slepcver}")
  set(SLEPC_VERSION_MAJOR ${CMAKE_MATCH_1})
  string(REGEX MATCH "define[ ]+SLEPC_VERSION_MINOR[ ]+([0-9]+)" _ "${_slepcver}")
  set(SLEPC_VERSION_MINOR ${CMAKE_MATCH_1})
  string(REGEX MATCH "define[ ]+SLEPC_VERSION_SUBMINOR[ ]+([0-9]+)" _ "${_slepcver}")
  set(SLEPC_VERSION_PATCH ${CMAKE_MATCH_1})
  set(SLEPC_VERSION "${SLEPC_VERSION_MAJOR}.${SLEPC_VERSION_MINOR}.${SLEPC_VERSION_PATCH}")
endif()

if (NOT SLEPC_INCLUDE_DIR)
  message(FATAL_ERROR "PETSc has not been configured with SLEPc support.")
endif()
if (NOT SLEPC_LIBRARY)
  message(FATAL_ERROR "PETSc has not been configured with SLEPc support.")
endif()
set(SLEPc_FOUND TRUE)

list(APPEND PETSC_LIBRARY ${SLEPC_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    PETSc
    REQUIRED_VARS PETSC_LIBRARY PETSC_INCLUDE_DIR
    VERSION_VAR PETSC_VERSION
)
