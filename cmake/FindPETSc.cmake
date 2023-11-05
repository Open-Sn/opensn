# Find PETSc
#
# Once done this will define
#  PETSC_FOUND - System has PETSc
#  PETSC_INCLUDE_DIR - The PETSc include directory
#  PETSC_LIBRARY - The PETSc library
#  PETSC_VERSION - The PETSc version

find_path(
    PETSC_INCLUDE_DIR
        petsc.h
    PATHS
        $ENV{PETSC_DIR}/include
        $ENV{PETSC_ROOT}/include
)

find_library(
    PETSC_LIBRARY
        petsc
    PATHS
        $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
        $ENV{PETSC_DIR}/lib
        $ENV{PETSC_ROOT}/$ENV{PETSC_ARCH}/lib
        $ENV{PETSC_ROOT}/lib
)

set(PETSC_VERSION "unknown")
find_file(PETSCVERSION_H petscversion.h
    PATHS
        $ENV{PETSC_DIR}/include
        $ENV{PETSC_ROOT}/include
)
mark_as_advanced(PETSCVERSION_H)
if (PETSCVERSION_H)
    file(READ ${PETSCVERSION_H} PETSC_VERSION_FILE)
    string(REGEX MATCH "define[ ]+PETSC_VERSION_MAJOR[ ]+([0-9]+)" TMP "${PETSC_VERSION_FILE}")
    set(PETSC_VERSION_MAJOR ${CMAKE_MATCH_1})
    string(REGEX MATCH "define[ ]+PETSC_VERSION_MINOR[ ]+([0-9]+)" TMP "${PETSC_VERSION_FILE}")
    set(PETSC_VERSION_MINOR ${CMAKE_MATCH_1})
    string(REGEX MATCH "define[ ]+PETSC_VERSION_SUBMINOR[ ]+([0-9]+)" TMP "${PETSC_VERSION_FILE}")
    set(PETSC_VERSION_PATCH ${CMAKE_MATCH_1})
    set(PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_PATCH}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    PETSc
    REQUIRED_VARS PETSC_LIBRARY PETSC_INCLUDE_DIR
    VERSION_VAR PETSC_VERSION
)
