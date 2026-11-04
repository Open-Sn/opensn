# Find PETSc
#
# Once done this will define
#  PETSC_FOUND - System has PETSc
#  PETSC_INCLUDE_DIR - The PETSc include directory
#  PETSC_LIBRARY - The PETSc library
#  PETSC_VERSION - The PETSc version
#  PETSC_BLASLAPACK_LIBRARIES - The BLAS/LAPACK libraries PETSc was built with, from
#                               petscvariables, or else those of a system LAPACK
#  PETSc::PETSc - Imported target for the PETSc library and include directory
#  PETSc::BLASLAPACK - Imported target for PETSC_BLASLAPACK_LIBRARIES, defined when they are found

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

# BLAS/LAPACK that PETSc was built with. PETSc keeps it private, so a client that calls LAPACK
# through petscblaslapack.h must link it. petscvariables records it as BLASLAPACK_LIB; each -l entry
# is resolved to a full path in the -L directories that precede it, run-path flags are dropped, and
# every other item (absolute paths, -framework, -Wl,--start-group, -qmkl, ...) is kept in order.
set(PETSC_BLASLAPACK_LIBRARIES "")
if(PETSC_LIBRARY)
  get_filename_component(_petsc_lib_dir "${PETSC_LIBRARY}" DIRECTORY)
  set(_petsc_variables "${_petsc_lib_dir}/petsc/conf/petscvariables")
  if(EXISTS "${_petsc_variables}")
    file(STRINGS "${_petsc_variables}" _petsc_blaslapack REGEX "^BLASLAPACK_LIB =")
    string(REGEX REPLACE "^BLASLAPACK_LIB =[ ]*" "" _petsc_blaslapack "${_petsc_blaslapack}")
    separate_arguments(_petsc_blaslapack UNIX_COMMAND "${_petsc_blaslapack}")
    set(_petsc_blaslapack_dirs "")
    set(_petsc_blaslapack_framework FALSE)
    foreach(_item IN LISTS _petsc_blaslapack)
      if(_petsc_blaslapack_framework)
        list(APPEND PETSC_BLASLAPACK_LIBRARIES "-framework ${_item}")
        set(_petsc_blaslapack_framework FALSE)
      elseif(_item STREQUAL "-framework")
        set(_petsc_blaslapack_framework TRUE)
      elseif(_item MATCHES "^-Wl,-rpath")
        # Run paths come from the libraries' directories.
      elseif(_item MATCHES "^-L(.+)")
        list(APPEND _petsc_blaslapack_dirs "${CMAKE_MATCH_1}")
      elseif(_item MATCHES "^-l(.+)")
        set(_petsc_blaslapack_name "${CMAKE_MATCH_1}")
        unset(_petsc_blaslapack_path CACHE)
        find_library(_petsc_blaslapack_path NAMES ${_petsc_blaslapack_name}
                     PATHS ${_petsc_blaslapack_dirs} NO_DEFAULT_PATH)
        if(NOT _petsc_blaslapack_path)
          find_library(_petsc_blaslapack_path NAMES ${_petsc_blaslapack_name})
        endif()
        if(_petsc_blaslapack_path)
          list(APPEND PETSC_BLASLAPACK_LIBRARIES "${_petsc_blaslapack_path}")
        else()
          list(APPEND PETSC_BLASLAPACK_LIBRARIES "${_item}")
        endif()
        unset(_petsc_blaslapack_path CACHE)
      else()
        list(APPEND PETSC_BLASLAPACK_LIBRARIES "${_item}")
      endif()
    endforeach()
  endif()
endif()

if(PETSC_FOUND AND NOT TARGET PETSc::PETSc)
    add_library(PETSc::PETSc UNKNOWN IMPORTED)
    set_target_properties(PETSc::PETSc PROPERTIES
        IMPORTED_LOCATION "${PETSC_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${PETSC_INCLUDE_DIR}"
    )
endif()

# Without a recorded BLASLAPACK_LIB, fall back to a system LAPACK with PETSc's BLAS integer width,
# the width of PetscBLASInt. QUIET because clients of the installed module that do not call LAPACK
# need neither.
if(PETSC_FOUND AND NOT TARGET PETSc::BLASLAPACK)
  if(PETSC_BLASLAPACK_LIBRARIES)
    add_library(PETSc::BLASLAPACK INTERFACE IMPORTED)
    set_target_properties(PETSc::BLASLAPACK PROPERTIES
        INTERFACE_LINK_LIBRARIES "${PETSC_BLASLAPACK_LIBRARIES}"
    )
  else()
    set(_petsc_bla_sizeof_integer 4)
    if(EXISTS "${PETSC_INCLUDE_DIR}/petscconf.h")
      file(STRINGS "${PETSC_INCLUDE_DIR}/petscconf.h" _petsc_64bit_blas
           REGEX "^#define PETSC_HAVE_64BIT_BLAS_INDICES 1")
      if(_petsc_64bit_blas)
        set(_petsc_bla_sizeof_integer 8)
      endif()
    endif()
    if(DEFINED BLA_SIZEOF_INTEGER)
      set(_petsc_saved_bla_sizeof_integer "${BLA_SIZEOF_INTEGER}")
    endif()
    set(BLA_SIZEOF_INTEGER ${_petsc_bla_sizeof_integer})
    find_package(LAPACK QUIET)
    if(DEFINED _petsc_saved_bla_sizeof_integer)
      set(BLA_SIZEOF_INTEGER "${_petsc_saved_bla_sizeof_integer}")
      unset(_petsc_saved_bla_sizeof_integer)
    else()
      unset(BLA_SIZEOF_INTEGER)
    endif()
    if(LAPACK_FOUND)
      add_library(PETSc::BLASLAPACK INTERFACE IMPORTED)
      set_target_properties(PETSc::BLASLAPACK PROPERTIES
          INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}"
      )
      set(PETSC_BLASLAPACK_LIBRARIES ${LAPACK_LIBRARIES})
    endif()
  endif()
endif()
