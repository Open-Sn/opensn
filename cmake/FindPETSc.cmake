# FindPETSc.cmake
#
# Once done this will define
#  PETSc_FOUND - System has PETSc and SLEPc
#  PETSC_FOUND - Compatibility alias for PETSc_FOUND
#  PETSC_INCLUDE_DIR - PETSc/SLEPc include directories
#  PETSC_LIBRARY - PETSc/SLEPc libraries or imported targets
#  PETSC_VERSION - The PETSc version
#  PETSC_CONF_HEADER - Full path to petscconf.h
#  SLEPC_FOUND - System has SLEPc
#  SLEPc_FOUND - Compatibility alias for SLEPC_FOUND
#  SLEPC_INCLUDE_DIR - The SLEPc include directory
#  SLEPC_LIBRARY - The SLEPc library or imported target
#  SLEPC_VERSION - The SLEPc version

include(FindPackageHandleStandardArgs)

set(_petsc_roots)
foreach(_var IN ITEMS PETSC_DIR PETSC_ROOT)
  if(DEFINED ENV{${_var}} AND NOT "$ENV{${_var}}" STREQUAL "")
    list(APPEND _petsc_roots "$ENV{${_var}}")
  endif()
endforeach()
list(REMOVE_DUPLICATES _petsc_roots)

set(_slepc_roots)
foreach(_var IN ITEMS SLEPC_DIR SLEPC_ROOT)
  if(DEFINED ENV{${_var}} AND NOT "$ENV{${_var}}" STREQUAL "")
    list(APPEND _slepc_roots "$ENV{${_var}}")
  endif()
endforeach()
list(REMOVE_DUPLICATES _slepc_roots)

set(_petsc_config_paths)
set(_petsc_include_hints)
set(_petsc_library_hints)
foreach(_root IN LISTS _petsc_roots)
  list(APPEND _petsc_config_paths "${_root}/lib/cmake/PETSc")
  list(APPEND _petsc_include_hints "${_root}/include")
  list(APPEND _petsc_library_hints "${_root}/lib")
  if(DEFINED ENV{PETSC_ARCH} AND NOT "$ENV{PETSC_ARCH}" STREQUAL "")
    list(APPEND _petsc_config_paths "${_root}/$ENV{PETSC_ARCH}/lib/cmake/PETSc")
    list(APPEND _petsc_include_hints "${_root}/$ENV{PETSC_ARCH}/include")
    list(APPEND _petsc_library_hints "${_root}/$ENV{PETSC_ARCH}/lib")
  endif()
endforeach()

if(DEFINED PETSC_INCLUDE_DIR AND NOT "${PETSC_INCLUDE_DIR}" MATCHES "-NOTFOUND$")
  list(APPEND _petsc_include_hints ${PETSC_INCLUDE_DIR})
endif()
if(DEFINED PETSC_LIBRARY AND NOT "${PETSC_LIBRARY}" MATCHES "-NOTFOUND$")
  foreach(_library IN LISTS PETSC_LIBRARY)
    if(EXISTS "${_library}")
      get_filename_component(_library_dir "${_library}" DIRECTORY)
      get_filename_component(_prefix_dir "${_library_dir}" DIRECTORY)
      list(APPEND _petsc_library_hints "${_library_dir}")
      list(APPEND _petsc_include_hints "${_prefix_dir}/include")
    endif()
  endforeach()
endif()
list(REMOVE_DUPLICATES _petsc_include_hints)
list(REMOVE_DUPLICATES _petsc_library_hints)

set(_slepc_config_paths)
set(_slepc_include_hints)
set(_slepc_library_hints)
foreach(_root IN LISTS _slepc_roots)
  list(APPEND _slepc_config_paths "${_root}/lib/cmake/SLEPc")
  list(APPEND _slepc_include_hints "${_root}/include")
  list(APPEND _slepc_library_hints "${_root}/lib")
  if(DEFINED ENV{SLEPC_ARCH} AND NOT "$ENV{SLEPC_ARCH}" STREQUAL "")
    list(APPEND _slepc_config_paths "${_root}/$ENV{SLEPC_ARCH}/lib/cmake/SLEPc")
    list(APPEND _slepc_include_hints "${_root}/$ENV{SLEPC_ARCH}/include")
    list(APPEND _slepc_library_hints "${_root}/$ENV{SLEPC_ARCH}/lib")
  endif()
endforeach()
foreach(_root IN LISTS _petsc_roots)
  list(APPEND _slepc_include_hints "${_root}/include/slepc")
  list(APPEND _slepc_library_hints "${_root}/lib")
  if(DEFINED ENV{PETSC_ARCH} AND NOT "$ENV{PETSC_ARCH}" STREQUAL "")
    list(APPEND _slepc_include_hints "${_root}/$ENV{PETSC_ARCH}/include/slepc")
    list(APPEND _slepc_library_hints "${_root}/$ENV{PETSC_ARCH}/lib")
  endif()
endforeach()

if(DEFINED SLEPC_INCLUDE_DIR AND NOT "${SLEPC_INCLUDE_DIR}" MATCHES "-NOTFOUND$")
  list(APPEND _slepc_include_hints ${SLEPC_INCLUDE_DIR})
endif()
if(DEFINED SLEPC_LIBRARY AND NOT "${SLEPC_LIBRARY}" MATCHES "-NOTFOUND$")
  foreach(_library IN LISTS SLEPC_LIBRARY)
    if(EXISTS "${_library}")
      get_filename_component(_library_dir "${_library}" DIRECTORY)
      get_filename_component(_prefix_dir "${_library_dir}" DIRECTORY)
      list(APPEND _slepc_library_hints "${_library_dir}")
      list(APPEND _slepc_include_hints "${_prefix_dir}/include")
    endif()
  endforeach()
endif()
list(REMOVE_DUPLICATES _slepc_include_hints)
list(REMOVE_DUPLICATES _slepc_library_hints)

find_package(PETSc CONFIG QUIET PATHS ${_petsc_config_paths})
find_package(SLEPc CONFIG QUIET PATHS ${_slepc_config_paths})
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(PC_PETSC QUIET IMPORTED_TARGET PETSc)
  pkg_check_modules(PC_SLEPC QUIET IMPORTED_TARGET SLEPc)
  if(PC_PETSC_FOUND)
    list(APPEND _petsc_include_hints ${PC_PETSC_INCLUDE_DIRS})
    list(APPEND _petsc_library_hints ${PC_PETSC_LIBRARY_DIRS})
  endif()
  if(PC_SLEPC_FOUND)
    list(APPEND _slepc_include_hints ${PC_SLEPC_INCLUDE_DIRS})
    list(APPEND _slepc_library_hints ${PC_SLEPC_LIBRARY_DIRS})
  endif()
endif()

function(_opensn_append_target_property output target property)
  if(TARGET "${target}")
    get_target_property(_property_value "${target}" "${property}")
    if(_property_value AND NOT _property_value STREQUAL "_property_value-NOTFOUND")
      set(${output} ${${output}} ${_property_value} PARENT_SCOPE)
    endif()
  endif()
endfunction()

function(_opensn_find_imported_target output)
  foreach(_target IN LISTS ARGN)
    if(TARGET "${_target}")
      set(${output} "${_target}" PARENT_SCOPE)
      return()
    endif()
  endforeach()
  set(${output} "" PARENT_SCOPE)
endfunction()

_opensn_find_imported_target(_petsc_target PETSc::PETSc petsc PkgConfig::PC_PETSC)
_opensn_find_imported_target(_slepc_target SLEPc::SLEPc SLEPc::slepc slepc PkgConfig::PC_SLEPC)

if(_petsc_target)
  set(PETSC_LIBRARY "${_petsc_target}")
  _opensn_append_target_property(PETSC_INCLUDE_DIR "${_petsc_target}" INTERFACE_INCLUDE_DIRECTORIES)
  _opensn_append_target_property(PETSC_LIBRARY "${_petsc_target}" INTERFACE_LINK_LIBRARIES)
  get_target_property(_petsc_version "${_petsc_target}" INTERFACE_VERSION)
  if(_petsc_version AND NOT _petsc_version STREQUAL "_petsc_version-NOTFOUND")
    set(PETSC_VERSION "${_petsc_version}")
  endif()
else()
  find_path(PETSC_BASE_INCLUDE_DIR
    NAMES petsc.h
    HINTS ${_petsc_include_hints}
    DOC "PETSc base include directory"
  )
  find_library(PETSC_BASE_LIBRARY
    NAMES petsc
    HINTS ${_petsc_library_hints}
    DOC "PETSc library"
  )
  set(PETSC_INCLUDE_DIR ${PETSC_BASE_INCLUDE_DIR})
  set(PETSC_LIBRARY ${PETSC_BASE_LIBRARY})
endif()

foreach(_include_dir IN LISTS PETSC_INCLUDE_DIR)
  if(IS_DIRECTORY "${_include_dir}")
    list(APPEND _petsc_include_hints "${_include_dir}")
    list(APPEND _slepc_include_hints "${_include_dir}")
    if(IS_DIRECTORY "${_include_dir}/slepc")
      list(APPEND _slepc_include_hints "${_include_dir}/slepc")
    endif()
  endif()
endforeach()
list(REMOVE_DUPLICATES _petsc_include_hints)
list(REMOVE_DUPLICATES _slepc_include_hints)

find_path(PETSC_ARCH_INCLUDE_DIR
  NAMES petscconf.h
  HINTS ${_petsc_include_hints}
  DOC "PETSc architecture include directory"
)
find_file(PETSC_CONF_HEADER
  NAMES petscconf.h
  HINTS ${_petsc_include_hints}
  DOC "PETSc configuration header"
)
find_file(PETSCVERSION_H
  NAMES petscversion.h
  HINTS ${_petsc_include_hints}
  DOC "PETSc version header"
)

list(APPEND PETSC_INCLUDE_DIR ${PETSC_ARCH_INCLUDE_DIR})

if(NOT DEFINED PETSC_VERSION)
  set(PETSC_VERSION "unknown")
  if(PETSCVERSION_H)
    file(READ "${PETSCVERSION_H}" _petscver)
    string(REGEX MATCH "#[ \t]*define[ \t]+PETSC_VERSION_MAJOR[ \t]+([0-9]+)" _major_match
                 "${_petscver}")
    set(_petsc_major "${CMAKE_MATCH_1}")
    string(REGEX MATCH "#[ \t]*define[ \t]+PETSC_VERSION_MINOR[ \t]+([0-9]+)" _minor_match
                 "${_petscver}")
    set(_petsc_minor "${CMAKE_MATCH_1}")
    string(REGEX MATCH "#[ \t]*define[ \t]+PETSC_VERSION_SUBMINOR[ \t]+([0-9]+)" _patch_match
                 "${_petscver}")
    set(_petsc_patch "${CMAKE_MATCH_1}")
    if(_major_match AND _minor_match AND _patch_match)
      set(PETSC_VERSION "${_petsc_major}.${_petsc_minor}.${_petsc_patch}")
    endif()
  endif()
endif()

foreach(_library IN LISTS PETSC_LIBRARY)
  if(EXISTS "${_library}")
    get_filename_component(_library_dir "${_library}" DIRECTORY)
    get_filename_component(_prefix_dir "${_library_dir}" DIRECTORY)
    list(APPEND _slepc_library_hints "${_library_dir}")
    list(APPEND _slepc_include_hints "${_prefix_dir}/include")
    if(IS_DIRECTORY "${_prefix_dir}/include/slepc")
      list(APPEND _slepc_include_hints "${_prefix_dir}/include/slepc")
    endif()
  endif()
endforeach()
list(REMOVE_DUPLICATES _slepc_include_hints)
list(REMOVE_DUPLICATES _slepc_library_hints)

if(_slepc_target)
  set(SLEPC_LIBRARY "${_slepc_target}")
  _opensn_append_target_property(SLEPC_INCLUDE_DIR "${_slepc_target}" INTERFACE_INCLUDE_DIRECTORIES)
  _opensn_append_target_property(SLEPC_LIBRARY "${_slepc_target}" INTERFACE_LINK_LIBRARIES)
else()
  find_path(SLEPC_INCLUDE_DIR
    NAMES slepc.h
    HINTS ${_slepc_include_hints}
    DOC "SLEPc include directory"
  )
  find_library(SLEPC_LIBRARY
    NAMES slepc
    HINTS ${_slepc_library_hints}
    DOC "SLEPc library"
  )
endif()

find_file(SLEPCVERSION_H
  NAMES slepcversion.h
  HINTS ${_slepc_include_hints}
  DOC "SLEPc version header"
)

set(SLEPC_VERSION "unknown")
if(SLEPCVERSION_H)
  file(READ "${SLEPCVERSION_H}" _slepcver)
  string(REGEX MATCH "#[ \t]*define[ \t]+SLEPC_VERSION_MAJOR[ \t]+([0-9]+)" _major_match
               "${_slepcver}")
  set(_slepc_major "${CMAKE_MATCH_1}")
  string(REGEX MATCH "#[ \t]*define[ \t]+SLEPC_VERSION_MINOR[ \t]+([0-9]+)" _minor_match
               "${_slepcver}")
  set(_slepc_minor "${CMAKE_MATCH_1}")
  string(REGEX MATCH "#[ \t]*define[ \t]+SLEPC_VERSION_SUBMINOR[ \t]+([0-9]+)" _patch_match
               "${_slepcver}")
  set(_slepc_patch "${CMAKE_MATCH_1}")
  if(_major_match AND _minor_match AND _patch_match)
    set(SLEPC_VERSION "${_slepc_major}.${_slepc_minor}.${_slepc_patch}")
  endif()
endif()

list(APPEND PETSC_INCLUDE_DIR ${SLEPC_INCLUDE_DIR})
list(APPEND PETSC_LIBRARY ${SLEPC_LIBRARY})
list(REMOVE_DUPLICATES PETSC_INCLUDE_DIR)
list(REMOVE_DUPLICATES PETSC_LIBRARY)

find_package_handle_standard_args(
  SLEPc
  REQUIRED_VARS SLEPC_LIBRARY SLEPC_INCLUDE_DIR
  VERSION_VAR SLEPC_VERSION
  NAME_MISMATCHED
)
set(SLEPC_FOUND ${SLEPc_FOUND})

if(NOT SLEPc_FOUND)
  # OpenSn requires SLEPc wherever PETSc is required.
  set(PETSC_INCLUDE_DIR PETSC_INCLUDE_DIR-NOTFOUND)
  set(PETSC_LIBRARY PETSC_LIBRARY-NOTFOUND)
  if(PETSc_FIND_REQUIRED)
    message(FATAL_ERROR "PETSc/SLEPc support is required, but SLEPc was not found.")
  endif()
endif()

find_package_handle_standard_args(
  PETSc
  REQUIRED_VARS PETSC_LIBRARY PETSC_INCLUDE_DIR PETSC_CONF_HEADER
  VERSION_VAR PETSC_VERSION
)
set(PETSC_FOUND ${PETSc_FOUND})
