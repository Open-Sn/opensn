# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

# Checks that each dependency floor in the Spack recipe is at least the floor in
# cmake/OpenSnDependencyVersions.cmake, so Spack cannot pick a version that CMake rejects.
#
#   cmake -P tools/developer/check-spack-floors.cmake

cmake_minimum_required(VERSION 3.20)

set(_root "${CMAKE_CURRENT_LIST_DIR}/../..")
include("${_root}/cmake/OpenSnDependencyVersions.cmake")
set(_recipe distribution/spack/packages/opensn/package.py)
file(READ "${_root}/${_recipe}" _recipe_text)
# CMake's ^ matches only at the start of the text, so a declaration is matched after a newline,
# which skips commented-out ones.
string(PREPEND _recipe_text "\n")

# Spack package name and the floor variable in OpenSnDependencyVersions.cmake. preCICE has no
# Spack dependency.
set(_pairs
  mpicpp-lite=OPENSN_MIN_MPICPP_LITE_VERSION
  boost=OPENSN_MIN_BOOST_VERSION
  petsc=OPENSN_MIN_PETSC_VERSION
  hdf5=OPENSN_MIN_HDF5_VERSION
  vtk=OPENSN_MIN_VTK_VERSION
  caliper=OPENSN_MIN_CALIPER_VERSION
  cuda=OPENSN_MIN_CUDA_VERSION
  python=OPENSN_MIN_PYTHON_VERSION
)

set(_failed FALSE)
foreach(_pair IN LISTS _pairs)
  string(REPLACE "=" ";" _pair "${_pair}")
  list(GET _pair 0 _package)
  list(GET _pair 1 _floor_var)
  set(_floor "${${_floor_var}}")
  if(NOT _recipe_text MATCHES "\n[ \t]*depends_on\\([ \t\r\n]*\"${_package}@([0-9.]+):")
    message(SEND_ERROR "${_recipe} has no depends_on(\"${_package}@<version>:\")")
    set(_failed TRUE)
    continue()
  endif()
  set(_spack_floor "${CMAKE_MATCH_1}")
  if(_spack_floor VERSION_LESS _floor)
    message(SEND_ERROR
      "Spack accepts ${_package}@${_spack_floor}, but OpenSn requires ${_floor} or newer "
      "(${_floor_var}). Raise the floor in ${_recipe}.")
    set(_failed TRUE)
  else()
    message(STATUS "${_package}: Spack ${_spack_floor}, OpenSn ${_floor}")
  endif()
endforeach()

if(_failed)
  message(FATAL_ERROR "${_recipe} does not match cmake/OpenSnDependencyVersions.cmake")
endif()
