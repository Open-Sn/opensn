#!/bin/sh

# 
# SYNOPSIS:
#   ./configure.sh [clean] [Debug|Release|RelWithDebInfo|MinSizeRel] [...]
#
# DESCRIPTION:
#   This script generates the CMake build scripts. The precise behavior of this
#   script depends on the arguments passed:
#   
#   1. If called without arguments, CMake scripts for a release (optimized) build
#      of ChiTech will be generated in the directory chi_build. Additional
#      arguments will be passed to and interpreted by the CMake generator.
#
#   2. If called with the single argument from the following list:
#         Debug, Release, RelWithDebInfo, MinSizeRel,
#      scripts for the corresponding build of ChiTech will be generated in
#      chi_build. See the CMake documentation (CMAKE_BUILD_TYPE) for more
#      information about the available build types.
#
#   3. If called with the single argument 'clean', the CMake build directory
#      (chi_build) will be removed and the script will stop. If additional
#      arguments are supplied after 'clean', after removing the build directory, 
#      the script will behave as in case 2.
#
# EXAMPLES:
#   ./configure.sh clean Release
#      regenerate the CMake scripts for a Release build of ChiTech from scratch
#      (invalidating CMake cache that stores previously used settings)
#
#   ./configure.sh -DCMAKE_VERBOSE_MAKEFILE=1
#      generate the CMake scripts for a Release build of ChiTech (using CMake
#      cache if available); the 'make' command will be verbose_info (showing, e.g.,
#      what compiler flags are being used)
#
#   ./configure.sh Debug -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=$PWD/bin/Debug
#      generate the CMake scripts for a Debug build of ChiTech (using CMake
#      cache if available); the 'make' command will place the ChiTech executable
#      into the specified directory (instead of the default 'bin' subdirectory of
#      ChiTech root directory).
#

is_cmake_build_type() {
  for bt in Debug Release RelWithDebInfo MinSizeRel; do
    if [ "$1" = $bt ]; then
      return 0
    fi
  done
  return 1
}

CMAKE_ARGS=                 # additional CMake arguments
DO_CLEAN="No"               # if yes, remove the chi_build directory before generating CMake
DO_CMAKE="Yes"              # by default, run CMake generator after configuring
CMAKE_BUILD_TYPE="Release"  # by default, configure the Release build

# Go through all arguments
for arg in "$@"; do

  if [ "$arg" = "clean" ]; then
    # The build directory will be cleaned
    DO_CLEAN="Yes"
    
    if [ $# -eq 1 ]; then
      # If "clean" is the only argument passed
      DO_CMAKE="No"
    fi

  elif is_cmake_build_type "$arg"; then
    CMAKE_BUILD_TYPE="$arg"

  else
    # Additional CMake arguments (we handle -DCMAKE_BUILD_TYPE ourselves)
    case "$arg" in
      -DCMAKE_BUILD_TYPE*) CMAKE_BUILD_TYPE=$(cut -d'=' -f1 "$arg") ;;
      *)  CMAKE_ARGS="${CMAKE_ARGS} ${arg}" ;;
    esac
  fi

done

if [ $DO_CLEAN = "Yes" ]; then
  # If chi_build directory exists, remove it
  if [ -d "chi_build" ]; then
    rm -r chi_build
  fi
fi

if [ $DO_CMAKE = "Yes" ]; then
  mkdir -p chi_build && cd chi_build && \
  cmake -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} ${CMAKE_ARGS} ../ && cd ..
fi
