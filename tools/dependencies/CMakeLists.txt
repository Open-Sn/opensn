cmake_minimum_required(VERSION 3.20.2)
cmake_policy(SET CMP0135 NEW)
project(OpenSnDependencies)
include(ExternalProject)
include(ProcessorCount)
include(CheckSymbolExists)
ProcessorCount(JOBS)
if(NOT JOBS OR JOBS EQUAL 0)
  set(JOBS 4)
endif()

if(NOT DEFINED CMAKE_INSTALL_PREFIX OR CMAKE_INSTALL_PREFIX STREQUAL "")
  message(FATAL_ERROR "You must set CMAKE_INSTALL_PREFIX explicitly. For example:\n"
                      "  cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ..")
endif()
message(STATUS "CMAKE_INSTALL_PREFIX is set to: ${CMAKE_INSTALL_PREFIX}")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../../cmake")

set(OPENSN_ENV_DIR "${CMAKE_INSTALL_PREFIX}/bin")
set(ENV_SCRIPT "${OPENSN_ENV_DIR}/set_opensn_env.sh")

find_package(MPI REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(mpicpp-lite 2.2 QUIET)

option(ENABLE_BOOST   "Install Boost"   ON)
option(ENABLE_PETSC   "Install PETSc"   ON)
option(ENABLE_HDF5    "Install HDF5"    ON)
option(ENABLE_VTK     "Install VTK"     ON)
option(ENABLE_CALIPER "Install Caliper" ON)
if(ENABLE_VTK AND NOT ENABLE_HDF5)
  message(STATUS "ENABLE_VTK=ON requires HDF5; enabling HDF5")
  set(ENABLE_HDF5 ON CACHE BOOL "" FORCE)
endif()

set(INSTALLED_DEP_COUNT 0)
set(ALL_DEP_TARGETS)

# mpicpp-lite install
if (TARGET mpicpp-lite::mpicpp-lite)
  message(STATUS "Found mpicpp-lite: " ${MPICPP_LITE_VERSION})
  if(NOT TARGET mpicpp-lite)
    add_custom_target(mpicpp-lite)
  endif()
  list(APPEND ALL_DEP_TARGETS mpicpp-lite)
else()
  message(STATUS "No")
  ExternalProject_Add(mpicpp-lite
    URL https://github.com/andrsd/mpicpp-lite/archive/refs/tags/v2.2.0.tar.gz
    CONFIGURE_COMMAND ${CMAKE_COMMAND}
      -S <SOURCE_DIR>
      -B <BINARY_DIR>
      -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      -DCMAKE_BUILD_TYPE=Release
      <SOURCE_DIR>
    BUILD_COMMAND ""
    INSTALL_COMMAND
      ${CMAKE_COMMAND} --build <BINARY_DIR> --target install --
    LOG_DOWNLOAD  1
    LOG_UPDATE    1
    LOG_CONFIGURE 1
    LOG_BUILD     1
    LOG_INSTALL   1
    LOG_TEST      1
    LOG_OUTPUT_ON_FAILURE 1
  )
  list(APPEND ALL_DEP_TARGETS mpicpp-lite)
  math(EXPR INSTALLED_DEP_COUNT "${INSTALLED_DEP_COUNT}+1")
endif()

# Boost install (header-only)
if(ENABLE_BOOST)
  find_package(Boost CONFIG QUIET)
  if (NOT Boost_FOUND)
    find_package(Boost QUIET)
  endif()
  if(Boost_FOUND)
    message(STATUS "Found Boost: ${Boost_DIR}")
    if(NOT TARGET boost)
      add_custom_target(boost)
    endif()
    list(APPEND ALL_DEP_TARGETS boost)
  else()
    ExternalProject_Add(boost
      URL https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.gz
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory <SOURCE_DIR>/boost ${CMAKE_INSTALL_PREFIX}/include/boost
      LOG_DOWNLOAD  1
      LOG_UPDATE    1
      LOG_CONFIGURE 1
      LOG_BUILD     1
      LOG_INSTALL   1
      LOG_TEST      1
      LOG_OUTPUT_ON_FAILURE 1
    )
    list(APPEND ALL_DEP_TARGETS boost)
    math(EXPR INSTALLED_DEP_COUNT "${INSTALLED_DEP_COUNT}+1")
  endif()
endif()

# PETSc
if(ENABLE_PETSC)
  find_package(PETSc 3.17 QUIET)
  if (PETSC_FOUND)
    check_symbol_exists(PETSC_USE_64BIT_INDICES "${PETSC_INCLUDE_DIR}/petscconf.h" PETSC_USE_64BIT_INDICES)
    check_symbol_exists(PETSC_HAVE_PTSCOTCH "${PETSC_INCLUDE_DIR}/petscconf.h" PETSC_HAVE_PTSCOTCH)
    if(PETSC_USE_64BIT_INDICES AND PETSC_HAVE_PTSCOTCH)
      message(STATUS "Found PETSc: ${PETSC_INCLUDE_DIR}")
      set(PETSC_CONFIG_OK TRUE)
      if(NOT TARGET petsc)
        add_custom_target(petsc)
      endif()
      list(APPEND ALL_DEP_TARGETS petsc)
    endif()
  endif()
  if (NOT PETSC_FOUND OR NOT PETSC_CONFIG_OK)
    ExternalProject_Add(petsc
      URL https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.23.0.tar.gz
      SOURCE_DIR ${CMAKE_BINARY_DIR}/petsc-src
      CONFIGURE_COMMAND ${CMAKE_COMMAND} -E chdir <SOURCE_DIR> ./configure --prefix=${CMAKE_INSTALL_PREFIX}
                        --download-hypre=yes --with-ssl=no --with-debugging=no --with-pic=yes --with-fc=0
                        --with-shared-libraries=yes --download-bison=yes --download-f2cblaslapack=yes
                        --download-metis=yes --download-parmetis=yes --download-superlu_dist=yes
                        --download-ptscotch=yes --download-slepc=yes --with-cxx-dialect=C++11
                        --with-64-bit-indices=yes --CC=${MPI_C_COMPILER} --CXX=${MPI_CXX_COMPILER}
                        --CFLAGS+=-O3 --CXXFLAGS+=-O3
      BUILD_COMMAND ${CMAKE_COMMAND} -E chdir <SOURCE_DIR> make -j${JOBS}
      INSTALL_COMMAND ${CMAKE_COMMAND} -E chdir <SOURCE_DIR> make install
      LOG_DOWNLOAD  1
      LOG_UPDATE    1
      LOG_CONFIGURE 1
      LOG_BUILD     1
      LOG_INSTALL   1
      LOG_TEST      1
      LOG_OUTPUT_ON_FAILURE 1
    )
    list(APPEND ALL_DEP_TARGETS petsc)
    math(EXPR INSTALLED_DEP_COUNT "${INSTALLED_DEP_COUNT}+1")
  endif()
endif()

# HDF5
if(ENABLE_HDF5)
  find_package(HDF5 COMPONENTS C HL QUIET)
  if (HDF5_FOUND)
    message(STATUS "Found HDF5: ${HDF5_VERSION}")
    if(NOT TARGET hdf5)
      add_custom_target(hdf5)
    endif()
    list(APPEND ALL_DEP_TARGETS hdf5)
  else()
    ExternalProject_Add(hdf5
      URL https://support.hdfgroup.org/releases/hdf5/v1_14/v1_14_6/downloads/hdf5-1.14.6.tar.gz
      SOURCE_DIR ${CMAKE_BINARY_DIR}/hdf5-src
      BINARY_DIR ${CMAKE_BINARY_DIR}/hdf5-build
      CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      BUILD_COMMAND
        ${CMAKE_COMMAND} --build <BINARY_DIR> -- -j${JOBS}
      INSTALL_COMMAND
        ${CMAKE_COMMAND} --build <BINARY_DIR> --target install
      LOG_DOWNLOAD  1
      LOG_UPDATE    1
      LOG_CONFIGURE 1
      LOG_BUILD     1
      LOG_INSTALL   1
      LOG_TEST      1
      LOG_OUTPUT_ON_FAILURE 1
    )
    list(APPEND ALL_DEP_TARGETS hdf5)
    math(EXPR INSTALLED_DEP_COUNT "${INSTALLED_DEP_COUNT}+1")
  endif()
endif()

# VTK
if(ENABLE_VTK)
  find_package(VTK COMPONENTS
               CommonCore
               CommonDataModel
               IOLegacy
               IOCore
               IOXML
               ParallelCore
               IOParallelXML
               FiltersCore
               IOEnSight
               IOExodus
               QUIET
              )
  if (VTK_FOUND)
    message(STATUS "Found VTK: ${VTK_DIR}")
    if(NOT TARGET vtk)
      add_custom_target(vtk)
    endif()
    list(APPEND ALL_DEP_TARGETS vtk)
  else()
    ExternalProject_Add(vtk
      URL https://www.vtk.org/files/release/9.3/VTK-9.3.0.tar.gz
      SOURCE_DIR ${CMAKE_BINARY_DIR}/vtk-src
      CONFIGURE_COMMAND ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                        -DBUILD_SHARED_LIBS=ON -DVTK_USE_MPI=ON
                        -DVTK_GROUP_ENABLE_StandAlone=WANT -DVTK_GROUP_ENABLE_Rendering=DONT_WANT
                        -DVTK_GROUP_ENABLE_Imaging=DONT_WANT -DVTK_GROUP_ENABLE_Web=DONT_WANT
                        -DVTK_GROUP_ENABLE_Qt=DONT_WANT -DVTK_MODULE_USE_EXTERNAL_VTK_hdf5=ON
                        -DCMAKE_BUILD_TYPE=Release <SOURCE_DIR>
      BUILD_COMMAND make -j${JOBS}
      INSTALL_COMMAND make install
      DEPENDS hdf5
      LOG_DOWNLOAD  1
      LOG_UPDATE    1
      LOG_CONFIGURE 1
      LOG_BUILD     1
      LOG_INSTALL   1
      LOG_TEST      1
      LOG_OUTPUT_ON_FAILURE 1
    )
    list(APPEND ALL_DEP_TARGETS vtk)
    math(EXPR INSTALLED_DEP_COUNT "${INSTALLED_DEP_COUNT}+1")
  endif()
endif()

# Caliper
if(ENABLE_CALIPER)
  find_package(caliper QUIET)
  if(caliper_FOUND)
    message(STATUS "Found Caliper: ${caliper_DIR}")
    if(NOT TARGET caliper)
      add_custom_target(caliper)
    endif()
    list(APPEND ALL_DEP_TARGETS caliper)
  else()
    ExternalProject_Add(caliper
      URL https://github.com/LLNL/Caliper/archive/refs/tags/v2.13.0.tar.gz
      SOURCE_DIR ${CMAKE_BINARY_DIR}/caliper-src
      CONFIGURE_COMMAND ${CMAKE_COMMAND}
        -S <SOURCE_DIR>
        -B <BINARY_DIR>
        -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
        -DWITH_MPI=ON -DWITH_KOKKOS=OFF -DWITH_GOTCHA=OFF
      BUILD_COMMAND   ${CMAKE_COMMAND} --build <BINARY_DIR> -- -j${JOBS}
      INSTALL_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --target install --
      LOG_DOWNLOAD  1
      LOG_UPDATE    1
      LOG_CONFIGURE 1
      LOG_BUILD     1
      LOG_INSTALL   1
      LOG_TEST      1
    )
    list(APPEND ALL_DEP_TARGETS caliper)
    math(EXPR INSTALLED_DEP_COUNT "${INSTALLED_DEP_COUNT}+1")
  endif()
endif()

if(ALL_DEP_TARGETS)
  add_custom_target(all_deps_done ALL DEPENDS ${ALL_DEP_TARGETS})
else()
  add_custom_target(all_deps_done ALL)
endif()

if(INSTALLED_DEP_COUNT EQUAL 0)
  message(STATUS "All dependencies satisfied; no package installations necessary.")
endif()

add_custom_target(finalize ALL
  DEPENDS all_deps_done
  COMMAND ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} -P "${CMAKE_SOURCE_DIR}/env_script.cmake"
)
