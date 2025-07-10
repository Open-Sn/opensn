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

option(ENABLE_BOOST   "Install Boost"   ON)
option(ENABLE_PETSC   "Install PETSc"   ON)
option(ENABLE_HDF5    "Install HDF5"    ON)
option(ENABLE_VTK     "Install VTK"     ON)
option(ENABLE_CALIPER "Install Caliper" ON)

set(INSTALLED_DEP_COUNT 0)

# Boost install (header-only)
if(ENABLE_BOOST)
  find_package(Boost CONFIG QUIET)
  if (NOT Boost_FOUND)
    find_package(Boost QUIET)
  endif()
  if(Boost_FOUND)
    message(STATUS "Found Boost: ${Boost_DIR}")
  else()
    ExternalProject_Add(boost
      URL https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.gz
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory <SOURCE_DIR>/boost ${CMAKE_INSTALL_PREFIX}/include/boost
    )
    math(EXPR INSTALLED_DEP_COUNT "${INSTALLED_DEP_COUNT}+1")
  endif()
endif()

# PETSc
if(ENABLE_PETSC)
  find_package(PETSc 3.17 QUIET)
  if (PETSc_FOUND)
    check_symbol_exists(PETSC_USE_64BIT_INDICES "${PETSC_INCLUDE_DIR}/petscconf.h" PETSC_USE_64BIT_INDICES)
    check_symbol_exists(PETSC_HAVE_PTSCOTCH "${PETSC_INCLUDE_DIR}/petscconf.h" PETSC_HAVE_PTSCOTCH)
    if(PETSC_USE_64BIT_INDICES AND PETSC_HAVE_PTSCOTCH)
      message(STATUS "Found PETSc: ${PETSC_INCLUDE_DIR}")
      set(PETSC_CONFIG_OK TRUE)
    endif()
  endif()
  if (NOT PETSc_FOUND OR NOT PETSC_CONFIG_OK)
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
    )
    math(EXPR INSTALLED_DEP_COUNT "${INSTALLED_DEP_COUNT}+1")
  endif()
endif()

# HDF5
if(ENABLE_HDF5)
  find_package(HDF5 COMPONENTS C HL QUIET)
  if (HDF5_FOUND)
    message(STATUS "Found HDF5: ${HDF5_DIR}")
  else()
    ExternalProject_Add(hdf5
      URL https://support.hdfgroup.org/releases/hdf5/v1_14/v1_14_6/downloads/hdf5-1.14.6.tar.gz
      SOURCE_DIR ${CMAKE_BINARY_DIR}/hdf5-src
      CONFIGURE_COMMAND ${CMAKE_COMMAND} -E make_directory <SOURCE_DIR>/build && cd <SOURCE_DIR>/build &&
                        ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} ..
      BUILD_COMMAND ${CMAKE_COMMAND} -E chdir <SOURCE_DIR>/build make -j${JOBS}
      INSTALL_COMMAND ${CMAKE_COMMAND} -E chdir <SOURCE_DIR>/build make install
    )
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
    )
    math(EXPR INSTALLED_DEP_COUNT "${INSTALLED_DEP_COUNT}+1")
  endif()
endif()

# Caliper
if(ENABLE_CALIPER)
  find_package(caliper QUIET)
  if(caliper_FOUND)
    message(STATUS "Found Caliper: ${caliper_DIR}")
  else()
    ExternalProject_Add(caliper
      URL https://github.com/LLNL/Caliper/archive/refs/tags/v2.13.0.tar.gz
      SOURCE_DIR ${CMAKE_BINARY_DIR}/caliper-src
      CONFIGURE_COMMAND cmake -E make_directory <SOURCE_DIR>/build && cd <SOURCE_DIR>/build &&
                        cmake -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                        -DWITH_MPI=ON -DWITH_KOKKOS=OFF -DWITH_GOTCHA=OFF ..
      BUILD_COMMAND ${CMAKE_COMMAND} -E chdir <SOURCE_DIR>/build make -j${JOBS}
      INSTALL_COMMAND ${CMAKE_COMMAND} -E chdir <SOURCE_DIR>/build make install
    )
    math(EXPR INSTALLED_DEP_COUNT "${INSTALLED_DEP_COUNT}+1")
  endif()
endif()

if(INSTALLED_DEP_COUNT EQUAL 0)
  message(STATUS "All dependencies satisfied; no package installations necessary.")
else()
  add_custom_target(finalize ALL
    COMMAND ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} -P "${CMAKE_SOURCE_DIR}/env_script.cmake"
  )
endif()
