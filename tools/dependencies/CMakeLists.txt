cmake_minimum_required(VERSION 3.20.2)
cmake_policy(SET CMP0135 NEW)
project(OpenSnDependencies)
include(ExternalProject)
include(ProcessorCount)
ProcessorCount(JOBS)
if(NOT JOBS OR JOBS EQUAL 0)
  set(JOBS 4)
endif()

if(NOT DEFINED CMAKE_INSTALL_PREFIX OR CMAKE_INSTALL_PREFIX STREQUAL "")
  message(FATAL_ERROR "You must set CMAKE_INSTALL_PREFIX explicitly. For example:\n"
                      "  cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ..")
endif()
message(STATUS "CMAKE_INSTALL_PREFIX is set to: ${CMAKE_INSTALL_PREFIX}")

set(OPENSN_ENV_DIR "${CMAKE_INSTALL_PREFIX}/bin")
set(ENV_SCRIPT "${OPENSN_ENV_DIR}/set_opensn_env.sh")

find_package(MPI REQUIRED)

option(ENABLE_BOOST   "Install Boost"   ON)
option(ENABLE_LUA     "Install Lua"     ON)
option(ENABLE_PETSC   "Install PETSc"   ON)
option(ENABLE_HDF5    "Install HDF5"    ON)
option(ENABLE_VTK     "Install VTK"     ON)
option(ENABLE_CALIPER "Install Caliper" ON)

# Only needed for Lua
if(APPLE)
  set(OS_TAG macosx)
else()
  set(OS_TAG linux)
endif()

# Boost install (header-only)
if(ENABLE_BOOST)
  find_package(Boost CONFIG QUIET)
  if (NOT Boost_FOUND)
    find_package(Boost QUIET)
  endif()

  if(Boost_FOUND)
    message(STATUS "Boost version: ${Boost_VERSION}")
    if(DEFINED Boost_INCLUDE_DIRS)
      message(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")
    elseif(DEFINED Boost_INCLUDE_DIR)
      message(STATUS "Boost include dir: ${Boost_INCLUDE_DIR}")
    endif()
  else()
    ExternalProject_Add(boost
      URL https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.gz
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory <SOURCE_DIR>/boost ${CMAKE_INSTALL_PREFIX}/include/boost
    )
  endif()
endif()

# Lua
if(ENABLE_LUA)
  ExternalProject_Add(lua
    URL https://www.lua.org/ftp/lua-5.4.6.tar.gz
    SOURCE_DIR ${CMAKE_BINARY_DIR}/lua-src
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ${CMAKE_COMMAND} -E chdir <SOURCE_DIR> make ${OS_TAG} MYCFLAGS=-fPIC MYLIBS=-lncurses -j${JOBS}
    INSTALL_COMMAND ${CMAKE_COMMAND} -E chdir <SOURCE_DIR> make install INSTALL_TOP=${CMAKE_INSTALL_PREFIX} 
  )
endif()

# PETSc
if(ENABLE_PETSC)
  ExternalProject_Add(petsc
    URL https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.23.0.tar.gz
    SOURCE_DIR ${CMAKE_BINARY_DIR}/petsc-src
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E chdir <SOURCE_DIR> ./configure --prefix=${CMAKE_INSTALL_PREFIX}
                      --download-hypre=1 --with-ssl=0 --with-debugging=0 --with-pic=1 --with-fc=0
                      --with-shared-libraries=1 --download-bison=1 --download-f2cblaslapack=1
                      --download-metis=1 --download-parmetis=1 --download-superlu_dist=1
                      --download-ptscotch=1 --with-cxx-dialect=C++11 --with-64-bit-indices
                      CC=${MPI_C_COMPILER} CXX=${MPI_CXX_COMPILER} COPTFLAGS="-O3" CXXOPTFLAGS="-O3"
    BUILD_COMMAND ${CMAKE_COMMAND} -E chdir <SOURCE_DIR> make -j${JOBS}
    INSTALL_COMMAND ${CMAKE_COMMAND} -E chdir <SOURCE_DIR> make install
  )
endif()

# HDF5
if(ENABLE_HDF5)
  ExternalProject_Add(hdf5
    URL https://support.hdfgroup.org/releases/hdf5/v1_14/v1_14_6/downloads/hdf5-1.14.6.tar.gz
    SOURCE_DIR ${CMAKE_BINARY_DIR}/hdf5-src
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E make_directory <SOURCE_DIR>/build && cd <SOURCE_DIR>/build &&
		      ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} ..
    BUILD_COMMAND ${CMAKE_COMMAND} -E chdir <SOURCE_DIR>/build make -j${JOBS}
    INSTALL_COMMAND ${CMAKE_COMMAND} -E chdir <SOURCE_DIR>/build make install
  )
endif()

# VTK
if(ENABLE_VTK)
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
endif()

# Caliper
if(ENABLE_CALIPER)
  ExternalProject_Add(caliper
    URL https://github.com/LLNL/Caliper/archive/refs/tags/v2.10.0.tar.gz
    SOURCE_DIR ${CMAKE_BINARY_DIR}/caliper-src
    CONFIGURE_COMMAND cmake -E make_directory <SOURCE_DIR>/build && cd <SOURCE_DIR>/build &&
                      cmake -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                      -DWITH_MPI=ON -DWITH_KOKKOS=OFF -DWITH_GOTCHA=OFF ..
    BUILD_COMMAND ${CMAKE_COMMAND} -E chdir <SOURCE_DIR>/build make -j${JOBS}
    INSTALL_COMMAND ${CMAKE_COMMAND} -E chdir <SOURCE_DIR>/build make install
  )
endif()

file(WRITE "${ENV_SCRIPT}" "export CMAKE_PREFIX_PATH=\"${CMAKE_INSTALL_PREFIX}\"\${CMAKE_PREFIX_PATH:+:\${CMAKE_PREFIX_PATH}}\n")
file(CHMOD ${ENV_SCRIPT} FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ)
