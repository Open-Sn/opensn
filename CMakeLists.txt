cmake_minimum_required(VERSION 3.12)

project(opensn
    VERSION 0.0.1
    LANGUAGES C CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

include(GNUInstallDirs)
include(CheckTypeSize)
include(CheckSymbolExists)

option(OPENSN_WITH_LUA "Build with lua support" ON)

# dependencies
find_package(MPI REQUIRED)

if(OPENSN_WITH_LUA)
    find_package(Lua 5.4 REQUIRED)
endif()

find_package(VTK QUIET)
if(VTK_VERSION VERSION_LESS "8.90.0")
    find_package(VTK COMPONENTS
        vtkCommonCore vtkCommonDataModel
        vtkIOLegacy vtkIOCore
        vtkIOXML vtkParallelCore vtkIOParallelXML
        vtkFiltersCore
        vtkIOEnSight
        vtkIOExodus
        REQUIRED)
else()
    find_package(VTK COMPONENTS
        CommonCore CommonDataModel
        IOLegacy IOCore
        IOXML ParallelCore IOParallelXML
        FiltersCore
        IOEnSight
        IOExodus
        REQUIRED)
endif()

find_package(PETSc 3.17 REQUIRED)
check_symbol_exists(PETSC_USE_64BIT_INDICES "${PETSC_INCLUDE_DIR}/petscconf.h" PETSC_USE_64BIT_INDICES)
if(NOT ${PETSC_USE_64BIT_INDICES} MATCHES 1)
    message(FATAL_ERROR "PETSc has not been configured with the flag --with-64-bit-indices\n")
endif()

# libopensn
file(GLOB_RECURSE LIBOPENSN_SRCS
    framework/*.cc
    modules/*.cc
)
add_library(libopensn SHARED ${LIBOPENSN_SRCS})

target_include_directories(libopensn
    PRIVATE
    $<INSTALL_INTERFACE:include/opensn>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}
    PUBLIC
    ${LUA_INCLUDE_DIR}
    ${PETSC_INCLUDE_DIR}
)

target_link_libraries(libopensn
    PRIVATE
    ${PETSC_LIBRARY}
    ${VTK_LIBRARIES}
    MPI::MPI_CXX
)
if(OPENSN_WITH_LUA)
    target_link_libraries(libopensn PRIVATE ${LUA_LIBRARIES})
    target_compile_definitions(libopensn PRIVATE OPENSN_WITH_LUA)
endif()

if(NOT MSVC)
    set_target_properties(libopensn PROPERTIES OUTPUT_NAME opensn)
endif()


if(OPENSN_WITH_LUA)
    add_subdirectory(lua)
    add_subdirectory(test)
endif()

#
configure_file(config.h.in config.h)
