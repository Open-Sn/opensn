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

# dependencies
find_package(MPI REQUIRED)

find_package(Lua 5.4 REQUIRED)

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
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}
    PUBLIC
    ${LUA_INCLUDE_DIR}
    ${PETSC_INCLUDE_DIR}
)

target_link_libraries(libopensn
    PRIVATE
    ${LUA_LIBRARIES}
    ${PETSC_LIBRARY}
    ${VTK_LIBRARIES}
    MPI::MPI_CXX
)

if(NOT MSVC)
    set_target_properties(libopensn PROPERTIES OUTPUT_NAME opensn)
endif()

add_subdirectory(lua)
add_subdirectory(test)

#
configure_file(config.h.in include/opensn/config.h)

# Copy header files into "CMAKE_CURRENT_BINARY_DIR/opensn"
file(GLOB_RECURSE OPENSN_HEADERS framework/*.h modules/*.h)
foreach(f ${OPENSN_HEADERS})
    get_filename_component(destdir ${f} DIRECTORY)
    string(REGEX REPLACE "^${CMAKE_CURRENT_SOURCE_DIR}/" "" destdir ${destdir})
    get_filename_component(destfile ${f} NAME)
    configure_file(${f} include/opensn/${destdir}/${destfile} COPYONLY)
endforeach()
