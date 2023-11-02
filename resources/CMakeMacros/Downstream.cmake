message (STATUS "Loading Downstream.cmake")

if(UNIX AND NOT APPLE)
    add_definitions(-DUNIX_ENV)
elseif(APPLE)
    add_definitions(-DAPPLE_ENV)
    add_definitions(-DUNIX_ENV)
else()
    add_definitions(-DWINDOWS_ENV)
endif()

#------------------------------------------------ DEPENDENCIES
if (NOT DEFINED CHI_TECH_DIR)
    if (NOT (DEFINED ENV{CHI_TECH_DIR}))
        message(FATAL_ERROR "***** CHI_TECH_DIR is not set *****")
    else()
        set(CHI_TECH_DIR    "$ENV{CHI_TECH_DIR}")
    endif()
endif()
message(STATUS "CHI_TECH_DIR set to ${CHI_TECH_DIR}")

include("${CHI_TECH_DIR}/bin/config.cmake")

if (NOT DEFINED PETSC_ROOT)
    if (NOT (DEFINED ENV{PETSC_ROOT}))
        message(FATAL_ERROR "***** PETSC_ROOT is not set *****")
    else()
        set(PETSC_ROOT "$ENV{PETSC_ROOT}")
    endif()
endif()
message(STATUS "PETSC_ROOT set to ${PETSC_ROOT}")

if (NOT DEFINED LUA_ROOT)
    if (NOT (DEFINED ENV{LUA_ROOT}))
        message(FATAL_ERROR "***** LUA_ROOT is not set *****")
    else()
        set(LUA_ROOT    "$ENV{LUA_ROOT}")
    endif()
endif()
message(STATUS "LUA_ROOT set to ${LUA_ROOT}")

if (NOT DEFINED VTK_DIR)
    if (NOT (DEFINED ENV{VTK_DIR}))
        message(FATAL_ERROR "***** VTK_DIR is not set *****")
    else()
        set(VTK_DIR "$ENV{VTK_DIR}")
    endif()
endif()
message(STATUS "VTK_DIR set to ${VTK_DIR}")

find_package(MPI)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    "${CHI_TECH_DIR}/resources/CMakeMacros")

#================================================ Include macros
include(GNUInstallDirs)
include(Filter)
include(Checks)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

#================================================ Include directories
include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})
include_directories(SYSTEM "${LUA_ROOT}/include")
include_directories(SYSTEM "${PETSC_ROOT}/include")

include_directories(SYSTEM "${CHI_TECH_DIR}/ThirdParty/paraview")

include_directories("${CHI_TECH_DIR}/framework")
include_directories("${CHI_TECH_DIR}/framework/lua")
include_directories("${CHI_TECH_DIR}/framework/mpi")
include_directories("${CHI_TECH_DIR}/framework/logging")
include_directories("${CHI_TECH_DIR}/framework/utils")
include_directories("${CHI_TECH_DIR}/resources")
include_directories("${CHI_TECH_DIR}/modules")
include_directories("${CHI_TECH_DIR}/modules/LinearBoltzmannSolvers")

#================================================ Library directories
link_directories("${LUA_ROOT}/lib")
link_directories("${PETSC_ROOT}/lib")
link_directories("${CHI_TECH_DIR}/lib")

# --------------------------- VTK
find_package(VTK PATHS ${VTK_DIR} QUIET)
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
if (NOT VTK_FOUND)
    message(FATAL_ERROR "VTK not found: ${VTK_NOT_FOUND_MESSAGE}")
endif()

if (VTK_VERSION VERSION_LESS "8.90.0")
    find_package(VTK COMPONENTS
            vtkCommonCore vtkCommonDataModel
            vtkIOLegacy vtkIOCore
            vtkIOXML vtkParallelCore vtkIOParallelXML
            vtkFiltersCore
            vtkIOEnSight
            vtkIOExodus
            REQUIRED PATHS ${VTK_DIR})
    # old system
    include(${VTK_USE_FILE})
    include_directories(SYSTEM ${VTK_INCLUDE_DIRS})
else ()
    find_package(VTK COMPONENTS
            CommonCore CommonDataModel
            IOLegacy IOCore
            IOXML ParallelCore IOParallelXML
            FiltersCore
            IOEnSight
            IOExodus
            REQUIRED PATHS ${VTK_DIR})
    # vtk_module_autoinit is needed
    vtk_module_autoinit(TARGETS ${TARGET} MODULES ${VTK_LIBRARIES})
endif()

set(CHI_LIBS stdc++ lua m dl ${MPI_CXX_LIBRARIES} petsc ${VTK_LIBRARIES})
set(CHI_LIBS ${CHI_LIBS} external)

#================================================ Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")

#================================================ Linker flags
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")
