cmake_minimum_required(VERSION 3.12)

message(STATUS "Configuring ${CMAKE_BUILD_TYPE} build of ChiTech.")

if(UNIX AND NOT APPLE)
    add_definitions(-DUNIX_ENV)
elseif(APPLE)
    add_definitions(-DAPPLE_ENV)
    add_definitions(-DUNIX_ENV)
else()
    add_definitions(-DWINDOWS_ENV)
endif()

#------------------------------------------------ DEPENDENCIES
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
        set(LUA_ROOT "$ENV{LUA_ROOT}")
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

#================================================ Set project targets
SET(TARGET ChiTech)
project(${TARGET} LANGUAGES C CXX VERSION 1.0.2)

#================================================ Set variables
set(CHI_TECH_DIR "${PROJECT_SOURCE_DIR}")
set(CHI_TECH_MOD "modules")

configure_file("resources/chi_configuration.h.in"
               "${CHI_TECH_DIR}/framework/chi_configuration.h")

#================================================ Set cmake variables
find_package(MPI)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/resources/CMakeMacros")

if (NOT DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

    set(CONFIG_FILE_NAME "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/config.cmake")
    file(WRITE "${CONFIG_FILE_NAME}"
        "# Environment configuration file for ChiTech.\n")
endif()

#================================================ Export vars to config file
file(APPEND "${CONFIG_FILE_NAME}" "set(PETSC_ROOT ${PETSC_ROOT})\n")
file(APPEND "${CONFIG_FILE_NAME}" "set(LUA_ROOT ${LUA_ROOT})\n")
file(APPEND "${CONFIG_FILE_NAME}" "set(VTK_DIR ${VTK_DIR})\n")

#================================================ Include macros
include(GNUInstallDirs)
include(Filter)
include(Checks)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

#================================================ Default include directories
include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})
include_directories(SYSTEM "${LUA_ROOT}/include")
include_directories(SYSTEM "${PETSC_ROOT}/include")

include_directories(SYSTEM "${CHI_TECH_DIR}/external/paraview")

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

#================================================ Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
# Force -O3 in release builds (some OS's might downgrade it)
string(REPLACE "-O2" "-O3" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
string(REPLACE "-O2" "-O3" CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
string(REPLACE "-O2" "-O3" CMAKE_CXX_FLAGS_MINSIZEREL ${CMAKE_CXX_FLAGS_MINSIZEREL})

set(CHI_CXX_FLAGS)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    list(APPEND CHI_CXX_FLAGS "-pedantic")
    list(APPEND CHI_CXX_FLAGS "-Wall")
    list(APPEND CHI_CXX_FLAGS "-Wno-unused-variable")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")
    list(APPEND CHI_CXX_FLAGS "-pedantic")
    list(APPEND CHI_CXX_FLAGS "-Wall")
    list(APPEND CHI_CXX_FLAGS "-Wno-unused-variable")
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    list(APPEND CHI_CXX_FLAGS "-pedantic")
    list(APPEND CHI_CXX_FLAGS "-Wall")
    list(APPEND CHI_CXX_FLAGS "-Wno-unused-variable")
    list(APPEND CHI_CXX_FLAGS "-Wno-sign-compare")
    list(APPEND CHI_CXX_FLAGS "-Wno-psabi")
#elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
#elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
else()
    message(WARNING "Untested CMAKE_CXX_COMPILER_ID : ${CMAKE_CXX_COMPILER_ID}")
endif()

#================================================ Linker flags
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")

#================================================ Third-party library
add_subdirectory("${CHI_TECH_DIR}/external")

#================================================ Define source directories
set(SOURCES "${CHI_TECH_DIR}/framework/chi_runtime.cc"
            "${CHI_TECH_DIR}/framework/ChiObject.cc"
            "${CHI_TECH_DIR}/framework/ChiObjectFactory.cc")
add_subdirectory("${CHI_TECH_DIR}/framework/lua")
add_subdirectory("${CHI_TECH_DIR}/framework/console")
add_subdirectory("${CHI_TECH_DIR}/framework/math")
add_subdirectory("${CHI_TECH_DIR}/framework/physics")
add_subdirectory("${CHI_TECH_DIR}/framework/graphs")
add_subdirectory("${CHI_TECH_DIR}/framework/materials")
add_subdirectory("${CHI_TECH_DIR}/framework/parameters")
add_subdirectory("${CHI_TECH_DIR}/framework/mesh")
add_subdirectory("${CHI_TECH_DIR}/framework/mpi")
add_subdirectory("${CHI_TECH_DIR}/framework/logging")
add_subdirectory("${CHI_TECH_DIR}/framework/data_types")
add_subdirectory("${CHI_TECH_DIR}/framework/utils")
add_subdirectory("${CHI_TECH_DIR}/framework/interfaces")
add_subdirectory("${CHI_TECH_DIR}/framework/event_system")
add_subdirectory("${CHI_TECH_DIR}/framework/post_processors")

add_subdirectory("${CHI_TECH_MOD}")

add_subdirectory("${CHI_TECH_DIR}/test")

#================================================ Define targets
add_library(ChiLib SHARED ${SOURCES})
target_link_libraries(ChiLib ${CHI_LIBS} external)
target_compile_options(ChiLib PUBLIC ${CHI_CXX_FLAGS})

add_executable(${TARGET} "${CHI_TECH_DIR}/framework/chi_tech_main.cc")
add_executable(${TARGET}_test ${TEST_SOURCES})
#if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")

#else()
#    target_link_libraries(${TARGET} -Wl,--whole-archive ChiLib -Wl,--no-whole-archive )
#endif()

if(UNIX AND NOT APPLE)
    target_link_libraries(${TARGET} -Wl,--whole-archive ChiLib -Wl,--no-whole-archive )
    target_link_libraries(${TARGET}_test -Wl,--whole-archive ChiLib -Wl,--no-whole-archive )
elseif(APPLE)
    target_link_libraries(${TARGET} -Wl,-all_load ChiLib )
    target_link_libraries(${TARGET}_test -Wl,-all_load ChiLib )
endif()

target_compile_options(${TARGET} PUBLIC ${CHI_CXX_FLAGS})
target_compile_options(${TARGET}_test PUBLIC ${CHI_CXX_FLAGS})

set_target_properties( ${TARGET}_test
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/test/lib"
        LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/test/lib"
        RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/test/bin"
        )

# |------------ Write Makefile to root directory
file(WRITE ${PROJECT_SOURCE_DIR}/Makefile "subsystem:\n" "\t$(MAKE) -C chi_build \n\n"
        "clean:\n\t$(MAKE) -C chi_build clean\n")
