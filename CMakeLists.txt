cmake_minimum_required(VERSION 3.25.2)

project(opensn
    VERSION 0.0.1
    LANGUAGES C CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
include(CheckTypeSize)
include(CheckSymbolExists)
include(CheckCXXCompilerFlag)
include(CheckLinkerFlag)
include(FetchContent)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O3 -g -DNDEBUG" CACHE STRING "" FORCE)
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -DNDEBUG" CACHE STRING "" FORCE)
endif()

set(BuildValues "Release;Debug;RelWithDebInfo;MinSizeRel")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
  message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${BuildValues})

function(_cxx_flag_supported OUTVAR FLAG)
  set(_save_req "${CMAKE_REQUIRED_FLAGS}")
  if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_REQUIRED_FLAGS "${_save_req} -Werror -Werror=unknown-warning-option -Werror=unused-command-line-argument")
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_REQUIRED_FLAGS "${_save_req} -Werror")
  elseif (MSVC)
    set(CMAKE_REQUIRED_FLAGS "${_save_req} /WX")
  endif()

  check_cxx_compiler_flag("${FLAG}" _ok)
  set(CMAKE_REQUIRED_FLAGS "${_save_req}")
  set(${OUTVAR} "${_ok}" PARENT_SCOPE)
endfunction()

function(_ld_flag_supported OUTVAR FLAG)
  set(_save_req_link "${CMAKE_REQUIRED_LINK_OPTIONS}")
  if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_REQUIRED_LINK_OPTIONS "${_save_req_link};-Wl,-fatal_warnings")
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_REQUIRED_LINK_OPTIONS "${_save_req_link}")
  elseif (MSVC)
    set(CMAKE_REQUIRED_LINK_OPTIONS "${_save_req_link};/WX")
  endif()

  check_linker_flag(CXX "${FLAG}" _ok)
  set(CMAKE_REQUIRED_LINK_OPTIONS "${_save_req_link}")
  set(${OUTVAR} "${_ok}" PARENT_SCOPE)
endfunction()

function(_add_cxx_cfg CFG FLAG)
  _cxx_flag_supported(_ok "${FLAG}")
  if (_ok)
    add_compile_options($<$<CONFIG:${CFG}>:${FLAG}>)
  endif()
endfunction()

function(_add_ld_cfg CFG FLAG)
  _ld_flag_supported(_ok "${FLAG}")
  if (_ok)
    add_link_options($<$<CONFIG:${CFG}>:${FLAG}>)
  endif()
endfunction()

_add_cxx_cfg(Debug "-O0")
_add_cxx_cfg(Debug "-g")
_add_ld_cfg (Debug "-g")

foreach(cfg IN ITEMS Release RelWithDebInfo)
  foreach(f IN ITEMS "-march=native" "-mtune=native")
    _add_cxx_cfg(${cfg} "${f}")
  endforeach()

  foreach(f IN ITEMS "-fno-math-errno" "-fno-trapping-math" "-ffp-contract=fast")
    _add_cxx_cfg(${cfg} "${f}")
  endforeach()

  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    foreach(f IN ITEMS "-frename-registers" "-fno-semantic-interposition")
      _add_cxx_cfg(${cfg} "${f}")
    endforeach()
  endif()
endforeach()

option(OPENSN_WITH_CUDA "Enable CUDA support" OFF)
if (OPENSN_WITH_CUDA)
  set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
  find_package(CUDAToolkit 12.0 REQUIRED)
  message(STATUS "CUDA support enabled")
  message(STATUS "CUDA architectures detected: ${CMAKE_CUDA_ARCHITECTURES}")
  enable_language(CUDA)
endif()
option(OPENSN_WITH_PYTHON "Build with python support" ON)
option(OPENSN_WITH_PYTHON_MODULE "Build with python module" OFF)

if(OPENSN_WITH_PYTHON_MODULE)
  set(OPENSN_WITH_PYTHON ON CACHE BOOL "Python support requested for building module" FORCE)
endif()

# dependencies
find_package(MPI REQUIRED)

find_package(HDF5 REQUIRED COMPONENTS C HL)
find_package(Boost CONFIG)
if (NOT Boost_FOUND)
  message(STATUS "Boost CONFIG not found, falling back to module mode")
  find_package(Boost REQUIRED)
endif()

find_package(VTK QUIET)
if(VTK_VERSION VERSION_GREATER_EQUAL "9.0.0")
  find_package(VTK
        COMPONENTS
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
        REQUIRED
    )
  message(STATUS "Found VTK: ${VTK_DIR} (found version \"${VTK_VERSION}\")")
elseif(VTK_VERSION VERSION_GREATER_EQUAL "8.0.0")
  find_package(VTK
        COMPONENTS
            vtkCommonCore
            vtkCommonDataModel
            vtkIOLegacy
            vtkIOCore
            vtkIOXML
            vtkParallelCore
            vtkIOParallelXML
            vtkFiltersCore
            vtkIOEnSight
            vtkIOExodus
        REQUIRED
    )
  message(STATUS "Found VTK: ${VTK_DIR} (found version \"${VTK_VERSION}\")")
else()
  if (${VTK_DIR} STREQUAL VTK_DIR-NOTFOUND)
    message(FATAL_ERROR "Could not find VTK")
  else()
    message(FATAL_ERROR "Found unsupported version of VTK: ${VTK_DIR}, version ${VTK_VERSION}")
  endif()
endif()

find_package(PETSc 3.17 REQUIRED)
check_symbol_exists(PETSC_USE_64BIT_INDICES "${PETSC_INCLUDE_DIR}/petscconf.h" PETSC_USE_64BIT_INDICES)
if(NOT ${PETSC_USE_64BIT_INDICES} MATCHES 1)
  message(FATAL_ERROR "PETSc has not been configured with the flag --with-64-bit-indices\n")
endif()
check_symbol_exists(PETSC_HAVE_PTSCOTCH "${PETSC_INCLUDE_DIR}/petscconf.h" PETSC_HAVE_PTSCOTCH)
if(NOT ${PETSC_HAVE_PTSCOTCH} MATCHES 1)
  message(FATAL_ERROR "PETSc has not been configured with PTSCOTCH\n")
endif()

find_package(caliper REQUIRED)
if(caliper_FOUND)
  message(STATUS "Found Caliper: ${caliper_DIR}")
endif()

# google test
find_package(GTest QUIET)
if(NOT TARGET GTest::gtest_main)
  message(STATUS "System GoogleTest not found. Fetching GoogleTest via FetchContent.")
  FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.15.2
    )
  FetchContent_MakeAvailable(googletest)
endif()

# compile options
set(OPENSN_CXX_FLAGS)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  if (NOT OPENSN_WITH_CUDA)
    list(APPEND OPENSN_CXX_FLAGS "-pedantic")
  endif()
  list(APPEND OPENSN_CXX_FLAGS "-Wall")
  list(APPEND OPENSN_CXX_FLAGS "-Wno-unused-variable")
  list(APPEND OPENSN_CXX_FLAGS "-Wno-c11-extensions")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")
  list(APPEND OPENSN_CXX_FLAGS "-pedantic")
  list(APPEND OPENSN_CXX_FLAGS "-Wall")
  list(APPEND OPENSN_CXX_FLAGS "-Wno-unused-variable")
  list(APPEND OPENSN_CXX_FLAGS "-Wno-c11-extensions")
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
  if (NOT OPENSN_WITH_CUDA)
    list(APPEND OPENSN_CXX_FLAGS "-pedantic")
  endif()
  list(APPEND OPENSN_CXX_FLAGS "-Wall")
  list(APPEND OPENSN_CXX_FLAGS "-Wno-unused-variable")
  list(APPEND OPENSN_CXX_FLAGS "-Wno-sign-compare")
  list(APPEND OPENSN_CXX_FLAGS "-Wno-psabi")
  list(APPEND OPENSN_CXX_FLAGS "-Wno-deprecated-declarations")
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
  list(APPEND OPENSN_CXX_FLAGS "-fp-model=precise")
else()
  message(WARNING "Untested compiler : ${CMAKE_CXX_COMPILER_ID}")
endif()

# libopensn
file(GLOB_RECURSE LIBOPENSN_SRCS CONFIGURE_DEPENDS
    framework/*.cc
    modules/*.cc
)

if (OPENSN_WITH_CUDA)
  file(GLOB_RECURSE LIBOPENSN_CUDA_SRCS CONFIGURE_DEPENDS modules/*.cu)
  add_library(libopensncuda STATIC ${LIBOPENSN_CUDA_SRCS})
endif()

add_library(libopensn SHARED ${LIBOPENSN_SRCS})

if (OPENSN_WITH_CUDA)
  target_compile_options(libopensncuda PRIVATE --expt-extended-lambda --expt-relaxed-constexpr)
  target_compile_features(libopensncuda PUBLIC cuda_std_20)
  target_compile_definitions(libopensncuda PRIVATE __OPENSN_USE_CUDA__)
  target_compile_definitions(libopensn PRIVATE __OPENSN_USE_CUDA__)
  set_property(TARGET libopensncuda PROPERTY POSITION_INDEPENDENT_CODE ON)
  set_property(TARGET libopensncuda PROPERTY OUTPUT_NAME opensncuda)
endif()

target_include_directories(libopensn
    PRIVATE
    $<INSTALL_INTERFACE:include/opensn>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/external
    PUBLIC
    $<BUILD_INTERFACE:${HDF5_INCLUDE_DIRS}>
    $<BUILD_INTERFACE:${caliper_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${PETSC_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${Boost_INCLUDE_DIRS}>
)
if (OPENSN_WITH_CUDA)
  target_include_directories(libopensncuda PUBLIC $<BUILD_INTERFACE:${MPI_CXX_INCLUDE_DIRS}>)
  target_compile_options(libopensncuda PRIVATE $<BUILD_INTERFACE:${OPENSN_CXX_FLAGS}>)
  target_compile_options(libopensncuda PRIVATE $<BUILD_INTERFACE:-diag-suppress 1301>)
  target_compile_options(libopensncuda PRIVATE
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<COMPILE_LANGUAGE:CUDA>>:
    $<BUILD_INTERFACE:-Wno-deprecated-declarations -Wno-non-template-friend -Wno-parentheses>>
  )
  target_compile_options(libopensncuda PRIVATE
    $<$<AND:$<CXX_COMPILER_ID:Clang>,$<COMPILE_LANGUAGE:CUDA>>:
    $<BUILD_INTERFACE:-Wno-unused-function -Wno-constant-logical-operand>>
  )
  target_include_directories(libopensncuda
      PRIVATE
      $<INSTALL_INTERFACE:include/opensn>
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
      $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${PROJECT_SOURCE_DIR}
      ${PROJECT_SOURCE_DIR}/external
      PUBLIC
      $<BUILD_INTERFACE:${HDF5_INCLUDE_DIRS}>
      $<BUILD_INTERFACE:${caliper_INCLUDE_DIR}>
      $<BUILD_INTERFACE:${PETSC_INCLUDE_DIR}>
      $<BUILD_INTERFACE:${Boost_INCLUDE_DIRS}>
  )
  target_link_libraries(libopensncuda PRIVATE MPI::MPI_CXX)
endif()

target_link_libraries(libopensn
    PRIVATE
    ${PETSC_LIBRARY}
    ${VTK_LIBRARIES}
    caliper
    ${HDF5_LIBRARIES}
    MPI::MPI_CXX
)
if (OPENSN_WITH_CUDA)
  target_link_libraries(libopensncuda PRIVATE ${CUDA_LIBRARIES} CUDA::cublas)
  target_link_libraries(libopensn INTERFACE libopensncuda)
endif()

target_compile_options(libopensn PRIVATE ${OPENSN_CXX_FLAGS})

set_target_properties(
    libopensn
    PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
        OUTPUT_NAME opensn
)

if(OPENSN_WITH_PYTHON)
  add_subdirectory(python)
  add_subdirectory(test)

  if(OPENSN_WITH_PYTHON_MODULE)
    add_subdirectory(pyopensn)
  endif()

  add_custom_target(test
      COMMAND
      ${CMAKE_BINARY_DIR}/test/opensn-unit &&
              ${CMAKE_SOURCE_DIR}/test/run_tests
              -d ${CMAKE_SOURCE_DIR}/test/python
              --exe ${CMAKE_BINARY_DIR}/test/python/opensn-test
              -w3
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
endif()

configure_file(config.h.in config.h)

# install
configure_package_config_file(
    ${PROJECT_SOURCE_DIR}/cmake/OpenSnConfig.cmake.in
    OpenSnConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/opensn
    NO_SET_AND_CHECK_MACRO
)

write_basic_package_version_file(
    OpenSnConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(
    TARGETS libopensn
    EXPORT openSnTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

if(OPENSN_WITH_CUDA)
  install(
        TARGETS libopensncuda
        EXPORT openSnTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )
endif()

install(
    DIRECTORY ${CMAKE_SOURCE_DIR}/framework
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/opensn
    FILES_MATCHING PATTERN "*.h"
)

install(
    DIRECTORY ${CMAKE_SOURCE_DIR}/modules
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/opensn
    FILES_MATCHING PATTERN "*.h"
)

install(
    DIRECTORY ${CMAKE_SOURCE_DIR}/external/mpicpp-lite
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/opensn
    FILES_MATCHING PATTERN "*.h"
)

install(
    EXPORT openSnTargets
    FILE OpenSnTargets.cmake
    NAMESPACE opensn::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/opensn
)

install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/OpenSnConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/OpenSnConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/opensn
)

if(OPENSN_WITH_PYTHON)
  install(
        TARGETS opensn libopensnpy
        EXPORT openSnTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
    )

  install(
        DIRECTORY ${CMAKE_SOURCE_DIR}/python
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/opensn
        FILES_MATCHING PATTERN "*.h"
    )
endif()
