# Find Python3
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# Get pybind11 CMake directory and append it to CMake prefix path
execute_process(
    COMMAND ${Python3_EXECUTABLE} -m pybind11 --cmakedir
    OUTPUT_VARIABLE PYBIND11_CMAKE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
list(APPEND CMAKE_PREFIX_PATH "${PYBIND11_CMAKE_DIR}")

# Locate Pybind11
find_package(pybind11 CONFIG REQUIRED)

# libopensnpython
file(GLOB_RECURSE LIBOPENSN_PY_SRCS CONFIGURE_DEPENDS lib/*.cc)
add_library(libopensnpy STATIC ${LIBOPENSN_PY_SRCS})
target_include_directories(libopensnpy
    PRIVATE
    $<INSTALL_INTERFACE:include/opensn>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/external
)
target_link_libraries(libopensnpy
    PUBLIC
    libopensn
    ${PETSC_LIBRARY}
    MPI::MPI_CXX
    caliper
    pybind11::embed
)
target_compile_options(libopensnpy PRIVATE ${OPENSN_CXX_FLAGS})
if(NOT MSVC)
    set_target_properties(libopensnpy PROPERTIES OUTPUT_NAME opensnpy)
endif()
set_target_properties(
    libopensnpy
    PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
        POSITION_INDEPENDENT_CODE ON
        INTERPROCEDURAL_OPTIMIZATION OFF
)

# opensn binary
add_executable(opensn "main.cc")
target_include_directories(opensn
    PRIVATE
    $<INSTALL_INTERFACE:include/opensn>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
    ${PROJECT_SOURCE_DIR}/external
)
target_link_libraries(opensn
    PRIVATE
    libopensnpy
)
target_compile_options(opensn PRIVATE ${OPENSN_CXX_FLAGS})
