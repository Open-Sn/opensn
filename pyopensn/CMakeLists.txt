# Find Python3
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# Get pybind11 CMake directory and append it to CMake prefix path
execute_process(
    COMMAND ${Python3_EXECUTABLE} -m pybind11 --cmakedir
    OUTPUT_VARIABLE PYBIND11_CMAKE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
list(APPEND CMAKE_PREFIX_PATH "${PYBIND11_CMAKE_DIR}")

# Execute Python to get the mpi4py include directory.
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import mpi4py; print(mpi4py.get_include())"
    OUTPUT_VARIABLE MPI4PY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Locate Pybind11
find_package(pybind11 CONFIG)

# Construct module
file(GLOB_RECURSE PY_INTERFACE_FILE_LIST ${CMAKE_CURRENT_SOURCE_DIR}/*.cc)
pybind11_add_module(__init__ NO_EXTRAS ${PY_INTERFACE_FILE_LIST})

# Module properties
target_include_directories(__init__
    PRIVATE
    $<INSTALL_INTERFACE:include/opensn>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/external
    ${MPI4PY_INCLUDE_DIR}
)
target_link_libraries(__init__
    PRIVATE
    libopensn
    libopensnpy
    ${PETSC_LIBRARY}
    MPI::MPI_CXX
    caliper
)
