# libopensnpython
file(GLOB_RECURSE LIBOPENSN_PY_SRCS CONFIGURE_DEPENDS lib/*.cc)
add_library(libopensnpy SHARED ${LIBOPENSN_PY_SRCS})

target_include_directories(libopensnpy
    PRIVATE
    $<INSTALL_INTERFACE:include/opensn>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/external
    pybind11::embed
)

target_link_libraries(libopensnpy
    PRIVATE
    libopensn
    ${PETSC_LIBRARY}
    MPI::MPI_CXX
    caliper
)

target_compile_definitions(libopensnpy PRIVATE OPENSN_WITH_PYTHON)

target_compile_options(libopensnpy PRIVATE ${OPENSN_CXX_FLAGS})

if(NOT MSVC)
    set_target_properties(libopensnpy PROPERTIES OUTPUT_NAME openspy)
endif()

set_target_properties(
    libopensnpy
    PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
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
    libopensn
    libopensnpy
    MPI::MPI_CXX
)

target_compile_definitions(opensn PRIVATE OPENSN_WITH_PYTHON)

target_compile_options(opensn PRIVATE ${OPENSN_CXX_FLAGS})
