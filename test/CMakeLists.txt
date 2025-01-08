# test binary
file(GLOB_RECURSE TEST_SRCS CONFIGURE_DEPENDS framework/*.cc modules/*.cc lua/*.cc src/*.cc)

add_executable(opensn-test ${TEST_SRCS})

target_include_directories(opensn-test
    PRIVATE
    $<INSTALL_INTERFACE:include/opensn>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/external
)

target_link_libraries(opensn-test
    PRIVATE
    libopensn
    libopensnlua
    ${LUA_LIBRARIES}
    ${PETSC_LIBRARY}
    ${HDF5_LIBRARIES}
    caliper
    MPI::MPI_CXX
)

target_compile_definitions(opensn-test PRIVATE OPENSN_WITH_LUA)

target_compile_options(opensn-test PRIVATE ${OPENSN_CXX_FLAGS})

# unit test binary
file(GLOB_RECURSE UNIT_TEST_SRCS CONFIGURE_DEPENDS unit/*.cc)
add_executable(opensn-unit ${UNIT_TEST_SRCS})

include(GoogleTest)

target_include_directories(opensn-unit
    PRIVATE
    $<INSTALL_INTERFACE:include/opensn>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/external
)

target_link_libraries(opensn-unit
    PRIVATE
    libopensn
    libopensnlua
    ${LUA_LIBRARIES}
    ${PETSC_LIBRARY}
    ${HDF5_LIBRARIES}
    caliper
    MPI::MPI_CXX
    gmock_main
    gtest_main
)

target_compile_definitions(opensn-unit PRIVATE OPENSN_WITH_LUA)

target_compile_options(opensn-unit PRIVATE ${OPENSN_CXX_FLAGS})
