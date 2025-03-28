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
    ${PETSC_LIBRARY}
    ${HDF5_LIBRARIES}
    caliper
    MPI::MPI_CXX
    GTest::gmock_main
    GTest::gtest_main
)
if(OPENSN_WITH_LUA)
    add_subdirectory(lua)
    target_link_libraries(opensn-unit PRIVATE libopensnlua)
    target_link_libraries(opensn-unit PRIVATE ${LUA_LIBRARIES})
    target_compile_definitions(opensn-unit PRIVATE OPENSN_WITH_LUA)
endif()
if(OPENSN_WITH_PYTHON)
    target_link_libraries(opensn-unit PRIVATE libopensnpy)
    target_compile_definitions(opensn-unit PRIVATE OPENSN_WITH_PYTHON)
endif()

target_compile_options(opensn-unit PRIVATE ${OPENSN_CXX_FLAGS})
