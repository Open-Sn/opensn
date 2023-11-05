# test binary
file(GLOB_RECURSE TEST_SRCS *.cc)

add_executable(opensn-test ${TEST_SRCS})

target_include_directories(opensn-test
    PRIVATE
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    ${PROJECT_SOURCE_DIR}
)

target_link_libraries(opensn-test
    PRIVATE
    libopensn
    libopensnlua
    ${LUA_LIBRARIES}
    ${PETSC_LIBRARY}
    MPI::MPI_CXX
)

target_compile_definitions(opensn-test PRIVATE OPENSN_WITH_LUA)
