# test binary
file(GLOB_RECURSE TEST_SRCS *.cc)

add_executable(opensn-test ${TEST_SRCS})
target_include_directories(opensn-test
    PRIVATE
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
)
target_link_libraries(opensn-test
    PRIVATE
    libopensn
    libopensnlua
    ${LUA_LIBRARIES}
    ${PETSC_LIBRARY}
    MPI::MPI_CXX
)
