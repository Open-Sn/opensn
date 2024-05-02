# test binary
file(GLOB_RECURSE TEST_SRCS CONFIGURE_DEPENDS *.cc)

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
)

target_compile_definitions(opensn-test PRIVATE OPENSN_WITH_LUA)

target_compile_options(opensn-test PRIVATE ${OPENSN_CXX_FLAGS})
