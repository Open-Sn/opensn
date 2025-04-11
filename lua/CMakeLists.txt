# libopensnlua
file(GLOB_RECURSE LIBOPENSN_LUA_SRCS CONFIGURE_DEPENDS lib/*.cc)
add_library(libopensnlua SHARED ${LIBOPENSN_LUA_SRCS})

target_include_directories(libopensnlua
    PRIVATE
    $<INSTALL_INTERFACE:include/opensn>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/external
    PUBLIC
    ${LUA_INCLUDE_DIR}
)

target_link_libraries(libopensnlua
    PRIVATE
    libopensn
    ${PETSC_LIBRARY}
    MPI::MPI_CXX
    caliper
    PUBLIC
    ${LUA_LIBRARIES}
)

target_compile_definitions(libopensnlua PRIVATE OPENSN_WITH_LUA)

target_compile_options(libopensnlua PRIVATE ${OPENSN_CXX_FLAGS})

if(NOT MSVC)
    set_target_properties(libopensnlua PROPERTIES OUTPUT_NAME opensnlua)
endif()

set_target_properties(
    libopensnlua
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
    ${LUA_LIBRARIES}
    libopensn
    libopensnlua
    MPI::MPI_CXX
)

target_compile_definitions(opensn PRIVATE OPENSN_WITH_LUA)

target_compile_options(opensn PRIVATE ${OPENSN_CXX_FLAGS})
