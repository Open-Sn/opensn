add_subdirectory("paraview")

set(TARGET external)
add_library(${TARGET} STATIC
    $<TARGET_OBJECTS:ThirdPartyParaView>)
