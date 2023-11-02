file (GLOB_RECURSE MORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")

FOREACH(item ${MORE_SOURCES})
    IF(${item} MATCHES "lua/diffusion_maindoc_page.cc")
        LIST(REMOVE_ITEM MORE_SOURCES ${item})
    ENDIF(${item} MATCHES "lua/diffusion_maindoc_page.cc")
ENDFOREACH(item)

set(SOURCES ${SOURCES} ${MORE_SOURCES} PARENT_SCOPE)

#set_source_files_properties(chi_console_00_constrdestr.cc PROPERTIES COMPILE_FLAGS -Wno-effc++)