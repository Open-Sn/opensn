file (GLOB_RECURSE MORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")

set(SOURCES ${SOURCES} ${MORE_SOURCES} PARENT_SCOPE)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/Groupset")

#set_source_files_properties(chi_console_00_constrdestr.cc PROPERTIES COMPILE_FLAGS -Wno-effc++)