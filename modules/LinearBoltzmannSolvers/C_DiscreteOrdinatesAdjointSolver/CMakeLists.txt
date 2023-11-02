# The command below recursively looks for *.cc files and
# adds them to the list "MORE_SOURCES".
file (GLOB_RECURSE MORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")

# Example of how to remove certain source files from
# the "MORE_SOURCES" list. These kind of file are usually those
# that form part of doxygen documentation.
#
#    FOREACH(item ${MORE_SOURCES})
#        IF(${item} MATCHES "lua/mesh_maindoc_page.cc")
#            LIST(REMOVE_ITEM MORE_SOURCES ${item})
#        ENDIF(${item} MATCHES "lua/mesh_maindoc_page.cc")
#    ENDFOREACH(item)

# The command below adds the list "MORE_SOURCES" to
# existing "SOURCES" list.
set(SOURCES ${SOURCES} ${MORE_SOURCES} PARENT_SCOPE)

# Example of how to set custom compile options_ for specific source files.
#
# set_source_files_properties(chi_console_00_constrdestr.cc PROPERTIES COMPILE_FLAGS -Wno-effc++)