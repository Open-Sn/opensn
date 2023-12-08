find_package(Doxygen QUIET)
find_program(SPHINX_BUILD sphinx-build)
mark_as_advanced(FORCE SPHINX_BUILD)

if(DOXYGEN_FOUND AND NOT ${SPHINX_BUILD} STREQUAL SPHINX_BUILD-NOTFOUND)
    macro(lua2md src dest)
        add_custom_command(
            OUTPUT ${dest}
            COMMAND python ${CMAKE_SOURCE_DIR}/doc/scripts/lua2md.py ${src} ${dest}
            DEPENDS ${src} ${CMAKE_SOURCE_DIR}/doc/scripts/lua2md.py
        )
        list(APPEND RST_FILES ${dest})
    endmacro()

    configure_file(Doxyfile.in Doxyfile)
    configure_file(conf.py.in conf.py)

    add_custom_target(doc DEPENDS ${PROJECT_BINARY_DIR}/html/index.html)

    configure_file(index.rst src/index.rst)

    file(GLOB_RECURSE CONFIGURE_DEPENDS RST_FILES
        ${PROJECT_BINARY_DIR}/doc/*.rst
        ${PROJECT_BINARY_DIR}/doc/*.md
    )
    file(GLOB_RECURSE HEADER_FILES
        ${PROJECT_SOURCE_DIR}/framework/*.h
        ${PROJECT_SOURCE_DIR}/modules/*.h
    )

    add_custom_command(
        OUTPUT
            ${PROJECT_BINARY_DIR}/html/index.html
        COMMAND
            ${SPHINX_BUILD} -c ${PROJECT_BINARY_DIR}/doc -b html -d ${PROJECT_BINARY_DIR}/doc/.doctrees ${PROJECT_BINARY_DIR}/doc/src html
        DEPENDS
            ${PROJECT_BINARY_DIR}/doc/conf.py
            ${PROJECT_BINARY_DIR}/doc/xml/index.xml
            ${RST_FILES}
    )

    add_custom_command(
        OUTPUT
            ${PROJECT_BINARY_DIR}/doc/xml/index.xml
        COMMAND
            doxygen Doxyfile
        DEPENDS
            ${PROJECT_BINARY_DIR}/doc/Doxyfile
            ${HEADER_FILES}
    )

    add_custom_command(
        TARGET doc
        POST_BUILD
        COMMAND echo "Open ${PROJECT_BINARY_DIR}/doc/html/index.html in your browser to view the documentation."
    )

else()
    add_custom_target(doc)

    add_custom_command(
        TARGET doc
        COMMAND echo "Unable to generate documentation:"
    )
    if (NOT DOXYGEN_FOUND)
        add_custom_command(
            TARGET doc
            COMMAND echo "- 'doxygen' is not installed."
        )
    endif()

    if (${SPHINX_BUILD} STREQUAL SPHINX_BUILD-NOTFOUND)
        add_custom_command(
            TARGET doc
            COMMAND echo "- 'sphinx-build' is not installed."
        )
    endif()
endif()
