find_package(Python3 REQUIRED COMPONENTS Interpreter)
find_package(Doxygen QUIET)
find_program(SPHINX_BUILD sphinx-build)
mark_as_advanced(FORCE SPHINX_BUILD)

if(DOXYGEN_FOUND AND NOT ${SPHINX_BUILD} STREQUAL SPHINX_BUILD-NOTFOUND)
    macro(copy_files ext srcDir destDir)
        message(STATUS "Copying .${ext} files to ${destDir}")
        file(MAKE_DIRECTORY ${destDir})

        file(GLOB_RECURSE templateFiles RELATIVE ${srcDir} "${srcDir}/*.${ext}")
        foreach(templateFile ${templateFiles})
            set(srcTemplatePath ${srcDir}/${templateFile})
            if(NOT IS_DIRECTORY ${srcTemplatePath})
                configure_file(${srcTemplatePath} ${destDir}/${templateFile} COPYONLY)
            endif()
        endforeach()
    endmacro()

    macro(lua2md src dest)
        add_custom_command(
            OUTPUT ${dest}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/doc/scripts/lua2md.py ${src} ${dest}
            DEPENDS ${src} ${CMAKE_SOURCE_DIR}/doc/scripts/lua2md.py
        )
        list(APPEND RST_FILES ${dest})
    endmacro()

    configure_file(Doxyfile.in Doxyfile)
    configure_file(conf.py.in conf.py)

    add_custom_target(doc DEPENDS ${PROJECT_BINARY_DIR}/html/index.html)

    configure_file(index.rst src/index.rst)
    configure_file(tutorials/old/index.rst src/tutorials/old/index.rst)
    configure_file(tutorials/lbs/index.rst src/tutorials/lbs/index.rst)
    configure_file(tutorials/lbs/meshing/index.rst src/tutorials/lbs/meshing/index.rst)
    configure_file(tutorials/lbs/material/index.rst src/tutorials/lbs/material/index.rst)
    copy_files("md" ${CMAKE_SOURCE_DIR}/doc ${CMAKE_BINARY_DIR}/doc/src)
    copy_files("png" ${CMAKE_SOURCE_DIR}/doc ${CMAKE_BINARY_DIR}/doc/src)
    copy_files("gif" ${CMAKE_SOURCE_DIR}/doc ${CMAKE_BINARY_DIR}/doc/src)

    file(GLOB_RECURSE CONFIGURE_DEPENDS RST_FILES
        ${PROJECT_BINARY_DIR}/doc/*.rst
        ${PROJECT_BINARY_DIR}/doc/*.md
    )
    file(GLOB_RECURSE HEADER_FILES
        ${PROJECT_SOURCE_DIR}/framework/*.h
        ${PROJECT_SOURCE_DIR}/modules/*.h
    )

    lua2md(${CMAKE_SOURCE_DIR}/doc/tutorials/lbs/first/first_example.lua
            ${CMAKE_BINARY_DIR}/doc/src/tutorials/lbs/first/first_example.md)
    lua2md(${CMAKE_SOURCE_DIR}/doc/tutorials/lbs/material/mg_xs.lua
            ${CMAKE_BINARY_DIR}/doc/src/tutorials/lbs/material/mg_xs.md)
    lua2md(${CMAKE_SOURCE_DIR}/doc/tutorials/lbs/meshing/ortho_2D.lua
            ${CMAKE_BINARY_DIR}/doc/src/tutorials/lbs/meshing/ortho_2D.md)
    lua2md(${CMAKE_SOURCE_DIR}/doc/tutorials/lbs/meshing/ortho_3D.lua
            ${CMAKE_BINARY_DIR}/doc/src/tutorials/lbs/meshing/ortho_3D.md)
    lua2md(${CMAKE_SOURCE_DIR}/doc/tutorials/lbs/meshing/read_2D_obj.lua
            ${CMAKE_BINARY_DIR}/doc/src/tutorials/lbs/meshing/read_2D_obj.md)
    lua2md(${CMAKE_SOURCE_DIR}/doc/tutorials/lbs/meshing/transport_simulation_part_2D.lua
            ${CMAKE_BINARY_DIR}/doc/src/tutorials/lbs/meshing/transport_simulation_part_2D.md)
    lua2md(${CMAKE_SOURCE_DIR}/doc/tutorials/lbs/meshing/extrusion.lua
            ${CMAKE_BINARY_DIR}/doc/src/tutorials/lbs/meshing/extrusion.md)
    lua2md(${CMAKE_SOURCE_DIR}/doc/tutorials/lbs/meshing/transport_simulation_part_3D.lua
            ${CMAKE_BINARY_DIR}/doc/src/tutorials/lbs/meshing/transport_simulation_part_3D.md)
    lua2md(${CMAKE_SOURCE_DIR}/doc/tutorials/lbs/meshing/read_3D_vtu.lua
            ${CMAKE_BINARY_DIR}/doc/src/tutorials/lbs/meshing/read_3D_vtu.md)
    lua2md(${CMAKE_SOURCE_DIR}/doc/tutorials/lbs/meshing/transport_simulation_part_3D_single_agg.lua
            ${CMAKE_BINARY_DIR}/doc/src/tutorials/lbs/meshing/transport_simulation_part_3D_single_agg.md)

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
