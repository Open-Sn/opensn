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
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/doc/scripts/lua2md.py ${src} ${dest} -d ${CMAKE_SOURCE_DIR}
            DEPENDS ${src} ${CMAKE_SOURCE_DIR}/doc/scripts/lua2md.py
        )
        list(APPEND RST_FILES ${dest})
    endmacro()

    macro(generate_tutorial srcDir)
        message(STATUS "Generating tutorial ${srcDir}")
        file(RELATIVE_PATH tutPath ${CMAKE_SOURCE_DIR}/tutorials ${srcDir})
        set(destDir ${CMAKE_BINARY_DIR}/doc/src/tutorials)
        # copy images
        file(MAKE_DIRECTORY ${destDir}/${tutPath})
        copy_files("png" ${srcDir} ${destDir}/${tutPath})
        copy_files("jpg" ${srcDir} ${destDir}/${tutPath})
        copy_files("gif" ${srcDir} ${destDir}/${tutPath})
        # process lua files
        file(GLOB_RECURSE luaFiles RELATIVE ${srcDir} "${srcDir}/*.lua")
        foreach(luaFile ${luaFiles})
            set(srcLuaPath ${srcDir}/${luaFile})
            if(NOT IS_DIRECTORY ${srcLuaPath})
                get_filename_component(baseName ${luaFile} NAME_WE)
                lua2md(${srcLuaPath} ${destDir}/${tutPath}/${baseName}.md)
            endif()
        endforeach()
    endmacro()

    configure_file(Doxyfile.in Doxyfile)
    configure_file(conf.py.in conf.py)

    add_custom_target(doc DEPENDS ${PROJECT_BINARY_DIR}/html/index.html)

    configure_file(index.rst src/index.rst)
    configure_file(install/index.rst src/install/index.rst COPYONLY)
    configure_file(install/install_linux.rst src/install/install_linux.rst COPYONLY)
    configure_file(install/install_linux_easy.rst src/install/install_linux_easy.rst COPYONLY)
    configure_file(install/install_macos.rst src/install/install_macos.rst COPYONLY)
    configure_file(install/install_macos_easy.rst src/install/install_macos_easy.rst COPYONLY)
    configure_file(theory/index.rst src/theory/index.rst)
    configure_file(theory/background.rst src/theory/background.rst)
    configure_file(theory/discretization.rst src/theory/discretization.rst)
    configure_file(theory/cross_sections.rst src/theory/cross_sections.rst)
    configure_file(theory/outcome.rst src/theory/outcome.rst)
    configure_file(theory/iterative.rst src/theory/iterative.rst)
    configure_file(theory/parallelization.rst src/theory/parallelization.rst)
    configure_file(theory/adjoint.rst src/theory/adjoint.rst)
    configure_file(devguide/index.rst src/devguide/index.rst)
    configure_file(devguide/coding_standard.md src/devguide/coding_standard.md)
    configure_file(devguide/workflow.md src/devguide/workflow.md)
    configure_file(devguide/lua_bindings.md src/devguide/lua_bindings.md)
    configure_file(devguide/py_bindings.rst src/devguide/py_bindings.rst COPYONLY)
    configure_file(pyapi/index.rst src/pyapi/index.rst COPYONLY)
    file(MAKE_DIRECTORY src/pyapi/generated)
    configure_file(_static/pyopensn.css src/_static/pyopensn.css COPYONLY)
    configure_file(_templates/function.rst src/_templates/function.rst COPYONLY)
    configure_file(_templates/noinit.rst src/_templates/noinit.rst COPYONLY)
    configure_file(_templates/python.rst src/_templates/python.rst COPYONLY)
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/lua/index.rst src/tutorials/index.rst)
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/lua/lbs/index.rst src/tutorials/lbs/index.rst)
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/lua/meshing/index.rst src/tutorials/meshing/index.rst)
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/lua/material/index.rst src/tutorials/material/index.rst)

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

    generate_tutorial(${CMAKE_SOURCE_DIR}/tutorials/lua/meshing)
    copy_files("md" ${CMAKE_SOURCE_DIR}/tutorials/lua/meshing ${CMAKE_BINARY_DIR}/doc/src/tutorials/meshing)
    generate_tutorial(${CMAKE_SOURCE_DIR}/tutorials/lua/material)
    generate_tutorial(${CMAKE_SOURCE_DIR}/tutorials/lua/lbs/first)

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
