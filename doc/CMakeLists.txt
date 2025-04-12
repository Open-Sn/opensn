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

    macro(generate_tutorial srcDir)
        message(STATUS "Generating tutorial ${srcDir}")
        file(RELATIVE_PATH tutPath ${CMAKE_SOURCE_DIR}/tutorials/python ${srcDir})
        set(destDir ${CMAKE_BINARY_DIR}/doc/src/tutorials)
        # copy images
        file(MAKE_DIRECTORY ${destDir}/${tutPath})
        copy_files("png" ${srcDir} ${destDir}/${tutPath})
        copy_files("jpg" ${srcDir} ${destDir}/${tutPath})
        copy_files("gif" ${srcDir} ${destDir}/${tutPath})
    endmacro()

    configure_file(Doxyfile.in Doxyfile)
    configure_file(conf.py.in conf.py)

    add_custom_target(doc DEPENDS ${PROJECT_BINARY_DIR}/html/index.html)

    configure_file(index.rst src/index.rst)
    # installation guide
    configure_file(install/index.rst src/install/index.rst COPYONLY)
    configure_file(install/install_linux.rst src/install/install_linux.rst COPYONLY)
    configure_file(install/install_linux_easy.rst src/install/install_linux_easy.rst COPYONLY)
    configure_file(install/install_macos.rst src/install/install_macos.rst COPYONLY)
    configure_file(install/install_macos_easy.rst src/install/install_macos_easy.rst COPYONLY)
    # theory guide
    configure_file(theory/index.rst src/theory/index.rst)
    configure_file(theory/background.rst src/theory/background.rst)
    configure_file(theory/discretization.rst src/theory/discretization.rst)
    configure_file(theory/cross_sections.rst src/theory/cross_sections.rst)
    configure_file(theory/outcome.rst src/theory/outcome.rst)
    configure_file(theory/iterative.rst src/theory/iterative.rst)
    configure_file(theory/parallelization.rst src/theory/parallelization.rst)
    configure_file(theory/adjoint.rst src/theory/adjoint.rst)
    # developers guide
    configure_file(devguide/index.rst src/devguide/index.rst)
    configure_file(devguide/coding_standard.md src/devguide/coding_standard.md)
    configure_file(devguide/workflow.md src/devguide/workflow.md)
    configure_file(devguide/lua_bindings.md src/devguide/lua_bindings.md)
    configure_file(devguide/py_bindings.rst src/devguide/py_bindings.rst COPYONLY)
    # Python API
    configure_file(pyapi/index.rst src/pyapi/index.rst COPYONLY)
    file(MAKE_DIRECTORY src/pyapi/generated)
    configure_file(_static/pyopensn.css src/_static/pyopensn.css COPYONLY)
    configure_file(_templates/function.rst src/_templates/function.rst COPYONLY)
    configure_file(_templates/noinit.rst src/_templates/noinit.rst COPYONLY)
    configure_file(_templates/python.rst src/_templates/python.rst COPYONLY)
    # Tutorials TOC
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/python/index.rst src/tutorials/index.rst COPYONLY)
    # Tutorials LBS
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/python/lbs/index.rst src/tutorials/lbs/index.rst COPYONLY)
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/python/lbs/first/first_example.ipynb
                   src/tutorials/lbs/first/first_example.ipynb COPYONLY)
    # Tutorials mesh
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/python/mesh/index.rst
                   src/tutorials/mesh/index.rst COPYONLY)
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/python/mesh/introduction.rst
                   src/tutorials/mesh/introduction.rst COPYONLY)
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/python/mesh/mesh_ortho_2D.ipynb
                   src/tutorials/mesh/mesh_ortho_2D.ipynb COPYONLY)
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/python/mesh/mesh_ortho_3D.ipynb
            src/tutorials/mesh/mesh_ortho_3D.ipynb COPYONLY)
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/python/mesh/mesh_read_gmsh_2D.ipynb
            src/tutorials/mesh/mesh_read_gmsh_2D.ipynb COPYONLY)
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/python/mesh/mesh_read_gmsh_3D.ipynb
            src/tutorials/mesh/mesh_read_gmsh_3D.ipynb COPYONLY)
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/python/mesh/mesh_read_obj_2D.ipynb
            src/tutorials/mesh/mesh_read_obj_2D.ipynb COPYONLY)
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/python/mesh/mesh_read_vtu_3D.ipynb
            src/tutorials/mesh/mesh_read_vtu_3D.ipynb COPYONLY)
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/python/mesh/mesh_extrude.ipynb
            src/tutorials/mesh/mesh_extrude.ipynb COPYONLY)
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/python/mesh/mesh_ortho_3D_blockIDs.ipynb
            src/tutorials/mesh/mesh_ortho_3D_blockIDs.ipynb COPYONLY)
    # Tutorials XS
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/python/xs/index.rst
            src/tutorials/xs/index.rst COPYONLY)
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/python/xs/introduction.rst
            src/tutorials/xs/introduction.rst COPYONLY)
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/python/xs/1g_xs_simple.ipynb
            src/tutorials/xs/1g_xs_simple.ipynb COPYONLY)
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/python/xs/1g_xs_read_file.ipynb
            src/tutorials/xs/1g_xs_read_file.ipynb COPYONLY)
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/python/xs/mg_xs_read_file.ipynb
            src/tutorials/xs/mg_xs_read_file.ipynb COPYONLY)
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/python/xs/mg_xs_read_openmc.ipynb
            src/tutorials/xs/mg_xs_read_openmc.ipynb COPYONLY)
    # Tutorials aquad
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/python/aquad/index.rst
            src/tutorials/aquad/index.rst COPYONLY)
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/python/aquad/introduction.rst
            src/tutorials/aquad/introduction.rst COPYONLY)
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/python/aquad/product_GLC.ipynb
            src/tutorials/aquad/product_GLC.ipynb COPYONLY)
    configure_file(${CMAKE_SOURCE_DIR}/tutorials/python/aquad/sLDFE.ipynb
            src/tutorials/aquad/sLDFE.ipynb COPYONLY)

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

    generate_tutorial(${CMAKE_SOURCE_DIR}/tutorials/python/lbs/first)
    generate_tutorial(${CMAKE_SOURCE_DIR}/tutorials/python/mesh)
    generate_tutorial(${CMAKE_SOURCE_DIR}/tutorials/python/aquad)

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
