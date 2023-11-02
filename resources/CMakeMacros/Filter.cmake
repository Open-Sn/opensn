if(__FILTER_INCLUDED)
    return()
endif()
set(__FILTER_INCLUDED TRUE)

set (FILTER_CMAKE_MODULE_VERSION "0.2")

#function (filterItems varFilter regFilter)
#    #----------------------------------------------------- Checks each item in the list
#    foreach (fItem ${${varFilter}})
#        #------------------------------------------------- Checks for matches
#        if ("${fItem}" MATCHES ${regFilter})
#            #--------------------------------------------- Removes the match
#            list (REMOVE_ITEM ${varFilter} ${fItem})
#        endif ("${fItem}" MATCHES ${regFilter})
#    endforeach(fItem)
#    #----------------------------------------------------- Delivers refunded value
#    set(${varFilter} ${${varFilter}} PARENT_SCOPE)
#endfunction (filterItems)

FUNCTION (filterDirectories _InFileList _excludeDirName _verbose)
    foreach (ITR ${_InFileList})
        if ("${_verbose}")
            message(STATUS "ITR=${ITR}")
        endif ("${_verbose}")
        #-------------------------------------------------------------------------------- Check if the item matches the directory name in _excludeDirName
        if ("${ITR}" MATCHES "(.*)${_excludeDirName}(.*)")
            #---------------------------------------------------------------------------- Remove the item from the list
            #message(STATUS "Remove Item from List:${ITR}")
            list (REMOVE_ITEM _InFileList ${ITR})
        endif ("${ITR}" MATCHES "(.*)${_excludeDirName}(.*)")
    endforeach(ITR)
    #------------------------------------------------------------------------------------ Return the SOURCE_FILES variable to the calling parent
    set(SOURCE_FILES ${_InFileList} PARENT_SCOPE)
ENDFUNCTION (filterDirectories)

message (STATUS "Filter loaded.")
