#pragma once

#include "framework/chi_lua.h"

/**Gets a field-function handle by name.
 * \param FFname string Name of the field function.
 *
 * \return handle If the field-function was found and a handle identified the valid
 *                handle will be returned (i.e., a natural number >= 0). If the
 *                field-function by the given name was not found then the function
 *                will return null.
 *
 * \ingroup LuaFieldFunc
 */
int chiGetFieldFunctionHandleByName(lua_State* L);

/** Exports a field function to VTK format.
 *
 * \param FFHandle int Global handle to the field function.
 * \param BaseName char Base name for the exported file.
 *
 * \ingroup LuaFieldFunc
 * \author Jan
 */
int chiExportFieldFunctionToVTK(lua_State* L);

/** Exports all the field functions in a list to VTK format.
 *
 * \param listFFHandles table Global handles or names to the field functions
 * \param BaseName char Base name for the exported file.
 *
 * \ingroup LuaFieldFunc
 * \author Jan
 */
int chiExportMultiFieldFunctionToVTK(lua_State* L);
