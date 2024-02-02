#pragma once

#include "framework/lua.h"

namespace opensnlua
{

/**Executes a field function operation.
 * \param handle int Handle to the field function operation object.*/
int FieldOperationExecute(lua_State* L);

} // namespace opensnlua
