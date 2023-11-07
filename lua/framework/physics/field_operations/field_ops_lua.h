#pragma once

#include "framework/lua.h"

namespace chi_physics::field_operations::lua_utils
{

/**Executes a field function operation.
 * \param handle int Handle to the field function operation object.*/
int chiFieldOperationExecute(lua_State* L);

} // namespace chi_physics::field_operations::lua_utils
