#ifndef CHITECH_LUA_MESHMODIFIERS_H
#define CHITECH_LUA_MESHMODIFIERS_H

#include "parameters/input_parameters.h"

namespace chi_mesh::lua_utils
{
chi::InputParameters MeshModifiersApply_Syntax();
chi::ParameterBlock MeshModifiersApply(const chi::InputParameters& params);
} // namespace chi_mesh::lua_utils

#endif // CHITECH_LUA_MESHMODIFIERS_H
