#pragma once

#include "chi_lua.h"

namespace prk::lua_utils
{
int chiPRKGetParam(lua_State* L);
int chiPRKSetParam(lua_State* L);

chi::InputParameters GetSyntax_SetParam();
chi::ParameterBlock SetParam(const chi::InputParameters& params);

chi::InputParameters GetParamSyntax();
chi::ParameterBlock GetParam(const chi::InputParameters& params);
} // namespace prk::lua_utils

