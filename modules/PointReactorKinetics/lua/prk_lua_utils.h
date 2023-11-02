#ifndef CHITECH_POINT_REACTOR_KINETICS_LUA_UTILS_H
#define CHITECH_POINT_REACTOR_KINETICS_LUA_UTILS_H

#include "chi_lua.h"

namespace prk::lua_utils
{
int chiPRKGetParam(lua_State* L);
int chiPRKSetParam(lua_State* L);

chi::InputParameters GetSyntax_SetParam();
chi::ParameterBlock
SetParam(const chi::InputParameters& params);

chi::InputParameters GetParamSyntax();
chi::ParameterBlock
GetParam(const chi::InputParameters& params);
} // namespace prk::lua_utils

#endif