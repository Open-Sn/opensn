#ifndef LBS_LUA_UTILS_H
#define LBS_LUA_UTILS_H

#include "chi_lua.h"

namespace lbs::disc_ord_lua_utils
{
int chiLBSComputeBalance(lua_State* L);
int chiLBSComputeLeakage(lua_State* L);
} // namespace lbs::disc_ord_lua_utils

#endif
