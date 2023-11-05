#pragma once

#include "opensn/framework/chi_lua.h"

namespace lbs::disc_ord_lua_utils
{
int chiLBSComputeBalance(lua_State* L);
int chiLBSComputeLeakage(lua_State* L);
} // namespace lbs::disc_ord_lua_utils
