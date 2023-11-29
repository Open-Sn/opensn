#include "lbts_lua_utils.h"

#if 0
#define LUA_FMACRO1(x) lua_register(L, #x, lbs::lbts_lua_utils::x)
#define LUA_CMACRO1(x, y)                                                                          \
  lua_pushnumber(L, y);                                                                            \
  lua_setglobal(L, #x)

#define LUA_CTABLE1(x)                                                                             \
  lua_newtable(L);                                                                                 \
  lua_setglobal(L, #x)

#define LUA_CADDCONST_VALUE_TO_TABLE1(const_name, const_value, namespace_name)                     \
  lua_getglobal(L, #namespace_name);                                                               \
  lua_pushstring(L, #const_name);                                                                  \
  lua_pushnumber(L, const_value);                                                                  \
  lua_settable(L, -3);                                                                             \
  lua_pop(L, 1)

void
lbs::lbts_lua_utils::RegisterLuaEntities(lua_State* L)
{
  LUA_FMACRO1(LBSCreateTransientSolver);
  LUA_FMACRO1(LBTSSetProperty);
  LUA_FMACRO1(LBTSGetProperty);
  LUA_FMACRO1(chiLBTSAdvanceTimeData);
}
#endif
