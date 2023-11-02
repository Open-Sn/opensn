#ifndef CHITECH_QUADRATURES_LUA_H
#define CHITECH_QUADRATURES_LUA_H

int chiCreateCustomAngularQuadrature(lua_State *L);
int chiCreateCylindricalProductQuadrature(lua_State *L);
int chiCreateSphericalProductQuadrature(lua_State *L);
int chiCreateProductQuadrature(lua_State *L);
int chiCreateLineQuadrature(lua_State *L);
int chiGetProductQuadrature(lua_State *L);
int chiOptimizeAngularQuadratureForPolarSymmetry(lua_State* L);

#endif //CHITECH_QUADRATURES_LUA_H
