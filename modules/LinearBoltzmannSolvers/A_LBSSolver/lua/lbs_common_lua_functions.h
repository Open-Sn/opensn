#ifndef CHITECH_LBS_COMMON_LUA_FUNCTIONS_H
#define CHITECH_LBS_COMMON_LUA_FUNCTIONS_H

#include "chi_lua.h"

namespace lbs::common_lua_utils
{
  int chiLBSSetProperty(lua_State *L);

  int chiLBSCreateGroupset(lua_State *L);
  int chiLBSCreateGroup(lua_State *L);
  int chiLBSGroupsetAddGroups(lua_State *L);
  int chiLBSGroupsetSetQuadrature(lua_State *L);
  int chiLBSGroupsetSetAngleAggregationType(lua_State *L);
  int chiLBSGroupsetSetAngleAggDiv(lua_State *L);
  int chiLBSGroupsetSetGroupSubsets(lua_State *L);
  int chiLBSGroupsetSetIterativeMethod(lua_State *L);
  int chiLBSGroupsetSetResidualTolerance(lua_State *L);
  int chiLBSGroupsetSetMaxIterations(lua_State *L);
  int chiLBSGroupsetSetGMRESRestartIntvl(lua_State *L);
  int chiLBSGroupsetSetEnableSweepLog(lua_State *L);
  int chiLBSGroupsetSetWGDSA(lua_State *L);
  int chiLBSGroupsetSetTGDSA(lua_State *L);

  int chiLBSGetScalarFieldFunctionList(lua_State *L);

  int chiLBSWriteGroupsetAngularFlux(lua_State *L);
  int chiLBSReadGroupsetAngularFlux(lua_State *L);

  int chiLBSWriteFluxMoments(lua_State *L);
  int chiLBSCreateAndWriteSourceMoments(lua_State *L);
  int chiLBSReadFluxMomentsAndMakeSourceMoments(lua_State *L);
  int chiLBSReadSourceMoments(lua_State *L);
  int chiLBSReadFluxMoments(lua_State *L);

  int chiLBSComputeFissionRate(lua_State *L);
  int chiLBSInitializeMaterials(lua_State* L);

  int chiLBSAddPointSource(lua_State *L);
  int chiLBSClearPointSources(lua_State *L);
  int chiLBSInitializePointSources(lua_State *L);
}

#endif //CHITECH_LBS_COMMON_LUA_FUNCTIONS_H
