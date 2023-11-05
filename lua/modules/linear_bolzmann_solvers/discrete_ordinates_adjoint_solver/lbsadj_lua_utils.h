#pragma once

#include "opensn/framework/chi_lua.h"

namespace lbs::adjoint_lua_utils
{
int chiAdjointSolverCreate(lua_State* L);
int chiAdjointSolverAddResponseFunction(lua_State* L);
int chiAdjointSolverMakeExpRepFromP1Moments(lua_State* L);
int chiAdjointSolverExportImportanceMapBinary(lua_State* L);
int chiAdjointSolverComputeInnerProduct(lua_State* L);

int chiAdjointSolverReadFluxMomentsToBuffer(lua_State* L);
int chiAdjointSolverApplyFluxMomentBuffer(lua_State* L);
} // namespace lbs::adjoint_lua_utils
