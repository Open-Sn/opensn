#pragma once

#include "framework/lua.h"

namespace opensnlua::lbs
{
int AdjointSolverCreate(lua_State* L);
int AdjointSolverAddResponseFunction(lua_State* L);
int AdjointSolverMakeExpRepFromP1Moments(lua_State* L);
int AdjointSolverExportImportanceMapBinary(lua_State* L);
int AdjointSolverComputeInnerProduct(lua_State* L);

/**Reads flux-moments file to a buffer and returns a handle to that buffer.
 *
 * \param SolverHandle int Handle to the relevant solver.
 * \param FileBaseName string The base-name of the file(s) from which to read
 *                            the flux moments.
 *
 * \return handle int A handle that can be used with
 *                    `AdjointSolverApplyFluxMomentBuffer`.
 */
int AdjointSolverReadFluxMomentsToBuffer(lua_State* L);

/**Applies buffered flux-moments data to the current phi-old.
 *
 * \param SolverHandle int Handle to the relevant solver.
 * \param BufferHandle int The handle to the buffer-position to be applied.
 */
int AdjointSolverApplyFluxMomentBuffer(lua_State* L);
} // namespace opensnlua::lbs
