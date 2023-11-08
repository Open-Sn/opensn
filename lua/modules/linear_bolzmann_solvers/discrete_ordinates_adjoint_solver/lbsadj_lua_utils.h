#pragma once

#include "framework/chi_lua.h"

namespace lbs::adjoint_lua_utils
{
int chiAdjointSolverCreate(lua_State* L);
int chiAdjointSolverAddResponseFunction(lua_State* L);
int chiAdjointSolverMakeExpRepFromP1Moments(lua_State* L);
int chiAdjointSolverExportImportanceMapBinary(lua_State* L);
int chiAdjointSolverComputeInnerProduct(lua_State* L);

/**Reads flux-moments file to a buffer and returns a handle to that buffer.
 *
 * \param SolverHandle int Handle to the relevant solver.
 * \param FileBaseName string The base-name of the file(s) from which to read
 *                            the flux moments.
 *
 * \return handle int A handle that can be used with
 *                    `chiAdjointSolverApplyFluxMomentBuffer`.
 */
int chiAdjointSolverReadFluxMomentsToBuffer(lua_State* L);

/**Applies buffered flux-moments data to the current phi-old.
 *
 * \param SolverHandle int Handle to the relevant solver.
 * \param BufferHandle int The handle to the buffer-position to be applied.
 */
int chiAdjointSolverApplyFluxMomentBuffer(lua_State* L);
} // namespace lbs::adjoint_lua_utils
