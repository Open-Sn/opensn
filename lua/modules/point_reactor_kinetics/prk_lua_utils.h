// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "lua/framework/lua.h"

namespace opensnlua::prk
{
/**Gets a parameter from the prk::TransientSolver.
 *
 * \param handle int Handle of the solver.
 * \param param_name  string Name of the parameter to retrieve.
 * \return Varying
 */
int PRKGetParam(lua_State* L);

/**Gets a parameter from the prk::TransientSolver.
 *
 * \param handle int Handle of the solver.
 * \param param_name  string Name of the parameter to retrieve.
 * \param value Varying The value to be set to the parameter.
 * \return Varying
 */
int PRKSetParam(lua_State* L);

opensn::InputParameters GetSyntax_SetParam();
opensn::ParameterBlock SetParam(const opensn::InputParameters& params);

opensn::InputParameters GetParamSyntax();
opensn::ParameterBlock GetParam(const opensn::InputParameters& params);

} // namespace opensnlua::prk
