#pragma once

#include "framework/lua.h"

namespace prk::lua_utils
{
/**Gets a parameter from the prk::TransientSolver.
 *
 * \param handle int Handle of the solver.
 * \param param_name  string Name of the parameter to retrieve.
 * \return Varying
 */
int chiPRKGetParam(lua_State* L);

/**Gets a parameter from the prk::TransientSolver.
 *
 * \param handle int Handle of the solver.
 * \param param_name  string Name of the parameter to retrieve.
 * \param value Varying The value to be set to the parameter.
 * \return Varying
 */
int chiPRKSetParam(lua_State* L);

chi::InputParameters GetSyntax_SetParam();
chi::ParameterBlock SetParam(const chi::InputParameters& params);

chi::InputParameters GetParamSyntax();
chi::ParameterBlock GetParam(const chi::InputParameters& params);
} // namespace prk::lua_utils
