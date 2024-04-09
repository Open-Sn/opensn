// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "lua/framework/lua.h"

namespace opensn
{
class InputParameters;
class ParameterBlock;
} // namespace opensn

namespace opensnlua::lbs
{

using namespace opensn;

/**
 * Removes all forward sources currently set within a response evaluator.
 *
 * \param ResponseEvaluatorIndex A handle to a response evaluator
 */
int ClearResponseSources(lua_State* L);

InputParameters GetResponseBufferSyntax();
ParameterBlock AddResponseBuffers(const InputParameters& params);

InputParameters GetResponseSourceSyntax();
ParameterBlock AddResponseSources(const InputParameters& params);

int EvaluateResponse(lua_State* L);

} // namespace opensnlua::lbs
