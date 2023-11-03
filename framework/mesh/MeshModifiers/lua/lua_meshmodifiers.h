#pragma once

#include "parameters/input_parameters.h"

namespace chi_mesh::lua_utils
{
chi::InputParameters MeshModifiersApply_Syntax();
chi::ParameterBlock MeshModifiersApply(const chi::InputParameters& params);
} // namespace chi_mesh::lua_utils


