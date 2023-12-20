#pragma once

#include "framework/parameters/input_parameters.h"

namespace opensnlua
{
opensn::InputParameters MeshModifiersApply_Syntax();
opensn::ParameterBlock MeshModifiersApply(const opensn::InputParameters& params);
} // namespace opensnlua
