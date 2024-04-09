// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/functions/spatial_material_function.h"

namespace opensn
{

InputParameters
SpatialMaterialFunction::GetInputParameters()
{
  InputParameters params = Function::GetInputParameters();
  return params;
}

SpatialMaterialFunction::SpatialMaterialFunction(const InputParameters& params) : Function(params)
{
}

} // namespace opensn
