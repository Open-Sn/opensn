// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/functions/scalar_spatial_material_function.h"
#include "framework/object_factory.h"

namespace opensnlua
{

class LuaScalarSpatialMaterialFunction : public opensn::ScalarSpatialMaterialFunction
{
public:
  explicit LuaScalarSpatialMaterialFunction(const opensn::InputParameters& params);
  double Evaluate(int mat_id, const opensn::Vector3& xyz) const override;

private:
  const std::string function_name_;

public:
  static opensn::InputParameters GetInputParameters();
  static std::shared_ptr<LuaScalarSpatialMaterialFunction>
  Create(const opensn::ParameterBlock& params);
};

} // namespace opensnlua
