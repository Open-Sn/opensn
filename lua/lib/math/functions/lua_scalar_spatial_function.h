// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/functions/scalar_spatial_function.h"
#include "framework/object_factory.h"

namespace opensnlua
{

class LuaScalarSpatialFunction : public opensn::ScalarSpatialFunction
{
public:
  explicit LuaScalarSpatialFunction(const opensn::InputParameters& params);
  double Evaluate(const opensn::Vector3& xyz) const override;

private:
  const std::string function_name_;

public:
  static opensn::InputParameters GetInputParameters();
  static std::shared_ptr<LuaScalarSpatialFunction> Create(const opensn::ParameterBlock& params);
};

} // namespace opensnlua
