// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/functions/scalar_material_function.h"
#include "framework/object_factory.h"

namespace opensnlua
{

class LuaScalarMaterialFunction : public opensn::ScalarMaterialFunction
{
public:
  explicit LuaScalarMaterialFunction(const opensn::InputParameters& params);
  double Evaluate(double val, int mat_id) const override;

private:
  const std::string function_name_;

public:
  static opensn::InputParameters GetInputParameters();
  static std::shared_ptr<LuaScalarMaterialFunction> Create(const opensn::ParameterBlock& params);
};

} // namespace opensnlua
