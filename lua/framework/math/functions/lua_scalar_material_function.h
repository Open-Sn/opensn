// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/functions/scalar_material_function.h"

namespace opensnlua
{

class LuaScalarMaterialFunction : public opensn::ScalarMaterialFunction
{
public:
  static opensn::InputParameters GetInputParameters();
  explicit LuaScalarMaterialFunction(const opensn::InputParameters& params);
  double Evaluate(double val, int mat_id) const override;

private:
  const std::string lua_function_name_;
};

} // namespace opensnlua
