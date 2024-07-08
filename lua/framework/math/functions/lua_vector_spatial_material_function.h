// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/functions/vector_spatial_material_function.h"
#include "framework/mesh/mesh_vector.h"

namespace opensnlua
{

class LuaVectorSpatialMaterialFunction : public opensn::VectorSpatialMaterialFunction
{
public:
  static opensn::InputParameters GetInputParameters();
  explicit LuaVectorSpatialMaterialFunction(const opensn::InputParameters& params);

  std::vector<double>
  Evaluate(const opensn::Vector3& xyz, int mat_id, int num_components) const override;

private:
  const std::string lua_function_name_;
};

} // namespace opensnlua
