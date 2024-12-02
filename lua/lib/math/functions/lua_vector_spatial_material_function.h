// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/functions/vector_spatial_material_function.h"
#include "framework/mesh/mesh_vector.h"
#include "framework/object_factory.h"

namespace opensnlua
{

class LuaVectorSpatialMaterialFunction : public opensn::VectorSpatialMaterialFunction
{
public:
  explicit LuaVectorSpatialMaterialFunction(const opensn::InputParameters& params);

  std::vector<double>
  Evaluate(const opensn::Vector3& xyz, int mat_id, int num_components) const override;

private:
  const std::string function_name_;

public:
  static opensn::InputParameters GetInputParameters();
  static std::shared_ptr<LuaVectorSpatialMaterialFunction>
  Create(const opensn::ParameterBlock& params);
};

} // namespace opensnlua
