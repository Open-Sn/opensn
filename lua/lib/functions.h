// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/functions/scalar_material_function.h"
#include "framework/math/functions/scalar_spatial_function.h"
#include "framework/math/functions/scalar_spatial_material_function.h"
#include "framework/math/functions/vector_spatial_function.h"
#include "framework/math/functions/vector_spatial_material_function.h"
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

//

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

//

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

//

class LuaVectorSpatialFunction : public opensn::VectorSpatialFunction
{
public:
  explicit LuaVectorSpatialFunction(const opensn::InputParameters& params);

  std::vector<double> Evaluate(const opensn::Vector3& xyz, int num_components) const override;

private:
  const std::string function_name_;

public:
  static opensn::InputParameters GetInputParameters();
  static std::shared_ptr<LuaVectorSpatialFunction> Create(const opensn::ParameterBlock& params);
};

//

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
