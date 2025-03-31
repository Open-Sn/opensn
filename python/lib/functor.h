// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

// This is only a temporary solution!
// Once the Lua API is eradicated and opensn::Function is replace by std::function
// this API is no longer needed.

#pragma once

#include "framework/math/functions/function.h"
#include "framework/math/functions/scalar_material_function.h"
#include "framework/math/functions/scalar_spatial_function.h"
#include "framework/math/functions/scalar_spatial_material_function.h"
#include "framework/math/functions/vector_spatial_function.h"
#include <functional>

namespace opensn
{

/// Bind class for scalar material function
class PySMFunction : public ScalarMaterialFunction
{
public:
  PySMFunction(const std::function<double(double, int)>& func)
    : ScalarMaterialFunction(), func_(func)
  {
  }

  inline double Evaluate(double val, int mat_id) const override { return this->func_(val, mat_id); }

protected:
  std::function<double(double, int)> func_;
};

/// Bind class for scalar spatial function
class PySSFunction : public ScalarSpatialFunction
{
public:
  PySSFunction(const std::function<double(const Vector3&)>& func)
    : ScalarSpatialFunction(), func_(func)
  {
  }

  inline double Evaluate(const Vector3& xyz) const override { return this->func_(xyz); }

protected:
  std::function<double(const Vector3&)> func_;
};

/// Bind class for scalar spatial material function
class PySSMFunction : public ScalarSpatialMaterialFunction
{
public:
  PySSMFunction(const std::function<double(int, const Vector3&)>& func)
    : ScalarSpatialMaterialFunction(), func_(func)
  {
  }

  inline double Evaluate(int mat_id, const Vector3& xyz) const override
  {
    return this->func_(mat_id, xyz);
  }

protected:
  std::function<double(int, const Vector3&)> func_;
};

/// Bind class for vector spatial function
class PyVSFunction : public VectorSpatialFunction
{
public:
  PyVSFunction(const std::function<std::vector<double>(const Vector3&, int)>& func)
    : VectorSpatialFunction(), func_(func)
  {
  }

  inline std::vector<double> Evaluate(const Vector3& xyz, int num_components) const override
  {
    return this->func_(xyz, num_components);
  }

protected:
  std::function<std::vector<double>(const Vector3&, int)> func_;
};

} // namespace opensn
