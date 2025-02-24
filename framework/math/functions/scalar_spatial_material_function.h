// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/functions/function.h"
#include "framework/math/vector3.h"

namespace opensn
{

/// Base class for evaluating functions given material ID and spatial location
class ScalarSpatialMaterialFunction : public Function
{
public:
  static InputParameters GetInputParameters() { return Function::GetInputParameters(); }
  explicit ScalarSpatialMaterialFunction(const InputParameters& params) : Function(params) {}

  /**
   * Evaluate this function
   *
   * \param mat_id The material ID of the cell
   * \param xyz The xyz coordinates of the point where the function is called.
   * \return Function value
   */
  virtual double Evaluate(int mat_id, const Vector3& xyz) const = 0;
};

} // namespace opensn
