// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/functions/function.h"
#include "framework/mesh/mesh_vector.h"

namespace opensn
{

/**
 * Base class for evaluating spatial material functions given a coordinate.
 */
class VectorSpatialFunction : public Function
{
public:
  static InputParameters GetInputParameters() { return Function::GetInputParameters(); }
  explicit VectorSpatialFunction(const InputParameters& params) : Function(params) {}

  /**
   * Evaluate the function.
   *
   * \param xyz The xyz coordinates of the point where the function is evaluated.
   * \param num_components The number of components
   * \return A vector with the function evaluation (should have `num_groups` entries)
   */
  virtual std::vector<double> Evaluate(const Vector3& xyz, int num_components) const = 0;
};

} // namespace opensn
