// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/functions/function.h"
#include "framework/math/vector3.h"

namespace opensn
{

/// Base class for evaluating functions given material ID and single value
class ScalarMaterialFunction : public Function
{
public:
  ScalarMaterialFunction() = default;
  static InputParameters GetInputParameters() { return Function::GetInputParameters(); }
  explicit ScalarMaterialFunction(const InputParameters& params) : Function(params) {}

  /**
   * Evaluate this function
   *
   * \param val The scalar value (for example, a field function value)
   * \param mat_id The material ID of the cell
   * \return Function value
   */
  virtual double Evaluate(double val, int mat_id) const = 0;
};

} // namespace opensn
