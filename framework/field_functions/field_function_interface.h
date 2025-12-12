// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/parameters/input_parameters.h"

namespace opensn
{

class FieldFunction;

/**
 * Interface class to add a dependency on a logical volume. Two things need to
 * be done to use this interface. 1) Derive from it. 2) Add its parameters to
 * the child class. Now it will require a handle to a FieldFunction in
 * the input language.
 */
class FieldFunctionInterface
{
protected:
  explicit FieldFunctionInterface(const InputParameters& params);

  std::shared_ptr<FieldFunction> GetFieldFunction() const;

private:
  std::shared_ptr<FieldFunction> field_function_;

protected:
  static InputParameters GetInputParameters();
};

} // namespace opensn
