// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/object.h"

namespace opensn
{

/**The base field operation class.*/
class FieldOperation : public Object
{
public:
  /**Returns the input parameters.*/
  static InputParameters GetInputParameters();

  /**Constructor.*/
  explicit FieldOperation(const InputParameters& params);

  virtual void Execute() = 0;

  virtual ~FieldOperation() = default;
};

} // namespace opensn
