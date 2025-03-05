// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/object.h"

namespace opensn
{

/// Base class for functions
class Function : public Object
{
public:
  Function() = default;
  static InputParameters GetInputParameters() { return Object::GetInputParameters(); }

protected:
  explicit Function(const InputParameters& params) : Object(params) {}
};

} // namespace opensn
