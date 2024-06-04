// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/functions/function.h"

namespace opensn
{

InputParameters
Function::GetInputParameters()
{
  InputParameters params = Object::GetInputParameters();
  return params;
}

Function::Function(const InputParameters& params) : Object(params)
{
}

} // namespace opensn
