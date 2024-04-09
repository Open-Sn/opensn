// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/field_functions/operations/field_operation.h"

namespace opensn
{

// Since there are no input parameters we will not register this object

InputParameters
FieldOperation::GetInputParameters()
{
  return Object::GetInputParameters();
}

FieldOperation::FieldOperation(const InputParameters& params) : Object(params)
{
}

} // namespace opensn
