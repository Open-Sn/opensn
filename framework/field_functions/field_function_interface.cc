// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/field_functions/field_function_interface.h"
#include "framework/field_functions/field_function.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

namespace opensn
{

InputParameters
FieldFunctionInterface::GetInputParameters()
{
  InputParameters params;

  params.AddRequiredParameter<std::shared_ptr<FieldFunction>>("field_function", "Field function.");
  params.SetParameterTypeMismatchAllowed("field_function");

  return params;
}

FieldFunctionInterface::FieldFunctionInterface(const InputParameters& params)
  : field_function_(params.GetSharedPtrParam<FieldFunction>("field_function"))
{
}

std::shared_ptr<FieldFunction>
FieldFunctionInterface::GetFieldFunction() const
{
  return field_function_;
}

} // namespace opensn
