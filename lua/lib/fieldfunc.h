// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/field_functions/field_function.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/field_functions/interpolation/ffinterpolation.h"

namespace opensnlua
{

void ExportFieldFunctionToVTK(std::shared_ptr<opensn::FieldFunctionGridBased> ff,
                              const std::string& base_name);
void
ExportMultiFieldFunctionToVTK(const std::vector<std::shared_ptr<opensn::FieldFunction>>& ff_handles,
                              const std::string& base_name);
void FFInterpolationExportToCSV(const opensn::FieldFunctionInterpolation& ffi);

std::shared_ptr<opensn::FieldFunction> GetFieldFunctionHandleByName(const std::string& ff_name);

} // namespace opensnlua
