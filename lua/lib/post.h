// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/parameters/input_parameters.h"
#include "framework/post_processors/post_processor.h"

namespace opensnlua
{

void PostProcessorPrinterSetOptions(const opensn::InputParameters& params);

void PrintPostProcessors(
  const std::vector<std::shared_ptr<opensn::PostProcessor>>& postprocessor_ptr_list);

void ExecutePostProcessors(
  const std::vector<std::shared_ptr<opensn::PostProcessor>>& postprocessor_ptr_list);

} // namespace opensnlua
