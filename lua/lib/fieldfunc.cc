// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/lib/fieldfunc.h"
#include "framework/field_functions/interpolation/ffinterpolation.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "lua/lib/console.h"
#include <iostream>

using namespace opensn;

namespace opensnlua
{

void
ExportFieldFunctionToVTK(std::shared_ptr<FieldFunctionGridBased> ff, const std::string& base_name)
{
  FieldFunctionGridBased::ExportMultipleToVTK(base_name, {ff});
}

void
ExportMultiFieldFunctionToVTK(const std::vector<std::shared_ptr<FieldFunction>>& ff_handles,
                              const std::string& base_name)
{
  std::vector<std::shared_ptr<const FieldFunctionGridBased>> ffs;
  ffs.reserve(ff_handles.size());
  for (std::size_t i = 0; i < ff_handles.size(); ++i)
  {
    auto ff = std::dynamic_pointer_cast<FieldFunctionGridBased>(ff_handles[i]);
    OpenSnLogicalErrorIf(not ff, "Only grid-based field functions can be exported");
    ffs.push_back(ff);
  }

  FieldFunctionGridBased::ExportMultipleToVTK(base_name, ffs);
}

void
FFInterpolationExportToCSV(const FieldFunctionInterpolation& ffi)
{
  ffi.ExportToCSV(opensn::input_path.stem());
}

std::shared_ptr<FieldFunction>
GetFieldFunctionHandleByName(const std::string& ff_name)
{
  const std::string fname = "fieldfunc.GetHandleByName";
  std::vector<std::shared_ptr<FieldFunction>> matched_ff;
  for (const auto& pff : opensn::field_function_stack)
  {
    if (pff->GetName() == ff_name)
      matched_ff.push_back(pff);
  }

  auto num_ffuncs = matched_ff.size();

  if (num_ffuncs == 0)
  {
    opensn::log.Log0Warning() << fname << ": No field-functions were found that "
                              << "matched the requested name:\"" << ff_name
                              << "\". A null handle will "
                              << "be returned." << std::endl;

    return nullptr;
  }

  if (num_ffuncs > 1)
    opensn::log.Log0Warning() << fname << ": A total of " << num_ffuncs
                              << " field-functions were found that matched the "
                              << " requested name. Only the first match will be "
                              << " returned.";

  return matched_ff.front();
}

} // namespace opensnlua
