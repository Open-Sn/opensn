// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/post_processors/post_processor.h"
#include "framework/field_functions/grid_based_field_function_interface.h"
#include "framework/mesh/logical_volume/logical_volume_interface.h"

namespace opensn
{
class LogicalVolume;
class FieldFunctionGridBased;

class CellVolumeIntegralPostProcessor : public PostProcessor,
                                        public GridBasedFieldFunctionInterface,
                                        public LogicalVolumeInterface
{
public:
  explicit CellVolumeIntegralPostProcessor(const InputParameters& params);

  void Execute(const Event& event_context) override;

protected:
  void Initialize();

  const bool compute_volume_average_;
  bool initialized_ = false;
  std::vector<uint64_t> cell_local_ids_;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<CellVolumeIntegralPostProcessor> Create(const ParameterBlock& params);
};

} // namespace opensn
