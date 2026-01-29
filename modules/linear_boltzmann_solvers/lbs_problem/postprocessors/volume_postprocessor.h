// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/parameters/input_parameters.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include <memory>

namespace opensn
{

class VolumePostprocessor
{
public:
  /// Input parameters based construction.
  explicit VolumePostprocessor(const InputParameters& params);

  void Initialize();

  void Execute();

  std::vector<double> GetValue() const;

private:
  std::shared_ptr<LBSProblem> lbs_problem_;
  /// Logical volume associated with this PPS (can be null)
  std::shared_ptr<LogicalVolume> logical_volume_;
  /// Local cell IDs
  std::vector<std::uint32_t> cell_local_ids_;
  /// Groups
  std::vector<unsigned int> groups_;
  /// Computed postprocessed values
  std::vector<double> values_;

public:
  /// Returns the input parameters for this object.
  static InputParameters GetInputParameters();
  static std::shared_ptr<VolumePostprocessor> Create(const ParameterBlock& params);
};

} // namespace opensn
