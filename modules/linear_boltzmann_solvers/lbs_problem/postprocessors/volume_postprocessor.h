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
  enum class ValueType
  {
    INTEGRAL,
    MAX,
    MIN,
    AVERAGE
  };

  /// Input parameters based construction.
  explicit VolumePostprocessor(const InputParameters& params);

  void Execute();

  std::vector<std::vector<double>> GetValue() const;

private:
  void CreateSpatialRestriction();
  void CreateEnergyRestriction();
  std::vector<std::uint32_t> GetLogivalVolumeCellIDs(std::shared_ptr<LogicalVolume> log_vol);

  std::shared_ptr<LBSProblem> lbs_problem_;
  /// Block IDs this postprocessor is restricted to
  std::vector<int> block_ids_;
  /// Logical volume associated with this PPS (can be null)
  std::vector<std::shared_ptr<LogicalVolume>> logical_volumes_;
  /// Local cell IDs
  std::vector<std::vector<std::uint32_t>> cell_local_ids_;
  /// Groups
  std::vector<unsigned int> groups_;
  /// Type of value to compute
  ValueType value_type_;
  /// Computed postprocessed values
  std::vector<std::vector<double>> values_;
  /// Selected group (-1 means not selected)
  std::optional<unsigned int> selected_group_;
  /// Selected groupset (-1 means not selected)
  std::optional<unsigned int> selected_groupset_;

  // Helper functions for different computation types
  std::vector<double> ComputeIntegral(const std::vector<uint32_t>& cell_local_ids);
  std::vector<double> ComputeMax(const std::vector<uint32_t>& cell_local_ids);
  std::vector<double> ComputeMin(const std::vector<uint32_t>& cell_local_ids);
  std::vector<double> ComputeVolumeWeightedAverage(const std::vector<uint32_t>& cell_local_ids);

public:
  /// Returns the input parameters for this object.
  static InputParameters GetInputParameters();
  static std::shared_ptr<VolumePostprocessor> Create(const ParameterBlock& params);
};

} // namespace opensn
