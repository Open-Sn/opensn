// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/parameters/input_parameters.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include <memory>

namespace opensn
{

class SurfacePostprocessor
{
public:
  enum class ValueType
  {
    INTEGRAL,
    MAX,
    MIN,
    AVERAGE
  };

  enum class CurrentType
  {
    NET,
    INCOMING,
    OUTGOING
  };

  /// Input parameters based construction.
  explicit SurfacePostprocessor(const InputParameters& params);

  void Execute();

  std::vector<std::vector<double>> GetValue() const;

private:
  struct SideFace
  {
    uint64_t cell_local_index;
    int local_face;
  };

  void CreateSpatialRestriction();
  void CreateEnergyRestriction();
  std::vector<SurfacePostprocessor::SideFace>
  GetLogicalVolumeSides(const std::vector<SideFace>& all_side_faces,
                        std::shared_ptr<LogicalVolume> log_vol);
  std::vector<double> ComputeIntegral(const std::vector<SideFace>& side_faces);
  std::vector<double> ComputeMin(const std::vector<SideFace>& side_faces);
  std::vector<double> ComputeMax(const std::vector<SideFace>& side_faces);
  std::vector<double> ComputeAvg(const std::vector<SideFace>& side_faces);

  std::shared_ptr<DiscreteOrdinatesProblem> do_problem_;
  ///
  std::vector<std::string> boundary_names_;
  /// Boundary IDs this postprocessor is restricted to
  std::vector<uint64_t> boundary_ids_;
  /// Logical volume associated with this PPS (can be null)
  std::vector<std::shared_ptr<LogicalVolume>> logical_volumes_;
  ///
  std::vector<std::vector<SideFace>> side_faces_;
  /// Groups
  std::vector<unsigned int> groups_;
  /// Groupset IDs
  std::vector<unsigned int> groupset_ids_;
  /// Type of value to compute
  ValueType value_type_;
  /// Type of current
  CurrentType current_type_;
  /// Computed postprocessed values
  std::vector<std::vector<double>> values_;
  /// Selected group (-1 means not selected)
  std::optional<unsigned int> selected_group_;
  /// Selected groupset (-1 means not selected)
  std::optional<unsigned int> selected_groupset_;

public:
  /// Returns the input parameters for this object.
  static InputParameters GetInputParameters();
  static std::shared_ptr<SurfacePostprocessor> Create(const ParameterBlock& params);
};

} // namespace opensn
