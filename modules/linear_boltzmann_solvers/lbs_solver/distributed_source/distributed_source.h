// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/object.h"

#include <vector>
#include <iostream>

namespace opensn
{
struct Vector3;
class Cell;
class LogicalVolume;
class SpatialMaterialFunction;

namespace lbs
{
class LBSSolver;

/**
 * A class for multi-group isotropic distributed sources.
 *
 * This class differs from the standard material property sources in that it can be specified
 * for an arbitrary logical volume that may span multiple material regions and with spatial,
 * material id, and group-wise behaviors via a SpatialMaterialFunction.
 *
 * The flexibility of this object allows for its use as a standard distributed source or as
 * a volumetric response function for adjoint calculations.
 */
class DistributedSource : public Object
{
public:
  static InputParameters GetInputParameters();
  explicit DistributedSource(const InputParameters& params);

  void Initialize(const LBSSolver& lbs_solver);

  /**
   * Evaluate the distributed source for a given cell node for all groups
   *
   * If the cell does not belong to the logical volume tied to this source,
   * a vector of zeros are returned.
   */
  std::vector<double> operator()(const Cell& cell, const Vector3& xyz, const int num_groups) const;

  size_t NumLocalSubscribers() const { return num_local_subsribers_; }
  size_t NumGlobalSubsribers() const { return num_global_subscribers_; }

  const std::vector<uint64_t>& Subscribers() const { return subscribers_; }
  const std::shared_ptr<opensn::LogicalVolume>& LogicalVolume() const
  {
    return logical_volume_ptr_;
  }

private:
  const std::shared_ptr<opensn::LogicalVolume> logical_volume_ptr_;
  const std::shared_ptr<SpatialMaterialFunction> function_;

  size_t num_local_subsribers_ = 0;
  size_t num_global_subscribers_ = 0;

  std::vector<uint64_t> subscribers_;
};

} // namespace lbs
} // namespace opensn
