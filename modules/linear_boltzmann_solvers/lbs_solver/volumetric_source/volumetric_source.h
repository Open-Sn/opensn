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
class VectorSpatialFunction;
class LBSSolver;

/**
 * A class for multi-group isotropic volumetric sources.
 *
 * This class differs from the standard material property sources in that it can be specified
 * for an arbitrary logical volume that may span multiple material regions and with spatial,
 * material id, and group-wise behaviors via a VectorSpatialMaterialFunction.
 *
 * The flexibility of this object allows for its use as a standard volumetric source or as
 * a volumetric response function for adjoint calculations.
 */
class VolumetricSource : public Object
{
public:
  static InputParameters GetInputParameters();
  explicit VolumetricSource(const InputParameters& params);

  void Initialize(const LBSSolver& lbs_solver);

  /**
   * Evaluate the distributed source for a given cell node for all groups
   *
   * If the cell does not belong to the logical volume tied to this source,
   * a vector of zeros are returned.
   */
  std::vector<double> operator()(const Cell& cell, const Vector3& xyz, int num_groups) const;

  size_t NumLocalSubscribers() const { return num_local_subsribers_; }
  size_t NumGlobalSubsribers() const { return num_global_subscribers_; }

  const std::vector<uint64_t>& Subscribers() const { return subscribers_; }
  std::shared_ptr<LogicalVolume> GetLogicalVolume() const { return logvol_; }
  const std::vector<int>& BlockIDs() const { return block_ids_; }

private:
  std::vector<int> block_ids_;
  const std::shared_ptr<LogicalVolume> logvol_;

  std::vector<double> strength_;
  const std::shared_ptr<VectorSpatialFunction> function_;

  size_t num_local_subsribers_ = 0;
  size_t num_global_subscribers_ = 0;

  std::vector<uint64_t> subscribers_;
};

} // namespace opensn
