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
class LBSProblem;

/**
 * A class for multi-group isotropic volumetric sources.
 *
 * This class differs from the standard material property sources in that it can be specified
 * for an arbitrary logical volume that may span multiple material regions and with spatial,
 * material id, and group-wise behaviors via a VectorSpatialFunction.
 *
 * The flexibility of this object allows for its use as a standard volumetric source or as
 * a volumetric response function for adjoint calculations.
 */
class VolumetricSource : public Object
{
public:
  explicit VolumetricSource(const InputParameters& params);

  void Initialize(const LBSProblem& lbs_problem);

  /**
   * Evaluate the distributed source for a given cell node for all groups
   *
   * If the cell does not belong to the logical volume tied to this source,
   * a vector of zeros are returned.
   */
  std::vector<double> operator()(const Cell& cell, const Vector3& xyz, int num_groups) const;

  size_t GetNumLocalSubscribers() const { return num_local_subsribers_; }
  size_t GetNumGlobalSubsribers() const { return num_global_subscribers_; }

  const std::vector<uint64_t>& GetSubscribers() const { return subscribers_; }
  std::shared_ptr<LogicalVolume> GetLogicalVolume() const { return logvol_; }
  const std::vector<int>& GetBlockIDs() const { return block_ids_; }

private:
  int id_;

  std::vector<int> block_ids_;
  const std::shared_ptr<LogicalVolume> logvol_;

  std::vector<double> strength_;
  const std::shared_ptr<VectorSpatialFunction> function_;

  size_t num_local_subsribers_ = 0;
  size_t num_global_subscribers_ = 0;

  std::vector<uint64_t> subscribers_;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<VolumetricSource> Create(const ParameterBlock& params);

private:
  static int next_id_;
};

} // namespace opensn
