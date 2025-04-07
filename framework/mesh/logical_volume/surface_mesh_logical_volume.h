// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/logical_volume/logical_volume.h"

namespace opensn
{

class SurfaceMesh;

/// SurfaceMesh volume
class SurfaceMeshLogicalVolume : public LogicalVolume
{
public:
  explicit SurfaceMeshLogicalVolume(const InputParameters& params);

  bool Inside(const Vector3& point) const override;

private:
  const std::shared_ptr<SurfaceMesh> surf_mesh_ = nullptr;
  std::array<double, 2> xbounds_;
  std::array<double, 2> ybounds_;
  std::array<double, 2> zbounds_;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<SurfaceMeshLogicalVolume> Create(const ParameterBlock& params);
};

} // namespace opensn
