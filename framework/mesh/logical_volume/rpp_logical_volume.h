// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/logical_volume/logical_volume.h"

namespace opensn
{

/**Rectangular Parallel Piped (RPP) logical volume*/
class RPPLogicalVolume : public LogicalVolume
{
  enum CornerName
  {
    XMAX = 0,
    XMIN = 1,
    YMAX = 2,
    YMIN = 3,
    ZMAX = 4,
    ZMIN = 5
  };

public:
  static InputParameters GetInputParameters();
  explicit RPPLogicalVolume(const InputParameters& params);

  bool Inside(const Vector3& point) const override;

protected:
  double xmin_, xmax_;
  double ymin_, ymax_;
  double zmin_, zmax_;
  bool infx_, infy_, infz_;
};

} // namespace opensn
