// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/spatial_weight_function.h"

namespace opensn
{

std::shared_ptr<SpatialWeightFunction>
SpatialWeightFunction::FromGeometryType(GeometryType geometry_type)
{
  if (geometry_type == GeometryType::ONED_SPHERICAL)
    return std::make_shared<SphericalWeightFunction>();
  else if (geometry_type == GeometryType::TWOD_CYLINDRICAL)
    return std::make_shared<CylindricalWeightFunction>();
  else
    return std::make_shared<SpatialWeightFunction>();
}

} // namespace opensn
