// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/spatial_weight_function.h"
#include "framework/mesh/mesh.h"
#include <stdexcept>

namespace opensn
{

std::shared_ptr<SpatialWeightFunction>
SpatialWeightFunction::FromCoordinateType(CoordinateSystemType coord_sys)
{
  if (coord_sys == CoordinateSystemType::CARTESIAN)
    return std::make_shared<CartesianSpatialWeightFunction>();
  else if (coord_sys == CoordinateSystemType::SPHERICAL)
    return std::make_shared<SphericalSpatialWeightFunction>();
  else if (coord_sys == CoordinateSystemType::CYLINDRICAL)
    return std::make_shared<CylindricalSpatialWeightFunction>();
  else
    throw std::runtime_error("Undefined coordinate system type");
}

} // namespace opensn
