// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/spatial_discretization/cell_mappings/finite_volume/finite_volume_mapping.h"
#include "framework/math/spatial_discretization/finite_element/finite_element_data.h"

namespace opensn
{

VolumetricFiniteElementData
FiniteVolumeMapping::MakeVolumetricFiniteElementData() const
{
  return VolumetricFiniteElementData({0},
                                     {{cell_.centroid_}},
                                     {{1.0}},
                                     {{Vector3(0, 0, 0)}},
                                     {volume_},
                                     face_node_mappings_,
                                     num_nodes_);
}

SurfaceFiniteElementData
FiniteVolumeMapping::MakeSurfaceFiniteElementData(size_t face_index) const
{
  return SurfaceFiniteElementData({0},
                                  {{Vector3(0, 0, 0)}},
                                  {{1.0}},
                                  {{Vector3(0, 0, 0)}},
                                  {areas_[face_index]},
                                  {{Vector3(0, 0, 0)}},
                                  face_node_mappings_,
                                  1);
}

} // namespace opensn
