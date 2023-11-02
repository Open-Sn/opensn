#include "math/SpatialDiscretization/CellMappings/FiniteVolumeMapping.h"

#include "math/SpatialDiscretization/FiniteElement/QuadraturePointData.h"

namespace chi_math::cell_mapping
{

finite_element::VolumetricQuadraturePointData
FiniteVolumeMapping::MakeVolumetricQuadraturePointData() const
{
  return finite_element::VolumetricQuadraturePointData(
    {0},
    {{cell_.centroid_}},
    {{1.0}},
    {{chi_mesh::Vector3(0, 0, 0)}},
    {volume_},
    face_node_mappings_,
    num_nodes_);
}

finite_element::SurfaceQuadraturePointData
FiniteVolumeMapping::MakeSurfaceQuadraturePointData(size_t face_index) const
{
  return finite_element::SurfaceQuadraturePointData({0},
                                                 {{chi_mesh::Vector3(0, 0, 0)}},
                                                 {{1.0}},
                                                 {{chi_mesh::Vector3(0, 0, 0)}},
                                                 {areas_[face_index]},
                                                 {{chi_mesh::Vector3(0, 0, 0)}},
                                                 face_node_mappings_,
                                                 1);
}

} // namespace chi_math
