#include "framework/math/spatial_discretization/cell_mappings/finite_volume_mapping.h"

#include "framework/math/spatial_discretization/finite_element/quadrature_point_data.h"

namespace opensn
{

VolumetricQuadraturePointData
FiniteVolumeMapping::MakeVolumetricQuadraturePointData() const
{
  return VolumetricQuadraturePointData({0},
                                       {{cell_.centroid_}},
                                       {{1.0}},
                                       {{Vector3(0, 0, 0)}},
                                       {volume_},
                                       face_node_mappings_,
                                       num_nodes_);
}

SurfaceQuadraturePointData
FiniteVolumeMapping::MakeSurfaceQuadraturePointData(size_t face_index) const
{
  return SurfaceQuadraturePointData({0},
                                    {{Vector3(0, 0, 0)}},
                                    {{1.0}},
                                    {{Vector3(0, 0, 0)}},
                                    {areas_[face_index]},
                                    {{Vector3(0, 0, 0)}},
                                    face_node_mappings_,
                                    1);
}

} // namespace opensn
