// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/spatial_discretization/cell_mappings/cell_mapping.h"
#include "framework/mesh/cell/cell.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include <utility>

namespace opensn
{

/**
 * Cell mapping for a finite volume representation of a cell.
 *
 * \ingroup doc_CellMappings
 */
class FiniteVolumeMapping : public CellMapping
{
public:
  explicit FiniteVolumeMapping(const std::shared_ptr<MeshContinuum> grid,
                               const Cell& cell,
                               const Vector3& cc,
                               std::vector<std::vector<int>> face_node_mappings)
    : CellMapping(grid, cell, 1, {cell.centroid}, std::move(face_node_mappings))
  {
  }

  double ShapeValue(size_t i, const Vector3& xyz) const override
  {
    return grid_->CheckPointInsideCell(cell_, xyz) ? 1.0 : 0.0;
  }

  void ShapeValues(const Vector3& xyz, Vector<double>& shape_values) const override
  {
    if (grid_->CheckPointInsideCell(cell_, xyz))
      shape_values = Vector<double>(num_nodes_, 1.0);
    else
      shape_values = Vector<double>(num_nodes_, 0.0);
  }

  Vector3 GradShapeValue(size_t i, const Vector3& xyz) const override
  {
    return Vector3(0.0, 0.0, 0.0);
  }

  void GradShapeValues(const Vector3& xyz, std::vector<Vector3>& gradshape_values) const override
  {
    gradshape_values.assign(num_nodes_, Vector3(0, 0, 0));
  }

  VolumetricFiniteElementData MakeVolumetricFiniteElementData() const override;

  SurfaceFiniteElementData MakeSurfaceFiniteElementData(size_t face_index) const override;
};

} // namespace opensn
