#pragma once

#include <utility>

#include "CellMapping.h"

#include "mesh/Cell/cell.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

// ######################################################### Class def
namespace chi_math::cell_mapping
{

/**Cell mapping for a finite volume representation of a cell.
 * \ingroup doc_CellMappings*/
class FiniteVolumeMapping : public CellMapping
{
public:
  explicit FiniteVolumeMapping(const chi_mesh::MeshContinuum& grid,
                               const chi_mesh::Cell& cell,
                               const chi_mesh::Vector3& cc,
                               std::vector<std::vector<int>> face_node_mappings)
    : CellMapping(grid,
                  cell,
                  1,
                  {cell.centroid_},
                  std::move(face_node_mappings),
                  &CellMapping::ComputeCellVolumeAndAreas)
  {
  }

  // 02 Shapefuncs
  double ShapeValue(int i, const chi_mesh::Vector3& xyz) const override
  {
    if (ref_grid_.CheckPointInsideCell(cell_, xyz)) return 1.0;
    else
      return 0.0;
  }
  void ShapeValues(const chi_mesh::Vector3& xyz, std::vector<double>& shape_values) const override
  {
    if (ref_grid_.CheckPointInsideCell(cell_, xyz)) shape_values.assign(num_nodes_, 1.0);
    else
      shape_values.assign(num_nodes_, 0.0);
  }
  chi_mesh::Vector3 GradShapeValue(int i, const chi_mesh::Vector3& xyz) const override
  {
    return chi_mesh::Vector3(0.0, 0.0, 0.0);
  }
  void GradShapeValues(const chi_mesh::Vector3& xyz,
                       std::vector<chi_mesh::Vector3>& gradshape_values) const override
  {
    gradshape_values.assign(num_nodes_, chi_mesh::Vector3(0, 0, 0));
  }

  // 03 Quadrature
  finite_element::VolumetricQuadraturePointData MakeVolumetricQuadraturePointData() const override;

  finite_element::SurfaceQuadraturePointData
  MakeSurfaceQuadraturePointData(size_t face_index) const override;
};

} // namespace chi_math::cell_mapping

