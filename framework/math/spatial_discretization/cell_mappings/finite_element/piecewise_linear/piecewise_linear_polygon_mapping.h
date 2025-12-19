// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/spatial_discretization/cell_mappings/finite_element/piecewise_linear/piecewise_linear_base_mapping.h"
#include "framework/data_types/matrix3x3.h"
#include <array>

namespace opensn
{

class Cell;
class TriangleQuadrature;
class LineQuadrature;

/**
 * Object for handling polygon shaped 2D cells.
 *
 * \ingroup doc_CellMappings
 */
class PieceWiseLinearPolygonMapping : public PieceWiseLinearBaseMapping
{
public:
  PieceWiseLinearPolygonMapping(const Cell& poly_cell,
                                std::shared_ptr<MeshContinuum> ref_grid,
                                const TriangleQuadrature& volume_quadrature,
                                const LineQuadrature& surface_quadrature);

  VolumetricFiniteElementData MakeVolumetricFiniteElementData() const override;

  SurfaceFiniteElementData MakeSurfaceFiniteElementData(size_t face_index) const override;

  /// Pre-computation of the partial derivative along x of the shape function at a quadrature point.
  double SideGradShape_x(size_t side, size_t i) const;

  /// Pre-computation of the partial derivative along y of the shape function at a quadrature point.
  double SideGradShape_y(size_t side, size_t i) const;

  double ShapeValue(size_t i, const Vector3& xyz) const override;

  Vector3 GradShapeValue(size_t i, const Vector3& xyz) const override;

  void ShapeValues(const Vector3& xyz, Vector<double>& shape_values) const override;

  void GradShapeValues(const Vector3& xyz, std::vector<Vector3>& gradshape_values) const override;

private:
  /// Precomputation of the shape function at a quadrature point.
  double SideShape(size_t side, size_t i, const Vector3& qpoint, bool on_surface = false) const;

  /// This structure goes into sides
  struct FEside_data2d
  {
    double detJ{};
    double detJ_surf{};
    std::array<uint64_t, 2> v_index{};
    Vector3 v0;
    Matrix3x3 J;
    Matrix3x3 Jinv;
    Matrix3x3 JTinv;
    Vector3 normal;
  };

  std::vector<FEside_data2d> sides_;
  const TriangleQuadrature& volume_quadrature_;
  const LineQuadrature& surface_quadrature_;

  std::size_t num_of_subtris_;
  double beta_;
  Vector3 vc_;
  std::vector<std::vector<int>> node_to_side_map_;

private:
  /// Define standard triangle linear shape functions
  static double TriShape(uint32_t index, const Vector3& qpoint, bool on_surface = false);
};

} // namespace opensn
