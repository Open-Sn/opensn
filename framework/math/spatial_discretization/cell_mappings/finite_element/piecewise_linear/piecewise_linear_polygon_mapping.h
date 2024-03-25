#pragma once

#include "framework/math/spatial_discretization/cell_mappings/finite_element/piecewise_linear/piecewise_linear_base_mapping.h"
#include "framework/math/quadratures/spatial/line_quadrature.h"
#include "framework/math/quadratures/spatial/triangle_quadrature.h"
#include "framework/mesh/cell/cell.h"
#include <array>

namespace opensn
{

/**
 * Object for handling polygon shaped 2D cells.
 * \ingroup doc_CellMappings
 */
class PieceWiseLinearPolygonMapping : public PieceWiseLinearBaseMapping
{
public:
  PieceWiseLinearPolygonMapping(const Cell& poly_cell,
                                const MeshContinuum& ref_grid,
                                const QuadratureTriangle& volume_quadrature,
                                const QuadratureLine& surface_quadrature);

  VolumetricFiniteElementData MakeVolumetricFiniteElementData() const override;

  SurfaceFiniteElementData MakeSurfaceFiniteElementData(size_t face_index) const override;

  /**
   * Precomputation of the partial derivative along x of the
   * shape function at a quadrature point.
   */
  double SideGradShape_x(uint32_t side, uint32_t i) const;

  /**
   * Precomputation of the partial derivative along y of the shape function at a quadrature point.
   */
  double SideGradShape_y(uint32_t side, uint32_t i) const;

  double ShapeValue(int i, const Vector3& xyz) const override;

  Vector3 GradShapeValue(int i, const Vector3& xyz) const override;

  void ShapeValues(const Vector3& xyz, std::vector<double>& shape_values) const override;

  void GradShapeValues(const Vector3& xyz, std::vector<Vector3>& gradshape_values) const override;

private:
  /**
   * Define standard triangle linear shape functions
   */
  static double TriShape(uint32_t index, const Vector3& qpoint, bool on_surface = false);

  /**
   * Precomputation of the shape function at a quadrature point.
   */
  double SideShape(uint32_t side, uint32_t i, const Vector3& qpoint, bool on_surface = false) const;

  // This structure goes into sides
  struct FEside_data2d
  {
    double detJ;
    double detJ_surf;
    std::array<uint64_t, 2> v_index;
    Vector3 v0;
    Matrix3x3 J;
    Matrix3x3 Jinv;
    Matrix3x3 JTinv;
    Vector3 normal;
  };

  std::vector<FEside_data2d> sides_;
  const QuadratureTriangle& volume_quadrature_;
  const QuadratureLine& surface_quadrature_;

  int num_of_subtris_;
  double beta_;
  Vertex vc_;
  std::vector<std::vector<int>> node_to_side_map_;
};

} // namespace opensn
