#ifndef PWL_POLYGON_VALUES_H
#define PWL_POLYGON_VALUES_H

#include "math/SpatialDiscretization/CellMappings/PieceWiseLinearBaseMapping.h"
#include "math/Quadratures/quadrature_line.h"
#include "math/Quadratures/quadrature_triangle.h"
#include "mesh/Cell/cell.h"

#include <array>

// ###################################################################
namespace chi_math::cell_mapping
{
/** Object for handling polygon shaped 2D cells.
* \ingroup doc_CellMappings*/
class PieceWiseLinearPolygonMapping : public PieceWiseLinearBaseMapping
{
public:
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constructor
  PieceWiseLinearPolygonMapping(const chi_mesh::Cell& poly_cell,
                       const chi_mesh::MeshContinuum& ref_grid,
                       const QuadratureTriangle& volume_quadrature,
                       const QuadratureLine& surface_quadrature);

  finite_element::VolumetricQuadraturePointData
  MakeVolumetricQuadraturePointData() const override;

  finite_element::SurfaceQuadraturePointData
  MakeSurfaceQuadraturePointData(size_t face_index) const override;

  double SideGradShape_x(uint32_t side, uint32_t i) const;
  double SideGradShape_y(uint32_t side, uint32_t i) const;

  double ShapeValue(int i, const chi_mesh::Vector3& xyz) const override;

  chi_mesh::Vector3 GradShapeValue(int i,
                                   const chi_mesh::Vector3& xyz) const override;

  void ShapeValues(const chi_mesh::Vector3& xyz,
                   std::vector<double>& shape_values) const override;

  void GradShapeValues(
    const chi_mesh::Vector3& xyz,
    std::vector<chi_mesh::Vector3>& gradshape_values) const override;

private:
  // ################################################## Define standard
  //                                                    triangle linear shape
  //                                                    functions
  static double TriShape(uint32_t index,
                         const chi_mesh::Vector3& qpoint,
                         bool on_surface = false);

  // ############################################### Shape functions per side
  double SideShape(uint32_t side,
                   uint32_t i,
                   const chi_mesh::Vector3& qpoint,
                   bool on_surface = false) const;

  // This structure goes into sides
  struct FEside_data2d
  {
    double detJ;
    double detJ_surf;
    std::array<uint64_t, 2> v_index;
    chi_mesh::Vector3 v0;
    chi_mesh::Matrix3x3 J;
    chi_mesh::Matrix3x3 Jinv;
    chi_mesh::Matrix3x3 JTinv;
    chi_mesh::Vector3 normal;
  };

  std::vector<FEside_data2d> sides_;
  const QuadratureTriangle& volume_quadrature_;
  const QuadratureLine& surface_quadrature_;

  int num_of_subtris_;
  double beta_;
  chi_mesh::Vertex vc_;
  std::vector<std::vector<int>> node_to_side_map_;
};
} // namespace chi_math::cell_mapping

#endif
