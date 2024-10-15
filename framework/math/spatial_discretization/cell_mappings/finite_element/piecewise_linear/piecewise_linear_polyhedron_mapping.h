// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/spatial_discretization/cell_mappings/finite_element/piecewise_linear/piecewise_linear_base_mapping.h"
#include "framework/math/quadratures/spatial/tetrahedra_quadrature.h"
#include "framework/math/quadratures/spatial/triangle_quadrature.h"
#include "framework/mesh/cell/cell.h"

namespace opensn
{

/**
 * Object for handling piecewise linear shape functions on polyhedron shaped 3D cells.
 *
 * \ingroup doc_CellMappings
 */
class PieceWiseLinearPolyhedronMapping : public PieceWiseLinearBaseMapping
{
public:
  /// Constructor for the Piecewise Linear Polyhedron cell finite element view.
  PieceWiseLinearPolyhedronMapping(const Cell& polyh_cell,
                                   const MeshContinuum& ref_grid,
                                   const TetrahedraQuadrature& volume_quadrature,
                                   const TriangleQuadrature& surface_quadrature);

  VolumetricFiniteElementData MakeVolumetricFiniteElementData() const override;

  SurfaceFiniteElementData MakeSurfaceFiniteElementData(size_t face_index) const override;

  /// Actual shape functions as function of cartesian coordinates
  double ShapeValue(int i, const Vector3& xyz) const override;

  Vector3 GradShapeValue(int i, const Vector3& xyz) const override;

  void ShapeValues(const Vector3& xyz, std::vector<double>& shape_values) const override;

  void GradShapeValues(const Vector3& xyz, std::vector<Vector3>& gradshape_values) const override;

private:
  /// Define standard tetrahedron linear shape functions
  static double TetShape(int index, const Vector3& qpoint, bool on_surface = false);

  static double TetGradShape_x(int index);
  static double TetGradShape_y(int index);
  static double TetGradShape_z(int index);

  /// Precomputes the shape function values of a face-side pair at a quadrature point
  double FaceSideShape(uint32_t face_index,
                       uint32_t side_index,
                       uint32_t i,
                       const Vector3& qpoint,
                       bool on_surface = false) const;

  /// Precomputes the gradx-shape function values of a face-side pair at a quadrature point
  double FaceSideGradShape_x(uint32_t face_index, uint32_t side_index, uint32_t i) const;

  /// Precomputes the grady-shape function values of a face-side pair at a quadrature point
  double FaceSideGradShape_y(uint32_t face_index, uint32_t side_index, uint32_t i) const;

  /// Precomputes the gradz-shape function values of a face-side pair at a quadrature point
  double FaceSideGradShape_z(uint32_t face_index, uint32_t side_index, uint32_t i) const;

  /// Stores the data for each side's tetrahedron.
  struct FEside_data3d
  {
    double detJ = 0.0;
    double detJ_surf = 0.0;
    std::vector<uint64_t> v_index;
    Vector3 v0;
    Matrix3x3 J;
    Matrix3x3 Jinv;
    Matrix3x3 JTinv;
  };

  /// Stores data for each face.
  struct FEface_data
  {
    std::vector<FEside_data3d> sides;
    Vector3 normal;
  };

  /// Lowest level of mapping dof i.
  struct FEnodeSideMap
  {
    int index = -1;
    bool part_of_face = false;
  };

  /// Intermediate level of mapping.
  struct FEnodeFaceMap
  {
    std::vector<FEnodeSideMap> side_map;
  };

  /// Node map per face.
  struct FEnodeMap
  {
    std::vector<FEnodeFaceMap> face_map;
  };
  // Goes into node_maps
  //  node n
  //  face f
  //  side s
  //  node_maps[n]->face_map[f]->side_map[s]

  /// Face Beta-factor.
  std::vector<double> face_betaf_;
  /// Cell alpha-factor.
  double alphac_;
  /// Holds determinants and data tet-by-tet.
  std::vector<FEface_data> face_data_;
  /// Maps nodes to side tets.
  std::vector<FEnodeMap> node_side_maps_;

  const TetrahedraQuadrature& volume_quadrature_;
  const TriangleQuadrature& surface_quadrature_;
};

} // namespace opensn
