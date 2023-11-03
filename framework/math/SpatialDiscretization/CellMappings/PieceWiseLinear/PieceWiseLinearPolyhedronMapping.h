#pragma once

#include "math/SpatialDiscretization/CellMappings/PieceWiseLinearBaseMapping.h"
#include "math/Quadratures/quadrature_tetrahedron.h"
#include "math/Quadratures/quadrature_triangle.h"
#include "mesh/Cell/cell.h"

// ###################################################################
namespace chi_math::cell_mapping
{
/**Object for handling piecewise linear
 * shape functions on polyhedron shaped 3D cells.
 * \ingroup doc_CellMappings*/
class PieceWiseLinearPolyhedronMapping : public PieceWiseLinearBaseMapping
{
public:
  // 00_constrdestr.cc
  PieceWiseLinearPolyhedronMapping(const chi_mesh::Cell& polyh_cell,
                                   const chi_mesh::MeshContinuum& ref_grid,
                                   const QuadratureTetrahedron& volume_quadrature,
                                   const QuadratureTriangle& surface_quadrature);

  finite_element::VolumetricQuadraturePointData MakeVolumetricQuadraturePointData() const override;

  finite_element::SurfaceQuadraturePointData
  MakeSurfaceQuadraturePointData(size_t face_index) const override;

  // ############################################### Actual shape functions
  //                                                 as function of cartesian
  //                                                 coordinates
  double ShapeValue(int i, const chi_mesh::Vector3& xyz) const override;

  chi_mesh::Vector3 GradShapeValue(int i, const chi_mesh::Vector3& xyz) const override;

  void ShapeValues(const chi_mesh::Vector3& xyz, std::vector<double>& shape_values) const override;

  void GradShapeValues(const chi_mesh::Vector3& xyz,
                       std::vector<chi_mesh::Vector3>& gradshape_values) const override;

private:
  // ################################################## Define standard
  //                                                    tetrahedron linear shape
  //                                                    functions
  // 01a_reftet.cc
  static double TetShape(uint32_t index, const chi_mesh::Vector3& qpoint, bool on_surface = false);

  static double TetGradShape_x(uint32_t index);
  static double TetGradShape_y(uint32_t index);
  static double TetGradShape_z(uint32_t index);

  // ################################################## Shape functions per
  // face-side 01b_sidevalues.cc
  double FaceSideShape(uint32_t face_index,
                       uint32_t side_index,
                       uint32_t i,
                       const chi_mesh::Vector3& qpoint,
                       bool on_surface = false) const;

  double FaceSideGradShape_x(uint32_t face_index, uint32_t side_index, uint32_t i) const;

  double FaceSideGradShape_y(uint32_t face_index, uint32_t side_index, uint32_t i) const;

  double FaceSideGradShape_z(uint32_t face_index, uint32_t side_index, uint32_t i) const;
  /**Stores the data for each side's tetrahedron. */
  struct FEside_data3d
  {
    double detJ = 0.0;
    double detJ_surf = 0.0;
    std::vector<uint64_t> v_index;
    chi_mesh::Vector3 v0;
    chi_mesh::Matrix3x3 J;
    chi_mesh::Matrix3x3 Jinv;
    chi_mesh::Matrix3x3 JTinv;
  };
  // Goes into
  /**Stores data for each face.*/
  struct FEface_data
  {
    std::vector<FEside_data3d> sides;
    chi_mesh::Vector3 normal;
  };
  // Goes int face_data

  /**Lowest level of mapping dof i.*/
  struct FEnodeSideMap
  {
    int index = -1;
    bool part_of_face = false;
  };
  // Goes into
  /**Intermediate level of mapping.*/
  struct FEnodeFaceMap
  {
    std::vector<FEnodeSideMap> side_map;
  };
  // Goes into
  /**Node map per face.*/
  struct FEnodeMap
  {
    std::vector<FEnodeFaceMap> face_map;
  };
  // Goes into node_maps
  //  node n
  //  face f
  //  side s
  //  node_maps[n]->face_map[f]->side_map[s]

  std::vector<double> face_betaf_; ///< Face Beta-factor.
  double alphac_;                  ///< Cell alpha-factor.

  std::vector<FEface_data> face_data_;    ///< Holds determinants and data tet-by-tet.
  std::vector<FEnodeMap> node_side_maps_; ///< Maps nodes to side tets.

  const QuadratureTetrahedron& volume_quadrature_;
  const QuadratureTriangle& surface_quadrature_;
};

} // namespace chi_math::cell_mapping
