#include "framework/math/SpatialDiscretization/CellMappings/PieceWiseLinear/PieceWiseLinearSlabMapping.h"
#include "framework/mesh/MeshContinuum/chi_meshcontinuum.h"
#include "framework/math/SpatialDiscretization/CellMappings/PieceWiseLinearBaseMapping.h"
#include "framework/math/SpatialDiscretization/FiniteElement/QuadraturePointData.h"

namespace chi_math::cell_mapping
{

PieceWiseLinearSlabMapping::PieceWiseLinearSlabMapping(
  const chi_mesh::Cell& slab_cell,
  const chi_mesh::MeshContinuum& ref_grid,
  const chi_math::QuadratureLine& volume_quadrature)
  : PieceWiseLinearBaseMapping(ref_grid,
                               slab_cell,
                               2, // num_nodes
                               MakeFaceNodeMapping(slab_cell)),
    volume_quadrature_(volume_quadrature)
{
  v0i_ = slab_cell.vertex_ids_[0];
  v1i_ = slab_cell.vertex_ids_[1];
  v0_ = ref_grid_.vertices[v0i_];
  const auto& v1 = ref_grid_.vertices[v1i_];

  chi_mesh::Vector3 v01 = v1 - v0_;
  h_ = v01.Norm();

  normals_[0] = slab_cell.faces_[0].normal_;
  normals_[1] = slab_cell.faces_[1].normal_;
}

double
PieceWiseLinearSlabMapping::SlabShape(uint32_t index,
                                      const chi_mesh::Vector3& qpoint,
                                      bool on_surface,
                                      uint32_t edge) const
{
  double xi = 0.0;
  if (!on_surface) xi = qpoint.x;
  else
    xi = static_cast<double>(edge);

  double value = 0.0;
  if (index == 0) value = 1.0 - xi;
  else if (index == 1)
    value = xi;

  return value;
}

double
PieceWiseLinearSlabMapping::SlabGradShape(uint32_t index) const
{
  double value = 0.0;

  if (index == 0) value = -1.0 / h_;
  else
    value = 1.0 / h_;

  return value;
}

double
PieceWiseLinearSlabMapping::ShapeValue(const int i, const chi_mesh::Vector3& xyz) const
{
  const auto& p0 = ref_grid_.vertices[v0i_];
  const auto& p1 = ref_grid_.vertices[v1i_];
  chi_mesh::Vector3 xyz_ref = xyz - p0;

  chi_mesh::Vector3 v01 = p1 - p0;

  double xi = v01.Dot(xyz_ref) / v01.Norm() / h_;

  if ((xi >= -1.0e-6) and (xi <= 1.0 + 1.0e-6))
  {
    if (i == 0) return 1.0 - xi;
    else
      return xi;
  } // if in cell

  return 0.0;
}

void
PieceWiseLinearSlabMapping::ShapeValues(const chi_mesh::Vector3& xyz,
                                        std::vector<double>& shape_values) const
{
  shape_values.resize(num_nodes_, 0.0);
  const auto& p0 = ref_grid_.vertices[v0i_];
  const auto& p1 = ref_grid_.vertices[v1i_];
  chi_mesh::Vector3 xyz_ref = xyz - p0;

  chi_mesh::Vector3 v01 = p1 - p0;

  double xi = v01.Dot(xyz_ref) / v01.Norm() / h_;

  if ((xi >= -1.0e-6) and (xi <= 1.0 + 1.0e-6))
  {
    for (int i = 0; i < num_nodes_; i++)
    {
      if (i == 0) shape_values[i] = 1.0 - xi;
      else
        shape_values[i] = xi;
    } // for dof

    return;
  } // if in cell
}

chi_mesh::Vector3
PieceWiseLinearSlabMapping::GradShapeValue(const int i, const chi_mesh::Vector3& xyz) const
{
  if (i == 0) return chi_mesh::Vector3(0.0, 0.0, -1.0 / h_);
  else
    return chi_mesh::Vector3(0.0, 0.0, 1.0 / h_);
}

void
PieceWiseLinearSlabMapping::GradShapeValues(const chi_mesh::Vector3& xyz,
                                            std::vector<chi_mesh::Vector3>& gradshape_values) const
{
  gradshape_values.clear();
  gradshape_values.emplace_back(GradShapeValue(0, xyz));
  gradshape_values.emplace_back(GradShapeValue(1, xyz));
}

finite_element::VolumetricQuadraturePointData
PieceWiseLinearSlabMapping::MakeVolumetricQuadraturePointData() const
{
  //=================================== Determine number of internal qpoints
  size_t ttl_num_vol_qpoints = volume_quadrature_.qpoints_.size();

  //=================================== Declare necessary vars
  std::vector<unsigned int> V_quadrature_point_indices;
  VecVec3 V_qpoints_xyz;
  std::vector<VecDbl> V_shape_value;
  std::vector<VecVec3> V_shape_grad;
  VecDbl V_JxW;
  size_t V_num_nodes;

  //=================================== Init volumetric quadrature
  V_quadrature_point_indices.reserve(ttl_num_vol_qpoints);
  for (unsigned int qp = 0; qp < ttl_num_vol_qpoints; ++qp)
    V_quadrature_point_indices.push_back(qp);

  V_shape_value.reserve(num_nodes_);
  V_shape_grad.reserve(num_nodes_);
  for (size_t i = 0; i < num_nodes_; i++)
  {
    VecDbl node_shape_value;
    VecVec3 node_shape_grad;

    node_shape_value.reserve(ttl_num_vol_qpoints);
    node_shape_grad.reserve(ttl_num_vol_qpoints);

    for (const auto& qpoint : volume_quadrature_.qpoints_)
    {
      node_shape_value.push_back(SlabShape(i, qpoint));
      node_shape_grad.emplace_back(0.0,               // x
                                   0.0,               // y
                                   SlabGradShape(i)); // z
    }                                                 // for qp

    V_shape_value.push_back(node_shape_value);
    V_shape_grad.push_back(node_shape_grad);
  } // for i

  V_JxW.reserve(ttl_num_vol_qpoints);
  V_qpoints_xyz.reserve(ttl_num_vol_qpoints);
  const double J = h_;
  for (size_t qp = 0; qp < ttl_num_vol_qpoints; ++qp)
  {
    const double w = volume_quadrature_.weights_[qp];
    V_JxW.push_back(J * w);

    const double qp_xyz_tilde = volume_quadrature_.qpoints_[qp][0];
    V_qpoints_xyz.push_back(v0_ + J * chi_mesh::Vector3(0.0, 0.0, qp_xyz_tilde));
  } // for qp

  V_num_nodes = num_nodes_;

  return finite_element::VolumetricQuadraturePointData(V_quadrature_point_indices,
                                                       V_qpoints_xyz,
                                                       V_shape_value,
                                                       V_shape_grad,
                                                       V_JxW,
                                                       face_node_mappings_,
                                                       V_num_nodes);
}

finite_element::SurfaceQuadraturePointData
PieceWiseLinearSlabMapping::MakeSurfaceQuadraturePointData(size_t face_index) const
{
  const bool ON_SURFACE = true;

  //=================================== Init surface quadrature
  size_t num_srf_qpoints = 1;

  const unsigned int f = face_index;

  //=================================== Declare necessary vars
  std::vector<unsigned int> F_quadrature_point_indices;
  VecVec3 F_qpoints_xyz;
  std::vector<VecDbl> F_shape_value;
  std::vector<VecVec3> F_shape_grad;
  VecDbl F_JxW;
  VecVec3 F_normals;
  size_t F_num_nodes;

  size_t ttl_num_face_qpoints = num_srf_qpoints;

  F_quadrature_point_indices.reserve(ttl_num_face_qpoints);
  for (unsigned int qp = 0; qp < ttl_num_face_qpoints; ++qp)
    F_quadrature_point_indices.push_back(qp);

  F_normals.reserve(ttl_num_face_qpoints);
  for (size_t qp = 0; qp < ttl_num_face_qpoints; ++qp)
    F_normals.push_back(normals_[f]);

  F_shape_value.reserve(num_nodes_);
  F_shape_grad.reserve(num_nodes_);
  for (size_t i = 0; i < num_nodes_; i++)
  {
    VecDbl node_shape_value;
    VecVec3 node_shape_grad;

    node_shape_value.reserve(ttl_num_face_qpoints);
    node_shape_grad.reserve(ttl_num_face_qpoints);

    for (const auto& qpoint : {chi_mesh::Vector3(0.0, 0.0, 0.0)})
    {
      node_shape_value.push_back(SlabShape(i, qpoint, ON_SURFACE, f));
      node_shape_grad.emplace_back(0.0,               // x
                                   0.0,               // y
                                   SlabGradShape(i)); // z
    }                                                 // for qp
    F_shape_value.push_back(node_shape_value);
    F_shape_grad.push_back(node_shape_grad);
  } // for i

  F_JxW.reserve(ttl_num_face_qpoints);
  F_qpoints_xyz.reserve(ttl_num_face_qpoints);
  const double JxW = 1.0;
  for (size_t qp = 0; qp < num_srf_qpoints; ++qp)
  {
    F_JxW.push_back(JxW);

    F_qpoints_xyz.push_back(chi_mesh::Vector3(0.0, 0.0, f));
  }

  F_num_nodes = 1;

  return finite_element::SurfaceQuadraturePointData(F_quadrature_point_indices,
                                                    F_qpoints_xyz,
                                                    F_shape_value,
                                                    F_shape_grad,
                                                    F_JxW,
                                                    F_normals,
                                                    face_node_mappings_,
                                                    F_num_nodes);
}

} // namespace chi_math::cell_mapping
