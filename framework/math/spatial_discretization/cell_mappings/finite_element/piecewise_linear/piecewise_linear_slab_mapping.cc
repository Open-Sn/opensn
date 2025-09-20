// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/spatial_discretization/cell_mappings/finite_element/piecewise_linear/piecewise_linear_slab_mapping.h"
#include "framework/math/spatial_discretization/finite_element/finite_element_data.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"

namespace opensn
{

PieceWiseLinearSlabMapping::PieceWiseLinearSlabMapping(
  const Cell& slab_cell,
  const std::shared_ptr<MeshContinuum> ref_grid,
  const LineQuadrature& volume_quadrature)
  : PieceWiseLinearBaseMapping(ref_grid, slab_cell, 2, MakeFaceNodeMapping(slab_cell)),
    v0i_(slab_cell.vertex_ids[0]),
    v1i_(slab_cell.vertex_ids[1]),
    volume_quadrature_(volume_quadrature)
{
  v0_ = grid_->vertices[v0i_];
  const auto& v1 = grid_->vertices[v1i_];

  Vector3 v01 = v1 - v0_;
  h_ = v01.Norm();

  normals_[0] = slab_cell.faces[0].normal;
  normals_[1] = slab_cell.faces[1].normal;
}

double
PieceWiseLinearSlabMapping::SlabShape(uint32_t index,
                                      const Vector3& qpoint,
                                      bool on_surface,
                                      uint32_t edge) const
{
  double xi = 0.0;
  if (not on_surface)
    xi = qpoint.x;
  else
    xi = static_cast<double>(edge);

  double value = 0.0;
  if (index == 0)
    value = 1.0 - xi;
  else if (index == 1)
    value = xi;

  return value;
}

double
PieceWiseLinearSlabMapping::SlabGradShape(uint32_t index) const
{
  double value = 0.0;

  if (index == 0)
    value = -1.0 / h_;
  else
    value = 1.0 / h_;

  return value;
}

double
PieceWiseLinearSlabMapping::ShapeValue(const int i, const Vector3& xyz) const
{
  const auto& p0 = grid_->vertices[v0i_];
  const auto& p1 = grid_->vertices[v1i_];
  Vector3 xyz_ref = xyz - p0;

  Vector3 v01 = p1 - p0;

  double xi = v01.Dot(xyz_ref) / v01.Norm() / h_;

  if ((xi >= -1.0e-6) and (xi <= 1.0 + 1.0e-6))
  {
    if (i == 0)
      return 1.0 - xi;
    else
      return xi;
  } // if in cell

  return 0.0;
}

void
PieceWiseLinearSlabMapping::ShapeValues(const Vector3& xyz, Vector<double>& shape_values) const
{
  shape_values.Resize(num_nodes_, 0.0);
  const auto& p0 = grid_->vertices[v0i_];
  const auto& p1 = grid_->vertices[v1i_];
  Vector3 xyz_ref = xyz - p0;

  Vector3 v01 = p1 - p0;

  double xi = v01.Dot(xyz_ref) / v01.Norm() / h_;

  if ((xi >= -1.0e-6) and (xi <= 1.0 + 1.0e-6))
  {
    for (size_t i = 0; i < num_nodes_; ++i)
    {
      if (i == 0)
        shape_values(i) = 1.0 - xi;
      else
        shape_values(i) = xi;
    } // for dof

    return;
  } // if in cell
}

Vector3
PieceWiseLinearSlabMapping::GradShapeValue(const int i, const Vector3& xyz) const
{
  if (i == 0)
    return Vector3(0.0, 0.0, -1.0 / h_);
  else
    return Vector3(0.0, 0.0, 1.0 / h_);
}

void
PieceWiseLinearSlabMapping::GradShapeValues(const Vector3& xyz,
                                            std::vector<Vector3>& gradshape_values) const
{
  gradshape_values.clear();
  gradshape_values.emplace_back(GradShapeValue(0, xyz));
  gradshape_values.emplace_back(GradShapeValue(1, xyz));
}

VolumetricFiniteElementData
PieceWiseLinearSlabMapping::MakeVolumetricFiniteElementData() const
{
  // Determine number of internal qpoints
  size_t ttl_num_vol_qpoints = volume_quadrature_.qpoints.size();

  // Declare necessary vars
  std::vector<unsigned int> V_quadrature_point_indices;
  std::vector<Vector3> V_qpoints_xyz;
  std::vector<std::vector<double>> V_shape_value;
  std::vector<std::vector<Vector3>> V_shape_grad;
  std::vector<double> V_JxW;
  size_t V_num_nodes;

  // Init volumetric quadrature
  V_quadrature_point_indices.reserve(ttl_num_vol_qpoints);
  for (unsigned int qp = 0; qp < ttl_num_vol_qpoints; ++qp)
    V_quadrature_point_indices.push_back(qp);

  V_shape_value.reserve(num_nodes_);
  V_shape_grad.reserve(num_nodes_);
  for (size_t i = 0; i < num_nodes_; ++i)
  {
    std::vector<double> node_shape_value;
    std::vector<Vector3> node_shape_grad;

    node_shape_value.reserve(ttl_num_vol_qpoints);
    node_shape_grad.reserve(ttl_num_vol_qpoints);

    for (const auto& qpoint : volume_quadrature_.qpoints)
    {
      node_shape_value.push_back(SlabShape(i, qpoint));
      node_shape_grad.emplace_back(0.0, 0.0, SlabGradShape(i));
    } // for qp

    V_shape_value.push_back(node_shape_value);
    V_shape_grad.push_back(node_shape_grad);
  } // for i

  V_JxW.reserve(ttl_num_vol_qpoints);
  V_qpoints_xyz.reserve(ttl_num_vol_qpoints);
  const double J = h_;
  for (size_t qp = 0; qp < ttl_num_vol_qpoints; ++qp)
  {
    const double w = volume_quadrature_.weights[qp];
    V_JxW.push_back(J * w);

    const double qp_xyz_tilde = volume_quadrature_.qpoints[qp][0];
    V_qpoints_xyz.push_back(v0_ + J * Vector3(0.0, 0.0, qp_xyz_tilde));
  } // for qp

  V_num_nodes = num_nodes_;

  return VolumetricFiniteElementData(V_quadrature_point_indices,
                                     V_qpoints_xyz,
                                     V_shape_value,
                                     V_shape_grad,
                                     V_JxW,
                                     face_node_mappings_,
                                     V_num_nodes);
}

SurfaceFiniteElementData
PieceWiseLinearSlabMapping::MakeSurfaceFiniteElementData(size_t face_index) const
{
  const bool ON_SURFACE = true;

  // Init surface quadrature
  size_t num_srf_qpoints = 1;

  const unsigned int f = face_index;

  // Declare necessary vars
  std::vector<unsigned int> F_quadrature_point_indices;
  std::vector<Vector3> F_qpoints_xyz;
  std::vector<std::vector<double>> F_shape_value;
  std::vector<std::vector<Vector3>> F_shape_grad;
  std::vector<double> F_JxW;
  std::vector<Vector3> F_normals;
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
  for (size_t i = 0; i < num_nodes_; ++i)
  {
    std::vector<double> node_shape_value;
    std::vector<Vector3> node_shape_grad;

    node_shape_value.reserve(ttl_num_face_qpoints);
    node_shape_grad.reserve(ttl_num_face_qpoints);

    for (const auto& qpoint : {Vector3(0.0, 0.0, 0.0)})
    {
      node_shape_value.push_back(SlabShape(i, qpoint, ON_SURFACE, f));
      node_shape_grad.emplace_back(0.0, 0.0, SlabGradShape(i));
    } // for qp
    F_shape_value.push_back(node_shape_value);
    F_shape_grad.push_back(node_shape_grad);
  } // for i

  F_JxW.reserve(ttl_num_face_qpoints);
  F_qpoints_xyz.reserve(ttl_num_face_qpoints);
  const double JxW = 1.0;
  for (size_t qp = 0; qp < num_srf_qpoints; ++qp)
  {
    F_JxW.push_back(JxW);

    F_qpoints_xyz.emplace_back(0.0, 0.0, f);
  }

  F_num_nodes = 1;

  return SurfaceFiniteElementData(F_quadrature_point_indices,
                                  F_qpoints_xyz,
                                  F_shape_value,
                                  F_shape_grad,
                                  F_JxW,
                                  F_normals,
                                  face_node_mappings_,
                                  F_num_nodes);
}

} // namespace opensn
