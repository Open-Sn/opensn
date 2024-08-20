// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/spatial_discretization/cell_mappings/finite_element/piecewise_linear/piecewise_linear_polygon_mapping.h"
#include "framework/math/spatial_discretization/finite_element/finite_element_data.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"

namespace opensn
{

PieceWiseLinearPolygonMapping::PieceWiseLinearPolygonMapping(
  const Cell& poly_cell,
  const MeshContinuum& ref_grid,
  const TriangleQuadrature& volume_quadrature,
  const LineQuadrature& surface_quadrature)
  : PieceWiseLinearBaseMapping(
      ref_grid, poly_cell, poly_cell.vertex_ids.size(), MakeFaceNodeMapping(poly_cell)),
    volume_quadrature_(volume_quadrature),
    surface_quadrature_(surface_quadrature)
{
  num_of_subtris_ = static_cast<int>(poly_cell.faces.size());
  beta_ = 1.0 / num_of_subtris_;

  // Get raw vertices
  vc_ = poly_cell.centroid;

  // Calculate legs and determinants
  for (int side = 0; side < num_of_subtris_; ++side)
  {
    const CellFace& face = poly_cell.faces[side];

    const auto& v0 = ref_grid_.vertices[face.vertex_ids[0]];
    const auto& v1 = ref_grid_.vertices[face.vertex_ids[1]];
    Vector3 v2 = vc_;

    Vector3 sidev01 = v1 - v0;
    Vector3 sidev02 = v2 - v0;

    double sidedetJ = ((sidev01.x) * (sidev02.y) - (sidev02.x) * (sidev01.y));

    FEside_data2d triangle_data;
    triangle_data.detJ = sidedetJ;
    triangle_data.detJ_surf = sidev01.Norm();

    triangle_data.v_index[0] = face.vertex_ids[0];
    triangle_data.v_index[1] = face.vertex_ids[1];

    triangle_data.v0 = v0;

    // Set Jacobian
    triangle_data.J.SetIJ(0, 0, sidev01.x);
    triangle_data.J.SetIJ(1, 0, sidev01.y);
    triangle_data.J.SetIJ(0, 1, sidev02.x);
    triangle_data.J.SetIJ(1, 1, sidev02.y);
    triangle_data.J.SetIJ(2, 2, 0.0);

    // Set Jacobian inverse
    triangle_data.Jinv.SetIJ(0, 0, sidev02.y / sidedetJ);
    triangle_data.Jinv.SetIJ(1, 0, -sidev01.y / sidedetJ);
    triangle_data.Jinv.SetIJ(0, 1, -sidev02.x / sidedetJ);
    triangle_data.Jinv.SetIJ(1, 1, sidev01.x / sidedetJ);
    triangle_data.Jinv.SetIJ(2, 2, 0.0);

    // Set Jacobian-Transpose inverse
    triangle_data.JTinv.SetIJ(0, 0, sidev02.y / sidedetJ);
    triangle_data.JTinv.SetIJ(1, 0, -sidev02.x / sidedetJ);
    triangle_data.JTinv.SetIJ(0, 1, -sidev01.y / sidedetJ);
    triangle_data.JTinv.SetIJ(1, 1, sidev01.x / sidedetJ);
    triangle_data.JTinv.SetIJ(2, 2, 0.0);

    // Set face normal
    triangle_data.normal = face.normal;

    sides_.push_back(triangle_data);
  }

  // Compute node to side mapping
  for (int v = 0; v < poly_cell.vertex_ids.size(); ++v)
  {
    const uint64_t vindex = poly_cell.vertex_ids[v];
    std::vector<int> side_mapping(num_of_subtris_);
    for (int side = 0; side < num_of_subtris_; ++side)
    {
      side_mapping[side] = -1;

      const CellFace& face = poly_cell.faces[side];
      if (face.vertex_ids[0] == vindex)
      {
        side_mapping[side] = 0;
      }
      if (face.vertex_ids[1] == vindex)
      {
        side_mapping[side] = 1;
      }
    }
    node_to_side_map_.push_back(side_mapping);
  }
}

double
PieceWiseLinearPolygonMapping::TriShape(uint32_t index, const Vector3& qpoint, bool on_surface)
{
  double xi;
  double eta;
  if (not on_surface)
  {
    xi = qpoint.x;
    eta = qpoint.y;
  }
  else
  {
    xi = qpoint.x;
    eta = 0.0;
  }

  double value = 0.0;
  if (index == 0)
    value = 1.0 - xi - eta;
  else if (index == 1)
    value = xi;
  else if (index == 2)
    value = eta;

  return value;
}

double
PieceWiseLinearPolygonMapping::SideShape(uint32_t side,
                                         uint32_t i,
                                         const Vector3& qpoint,
                                         bool on_surface) const
{
  int index = node_to_side_map_[i][side];
  double value = 0.0;
  if (index == 0 or index == 1)
    value = TriShape(index, qpoint, on_surface);

  value += beta_ * TriShape(2, qpoint, on_surface);

  return value;
}

double
PieceWiseLinearPolygonMapping::SideGradShape_x(uint32_t side, uint32_t i) const
{
  int index = node_to_side_map_[i][side];
  double value = 0;
  if (index == 0)
  {

    value = sides_[side].JTinv.GetIJ(0, 0) * -1.0 + sides_[side].JTinv.GetIJ(0, 1) * -1.0;
  }
  if (index == 1)
  {

    value = sides_[side].JTinv.GetIJ(0, 0) * 1.0 + sides_[side].JTinv.GetIJ(0, 1) * 0.0;
  }

  value += beta_ * (sides_[side].JTinv.GetIJ(0, 0) * 0.0 + sides_[side].JTinv.GetIJ(0, 1) * 1.0);

  return value;
}

double
PieceWiseLinearPolygonMapping::SideGradShape_y(uint32_t side, uint32_t i) const
{
  int index = node_to_side_map_[i][side];
  double value = 0;
  if (index == 0)
  {

    value = sides_[side].JTinv.GetIJ(1, 0) * -1.0 + sides_[side].JTinv.GetIJ(1, 1) * -1.0;
  }
  if (index == 1)
  {

    value = sides_[side].JTinv.GetIJ(1, 0) * 1.0 + sides_[side].JTinv.GetIJ(1, 1) * 0.0;
  }

  value += beta_ * (sides_[side].JTinv.GetIJ(1, 0) * 0.0 + sides_[side].JTinv.GetIJ(1, 1) * 1.0);

  return value;
}

double
PieceWiseLinearPolygonMapping::ShapeValue(const int i, const Vector3& xyz) const
{
  for (int s = 0; s < num_of_subtris_; ++s)
  {
    const auto& p0 = ref_grid_.vertices[sides_[s].v_index[0]];
    Vector3 xyz_ref = xyz - p0;

    Vector3 xi_eta_zeta = sides_[s].Jinv * xyz_ref;

    double xi = xi_eta_zeta.x;
    double eta = xi_eta_zeta.y;

    // Determine if inside tet
    if ((xi >= -1.0e-12) and (eta >= -1.0e-12) and ((xi + eta) <= (1.0 + 1.0e-12)))
    {
      int index = node_to_side_map_[i][s];
      double value = 0.0;

      if (index == 0)
      {
        value = 1.0 - xi - eta;
      }
      if (index == 1)
      {
        value = xi;
      }

      value += beta_ * eta;

      return value;
    }
  }

  return 0.0;
}

void
PieceWiseLinearPolygonMapping::ShapeValues(const Vector3& xyz, Vector<double>& shape_values) const
{
  shape_values.Resize(num_nodes_, 0.0);
  for (int s = 0; s < num_of_subtris_; ++s)
  {
    const auto& p0 = ref_grid_.vertices[sides_[s].v_index[0]];
    Vector3 xi_eta_zeta = sides_[s].Jinv * (xyz - p0);

    double xi = xi_eta_zeta.x;
    double eta = xi_eta_zeta.y;

    // Determine if inside tet
    if ((xi >= -1.0e-12) and (eta >= -1.0e-12) and ((xi + eta) <= (1.0 + 1.0e-12)))
    {
      for (int i = 0; i < num_nodes_; ++i)
      {
        int index = node_to_side_map_[i][s];
        double value = 0.0;

        if (index == 0)
        {
          value = 1.0 - xi - eta;
        }
        if (index == 1)
        {
          value = xi;
        }

        value += beta_ * eta;

        shape_values(i) = value;
      }
      return;
    } // if in triangle
  }   // for side
}

Vector3
PieceWiseLinearPolygonMapping::GradShapeValue(const int i, const Vector3& xyz) const
{
  Vector3 grad_r;
  Vector3 grad;

  for (int e = 0; e < num_of_subtris_; ++e)
  {
    const auto& p0 = ref_grid_.vertices[sides_[e].v_index[0]];
    Vector3 xyz_ref = xyz - p0;

    Vector3 xi_eta_zeta = sides_[e].Jinv * xyz_ref;

    double xi = xi_eta_zeta.x;
    double eta = xi_eta_zeta.y;

    if ((xi >= -1.0e-12) and (eta >= -1.0e-12) and ((xi + eta) <= (1.0 + 1.0e-12)))
    {
      int index = node_to_side_map_[i][e];

      if (index == 0)
      {
        grad_r.x += -1.0;
        grad_r.y += -1.0;
      }
      if (index == 1)
      {
        grad_r.x += 1.0;
        grad_r.y += 0.0;
      }

      grad_r.y += beta_ * 1.0;

      grad = sides_[e].JTinv * grad_r;

      return grad;
    }
  }

  return grad;
}

void
PieceWiseLinearPolygonMapping::GradShapeValues(const Vector3& xyz,
                                               std::vector<Vector3>& gradshape_values) const
{
  gradshape_values.clear();
  for (int i = 0; i < num_nodes_; ++i)
    gradshape_values.emplace_back(GradShapeValue(i, xyz));
}

VolumetricFiniteElementData
PieceWiseLinearPolygonMapping::MakeVolumetricFiniteElementData() const
{
  // Determine number of internal qpoints
  size_t num_tris = sides_.size();
  size_t num_vol_qpoints = volume_quadrature_.qpoints.size();
  size_t ttl_num_vol_qpoints = num_tris * num_vol_qpoints;

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

    for (size_t s = 0; s < sides_.size(); ++s)
    {
      for (const auto& qpoint : volume_quadrature_.qpoints)
      {
        node_shape_value.push_back(SideShape(s, i, qpoint));
        node_shape_grad.emplace_back(SideGradShape_x(s, i), SideGradShape_y(s, i), 0.0);
      } // for qp
    }   // for side

    V_shape_value.push_back(node_shape_value);
    V_shape_grad.push_back(node_shape_grad);
  } // for i

  V_JxW.reserve(ttl_num_vol_qpoints);
  V_qpoints_xyz.reserve(ttl_num_vol_qpoints);
  for (const auto& side : sides_)
  {
    for (size_t qp = 0; qp < num_vol_qpoints; ++qp)
    {
      const auto w = volume_quadrature_.weights[qp];
      V_JxW.push_back(side.detJ * w);

      const auto& qp_xyz_tilde = volume_quadrature_.qpoints[qp];
      V_qpoints_xyz.push_back(side.v0 + side.J * qp_xyz_tilde);
    } // for qp
  }   // for side

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
PieceWiseLinearPolygonMapping::MakeSurfaceFiniteElementData(size_t face_index) const
{
  const bool ON_SURFACE = true;

  // Init surface quadrature
  size_t num_srf_qpoints = surface_quadrature_.qpoints.size();

  unsigned int s = face_index;
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
    F_normals.push_back(sides_[s].normal);

  F_shape_value.reserve(num_nodes_);
  F_shape_grad.reserve(num_nodes_);
  for (size_t i = 0; i < num_nodes_; ++i)
  {
    std::vector<double> node_shape_value;
    std::vector<Vector3> node_shape_grad;

    node_shape_value.reserve(ttl_num_face_qpoints);
    node_shape_grad.reserve(ttl_num_face_qpoints);

    for (const auto& qpoint : surface_quadrature_.qpoints)
    {
      node_shape_value.push_back(SideShape(s, i, qpoint, ON_SURFACE));
      node_shape_grad.emplace_back(SideGradShape_x(s, i), SideGradShape_y(s, i), 0.0);
    } // for qp
    F_shape_value.push_back(node_shape_value);
    F_shape_grad.push_back(node_shape_grad);
  } // for i

  F_JxW.reserve(ttl_num_face_qpoints);
  F_qpoints_xyz.reserve(ttl_num_face_qpoints);
  for (size_t qp = 0; qp < num_srf_qpoints; ++qp)
  {
    const auto w = surface_quadrature_.weights[qp];
    F_JxW.push_back(sides_[s].detJ_surf * w);

    const auto& qp_xyz_tilde = surface_quadrature_.qpoints[qp];
    F_qpoints_xyz.push_back(sides_[s].v0 + sides_[s].J * qp_xyz_tilde);
  }

  F_num_nodes = 2;

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
