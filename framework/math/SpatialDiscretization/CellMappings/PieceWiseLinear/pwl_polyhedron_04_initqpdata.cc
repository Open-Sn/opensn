#include "PieceWiseLinearPolyhedronMapping.h"

#include "math/SpatialDiscretization/FiniteElement/QuadraturePointData.h"

namespace chi_math::cell_mapping
{

finite_element::VolumetricQuadraturePointData
PieceWiseLinearPolyhedronMapping::MakeVolumetricQuadraturePointData() const
{
  //=================================== Determine number of internal qpoints
  size_t num_tets = 0;
  for (auto& face : face_data_)
    num_tets += face.sides.size();

  size_t num_vol_qpoints = volume_quadrature_.qpoints_.size();
  size_t ttl_num_vol_qpoints = num_tets * num_vol_qpoints;

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

    for (size_t f = 0; f < face_data_.size(); f++)
    {
      for (size_t s = 0; s < face_data_[f].sides.size(); s++)
      {
        for (const auto& qpoint : volume_quadrature_.qpoints_)
        {
          node_shape_value.push_back(FaceSideShape(f, s, i, qpoint));
          node_shape_grad.emplace_back(FaceSideGradShape_x(f, s, i),  // x
                                       FaceSideGradShape_y(f, s, i),  // y
                                       FaceSideGradShape_z(f, s, i)); // z
        }                                                             // for qp
      }                                                               // for side
    }                                                                 // for face

    V_shape_value.push_back(node_shape_value);
    V_shape_grad.push_back(node_shape_grad);
  } // for i

  V_JxW.reserve(ttl_num_vol_qpoints);
  V_qpoints_xyz.reserve(ttl_num_vol_qpoints);
  for (const auto& face : face_data_)
  {
    for (const auto& side : face.sides)
    {
      for (size_t qp = 0; qp < num_vol_qpoints; ++qp)
      {
        const auto w = volume_quadrature_.weights_[qp];
        V_JxW.push_back(side.detJ * w);

        const auto& qp_xyz_tilde = volume_quadrature_.qpoints_[qp];
        V_qpoints_xyz.push_back(side.v0 + side.J * qp_xyz_tilde);
      } // for qp
    }   // for side
  }     // for face

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
PieceWiseLinearPolyhedronMapping::MakeSurfaceQuadraturePointData(size_t face_index) const
{
  const bool ON_SURFACE = true;

  //=================================== Init surface quadrature
  size_t num_srf_qpoints = surface_quadrature_.qpoints_.size();

  unsigned int f = face_index;
  //=================================== Declare necessary vars
  std::vector<unsigned int> F_quadrature_point_indices;
  VecVec3 F_qpoints_xyz;
  std::vector<VecDbl> F_shape_value;
  std::vector<VecVec3> F_shape_grad;
  VecDbl F_JxW;
  VecVec3 F_normals;
  size_t F_num_nodes;

  size_t num_tris = face_data_[f].sides.size();
  size_t ttl_num_face_qpoints = num_tris * num_srf_qpoints;

  F_quadrature_point_indices.reserve(ttl_num_face_qpoints);
  for (unsigned int qp = 0; qp < ttl_num_face_qpoints; ++qp)
    F_quadrature_point_indices.push_back(qp);

  F_normals.reserve(ttl_num_face_qpoints);
  for (size_t qp = 0; qp < ttl_num_face_qpoints; ++qp)
    F_normals.push_back(face_data_[f].normal);

  F_shape_value.reserve(num_nodes_);
  F_shape_grad.reserve(num_nodes_);
  for (size_t i = 0; i < num_nodes_; i++)
  {
    VecDbl node_shape_value;
    VecVec3 node_shape_grad;

    node_shape_value.reserve(ttl_num_face_qpoints);
    node_shape_grad.reserve(ttl_num_face_qpoints);

    for (size_t s = 0; s < face_data_[f].sides.size(); s++)
    {
      for (const auto& qpoint : surface_quadrature_.qpoints_)
      {
        node_shape_value.push_back(FaceSideShape(f, s, i, qpoint, ON_SURFACE));
        node_shape_grad.emplace_back(FaceSideGradShape_x(f, s, i),  // x
                                     FaceSideGradShape_y(f, s, i),  // y
                                     FaceSideGradShape_z(f, s, i)); // z
      }                                                             // for qp
    }                                                               // for s
    F_shape_value.push_back(node_shape_value);
    F_shape_grad.push_back(node_shape_grad);
  } // for i

  F_JxW.reserve(ttl_num_face_qpoints);
  F_qpoints_xyz.reserve(ttl_num_face_qpoints);
  for (const auto& side : face_data_[f].sides)
    for (size_t qp = 0; qp < num_srf_qpoints; ++qp)
    {
      const auto w = surface_quadrature_.weights_[qp];
      F_JxW.push_back(side.detJ_surf * w);

      const auto& qp_xyz_tilde = surface_quadrature_.qpoints_[qp];
      F_qpoints_xyz.push_back(side.v0 + side.J * qp_xyz_tilde);
    }

  F_num_nodes = face_data_[f].sides.size();

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
