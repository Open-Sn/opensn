#include "QuadraturePointData.h"

namespace chi_math::finite_element
{
VolumetricQuadraturePointData::VolumetricQuadraturePointData() {}

VolumetricQuadraturePointData::VolumetricQuadraturePointData(
  std::vector<unsigned int> quadrature_point_indices,
  VecVec3 qpoints_xyz,
  std::vector<VecDbl> shape_value,
  std::vector<VecVec3> shape_grad,
  VecDbl JxW,
  std::vector<std::vector<int>> face_dof_mappings,
  size_t num_nodes)
  : quadrature_point_indices_(std::move(quadrature_point_indices)),
    qpoints_xyz_(std::move(qpoints_xyz)),
    shape_value_(std::move(shape_value)),
    shape_grad_(std::move(shape_grad)),
    JxW_(std::move(JxW)),
    face_dof_mappings_(std::move(face_dof_mappings)),
    num_nodes_(num_nodes)
{
}

const std::vector<unsigned int>&
VolumetricQuadraturePointData::QuadraturePointIndices() const
{
  return quadrature_point_indices_;
}
chi_mesh::Vector3
VolumetricQuadraturePointData::QPointXYZ(unsigned int qp) const
{
  return qpoints_xyz_.at(qp);
}
double VolumetricQuadraturePointData::ShapeValue(unsigned int i,
                                               unsigned int qp) const
{
  auto& qp_data = shape_value_.at(i);
  return qp_data.at(qp);
}
chi_mesh::Vector3
VolumetricQuadraturePointData::ShapeGrad(unsigned int i,
                                                         unsigned int qp) const
{
  auto& qp_data = shape_grad_.at(i);
  return qp_data.at(qp);
}

const VecVec3& VolumetricQuadraturePointData::QPointsXYZ() const
{
  return qpoints_xyz_;
}

const std::vector<VecDbl>& VolumetricQuadraturePointData::ShapeValues() const
{
  return shape_value_;
}
const std::vector<VecVec3>&
VolumetricQuadraturePointData::ShapeGradValues() const
{
  return shape_grad_;
}
const std::vector<double>& VolumetricQuadraturePointData::JxW_Values() const
{
  return JxW_;
}
double VolumetricQuadraturePointData::JxW(unsigned int qp) const
{
  return JxW_.at(qp);
}
int VolumetricQuadraturePointData::FaceDofMapping(size_t face,
                                                size_t face_node_index) const
{
  auto& face_data = face_dof_mappings_.at(face);
  return face_data.at(face_node_index);
}
size_t VolumetricQuadraturePointData::NumNodes() const { return num_nodes_; }

SurfaceQuadraturePointData::SurfaceQuadraturePointData()
{}

SurfaceQuadraturePointData::SurfaceQuadraturePointData(
  std::vector<unsigned int> quadrature_point_indices,
  VecVec3 qpoints_xyz,
  std::vector<VecDbl> shape_value,
  std::vector<VecVec3> shape_grad,
  VecDbl JxW,
  VecVec3 normals,
  std::vector<std::vector<int>> face_dof_mappings,
  size_t num_nodes)
  : VolumetricQuadraturePointData(std::move(quadrature_point_indices),
                                std::move(qpoints_xyz),
                                std::move(shape_value),
                                std::move(shape_grad),
                                std::move(JxW),
                                std::move(face_dof_mappings),
                                num_nodes),
    normals_(std::move(normals))
{
}

chi_mesh::Vector3 SurfaceQuadraturePointData::Normal(unsigned int qp) const
{
  return normals_.at(qp);
}

const VecVec3& SurfaceQuadraturePointData::Normals() const { return normals_; }

} // namespace chi_math::finite_element
