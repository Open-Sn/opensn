#include "framework/math/spatial_discretization/finite_element/finite_element_data.h"

namespace opensn
{

const std::vector<unsigned int>&
VolumetricFiniteElementData::QuadraturePointIndices() const
{
  return quadrature_point_indices_;
}

Vector3
VolumetricFiniteElementData::QPointXYZ(unsigned int qp) const
{
  return qpoints_xyz_.at(qp);
}

double
VolumetricFiniteElementData::ShapeValue(unsigned int i, unsigned int qp) const
{
  auto& qp_data = shape_value_.at(i);
  return qp_data.at(qp);
}

Vector3
VolumetricFiniteElementData::ShapeGrad(unsigned int i, unsigned int qp) const
{
  auto& qp_data = shape_grad_.at(i);
  return qp_data.at(qp);
}

const std::vector<Vector3>&
VolumetricFiniteElementData::QPointsXYZ() const
{
  return qpoints_xyz_;
}

const std::vector<VecDbl>&
VolumetricFiniteElementData::ShapeValues() const
{
  return shape_value_;
}

const std::vector<std::vector<Vector3>>&
VolumetricFiniteElementData::ShapeGradValues() const
{
  return shape_grad_;
}

const std::vector<double>&
VolumetricFiniteElementData::JxW_Values() const
{
  return JxW_;
}

double
VolumetricFiniteElementData::JxW(unsigned int qp) const
{
  return JxW_.at(qp);
}

int
VolumetricFiniteElementData::FaceDofMapping(size_t face, size_t face_node_index) const
{
  auto& face_data = face_dof_mappings_.at(face);
  return face_data.at(face_node_index);
}

size_t
VolumetricFiniteElementData::NumNodes() const
{
  return num_nodes_;
}

Vector3
SurfaceFiniteElementData::Normal(unsigned int qp) const
{
  return normals_.at(qp);
}

const std::vector<Vector3>&
SurfaceFiniteElementData::Normals() const
{
  return normals_;
}

} // namespace opensn
