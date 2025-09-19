// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/spatial_discretization/finite_element/finite_element_data.h"
#include <cassert>

namespace opensn
{

const std::vector<unsigned int>&
VolumetricFiniteElementData::GetQuadraturePointIndices() const
{
  return quadrature_point_indices_;
}

Vector3
VolumetricFiniteElementData::QPointXYZ(size_t qp) const
{
  assert(qp < qpoints_xyz_.size());
  return qpoints_xyz_[qp];
}

double
VolumetricFiniteElementData::ShapeValue(size_t i, size_t qp) const
{
  assert(i < shape_value_.size() and qp < shape_value_[i].size());
  return shape_value_[i][qp];
}

Vector3
VolumetricFiniteElementData::ShapeGrad(size_t i, size_t qp) const
{
  assert(i < shape_grad_.size() and qp < shape_grad_[i].size());
  return shape_grad_[i][qp];
}

const std::vector<Vector3>&
VolumetricFiniteElementData::GetQPointsXYZ() const
{
  return qpoints_xyz_;
}

const std::vector<std::vector<double>>&
VolumetricFiniteElementData::GetShapeValues() const
{
  return shape_value_;
}

const std::vector<std::vector<Vector3>>&
VolumetricFiniteElementData::GetShapeGradValues() const
{
  return shape_grad_;
}

const std::vector<double>&
VolumetricFiniteElementData::GetJxWValues() const
{
  return JxW_;
}

double
VolumetricFiniteElementData::JxW(size_t qp) const
{
  assert(qp < JxW_.size());
  return JxW_[qp];
}

int
VolumetricFiniteElementData::FaceDofMapping(size_t face, size_t face_node_index) const
{
  assert(face < face_dof_mappings_.size() and face_node_index < face_dof_mappings_[face].size());
  return face_dof_mappings_[face][face_node_index];
}

size_t
VolumetricFiniteElementData::GetNumNodes() const
{
  return num_nodes_;
}

Vector3
SurfaceFiniteElementData::Normal(size_t qp) const
{
  assert(qp < normals_.size());
  return normals_[qp];
}

const std::vector<Vector3>&
SurfaceFiniteElementData::GetNormals() const
{
  return normals_;
}

} // namespace opensn
