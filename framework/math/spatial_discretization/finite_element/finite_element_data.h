// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/math.h"

namespace opensn
{
/**Stores finite element data for volumetric integrals.*/
class VolumetricFiniteElementData
{
public:
  VolumetricFiniteElementData(std::vector<unsigned int> quadrature_point_indices,
                              std::vector<Vector3> qpoints_xyz,
                              std::vector<std::vector<double>> shape_value,
                              std::vector<std::vector<Vector3>> shape_grad,
                              std::vector<double> JxW,
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

  const std::vector<unsigned int>& QuadraturePointIndices() const;

  Vector3 QPointXYZ(unsigned int qp) const;

  double ShapeValue(unsigned int i, unsigned int qp) const;

  Vector3 ShapeGrad(unsigned int i, unsigned int qp) const;

  const std::vector<Vector3>& QPointsXYZ() const;

  const std::vector<std::vector<double>>& ShapeValues() const;

  const std::vector<std::vector<Vector3>>& ShapeGradValues() const;

  const std::vector<double>& JxW_Values() const;

  double JxW(unsigned int qp) const;

  int FaceDofMapping(size_t face, size_t face_node_index) const;

  size_t NumNodes() const;

private:
  std::vector<unsigned int> quadrature_point_indices_; ///< qp index only
  std::vector<Vector3> qpoints_xyz_;                   ///< qp index only
  std::vector<std::vector<double>> shape_value_;       ///< Node i, then qp
  std::vector<std::vector<Vector3>> shape_grad_;       ///< Node i, then qp
  std::vector<double> JxW_;                            ///< qp index only
  std::vector<std::vector<int>> face_dof_mappings_;    ///< Face f,then fi
  size_t num_nodes_;
};

/**Stores finite element information for surface integrals.*/
class SurfaceFiniteElementData : public VolumetricFiniteElementData
{
public:
  SurfaceFiniteElementData(std::vector<unsigned int> quadrature_point_indices,
                           std::vector<Vector3> qpoints_xyz,
                           std::vector<std::vector<double>> shape_value,
                           std::vector<std::vector<Vector3>> shape_grad,
                           std::vector<double> JxW,
                           std::vector<Vector3> normals,
                           std::vector<std::vector<int>> face_dof_mappings,
                           size_t num_nodes)
    : VolumetricFiniteElementData(std::move(quadrature_point_indices),
                                  std::move(qpoints_xyz),
                                  std::move(shape_value),
                                  std::move(shape_grad),
                                  std::move(JxW),
                                  std::move(face_dof_mappings),
                                  num_nodes),
      normals_(std::move(normals))
  {
  }

  Vector3 Normal(unsigned int qp) const;

  const std::vector<Vector3>& Normals() const;

private:
  std::vector<Vector3> normals_; ///< node i, then qp
};

} // namespace opensn
