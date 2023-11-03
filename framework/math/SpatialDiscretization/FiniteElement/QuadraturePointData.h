#pragma once

#include "math/chi_math.h"

namespace chi_math::finite_element
{
typedef std::vector<chi_mesh::Vector3> VecVec3;

// #############################################
/**Stored relevant quadrature point information
 * for volumetric integrals.*/
class VolumetricQuadraturePointData
{
public:
  VolumetricQuadraturePointData();
  VolumetricQuadraturePointData(std::vector<unsigned int> quadrature_point_indices,
                                VecVec3 qpoints_xyz,
                                std::vector<VecDbl> shape_value,
                                std::vector<VecVec3> shape_grad,
                                VecDbl JxW,
                                std::vector<std::vector<int>> face_dof_mappings,
                                size_t num_nodes);

  const std::vector<unsigned int>& QuadraturePointIndices() const;
  chi_mesh::Vector3 QPointXYZ(unsigned int qp) const;
  double ShapeValue(unsigned int i, unsigned int qp) const;
  chi_mesh::Vector3 ShapeGrad(unsigned int i, unsigned int qp) const;
  const VecVec3& QPointsXYZ() const;
  const std::vector<VecDbl>& ShapeValues() const;
  const std::vector<VecVec3>& ShapeGradValues() const;
  const std::vector<double>& JxW_Values() const;

  double JxW(unsigned int qp) const;
  int FaceDofMapping(size_t face, size_t face_node_index) const;
  size_t NumNodes() const;

protected:
  std::vector<unsigned int> quadrature_point_indices_; ///< qp index only
  VecVec3 qpoints_xyz_;                                ///< qp index only
  std::vector<VecDbl> shape_value_;                    ///< Node i, then qp
  std::vector<VecVec3> shape_grad_;                    ///< Node i, then qp
  VecDbl JxW_;                                         ///< qp index only
  std::vector<std::vector<int>> face_dof_mappings_;    ///< Face f,then fi
  size_t num_nodes_ = 0;
};

// #############################################
/**Stores relevant quadrature point information
 * for surface integrals.*/
class SurfaceQuadraturePointData : public VolumetricQuadraturePointData
{
public:
  SurfaceQuadraturePointData();
  SurfaceQuadraturePointData(std::vector<unsigned int> quadrature_point_indices,
                             VecVec3 qpoints_xyz,
                             std::vector<VecDbl> shape_value,
                             std::vector<VecVec3> shape_grad,
                             VecDbl JxW,
                             VecVec3 normals,
                             std::vector<std::vector<int>> face_dof_mappings,
                             size_t num_nodes);
  chi_mesh::Vector3 Normal(unsigned int qp) const;
  const VecVec3& Normals() const;

protected:
  VecVec3 normals_; ///< node i, then qp
};
} // namespace chi_math::finite_element
