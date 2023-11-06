#include "framework/math/SpatialDiscretization/CellMappings/Lagrange/LagrangeTetMapping.h"

#include "framework/mesh/Cell/cell.h"

#include "framework/logging/chi_log_exceptions.h"

namespace chi_math::cell_mapping
{

LagrangeTetMapping::LagrangeTetMapping(const chi_mesh::MeshContinuum& grid,
                                       const chi_mesh::Cell& cell,
                                       const Quadrature& volume_quadrature,
                                       const Quadrature& surface_quadrature)
  : LagrangeBaseMapping(
      grid, cell, 4, MakeFaceNodeMapping(cell), volume_quadrature, surface_quadrature)
{
}

double
LagrangeTetMapping::RefShape(uint32_t i, const Vec3& qpoint) const
{
  const double x = qpoint.x;
  const double y = qpoint.y;
  const double z = qpoint.z;

  // clang-format off
  if (i == 0) return 1.0 - x - y - z;
  if (i == 1) return x;
  if (i == 2) return y;
  if (i == 3) return z;
  // clang-format on

  ChiLogicalError("Invalid shapefunction index " + std::to_string(i));
}

chi_mesh::Vector3
LagrangeTetMapping::RefGradShape(uint32_t i, const Vec3& qpoint) const
{
  // clang-format off
  if (i == 0) return Vec3( -1.0, -1.0, -1.0);
  if (i == 1) return Vec3(  1.0,  0.0,  0.0);
  if (i == 2) return Vec3(  0.0,  1.0,  0.0);
  if (i == 3) return Vec3(  0.0,  0.0,  1.0);
  // clang-format on

  ChiLogicalError("Invalid shapefunction index " + std::to_string(i));
}

LagrangeBaseMapping::MatDbl
LagrangeTetMapping::RefJacobian(const Vec3& qpoint) const
{
  MatDbl J(3, VecDbl(3, 0.0));
  for (size_t i = 0; i < num_nodes_; ++i)
  {
    const Vec3 grad_shape_i = RefGradShape(i, qpoint);
    const auto& node_i = node_locations_[i];
    const double x_i = node_i.x;
    const double y_i = node_i.y;
    const double z_i = node_i.z;

    // Loops unrolled for performance
    J[0][0] += grad_shape_i.x * x_i;
    J[0][1] += grad_shape_i.y * x_i;
    J[0][2] += grad_shape_i.z * x_i;

    J[1][0] += grad_shape_i.x * y_i;
    J[1][1] += grad_shape_i.y * y_i;
    J[1][2] += grad_shape_i.z * y_i;

    J[2][0] += grad_shape_i.x * z_i;
    J[2][1] += grad_shape_i.y * z_i;
    J[2][2] += grad_shape_i.z * z_i;
  }

  return J;
}

std::pair<double, LagrangeBaseMapping::Vec3>
LagrangeTetMapping::RefFaceJacobianDeterminantAndNormal(size_t face_index,
                                                        const Vec3& qpoint_face) const
{
  const auto& x0 = node_locations_[face_node_mappings_[face_index][0]];
  const auto& x1 = node_locations_[face_node_mappings_[face_index][1]];
  const auto& x2 = node_locations_[face_node_mappings_[face_index][2]];

  Vec3 dx_dxbar;
  Vec3 dx_dybar;
  for (size_t i = 0; i < 3; ++i)
  {
    double dN_dxbar = 0.0, dN_dybar = 0.0;
    Vec3 xi;
    // clang-format off
    if (i == 0) {dN_dxbar = -1.0; dN_dybar=-1.0; xi = x0;}
    if (i == 1) {dN_dxbar =  1.0; dN_dybar= 0.0; xi = x1;}
    if (i == 2) {dN_dxbar =  0.0; dN_dybar= 1.0; xi = x2;}
    // clang-format on

    dx_dxbar += dN_dxbar * xi;
    dx_dybar += dN_dybar * xi;
  }

  const auto cross = dx_dxbar.Cross(dx_dybar);
  const double detJ = cross.Norm();

  return {detJ, cross / detJ};
}

LagrangeBaseMapping::Vec3
LagrangeTetMapping::FaceToElementQPointConversion(size_t face_index, const Vec3& qpoint_face) const
{
  const double x = qpoint_face.x;
  const double y = qpoint_face.y;

  // clang-format off
  if (face_index == 0/*Base*/)    return Vec3(x, y, 0);
  if (face_index == 1/*XZ*/)      return Vec3(x, 0, y);
  if (face_index == 2/*YZ*/)      return Vec3(0, x, y);
  if (face_index == 3/*Slanted*/) return Vec3(x, y, 1.0 - x - y);
  // clang-format on

  ChiLogicalError("Invalid face index " + std::to_string(face_index));
}

} // namespace chi_math::cell_mapping
