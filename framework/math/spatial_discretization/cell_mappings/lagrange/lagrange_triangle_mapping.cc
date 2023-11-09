#include "framework/math/spatial_discretization/cell_mappings/lagrange/lagrange_triangle_mapping.h"

#include "framework/mesh/cell/cell.h"

#include "framework/logging/log_exceptions.h"

namespace chi_math::cell_mapping
{

LagrangeTriangleMapping::LagrangeTriangleMapping(const chi_mesh::MeshContinuum& grid,
                                                 const chi_mesh::Cell& cell,
                                                 const Quadrature& volume_quadrature,
                                                 const Quadrature& surface_quadrature)
  : LagrangeBaseMapping(
      grid, cell, 3, MakeFaceNodeMapping(cell), volume_quadrature, surface_quadrature)
{
}

double
LagrangeTriangleMapping::RefShape(uint32_t i, const Vec3& qpoint) const
{
  // clang-format off
  if (i == 0) return 1.0 - qpoint.x - qpoint.y;
  if (i == 1) return qpoint.x;
  if (i == 2) return qpoint.y;
  // clang-format on

  ChiLogicalError("Invalid shapefunction index " + std::to_string(i));
}

chi_mesh::Vector3
LagrangeTriangleMapping::RefGradShape(uint32_t i, const Vec3& qpoint) const
{
  ChiLogicalErrorIf(i >= 3, "Invalid shapefunction index " + std::to_string(i));

  // clang-format off
  if (i == 0) return Vec3(-1.0, -1.0, 0.0);
  if (i == 1) return Vec3( 1.0,  0.0, 0.0);
  if (i == 2) return Vec3( 0.0,  1.0, 0.0);
  // clang-format on

  ChiLogicalError("Invalid shapefunction index " + std::to_string(i));
}

LagrangeBaseMapping::MatDbl
LagrangeTriangleMapping::RefJacobian(const Vec3& qpoint) const
{
  MatDbl J(3, VecDbl(3, 0.0));
  for (size_t i = 0; i < num_nodes_; ++i)
  {
    const Vec3 grad_shape_i = RefGradShape(i, qpoint);
    const auto& node_i = node_locations_[i];
    const double x_i = node_i.x;
    const double y_i = node_i.y;

    // Loops unrolled for performance
    J[0][0] += grad_shape_i.x * x_i;
    J[0][1] += grad_shape_i.y * x_i;
    J[1][0] += grad_shape_i.x * y_i;
    J[1][1] += grad_shape_i.y * y_i;
  }
  J[2][2] = 1.0;

  return J;
}

std::pair<double, LagrangeBaseMapping::Vec3>
LagrangeTriangleMapping::RefFaceJacobianDeterminantAndNormal(size_t face_index,
                                                             const Vec3& qpoint_face) const
{
  // x = sum_i N_i x_i
  const auto& x0 = node_locations_[face_node_mappings_[face_index][0]];
  const auto& x1 = node_locations_[face_node_mappings_[face_index][1]];
  Vec3 dx_dxbar;
  for (size_t i = 0; i < 2; ++i)
  {
    if (i == 0) dx_dxbar += -0.5 * x0;
    if (i == 1) dx_dxbar += 0.5 * x1;
  }

  const auto cross = dx_dxbar.Cross(Vec3(0.0, 0.0, 1.0));

  return {dx_dxbar.Norm(), cross.Normalized()};
}

LagrangeBaseMapping::Vec3
LagrangeTriangleMapping::FaceToElementQPointConversion(size_t face_index,
                                                       const Vec3& qpoint_face) const
{
  // The quadrature for the face of a triangle is a line-quadrature, thus
  // only x-component has a value.
  const double x = 0.5 * (qpoint_face.x + 1.0);
  // clang-format off
  if (face_index == 0) return Vec3(  x,     0.0, 0.0);
  if (face_index == 1) return Vec3(  x, 1.0 - x, 0.0);
  if (face_index == 2) return Vec3(0.0,       x, 0.0);
  // clang-format on

  ChiLogicalError("Invalid face index " + std::to_string(face_index));
}

} // namespace chi_math::cell_mapping
