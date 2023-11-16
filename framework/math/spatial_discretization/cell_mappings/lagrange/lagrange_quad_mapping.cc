#include "framework/math/spatial_discretization/cell_mappings/lagrange/lagrange_quad_mapping.h"

#include "framework/mesh/cell/cell.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"

#include "framework/logging/log_exceptions.h"

namespace opensn
{

LagrangeQuadMapping::LagrangeQuadMapping(const MeshContinuum& grid,
                                         const Cell& cell,
                                         const Quadrature& volume_quadrature,
                                         const Quadrature& surface_quadrature)
  : LagrangeBaseMapping(
      grid, cell, 4, MakeFaceNodeMapping(cell), volume_quadrature, surface_quadrature)
{
}

double
LagrangeQuadMapping::RefShape(uint32_t i, const Vec3& qpoint) const
{
  ChiLogicalErrorIf(i >= 4, "Invalid shapefunction index " + std::to_string(i));

  double a, b;
  // clang-format off
  if (i == 0) {a = -1.0; b=-1.0;}
  if (i == 1) {a =  1.0; b=-1.0;}
  if (i == 2) {a =  1.0; b= 1.0;}
  if (i == 3) {a = -1.0; b= 1.0;}
  // clang-format on

  return 0.25 * (1.0 + a * qpoint.x) * (1.0 + b * qpoint.y);
}

Vector3
LagrangeQuadMapping::RefGradShape(uint32_t i, const Vec3& qpoint) const
{
  ChiLogicalErrorIf(i >= 4, "Invalid shapefunction index " + std::to_string(i));

  double a, b;
  // clang-format off
  if (i == 0) {a = -1.0; b=-1.0;}
  if (i == 1) {a =  1.0; b=-1.0;}
  if (i == 2) {a =  1.0; b= 1.0;}
  if (i == 3) {a = -1.0; b= 1.0;}
  // clang-format on

  const double ab = a * b;
  return Vec3(0.25 * (a + ab * qpoint.y), 0.25 * (b + ab * qpoint.x), 0.0);
}

LagrangeBaseMapping::MatDbl
LagrangeQuadMapping::RefJacobian(const Vec3& qpoint) const
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
LagrangeQuadMapping::RefFaceJacobianDeterminantAndNormal(size_t face_index,
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
LagrangeQuadMapping::FaceToElementQPointConversion(size_t face_index, const Vec3& qpoint_face) const
{
  // The quadrature for the face of a quad is a line-quadrature, thus
  // only x-component has a value.
  const double x = qpoint_face.x + 1.0;
  if (face_index == 0) return Vec3(-1.0 + x, -1.0, 0.0);
  if (face_index == 1) return Vec3(1.0, -1.0 + x, 0.0);
  if (face_index == 2) return Vec3(1.0 - x, 1.0, 0.0);
  if (face_index == 3) return Vec3(-1.0, 1.0 - x, 0.0);

  ChiLogicalError("Invalid face index " + std::to_string(face_index));
}

} // namespace opensn
