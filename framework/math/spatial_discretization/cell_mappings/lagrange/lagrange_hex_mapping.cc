#include "framework/math/spatial_discretization/cell_mappings/lagrange/lagrange_hex_mapping.h"

#include "framework/mesh/cell/cell.h"

#include "framework/logging/log.h"

namespace opensn
{

LagrangeHexMapping::LagrangeHexMapping(const MeshContinuum& grid,
                                       const Cell& cell,
                                       const Quadrature& volume_quadrature,
                                       const Quadrature& surface_quadrature)
  : LagrangeBaseMapping(
      grid, cell, 8, MakeFaceNodeMapping(cell), volume_quadrature, surface_quadrature)
{
}

double
LagrangeHexMapping::RefShape(uint32_t i, const Vec3& qpoint) const
{
  ChiLogicalErrorIf(i >= 8, "Invalid shapefunction index " + std::to_string(i));

  double a, b, c;
  // clang-format off
  if (i == 0) {a = -1.0; b=-1.0; c=-1.0;}
  if (i == 1) {a =  1.0; b=-1.0; c=-1.0;}
  if (i == 2) {a =  1.0; b= 1.0; c=-1.0;}
  if (i == 3) {a = -1.0; b= 1.0; c=-1.0;}

  if (i == 4) {a = -1.0; b=-1.0; c= 1.0;}
  if (i == 5) {a =  1.0; b=-1.0; c= 1.0;}
  if (i == 6) {a =  1.0; b= 1.0; c= 1.0;}
  if (i == 7) {a = -1.0; b= 1.0; c= 1.0;}

  return 0.125 *
    (1.0 + a * qpoint.x) * (1.0 + b * qpoint.y) * (1.0 + c * qpoint.z);
  // clang-format on
}

Vector3
LagrangeHexMapping::RefGradShape(uint32_t i, const Vec3& qpoint) const
{
  ChiLogicalErrorIf(i >= 8, "Invalid shapefunction index " + std::to_string(i));

  double a, b, c;
  // clang-format off
  if (i == 0) {a = -1.0; b=-1.0; c=-1.0;}
  if (i == 1) {a =  1.0; b=-1.0; c=-1.0;}
  if (i == 2) {a =  1.0; b= 1.0; c=-1.0;}
  if (i == 3) {a = -1.0; b= 1.0; c=-1.0;}

  if (i == 4) {a = -1.0; b=-1.0; c= 1.0;}
  if (i == 5) {a =  1.0; b=-1.0; c= 1.0;}
  if (i == 6) {a =  1.0; b= 1.0; c= 1.0;}
  if (i == 7) {a = -1.0; b= 1.0; c= 1.0;}
  // clang-format on

  const double x = qpoint.x;
  const double y = qpoint.y;
  const double z = qpoint.z;

  return Vec3(0.125 * a * (1.0 + b * y) * (1.0 + c * z),
              0.125 * b * (1.0 + a * x) * (1.0 + c * z),
              0.125 * c * (1.0 + a * x) * (1.0 + b * y));
}

LagrangeBaseMapping::MatDbl
LagrangeHexMapping::RefJacobian(const Vec3& qpoint) const
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
LagrangeHexMapping::RefFaceJacobianDeterminantAndNormal(size_t face_index,
                                                        const Vec3& qpoint_face) const
{
  const auto& x0 = node_locations_[face_node_mappings_[face_index][0]];
  const auto& x1 = node_locations_[face_node_mappings_[face_index][1]];
  const auto& x2 = node_locations_[face_node_mappings_[face_index][2]];
  const auto& x3 = node_locations_[face_node_mappings_[face_index][3]];

  const double xbar = qpoint_face.x;
  const double ybar = qpoint_face.y;
  Vec3 dx_dxbar;
  Vec3 dx_dybar;
  for (size_t i = 0; i < 4; ++i)
  {
    double a, b;
    Vec3 xi;
    // clang-format off
    if (i == 0) {a = -1.0; b=-1.0; xi = x0;}
    if (i == 1) {a =  1.0; b=-1.0; xi = x1;}
    if (i == 2) {a =  1.0; b= 1.0; xi = x2;}
    if (i == 3) {a = -1.0; b= 1.0; xi = x3;}
    // clang-format on

    const double ab = a * b;
    dx_dxbar += 0.25 * (a + ab * ybar) * xi;
    dx_dybar += 0.25 * (b + ab * xbar) * xi;
  }

  const auto cross = dx_dxbar.Cross(dx_dybar);
  const double detJ = cross.Norm();

  return {detJ, cross / detJ};
}

LagrangeBaseMapping::Vec3
LagrangeHexMapping::FaceToElementQPointConversion(size_t face_index, const Vec3& qpoint_face) const
{
  const double x = qpoint_face.x + 1.0;
  const double y = qpoint_face.y + 1.0;
  // clang-format off
  if (face_index == 0/*XMAX*/) return Vec3( 1.0, -1.0, -1.0) + Vec3(0.0,  x,  y);
  if (face_index == 1/*XMIN*/) return Vec3(-1.0,  1.0, -1.0) + Vec3(0.0, -x,  y);
  if (face_index == 2/*YMAX*/) return Vec3( 1.0,  1.0, -1.0) + Vec3( -x,0.0,  y);
  if (face_index == 3/*YMIN*/) return Vec3(-1.0, -1.0, -1.0) + Vec3(  x,0.0,  y);
  if (face_index == 4/*ZMAX*/) return Vec3(-1.0, -1.0,  1.0) + Vec3(  x,  y,0.0);
  if (face_index == 5/*ZMIN*/) return Vec3(-1.0,  1.0, -1.0) + Vec3(  x, -y,0.0);
  // clang-format on

  ChiLogicalError("Invalid face index " + std::to_string(face_index));
}

} // namespace opensn
