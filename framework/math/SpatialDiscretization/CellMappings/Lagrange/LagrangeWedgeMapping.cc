#include "framework/math/SpatialDiscretization/CellMappings/Lagrange/LagrangeWedgeMapping.h"

#include "framework/mesh/Cell/cell.h"

#include "framework/logging/chi_log_exceptions.h"

namespace chi_math::cell_mapping
{

LagrangeWedgeMapping::LagrangeWedgeMapping(const chi_mesh::MeshContinuum& grid,
                                           const chi_mesh::Cell& cell,
                                           const Quadrature& volume_quadrature,
                                           const Quadrature& surface_quadrature,
                                           const Quadrature& aux_surface_quadrature)
  : LagrangeBaseMapping(
      grid, cell, 6, MakeFaceNodeMapping(cell), volume_quadrature, surface_quadrature),
    aux_surface_quadrature_(aux_surface_quadrature)
{
}

double
LagrangeWedgeMapping::RefShape(uint32_t i, const Vec3& qpoint) const
{
  const double x = qpoint.x;
  const double y = qpoint.y;
  const double z = qpoint.z;

  // clang-format off
  if (i == 0) return (1.0 - x - y) * 0.5 * (1.0 - z);
  if (i == 1) return x * 0.5 * (1.0 - z);
  if (i == 2) return y * 0.5 * (1.0 - z);

  if (i == 3) return (1.0 - x - y) * 0.5 * (1.0 + z);
  if (i == 4) return x * 0.5 * (1.0 + z);
  if (i == 5) return y * 0.5 * (1.0 + z);
  // clang-format on

  ChiLogicalError("Invalid shapefunction index " + std::to_string(i));
}

chi_mesh::Vector3
LagrangeWedgeMapping::RefGradShape(uint32_t i, const Vec3& qpoint) const
{

  const double x = qpoint.x;
  const double y = qpoint.y;
  const double z = qpoint.z;

  // clang-format off
  if (i == 0) return Vec3(-0.5 * (1.0 - z), -0.5 * (1.0 - z), -0.5 * (1.0 - x - y));
  if (i == 1) return Vec3( 0.5 * (1.0 - z),              0.0, -0.5 * x);
  if (i == 2) return Vec3(             0.0,  0.5 * (1.0 - z), -0.5 * y);

  if (i == 3) return Vec3(-0.5 * (1.0 + z), -0.5 * (1.0 + z), 0.5 * (1.0 - x - y));
  if (i == 4) return Vec3( 0.5 * (1.0 + z),              0.0, 0.5 * x);
  if (i == 5) return Vec3(             0.0,  0.5 * (1.0 + z), 0.5 * y);
  // clang-format on

  ChiLogicalError("Invalid shapefunction index " + std::to_string(i));
}

LagrangeBaseMapping::MatDbl
LagrangeWedgeMapping::RefJacobian(const Vec3& qpoint) const
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
LagrangeWedgeMapping::RefFaceJacobianDeterminantAndNormal(size_t face_index,
                                                          const Vec3& qpoint_face) const
{
  if (face_index <= 2)
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

    return {cross.Norm(), cross.Normalized()};
  }
  else
  {
    const auto& x0 = node_locations_[face_node_mappings_[face_index][0]];
    const auto& x1 = node_locations_[face_node_mappings_[face_index][1]];
    const auto& x2 = node_locations_[face_node_mappings_[face_index][2]];

    Vec3 dx_dxbar;
    Vec3 dx_dybar;
    for (size_t i = 0; i < 3; ++i)
    {
      double dN_dxbar, dN_dybar;
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
}

LagrangeBaseMapping::Vec3
LagrangeWedgeMapping::FaceToElementQPointConversion(size_t face_index,
                                                    const Vec3& qpoint_face) const
{
  if (face_index <= 2)
  {

    const double x = 0.5 * (qpoint_face.x + 1.0);
    const double y = qpoint_face.y + 1.0;

    // clang-format off
    if (face_index == 0/*Side0*/) return Vec3(    x,     0.0, -1.0 + y);
    if (face_index == 1/*Side1*/) return Vec3(1.0-x,       x, -1.0 + y);
    if (face_index == 2/*Side2*/) return Vec3(  0.0,   1.0-x, -1.0 + y);
    // clang-format on
  }
  else
  {
    const double x = qpoint_face.x;
    const double y = qpoint_face.y;
    // clang-format off
    if (face_index == 3/*ZMAX*/) return Vec3( x, y,  1.0);
    if (face_index == 4/*ZMIN*/) return Vec3( x, y, -1.0);
    // clang-format on
  }

  ChiLogicalError("Invalid face index " + std::to_string(face_index));
}

const Quadrature&
LagrangeWedgeMapping::GetSurfaceQuadrature(size_t face_index) const
{
  if (face_index <= 2) return surface_quadrature_;
  else
    return aux_surface_quadrature_;
}

} // namespace chi_math::cell_mapping
