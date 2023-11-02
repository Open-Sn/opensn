#include "LagrangeSlabMapping.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log_exceptions.h"

namespace chi_math::cell_mapping
{

LagrangeSlabMapping::LagrangeSlabMapping(const chi_mesh::MeshContinuum& grid,
                                         const chi_mesh::Cell& cell,
                                         const Quadrature& volume_quadrature,
                                         const Quadrature& surface_quadrature)
  : LagrangeBaseMapping(
      grid, cell, 2, MakeFaceNodeMapping(cell), volume_quadrature, surface_quadrature)

{
}

double
LagrangeSlabMapping::RefShape(uint32_t i, const Vec3& qpoint) const
{
  if (i == 0) return 0.5 * (1.0 - qpoint.x);
  if (i == 1) return 0.5 * (1.0 + qpoint.x);

  ChiLogicalError("Invalid shapefunction index " + std::to_string(i));
}

chi_mesh::Vector3
LagrangeSlabMapping::RefGradShape(uint32_t i, const Vec3& qpoint) const
{
  if (i == 0) return Vec3(0.0, 0.0, -0.5);
  if (i == 1) return Vec3(0.0, 0.0, 0.5);

  ChiLogicalError("Invalid shapefunction index " + std::to_string(i));
}

LagrangeBaseMapping::MatDbl
LagrangeSlabMapping::RefJacobian(const Vec3& qpoint) const
{
  const double zz = 0.5 * (node_locations_[1] - node_locations_[0]).Norm();
  return {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, zz}};
}

std::pair<double, LagrangeBaseMapping::Vec3>
LagrangeSlabMapping::RefFaceJacobianDeterminantAndNormal(size_t face_index,
                                                         const Vec3& qpoint_face) const
{
  return {1.0, cell_.faces_[face_index].normal_};
}

LagrangeBaseMapping::Vec3
LagrangeSlabMapping::FaceToElementQPointConversion(size_t face_index, const Vec3& qpoint_face) const
{
  if (face_index == 0) return Vec3(-1.0, 0.0, 0.0);
  if (face_index == 1) return Vec3(1.0, 0.0, 0.0);

  ChiLogicalError("Invalid face index " + std::to_string(face_index));
}

} // namespace chi_math::cell_mapping