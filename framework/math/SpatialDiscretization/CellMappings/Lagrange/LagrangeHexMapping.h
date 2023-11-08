#pragma once

#include "framework/math/SpatialDiscretization/CellMappings/LagrangeBaseMapping.h"

namespace chi_math::cell_mapping
{

/**Lagrange element mapping for a hexahedron.
 * \ingroup doc_CellMappings*/
class LagrangeHexMapping : public LagrangeBaseMapping
{
public:
  LagrangeHexMapping(const chi_mesh::MeshContinuum& grid,
                     const chi_mesh::Cell& cell,
                     const Quadrature& volume_quadrature,
                     const Quadrature& surface_quadrature);

protected:
  double RefShape(uint32_t i, const Vec3& qpoint) const override;
  Vec3 RefGradShape(uint32_t i, const Vec3& qpoint) const override;

  MatDbl RefJacobian(const Vec3& qpoint) const override;

  std::pair<double, Vec3>
  RefFaceJacobianDeterminantAndNormal(size_t face_index, const Vec3& qpoint_face) const override;

  Vec3 FaceToElementQPointConversion(size_t face_index, const Vec3& qpoint_face) const override;
};

} // namespace chi_math::cell_mapping
