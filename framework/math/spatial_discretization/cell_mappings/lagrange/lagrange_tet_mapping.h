#pragma once

#include "framework/math/spatial_discretization/cell_mappings/lagrange_base_mapping.h"

namespace opensn
{

/**Lagrange element mapping for a tetrahedron.
 * \ingroup doc_CellMappings*/
class LagrangeTetMapping : public LagrangeBaseMapping
{
public:
  LagrangeTetMapping(const MeshContinuum& grid,
                     const Cell& cell,
                     const Quadrature& volume_quadrature,
                     const Quadrature& surface_quadrature);

protected:
  double RefShape(uint32_t i, const Vec3& qpoint) const override;
  Vec3 RefGradShape(uint32_t i, const Vec3& qpoint) const override;

  MatDbl RefJacobian(const Vec3& qpoint) const override;

  Vec3 FaceToElementQPointConversion(size_t face_index, const Vec3& qpoint_face) const override;

  std::pair<double, Vec3>
  RefFaceJacobianDeterminantAndNormal(size_t face_index, const Vec3& qpoint_face) const override;
};

} // namespace opensn
