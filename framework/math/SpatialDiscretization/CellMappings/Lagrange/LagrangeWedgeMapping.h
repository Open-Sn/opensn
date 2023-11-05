#pragma once

#include "opensn/framework/math/SpatialDiscretization/CellMappings/LagrangeBaseMapping.h"

namespace chi_math::cell_mapping
{

/**Lagrange element mapping for a wedge (extruded triangle).
 * \ingroup doc_CellMappings*/
class LagrangeWedgeMapping : public LagrangeBaseMapping
{
public:
  LagrangeWedgeMapping(const chi_mesh::MeshContinuum& grid,
                       const chi_mesh::Cell& cell,
                       const Quadrature& volume_quadrature,
                       const Quadrature& surface_quadrature,
                       const Quadrature& aux_surface_quadrature);

protected:
  double RefShape(uint32_t i, const Vec3& qpoint) const override;
  Vec3 RefGradShape(uint32_t i, const Vec3& qpoint) const override;

  MatDbl RefJacobian(const Vec3& qpoint) const override;

  Vec3 FaceToElementQPointConversion(size_t face_index, const Vec3& qpoint_face) const override;

  std::pair<double, Vec3>
  RefFaceJacobianDeterminantAndNormal(size_t face_index, const Vec3& qpoint_face) const override;

  const Quadrature& GetSurfaceQuadrature(size_t face_index) const override;

  // for use on the bottom and top triangular faces
  const Quadrature& aux_surface_quadrature_;
};

} // namespace chi_math::cell_mapping
