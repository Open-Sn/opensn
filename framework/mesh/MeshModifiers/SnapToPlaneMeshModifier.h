#pragma once

#include "MeshModifier.h"

#include "mesh/chi_mesh.h"

namespace chi_mesh
{

/**Object that will modify boundary elements that are close to the specified
 * direction to become perfectly aligned with that direction.*/
class SnapToPlaneMeshModifier : public MeshModifier
{
public:
  static chi::InputParameters GetInputParameters();
  explicit SnapToPlaneMeshModifier(const chi::InputParameters& params);

  void Apply() override;

protected:
  const chi_mesh::Vector3 normal_;
  const chi_mesh::Vector3 point_;
  const bool boundary_nodes_only_ = true;
  const bool check_face_alignment_ = false;
  const double tol_ = 1.0e-5;
};

} // namespace chi_mesh


