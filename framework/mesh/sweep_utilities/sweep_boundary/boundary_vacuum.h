#pragma once

#include "framework/mesh/mesh.h"
#include "framework/math/math.h"
#include "framework/mesh/sweep_utilities/sweep_boundary/sweep_boundary.h"
#include <vector>
#include <limits>

namespace chi_mesh::sweep_management
{

/**
 * Zero fluxes homogenous on a boundary and in angle.
 */
class BoundaryVaccuum : public SweepBoundary
{
private:
  std::vector<double> boundary_flux_;

public:
  explicit BoundaryVaccuum(
    size_t in_num_groups,
    chi_math::CoordinateSystemType coord_type = chi_math::CoordinateSystemType::CARTESIAN)
    : SweepBoundary(BoundaryType::INCIDENT_VACCUUM, in_num_groups, coord_type),
      boundary_flux_(in_num_groups, 0.0)
  {
  }

  double* HeterogeneousPsiIncoming(uint64_t cell_local_id,
                                   unsigned int face_num,
                                   unsigned int fi,
                                   unsigned int angle_num,
                                   int group_num,
                                   size_t gs_ss_begin) override;
};

} // namespace chi_mesh::sweep_management
