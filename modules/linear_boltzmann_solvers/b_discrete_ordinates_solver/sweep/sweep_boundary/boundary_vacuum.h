#pragma once

#include "framework/mesh/mesh.h"
#include "framework/math/math.h"
#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep/sweep_boundary/sweep_boundary.h"
#include <vector>
#include <limits>

namespace opensn
{

/**
 * Zero fluxes homogenous on a boundary and in angle.
 */
class BoundaryVaccuum : public SweepBoundary
{
private:
  std::vector<double> boundary_flux_;

public:
  explicit BoundaryVaccuum(size_t in_num_groups,
                           CoordinateSystemType coord_type = CoordinateSystemType::CARTESIAN)
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

} // namespace opensn
