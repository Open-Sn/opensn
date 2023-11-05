#pragma once

#include "mesh/chi_mesh.h"
#include "math/chi_math.h"
#include "sweep_boundary.h"
#include <vector>
#include <limits>

namespace chi_mesh::sweep_management
{

/**
 * Specified incident fluxes homogenous on a boundary.
 */
class BoundaryIsotropicHomogenous : public SweepBoundary
{
private:
  std::vector<double> boundary_flux;

public:
  explicit BoundaryIsotropicHomogenous(
    size_t in_num_groups,
    std::vector<double> ref_boundary_flux,
    chi_math::CoordinateSystemType coord_type = chi_math::CoordinateSystemType::CARTESIAN)
    : SweepBoundary(BoundaryType::INCIDENT_ISOTROPIC_HOMOGENOUS, in_num_groups, coord_type),
      boundary_flux(std::move(ref_boundary_flux))
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
