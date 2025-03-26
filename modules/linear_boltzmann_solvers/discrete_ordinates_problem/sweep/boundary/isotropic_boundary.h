// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include "framework/mesh/mesh.h"
#include "framework/math/math.h"
#include <vector>
#include <limits>

namespace opensn
{

/// Specified isotropic incident fluxes on a boundary.
class IsotropicBoundary : public SweepBoundary
{
private:
  std::vector<double> boundary_flux_;

public:
  explicit IsotropicBoundary(size_t num_groups,
                             std::vector<double> boundary_flux,
                             CoordinateSystemType coord_type = CoordinateSystemType::CARTESIAN)
    : SweepBoundary(LBSBoundaryType::ISOTROPIC, num_groups, coord_type),
      boundary_flux_(std::move(boundary_flux))
  {
  }

  double* PsiIncoming(uint64_t cell_local_id,
                      unsigned int face_num,
                      unsigned int fi,
                      unsigned int angle_num,
                      int group_num) override
  {
    return &boundary_flux_[group_num];
  }
};

} // namespace opensn
