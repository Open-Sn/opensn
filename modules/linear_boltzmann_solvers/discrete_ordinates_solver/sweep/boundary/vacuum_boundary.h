// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/boundary/sweep_boundary.h"
#include "framework/mesh/mesh.h"
#include "framework/math/math.h"
#include <vector>
#include <limits>

namespace opensn
{
namespace lbs
{

class VacuumBoundary : public SweepBoundary
{
private:
  std::vector<double> boundary_flux_;

public:
  explicit VacuumBoundary(size_t num_groups,
                          CoordinateSystemType coord_type = CoordinateSystemType::CARTESIAN)
    : SweepBoundary(BoundaryType::VACUUM, num_groups, coord_type), boundary_flux_(num_groups, 0.0)
  {
  }

  double* PsiIncoming(uint64_t cell_local_id,
                      unsigned int face_num,
                      unsigned int fi,
                      unsigned int angle_num,
                      int group_num,
                      size_t gs_ss_begin) override
  {
    return &boundary_flux_[group_num];
  }
};

} // namespace lbs
} // namespace opensn
