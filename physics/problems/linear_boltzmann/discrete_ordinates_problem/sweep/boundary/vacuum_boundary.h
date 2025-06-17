// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "physics/problems/linear_boltzmann/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include "framework/mesh/mesh.h"
#include "framework/math/math.h"
#include <vector>
#include <limits>

namespace opensn
{

class VacuumBoundary : public SweepBoundary
{
private:
  std::vector<double> boundary_flux_;

public:
  explicit VacuumBoundary(size_t num_groups,
                          CoordinateSystemType coord_type = CoordinateSystemType::CARTESIAN)
    : SweepBoundary(LBSBoundaryType::VACUUM, num_groups, coord_type),
      boundary_flux_(num_groups, 0.0)
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
