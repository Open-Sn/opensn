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

class VacuumBoundary : public SweepBoundary
{
public:
  explicit VacuumBoundary(unsigned int num_groups,
                          CoordinateSystemType coord_type = CoordinateSystemType::CARTESIAN)
    : SweepBoundary(LBSBoundaryType::VACUUM, num_groups, coord_type)
  {
  }

  double* PsiIncoming(std::uint32_t cell_local_id,
                      unsigned int face_num,
                      unsigned int fi,
                      unsigned int angle_num,
                      unsigned int group_num) override
  {
    return ZeroFlux(group_num);
  }

private:
};

} // namespace opensn
