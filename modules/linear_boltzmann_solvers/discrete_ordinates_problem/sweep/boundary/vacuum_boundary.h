// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include "framework/mesh/mesh.h"
#include "framework/math/math.h"
#include <algorithm>
#include <vector>

namespace opensn
{

class VacuumBoundary : public SweepBoundary
{
public:
  explicit VacuumBoundary(BoundaryBank& bank) : SweepBoundary(bank, LBSBoundaryType::VACUUM)
  {
    std::fill(offset_.begin(), offset_.end(), 0);
  }

  double* PsiIncoming(std::uint32_t cell_local_id,
                      unsigned int face_num,
                      unsigned int fi,
                      unsigned int angle_num,
                      int groupset_id,
                      unsigned int group_idx) override
  {
    return ZeroFlux(groupset_id, group_idx);
  }

  std::uint64_t
  GetOffsetToAngleset(const FaceNode& face_node, AngleSet& anglset, bool is_outgoing) override
  {
    if (is_outgoing)
      throw std::logic_error("VacuumBoundary does not support outgoing flux.");
    return 0;
  }
};

} // namespace opensn
