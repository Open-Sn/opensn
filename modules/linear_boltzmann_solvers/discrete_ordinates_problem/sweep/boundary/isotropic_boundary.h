// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include "framework/mesh/mesh.h"
#include "framework/math/math.h"
#include <vector>

namespace opensn
{

/// Specified isotropic incident fluxes on a boundary.
class IsotropicBoundary : public SweepBoundary
{
public:
  explicit IsotropicBoundary(BoundaryBank& bank,
                             const std::vector<LBSGroupset>& groupsets,
                             std::vector<double>& boundary_flux);

  double* PsiIncoming(std::uint32_t cell_local_id,
                      unsigned int face_num,
                      unsigned int fi,
                      unsigned int angle_num,
                      int groupset_id,
                      unsigned int group_idx) override
  {
    return GetBoundaryFlux(groupset_id) + group_idx;
  }

  std::uint64_t
  GetOffsetToAngleset(const FaceNode& face_node, AngleSet& anglset, bool is_outgoing) override;
};

} // namespace opensn
