// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/isotropic_boundary.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/mesh/cell/cell.h"

namespace opensn
{

IsotropicBoundary::IsotropicBoundary(BoundaryBank& bank,
                                     const std::vector<LBSGroupset>& groupsets,
                                     std::vector<double>& boundary_flux)
  : SweepBoundary(bank, LBSBoundaryType::ISOTROPIC)
{
  for (const auto& groupset : groupsets)
  {
    offset_[groupset.id] = bank_[groupset.id].counter++;
    double* flux = bank_.ExtendBoundaryFlux(groupset.id, groupset.GetNumGroups());
    for (unsigned int g = 0; g < groupset.GetNumGroups(); ++g)
    {
      flux[g] = boundary_flux[groupset.first_group + g];
    }
  }
}

std::uint64_t
IsotropicBoundary::GetOffsetToAngleset(const FaceNode& face_node,
                                       AngleSet& anglset,
                                       bool is_outgoing)
{
  if (is_outgoing)
    throw std::logic_error("IsotropicBoundary does not support outgoing flux.");
  int groupset_id = anglset.GetGroupsetID();
  return offset_[groupset_id] * bank_[groupset_id].groupset_size;
}

} // namespace opensn
