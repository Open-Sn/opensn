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
                                     const std::vector<double>& boundary_flux,
                                     double start_time,
                                     double end_time)
  : SweepBoundary(bank, LBSBoundaryType::ISOTROPIC),
    active_boundary_flux_(boundary_flux),
    start_time_(start_time),
    end_time_(end_time)
{
  for (const auto& groupset : groupsets)
  {
    offset_[groupset.id] = bank_[groupset.id].counter++;
    bank_.ExtendBoundaryFlux(groupset.id, groupset.GetNumGroups());
  }
  IsotropicBoundary::UpdateBoundaryFlux(groupsets);
}

void
IsotropicBoundary::UpdateBoundaryFlux(const std::vector<LBSGroupset>& groupsets)
{
  bool is_active = (evaluation_time_ >= start_time_ and evaluation_time_ <= end_time_);
  if (current_state_ != is_active)
  {
    for (const auto& groupset : groupsets)
    {
      double* boundary_flux = GetBoundaryFlux(groupset.id);
      for (unsigned int g = 0; g < groupset.GetNumGroups(); ++g)
      {
        boundary_flux[g] = ((is_active) ? active_boundary_flux_[groupset.first_group + g] : 0.0);
      }
    }
    current_state_ = is_active;
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
