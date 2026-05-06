// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include <algorithm>

namespace opensn
{

AngleSet::AngleSet(size_t id,
                   const LBSGroupset& groupset,
                   const SPDS& spds,
                   std::shared_ptr<FLUDS>& fluds,
                   const std::vector<size_t>& angle_indices,
                   std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries)
  : id_(id),
    groupset_id_(groupset.id),
    num_groups_(groupset.GetNumGroups()),
    spds_(spds),
    fluds_(fluds),
    angles_(angle_indices.begin(), angle_indices.end()),
    boundaries_(boundaries)
{
}

bool
AngleSet::HasAngleIndex(std::uint32_t angle_index) const
{
  return std::find(angles_.begin(), angles_.end(), angle_index) != angles_.end();
}

void
AngleSet::UpdateSweepDependencies(std::set<AngleSet*>& following_angle_sets)
{
  following_angle_sets_.insert(
    following_angle_sets_.end(), following_angle_sets.begin(), following_angle_sets.end());
  for (auto* angleset : following_angle_sets_)
  {
    ++(angleset->num_dependencies_);
  }
}

} // namespace opensn
