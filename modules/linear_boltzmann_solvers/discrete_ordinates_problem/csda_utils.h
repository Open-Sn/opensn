// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_structs.h"
#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace opensn
{

/**
 * Returns the problem's charged-particle blocks: contiguous ranges of groups whose
 * stopping power is nonzero in at least one material, as half-open [begin, end) ranges
 * ordered from high to low energy.
 */
inline std::vector<std::pair<unsigned int, unsigned int>>
FindCSDAProblemChargedGroupRanges(const BlockID2XSMap& xs_map, const unsigned int num_groups)
{
  std::vector<bool> active(num_groups, false);
  for (const auto& [_, xs] : xs_map)
    for (const auto& [begin, end] : xs->GetStoppingPowerGroupRanges())
      for (auto g = begin; g < std::min(end, num_groups); ++g)
        active[g] = true;

  std::vector<std::pair<unsigned int, unsigned int>> ranges;
  unsigned int g = 0;
  while (g < active.size())
  {
    while (g < active.size() and not active[g])
      ++g;
    if (g >= active.size())
      break;

    const unsigned int begin = g;
    while (g < active.size() and active[g])
      ++g;
    ranges.emplace_back(begin, g);
  }
  return ranges;
}

/**
 * Returns the charge sign of particles in group g: +1 in the first problem-level
 * charged-particle block (electrons) and -1 in the second (positrons). Validation
 * limits a problem to at most two blocks.
 */
inline double
CSDAChargeSign(const std::vector<std::pair<unsigned int, unsigned int>>& problem_ranges,
               const unsigned int g)
{
  return (problem_ranges.empty() or g < problem_ranges.front().second) ? 1.0 : -1.0;
}

} // namespace opensn
