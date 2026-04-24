// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cmath>
#include <utility>
#include <vector>

namespace opensn
{

inline constexpr double kCSDATolerance = 1.0e-12;

struct CSDATerminalDepositionRates
{
  double particle;
  double charge;
};

inline CSDATerminalDepositionRates
ComputeCSDATerminalDepositionRates(const double terminal_rate,
                                   const std::size_t charged_block_index)
{
  return {terminal_rate, charged_block_index == 0 ? terminal_rate : -terminal_rate};
}

inline std::vector<std::pair<unsigned int, unsigned int>>
FindCSDAChargedGroupRanges(const std::vector<double>& stopping_power)
{
  std::vector<std::pair<unsigned int, unsigned int>> ranges;

  unsigned int g = 0;
  while (g < stopping_power.size())
  {
    while (g < stopping_power.size() and std::abs(stopping_power[g]) <= kCSDATolerance)
      ++g;
    if (g >= stopping_power.size())
      break;

    const unsigned int g_begin = g;
    while (g < stopping_power.size() and std::abs(stopping_power[g]) > kCSDATolerance)
      ++g;
    ranges.emplace_back(g_begin, g);
  }

  return ranges;
}

inline double
CSDAGroupCenterEnergy(const std::vector<double>& energy_bounds, const unsigned int g)
{
  return 0.5 * (energy_bounds.at(g) + energy_bounds.at(g + 1));
}

} // namespace opensn
