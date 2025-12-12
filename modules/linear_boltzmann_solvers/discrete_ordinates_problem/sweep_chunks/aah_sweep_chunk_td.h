// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include <cstdint>
#include <array>
#include <cmath>
#include <map>
#include <vector>

namespace opensn
{

class DiscreteOrdinatesProblem;

class AAHSweepChunkTD : public SweepChunk
{
public:
  AAHSweepChunkTD(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset);

  ~AAHSweepChunkTD() override;

  void Sweep(AngleSet& angle_set) override;

protected:
  void CPUSweep(AngleSet& angle_set);

  DiscreteOrdinatesProblem& problem_;
  size_t max_level_size_;
  void* level_vector_ = nullptr;
  const std::vector<double>& psi_old_;
  bool use_gpus_;
};

} // namespace opensn
