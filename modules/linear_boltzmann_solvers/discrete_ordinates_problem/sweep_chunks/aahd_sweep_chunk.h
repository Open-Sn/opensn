// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include <cstdint>
#include <array>
#include <cmath>

namespace opensn
{

class DiscreteOrdinatesProblem;

class AAHDSweepChunk : public SweepChunk
{
public:
  AAHDSweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset);

  DiscreteOrdinatesProblem& GetProblem() { return problem_; }
  MeshContinuum& GetGrid() { return *grid_; }
  const LBSGroupset& GetGroupset() const { return groupset_; }

  void Sweep(AngleSet& angle_set) override;

protected:
  DiscreteOrdinatesProblem& problem_;
};

} // namespace opensn
