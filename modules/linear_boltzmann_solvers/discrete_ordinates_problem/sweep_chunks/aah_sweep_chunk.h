// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/avx_sweep_chunk_utils.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include <cstdint>
#include <array>
#include <cmath>

namespace opensn
{

class DiscreteOrdinatesProblem;

class AAHSweepChunk : public SweepChunk
{
public:
  AAHSweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset);

  void Sweep(AngleSet& angle_set) override;

protected:
  DiscreteOrdinatesProblem& problem_;
  unsigned int group_block_size_;

private:
  using SweepFunc = void (AAHSweepChunk::*)(AngleSet&);
  SweepFunc sweep_impl_ = nullptr;

  void Sweep_Generic(AngleSet& angle_set);
  template <unsigned int NumNodes>
  void Sweep_FixedN(AngleSet& angle_set);
};

} // namespace opensn
