// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"

namespace opensn
{

class DiscreteOrdinatesProblem;

class AAHCSDASweepChunk : public SweepChunk
{
public:
  AAHCSDASweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset);

  void Sweep(AngleSet& angle_set) override;
  void ZeroDestinationPhi() override;

private:
  DiscreteOrdinatesProblem& problem_;
  std::vector<double>& destination_phi_e_;
};

} // namespace opensn
