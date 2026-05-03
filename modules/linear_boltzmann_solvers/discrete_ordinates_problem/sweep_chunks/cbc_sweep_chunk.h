// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_kernels.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"

namespace opensn
{

class DiscreteOrdinatesProblem;

/// Host CBC sweep chunk.
class CBCSweepChunk : public SweepChunk
{
public:
  CBCSweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset);

  void SetAngleSet(AngleSet& angle_set) override;

  void Sweep(AngleSet& angle_set) override;

protected:
  template <bool time_dependent, class SweepChunkT>
  friend void CBC_Sweep_Generic(SweepChunkT& sweep_chunk, AngleSet& angle_set);

  template <unsigned int NumNodes, bool time_dependent, class SweepChunkT>
  friend void CBC_Sweep_FixedN(SweepChunkT& sweep_chunk, AngleSet& angle_set);

  CBC_FLUDS* fluds_ = nullptr;
  CBC_AsynchronousCommunicator* async_comm_ = nullptr;

  unsigned int group_block_size_ = 0;
  CBCSweepWorkspace workspace_;

private:
  using SweepFunc = void (CBCSweepChunk::*)(AngleSet&);

  SweepFunc sweep_impl_ = nullptr;

  void Sweep_Generic(AngleSet& angle_set);

  template <unsigned int NumNodes>
  void Sweep_FixedN(AngleSet& angle_set);
};

} // namespace opensn
