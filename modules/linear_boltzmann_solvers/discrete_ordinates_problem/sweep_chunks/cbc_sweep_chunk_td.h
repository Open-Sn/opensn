// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_kernels.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"

namespace opensn
{

/// Host CBC sweep chunk for transient problems.
class CBCSweepChunkTD : public SweepChunk
{
public:
  CBCSweepChunkTD(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset);

  void SetAngleSet(AngleSet& angle_set) override;

  void Sweep(AngleSet& angle_set) override;

  bool IsTimeDependent() const override { return true; }

protected:
  template <bool time_dependent, class SweepChunkT>
  friend void CBC_Sweep_Generic(SweepChunkT& sweep_chunk, AngleSet& angle_set);

  template <unsigned int NumNodes, bool time_dependent, class SweepChunkT>
  friend void CBC_Sweep_FixedN(SweepChunkT& sweep_chunk, AngleSet& angle_set);

  using SweepFunc = void (CBCSweepChunkTD::*)(AngleSet&);

  void Sweep_Generic(AngleSet& angle_set);

  template <unsigned int NumNodes>
  void Sweep_FixedN(AngleSet& angle_set);

  DiscreteOrdinatesProblem& problem_;
  const std::vector<double>& psi_old_;

  CBC_FLUDS* fluds_ = nullptr;
  CBC_AsynchronousCommunicator* async_comm_ = nullptr;

  unsigned int group_block_size_ = 0;
  CBCSweepWorkspace workspace_;

private:
  SweepFunc sweep_impl_td_ = nullptr;
};

} // namespace opensn
