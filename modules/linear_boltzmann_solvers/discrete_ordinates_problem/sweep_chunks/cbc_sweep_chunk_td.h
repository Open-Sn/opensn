// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/avx_sweep_chunk_utils.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_chunk_shared.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"

namespace opensn
{

/// CBC sweep chunk for transient problems.
class CBCSweepChunkTD : public SweepChunk
{
public:
  CBCSweepChunkTD(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset);
  ~CBCSweepChunkTD() override = default;

  void SetAngleSet(AngleSet& angle_set) override;

  void SetCell(const Cell* cell_ptr, AngleSet& angle_set) override;

  void Sweep(AngleSet& angle_set) override;

  bool IsTimeDependent() const override { return true; }

protected:
  using SweepFunc = void (CBCSweepChunkTD::*)(AngleSet&);

  void Sweep_Generic(AngleSet& angle_set);

  template <unsigned int NumNodes>
  void Sweep_FixedN(AngleSet& angle_set);

  /// Owning discrete ordinates problem.
  DiscreteOrdinatesProblem& problem_;

  /// Previous-step angular flux vector.
  const std::vector<double>& psi_old_;

  /// Number of groups solved in one block.
  unsigned int group_block_size_ = 0;

  /// Reusable CBC sweep context.
  CBCSweepChunkContext ctx_;

private:
  /// Selected time-dependent sweep implementation.
  SweepFunc sweep_impl_td_ = nullptr;
};

} // namespace opensn
