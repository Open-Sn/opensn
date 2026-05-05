// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/avx_sweep_chunk_utils.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_chunk_shared.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"

namespace opensn
{

class CellMapping;
class DiscreteOrdinatesProblem;

/// CBC sweep chunk.
class CBCSweepChunk : public SweepChunk
{
public:
  CBCSweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset);

  void SetAngleSet(AngleSet& angle_set) override;

  void SetCell(const Cell* cell_ptr, AngleSet& angle_set) override;

  void Sweep(AngleSet& angle_set) override;

protected:
  /// Owning discrete ordinates problem.
  DiscreteOrdinatesProblem& problem_;

  /// Reusable CBC sweep context.
  CBCSweepChunkContext ctx_;

  /// Number of groups solved in one block.
  unsigned int group_block_size_ = 0;

private:
  using SweepFunc = void (CBCSweepChunk::*)(AngleSet&);

  /// Selected sweep implementation.
  SweepFunc sweep_impl_ = nullptr;

  void Sweep_Generic(AngleSet& angle_set);

  template <unsigned int NumNodes>
  void Sweep_FixedN(AngleSet& angle_set);
};

} // namespace opensn
