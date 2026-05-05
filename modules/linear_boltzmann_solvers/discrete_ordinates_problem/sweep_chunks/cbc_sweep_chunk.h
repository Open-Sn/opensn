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

/**
 * Host CBC sweep chunk.
 *
 * Dispatches between the generic and fixed-node CBC sweep kernels for the
 * currently bound angle set and cell.
 */
class CBCSweepChunk : public SweepChunk
{
public:
  /**
   * Construct one CBC sweep chunk for a groupset.
   *
   * \param problem Owning discrete-ordinates problem.
   * \param groupset Groupset swept by this chunk.
   */
  CBCSweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset);

  /**
   * Bind the current angle set.
   *
   * \param angle_set Angle set to bind.
   */
  void SetAngleSet(AngleSet& angle_set) override;

  /**
   * Bind the current cell to be swept.
   *
   * \param cell_ptr Cell to bind.
   * \param angle_set Owning angle set.
   */
  void SetCell(Cell const* cell_ptr, AngleSet& angle_set) override;

  /**
   * Sweep the currently bound cell for the provided angle set.
   *
   * Selects the fixed-node kernel when all local cells have the same node count
   * in the supported range, otherwise falls back to the generic CBC kernel.
   *
   * \param angle_set Angle set currently being advanced.
   */
  void Sweep(AngleSet& angle_set) override;

protected:
  /// Owning discrete-ordinates problem.
  DiscreteOrdinatesProblem& problem_;
  /// Cached per-cell and per-angleset context.
  CBCSweepChunkContext ctx_;
  /// Group block size for SIMD batch solves.
  unsigned int group_block_size_ = 0;
  /// Reusable scratch buffers for generic sweep chunk kernel.
  CBCGenericSweepScratch generic_scratch_;

private:
  /// Pointer-to-member for the selected sweep implementation (generic or fixed-node).
  using SweepFunc = void (CBCSweepChunk::*)(AngleSet&);
  /// Selected sweep function pointer (generic or fixed-node).
  SweepFunc sweep_impl_ = nullptr;

  /// Construct the aggregated sweep data for the current cell.
  CBCSweepData MakeSweepData(const std::vector<double>* psi_old);
  /// Sweep using the generic kernel.
  void Sweep_Generic(AngleSet& angle_set);
  /// Sweep using the fixed-node kernel.
  template <unsigned int NumNodes>
  void Sweep_FixedN(AngleSet& angle_set);
};

} // namespace opensn
