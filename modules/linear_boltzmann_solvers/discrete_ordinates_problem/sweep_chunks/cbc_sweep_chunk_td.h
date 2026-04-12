// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/avx_sweep_chunk_utils.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_chunk_shared.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"

namespace opensn
{

/**
 * Time-dependent host-side CBC sweep chunk.
 *
 * Identical to CBCSweepChunk but instantiates the Generic and FixedN kernels
 * with \c time_dependent=true, adding the \f$v_g^{-1}/(\theta\Delta t)\f$
 * time-absorption term and the previous-time-step angular flux source.
 */
class CBCSweepChunkTD : public SweepChunk
{
public:
  CBCSweepChunkTD(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset);
  ~CBCSweepChunkTD() override = default;

  /// Cache angle-set-level data and select the FixedN or Generic kernel.
  void SetAngleSet(AngleSet& angle_set) override;
  /// Cache cell-level data for the next Sweep call.
  void SetCell(const Cell* cell_ptr, AngleSet& angle_set) override;
  /// Sweep the current cell for all angles and groups (time-dependent).
  void Sweep(AngleSet& angle_set) override;
  /// Indicate this chunk uses the time-dependent kernel variant.
  bool IsTimeDependent() const override { return true; }

protected:
  /// Pointer-to-member for the selected sweep implementation.
  using SweepFunc = void (CBCSweepChunkTD::*)(AngleSet&);
  /// Construct the aggregated sweep data struct for the current cell.
  CBCSweepData MakeSweepData(const std::vector<double>* psi_old);
  /// Sweep using the generic (dynamic-size) kernel.
  void Sweep_Generic(AngleSet& angle_set);
  /// Sweep using the FixedN (compile-time node count) kernel.
  template <unsigned int NumNodes>
  void Sweep_FixedN(AngleSet& angle_set);

  /// Owning discrete ordinates problem.
  DiscreteOrdinatesProblem& problem_;
  /// Previous-time-step angular flux vector.
  const std::vector<double>& psi_old_;
  /// Energy group block size for SIMD batch solve.
  unsigned int group_block_size_ = 0;
  /// Cached per-cell and per-angle-set context.
  CBCSweepChunkContext ctx_;
  /// Reusable scratch buffers for the Generic kernel.
  CBCGenericSweepScratch generic_scratch_;

private:
  /// Selected sweep function pointer (Generic or FixedN).
  SweepFunc sweep_impl_td_ = nullptr;
};

} // namespace opensn