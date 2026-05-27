// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_kernels.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"
#include <vector>

namespace opensn
{

class DiscreteOrdinatesProblem;
class LBSGroupset;
class CurvilinearProductQuadrature;

/// Host CBC sweep chunk for 2D RZ curvilinear coordinates.
class CBCSweepChunkRZ : public SweepChunk
{
public:
  CBCSweepChunkRZ(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset);

  void SetAngleSet(AngleSet& angle_set) override;

  void Sweep(AngleSet& angle_set) override;

private:
  using SweepFunc = void (CBCSweepChunkRZ::*)(AngleSet&);

  void Sweep_Generic(AngleSet& angle_set);

  template <unsigned int NumNodes>
  void Sweep_FixedN(AngleSet& angle_set);

  /// Curvilinear angular quadrature used for RZ streaming factors.
  const CurvilinearProductQuadrature& curvilinear_quadrature_;
  /// Secondary spatial discretization cell matrices.
  const std::vector<UnitCellMatrices>& secondary_unit_cell_matrices_;
  /// Unknown manager for polar-level sweep dependency storage.
  UnknownManager unknown_manager_;
  /// Sweeping dependency angular intensity for each polar level.
  std::vector<double> psi_sweep_;
  /// Direction index to polar level mapping.
  std::vector<unsigned int> direction_polar_level_;
  /// Normal vector used to identify the point/axis of symmetry.
  Vector3 normal_vector_boundary_;

  CBC_FLUDS* fluds_ = nullptr;
  CBC_AsynchronousCommunicator* async_comm_ = nullptr;

  /// Number of groups solved in one block.
  unsigned int group_block_size_;
  CBCSweepWorkspace workspace_;
  /// Selected sweep implementation.
  SweepFunc sweep_impl_ = nullptr;
};

} // namespace opensn
