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

/// CBC sweep chunk in point-symmetric and axial-symmetric curvilinear coordinates.
class CBCSweepChunkRZ : public SweepChunk
{
public:
  CBCSweepChunkRZ(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset);

  /// Set the CBC FLUDS and communicator for the active angle set.
  void SetAngleSet(AngleSet& angle_set) override;

  /// Sweep the active cell for an angle set.
  void Sweep(AngleSet& angle_set) override;

private:
  using SweepFunc = void (CBCSweepChunkRZ::*)(AngleSet&);

  /// Dispatch a cell to the fixed-size or generic implementation.
  void Sweep_Dispatch(AngleSet& angle_set);

  /// Sweep a cell whose node count is not handled by a fixed-size implementation.
  void Sweep_Generic(AngleSet& angle_set);

  /// Sweep a cell with a compile-time node count.
  template <unsigned int NumNodes>
  void Sweep_FixedN(AngleSet& angle_set);

  /// Curvilinear angular quadrature used for RZ streaming factors.
  const CurvilinearProductQuadrature& curvilinear_quadrature_;
  /// Secondary spatial discretization cell matrices.
  const std::vector<UnitCellMatrices>& secondary_unit_cell_matrices_;
  /// Unknown manager.
  UnknownManager unknown_manager_;
  /// Sweeping dependency angular intensity for each polar level.
  std::vector<double> psi_sweep_;
  /// Direction index to polar level mapping.
  std::vector<unsigned int> direction_polar_level_;
  /// Normal vector used to determine symmetric boundary condition.
  Vector3 normal_vector_boundary_;
  /// CBC FLUDS for the active angle set.
  CBC_FLUDS* fluds_ = nullptr;
  /// CBC async communicator for the active angle set.
  CBC_AsynchronousCommunicator* async_comm_ = nullptr;

  /// Number of groups solved in one block.
  unsigned int group_block_size_;
  /// Reusable CBC auxiliary sweep data storage.
  CBCSweepWorkspace workspace_;
  /// Selected sweep implementation.
  SweepFunc sweep_impl_ = nullptr;
};

} // namespace opensn
