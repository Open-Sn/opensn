// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include <cstdint>
#include <array>
#include <cmath>

namespace opensn
{

// experimental, to be moved to a higher level header file
static constexpr size_t simd_width =
#if defined(__AVX512F__)
  8; // 8 lanes (512-bit, doubles)
#elif defined(__AVX2__)
  4; // 4 lanes (256-bit, doubles)
#else
  1; // scalar
#endif

class DiscreteOrdinatesProblem;

class AAHSweepChunk : public SweepChunk
{
public:
  AAHSweepChunk(const std::shared_ptr<MeshContinuum>& grid,
                const SpatialDiscretization& discretization,
                const std::vector<UnitCellMatrices>& unit_cell_matrices,
                std::vector<CellLBSView>& cell_transport_views,
                const std::vector<double>& densities,
                std::vector<double>& destination_phi,
                std::vector<double>& destination_psi,
                const std::vector<double>& source_moments,
                const LBSGroupset& groupset,
                const std::map<int, std::shared_ptr<MultiGroupXS>>& xs,
                int num_moments,
                int max_num_cell_dofs,
                int min_num_cell_dofs,
                DiscreteOrdinatesProblem& problem,
                size_t max_level_size,
                size_t max_groupset_size,
                size_t max_angleset_size,
                bool use_gpus);

  ~AAHSweepChunk() = default;

  void Sweep(AngleSet& angle_set) override;

protected:

  void CPUSweep(AngleSet& angle_set);
  void GPUSweep(AngleSet& angle_set);

  DiscreteOrdinatesProblem& problem_;
  size_t max_level_size_;
  size_t group_block_size_;
  bool use_gpus_;
  void* level_vector_ = nullptr;

private:
  using CpuSweepFunc = void (AAHSweepChunk::*)(AngleSet&);
  CpuSweepFunc cpu_sweep_impl_ = nullptr;

  void CPUSweep_Generic(AngleSet& angle_set);
  void CPUSweep_N4(AngleSet& angle_set);
};

} // namespace opensn
