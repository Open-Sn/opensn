// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep_runtime_builder.h"

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep_parallel_for.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "framework/utils/timer.h"

namespace opensn::detail
{

void
BuildAAHGPUFludsCommonData(SweepRuntime& runtime,
                           const SpatialDiscretization& discretization,
                           const std::vector<CellFaceNodalMapping>& grid_nodal_mappings)
{
  struct WorkItem
  {
    std::shared_ptr<AngularQuadrature> quadrature;
    std::shared_ptr<SPDS> spds;
  };
  std::vector<WorkItem> work;
  for (const auto& [quadrature, spds_list] : runtime.quadrature_spds_map)
    for (const auto& spds : spds_list)
      work.push_back({quadrature, spds});
  if (work.empty())
    return;

  const size_t nthreads = std::min<size_t>(std::max(1U, opensn_num_threads), work.size());
  log.Log() << program_timer.GetTimeString() << " FLUDS construction: " << work.size()
            << " angles (" << nthreads << " threads(s)).";
  std::vector<std::unique_ptr<AAHD_FLUDSCommonData>> result(work.size());
  opensn::ParallelFor(work.size(),
                      nthreads,
                      [&](size_t i)
                      {
                        result[i] = std::make_unique<AAHD_FLUDSCommonData>(
                          *work[i].spds, grid_nodal_mappings, discretization);
                      });
  for (size_t i = 0; i < work.size(); ++i)
    runtime.quadrature_fluds_commondata_map[work[i].quadrature].push_back(std::move(result[i]));
  log.Log() << program_timer.GetTimeString() << " FLUDS construction done.";
}

void
BuildCBCGPUFludsCommonData(SweepRuntime& runtime,
                           const SpatialDiscretization& discretization,
                           const std::vector<CellFaceNodalMapping>& grid_nodal_mappings)
{
  struct WorkItem
  {
    std::shared_ptr<AngularQuadrature> quadrature;
    std::shared_ptr<SPDS> spds;
  };
  std::vector<WorkItem> work;
  for (const auto& [quadrature, spds_list] : runtime.quadrature_spds_map)
    for (const auto& spds : spds_list)
      work.push_back({quadrature, spds});
  if (work.empty())
    return;

  const auto nthreads = std::min<std::size_t>(work.size(), std::max(1U, opensn_num_threads));
  log.Log() << program_timer.GetTimeString() << " FLUDS construction: " << work.size()
            << " angles (" << nthreads << " thread(s)).";
  std::vector<std::unique_ptr<CBCD_FLUDSCommonData>> result(work.size());
  ParallelFor(work.size(),
              nthreads,
              [&](size_t i)
              {
                result[i] = std::make_unique<CBCD_FLUDSCommonData>(
                  *work[i].spds, grid_nodal_mappings, discretization);
              });
  for (size_t i = 0; i < work.size(); ++i)
    runtime.quadrature_fluds_commondata_map[work[i].quadrature].push_back(std::move(result[i]));
  log.Log() << program_timer.GetTimeString() << " FLUDS construction done.";
}

} // namespace opensn::detail
