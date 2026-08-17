// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep_runtime_builder.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/aah_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbc_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aah_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aah_sweep_chunk_td.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_chunk_td.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/utils/error.h"
#include "framework/utils/caliper_scopes.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace opensn
{

SweepKind
ParseSweepKind(const std::string& sweep_type, const std::string& problem_name)
{
  if (sweep_type == "AAH")
    return SweepKind::AAH;
  if (sweep_type == "CBC")
    return SweepKind::CBC;
  OpenSnInvalidArgument(problem_name + ": Unsupported sweep type \"" + sweep_type + "\"");
}

void
DiscreteOrdinatesProblem::InitializeSweepDataStructures()
{
  CALI_CXX_MARK_SCOPE("SweepDataStructures");

  auto sweep_runtime = BuildSweepRuntime(
    GetName(), groupsets_, grid_, sweep_type_, use_gpus_, *discretization_, grid_nodal_mappings_);

  quadrature_unq_so_grouping_map_ = std::move(sweep_runtime.quadrature_unq_so_grouping_map);
  quadrature_spds_map_ = std::move(sweep_runtime.quadrature_spds_map);
  quadrature_fluds_commondata_map_ = std::move(sweep_runtime.quadrature_fluds_commondata_map);
}

#ifndef __OPENSN_WITH_GPU__
std::shared_ptr<FLUDS>
DiscreteOrdinatesProblem::CreateAAHD_FLUDS(unsigned int num_groups,
                                           std::size_t num_angles,
                                           const FLUDSCommonData& common_data)
{
  throw std::runtime_error(
    "DiscreteOrdinatesProblem::CreateAAHD_FLUDS : OPENSN_WITH_CUDA not enabled.");
  return {};
}

std::shared_ptr<AngleSet>
DiscreteOrdinatesProblem::CreateAAHD_AngleSet(
  size_t id,
  const LBSGroupset& groupset,
  const SPDS& spds,
  std::shared_ptr<FLUDS>& fluds,
  std::vector<size_t>& angle_indices,
  std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
  int maximum_message_size,
  const MPICommunicatorSet& in_comm_set)
{
  throw std::runtime_error(
    "DiscreteOrdinatesProblem::CreateAAHD_AngleSet : OPENSN_WITH_CUDA not enabled.");
  return {};
}

std::shared_ptr<SweepChunk>
DiscreteOrdinatesProblem::CreateAAHD_SweepChunk(LBSGroupset& groupset)
{
  throw std::runtime_error(
    "DiscreteOrdinatesProblem::CreateAAHD_SweepChunk : OPENSN_WITH_CUDA not enabled.");
  return {};
}

std::shared_ptr<FLUDS>
DiscreteOrdinatesProblem::CreateCBCD_FLUDS(std::size_t num_groups,
                                           std::size_t num_angles,
                                           std::size_t num_local_cells,
                                           const FLUDSCommonData& common_data,
                                           const UnknownManager& psi_uk_man,
                                           const SpatialDiscretization& sdm,
                                           bool save_angular_flux)
{
  throw std::runtime_error(
    "DiscreteOrdinatesProblem::CreateCBCD_FLUDS : OPENSN_WITH_CUDA not enabled.");
  return {};
}

std::shared_ptr<AngleSet>
DiscreteOrdinatesProblem::CreateCBCD_AngleSet(
  size_t id,
  const LBSGroupset& groupset,
  const SPDS& spds,
  std::shared_ptr<FLUDS>& fluds,
  std::vector<size_t>& angle_indices,
  std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
  const MPICommunicatorSet& in_comm_set)
{
  throw std::runtime_error(
    "DiscreteOrdinatesProblem::CreateCBCD_AngleSet : OPENSN_WITH_CUDA not enabled.");
  return {};
}

std::shared_ptr<SweepChunk>
DiscreteOrdinatesProblem::CreateCBCDSweepChunk(LBSGroupset& groupset)
{
  throw std::runtime_error(
    "DiscreteOrdinatesProblem::CreateCBCDSweepChunk : OPENSN_WITH_CUDA not enabled.");
  return {};
}
#endif

void
DiscreteOrdinatesProblem::InitFluxDataStructures(LBSGroupset& groupset)
{
  CALI_CXX_MARK_SCOPE("FluxDataStructures");

  const auto& quadrature_sweep_info = quadrature_unq_so_grouping_map_[groupset.quadrature];

  const auto& unique_so_groupings = quadrature_sweep_info.first;
  const auto& dir_id_to_so_map = quadrature_sweep_info.second;

  const size_t gs_num_grps = groupset.GetNumGroups();

  // Passing the sweep boundaries to the angle aggregation
  groupset.angle_agg =
    std::make_shared<AngleAggregation>(groupset, sweep_boundaries_, groupset.quadrature, grid_);

  size_t angle_set_id = 0;
  for (const auto& so_grouping : unique_so_groupings)
  {
    const size_t master_dir_id = so_grouping.front();
    const size_t so_id = dir_id_to_so_map.at(master_dir_id);

    const auto& sweep_ordering = quadrature_spds_map_[groupset.quadrature][so_id];
    const auto& fluds_common_data = *quadrature_fluds_commondata_map_[groupset.quadrature][so_id];

    std::vector<size_t> angle_indices(so_grouping.begin(), so_grouping.end());

    switch (ParseSweepKind(sweep_type_, GetName()))
    {
      case SweepKind::AAH:
      {
        std::shared_ptr<FLUDS> fluds;
        if (use_gpus_)
        {
          fluds = CreateAAHD_FLUDS(gs_num_grps, angle_indices.size(), fluds_common_data);
        }
        else
        {
          fluds = std::make_shared<AAH_FLUDS>(
            gs_num_grps,
            angle_indices.size(),
            dynamic_cast<const AAH_FLUDSCommonData&>(fluds_common_data));
        }

        std::shared_ptr<AngleSet> angle_set;
        if (use_gpus_)
        {
          angle_set = CreateAAHD_AngleSet(angle_set_id++,
                                          groupset,
                                          *sweep_ordering,
                                          fluds,
                                          angle_indices,
                                          sweep_boundaries_,
                                          options_.max_mpi_message_size,
                                          *grid_local_comm_set_);
        }
        else
        {
          angle_set = std::make_shared<AAH_AngleSet>(angle_set_id++,
                                                     groupset,
                                                     *sweep_ordering,
                                                     fluds,
                                                     angle_indices,
                                                     sweep_boundaries_,
                                                     options_.max_mpi_message_size,
                                                     *grid_local_comm_set_);
        }
        groupset.angle_agg->GetAngleSetGroups().push_back(angle_set);
        break;
      }
      case SweepKind::CBC:
      {
        std::shared_ptr<FLUDS> fluds;
        if (use_gpus_)
        {
          fluds = CreateCBCD_FLUDS(gs_num_grps,
                                   angle_indices.size(),
                                   grid_->local_cells.size(),
                                   fluds_common_data,
                                   groupset.psi_uk_man_,
                                   *discretization_,
                                   (not GetPsiNewLocal()[groupset.id].empty()));
        }
        else
        {
          fluds =
            std::make_shared<CBC_FLUDS>(gs_num_grps,
                                        angle_indices.size(),
                                        dynamic_cast<const CBC_FLUDSCommonData&>(fluds_common_data),
                                        groupset.psi_uk_man_,
                                        *discretization_);
        }

        std::shared_ptr<AngleSet> angle_set;
        if (use_gpus_)
        {
          angle_set = CreateCBCD_AngleSet(angle_set_id++,
                                          groupset,
                                          *sweep_ordering,
                                          fluds,
                                          angle_indices,
                                          sweep_boundaries_,
                                          *grid_local_comm_set_);
        }
        else
        {
          angle_set = std::make_shared<CBC_AngleSet>(angle_set_id++,
                                                     groupset,
                                                     *sweep_ordering,
                                                     fluds,
                                                     angle_indices,
                                                     sweep_boundaries_,
                                                     options_.max_mpi_message_size,
                                                     *grid_local_comm_set_);
        }

        groupset.angle_agg->GetAngleSetGroups().push_back(angle_set);
        break;
      }
    }
  } // for so_grouping

  groupset.angle_agg->BuildDirnumToAnglesetMap();

  opensn::mpi_comm.barrier();
}

std::shared_ptr<SweepChunk>
DiscreteOrdinatesProblem::SetSweepChunk(LBSGroupset& groupset)
{

  const auto mode = sweep_chunk_mode_.value_or(SweepChunkMode::DEFAULT);

  const bool use_time_dependent_chunk = (mode == SweepChunkMode::TIME_DEPENDENT);

  switch (ParseSweepKind(sweep_type_, GetName()))
  {
    case SweepKind::AAH:
      if (use_time_dependent_chunk)
        return std::make_shared<AAHSweepChunkTD>(*this, groupset);
      if (use_gpus_)
        return CreateAAHD_SweepChunk(groupset);
      return std::make_shared<AAHSweepChunk>(*this, groupset);
    case SweepKind::CBC:
      if (use_time_dependent_chunk)
        return std::make_shared<CBCSweepChunkTD>(*this, groupset);
      if (use_gpus_)
        return CreateCBCDSweepChunk(groupset);
      return std::make_shared<CBCSweepChunk>(*this, groupset);
  }
  OpenSnLogicalError("Unreachable: unknown SweepKind.");
}

} // namespace opensn
