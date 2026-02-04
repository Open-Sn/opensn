// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/scheduler/sweep_scheduler.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/aah.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/math/quadratures/angular/curvilinear_product_quadrature.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <algorithm>
#include <unordered_map>

namespace opensn
{

SweepScheduler::SweepScheduler(SchedulingAlgorithm scheduler_type,
                               AngleAggregation& angle_agg,
                               SweepChunk& sweep_chunk)
  : scheduler_type_(scheduler_type), angle_agg_(angle_agg), sweep_chunk_(sweep_chunk)
{
  CALI_CXX_MARK_SCOPE("SweepScheduler::SweepScheduler");

  angle_agg_.InitializeReflectingBCs();

  if (scheduler_type_ == SchedulingAlgorithm::DEPTH_OF_GRAPH)
    InitializeAlgoDOG();

  if (scheduler_type_ == SchedulingAlgorithm::ALL_AT_ONCE ||
      scheduler_type_ == SchedulingAlgorithm::DEPTH_OF_GRAPH)
  {
    angle_agg_.SetupAngleSetDependencies();
  }

  // Initialize delayed upstream data
  for (auto& angset : angle_agg_)
    angset->InitializeDelayedUpstreamData();

  if (scheduler_type_ == SchedulingAlgorithm::DEPTH_OF_GRAPH)
  {
    const auto* curvi_quad =
      dynamic_cast<const CurvilinearProductQuadrature*>(angle_agg_.GetQuadrature().get());
    if (curvi_quad && angle_agg_.GetQuadrature()->GetDimension() == 2 &&
        angle_agg_.GetCoordinateSystem() == CoordinateSystemType::CYLINDRICAL)
    {
      const auto& following_map = angle_agg_.GetFollowingAngleSetsMap();
      for (const auto& [from, to_set] : following_map)
        for (auto* to : to_set)
          preceding_angle_sets_[to].insert(from);
    }
  }

  // Get local max num messages accross anglesets
  int local_max_num_messages = 0;
  for (auto& angset : angle_agg_)
    local_max_num_messages = std::max(angset->GetMaxBufferMessages(), local_max_num_messages);

  // Reconcile all local maximums
  int global_max_num_messages = 0;
  mpi_comm.all_reduce(local_max_num_messages, global_max_num_messages, mpi::op::max<int>());

  // Propogate items back to sweep buffers
  for (auto& angset : angle_agg_)
    angset->SetMaxBufferMessages(global_max_num_messages);
}

void
SweepScheduler::InitializeAlgoDOG()
{
  CALI_CXX_MARK_SCOPE("SweepScheduler::InitializeAlgoDOG");

  const bool is_cylindrical = angle_agg_.GetQuadrature()->GetDimension() == 2 &&
                              angle_agg_.GetCoordinateSystem() == CoordinateSystemType::CYLINDRICAL;

  std::unordered_map<unsigned int, int> angle_order;
  const auto* curvi_quad =
    dynamic_cast<const CurvilinearProductQuadrature*>(angle_agg_.GetQuadrature().get());
  const auto* product_quad =
    dynamic_cast<const ProductQuadrature*>(angle_agg_.GetQuadrature().get());
  if (is_cylindrical && curvi_quad && product_quad)
  {
    int order = 0;
    for (const auto& dir_set : product_quad->GetDirectionMap())
      for (const auto dir_id : dir_set.second)
        angle_order.emplace(dir_id, order++);
  }
  // Load all anglesets in preperation for sorting
  size_t num_anglesets = angle_agg_.GetNumAngleSets();
  for (size_t as = 0; as < num_anglesets; ++as)
  {
    auto angleset = angle_agg_[as];
    const auto& spds = dynamic_cast<const AAH_SPDS&>(angleset->GetSPDS());

    const std::vector<STDG>& leveled_graph = spds.GetGlobalSweepPlanes();

    // Find location depth
    int loc_depth = -1;
    for (size_t level = 0; level < leveled_graph.size(); ++level)
    {
      for (size_t index = 0; index < leveled_graph[level].item_id.size(); ++index)
      {
        if (leveled_graph[level].item_id[index] == opensn::mpi_comm.rank())
        {
          loc_depth = static_cast<int>(leveled_graph.size() - level);
          break;
        }
      } // for locations in plane
    } // for sweep planes

    // Set up rule values
    if (loc_depth >= 0)
    {
      RuleValues new_rule_vals(angleset);
      new_rule_vals.depth_of_graph = loc_depth;
      new_rule_vals.set_index = as;

      const auto& omega = spds.GetOmega();
      new_rule_vals.sign_of_omegax = (omega.x >= 0) ? 2 : 1;
      new_rule_vals.sign_of_omegay = (omega.y >= 0) ? 2 : 1;
      new_rule_vals.sign_of_omegaz = (omega.z >= 0) ? 2 : 1;
      if (is_cylindrical && !angle_order.empty() && angleset->GetNumAngles() == 1)
      {
        const auto angle_idx = angleset->GetAngleIndices().front();
        const auto it = angle_order.find(angle_idx);
        if (it != angle_order.end())
          new_rule_vals.azimuthal_order = it->second;
      }

      rule_values_.push_back(new_rule_vals);
    }
    else
      throw std::runtime_error("InitializeAlgoDOG: Failed to find location depth");
  } // for anglesets

  if (is_cylindrical)
    SortRuleValuesDOGRZ();
  else
    SortRuleValuesDOGDefault();
}

void
SweepScheduler::SortRuleValuesDOGRZ()
{
  std::stable_sort(rule_values_.begin(),
                   rule_values_.end(),
                   [](const RuleValues& a, const RuleValues& b)
                   {
                     if (a.sign_of_omegax != b.sign_of_omegax)
                       return a.sign_of_omegax > b.sign_of_omegax;
                     if (a.depth_of_graph != b.depth_of_graph)
                       return a.depth_of_graph > b.depth_of_graph;
                     if (a.sign_of_omegax != b.sign_of_omegax)
                       return a.sign_of_omegax > b.sign_of_omegax;
                     if (a.sign_of_omegay != b.sign_of_omegay)
                       return a.sign_of_omegay > b.sign_of_omegay;
                     if (a.sign_of_omegaz != b.sign_of_omegaz)
                       return a.sign_of_omegaz > b.sign_of_omegaz;
                     if (a.azimuthal_order != b.azimuthal_order)
                       return a.azimuthal_order < b.azimuthal_order;
                     return a.set_index < b.set_index;
                   });
}

void
SweepScheduler::SortRuleValuesDOGDefault()
{
  std::stable_sort(rule_values_.begin(),
                   rule_values_.end(),
                   [](const RuleValues& a, const RuleValues& b)
                   {
                     if (a.depth_of_graph != b.depth_of_graph)
                       return a.depth_of_graph > b.depth_of_graph;
                     if (a.sign_of_omegax != b.sign_of_omegax)
                       return a.sign_of_omegax > b.sign_of_omegax;
                     if (a.sign_of_omegay != b.sign_of_omegay)
                       return a.sign_of_omegay > b.sign_of_omegay;
                     return a.sign_of_omegaz > b.sign_of_omegaz;
                   });
}

void
SweepScheduler::ScheduleAlgoDOG(SweepChunk& sweep_chunk)
{
  CALI_CXX_MARK_SCOPE("SweepScheduler::ScheduleAlgoDOG");

  const bool is_cylindrical = angle_agg_.GetQuadrature()->GetDimension() == 2 &&
                              angle_agg_.GetCoordinateSystem() == CoordinateSystemType::CYLINDRICAL;

  if (is_cylindrical)
    ScheduleAlgoDOGRZ(sweep_chunk);
  else
    ScheduleAlgoDOGDefault(sweep_chunk);

  // Receive delayed data
  opensn::mpi_comm.barrier();
  bool received_delayed_data = false;
  while (not received_delayed_data)
  {
    received_delayed_data = true;

    for (auto& angle_set : angle_agg_)
    {
      if (angle_set->FlushSendBuffers() == AngleSetStatus::MESSAGES_PENDING)
        received_delayed_data = false;

      if (not angle_set->ReceiveDelayedData())
        received_delayed_data = false;
    }
  }

  // Reset all
  for (auto& angle_set : angle_agg_)
    angle_set->ResetSweepBuffers();

  for (const auto& [bid, bndry] : angle_agg_.GetSimBoundaries())
    bndry->ResetAnglesReadyStatus();
}

void
SweepScheduler::ScheduleAlgoDOGRZ(SweepChunk& sweep_chunk)
{
  bool finished = false;
  std::unordered_set<AngleSet*> completed;
  while (not finished)
  {
    finished = true;
    for (auto& rule_value : rule_values_)
    {
      auto angleset = rule_value.angle_set;

      auto dep_it = preceding_angle_sets_.find(angleset.get());
      if (dep_it != preceding_angle_sets_.end())
      {
        bool deps_ready = true;
        for (auto* dep : dep_it->second)
        {
          if (completed.find(dep) == completed.end())
          {
            deps_ready = false;
            break;
          }
        }
        if (!deps_ready)
        {
          AngleSetStatus status =
            angleset->AngleSetAdvance(sweep_chunk, AngleSetStatus::READY_TO_EXECUTE);
          if (status != AngleSetStatus::FINISHED)
            finished = false;
          else
            completed.insert(angleset.get());
          continue;
        }
      }

      AngleSetStatus status = angleset->AngleSetAdvance(sweep_chunk, AngleSetStatus::EXECUTE);
      if (status == AngleSetStatus::FINISHED)
        completed.insert(angleset.get());
      if (status != AngleSetStatus::FINISHED)
        finished = false;
    }
  }
}

void
SweepScheduler::ScheduleAlgoDOGDefault(SweepChunk& sweep_chunk)
{
  bool finished = false;
  while (not finished)
  {
    finished = true;
    for (auto& rule_value : rule_values_)
    {
      auto angleset = rule_value.angle_set;
      AngleSetStatus status = angleset->AngleSetAdvance(sweep_chunk, AngleSetStatus::EXECUTE);
      if (status != AngleSetStatus::FINISHED)
        finished = false;
    }
  }
}

void
SweepScheduler::ScheduleAlgoFIFO(SweepChunk& sweep_chunk)
{
  CALI_CXX_MARK_SCOPE("SweepScheduler::ScheduleAlgoFIFO");

  // Loop over AngleSetGroups
  bool finished = false;
  while (not finished)
  {
    finished = true;

    for (auto& angle_set : angle_agg_)
    {
      AngleSetStatus status = angle_set->AngleSetAdvance(sweep_chunk, AngleSetStatus::EXECUTE);
      if (status != AngleSetStatus::FINISHED)
        finished = false;
    } // for angleset
  } // while not finished

  // Receive delayed data
  opensn::mpi_comm.barrier();
  bool received_delayed_data = false;
  while (not received_delayed_data)
  {
    received_delayed_data = true;

    for (auto& angle_set : angle_agg_)
    {
      if (angle_set->FlushSendBuffers() == AngleSetStatus::MESSAGES_PENDING)
        received_delayed_data = false;

      if (not angle_set->ReceiveDelayedData())
        received_delayed_data = false;
    }
  }

  // Reset all
  for (auto& angle_set : angle_agg_)
    angle_set->ResetSweepBuffers();

  for (const auto& [bid, bndry] : angle_agg_.GetSimBoundaries())
    bndry->ResetAnglesReadyStatus();
}

#ifndef __OPENSN_WITH_GPU__

void
SweepScheduler::ScheduleAlgoAAO(SweepChunk& sweep_chunk)
{
  throw std::runtime_error("SweepScheduler::ScheduleAlgoAAO: AAO scheduling is only "
                           "available for builds with GPU support.");
}

#endif // __OPENSN_WITH_GPU__

void
SweepScheduler::Sweep()
{
  CALI_CXX_MARK_SCOPE("SweepScheduler::Sweep");

  if (scheduler_type_ == SchedulingAlgorithm::FIRST_IN_FIRST_OUT)
    ScheduleAlgoFIFO(sweep_chunk_);
  else if (scheduler_type_ == SchedulingAlgorithm::ALL_AT_ONCE)
    ScheduleAlgoAAO(sweep_chunk_);
  else if (scheduler_type_ == SchedulingAlgorithm::DEPTH_OF_GRAPH)
    ScheduleAlgoDOG(sweep_chunk_);
}

void
SweepScheduler::PrepareForSweep(bool use_boundary_source, bool zero_incoming_delayed_psi)
{
  if (zero_incoming_delayed_psi)
    angle_agg_.ZeroIncomingDelayedPsi();
  angle_agg_.ZeroOutgoingDelayedPsi();
  sweep_chunk_.ZeroDestinationPsi();
  sweep_chunk_.ZeroDestinationPhi();
  sweep_chunk_.SetBoundarySourceActiveFlag(use_boundary_source);
}

} // namespace opensn
