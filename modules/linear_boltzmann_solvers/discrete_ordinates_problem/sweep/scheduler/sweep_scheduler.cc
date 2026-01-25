// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/scheduler/sweep_scheduler.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/aah.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <algorithm>
#include <sstream>

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

  if (scheduler_type_ == SchedulingAlgorithm::ALL_AT_ONCE)
    angle_agg_.SetupAngleSetDependencies();

  // Initialize delayed upstream data
  for (auto& angset : angle_agg_)
    angset->InitializeDelayedUpstreamData();

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

      rule_values_.push_back(new_rule_vals);
    }
    else
      throw std::runtime_error("InitializeAlgoDOG: Failed to find location depth");
  } // for anglesets

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

  // Loop till done
  bool finished = false;
  while (not finished)
  {
    finished = true;
    for (auto& rule_value : rule_values_)
    {
      auto angleset = rule_value.angle_set;

      // Query angleset status
      // Angleset status will be one of the following:
      //  - RECEIVING.
      //      Meaning it is either waiting for messages or actively receiving it
      //  - READY_TO_EXECUTE.
      //      Meaning it has received all upstream data and can be executed
      //  - FINISHED.
      //      Meaning the angleset has executed its sweep chunk
      AngleSetStatus status = angleset->AngleSetAdvance(sweep_chunk, AngleSetStatus::EXECUTE);

      if (status != AngleSetStatus::FINISHED)
        finished = false;
    } // for each angleset rule
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
